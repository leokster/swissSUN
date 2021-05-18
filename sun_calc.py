import os
import numpy as np
import pandas as pd
import math
import datetime
import tqdm
from pathlib import Path
import urllib

from pysolar import radiation
from astropy.coordinates import get_sun, AltAz, EarthLocation
from astropy.time import Time

import rasterio as rs
from rasterio.plot import show

from pyrocko import orthodrome
from pysolar.solar import get_altitude, get_azimuth
from pysolar.radiation import get_radiation_direct

from swissreframe import Coordinate, initialize_reframe

r = initialize_reframe()

def get_dict_of_tiffs(path="./geo_tiff_links.csv"):
    '''
    returns nested dict of such that 
    dict[lat][lon] = url to the geo_tiff
    '''
    trans = lambda x: [int(xx) for xx in x.split("/")[4].split("_")[-1].split("-")]
    with open(path, "r") as f:
        list_of_tiffs = f.read().split('\n')

    dict_of_tiffs = dict()
    for ll in list_of_tiffs:
        keys = trans(ll)
        if dict_of_tiffs.get(keys[1]) is None:
            dict_of_tiffs[keys[1]] = dict()

        dict_of_tiffs[keys[1]][keys[0]] = ll
    
    return dict_of_tiffs

def cache_geotif(url):
    '''
    returns the cached version of url
    '''
    cache_dir = "./geotif_cache"
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    filename = url.split("/")[-1]
    file_path = os.path.join(cache_dir, filename)
    if not os.path.isfile(file_path):
        urllib.request.urlretrieve(url, file_path)
    
    return file_path
    
def get_height(lat, lon):
    '''
    extracts for lat, lon coordinates the corresponding height from
    SWISSALTI3D Model
    '''
    ch_lon, ch_lat, _ = WGS84_to_LV95(lat, lon, 0)
    url = dict_of_tiffs.get(int(ch_lat/1000)).get(int(ch_lon/1000))
    dataset = rs.open(cache_geotif(url))
    return dataset.read(1)[dataset.index(ch_lon, ch_lat)]

def get_cartesian(lat, lon):
    '''
    converts (lat, lon) to cartesian coordinates on a unit sphere
    '''
    lat_rad = lat/180*math.pi
    lon_rad = lon/180*math.pi
    return np.array(
        [
            math.cos(lat_rad)*math.cos(lon_rad),
            math.cos(lat_rad)*math.sin(lon_rad),
            math.sin(lat_rad)
        ]
    )

def get_angle_between(lat1, lon1, lat2, lon2):
    '''
    for given coordinates (lat1, lon1) and (lat2, lon2) the function
    returns the angle between these points at the earth middle point.
    '''
    return math.acos(np.dot(get_cartesian(lat1, lon1), get_cartesian(lat2, lon2)))/math.pi*180

def WGS84_to_LV95(lat, lon, height):
    '''
    converts coordinates from the WGS84 system into LV95 Swiss coordinate
    system 
    '''
    return r.compute_gpsref((lon, lat, height), 'etrf93_gepgraphic_to_lv95')

def LV95_to_WGS84(lat, lon, height):
    '''
    converts coordinates from the LV95 Swiss coordinate
    system into the WGS84 system
    '''
    return r.compute_gpsref((lon, lat, height), 'lv95_to_etrf93_geographic')

def add_offset_to_coord(lat, lon, azimuth, distance):
    '''
    for given coordinates (lat, lon) the function returns new coordinates (lat, lon)
    given by "walking" the distance in the azimuth direction
    '''
    azimuth_rad = azimuth/180*math.pi
    east = np.sin(azimuth_rad)*distance
    north = np.cos(azimuth_rad)*distance
    return dict(zip(["lat", "lon"], orthodrome.ne_to_latlon(lat, lon, north, east)))


def get_solar_radiation(lat, lon, date):
    '''
    get all solar properties for a location (lat, lon) on a date
    '''
    sun_time = Time(date)
    loc = EarthLocation(lon=lon, lat=lat)
    altaz = AltAz(obstime=sun_time, location=loc)
    altitude_deg = get_sun(sun_time).transform_to(altaz).alt.value
    return {
        "azimuth": get_sun(sun_time).transform_to(altaz).az.value,
        "altitude": altitude_deg,
        "radiation": 0
    }

def get_max_height(r0, altitude, angle):
    '''
    for a starting point r0 (distance from sea level) and an angle discribing a second
    point on earth on the orthodrome (at r0 and into sun direction) the function returns
    the maximal height this second point on earth is allowed to have such that it is not 
    covering r0 from the sun
    '''
    sea_level = 6371146
    altitude_rad = altitude/180*math.pi
    angle_rad = angle/180*math.pi
    return (r0+sea_level)/math.sin(math.pi/2-altitude_rad-angle_rad)*math.sin(math.pi/2+altitude_rad)-sea_level

def get_distance_offset(i):
    '''
    just a function to get quadratic increasing distant values
    '''
    return i**2/75
        

def start_list_from_coord(lat, lon, height, date):
    '''
    initializing a list of point on the earth. the first point is the one for which
    we want to know whether there is sun or not. 
    '''
    return [{**{
        "lat": lat,
        "lon": lon, 
        "height": height,
        "date": date,
    }, **get_solar_radiation(lat,lon,date)}]

def get_next_point(list_of_points, distance):
    '''
    operates on a list. it takes the first element as a base point (the point for which we
    want to get information whether it gets sun or not) 
    '''
    lat, lon = add_offset_to_coord(
        list_of_points[-1].get("lat"),
        list_of_points[-1].get("lon"),
        list_of_points[-1].get("azimuth"),
        distance).values()
    alpha = get_angle_between(list_of_points[0]["lat"], list_of_points[0]["lon"], lat, lon)
    new_point = {**{
        "lat": lat,
        "lon": lon, 
        "height": get_height(lat, lon),
        "date": list_of_points[0]["date"],
        "max_height": get_max_height(list_of_points[0]["height"], list_of_points[0]["altitude"], alpha),
    }, **get_solar_radiation(lat,lon,list_of_points[0]["date"])}
    list_of_points.append(new_point)

def sun_visibility_check(lat, lon, height, date):
    '''
    run full visibility check for the coordinates (lat, lon, height) at date
    returns a dataframe describing the results. the first row is the point of interest.
    the following rows describe potential obstacles. if max_height > height, the corresponding
    obstacle-candidate is not covering the sun. otherwise it is covering the sun. 
    '''

    list_of_points = start_list_from_coord(lat,lon, height, date)

    for i in tqdm.tqdm(range(5, 200)):
        get_next_point(list_of_points, get_distance_offset(i))

        if list_of_points[-1]["max_height"]+2<list_of_points[-1]["height"]:
            break

        if list_of_points[-1]["max_height"]> 4000:
            break

    results = pd.DataFrame(list_of_points)
    return results

if __name__ == "__main__":
    '''
    Main idea of the script:

    For a given location on earth, described by x0 = (lat, lon, height) and a date we want to estimate whether
    the sun is visible at x0 at date or not. To achieve this we need to know where potential
    obstacles can be. All potential obstacles are on a so called great circle or orthodrome (these are circles on the
    earth sphere which center is the earth center) through the point x0 and the point where a straight line from 
    the sun to the center of earth intersects the earth sphere. The same great circle (or at least a partial arc of it) 
    can be found by walking infinitesimal steps in the azimuthal direction of the sun. Given this arc, we have an idea
    of where potential obstacles can be. With the swissalti3d model from swisstopo we then evaluate the height of each 
    obstacle-candidate. The position of the obstacle-candidate together with the altitude angle of the sun at x0 gives
    also a maximal height the obstacle is allowed to have to not cover the sun. 
    '''
    dict_of_tiffs = get_dict_of_tiffs()

    #the date / time is in UTC
    date = datetime.datetime(2021, 5, 17, 18, 10, 1, 130320, tzinfo=datetime.timezone.utc)
    lat = 47.34372
    lon = 8.53095
    height = get_height(lat, lon)

    results = sun_visibility_check(lat, lon, height, date)

    if (results["max_height"]-results["height"]).min()>0:
        print("congrats you see the sun")
    else:
        print("unfortunately there is no sun anymore. The sun is covered by Lat/Lon: {}, {}".format(results.tail(1)["lat"].values[0], results.tail(1)["lon"].values[0]))