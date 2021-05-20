import os
import math
import datetime
import tqdm
import urllib
import yaml
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.coordinates import get_sun, AltAz, EarthLocation
from astropy.time import Time
import rasterio as rs
from rasterio.errors import RasterioIOError
from swissreframe import initialize_reframe
from itertools import groupby
import logging
import requests

requests.packages.urllib3.disable_warnings()

logging.basicConfig(filename="log.log",
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.INFO)

logging.info("Started siwssSUN")


def valid_datetime(s):
    try:
        return datetime.datetime.strptime(s, "%Y-%m-%d %H:%M")
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise argparse.ArgumentTypeError(msg)

def valid_date(s):
    try:
        return datetime.datetime.strptime(s, "%Y-%m-%d")
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise argparse.ArgumentTypeError(msg)

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

def get_dict_of_tiffs(path):
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
    cache_dir = config.get("cache_folder")
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    filename = url.split("/")[-1]
    file_path = os.path.join(cache_dir, filename)
    if not os.path.isfile(file_path):
        with open(file_path, 'wb') as f:
            resp = requests.get(url, verify=False)
            f.write(resp.content)
        #urllib.request.urlretrieve(url, file_path)
    
    return file_path

def get_url_from_coordinates(lat, lon):
    ch_lon, ch_lat, _ = WGS84_to_LV95(lat, lon, 0)
    url = dict_of_tiffs.get(int(ch_lat/1000), dict()).get(int(ch_lon/1000))
    return url

def get_height_multi(coords, url):
    '''
    extracts for lat, lon coordinates the corresponding height from
    SWISSALTI3D Model
    '''
    ch_coords = [WGS84_to_LV95(lat, lon, 0) for lat, lon in coords]
    dataset = rs.open(cache_geotif(url))
    try:
        result = [dataset.read(1)[dataset.index(ch_lon, ch_lat)] for ch_lon, ch_lat, _ in ch_coords]
    except RasterioIOError as ex:
        logging.error(
            " ".join(
                ["RasterIO Error while processing", 
                "(",
                str(coords[0][0]),
                ",",
                str(coords[0][1]),
                ")"]
            ))
        return([0 for _ in ch_coords])
        
    dataset.close()
    return result

def get_height(lat, lon):
    '''
    extracts for lat, lon coordinates the corresponding height from
    SWISSALTI3D Model
    '''
    ch_lon, ch_lat, _ = WGS84_to_LV95(lat, lon, 0)
    url = get_url_from_coordinates(lat, lon)
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

def spheric_from_lat_lon(lat, lon, rad=1):
    '''
    returns spherical coordinates given lat, lon, [radius] 
    '''
    return np.array([rad, (90-lat)/180*math.pi, lon/180*math.pi])

def lat_lon_from_spheric(x0):
    '''
    returns lat, lon for an input x0 = (radius, polar, azimuth)
    '''
    return -x0[1]*180/math.pi+90, x0[2]*180/math.pi

def spheric_to_cartesian(x0):
    x = x0[0]*math.cos(x0[2])*math.sin(x0[1])
    y = x0[0]*math.sin(x0[2])*math.sin(x0[1])
    z = x0[0]*math.cos(x0[1])
    return np.array([x, y, z])

def cartesian_to_spheric(x0):
    r = math.sqrt(np.sum(x0**2))
    theta = math.acos(x0[2]/r)
    phi = math.atan2(x0[1], x0[0])
    return np.array([r, theta, phi])

def rotation_matrix_x(theta):
    return np.array(
        [[1,0,0],
        [0, math.cos(theta), -math.sin(theta)],
        [0, math.sin(theta), math.cos(theta)]]
    )

def rotation_matrix_z(phi):
    return np.array(
        [[math.cos(phi),-math.sin(phi),0],
        [math.sin(phi), math.cos(phi), 0],
        [0, 0, 1]]
    )

def reparameterize_sphere(x0):
    '''
    takes x no as input and returns two functions, which reparameterizes any 
    input point in spherical coordinates into a coordinate system where x0
    is the north pole. The first function transforms into the space, where
    x0 is the north pole, the second function transforms back.
    '''
    theta = x0[1]-0.00000000001
    phi = x0[2]-math.pi/2
    f1 = lambda x: cartesian_to_spheric(
        rotation_matrix_x(theta) @ rotation_matrix_z(-phi) @ spheric_to_cartesian(x)
    )
    f2 = lambda x: cartesian_to_spheric(
        rotation_matrix_z(phi) @ rotation_matrix_x(-theta) @ spheric_to_cartesian(x)
    )
    return f1, f2

def get_great_circle_at_x0(lat, lon, azimuth):
    x0 = spheric_from_lat_lon(lat, lon)
    f1, f2 = reparameterize_sphere(x0)
    great_circle = lambda t: lat_lon_from_spheric(f2(f1(x0)+np.array([0, t*2*math.pi, azimuth/180*math.pi])))
    return great_circle

def get_angle_between(lat1, lon1, lat2, lon2):
    '''
    for given coordinates (lat1, lon1) and (lat2, lon2) the function
    returns the angle between these points at the earth middle point.
    '''
    return math.acos(np.dot(get_cartesian(lat1, lon1), get_cartesian(lat2, lon2)))/math.pi*180

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
    sea_level = float(config.get("earth_radius"))
    altitude_rad = altitude/180*math.pi
    angle_rad = angle/180*math.pi
    return (r0+sea_level)/math.sin(math.pi/2-altitude_rad-angle_rad)*math.sin(math.pi/2+altitude_rad)-sea_level

def get_distance_offset(i, delta=50):
    '''
    just a function to get quadratic increasing distant values
    '''
    return i*1/float(config.get("earth_radius"))/2/math.pi*delta

def get_location_on_great_circle(great_circle, t, h0, altitude):
    lat, lon = great_circle(t)
    return {
        "lat": lat,
        "lon": lon,
        "max_height": get_max_height(h0, altitude, t*360),
        "alpha": t*360,
        "url": get_url_from_coordinates(lat, lon)
    }

def sun_visibility_check(lat, lon, height, date, silent=False, verbose=False):
    solar = get_solar_radiation(lat,lon,date)
    great_circle = get_great_circle_at_x0(lat, lon, solar.get("azimuth"))

    list_of_points = [get_location_on_great_circle(great_circle, get_distance_offset(t), height, solar.get("altitude"))
                   for t in range(5, 3000)]
    grouped = groupby(list_of_points, lambda x: x["url"])
    for k, v in tqdm.tqdm(grouped, disable=not verbose):
        vv = list(v)
        if k is None:
            #raise Warning("Point Lat/Lon {}, {} is not in Switzerland, be careful".format(vv[0].get("lat"), vv[0].get("lon")))
            logging.warning("Point Lat/Lon {}, {} is not in Switzerland, be careful".format(vv[0].get("lat"), vv[0].get("lon")))
            continue
        heights = np.array(get_height_multi([(ll["lat"], ll["lon"]) for ll in vv], k))
        max_heights = np.array([ll["max_height"] for ll in vv])
        if min(max_heights-heights)<0:
            if not silent:
                idx = (np.argwhere(max_heights-heights<0))[0][0]
                print("Lat/Lon: {}, {}   , height: {}  max_height:{}".format(vv[idx].get("lat"), vv[idx].get("lon"), heights[idx], vv[idx].get("max_height")))
        
            return False
        if min(max_heights-heights)>5000:
            break

    if verbose and not silent:
        print("The sun should be visible. Checked for obstacles in a range of {:.2f} km".format(list_of_points[-1].get("alpha")/180000*math.pi*float(config.get("earth_radius"))))
    elif not silent:
        print("The sun should be visible")
    
    return True

if __name__ == "__main__":
    '''
    Main idea of the script:s

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

    #load config file
    config = yaml.load(open('config.yml', 'r'), Loader=yaml.BaseLoader)

    #initialize the module to load GEO TIFFS
    r = initialize_reframe()

    #get the dict of all GEO TIFF URLS
    dict_of_tiffs = get_dict_of_tiffs(path = config.get("geo_tiff_links"))

    #Argument parsing
    parser = argparse.ArgumentParser(description='Determines wheter a point in Switzerland is sunny or not')
    parser.add_argument('lat', metavar='LATITUDE', type=float, nargs=None,
                        help='The latitude of the point of interest')
    parser.add_argument('lon', metavar='LONGITUDE', type=float, nargs=None,
                        help='The latitude of the point of interest')
    parser.add_argument("-r", 
                        "--height", 
                        help="The height from sea level. Default is the ground level +2 meters.", 
                        required=False, 
                        type=float)
    parser.add_argument("-d", 
                        "--date", 
                        help="The datetime in UTC in the format 'YYYY-MM-DD HH:MM' to request \
                        whether the sun in visible. The default value is now.", 
                        required=False, 
                        type=valid_datetime, default=datetime.datetime.now())
    parser.add_argument("-c", 
                        "--complete", 
                        help="The date in format 'YYYY-MM-DD' to request sunset and sunrise", 
                        required=False, 
                        type=valid_date)
    parser.add_argument('-v', '--verbose', action='store_true', help='Shows progressbar if set')
    args = parser.parse_args()
    date = args.date
    lat = args.lat
    lon = args.lon


    #If no height is set, then take the ground height + 2 meters
    if not args.height:
        height = get_height(lat, lon)+2
    else:
        height = args.height


    if args.complete:
        dt = args.complete
        from_dt =  dt.replace(hour=2, minute=0, second=0, microsecond=0)
        to_dt = dt+datetime.timedelta(hours=21)
        dates = pd.date_range(from_dt, to_dt, freq="10min")
        memory = False
        for dd in tqdm.tqdm(dates, disable=not args.verbose):
            if not memory and sun_visibility_check(lat, lon, height, dd, silent=True, verbose=False):
                print("Sunrise at: {}".format(dd))
                memory = True
            
            if memory and not sun_visibility_check(lat, lon, height, dd, silent=True, verbose=False):
                print("Sunset at: {}".format(dd))
                memory = False
    else:
        sun_visibility_check(lat, lon, height, date, verbose=args.verbose)