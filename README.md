# swissSUN

This is a simple piece of code to compute whether the sun is visible in Switzerland at a specific location/time or not. To install run the following code:
```
git clone https://github.com/leokster/swissSUN.git
cd swissSUN
python -m venv venv
source ./venv/bin/activate
pip install -r requirements.txt
```

To run the code, activate the venv
```
source ./venv/bin/activate
```

To compute if the sun is now visible at (47.36514, 8.54004) run
```
python sun_calc.py 47.36514 8.54004
```

To compute if the sun is visible at (47.36514, 8.54004) on March 2nd 2021, 16:00 run
```
python sun_calc.py 47.36514 8.54004 -d "2021-03-02 16:00"
```

To compute if the sun is visible at (47.36514, 8.54004) on March 2nd 2021, 16:00 1000m above sea ground run
```
python sun_calc.py 47.36514 8.54004 -d "2021-03-02 16:00" -r 1000
```

To compute the sunrise and sunset at (47.36514, 8.54004) on March 2nd 2021 run
```
python sun_calc.py 47.36514 8.54004 -c 2021-03-02
```

The verbose flag ```-v``` or ```--verbose``` can be set to show progress bars during calculation.



## Remark
The code might take a while the first few times, since it has to download and cache the GEO TIFFS


## Data Source
The model works only in Switzerland and uses the [swissALTI3D](https://www.swisstopo.admin.ch/de/geodata/height/alti3d.html) topographical model of Switzerland. To get the sun position we use the Python package [astropy](https://www.astropy.org). A recent version of the GEO TIFFS link list can be downloaded [here](https://ogd.swisstopo.admin.ch/resources/ch.swisstopo.swissalti3d-2eb0dsEH.csv).