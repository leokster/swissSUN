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




## Remark
The code might take a while the first few times, since it has to download and cache the GEO TIFFS