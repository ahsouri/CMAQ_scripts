from geoscf2bc.drivers import default
import numpy as np
import datetime
import warnings
import time
from joblib import Parallel, delayed


warnings.filterwarnings("ignore", category=RuntimeWarning)


def _daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + datetime.timedelta(n)

def run_bcon(date1,date2):
    try:
       bcpaths = default(
         GDNAM='CONUS_8km', gdpath='./GRIDDESC',
         SDATE=f'{date1}T00', EDATE=f'{date2}T00',m3path ='/discover/nobackup/asouri/MODELS/CMAQv5.5/data/mcip/CONUS_8km/METCRO3D_CONUS_8km_20240512.nc'
         )
    except:
       print("oopsy!")

    return None

datarange = _daterange(datetime.date(2024, 5, 1), datetime.date(2024, 6, 1))
datarange = list(datarange)
output_files = []
out = Parallel(n_jobs=6,verbose=10)(delayed(run_bcon)(
           datarange[k].strftime('%Y%m%d'),(datarange[k] + datetime.timedelta(days=1)).strftime('%Y%m%d')) for k in range(len(datarange)))

#icpaths = default(
#    GDNAM='CONUS_8km', gdpath='./GRIDDESC', ftype=1,
#    SDATE='2023-06-01T00', EDATE='2023-06-01T00',m3path ='/discover/nobackup/asouri/MODELS/CMAQv5.5/data/mcip/CONUS_8km/METCRO3D_CONUS_8km_20231030.nc'
#)
