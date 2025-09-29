from netCDF4 import Dataset
import numpy as np
from pathlib import Path
import pandas as pd
from scipy.ndimage import generic_filter

data_dir = Path('/discover/nobackup/asouri/MODELS/WPS/')
dates = pd.date_range('2023-10-01', '2024-01-01', freq='3H')

for dt in dates:
    fname = data_dir / f"met_em.d01.{dt.strftime('%Y-%m-%d_%H:%M:%S')}.nc"
    if not fname.exists():
        print(f"Missing {fname}")
        continue

    with Dataset(fname, 'r+') as ds:  # r+ allows overwrite
        if 'SKINTEMP' not in ds.variables:
            print("No TSK in this file")
            continue

        tsk_var = ds.variables['SKINTEMP']
        arr = tsk_var[:]  # read all at once (time,y,x) or (y,x)
        arr[0,440:460,570:591] = 273.01
        # write back in place
        tsk_var[:] = arr
        print(f"Fixed {fname}")
