from netCDF4 import Dataset
import os
import numpy as np
import scipy.io as sio
import netCDF4 as nc
import warnings
from glob import glob
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.interpolate import RBFInterpolator
from datetime import datetime,timedelta
import copy
import re
from scipy.io import savemat

warnings.filterwarnings("ignore", category=RuntimeWarning)

def _read_nc(filename, var):
    # reading nc files without a group
    nc_f = filename
    nc_fid = Dataset(nc_f, 'r')
    out = np.array(nc_fid.variables[var])
    nc_fid.close()
    return np.squeeze(out)

def doy_to_date(year, doy):
    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    return date.month, date.day

def extract_doy_from_filename(filename):
    parts = filename.split('_')
    year_doy_str = parts[-1].split('.')[0]  # Extract '2019121' from the filename
    year = int(year_doy_str[:4])  # Get the year (2019)
    doy = int(year_doy_str[4:])  # Get the DOY (121)
    return year, doy


for month in range(1,13):
    NEI2016_raw = sorted(glob("./NEI_2016_raw/*_2016" + f"{month:02d}" + "*.ncf"))
    dirunal_scales_species_data = {}
    NEI_EMIS = []
    with Dataset(NEI2016_raw[0], 'r') as dataset:
        # Iterate over all variables
        for var_name in dataset.variables:
            var = dataset.variables[var_name]
            # Check if variable has 2 dimensions
            if len(var.dimensions) == 4:
               NEI_EMIS.append(var_name)
    for sp in NEI_EMIS:
        print(f"Processing for {sp}")
        dirunal_scales_species_data[f"{sp}_weekday"] = []
        dirunal_scales_species_data[f"{sp}_weekend"] = []
    for f in NEI2016_raw:
        print("reading NEI for " + str(f))
        # Look for 8-digit number (YYYYMMDD)
        match = re.search(r"\d{8}", f)
        date_str = match.group()
        date_obj = datetime.strptime(date_str, "%Y%m%d").date()
        if date_obj.weekday() >= 5:
           weekend = True
        else:
           weekend = False
        for gas in NEI_EMIS:
            print(gas)
            emis = _read_nc(f,gas)
            #lat = _read_nc(f,'lat')
            #lon = _read_nc(f,'lon')
            # calculate the diurnal cycle
            mean_emis = np.mean(emis[:,:,:], axis=0).squeeze()
            diurnal_scaling_emis = emis[:,:,:]/mean_emis
            mass_factor = np.sum(emis[:,:,:],axis=0)/\
                          np.sum(emis[:,:,:]*diurnal_scaling_emis,axis=0)
            diurnal_scaling_emis = diurnal_scaling_emis*mass_factor

            # remove bad diurnal scales
            diurnal_scaling_emis[np.isnan(diurnal_scaling_emis)]= 1.0
            diurnal_scaling_emis[np.isinf(diurnal_scaling_emis)]= 1.0
            if weekend:
               dirunal_scales_species_data[f"{gas}_weekend"].append(diurnal_scaling_emis)
            else:
               dirunal_scales_species_data[f"{gas}_weekday"].append(diurnal_scaling_emis)
    #dirunal_scales_species_data["lon"] = lon
    #dirunal_scales_species_data["lat"] = lat
    # Convert all float arrays to float32
    dirunal_scales_species_data = {
    key: value.astype(np.float32) if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, np.floating) else value
    for key, value in dirunal_scales_species_data.items()
    }
    for gas in NEI_EMIS:
        dirunal_scales_species_data[f"{gas}_weekend"] = np.nanmean(dirunal_scales_species_data[f"{gas}_weekend"],axis=0)
        dirunal_scales_species_data[f"{gas}_weekday"] = np.nanmean(dirunal_scales_species_data[f"{gas}_weekday"],axis=0)
    savemat("./mat_diurnal_scales/" + f"Scales_2016{month:02d}.mat", dirunal_scales_species_data)
