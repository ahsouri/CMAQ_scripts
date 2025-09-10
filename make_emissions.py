import numpy as np
import datetime
from netCDF4 import Dataset
import warnings
import time
from scipy.io import savemat, loadmat
from scipy.spatial import Delaunay
from scipy.interpolate import NearestNDInterpolator
from joblib import Parallel, delayed
import xarray as xr


warnings.filterwarnings("ignore", category=RuntimeWarning)


def _daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + datetime.timedelta(n)

def _cb06_mapping():
    cb06_species_info = {
    "ACET": {"molecular_mass": 58.08, "num_carbons": 3},
    "ACROLEIN": {"molecular_mass": 56.06, "num_carbons": 3},
    "ALD2": {"molecular_mass": 44.05, "num_carbons": 2},
    "ALD2_PRIMARY": {"molecular_mass": 44.05, "num_carbons": 2},
    "ALDX": {"molecular_mass": 72.11, "num_carbons": 4},
    "BENZ": {"molecular_mass": 78.11, "num_carbons": 6},
    "BUTADIENE13": {"molecular_mass": 54.09, "num_carbons": 4},
    "CH4": {"molecular_mass": 16.04, "num_carbons": 1},
    "CH4_INV": {"molecular_mass": 16.04, "num_carbons": 1},
    "CL2": {"molecular_mass": 70.90, "num_carbons": ""},
    "CO": {"molecular_mass": 28.01, "num_carbons": ""},
    "CO2_INV": {"molecular_mass": 44.01, "num_carbons": 1},
    "ETH": {"molecular_mass": 30.07, "num_carbons": 2},
    "ETHA": {"molecular_mass": 30.07, "num_carbons": 2},
    "ETHY": {"molecular_mass": 44.10, "num_carbons": 2},
    "ETOH": {"molecular_mass": 46.07, "num_carbons": 2},
    "FORM": {"molecular_mass": 30.03, "num_carbons": 1},
    "FORM_PRIMARY": {"molecular_mass": 30.03, "num_carbons": 1},
    "HCL": {"molecular_mass": 36.46, "num_carbons": ""},
    "HONO": {"molecular_mass": 47.01, "num_carbons": ""},
    "IOLE": {"molecular_mass": 68.12, "num_carbons": 5},  # Estimate
    "ISOP": {"molecular_mass": 68.12, "num_carbons": 5},
    "IVOC": {"molecular_mass": 120.0, "num_carbons": 8},  # Estimate
    "KET": {"molecular_mass": 72.11, "num_carbons": 4},
    "MEOH": {"molecular_mass": 32.04, "num_carbons": 1},
    "N2O_INV": {"molecular_mass": 44.01, "num_carbons": ""},
    "NAPH": {"molecular_mass": 128.17, "num_carbons": 10},
    "NH3": {"molecular_mass": 17.03, "num_carbons": ""},
    "NH3_FERT": {"molecular_mass": 17.03, "num_carbons": ""},
    "NO": {"molecular_mass": 30.01, "num_carbons": ""},
    "NO2": {"molecular_mass": 46.01, "num_carbons": ""},
    "NVOL": {"molecular_mass": 150.0, "num_carbons": 10},  # Estimate
    "OLE": {"molecular_mass": 56.11, "num_carbons": 4},
    "PAL": {"molecular_mass": 250.0, "num_carbons": 18},  # Palmitic acid proxy
    "PAR": {"molecular_mass": 14.0, "num_carbons": 1},
    "PCA": {"molecular_mass": 250.0, "num_carbons": ""},
    "PCL": {"molecular_mass": 35.45, "num_carbons": ""},
    "PEC": {"molecular_mass": 12.01, "num_carbons": 1},
    "PFE": {"molecular_mass": 55.85, "num_carbons": ""},
    "PH2O": {"molecular_mass": 18.02, "num_carbons": ""},
    "PK": {"molecular_mass": 39.10, "num_carbons": ""},
    "PMC": {"molecular_mass": 12.01, "num_carbons": 1},
    "PMG": {"molecular_mass": 24.31, "num_carbons": ""},
    "PMN": {"molecular_mass": 14.01, "num_carbons": ""},
    "PMOTHR": {"molecular_mass": 150.0, "num_carbons": ""},
    "PNA": {"molecular_mass": 22.99, "num_carbons": ""},
    "PNCOM": {"molecular_mass": 150.0, "num_carbons": 10},  # Estimate
    "PNH4": {"molecular_mass": 18.04, "num_carbons": ""},
    "PNO3": {"molecular_mass": 62.00, "num_carbons": ""},
    "POC": {"molecular_mass": 12.01, "num_carbons": 1},
    "PRPA": {"molecular_mass": 44.10, "num_carbons": 3},
    "PSI": {"molecular_mass": 28.09, "num_carbons": ""},
    "PSO4": {"molecular_mass": 96.06, "num_carbons": ""},
    "PTI": {"molecular_mass": 47.87, "num_carbons": ""},
    "SO2": {"molecular_mass": 64.07, "num_carbons": ""},
    "SOAALK": {"molecular_mass": 150.0, "num_carbons": 10},  # Estimate
    "SULF": {"molecular_mass": 96.06, "num_carbons": ""},
    "TERP": {"molecular_mass": 136.24, "num_carbons": 10},
    "TOG_INV": {"molecular_mass": 120.0, "num_carbons": 8},  # Estimate
    "TOL": {"molecular_mass": 92.14, "num_carbons": 7},
    "UNK": {"molecular_mass": 100.0, "num_carbons": ""},
    "UNR": {"molecular_mass": 100.0, "num_carbons": ""},
    "VOC_INV": {"molecular_mass": 120.0, "num_carbons": 8},  # Estimate
    "XYLMN": {"molecular_mass": 106.17, "num_carbons": 8}
    }
    return cb06_species_info

def create_emissions_template(source_file, output_file, grid_dims, data_dict, date_value1,date_value2):
    """
    Create template for SMOKE emissions with new grid dimensions and populate with actual data
    
    Parameters:
    - source_file: path to source NetCDF emissions file
    - output_file: path to output template file
    - grid_dims: dict with keys 'nrows', 'ncols', and optionally 'tsteps', 'nlays'
    - data_dict: dict with species names as keys and numpy arrays as values
                 e.g., {'CO': co_data_array, 'NO': no_data_array, ...}
    """
    
    # Open source file
    ds_source = xr.open_dataset(source_file)
    
    # Create new dataset with same global attributes
    ds_new = xr.Dataset(attrs=ds_source.attrs.copy())
    
    # Set up new dimensions
    new_nrows = grid_dims['nrows']
    new_ncols = grid_dims['ncols']
    new_tsteps = grid_dims.get('tsteps', 1)
    new_nlays = grid_dims.get('nlays', 1)
    new_nvars = len(data_dict)
    
    # Create dimensions
    ds_new = ds_new.expand_dims({
        'TSTEP': new_tsteps,
        'LAY': new_nlays,
        'ROW': new_nrows,
        'COL': new_ncols,
        'VAR': new_nvars,
        'DATE-TIME': 2
    })
    
    # Handle TFLAG variable
    if 'TFLAG' in ds_source:
        # Create TFLAG for new number of variables and time steps
        tflag_data = ds_source['TFLAG'][:]
        tflag_data[0:24,:,0] = date_value1
        tflag_data[24,:,0] = date_value2
        # You can populate TFLAG with actual time stamps here if needed
        # For example: tflag_data[:, :, 0] = your_julian_dates
        #             tflag_data[:, :, 1] = your_time_stamps
        
        ds_new['TFLAG'] = xr.DataArray(
            tflag_data,
            dims=['TSTEP', 'VAR', 'DATE-TIME'],
            attrs=ds_source['TFLAG'].attrs.copy()
        )
    
    # Add emission species variables with actual data
    for species_name, species_data in data_dict.items():
        # Get attributes from source file if the species exists there
        attrs = {}
        if species_name in ds_source:
            attrs = ds_source[species_name].attrs.copy()
        
        # Ensure data has the right shape: (TSTEP, LAY, ROW, COL)
        if species_data.ndim == 2:
            # If 2D data, add time and layer dimensions
            final_data = species_data[np.newaxis, np.newaxis, :, :]
        elif species_data.ndim == 3:
            # If 3D data, add one dimension (either time or layer)
            if species_data.shape[0] == new_tsteps:
                # Assume first dimension is time, add layer dimension
                final_data = species_data[:, np.newaxis, :, :]
            else:
                # Add time dimension
                final_data = species_data[np.newaxis, :, :, :]
        elif species_data.ndim == 4:
            # Data already has correct dimensions
            species_data[np.isnan(species_data)]=0.0
            final_data = species_data
        else:
            raise ValueError(f"Invalid data shape for {species_name}: {species_data.shape}")
        
        # Verify final shape matches expected dimensions
        expected_shape = (new_tsteps, new_nlays, new_nrows, new_ncols)
        if final_data.shape != expected_shape:
            raise ValueError(f"Data shape {final_data.shape} doesn't match expected {expected_shape} for {species_name}")
        
        ds_new[species_name] = xr.DataArray(
            final_data.astype(np.float32),
            dims=['TSTEP', 'LAY', 'ROW', 'COL'],
            attrs=attrs
        )
    
    # Update global attributes for new grid
    ds_new.attrs.update({
        'NCOLS': new_ncols,
        'NROWS': new_nrows,
        'NLAYS': new_nlays,
        'NVARS': new_nvars,
        'SDATE': int(date_value1),
        'P_ALP': 38.9799995422363,
        'P_BET': 38.9799995422363,
        'P_GAM': -100.212997436523,
        'XCENT': -100.212997436523,
        'YCENT': 38.5340003967285,
        'XORIG': -2399755,
        'YORIG': -1640765,
        'XCELL': 8000,
        'YCELL': 8000
    })


    ds_new.attrs.update({
        'NCOLS': new_ncols,
        'NROWS': new_nrows,
        'NLAYS': new_nlays,
        'NVARS': new_nvars
    })    
    # Save to file
    ds_new.to_netcdf(output_file)
    
    # Close datasets
    ds_source.close()
    ds_new.close()
    
    print(f"Emissions file created: {output_file}")
    print(f"Grid: {new_nrows} x {new_ncols}")
    print(f"Species: {list(data_dict.keys())}")

def _read_nc(filename, var):
    # reading nc files without a group
    nc_f = filename
    nc_fid = Dataset(nc_f, 'r')
    out = np.array(nc_fid.variables[var])
    nc_fid.close()
    return np.squeeze(out)

from pyproj import Geod

def grid_area(lat, lon, fill_value=0.0):
    """
    Calculate area (m²) of each cell in an irregular 2D lat/lon grid.
    Output is resized to match original shape with padding.

    Parameters:
        lat (2D array): Latitude in degrees
        lon (2D array): Longitude in degrees
        fill_value (float): Value to fill in padded edge cells (default: 0.0)
    
    Returns:
        area (2D array): Area of grid cells in m², shape (nlat, nlon)
    """
    geod = Geod(ellps='WGS84')
    nlat, nlon = lat.shape
    area_core = np.empty((nlat - 1, nlon - 1), dtype=np.float64)

    # Precompute points for all cells
    for i in range(nlat - 1):
        lat_t = lat[i:i+2, :]  # (2, nlon)
        lon_t = lon[i:i+2, :]  # (2, nlon)
        for j in range(nlon - 1):
            # Corner points of cell (i,j) in counter-clockwise order
            lons = [lon_t[0, j], lon_t[0, j+1], lon_t[1, j+1], lon_t[1, j]]
            lats = [lat_t[0, j], lat_t[0, j+1], lat_t[1, j+1], lat_t[1, j]]
            area_core[i, j], _ = geod.polygon_area_perimeter(lons, lats)

    # Pad to match original shape
    area_full = np.full((nlat, nlon), fill_value)
    area_full[:-1, :-1] = np.abs(area_core)

    return area_full

def NEI_extracter(tri1, tri2, emis, specie_info, date_i, lon_org, lat_org):
    # we need to merge NEI-2016 with QFED
    # NEI2016
    NEI_file = "/discover/nobackup/asouri/SHARED/NEI_2016/nei2016_monthly/2016fh_16j_merge_0pt1degree_month_" +\
        f"{date_i.month:02d}" + ".ncf"
    print("Reading NEI file from " + NEI_file)
    try:
       NEI_emis = _read_nc(NEI_file, emis)
       dataset = Dataset(NEI_file, 'r')
       unit_NEI_GC = getattr(dataset.variables[emis],'units')
       dataset.close()
    except:
       print(f"GC files does not have any values for {emis}")
       return np.zeros((np.shape(lon_org)[0],np.shape(lon_org)[1],25))
    lon_NEI = _read_nc(NEI_file, "lon")
    lat_NEI = _read_nc(NEI_file, "lat")
    lon_NEI, lat_NEI = np.meshgrid(lon_NEI, lat_NEI)
    #points = np.zeros((np.size(lon_NEI), 2))
    #points[:, 0] = lon_NEI.flatten()
    #points[:, 1] = lat_NEI.flatten()
    #tri = Delaunay(points)
    interpolator = NearestNDInterpolator(tri1, (NEI_emis[:, :]).flatten())
    NEI_emis_mapped = interpolator((lon_org, lat_org))
    # remove data outside of lon_NEI and lat_NEI max and mins
    inside_box = (
        (lat_org >= np.min(lat_NEI.flatten())) &
        (lat_org <= np.max(lat_NEI.flatten())) &
        (lon_org >= np.min(lon_NEI.flatten())) &
        (lon_org <= np.max(lon_NEI.flatten()))
    )
    # apply mask to NOx data: keep values inside the box, zero out others
    NEI_emis_mapped = np.where(inside_box, NEI_emis_mapped, 0.0)
    # apply the diurnal factor
    diurnal_scale_file = "/discover/nobackup/asouri/SHARED/NEI_2016/diurnal_scales/Scales_2016" +\
        f"{date_i.month:02d}.mat"
    print("Reading the scaling factor file from " + diurnal_scale_file)
    diurnal_scales = loadmat(diurnal_scale_file)
    # Check if it's a weekend
    if date_i.weekday() >= 5:  # 5 = Saturday, 6 = Sunday
        diurnal_scales = diurnal_scales[f"{emis}_weekend"]
    else:
        diurnal_scales = diurnal_scales[f"{emis}_weekday"]
    lat_scale = _read_nc(
        "/discover/nobackup/asouri/SHARED/NEI_2016/diurnal_scales/GRIDCRO2D_20190201.nc4", 'LAT')
    lon_scale = _read_nc(
        "/discover/nobackup/asouri/SHARED/NEI_2016/diurnal_scales/GRIDCRO2D_20190201.nc4", 'LON')
    #points = np.zeros((np.size(lat_scale), 2))
    #points[:, 0] = lon_scale.flatten()
    #points[:, 1] = lat_scale.flatten()
    #tri = Delaunay(points)
    emis_to_saved = np.zeros((np.shape(NEI_emis_mapped)[0],np.shape(NEI_emis_mapped)[1],25))
    for hour in range(0, 25):
        hour_t = hour # we need to have a 25th index in CMAQ
        if hour == 24:
           hour_t = 23
        interpolator = NearestNDInterpolator(
            tri2, (diurnal_scales[hour_t, :, :]).flatten())
        diurnal_scales_mapped = interpolator((lon_org, lat_org))
        # make the dirunal scales = 1.0 outside of the domain
        inside_box = (
            (lat_org >= np.min(lat_scale.flatten())) &
            (lat_org <= np.max(lat_scale.flatten())) &
            (lon_org >= np.min(lon_scale.flatten())) &
            (lon_org <= np.max(lon_scale.flatten()))
        )
        # apply the mask
        diurnal_scales_mapped = np.where(
            inside_box, diurnal_scales_mapped, 1.0)
        emis_to_saved[:,:,hour] = diurnal_scales_mapped * NEI_emis_mapped
        #area = grid_area(lon_org,lat_org)
        area = 8000.0**2
        emis_to_saved[:,:,hour] = emis_to_saved[:,:,hour]*area*1000.0 # this is now g/s

    print(unit_NEI_GC)
    if unit_NEI_GC == 'kg/m2/s':
            emis_to_saved = emis_to_saved/specie_info["molecular_mass"] # now this is moles/s
    elif unit_NEI_GC == 'kgNO2/m2/s':
            emis_to_saved = emis_to_saved*(30.0/46.0)/specie_info["molecular_mass"]
            print('hello')
    elif unit_NEI_GC == 'kgC/m2/s':
            emis_to_saved = emis_to_saved*specie_info["molecular_mass"]/(12.01*specie_info["num_carbons"])/specie_info["molecular_mass"]
    else:
        raise ValueError("The unit is undefined, please check the data")      
    return emis_to_saved

def QFED_extracter(emis,  date_i, lon_org, lat_org):
    # QFED emissions should be vertically allocated (65% within PBL and 35% in the free troposphere)
    pass

def process_emis(skeleton,tri1,tri2,CB06_map,date_i,lon_org,lat_org,output_file):

     with Dataset(skeleton, 'r') as dataset:
         # Iterate over all variables
         for var_name in dataset.variables:
             var = dataset.variables[var_name]
             # Check if variable has 2 dimensions
             if len(var.dimensions) == 4:
                    print(f"Processing 4D variable: {var_name}")
                    # read unit
                    unit = getattr(dataset.variables[var_name], 'units')
                    result =  NEI_extracter(tri1, tri2, var_name, CB06_map[var_name], date_i, lon_org, lat_org)
                    result = np.moveaxis(result, -1, 0)  # Now shape is (25, 487, 757)
                    data_output[var_name] = np.expand_dims(result, axis=1)
     # scaling due to the year of our target (2023) while using 2016 emissions
     data_output['NO'] = 1.15*data_output['NO']
     data_output['NO2'] = 1.15*data_output['NO2']
     create_emissions_template('./emis_mole_all_20171226_12US1_nobeis_norwc_WR413_MYR_2017.nc4',
                                  output_file,grid_info,data_output, date_i.strftime("%Y%j"),
                                  (date_i + datetime.timedelta(days=1)).strftime('%Y%j'))
     return None

if __name__ == "__main__":

    skeleton = "./emis_mole_all_20171226_12US1_nobeis_norwc_WR413_MYR_2017.nc4"
    lat_org = _read_nc('./GRIDCRO2D.nc','LAT')
    lon_org = _read_nc('./GRIDCRO2D.nc','LON')
    CB06_map = _cb06_mapping()
    data_output = {}
    grid_info = {
       'nrows': 440,
       'ncols': 710,
       'tsteps': 25,
       'nlays': 1
    }

    # to increase the speed, let's make tri
    NEI_file = "/discover/nobackup/asouri/SHARED/NEI_2016/nei2016_monthly/2016fh_16j_merge_0pt1degree_month_" +\
        "01" + ".ncf"
    lon_NEI = _read_nc(NEI_file, "lon")
    lat_NEI = _read_nc(NEI_file, "lat")
    lon_NEI, lat_NEI = np.meshgrid(lon_NEI, lat_NEI)
    points = np.zeros((np.size(lon_NEI), 2))
    points[:, 0] = lon_NEI.flatten()
    points[:, 1] = lat_NEI.flatten()
    tri1 = Delaunay(points)

    lat_scale = _read_nc(
        "/discover/nobackup/asouri/SHARED/NEI_2016/diurnal_scales/GRIDCRO2D_20190201.nc4", 'LAT')
    lon_scale = _read_nc(
        "/discover/nobackup/asouri/SHARED/NEI_2016/diurnal_scales/GRIDCRO2D_20190201.nc4", 'LON')
    points = np.zeros((np.size(lat_scale), 2))
    points[:, 0] = lon_scale.flatten()
    points[:, 1] = lat_scale.flatten()
    tri2 = Delaunay(points)
    # loop over whole days ranging from 2023 till the end of 2024
    datarange = _daterange(datetime.date(2023, 7, 29), datetime.date(2024, 8, 1))
    datarange = list(datarange)
    output_files = []
    print(len(datarange))
    for date_i in datarange:
       output_files.append(f"EMIS_NEI_2016_ScaledEPA_All_Anthro_OneLayer_{date_i.strftime('%Y%m%d')}")
    out = Parallel(n_jobs=10,verbose=10)(delayed(process_emis)(
           skeleton, tri1, tri2,CB06_map, datarange[k],lon_org,lat_org,output_files[k]) for k in range(len(datarange)))

