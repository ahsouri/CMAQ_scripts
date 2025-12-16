import numpy as np
from netCDF4 import Dataset
import xarray as xr
import glob
import datetime

def _read_nc(filename, var):
    # reading nc files without a group
    nc_f = filename
    nc_fid = Dataset(nc_f, 'r')
    out = np.array(nc_fid.variables[var])
    nc_fid.close()
    return np.squeeze(out)

def _daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + datetime.timedelta(n)

def _calculate_ctm_partial_column(deltap, profile):
    """Calculate CTM partial column."""
    Mair = 28.97e-3
    g = 9.80665
    N_A = 6.02214076e23
    return deltap[0:24,...] * profile[0:24,...]*1e3 / g / Mair * N_A * 1e-4 * 1e-15 * 100.0 * 1e-9

def CMAQ_PA_reader(fname_cro3d,fname_cro2d,fname_pa,date_str):

    print("Currently reading: " + fname_pa.split('/')[-1])
    # reading time and coordinates
    prs = _read_nc(fname_cro3d, 'PRES').astype('float32')/100.0  # hPa
    surf_prs = _read_nc(fname_cro2d, 'PRSFC').astype('float32')/100.0
    delp = prs.copy()
    # calculate delta pressure
    for i in range(0, np.shape(prs)[1]):
        if i == 0:  # the first layer
            delp[:, i, :, :] = surf_prs - 0.5*(prs[:, 0, :, :] + prs[:, 1, :, :])
        elif i == np.shape(delp)[1]-1:  # the last layer
            delp[:, i, :, :] = prs[:, i-1, :, :] - prs[:, i, :, :]
        else:  # the between
            delp[:, i, :, :] = (prs[:, i, :, :] + prs[:, i-1, :, :]) * \
                    0.5 - (prs[:, i+1, :, :] + prs[:, i, :, :])*0.5
    var_col={}
    with Dataset(fname_pa, 'r') as dataset:
         for var_name in dataset.variables:
             var = dataset.variables[var_name]
             # Check if variable has 2 dimensions
             if len(var.dimensions) == 4:
                print(f"Processing 4D variable: {var_name}")
                var_col[var_name] = np.sum(_calculate_ctm_partial_column(delp,np.array(var[:])),axis=1).squeeze()
    
    # Open source file
    ds_source = xr.open_dataset(fname_pa)
    
    # Create new dataset with same global attributes
    ds_new = xr.Dataset(attrs=ds_source.attrs.copy())
    
    grid_info = {
       'nrows': 440,
       'ncols': 710,
       'tsteps': 24,
       'nlays': 1
    }
    # Set up new dimensions
    new_nrows = grid_info['nrows']
    new_ncols = grid_info['ncols']
    new_tsteps = grid_info.get('tsteps', 24)
    new_nlays = grid_info.get('nlays', 1)
    new_nvars = len(var_col)
    
    # Create dimensions
    ds_new = ds_new.expand_dims({
        'TSTEP': new_tsteps,
        'LAY': new_nlays,
        'ROW': new_nrows,
        'COL': new_ncols,
        'VAR': new_nvars,
        'DATE-TIME': 2
    })

    ds_new['TFLAG'] = ds_source['TFLAG']
    
    # Add emission species variables with actual data
    for species_name, species_data in var_col.items():
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
    

    ds_new.attrs.update({
        'NCOLS': new_ncols,
        'NROWS': new_nrows,
        'NLAYS': new_nlays,
        'NVARS': new_nvars
    })    
    # Save to file
    ds_new.to_netcdf("./COLUMN_PA_" + date_str + ".nc")
    # Close datasets
    ds_source.close()
    ds_new.close()

if __name__ == "__main__":

    datarange = _daterange(datetime.date(2024, 5, 15), datetime.date(2024, 10, 1))
    datarange = list(datarange)
    mcip_dir = "/discover/nobackup/asouri/MODELS/CMAQv5.5/data/mcip/CONUS_8km"
    cctm_dir = "/discover/nobackup/asouri/MODELS/CMAQv5.5/data/output_CCTM_v55_intel_CONUS_8km"
    for date in datarange:
        met_cro2d = f"{mcip_dir}/METCRO2D_CONUS_8km_{date.strftime('%Y%m%d')}.nc"
        met_cro3d = f"{mcip_dir}/METCRO3D_CONUS_8km_{date.strftime('%Y%m%d')}.nc"
        pa_file = f"{cctm_dir}/CCTM_PA_1_v55_intel_CONUS_8km_{date.strftime('%Y%m%d')}.nc"
        CMAQ_PA_reader(met_cro3d,met_cro2d,pa_file,date.strftime('%Y%m%d'))
