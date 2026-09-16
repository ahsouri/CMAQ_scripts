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

def CMAQ_nox_reader(fname_pa,fname_conc,date_str):

    print("Currently reading: " + fname_pa.split('/')[-1])
    # reading time and coordinates
    var_pa_col={}
    var_col={}
    with Dataset(fname_pa, 'r') as dataset:
         for var_name in dataset.variables:
             var = dataset.variables[var_name]
             # Check if variable has 3 dimensions
             if len(var.dimensions) == 4:
                if 'VDIF' in var_name:
                   print(f"Processing 4D variable: {var_name}")
                   var_pa_col[var_name] = np.array(var[:])
    with Dataset(fname_conc, 'r') as dataset:
         for var_name in dataset.variables:
             var = dataset.variables[var_name]
             # Check if variable has 3 dimensions
             if len(var.dimensions) == 4:
                   print(f"Processing 4D variable: {var_name}")
                   var_col[var_name] = np.array(var[:])

    print(var_col.keys())
    # make emis_no2
    NOx = (var_col['NO2']+var_col['NO']+var_col['NO3']+2*var_col['N2O5']+var_col['HONO'])
    emis_tend = var_pa_col['VDIF_NO2']+var_pa_col['VDIF_NO']+var_pa_col['VDIF_NO3']+2*var_pa_col['VDIF_N2O5']+var_pa_col['VDIF_HONO']
    emis_no2 = var_col['NO2'][0:24,...]/NOx[0:24,...]*emis_tend
    var_col ={'EMIS_NO2':emis_no2}
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
    new_tsteps = grid_info.get('tsteps', 25)
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
    ds_new.to_netcdf("./COLUMN_EMIS_NO2_" + date_str + ".nc")
    # Close datasets
    ds_source.close()
    ds_new.close()

if __name__ == "__main__":

    datarange = _daterange(datetime.date(2023, 10, 5), datetime.date(2024, 10, 6))
    datarange = list(datarange)
    data_dir = "/discover/nobackup/asouri/GITS/CMAQ_scripts/"
    for date in datarange:
        pa_file = f"{data_dir}/COLUMN_PA_{date.strftime('%Y%m%d')}.nc"
        conc_file = f"{data_dir}/COLUMN_CONC_{date.strftime('%Y%m%d')}.nc"
        print(pa_file)
        print(conc_file)
        CMAQ_nox_reader(pa_file,conc_file,date.strftime('%Y%m%d'))
