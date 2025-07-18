import pyproj
import numpy as np
from netCDF4 import Dataset
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.interpolate import RBFInterpolator
from scipy import signal
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from scipy.spatial import cKDTree
import xarray as xr

def create_emissions_template(source_file, output_file, grid_dims, data_dict):
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

def _interpolosis(interpol_func, Z: np.array, X: np.array, Y: np.array, interpolator_type: int, dists: np.array, threshold: float) -> np.array:
    # to make the interpolator() shorter
    if interpolator_type == 1:
        interpolator = LinearNDInterpolator(
            interpol_func, (Z).flatten(), fill_value=np.nan)
        ZZ = interpolator((X, Y))
        ZZ[dists > threshold*3.0] = np.nan
    elif interpolator_type == 2:
        interpolator = NearestNDInterpolator(interpol_func, (Z).flatten())
        ZZ = interpolator((X, Y))
        ZZ[dists > threshold*3.0] = 0.0
    elif interpolator_type == 3:
        interpolator = RBFInterpolator(
            interpol_func, (Z).flatten(), neighbors=5)
        XX = np.stack([X.ravel(), Y.ravel()], -1)
        ZZ = interpolator(XX)
        ZZ = ZZ.reshape(np.shape(X))
        ZZ[dists > threshold*3.0] = np.nan
    else:
        raise Exception(
            "other type of interpolation methods has not been implemented yet")
    return ZZ

def _boxfilter(size_kernel_x, size_kernel_y) -> np.array:
    return np.ones((int(size_kernel_x), int(size_kernel_y)))/(size_kernel_x*size_kernel_y)

def get_latlon_old():
    # Define the Lambert Conformal Conic projection for your input data
    proj_lcc = pyproj.Proj(proj='lcc', lat_1=33, lat_2=45, lat_0=40, lon_0=-97, 
                       x_0=0, y_0=0, datum='NAD83', units='m')

    # Calculate grid coordinates
    xorig, yorig = -2952000, -2772000
    xcell, ycell = 1000, 1000
    ncols, nrows = 6192, 5328

    # Create coordinate arrays
    x = np.arange(ncols) * xcell + xorig + xcell/2  # Cell centers
    y = np.arange(nrows) * ycell + yorig + ycell/2

    X, Y = np.meshgrid(x, y)

    # Convert to lat/lon
    lon, lat = proj_lcc(X, Y, inverse=True)
    print(np.shape(lat))
    print(np.shape(lon))
    return lon,lat

def get_latlon():
    # Define the Lambert Conformal Conic projection for your input data
    # Parameters from your provided metadata
    proj_lcc = pyproj.Proj(proj='lcc', 
                          lat_1=33.0,     # P_ALP
                          lat_2=45.0,     # P_BET
                          lat_0=40.0,     # YCENT
                          lon_0=-97.0,    # P_GAM/XCENT
                          x_0=0, 
                          y_0=0, 
                          datum='NAD83', 
                          units='m')

    # Grid parameters from your metadata
    xorig, yorig = -2556000.0, -1728000.0  # XORIG, YORIG
    xcell, ycell = 12000.0, 12000.0        # XCELL, YCELL
    ncols, nrows = 459, 299                # NCOLS, NROWS

    # Create coordinate arrays (cell centers)
    x = np.arange(ncols) * xcell + xorig + xcell/2
    y = np.arange(nrows) * ycell + yorig + ycell/2

    # Create 2D grid of coordinates
    X, Y = np.meshgrid(x, y)

    # Convert projection coordinates to lat/lon
    lon, lat = proj_lcc(X, Y, inverse=True)
    
    print(f"Latitude shape: {np.shape(lat)}")
    print(f"Longitude shape: {np.shape(lon)}")
    
    return lon, lat

lon_input,lat_input = get_latlon()
lon_output = _read_nc('/discover/nobackup/asouri/MODELS/CMAQv5.5/data/mcip/CONUS_8km/GRIDCRO2D_CONUS_8km_20230618.nc','LON')
lat_output = _read_nc('/discover/nobackup/asouri/MODELS/CMAQv5.5/data/mcip/CONUS_8km/GRIDCRO2D_CONUS_8km_20230618.nc','LAT')
points = np.zeros((np.size(lon_input), 2))
points[:, 0] = lon_input.flatten()
points[:, 1] = lat_input.flatten()
tri = Delaunay(points)
# calculate distance to remove too-far estimates
tree = cKDTree(points)
grid = np.zeros((2, np.shape(lon_output)[0], np.shape(lon_output)[1]))
grid[0, :, :] = lon_output
grid[1, :, :] = lat_output
xi = _ndim_coords_from_arrays(tuple(grid), ndim=points.shape[1])
dists, _ = tree.query(xi)
data_output_dict = {}
with Dataset('./beis_norm_emis_12US1_2020ha2_cb6_20k.ncf', 'r') as dataset:
    # Iterate over all variables
    for var_name in dataset.variables:
        var = dataset.variables[var_name]
        # Check if variable has 2 dimensions
        if len(var.dimensions) == 4:
            print(f"Processing 2D variable: {var_name}")
            # Read data
            data = var[:,:,:,:].squeeze()
            # first convolve the 1km data with a 8x8 km low pass filter
            #kernel = _boxfilter(8,8)
            #data = signal.convolve2d(data, kernel, boundary='symm', mode='same')
            # now interpolate with NN
            data_output = _interpolosis(tri, data, lon_output, lat_output, 2, dists, 0.08)
            data_output_dict.update({var_name: data_output})

grid_info = {
    'nrows': 487,
    'ncols': 757,
    'tsteps': 1,
    'nlays': 1
}
create_emissions_template('./beis_norm_emis_12US1_2020ha2_cb6_20k.ncf', './BELD5_avg_output.ncf',
                         grid_info, data_output_dict)
