# 20200430 hwang adapted gga's TEMPO routines to MEASURES
# mapping based on cartopy
import sys
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.style as mplstyle
# comment the following out for publication quality 
# leave it in for faster rendering
mplstyle.use('fast')
#
####################################################################
# helper routine
####################################################################
#
#------
# Preparing for global plot
#------
def pre_mapping_plot(proj,central_longitude=0.,spnum=111,
                    extent=(-180.,180.,-90.,90.)):
    
    if (proj == 'Mollweide'):
       projection=ccrs.Mollweide(central_longitude=central_longitude)
    elif (proj == 'Robinson'):
       projection=ccrs.Robinson(central_longitude=central_longitude)
    elif (proj == 'Miller'):
       projection=ccrs.Miller(central_longitude=central_longitude)
    elif (proj == 'Mercator'):
       projection=ccrs.Mercator(central_longitude=central_longitude)
    else:
       projection=ccrs.PlateCarree(central_longitude=central_longitude)

    ax = plt.subplot(spnum,projection=projection)
    # set global extent as necessary
    if (extent == (-180.,180.,-90.,90.)): 
       ax.set_global()
    else:
       ax.set_extent(extent, crs=ccrs.PlateCarree())

    ax.stock_img()

    return ax, projection

#------
# add map details, e.g., gridline and labels
#------
def finish_map_detail(ax,extent,central_longitude=0.,
           title=' '):
    import cartopy.mpl.ticker as cticker

    print('   add ax title '+title)
    ax.set_title(title)

    print('   add coastline and borders')
    ax.add_feature(cfeature.COASTLINE) 
    ax.add_feature(cfeature.BORDERS) 
    ax.add_feature(cfeature.LAKES)

    # somehow cartopy gridliner only worked for some cases
    # the following is a work-around from the web
    # it only works if the projected map is rectangular
    # such as, PlateCarree, Mercator, Miller
    comm = '''
    print,('   add gridline')
 
    if (extent == (-180.,180.,-90.,90.)) :
       latlab = [-45., 0., 45.]
       if (central_longitude == 0.) :
           lonlab = [-180.,-90.,0.,90.,180.]
       else:
           lonlab = [0., 90., 180., 270., 360.]

       ax.set_xticks(lonlab,crs=ccrs.PlateCarree())
       ax.set_xticklabels(str(lonlab))
       ax.set_yticks(latlab,crs=ccrs.PlateCarree())
       ax.set_yticklabels(str(latlab)) 
       lon_formatter = cticker.LongitudeFormatter()
       lat_formatter = cticker.LatitudeFormatter()
       ax.xaxis.set_major_formatter(lon_formatter)
       ax.yaxis.set_major_formatter(lat_formatter)
       ax.grid(alpha=0.5,linestyle='--',color='white')
    '''
#------
# post plot action 
#------
def post_mapping_plot(fig,figtitle=' ',filename=None,block=False):
    fig.suptitle(figtitle)
    
    if filename:
       plt.savefig(filename)
       plt.close()
    else:
       plt.show(block=block)

#------
# transform and mask data
#------
def get_masked_data(data,lon,lat,central_longitude):
    import numpy.ma as ma

    pixsize = 1.0
    if (central_longitude == 0.):
       leftlimit = -180.+pixsize
       rightlimit = 180.-pixsize
       ind = (lon>180.)
       lon[ind] = lon[ind]-360.
    else:
       leftlimit = pixsize
       rightlimit = 360. - pixsize
       ind = (lon < 0.)
       lon[ind] = lon[ind]+360.

    mapmask = (data<0.) | (abs(lat)>85.) | \
              (lon<leftlimit) | (lon > rightlimit)

    z = ma.array(data,mask=mapmask)

    return z,lon

#------       
# transform and mask data, lon, lat
# for use with map_data2d_color 
#------
def transform_mask_datalonlat(data,lon,lat, projection, central_longitude):
    import numpy.ma as ma 
    pixsize =1.0
    transxy=projection.transform_points(ccrs.Geodetic(),lon,lat)

    xx = transxy[:,:,0]
    yy = transxy[:,:,1]

    if (central_longitude == 0):
        leftlimit = -180. + pixsize
        rightlimit = 180. - pixsize
        ind = (lon > 180.)
        lon[ind] = lon[ind] - 360.
    else:
        leftlimit = pixsize 
        rightlimit = 360. - pixsize
        ind = (lon < 0.)
        lon[ind] = lon[ind] + 360.

    mapmask1 = (data < 0.) | (abs(lat)>85.0) | \
               (lon < leftlimit) | (lon > rightlimit)

    zz = ma.array(data, mask=mapmask1)

    return  zz,xx,yy
#
####################################################
# General plotting 
####################################################
#
#------
# add 2D data on existing map
#------
def add_data2d_cartopy(ax,data,lon,lat,vmin=None,vmax=None,
        cmap='jet',alpha=0.75,Nlevels=21,method='pcolormesh'):
         
    orgproj = ccrs.PlateCarree()

    #contourf appears slower than pcolormesh
    #Note transform should be set to the "from" projection

    if (method == 'pcolormesh'):
        img = ax.pcolormesh(lon, lat, data, cmap=cmap,
             vmin=vmin,vmax=vmax,transform=orgproj,alpha=alpha)
    else:
        levels=np.linspace(vmin,vmax,Nlevels) 
        img = ax.contourf(lon,lat,data,cmap=cmap,levels=levels,
              vmin=vmin,vmax=vmax,transform=orgproj,alpha=alpha)

    return img

#------
# Plot 2D data on world map using various projections
#------
def map_data2d_cartopy(data,lon,lat,vmin=None,vmax=None,
        spnum=111,method='pcolormesh',proj='PlateCarree',
        Nlevels=21, cbar_legend=None,cmap='jet',
        central_longitude=0.,extent=(-180.,180.,-90.,90.),
        title=' ',filename=None):

    if (lat.shape != data.shape) or (lon.shape != data.shape) :
        print(lat.shape,lon.shape,data.shape)
        sys.exit('shape confliction')

    print('initiating plot') 
    fig = plt.figure()

    ax, projection=pre_mapping_plot(proj,extent=extent,spnum=spnum,
                  central_longitude=central_longitude)

    print('masking out invalid pixels')
    z, lon = get_masked_data(data,lon,lat,central_longitude)

    print('ploting data on map')
    img = add_data2d_cartopy(ax, z, lon, lat, cmap=cmap,
            vmin=vmin,vmax=vmax,alpha=0.75,method=method)

    cbar=plt.colorbar(img,ax=ax,orientation='horizontal',
                      fraction=0.03)
    cbar.set_label(cbar_legend)

    print('perform finishing touches')
    axtitle = title 
    finish_map_detail(ax,extent,title=axtitle)
  
    figtitle = proj
    post_mapping_plot(fig,figtitle=figtitle,filename=filename,block=False)

#------
# Plot 2D data/lon/lat on global map using various projections
# although works, it is not as elegant as map_data2d_cartopy
#------
def map_data2d_color(data,lon,lat,cmap='jet',spnum=111,
        vmin=None,vmax=None,cbar_legend=None, proj='Mollweide', 
        central_longitude=0.,extent=(-180.,180.,-90.,90.),
        title=' ',filename=None):


    if (lat.shape != data.shape) or (lon.shape != data.shape) :
        print(lat.shape,lon.shape,data.shape)
        sys.exit('shape confliction')

    print('initiating plot') 
    fig = plt.figure()
    ax, projection = pre_mapping_plot(proj,extent=extent,
           spnum=spnum,central_longitude=central_longitude)

    print('transforming and masking data,lon,lat' )
    zz,xx,yy = transform_mask_datalonlat(data,lon,lat,
                         projection,central_longitude)

    img = ax.pcolormesh(xx,yy,zz,
             vmin=vmin,vmax=vmax,alpha=0.75,cmap=cmap)
            
    axtitle = title
    finish_map_detail(ax,extent,title=axtitle)

    cbar=plt.colorbar(img,ax=ax,orientation='horizontal',
              fraction=0.03)
    cbar.set_label(cbar_legend)

    figtitle=proj
    post_mapping_plot(fig,filename=filename,figtitle=figtitle,block=False)

#
##########################################################
#   Unit Conversion
##########################################################
#
#------
# convert from molec/cm2 to mm for H2O
#------
def h2o_molecpercm2_mm(data, missing=-1.e30):
    if (missing > 0.) :
        missing = -1.e30

    factor = 2.98904e-22
    output = data * factor
    ibad = (data < 0.)
    output[ibad] = missing

    return output

#------
# convert from mm to molec/cm2 for H2O
#------
def h2o_mm_molecpercm2(data, missing=-1.e30):
    if (missing > 0.) :
        missing = -1.e30
    factor = 2.98904e-22
    output = data / factor
    ibad = data < 0.
    output[ibad] = missing

    return output
#
#########################################################
# Measures reading
#########################################################
#
#------
# read one L2 measures output file
#------
def read_measures_l2(fnm,varnm,missing=-1.e30):
    from netCDF4 import Dataset

    f = Dataset(fnm)

    mdqfl = f['key_science_data']['main_data_quality_flag'][:,:]
    latitude = f['geolocation']['latitude'][:,:]
    longitude = f['geolocation']['longitude'][:,:]

    if (varnm == 'column_amount'):
       var = f['key_science_data']['column_amount'][:,:]
    elif (varnm == 'column_uncertainty'):
       var = f['key_science_data']['column_uncertainty'][:,:]
    elif (varnm == 'fit_rms_residual'):
       var = f['qa_statistics']['fit_rms_residual'][:,:] 
    elif (varnm == 'amf'):
       var = f['support_data']['amf'][:,:]
       #flag can be used for filtering
       #flag = f['support_data']['amf_diagnostic_flag'][:,:]
    elif (varnm == 'cloud_fraction'):
       var = f['support_data']['cloud_fraction'][:,:]
    elif (varnm == 'cloud pressure'):
       var = f['support_data']['cloud_pressure'][:,:]
    elif (varnm == 'albedo'):
       var = f['support_data']['albedo'][:,:]
    elif (varnm == 'surface_pressure'):
       var = f['support_data']['surface_pressure'][:,:]
    else:
       raise('trouble reading '+varnm+' from '+fnm) 

    var[mdqfl != 0] = missing

    return var,longitude,latitude
#------
# figure out column_amount normalization factor
#------
def figure_out_normcol(var):
    import math

    capct = np.percentile(var, 95)
    thispower = math.floor(math.log(capct,10))
    normcol = 10.**thispower

    return normcol,capct

#
##########################################################
# Measures plotting
##########################################################
#
#------
# Map column_amount over the globe for one L2 file
#------
def measures_plot_vcd_1orbit(fnm,missing=-1.e30,vmin=0.,vmax=4.,
        proj='PlateCarree',central_longitude=0.,filename=None,
        spnum=111, method='pcolormesh',unit_conversion=True,
        extent=(-180.,180.,-90.,90.),usecart=True):

    print('read '+fnm) 
    c_a,lon,lat = read_measures_l2(fnm,'column_amount',
                                   missing=missing)

    normcol,capct = figure_out_normcol(c_a)
    data = np.divide( c_a , normcol)
    v_max = vmax 
    legend = 'x'+str(normcol)+' molec/cm2'

    if ('H2O' in fnm):
       # unit_conversion is for H2O only
       if (unit_conversion):
          print('convert H2O from molec/cm2 to mm')
          normcol = 1.
          v_max = 120.
          legend = 'H2O (mm)'  
          data = h2o_molecpercm2_mm(c_a)

    print('normcol=',normcol)
    print('legend=',legend)

    print('plot data on world map')
    abc0=fnm.split('/')
    abc = abc0[-1].split('_')
    title=abc[0] + ' '+abc[1]

    if (usecart):
        print('map_data2d_caropy')
        map_data2d_cartopy(data,lon,lat,vmin=vmin,vmax=v_max,
             central_longitude=central_longitude,title=title,
             filename=filename,cbar_legend=legend,proj=proj,
             method=method,extent=extent,spnum=spnum) 
    else:
        print('map_data2d_color')
        map_data2d_color(data,lon,lat,vmin=vmin,vmax=v_max,
             central_longitude=central_longitude,title=title,
             filename=filename,cbar_legend=legend,
             proj=proj,spnum=spnum)

    return 

#------
# Map column_amount over the globe for files in list 
#------
def measures_plot_vcd_norbits(list,vmin=0.,vmax=4.,filename=None,
        proj='PlateCarree',extent=(-180.,180.,-90.,90.),
        central_longitude=0.,title=' ',cmap='jet'):

    # create figure with specified projection
    fig = plt.figure()
    ax, projection = pre_mapping_plot(proj,extent=extent,
                     central_longitude=central_longitude)

    # get normalization factor and legend
    # if molecule is H2O, do unit conversion
    f0 = open(list,'r')
    line0 = f0.readline()
    f0.close()
    if ('H2O' in line0):
        normcol = 1.
        v_max = 120.
        legend = 'H2O (mm)'
        unit_conversion = True
    # figure out normal for other molecules
    else:
        fnm0 = line0.rstrip("\n'").lstrip("b'")
        c_a0,lon,lat = read_measures_l2(fnm0,'column_amount')
        normcol,capct = figure_out_normcol(c_a0)
        v_max = vmax 
        legend = 'x'+str(normcol)+' molec/cm2'
        unit_conversion = False

    # plot each file in list
    with open(list) as f:
         for i, line in enumerate(f):
             fnm = line.rstrip("\n'").lstrip("b'")
             print ('...processing '+fnm)
             c_a,lon,lat = read_measures_l2(fnm,'column_amount')

             if (unit_conversion):
                 data = h2o_molecpercm2_mm(c_a)
             else:
                 data = np.divide(c_a, normcol)

             # mask out border pixels for plotting
             z, lon = get_masked_data(data,lon,lat,
                           central_longitude)

             # add swath to map
             img=add_data2d_cartopy(ax,z,lon,lat,vmin=vmin,
                   vmax=v_max,method='pcolormesh',cmap=cmap)

    finish_map_detail(ax,extent)

    cbar = plt.colorbar(img,ax=ax,orientation='horizontal',
              fraction=0.03)
    cbar.set_label(legend)

    post_mapping_plot(fig,figtitle=title,filename=filename)
    return 

#------
# map measures variables on world map
#------
def measures_plot_var_1orbit(fnm,varnm,proj='PlateCarree',
        central_longitude=0.,extent=(-180.,180.,-90.,90.),
        cmap='jet',filename=None,title=' '):

    missing = -1.e30
    var,lon,lat = read_measures_l2(fnm,varnm,missing=missing) 

    if var.max() < 0:
          raise(varnm+' are all negative.')

    normcol,capct = figure_out_normcol(var)
    data = np.divide(var, normcol)
    vmax = capct / normcol * 1.5
    legend='*'+str(normcol)
 
    fig = plt.figure()
    ax, projection = pre_mapping_plot(proj,extent=extent,
          central_longitude=central_longitude)

    z, lon = get_masked_data(data,lon,lat,central_longitude)
 
    img = add_data2d_cartopy(ax, z, lon, lat, vmin=0., vmax=vmax,
          cmap=cmap,method='pcolormesh')

    cbar = plt.colorbar(img,ax=ax,orientation='horizontal',
           fraction=0.03)
    cbar.set_label(legend)

    abc0 = fnm.split('/')
    abc=abc0[-1].split('_')
    axtitle = abc[0]+' '+abc[1]+' '+varnm
    finish_map_detail(ax,extent,title=axtitle)

    figtitle=proj
    post_mapping_plot(fig,figtitle=figtitle,filename=filename)

    return 

#------
# plot histogram of variable for 1 orbit
#------
def measures_histogram_1orbit(fnm,varnm,vmin=0.,vmax=None,
        nbin=100,density=True):

    var,lon,lat=read_measures_l2(fnm,varnm)
    
    normcol,capt=figure_out_normcol(var)

    var = var / normcol
    capt = capt / normcol
 
    vmax = capt * 2.

    data = var[var>vmin]

    fig = plt.figure()

    vrange = (vmin,vmax)

    plt.hist(data,bins=nbin,range=vrange,alpha=0.3,
             histtype='bar',density=density)

    abc0=fnm.split('/')
    abc = abc0[-1].split('_')
    title=abc[0]+' '+ abc[1] +' '+varnm
    xlabel = '* '+str(normcol)

    plt.title(title)
    plt.xlabel(xlabel)

    plt.show(block=False)

#------
# overlay histograms of variable for each orbit in list
#------
def measures_histogram_norbits(list,varnm,nbin=100,
                     density=True):

    fig, ax = plt.subplots()
 
    vmin=0. 
    # get normalization factor
    f0 = open(list,'r')
    line0 = f0.readline()
    fnm0 = line0.rstrip("\n'").lstrip("b'")
    f0.close()
    var,lon,lat = read_measures_l2(fnm0,varnm)
    normcol,capt = figure_out_normcol(var)
    varnorm = var / normcol
    captnorm = capt / normcol
    vmax = captnorm * 2.
    vrange = (vmin,vmax)
    xlabel = '*'+str(normcol)
    print('normcol=',normcol)
    print('vrange=',vrange)

    with open(list) as f:
        for i,line in enumerate(f):
            fnm=line.rstrip("\n'").lstrip("b'")
            print('... processing '+fnm)
            var,lon,lat = read_measures_l2(fnm,varnm)

            varnorm = var / normcol
            data = varnorm[varnorm>vmin]
            hist,bin_edge=np.histogram(data,bins=nbin,
                         range=vrange,density=density)
            if (i == 0):
                x = (bin_edge[:-1]+bin_edge[1:])/2.

            plt.plot(x,hist)
 
    plt.title('Histogram of '+varnm)
    plt.xlabel(xlabel)
     
    plt.show(block=False)
