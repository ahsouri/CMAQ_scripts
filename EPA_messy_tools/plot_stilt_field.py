class ViewGeoObj(object):

    def __init__(self):
        
        import numpy as np
        from scipy.interpolate import CubicSpline
        
        # Define some parameters for zensun
        self.nday=np.array([  1.0,   6.0,  11.0,  16.0,  21.0,  26.0,  31.0,  36.0,  41.0,  46.0,\
                             51.0,  56.0,  61.0,  66.0,  71.0,  76.0,  81.0,  86.0,  91.0,  96.0,\
                            101.0, 106.0, 111.0, 116.0, 121.0, 126.0, 131.0, 136.0, 141.0, 146.0,\
                            151.0, 156.0, 161.0, 166.0, 171.0, 176.0, 181.0, 186.0, 191.0, 196.0,\
                            201.0, 206.0, 211.0, 216.0, 221.0, 226.0, 231.0, 236.0, 241.0, 246.0,\
                            251.0, 256.0, 261.0, 266.0, 271.0, 276.0, 281.0, 286.0, 291.0, 296.0,\
                            301.0, 306.0, 311.0, 316.0, 321.0, 326.0, 331.0, 336.0, 341.0, 346.0,\
                            351.0, 356.0, 361.0, 366.0])
        
        self.eqt=np.array([ -3.23, -5.49, -7.60, -9.48,-11.09,-12.39,-13.34,-13.95,-14.23,-14.19,\
                           -13.85,-13.22,-12.35,-11.26,-10.01, -8.64, -7.18, -5.67, -4.16, -2.69,\
                            -1.29, -0.02,  1.10,  2.05,  2.80,  3.33,  3.63,  3.68,  3.49,  3.09,\
                             2.48,  1.71,  0.79, -0.24, -1.33, -2.41, -3.45, -4.39, -5.20, -5.84,\
                            -6.28, -6.49, -6.44, -6.15, -5.60, -4.82, -3.81, -2.60, -1.19,  0.36,\
                             2.03,  3.76,  5.54,  7.31,  9.04, 10.69, 12.20, 13.53, 14.65, 15.52,\
                            16.12, 16.41, 16.36, 15.95, 15.19, 14.09, 12.67, 10.93,  8.93,  6.70,\
                             4.32,  1.86, -0.62, -3.23])

        self.dec=np.array([-23.06,-22.57,-21.91,-21.06,-20.05,-18.88,-17.57,-16.13,-14.57,-12.91,\
                           -11.16, -9.34, -7.46, -5.54, -3.59, -1.62,  0.36,  2.33,  4.28,  6.19,\
                             8.06,  9.88, 11.62, 13.29, 14.87, 16.34, 17.70, 18.94, 20.04, 21.00,\
                            21.81, 22.47, 22.95, 23.28, 23.43, 23.40, 23.21, 22.85, 22.32, 21.63,\
                            20.79, 19.80, 18.67, 17.42, 16.05, 14.57, 13.00, 11.33,  9.60,  7.80,\
                             5.95,  4.06,  2.13,  0.19, -1.75, -3.69, -5.62, -7.51, -9.36,-11.16,\
                           -12.88,-14.53,-16.07,-17.50,-18.81,-19.98,-20.99,-21.85,-22.52,-23.02,\
                           -23.33,-23.44,-23.35,-23.06])
        
        # Create spline objects
        self.cs_eqt = CubicSpline( self.nday, self.eqt )
        self.cs_dec = CubicSpline( self.nday, self.dec )
         
    def zensun(self,day,time,lat,lon,local=False):
        
        import numpy as np
        
        ####################################
        # compute the subsolar coordinates #
        ####################################
        
        # Fractional day number wirht 1 Jan 12 am = 1
        tt = ( ( np.fix(day)+time/24.-1 ) % 365.25 ) + 1.
        
        #
        eqtime = self.cs_eqt( tt ) / 60.0
        decang = self.cs_dec( tt )
        
        latsun = decang
        
        if( local ):
            #
            lonorm=((lon + 360 + 180 ) % 360 ) - 180.
            tzone=np.fix((lonorm+7.5)/15)
            index = np.where( lonorm < 0 )
            if( lonorm[index].shape[0] > 0 ):
                tzone[index] = np.fix((lonorm[index]-7.5)/15)
            ut=(time-tzone+24.) % 24.
            noon=tzone+12.-lonorm/15.
        else:
            # 
            ut=time
            noon=12.-lon/15.
        
        lonsun=-15.*(ut-12.+eqtime)

        # compute the solar zenith, azimuth and flux multiplier
        dtor = np.pi / 180.0
        t0=(90.-lat)*dtor                            # colatitude of point
        t1=(90.-latsun)*dtor                         # colatitude of sun

        p0=lon*dtor                                  # longitude of point
        p1=lonsun*dtor                               # longitude of sun

        zz=np.cos(t0)*np.cos(t1)+np.sin(t0)*np.sin(t1)*np.cos(p1-p0) # up          \
        xx=np.sin(t1)*np.sin(p1-p0)                                  # east-west    > rotated coor
        yy=np.sin(t0)*np.cos(t1)-np.cos(t0)*np.sin(t1)*np.cos(p1-p0) # north-south /
        
        
        self.saa =atan_idl(xx,yy)/dtor 
        self.sza=np.arccos(zz)/dtor                        # solar zenith

        rsun=1.-0.01673*np.cos(.9856*(tt-2.)*dtor)      # earth-sun distance in AU
        self.solfac=zz/np.power(rsun,2)                 # flux multiplier
        
        angsun = 0.
        arg=-(np.sin(angsun)+np.cos(t0)*np.cos(t1))/(np.sin(t0)*np.sin(t1))
        self.sunrise = np.zeros(arg.shape)  
        self.sunset  = np.ones(arg.shape)*24. 
        index = np.where(abs(arg) < 1.0)
        if( arg[index].shape[0] > 0 ):
            dtime=np.arccos(arg[index])/(dtor*15)
            self.sunrise[index]=noon-dtime-eqtime[index]
            self.sunset[index] =noon+dtime-eqtime[index]
            
    def zenview(self,lonp, latp, lonss, latss, h, re):
        
        import numpy as np
        
        # Degrees to radians
        DtoR = np.pi / 180.0
        
        # sin(angular radius of Earth at height h)
        srho = re/(re+h)
        
        deltaL = abs(lonss-lonp)
        cdeltaL = np.cos(deltaL * DtoR)
        clatss = np.cos(latss * DtoR)
        slatss = np.sin(latss * DtoR)
        clatp = np.cos(latp * DtoR)
        slatp = np.sin(latp * DtoR)
        
        # Find lambda, central angle of great circle arc connecting P and SSP.
        # use Law of Cosines for spherical triangle with vertices P, SSP, and North Pole.
        # sides are central angles (great circle arcs) 90-LatP, 90-LatSS, and Lambda.
        # Law of Cosines: cos(c) = cos(a) cos(b) +sin(a) sin(b) cos(C),
        # where a, b, c are the sides and C is the corner angle opposite side c.
        
        # cos(lambda)
        clambda = slatss * slatp + clatss * clatp * cdeltaL
        
        # sin(lambda) 
        slambda = np.sin(np.arccos(clambda))
        
        
        # cos phiv (Phi_V is azimuth of satellite measured from North at target P).
        # Use Law of Cosines on Spherical Triangle formed by P, North Pole, SSP.
        
        cphiv = (slatss - slatp * clambda) / ( clatp * slambda)
        if( cphiv > 1.0 ):
            cphiv = 1.0
        if( cphiv < -1.0 ):
            cphiv = -1.0
            
        phiv = np.arccos(cphiv)
        
        # tan eta
        
        taneta = srho * slambda / (1.0 - srho * clambda)
        eta = np.arctan(taneta)
        
        # cos epsilon
        
        ceps = np.sin(eta)/srho
        eps = np.arccos(ceps)
        
        self.vaa = phiv / DtoR
        if( lonp-lonss > 0.0 ):
            self.vaa = 360.0 - self.vaa
        self.vza = eps / DtoR
        
        # Check for spacecraft below horizon at target
        
        lambda1 = np.arccos(clambda) / DtoR
        lambda0 = np.arccos(srho) / DtoR
        if( lambda1 > lambda0 ):
            #print, 'WARNING: SPACECRAFT BELOW HORIZON AT TARGET'
            #print, 'Lambda  (Central Angle between Target and SSP) = ', lambda
            #print, 'Lambda0 (Central Angle Visible to Spacecraft)  = ', lambda0
            self.vza = -self.vza
        
        
        self.vza = 90.0 - self.vza

    def compute_geom(self, lon, lat, satlon, satlat, satalt, \
                         earthradius, jday, utc):
         
        import numpy as np
        
        # Zensun requires np arrays
        self.lon    = np.array( [lon]  )
        self.lat    = np.array( [lat]    )
        self.satlon = np.array( [satlon] )
        self.satlat = np.array( [satlat] )
        self.satalt = np.array( [satalt] )
        self.earthradius = np.array( [earthradius] )
        self.jday   = np.array([jday])
        self.utc    = np.array([utc] )
        
        # Compute viewing geom
        self.zenview( self.lon, self.lat, self.satlon, self.satlat, \
                      self.satalt, self.earthradius)
            
        # Compute solar geometry
        self.zensun(self.jday, self.utc, self.lat, self.lon)
            
        # Compute relative azimuth
        self.aza = self.vaa - (180. + self.saa)
        if( self.aza < 0.0 ):
            self.aza = self.aza + 360.0
        if( self.aza > 360.0 ):
            self.aza = self.aza - 360.0

def interpol_extend_dim(xmid,ymid,zmesh):

    from scipy.interpolate import RegularGridInterpolator
    import numpy as np

    dy = ymid[1]-ymid[0]
    xmid_e = np.zeros(xmid.shape[0]+2)
    xmid_e[0] = xmid[-1]-360.0
    xmid_e[1:xmid.shape[0]+1] = xmid
    xmid_e[-1] = xmid[0]+360.0
    ymid_e = np.zeros(ymid.shape[0]+2)
    ymid_e[0] = ymid[0]-dy
    ymid_e[1:ymid.shape[0]+1] = ymid
    ymid_e[-1] = ymid[-1]+dy
    
    zmesh_e = np.zeros((xmid.shape[0]+2,ymid.shape[0]+2))

    # Bottom Row
    zmesh_e[0,0] = zmesh[-1,0] ; zmesh_e[-1,0] = zmesh[0,0]
    zmesh_e[1:xmid.shape[0]+1,0] = zmesh[:,0]

    # Top Row
    zmesh_e[0,-1] = zmesh[-1,-1] ; zmesh_e[-1,-1] = zmesh[0,-1]
    zmesh_e[1:xmid.shape[0]+1,-1] = zmesh[:,-1]

    # Center
    zmesh_e[1:xmid.shape[0]+1,1:ymid.shape[0]+1] = zmesh

    # Left/Right Columns
    zmesh_e[0,1:ymid.shape[0]+1] = zmesh[-1,:]
    zmesh_e[-1,1:ymid.shape[0]+1] = zmesh[0,:]

    f = RegularGridInterpolator((xmid_e,ymid_e),zmesh_e)

    return f


def compute_edges(lon,lat):

    import numpy as np
    from matplotlib import pyplot as plt

    # First we need to extend the mesh
    imx = lon.shape[0] ; jmx = lon.shape[1]
    lon_ext = np.zeros((imx+2,jmx+2))
    lat_ext = np.zeros((imx+2,jmx+2)) 

    # Dumpy the original coordinates 1-box from the border
    lon_ext[1:imx+1,1:jmx+1] = lon 
    lat_ext[1:imx+1,1:jmx+1] = lat

    # Bottom Left Corner
    v00 = np.array([lon[0,0],lat[0,0]])
    v10 = np.array([lon[1,0],lat[1,0]])
    v01 = np.array([lon[0,0],lat[0,1]])
    dv_x = v10-v00 ; dv_y = v01 - v00
    lon_ext[0,0] = v00[0] - dv_x[0] - dv_y[0]
    lat_ext[0,0] = v00[1] - dv_x[1] - dv_y[1]

    # Bottom/Top Row
    for i in range(imx):
        v00 = np.array([lon[i,0],lat[i,0]])
        v01 = np.array([lon[i,0],lat[i,1]])
        dv_y = v01 - v00
        lon_ext[i+1,0] = v00[0] - dv_y[0]
        lat_ext[i+1,0] = v00[1] - dv_y[1]

        v00 = np.array([lon[i,jmx-1],lat[i,jmx-1]])
        v01 = np.array([lon[i,jmx-2],lat[i,jmx-2]])
        dv_y = v01 - v00
        lon_ext[i+1,jmx+1] = v00[0] - dv_y[0]
        lat_ext[i+1,jmx+1] = v00[1] - dv_y[1]

    # Bottom Right Corner
    v00 = np.array([lon[imx-1,0],lat[imx-1,0]])
    v10 = np.array([lon[imx-2,0],lat[imx-2,0]])
    v01 = np.array([lon[imx-1,0],lat[imx-1,1]])
    dv_x = v10-v00 ; dv_y = v01-v00
    lon_ext[imx+1,0] = v00[0] - dv_x[0] - dv_y[0] ; xx = lon_ext[imx+1,0] 
    lat_ext[imx+1,0] = v00[1] - dv_x[1] - dv_y[1] ; yy = lat_ext[imx+1,0]

    # Right/Left Column
    for j in range(jmx):
        v00 = np.array([lon[imx-1,j],lat[imx-1,j]])
        v10 = np.array([lon[imx-2,j],lat[imx-2,j]])
        dv_x = v10 - v00
        lon_ext[imx+1,j+1] = v00[0] - dv_x[0]
        lat_ext[imx+1,j+1] = v00[1] - dv_x[1]

        v00 = np.array([lon[0,j],lat[0,j]])
        v10 = np.array([lon[1,j],lat[1,j]])
        dv_x = v10 - v00
        lon_ext[0,j+1] = v00[0] - dv_x[0]
        lat_ext[0,j+1] = v00[1] - dv_x[1]

    
    # Top Right Corner
    v00 = np.array([lon[imx-1,jmx-1],lat[imx-1,jmx-1]])
    v10 = np.array([lon[imx-2,jmx-1],lat[imx-2,jmx-1]])
    v01 = np.array([lon[imx-1,jmx-2],lat[imx-1,jmx-2]])
    dv_x = v10-v00 ; dv_y = v01 - v00
    lon_ext[imx+1,jmx+1] = v00[0] - dv_x[0] - dv_y[0]
    lat_ext[imx+1,jmx+1] = v00[1] - dv_x[1] - dv_y[1]

    # Top Left Corner
    v00 = np.array([lon[0,jmx-1],lat[0,jmx-1]])
    v10 = np.array([lon[1,jmx-1],lat[1,jmx-1]])
    v01 = np.array([lon[0,jmx-2],lat[0,jmx-2]])
    dv_x = v10-v00 ; dv_y = v01 - v00
    lon_ext[0,jmx+1] = v00[0] - dv_x[0] - dv_y[0]
    lat_ext[0,jmx+1] = v00[1] - dv_x[1] - dv_y[1]

    # Array of centers
    lon_c = np.zeros((imx,jmx,4))
    lat_c = np.zeros((imx,jmx,4))

    # Bottom Left
    lon_c[:,:,0] = 0.25*( lon_ext[0:imx,0:jmx]+lon_ext[1:imx+1,0:jmx] \
                         +lon_ext[0:imx,1:jmx+1]+lon_ext[1:imx+1,1:jmx+1])
    lat_c[:,:,0] = 0.25*( lat_ext[0:imx,0:jmx]+lat_ext[1:imx+1,0:jmx] \
                         +lat_ext[0:imx,1:jmx+1]+lat_ext[1:imx+1,1:jmx+1])

    # Top Left
    lon_c[:,:,1] = 0.25*( lon_ext[0:imx,1:jmx+1]+lon_ext[1:imx+1,1:jmx+1] \
                         +lon_ext[0:imx,2:jmx+2]+lon_ext[1:imx+1,2:jmx+2])
    lat_c[:,:,1] = 0.25*( lat_ext[0:imx,1:jmx+1]+lat_ext[1:imx+1,1:jmx+1] \
                         +lat_ext[0:imx,2:jmx+2]+lat_ext[1:imx+1,2:jmx+2])

    # Top Right
    lon_c[:,:,2] = 0.25*( lon_ext[1:imx+1,1:jmx+1]+lon_ext[2:imx+2,1:jmx+1] \
                         +lon_ext[1:imx+1,2:jmx+2]+lon_ext[2:imx+2,2:jmx+2])
    lat_c[:,:,2] = 0.25*( lat_ext[1:imx+1,1:jmx+1]+lat_ext[2:imx+2,1:jmx+1] \
                         +lat_ext[1:imx+1,2:jmx+2]+lat_ext[2:imx+2,2:jmx+2])

    # Bottom Right
    lon_c[:,:,3] = 0.25*( lon_ext[1:imx+1,0:jmx]+lon_ext[2:imx+2,0:jmx] \
                         +lon_ext[1:imx+1,1:jmx+1]+lon_ext[2:imx+2,1:jmx+1])
    lat_c[:,:,3] = 0.25*( lat_ext[1:imx+1,0:jmx]+lat_ext[2:imx+2,0:jmx] \
                         +lat_ext[1:imx+1,1:jmx+1]+lat_ext[2:imx+2,1:jmx+1])
    return lon_c,lat_c


def read_csv(infile):

    import pysplat
    import numpy as np
    import csv

    with open(infile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        first = True
        output = {} ; keys = []
        for row in csv_reader:
            if(first):
                #print(f'Column names are {", ".join(row)}')
                for name in row:
                    output[name] = [] ; keys.append(name) ; nvar = len(keys)
                first = False
            else:
                for n in range(nvar):
                    output[keys[n]].append(row[n])
    return output

def read_orbit(infile):

    import numpy as np
    import pysplat
    d=read_csv(infile)
    fldname = ['lon_left', 'lon_center', 'lon_right', 'lat_left', 'lat_center', 'lat_right','subsat_lon', 'subsat_lat']
    output = {}
    #print(d.keys())
    for fld in fldname:
        ts = ['nan' if x == 'NA' else x for x in d[fld]] ; d[fld] = ts
        output[fld] = np.array(d[fld],dtype=np.float)

    # Process time strings (Example:'2017-05-30 12:30:14')
    nymd = [] ; nhms = []
    for tstr in d['utc']:
        frag = tstr.split(' ') ; ymd = frag[0].split('-') ; hms = frag[1].split(':')
        nymd.append(ymd[0]+ymd[1]+ymd[2]) ; nhms.append(hms[0]+hms[1]+hms[2])
    output['nymd'] = np.array(nymd,dtype=np.int)
    output['nhms'] = np.array(nhms,dtype=np.int)
    output['tau'] = [] ; output['hour'] = []
    for a,b in zip(output['nymd'],output['nhms']):
        output['tau'].append(pysplat.time.nymd2tau(a,b))
        output['hour'].append( output['tau'][-1] - pysplat.time.nymd2tau(a) )
    output['tau'] = np.array(output['tau'])
    output['hour'] = np.array(output['hour'])
    
    # Now get the pixel lons
    nxtrk = 200 ; nt = len(output['hour'])
    xedge = np.zeros((nxtrk+1,nt))
    yedge = np.zeros((nxtrk+1,nt))
    for trk in range(nxtrk):
        fld = 'lon_' + str(trk+1)
        ts = ['nan' if x == 'NA' else x for x in d[fld]] 
        xedge[trk,:] = ts
        fld = 'lat_' + str(trk+1)
        ts = ['nan' if x == 'NA' else x for x in d[fld]] 
        yedge[trk,:] = ts
    
    # Pixel Coordinates
    output['pix_lon_cnr'] = np.zeros((nxtrk,nt-1,4))
    output['pix_lon_cnr'][:,:,0] = xedge[0:nxtrk,0:nt-1]
    output['pix_lon_cnr'][:,:,1] = xedge[0:nxtrk,1:nt]
    output['pix_lon_cnr'][:,:,2] = xedge[1:nxtrk+1,1:nt]
    output['pix_lon_cnr'][:,:,3] = xedge[1:nxtrk+1,0:nt-1]
    output['pix_lat_cnr'] = np.zeros((nxtrk,nt-1,4))
    output['pix_lat_cnr'][:,:,0] = yedge[0:nxtrk,0:nt-1]
    output['pix_lat_cnr'][:,:,1] = yedge[0:nxtrk,1:nt]
    output['pix_lat_cnr'][:,:,2] = yedge[1:nxtrk+1,1:nt]
    output['pix_lat_cnr'][:,:,3] = yedge[1:nxtrk+1,0:nt-1]
    output['pix_lon'] = np.sum(output['pix_lon_cnr'],axis=2)/4.0
    output['pix_lat'] = np.sum(output['pix_lat_cnr'],axis=2)/4.0
    output['pix_hour'] = 0.5*(output['hour'][0:nt-1] + output['hour'][1:nt])

    return output

# def read_scan(infile):

#     import numpy as np
#     d = read_csv(infile)

#     return {'lon':np.array(d['lon'],dtype=np.float),\
#             'lat':np.array(d['lat'],dtype=np.float),\
#             'time':np.array(d['time'],dtype=np.float)}

def read_3d_fld(infile):

    # Read the csv
    d = read_csv(infile) ; print(d.keys())

    # Get the grid
    output = {}
    output['lon'] = np.sort(np.array(np.unique(d['longitude']),dtype=np.float)) ; imx = len(output['lon'])
    output['lat'] = np.sort(np.array(np.unique(d['latitude']),dtype=np.float)) ; jmx = len(output['lat'])
    output['lev'] = np.sort(np.array(np.unique(d['level']),dtype=np.int)) ; lmx = len(output['lev'])
    npt = len(d['longitude'])
    
   

    # Output Variables 
    fldname = ['agl_m','aglprecise_m','temperature', 'pres_hPa', 'sphu_kgkg', 'uwnd_ms', 'vwnd_ms', 'ch4_ppb']
    for name in fldname:
        
        # Replace NA with nans
        ts = ['nan' if x == 'NA' else x for x in d[name]]
        d[name] = ts

        # Convert to float and reshape
        output[name] = np.array(ts,dtype=np.float).reshape((lmx,jmx,imx)).T
        
    return output

from pylab import *
from carbontracker import carbontracker, carbontracker_co2
from vertical_regrid_weights import vertical_regrid_weights
from scipy import constants
from scipy.interpolate import interp1d
import pysplat

# ===============================================================================
# INPUTS
# ===============================================================================

# Satellite information
orbit_file = 'datasets/STILT/Indianapolis/Orbit_17_5Hz_Indianapolis.csv'
stilt_file = 'datasets/STILT/Indianapolis/Fields_3d_Indianapolis_HRRR.csv'

# RTM Vertical grid
lmx = 20
eta_a = np.zeros(lmx+1)
eta_b = np.linspace(0.0,1.0,lmx+1)[::-1] ; eta_b[-1]=0.001
l1_outfile = 'MethaneSAT_L1_Indianapolis_20170506.nc'

# Wavelength grid
wvl_o2 = np.arange(1246.0,1293.0,0.05)
wvl_ch4 = np.arange(1606.0,1689.0,0.063)

# ===============================================================================

# Constants
MW_air = 28.97e-3        # Molecular weight dry air [kg/mol]
MW_h2o = 18.01528e-3     # Molecular weight water [kg/mol]
R = constants.R          # Universal Gas Constant [J/K/mol]
Na = constants.Avogadro  # Avogadros Number [molec/mol]
g = constants.g          # Surface Gravity [m/s^2]

# Concnentrations o
xn2 = 0.78084
xo2 = 0.20946
xar = 0.09340

# Read Datasets
orbit = read_orbit(orbit_file)
fld = read_3d_fld(stilt_file)

xmesh,ymesh = np.meshgrid( fld['lon'], fld['lat'], indexing='ij')

# Get Cross Track / Along Track indices
xtrk_id,atrk_id = np.meshgrid(np.arange(0, orbit['pix_lon'].shape[0],dtype=np.int),\
                              np.arange(0, orbit['pix_lon'].shape[1],dtype=np.int),\
                              indexing='ij')

# Get the field limits
lon0 = fld['lon'].min() ; lonf = fld['lon'].max()
lat0 = fld['lat'].min() ; latf = fld['lat'].max()
lon_c = fld['lon'].mean() ; lat_c = fld['lat'].mean() 

# Lets estimate the edges
lon_c,lat_c = compute_edges(orbit['pix_lon'],orbit['pix_lat'])

# Get pixels within range
pix_mask = np.logical_and.reduce([orbit['pix_lon']>=lon0,\
                                  orbit['pix_lon']<=lonf,\
                                  orbit['pix_lat']>=lat0,\
                                  orbit['pix_lat']<=latf])

# Get the overlapping indices
imin = xtrk_id[pix_mask].min() ; imax = xtrk_id[pix_mask].max() ; imx = imax-imin+1
jmin = atrk_id[pix_mask].min() ; jmax = atrk_id[pix_mask].max() ; jmx = jmax-jmin+1


fig = plt.figure()
plt.pcolormesh(xmesh,ymesh,fld['sphu_kgkg'][:,:,0])
plt.colorbar()
plt.scatter(orbit['pix_lon'][pix_mask],orbit['pix_lat'][pix_mask],alpha=0.2,marker='.',color='red')
fig.show()