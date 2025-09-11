import numpy as np
import datetime
from dataclasses import dataclass
from netCDF4 import Dataset
import warnings
import glob
from scipy import interpolate
from scipy.io import savemat
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning)


@dataclass
class aircraft_data:
    lat: np.ndarray
    lon: np.ndarray
    time: datetime.datetime
    HCHO_ppbv: np.ndarray
    NO2_ppbv: np.ndarray
    altp: np.ndarray
    profile_num: np.ndarray

@dataclass
class surface_data:
    lat: np.ndarray
    lon: np.ndarray
    time: datetime.datetime
    value: np.ndarray
    state: list


@dataclass
class ctm_model:
    latitude: np.ndarray
    longitude: np.ndarray
    time: list
    gas_profile_no2: np.ndarray
    gas_profile_hcho: np.ndarray
    pressure_mid: np.ndarray
    pressure_edge: np.ndarray
    ZL: np.ndarray
    PBLH: np.ndarray
    TEMP2: np.ndarray
    ctmtype: str


def _read_nc(filename, var):
    # reading nc files without a group
    nc_f = filename
    nc_fid = Dataset(nc_f, 'r')
    out = np.array(nc_fid.variables[var])
    nc_fid.close()
    return np.squeeze(out)


def _get_nc_attr(filename, var):
    # getting attributes
    nc_f = filename
    nc_fid = Dataset(nc_f, 'r')
    attr = {}
    for attrname in nc_fid.variables[var].ncattrs():
        attr[attrname] = getattr(nc_fid.variables[var], attrname)
    nc_fid.close()
    return attr


def aircraft_reader(filename: str, year: int, NO2_string: str, HCHO_string: str):
    # read the header
    with open(filename) as f:
        header = f.readline().split(',')
    # read the data
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    doy = data[:, header.index(' JDAY')]
    utc_time = data[:, header.index('  UTC')]
    lat = data[:, header.index(' LATITUDE')]
    lon = data[:, header.index(' LONGITUDE')]
    lon[lon > 180.0] = lon[lon > 180.0] - 360.0
    HCHO_ppbv = data[:, header.index(' ' + HCHO_string)]/1000.0
    HCHO_ppbv[HCHO_ppbv <= 0] = np.nan
    NO2_ppbv = data[:, header.index(' ' + NO2_string)]/1000.0
    NO2_ppbv[NO2_ppbv <= 0] = np.nan
    altp = data[:, header.index(' PRESSURE')]
    altp[altp <= 0] = np.nan
    profile_num = data[:, header.index(' ProfileNumber')]
    # convert doy and utc to datetime
    date = []
    for i in range(0, np.size(doy)):
        date.append(datetime.datetime.strptime(str(int(year)) + "-" + str(int(doy[i])), "%Y-%j") +
                    datetime.timedelta(seconds=int(utc_time[i])))
    # return a structure
    return aircraft_data(lat, lon, date, HCHO_ppbv, NO2_ppbv, altp, profile_num)

def AQS_reader(filename: str):
    # read the header
    data = pd.read_csv(filename)
    # read the data
    lat = np.array(data["Latitude"])
    lon = np.array(data["Longitude"])
    lon[lon > 180.0] = lon[lon > 180.0] - 360.0
    airtemp = np.array(data["Sample Measurement"]) 
    airtemp = (airtemp - 32.0)*(5.0/9.0) # far to deg
    # convert doy and utc to datetime
    data["DateTime_GTM"] = pd.to_datetime(data["Date_GTM"] + " " + data["Time_GTM"])
    date = pd.to_datetime(data["DateTime_GTM"]).to_list()
    state = data["State Code"].to_list()
    # return a structure
    return surface_data(lat, lon, date, airtemp, state)

def cmaq_reader_core(cmaq_target_file,met_file_3d_file,met_file_2d_file,grd_file_2d_file):
        
        print("Currently reading: " + cmaq_target_file.split('/')[-1])
        # reading time and coordinates
        lat = _read_nc(grd_file_2d_file, 'LAT')
        lon = _read_nc(grd_file_2d_file, 'LON')
        time_var = _read_nc(cmaq_target_file, 'TFLAG')
        # populating cmaq time
        time = []
        for t in range(0, np.shape(time_var)[0]):
            cmaq_date = datetime.datetime.strptime(
                str(time_var[t, 0, 0]), '%Y%j').date()
            time.append(datetime.datetime(int(cmaq_date.strftime('%Y')), int(cmaq_date.strftime('%m')),
                                      int(cmaq_date.strftime('%d')), int(time_var[t, 0, 1]/10000.0), 0, 0) +
                    datetime.timedelta(minutes=0))

        #prs = _read_nc(met_file_3d_file, 'PRES').astype('float32')/100.0  # hPa
        PBLH = _read_nc(met_file_2d_file, 'PBL').astype('float32')
        ZL = _read_nc(met_file_3d_file, 'ZH').astype('float32')
        TEMP2 = _read_nc(met_file_2d_file,'TEMP2').astype('float32')-273.15
        # read gas in ppbv
        #gas_hcho = _read_nc(cmaq_target_file, 'FORM')*1000.0  # ppb
        #gas_hcho = gas_hcho.astype('float32')
        #gas_no2 = _read_nc(cmaq_target_file, 'NO2')*1000.0  # ppb
        #gas_no2 = gas_no2.astype('float32')        
        # populate cmaq_data format
        cmaq_data = ctm_model(lat, lon, time, [], [], [], [], ZL, PBLH, TEMP2, 'CMAQ')
        return cmaq_data

def cmaq_reader(dir_mcip: str, dir_cmaq: str, YYYYMM: str) -> list:
    '''
       MINDS reader
       Inputs:
             product_dir [str]: the folder containing the GMI data
             YYYYMM [str]: the target month and year, e.g., 202005 (May 2020)
       Output:
             minds_fields [ctm_model]: a dataclass format
    '''
    # read meteorological and chemical fields
    cmaq_target_files = sorted(glob.glob(dir_cmaq + "/CCTM_CONC_*" + YYYYMM +  "*.nc"))
    grd_files_2d = sorted(
            glob.glob(dir_mcip + "/GRIDCRO2D_*" + \
        YYYYMM +  "*"))
    met_files_2d = sorted(
            glob.glob(dir_mcip + "/METCRO2D_*" + YYYYMM  + "*"))
    met_files_3d = sorted(
            glob.glob(dir_mcip + "/METCRO3D_*" + YYYYMM  + "*"))
    if len(cmaq_target_files) != len(met_files_3d):
            raise Exception(
                "the data are not consistent")

    # define gas profiles for saving
    outputs = []
    for k in range(len(met_files_3d)):
        outputs.append(cmaq_reader_core(cmaq_target_files[k],met_files_3d[k],met_files_2d[k],grd_files_2d[k]))

    return outputs

def colocate(ctmdata, airdata, date1, date2):
    '''
       Colocate the model and the aircraft based on the nearest neighbour
       Inputs:
             ctmdata: structure for the model
             airdata: structure for the aircraft dataset
       Output:
             a dict containing colocated model and aircraft
    '''
    print('Colocating Aircraft and CMAQ...')
    # list the time in ctm_data
    time_ctm = []
    time_ctm_datetype = []
    for ctm_granule in ctmdata:
        time_temp = ctm_granule.time
        for n in range(len(time_temp)):
            time_temp2 = time_temp[n].year*10000 + time_temp[n].month*100 +\
                time_temp[n].day + time_temp[n].hour/24.0 + \
                time_temp[n].minute/60.0/24.0 + time_temp[n].second/3600.0/24.0
            time_ctm.append(time_temp2)
        time_ctm_datetype.append(ctm_granule.time)

    time_ctm = np.array(time_ctm)

    time_surface = []
    for t in airdata.time:
        time_surface.append(t.year*10000 + t.month*100 +
                             t.day + t.hour/24.0 + t.minute /
                             60.0/24.0 + t.second/3600.0/24.0)

    time_surface = np.array(time_surface)
    # translating ctm data into aircraft location/time
    ctm_mapped_surface_value = []
    for t1 in range(0, np.size(airdata.time)):
        if np.isnan(airdata.lat[t1]):
           ctm_mapped_surface_value.append(-999.0)
           continue
        if (airdata.time[t1]<pd.Timestamp(date1) or airdata.time[t1]<pd.Timestamp(date2)):
            continue
        if (airdata.state[t1] != 'Tennessee'):
            continue
        # for t1 in range(0, 2500):
        # find the closest day
        closest_index = np.argmin(np.abs(time_surface[t1] - time_ctm))
        # find the closest hour (this only works for 3-hourly frequency)
        closest_index_day = int(np.floor(closest_index/25.0))
        closest_index_hour = int(closest_index % 25)
        print("The closest CMAQ file used for the aircraft at " + str(airdata.time[t1]) +
              " is at " + str(time_ctm_datetype[closest_index_day][closest_index_hour]))
        ctm_surface = ctmdata[closest_index_day].TEMP2[closest_index_hour, :, :].squeeze()
        # pinpoint the closest location based on NN
        cost = np.sqrt((airdata.lon[t1]-ctmdata[0].longitude)
                       ** 2 + (airdata.lat[t1]-ctmdata[0].latitude)**2)
        index_i, index_j = np.where(cost == min(cost.flatten()))
        # picking the right grid
        ctm_surface_new = ctm_surface[index_i, index_j].squeeze(
        )

        # in case we have more than one grid box cloes to the aircraft point
        if np.size(ctm_surface_new) > 1:
            ctm_surface_new = np.mean(ctm_surface_new)

        ctm_mapped_surface_value.append(ctm_surface_new)

    # converting the lists to numpy array
    ctm_mapped_surface_value = np.array(ctm_mapped_surface_value)


    # now colocating based on each spiral number
    output = {}
    output["MCIP_TEMP2"] = ctm_mapped_surface_value
    output["AQS TEMP2"] = airdata.value
    output["time"] = airdata.time
    output["lon"] = airdata.lon
    output["lat"] = airdata.lat
    output["state"] = airdata.lat.state
    savemat('./AQS_test.mat', output)
    return None


if __name__ == "__main__":
    AQS_files = sorted(glob.glob("./aeromma_data/A*.csv"))
    for file in AQS_files:
        try:
           aircraft_data1 = AQS_reader(file)
           split_file = file.split('_')
           cmaq_data = cmaq_reader('./cmaq_mcip/','./cmaq_mcip/', "202308")
           spiral_data = colocate(cmaq_data, aircraft_data1,"2023-08-01","2023-09-01")
           spiral_data = []
           cmaq_data = []
           aircraft_data1 = []
        except:
           print("something bad happened")