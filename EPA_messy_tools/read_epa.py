from read_csv import read_csv
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from time_epa import tau2yymmdd
from subsetter import subseter
from diurnal_cycle import diurnal_cycle
from matplotlib import cm
from trend_cal import trend_cal
from mda8 import mda8
from yearly_avg import yearly_avg
from scipy.io import savemat,loadmat
from average_target import average_target

subsets = list(['Site Num','Latitude','Longitude','Date Local','Time Local','Sample Measurement' , \
    'State Name', 'County Name'])

list_files = []
for yr in range(2000,2001):
    list_files.append('/Volumes/My Passport/SAO/Ozone_Trend_Louise/hourly_44201_' + str(yr) + '.csv')


    output_o3 = read_csv(list_files,subsets,"South Dakota")
#savemat("./ozone_south_dakota_2000_2020.mat", output_o3)
#exit()
output_o3 = loadmat("./ozone_south_dakota_2000_2020.mat")

date_start = 20000101
date_end = 20201230
num_year = np.floor(date_end/1e4) - np.floor(date_start/1e4) + 1
startyear = np.floor(date_start/1e4)
threshold_year = 10

# subset based on time and location
output_subsetter,num_station = subseter(output_o3,date_start,date_end,000000,240000,-90,-180,90,180)

x_all_stations,y_all_stations,std_all_stations = average_target(output_subsetter,num_station,num_year,
                                                                  startyear,threshold_year)


print(x_all_stations)
print(y_all_stations)


year = np.nanmean(x_all_stations,axis=1)
data_yearly_avg = np.nanmean(y_all_stations,axis=1)


trend_cal(year,data_yearly_avg,season=False)

exit()

data_hourly = diurnal_cycle(o3*1000.0,dates_f)

fig, ax = plt.subplots()
plt.style.use('seaborn')

bp = ax.boxplot(data_hourly,showfliers=False, patch_artist=True) 

cmap = cm.get_cmap('jet')
colors = []
for i in range(len(data_hourly)):
    c1 = (np.nanmean(data_hourly[i])/40)
    colors.append(c1)
colors = np.array(colors)
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(cm.jet(color))
# show plot 

ax.set_title('Surface Ozone')
ax.set_xlabel('Time [LST]')
ax.set_ylabel('[ppbv]')

plt.show()

exit()

plt.style.use('seaborn')

plt.plot_date(dates_july, o3_july*1000,'r-',linewidth=3,label= 'Jul')
#plt.plot_date(dates_jan, o3_jan*1000,'b-',linewidth=3,label= 'Jan')

plt.title('Surface Ozone')
plt.xlabel('Date')
plt.ylabel('[ppbv]')
plt.legend()
plt.show()



