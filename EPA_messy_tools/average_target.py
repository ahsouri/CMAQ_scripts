def average_target(variable,num_station,num_year,startyear,thr_yr):
    import numpy as np 
    from datetime import datetime,timedelta
    import pandas as pd
    from time_epa import tau2yymmdd
    from mda8 import mda8
    from yearly_avg import yearly_avg
    from monthly_avg import monthly_average
    import matplotlib.pyplot as plt
    # averaging over location
    x_all_stations = np.zeros((int(num_year),int(num_station)))*np.nan
    y_all_stations = np.zeros((int(num_year),int(num_station)))*np.nan
    std_all_stations = np.zeros((int(num_year),int(num_station)))*np.nan
    month_all_stations = np.zeros((int(12),int(num_station)))*np.nan
    counter = 0
    for nsta in range(num_station):
        print('reading station ' + str(nsta+1) + '/' + str(num_station))
        var_value = variable['station' + str(nsta) + '_O3']
        dates = variable['station' + str(nsta) + '_time']
        var_value[var_value<=0] = np.nan

        dates_f=[]
        for i in range(np.size(dates)):
            if np.isnan(dates[i]):
               dates_f.append(pd.NaT)
            else:
               dates_tmp = tau2yymmdd(dates[i])
               p=datetime(dates_tmp['Year'][0], dates_tmp['Month'][0], dates_tmp['Day'][0], dates_tmp['Hour'][0], \
                   dates_tmp['Minute'][0], dates_tmp['Second'][0])
               dates_f.append(p)


        xData,yData = mda8(dates_f,var_value*1000.0)

        monthly_mean_data = monthly_average(xData,yData)
        
        data_yearly_avg,data_yearly_std,year = yearly_avg(xData,yData,np.array(range(3,8)))

        if np.size(year) == 0:
           continue
        
        min_yr = np.floor(np.min(year))
        max_yr = np.floor(np.max(year))

        if max_yr - min_yr > thr_yr:
           
           counter = counter + 1
           month_all_stations[:,counter] = monthly_mean_data
           x_all_stations[int(min_yr-startyear):int(max_yr-startyear)+1,counter] = year
           y_all_stations[int(min_yr-startyear):int(max_yr-startyear)+1,counter] = data_yearly_avg
           std_all_stations[int(min_yr-startyear):int(max_yr-startyear)+1,counter] = data_yearly_std
    
    #plt.plot(month_all_stations, 'b-', label='data')
    #plt.show()
    return x_all_stations,y_all_stations,std_all_stations