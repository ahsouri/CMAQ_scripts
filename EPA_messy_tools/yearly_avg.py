def leap_year(y):
    if y % 400 == 0:
        return True
    if y % 100 == 0:
        return False
    if y % 4 == 0:
        return True
    else:
        return False

def yearly_avg(dates,data,months):

    import numpy as np
    import datetime


    yearfrc = np.array(dates)
    month = np.floor((yearfrc-np.floor(yearfrc))*12)+1

    min_yr = np.floor(np.min(yearfrc))
    max_yr = np.floor(np.max(yearfrc))


    month_flag = [n in months for n in month]
 
    data_yearly_avg = [] 
    data_yearly_std = []
    year = []
    for yr in range(int(min_yr),int(max_yr)):
        temp_yr = data [ (yearfrc>=yr) & (yearfrc<yr+1) & month_flag]
        temp_yr [temp_yr<=0] = np.nan
        data_yearly_avg.append(np.nanmean(temp_yr))
        data_yearly_std.append(np.nanstd(temp_yr))
        year.append(float(yr))
   
    return np.array(data_yearly_avg),np.array(data_yearly_std),np.array(year)