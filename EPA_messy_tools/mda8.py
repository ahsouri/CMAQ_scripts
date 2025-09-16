def leap_year(y):
    if y % 400 == 0:
        return True
    if y % 100 == 0:
        return False
    if y % 4 == 0:
        return True
    else:
        return False

def mda8(dates,o3):
    import numpy as np
    import datetime
    

    hours = []
    yeardoy = []

    for t in range(len(dates)):

        date_data = dates[t]
        doy = date_data.strftime('%j')
        year = date_data.year
        if leap_year(year):
            yeardoy.append(year + float(doy)/366)
        else:
            yeardoy.append(year + float(doy)/365)

        hours.append(date_data.hour)
    
    yeardoy = np.array(yeardoy)
    year = []
    hours = np.array(hours)
    MDA8 = []
    for t in range(np.shape(yeardoy)[0]):
        ozone_for_a_day = o3[yeardoy == yeardoy[t]]
        hours_for_a_day = hours[yeardoy == yeardoy[t]]
        avg_ozone = []
        for hr in range(np.size(hours_for_a_day)):
            if hr+8>np.size(hours_for_a_day): continue
            temp = ozone_for_a_day[hr:hr+8]
            temp[temp <= 0] = np.nan
            avg_ozone.append(np.mean(temp))
        if avg_ozone:
           MDA8.append(max(avg_ozone))
        else:
           MDA8.append(np.nan)
        year.append(yeardoy[t]/1000)
        
    return yeardoy,np.array(MDA8)


