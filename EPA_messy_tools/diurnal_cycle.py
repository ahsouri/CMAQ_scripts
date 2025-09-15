def diurnal_cycle(data,date_data):

    import numpy as np
     
    data_hourly = []

    for hr in range(0,24):
        tmp_data = []
        for i in range(np.size(data)):
            if date_data[i].hour == hr:
                tmp_data.append(data[i])
        tmp_data = np.array(tmp_data)
        tmp_data = tmp_data[tmp_data > 0 & ~np.isnan(tmp_data)]      
        data_hourly.append(tmp_data)
        
    return data_hourly
