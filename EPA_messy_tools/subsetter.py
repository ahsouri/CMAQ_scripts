def subseter(data,yyyy1mm1dd1,yyyy2mm2dd2,hh1mm1ss1,hh2mm2ss2,lat1,lon1,lat2,lon2,state=None):
    # start date: yyyy1mm1dd1hh1
    # end date: yyyy2mm2dd2hh2
    # convert start/end date to tau_data
    from time_epa import nymd2str,nymd2tau,tau2yymmdd
    import numpy as np
    start_date = nymd2tau(int(yyyy1mm1dd1),int(hh1mm1ss1))
    end_date = nymd2tau(int(yyyy2mm2dd2),int(hh2mm2ss2))
    data_subset = []
    date_subset = []
    lat_subset = []
    lon_subset = []

    for idx in range(len(data['Date Local'])):

        date_data = data['Date Local'][idx]
        time_date = data['Time Local'][idx]
        #date_data = date_data.replace("-","")
        date_data = int(date_data)
        time_date = int(time_date) 
        #time_date = time_date.replace(":","")

        #if len(time_date)==4:
        #    time_date = time_date + "00"
        #print(time_date)
        tau_data = nymd2tau(int(date_data),int(time_date))
        if tau_data>=start_date and tau_data<=end_date:
            if float(data['Latitude'][idx]) >= lat1 and float(data['Latitude'][idx]) <= lat2 and \
               float(data['Longitude'][idx]) >= lon1 and float(data['Longitude'][idx]) <= lon2:
               if not (state is None):
                  if data['State Name'][idx] == state:
                     data_subset.append(float(data['Sample Measurement'][idx]))
                     date_subset.append(tau_data)
                     lat_subset.append(float(data['Latitude'][idx]))
                     lon_subset.append(float(data['Longitude'][idx]))
               else:
                     data_subset.append(float(data['Sample Measurement'][idx]))
                     date_subset.append(tau_data)
                     lat_subset.append(float(data['Latitude'][idx]))
                     lon_subset.append(float(data['Longitude'][idx]))
    # uniqueness in time and space

    data_subset_tmp = np.array(data_subset)
    date_subset_tmp = np.array(date_subset).squeeze()
    lat_subset_tmp = np.array(lat_subset)
    lon_subset_tmp = np.array(lon_subset)

    sites_unique = np.unique(lat_subset)


    output_subsetter = {} ; keys = []
    num_station = len(sites_unique)

    for i in range(num_station):

        keys.append('station' + str(i) + '_O3')
        keys.append('station' + str(i) + '_time')
        keys.append('station' + str(i) + '_lat')
        keys.append('station' + str(i) + '_lon')

        mask = lat_subset_tmp == sites_unique[i] 

        output_subsetter[keys[i*4+0]] = data_subset_tmp[mask]
        output_subsetter[keys[i*4+1]] = date_subset_tmp[mask]
        output_subsetter[keys[i*4+2]] = np.nanmean(lat_subset_tmp[mask])
        output_subsetter[keys[i*4+3]] = np.nanmean(lon_subset_tmp[mask])
            

    return output_subsetter,num_station
