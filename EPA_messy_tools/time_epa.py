def nymd2str(nymd,nhms=None):

    import numpy as np
    # Ensure numpy array
    npymd = np.array(nymd) ; nphms = np.array(nhms)
    scalar_ymd = False
    if(len(npymd.shape) == 0):
        npymd = np.array([nymd])
        scalar_ymd = True
    Year = [] ; Month = []; Day = []
    for ymd in npymd:
        # Convert to string
        ymdstr = str(ymd).zfill(8)
        # Get YMD
         # Get YMD
        Year.append(int(ymdstr[0:4]))
        Month.append(int(ymdstr[4:6]))
        Day.append(int(ymdstr[6:8]))
    d = {}
    d['Year'] = np.array(Year) ; d['Month'] = np.array(Month) ; d['Day'] = np.array(Day)
    nt = len(Year)
    if(nhms is None):
        if(scalar_ymd):
            z = 0
        else:
            z = np.zeros(d['Year'].shape)
        d['Hour'] = z ; d['Minute'] = z ; d['Second'] = z
    else:
        scalar_hms = False
        if(len(nphms.shape) == 0):
            nphms = np.array([nhms])
            scalar_hms = True
        Hour = [] ; Minute = []; Second = []
        for hms in nphms:
            # Convert to string
            hmsstr = str(hms).zfill(6)

            # Get YMD
            Hour.append(int(hmsstr[0:2]))
            Minute.append(int(hmsstr[2:4]))
            Second.append(int(hmsstr[4:6]))

        if( len(Hour) < nt):
            d['Hour'] = np.ones(nt)*Hour[0]
            d['Minute'] = np.ones(nt)*Minute[0]
            d['Second'] = np.ones(nt)*Second[0]
        else:
            d['Hour'] = np.array(Hour[0:nt])
            d['Minute'] = np.array(Minute[0:nt])
            d['Second'] = np.array(Second[0:nt])
    return d

def nymd2tau(nymd,nhms=000000,nymd0=19850101,nhms0=000000):

    ''' Computes the value of TAU, which is the elapsed hours between
        the current date/time and the beginning of an epoch. This is
        the starting date for 3D model data assimilation.

        tau = nymd2tau(nymd,nhms,nymd0,nhms0)

        ARGS:
            nymd (long)  -> YYYY/MM/DD for this date (e.g. 19940101)
        OPTIONAL:
            nhms (long)  -> HH/MM/SS for this date (e.g. 060000)
                      will be defaulted to 000000
    
            nymd0 (long) -> YY/MM/DD for the start of the epoch
                      default is {19}850101 which is the GEOS-1 start
    
            nhms0 (long) -> HH/MM/SS for the start of the epoch
                      will be defaulted to 000000
        RETURN:
            Tau -> The function returns the TAU value as a 
;             double-precision number
    '''
    import numpy as np
    from datetime import date

    # Separate YMD 
    d0 = nymd2str(nymd0,nhms0)
    d  = nymd2str(nymd,nhms)
     # Set start date
    date0 = date(d0['Year'][0],d0['Month'][0],d0['Day'][0])

    # Compute fractional hour
    hr0 = d0['Hour'][0] + d0['Minute'][0]/60.0 + d0['Second'][0]/3600.0

    # Compute tau time
    tau = []
    for n in range(d0['Year'].shape[0]):

        # Compute date for nth time
        date1 = date(d['Year'][n],d['Month'][n],d['Day'][n])
        hr1 = d['Hour'][n] + d['Minute'][n]/60.0 + d['Second'][n]/3600.0

        # Compute difference in hours
        tau.append(float( (date1 - date0).days )*24.0+hr1-hr0)

    return np.array(tau) 
def tau2yymmdd(tau,nymd0=19850101,nhms0=000000):


    ''' Converts a tau value (elapsed hours between
        the current date/time and the beginning of an epoch) into a
        calendar date and time value.

        ARGS:
            tau -> the tau value to be converted (type long)
        OPTIONAL:
            nymd0 (long) -> YY/MM/DD for the start of the epoch
                   default is 19850101 which is the GEOS-1 start
            nhms0 (long) -> HH/MM/SS for the start of the epoch
                   will be defaulted to 000000
    '''

    import numpy as np
    from datetime import timedelta
    from datetime import datetime
    pass

    # Compute the date time for the start of the epoch
    d0 = nymd2str(nymd0,nhms0)
    time0 = datetime(d0['Year'][0],d0['Month'][0],d0['Day'][0],\
                     d0['Hour'][0],d0['Minute'][0],d0['Second'][0])

    # Check tau
    nptau = np.array(tau)
    if(len(nptau.shape) == 0):
        nptau = np.array([nptau])

    # Output
    d = {}

    Year = [] ; Month = [] ; Day = []
    Hour = [] ; Minute = [] ; Second  = []

    for n in range(nptau.shape[0]):
    # Compute the new date
        time = time0 + timedelta(seconds=nptau[n]*3600.0)
        Year.append(time.year)
        Month.append(time.month)
        Day.append(time.day)
        Hour.append(time.hour)
        Minute.append(time.minute)
        Second.append(time.second)

    d = {}
    d['Year'] = np.array(Year) ; d['Month'] = np.array(Month)
    d['Day'] = np.array(Day) ; d['Hour'] = np.array(Hour)
    d['Minute'] = np.array(Minute) ; d['Second'] = np.array(Second)

    return d