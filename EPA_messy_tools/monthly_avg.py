def monthly_average(xData,yData):
    import numpy as np 
    #list month names separate into bucket    
    months=["Jan","Feb","Mar","April","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    mda8_monthly_bucket = {}
    for name in months:
        mda8_monthly_bucket[name] = []
    
    for x in range(np.size(xData)):
        month = np.floor((xData[x]-np.floor(xData[x]))*12) 
        mda8_monthly_bucket[months[int(month)]].append(yData[x])

    mda8_monthly = np.zeros(12)
    for i in range(0,12):
        mda8_monthly[i] = np.nanmean(mda8_monthly_bucket[months[i]])

    return mda8_monthly
 



    
