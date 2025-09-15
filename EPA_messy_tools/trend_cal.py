def func_season(x, a, b, c, d):
    import numpy as np
    return a + b*(x-np.min(x)) + c*np.cos(2*np.pi*1*(x-d))

def func_year(x, a, b):
    import numpy as np
    return a + b*(x-np.min(x)) 
    
def trend_cal(xdata,ydata,season=False):
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    import numpy as np
    import pymannkendall as mk


    xdata[np.isnan(ydata)] = np.nan
    ydata[np.isnan(xdata)] = np.nan
    xdata = xdata[~np.isnan(xdata)]
    ydata = ydata[~np.isnan(ydata)]

    if season:
       popt, pcov = curve_fit(func_season, xdata, ydata)
       ydata_desea = ydata - popt[2]*np.cos(2*np.pi*1*(xdata-popt[3]))
    
       #mk_output = mk.original_test(ydata_desea, alpha=0.05)
       plt.plot(xdata, ydata_desea, 'b-', label='data')
       plt.plot(xdata, func_season(xdata, *popt), 'r--',label='fit: Avg=%5.3f, Trend=%5.3f, Seasonal=%5.3f, Phase=%5.3f' % tuple(popt))
       plt.legend(loc='best')
       plt.show()
    else:
       popt, pcov = curve_fit(func_year, xdata, ydata) 
       mk_output = mk.original_test(ydata, alpha=0.05)
       plt.plot(xdata, ydata, 'b-', label='data')
       plt.plot(xdata, func_year(xdata, *popt), 'r--',label='fit: Avg=%5.3f, Trend=%5.3f' % tuple(popt))
       plt.legend(loc='best')
       plt.show()
       print(mk_output)

    



