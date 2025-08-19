import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.colors import LogNorm
from scipy import odr
from netCDF4 import Dataset
import numpy as np
from scipy.stats import gaussian_kde
from scipy import stats
import math
from york_fit import bivariate_fit
from sklearn.linear_model import TheilSenRegressor,LinearRegression
# -------------------------------
# Example data
# -------------------------------
nc_f = './NO2_202307.nc'
nc_fid = Dataset(nc_f,'r')
no2_cmaq = nc_fid.variables["ctm_averaged_vcd_prior"]
no2_trop = nc_fid.variables["sat_averaged_vcd"]
no2_trop_std = nc_fid.variables["sat_averaged_error"]
x = np.array(no2_trop).flatten()
y = np.array(no2_cmaq).flatten()
mask = ~np.isnan(x) & ~np.isnan(y)
x = x[mask]
y = y[mask]
error = np.array(no2_trop_std).flatten()
error_x = error[mask]
# -------------------------------
# Metrics
# -------------------------------
rmse = np.sqrt(np.mean((y - x)**2))
mean_bias = np.mean(y - x)
mean_abs_bias = np.mean(np.abs(y - x))
index_of_agreement = 1 - (np.sum((y - x)**2) /
                          np.sum((np.abs(y - np.mean(x)) + np.abs(x - np.mean(x)))**2))

X = x.reshape(-1, 1)
y_model = y
linear_model = LinearRegression()
linear_model.fit(X, y_model)
slope = linear_model.coef_[0]
intercept = linear_model.intercept_
r_value, _ = stats.pearsonr(x, y)
r2 = r_value**2

# york fit
#error_y = 0.1*y
#intercept, slope, _, _ = bivariate_fit(x, y, error_x, error_y, ri=0.0, b0=1.0, maxIter=1e6)

# Subsample points for KDE to make it faster
#sub_idx = np.random.choice(len(x), size=5000, replace=False)
#xy_sub = np.vstack([x[sub_idx], y[sub_idx]])
#kde = gaussian_kde(xy_sub, bw_method=0.3)  # Larger bw = smoother & faster

# Only evaluate KDE on a subset or grid
#z = kde(np.vstack([x, y]))

# Sort by density for better plotting
#idx = z.argsort()
#x, y, density = x[idx], y[idx], z[idx]

# -------------------------------
# Plot
# -------------------------------
fig, ax = plt.subplots(figsize=(7, 7))

# Scatter plot with density coloring
#sc = ax.scatter(x, y, c=density, cmap='plasma', s=5)
ax.scatter(
    x, y,
    s=30,                 # small marker size
    alpha=0.005,           # transparency
    color="black",         # base color
)
# y = x reference line
lims = [0, max(x.max(), y.max())]
ax.plot(lims, lims, 'k--', alpha=0.7, label='y = x')
# linear fit line
fit_x = np.linspace(lims[0], lims[1], 100)
fit_y = slope * fit_x + intercept
ax.plot(fit_x, fit_y, 'r-', label=f'York fit')

# One-sigma bounds (approx.)
sigma = np.std(y - (slope*x + intercept))
ax.fill_between(fit_x, fit_y - sigma, fit_y + sigma, color='r', alpha=0.15, label='±1σ')

# Colorbar
#cbar = plt.colorbar(sc, ax=ax)
#cbar.set_label('Density')

# Labels
ax.set_xlabel("TROPOMI " + r'$[×10^{15} \ \mathrm{molec.cm^{-2}}]$')
ax.set_ylabel("WRF-CMAQ " + r'$[×10^{15} \ \mathrm{molec.cm^{-2}}]$')

# Stats text box
stats_text = (
    f"RMSE = {rmse:.2f}\n"
    f"Mean Bias = {mean_bias:.2f}\n"
    f"Mean Abs Bias = {mean_abs_bias:.2f}\n"
    f"Index of Agreement = {index_of_agreement:.2f}\n"
    f"$R^2$ = {r2:.2f}\n"
    f"y = {slope:.2f}x + {intercept:.2f}"
)
ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
        fontsize=12, verticalalignment='top',
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray'))

ax.legend(loc='lower right')
# -------------------------------
# Save as PNG with 300 dpi
# -------------------------------
plt.tight_layout()
plt.axis('equal')
plt.savefig("scatterplot.png", dpi=300)
plt.show()
