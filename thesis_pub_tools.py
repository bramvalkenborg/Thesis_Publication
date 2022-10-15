"""
------------------------------------------------------------------------------------------------------------------------
Base script containing important functions used in other functions
data_analysis_dimension     :   Loading the dimenstions (lon, lat & time) of a certain variable
data_analysis_variable      :   Loading the data of a certain variable
using_mpl_scatter_density   :   Creating a density plot of variables x and y
density_plots               :   Shortcut for creating a density plot
scatterplots                :   Shortcut for creating a scatterplot
scatter_axis                :   Shortcut for creating a scatterplot in a subplot
timeseries_2y_latlon        :   Shortcut for creating a timeseries with 2 y-axes
min_window                  :   Calculates the 4-week window around the minimum containing the minimum values and the
                                corresponding x and y values
max_window_wet              :   Calculates the 4-week window around the maximum containing the maximum values and the
                                corresponding x and y values
------------------------------------------------------------------------------------------------------------------------
"""

# Importing all the packages
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from SIF_tools.python.L2_tools import convert_time
import numpy as np
import mpl_scatter_density
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.axis as axis
from mpl_toolkits.basemap import Basemap
from scipy import stats
from scipy import optimize
import matplotlib as mpl
from pytesmo.time_series.anomaly import calc_climatology
from pytesmo.time_series.anomaly import calc_anomaly as calc_anom_pytesmo
import pandas as pd
import datetime
from validation_good_practice.ancillary.metrics import correct_n
from patsy import dmatrices
import statsmodels_adapted.api as sm

#-----------------------------------------------------------------------------------------------------------------------
# Important variables that are often used
path_gridding = '/data/leuven/336/vsc33653/OUTPUT/Gridding/'
# path_save = '/data/leuven/336/vsc33653/OUTPUT/'
save_graph = '/data/leuven/336/vsc33653/OUTPUT_pub/Graphs/'
save_maps = '/data/leuven/336/vsc33653/OUTPUT_pub/Maps/'

"""
Import data for data analysis: Ob South
SIF_file = 'Ob_south_SIF_8day_0.05deg_c30.nc'
SIFnorm_MERRA2_file = 'Ob_south_SIFnorm_MERRA2__omitnan_8day_0.05deg_c30.nc'    # Theoretical PAR
SIFnorm_CERES_file = 'Ob_south_SIFnorm_CERES_omitnan_8day_0.05deg_c30.nc'  # Observational PAR
fPAR_file = 'Ob_South_fPAR_8day_0.05deg.nc'
peat_file = 'Siberia_PeatMap_0.05deg.nc'
WTD_file = 'Ob_South_WTD_8day_0.05deg.nc'
treecover_file = 'Ob_South_TreeCover_0.05deg.nc'
"""

# Import data for data analysis: NH
SIF_file5 = 'NH_SIF_dN_omitnan_8day_0.05deg_c30.nc'
PAR_file5 = 'NH_PAR_8day_0.05deg.nc'
SIFnorm_MERRA2_file5 = 'NH_SIFnorm_MERRA2_omitnan_8day_0.05deg_c30.nc'    # Theoretical PAR
SIFnorm_CERES_file5 = 'NH_SIFnorm_CERES_omitnan_8day_0.05deg_c30.nc'  # Observational PAR
fPAR_file5 = 'NH_fPAR_8day_0.05deg.nc'
peat_file5 = 'NH_PeatMap_0.05deg.nc'
WTD_file5 = 'NH_WTD_8day_0.05deg.nc'
treecover_file5 = 'NH_TreeFraction_0.05deg.nc'
WTD10y_file5 = 'NH_WTD_10y_8day_0.5deg.nc'

SIF_file2 = 'NH_SIF_dN_omitnan_8day_0.2deg_c30.nc' # Contains SIFnorm (PAR based on CERES)
PAR_file2 = 'NH_PAR_8day_0.2deg.nc'
peat_file2 = 'NH_PeatMap_0.2deg.nc'
WTD_file2 = 'NH_WTD_8day_0.2deg.nc'
WTD10y_file2 = 'NH_WTD_10y_8day_0.2deg.nc'

SIF2021_file2 = 'NH_SIF_2021_dN_omitnan_8day_0.2deg_c30.nc'
WTD2021_file2 = 'NH_WTD_2021_8day_0.2deg.nc'
OW2021_file2 ='AMSRU_Mland_A_v03_8day_02.nc'
PAR2021_file2 = 'NH_PAR_2021_8day_0.2deg.nc'

SIF2021_file5 = 'NH_SIF_2021_dN_omitnan_8day_0.05deg_c30.nc'
WTD2021_file5 = 'NH_WTD_2021_8day_0.05deg.nc'
OW2021_file5 = 'AMSRU_Mland_A_v03_8day_005.nc'

SIF_unit = ' (Wm$^{-2}$sr$^{-1}$µm$^{-1}$)'
fPAR_unit = ' [-]'
WTD_unit = ' (m)'
SIFy_unit = ' (sr$^{-1}$mm$^{-1}$)'
PAR_unit = ' (Wm$^{-2}$)'
APAR_unit = ' (Wm$^{-2}$)'
slope_unit = ' (mm$^{-2}$sr$^{-1}$)'
i_unit = ' (m$^{-2}$sr$^{-1}$)'

#-----------------------------------------------------------------------------------------------------------------------
# Loading dimensions from the data (lon, lat and time)
def data_analysis_dimension(file, early_summer=True, late_summer=True, lat_lim=[], lon_lim=[]):
    nc = Dataset(path_gridding + file)
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    time = nc.variables['time'][:]
    time = convert_time(time, 'days since 1970-01-01', u"gregorian")
    month = time.astype('datetime64[M]').astype(int) % 12 + 1
    if lon_lim != [] and lat_lim != []:
        cond_lat = np.asarray((lat <= lat_lim[1]) & (lat >= lat_lim[0]))
        cond_lon = np.asarray((lon <= lon_lim[1]) & (lon >= lon_lim[0]))
        i_latmax = np.max(np.where(cond_lat))
        i_latmin = np.min(np.where(cond_lat))-1
        i_lonmax = np.max(np.where(cond_lon))
        i_lonmin = np.min(np.where(cond_lon))-1
        lon = lon[i_lonmin:i_lonmax]
        lat = lat[i_latmin:i_latmax]
    # if lonMin != '' and lonMax != '' and latMin != '' and latMax != '':
        # i_lonmin = (np.abs(lon - lonMin)).argmin()
        # i_lonmax = (np.abs(lon - lonMax)).argmin()
        # i_latmin = (np.abs(lat - latMin)).argmin()
        # i_latmax = (np.abs(lat - latMax)).argmin()
        # lat = lat[i_latmin:i_latmax]
        # lon = lon[i_lonmin:i_lonmax]
    # if start != '' and end != '':
        # i_start = (np.abs(time - start)).argmin()
        # i_end = (np.abs(time - end)).argmin()
        # time = time[i_start:i_end]
        # month = time.astype('datetime64[M]').astype(int) % 12 + 1
    if early_summer and late_summer:
        data = np.where(month == 6, time,
                        np.where(month == 7, time,
                                 np.where(month == 8, time,
                                          np.where(month == 9, time, np.nan))))
    elif early_summer:
        time = np.where(month == 6, time, np.where(month == 7, time, np.nan))
    elif late_summer:
        time = np.where(month == 8, time, np.where(month == 9, time, np.nan))
    return lon, lat, time

#-----------------------------------------------------------------------------------------------------------------------
# Loading variables from a certain netCDF file
def data_analysis_variable(file, variable, peat_min=0.5, early_summer=True, late_summer=True, lon_lim=[], lat_lim=[], SP_res=0.05, Time_2021=False):
    if Time_2021 and SP_res == 0.05:
        SIF_file = SIF2021_file5
        PAR_file = PAR_file5
        peat_file = peat_file5
        if file == WTD10y_file5:
            WTD_file = WTD10y_file5
        else:
            WTD_file = WTD2021_file5
        treecover_file = treecover_file5
    elif Time_2021 and SP_res == 0.2:
        SIF_file = SIF2021_file2
        peat_file = peat_file2
        if file == WTD10y_file2:
            WTD_file = WTD10y_file2
        else:
            WTD_file = WTD2021_file2
    elif not Time_2021 and SP_res == 0.05:
        SIF_file = SIF_file5
        PAR_file = PAR_file5
        peat_file = peat_file5
        if file == WTD10y_file5:
            WTD_file = WTD10y_file5
        else:
            WTD_file = WTD_file5
        treecover_file = treecover_file5
    else:
        SIF_file = SIF_file2
        peat_file = peat_file2
        if file == WTD10y_file2:
            WTD_file = WTD10y_file2
        else:
            WTD_file = WTD_file2
    nc = Dataset(path_gridding + file)
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    time = nc.variables['time'][:]
    nc_peat = Dataset(path_gridding + peat_file)
    nc_WTD = Dataset(path_gridding + WTD_file)
    if lon_lim != [] and lat_lim != []:
        cond_lat = np.asarray((lat <= lat_lim[1]) & (lat >= lat_lim[0]))
        cond_lon = np.asarray((lon <= lon_lim[1]) & (lon >= lon_lim[0]))
        i_latmax = np.max(np.where(cond_lat))
        i_latmin = np.min(np.where(cond_lat))-1
        i_lonmax = np.max(np.where(cond_lon))
        i_lonmin = np.min(np.where(cond_lon))-1
        peat = nc_peat.variables['peatmap'][i_latmin:i_latmax, i_lonmin:i_lonmax, :]
        peat = np.tile(peat, (1, 1, len(time)))
        tp1 = nc_WTD.variables['tp1'][i_latmin:i_latmax, i_lonmin:i_lonmax, 0:len(time)]
        data = nc.variables[variable][i_latmin:i_latmax, i_lonmin:i_lonmax, 0:len(time)]
        if (variable == 'sif_dc') | (variable == 'sif_norm'):
            n = nc.variables['n'][i_latmin:i_latmax, i_lonmin:i_lonmax, 0:len(time)]
        if variable == 'sif_norm':
            data = data*1000
    else:
        peat = nc_peat.variables['peatmap'][:, :, :]
        peat = np.tile(peat, (1, 1, len(time)))
        tp1 = nc_WTD.variables['tp1'][:, :, 0:len(time)]
        time = nc.variables['time'][0:len(time)]
        data = nc.variables[variable][:, :, 0:len(time)]
        if (variable == 'sif_dc') | (variable == 'sif_norm'):
            n = nc.variables['n'][:, :, 0:len(time)]
        if variable == 'sif_norm':
            data = data*1000
    time = convert_time(time, 'days since 1970-01-01', u"gregorian")
    month = time.astype('datetime64[M]').astype(int) % 12 + 1
    if early_summer and late_summer:
        data = np.where(month == 6, data,
                        np.where(month == 7, data,
                                 np.where(month == 8, data,
                                          np.where(month == 9, data, np.nan))))
    elif early_summer:
        data = np.where(month == 6, data, np.where(month == 7, data, np.nan))
    elif late_summer:
        data = np.where(month == 8, data, np.where(month == 9, data, np.nan))
    # Mask out wrong fPAR data: 29/09/2018 (38) for lat 60-70 and lon 30-65
    if variable == 'Fpar_500m':
        data[300:499, 4000:4700, 75] = np.nan
    # Selecting only the peatland cells based on peatmap
    data = np.where(peat >= peat_min, data, np.nan)
    # Correction for snow cover
    data = np.where(tp1 < 277.15, np.nan, data)
    # Change nan from Julia (-999) to nan of Python
    data = np.where(data == -999, np.nan, data)
    # Remove all the data where there is no WTD value
    data = np.where(np.isnan(tp1), np.nan, data)
    # Remove pixels with not sufficient SIF measurements
    if (variable == 'sif_dc') | (variable == 'sif_norm') :
        data = np.where(n > 1, data, np.nan)
    # Remove area in kazachstan:
    if lat_lim==[] and lon_lim==[]:
        data[10:35, 1075:1125, :] = np.nan
    # Negative SIF values
    # if variable == 'sif_dc':
    #     year = time.astype('datetime64[Y]').astype(int) + 1970
        # Select the data for the year 2018, 2019 or 2020
    #     dataJJA18 = np.where((year == 2018) & (month != 9), data, np.nan)
    #     dataJJA19 = np.where((year == 2019) & (month != 9), data, np.nan)
    #     dataJJA20 = np.where((year == 2020) & (month != 9), data, np.nan)
    #     dataS18 = np.where((year == 2018) & (month == 9), data, np.nan)
    #     dataS19 = np.where((year == 2019) & (month == 9), data, np.nan)
    #     dataS20 = np.where((year == 2020) & (month == 9), data, np.nan)
        # For each pixel find the index of the maximum value for each year
        # data18 = np.where(np.isnan(data18), -999, data18)
        # iMax18 = (np.argsort(data18, axis=2)[:, :, ::-1])[:, :, 0]
        # iMax18 = np.where(iMax18 == 0, np.nan, iMax18)
        # data19 = np.where(np.isnan(data19), -999, data19)
        # iMax19 = (np.argsort(data19, axis=2)[:, :, ::-1])[:, :, 0]
        # iMax19 = np.where(iMax19 == 0, np.nan, iMax19)
        # data20 = np.where(np.isnan(data20), -999, data20)
        # iMax20 = (np.argsort(data20, axis=2)[:, :, ::-1])[:, :, 0]
        # iMax20 = np.where(iMax20 == 0, np.nan, iMax20)
        # All the data that is after than the maximum is equaled to NaN
        # dataMax18 = data18
        # dataMax19 = data19
        # dataMax20 = data20
        # for x in range(iMax18.shape[0]):
            # for y in range(iMax18.shape[1]):
                # if not np.isnan(iMax18[x, y]):
                    # dataMax18[x, y, int(iMax18[x, y]):-1] = np.nan
                    # dataMax19[x, y, int(iMax19[x, y]):-1] = np.nan
                    # dataMax20[x, y, int(iMax20[x, y]):-1] = np.nan
        # dataMax18 = np.where(dataMax18 == -999, np.nan, dataMax18)
        # dataMax19 = np.where(dataMax19 == -999, np.nan, dataMax19)
        # dataMax20 = np.where(dataMax20 == -999, np.nan, dataMax20)
        # data1Max18 = np.maximum.accumulate((dataMax18 < 0)[:, :, ::-1], axis=2)[:, :, ::-1]
        # data1Max19 = np.maximum.accumulate((dataMax19 < 0)[:, :, ::-1], axis=2)[:, :, ::-1]
        # data1Max20 = np.maximum.accumulate((dataMax20 < 0)[:, :, ::-1], axis=2)[:, :, ::-1]
        # data18 = np.where(data18 == -999, np.nan, data18)
        # data19 = np.where(data19 == -999, np.nan, data19)
        # data20 = np.where(data20 == -999, np.nan, data20)
    #     dataJJA18 = np.maximum.accumulate((dataJJA18 < 0)[:, :, ::-1], axis=2)[:, :, ::-1]
    #     dataJJA19 = np.maximum.accumulate((dataJJA19 < 0)[:, :, ::-1], axis=2)[:, :, ::-1]
    #     dataJJA20 = np.maximum.accumulate((dataJJA20 < 0)[:, :, ::-1], axis=2)[:, :, ::-1]
    #     dataS18 = np.maximum.accumulate((dataS18 < 0)[:, :, ::-1], axis=2)[:, :, ::-1]
    #     dataS19 = np.maximum.accumulate((dataS19 < 0)[:, :, ::-1], axis=2)[:, :, ::-1]
    #     dataS20 = np.maximum.accumulate((dataS20 < 0)[:, :, ::-1], axis=2)[:, :, ::-1]
    #     dataJJA18 = np.where((year == 2017), False, dataJJA18)
    #     dataJJA19 = np.where((year==2017) | (year==2018), False, dataJJA19)
    #     dataJJA20 = np.where((year==2017) | (year==2018) | (year==2019), False, dataJJA20)
    #     dataS18 = np.where((year == 2017) | (year == 2019) | (year == 2020), False, dataS18)
    #     dataS19 = np.where((year == 2017) | (year == 2018) | (year == 2020), False, dataS19)
    #     dataS20 = np.where((year == 2017) | (year == 2018) | (year == 2019), False, dataS20)
    #     data = np.where(dataJJA18 | dataS18, np.nan, data)
    #     data = np.where(dataJJA19 | dataS19, np.nan, data)
    #     data = np.where(dataJJA20 | dataS20, np.nan, data)
    return data

#-----------------------------------------------------------------------------------------------------------------------
def using_mpl_scatter_density(fig, x, y):
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    density = ax.scatter_density(x, y, cmap=white_viridis)
    fig.colorbar(density, label='Number of pixels')

# "Viridis-like" colormap with white background
white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#ffffff'),
    (1e-20, '#440053'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#78d151'),
    (1, '#fde624'),
], N=256)

# "Viridis-like" colormap with white background
white_viridis_r = LinearSegmentedColormap.from_list('white_viridis_r', [
    (0, '#fde624'),
    (1e-20, '#78d151'),
    (0.2, '#21a784'),
    (0.4, '#2a788e'),
    (0.6, '#404388'),
    (0.8, '#440053'),
    (1, '#ffffff'),
], N=256)

#-----------------------------------------------------------------------------------------------------------------------
def density_plots(x, y, axis=[], title='', save='', xlabel='', ylabel='', grid=False):
    x1 = x[np.logical_not(np.isnan(x))]
    y1 = y[np.logical_not(np.isnan(x))]
    x2 = x1[np.logical_not(np.isnan(y1))]
    y2 = y1[np.logical_not(np.isnan(y1))]
    x3 = x2[np.logical_not(np.isinf(x2))]
    y3 = y2[np.logical_not(np.isinf(x2))]
    x4 = x3[np.logical_not(np.isinf(y3))]
    y4 = y3[np.logical_not(np.isinf(y3))]
    x = x4
    y = y4
    fig = plt.figure()
    using_mpl_scatter_density(fig, x, y)
    if axis != []:
        plt.axis(axis)
    if title != '':
        plt.title(title)
    if xlabel != '':
        plt.xlabel(xlabel)
    if ylabel != '':
        plt.ylabel(ylabel)
    if grid:
        plt.grid()
    if save != '':
        plt.savefig(save_graph+'DensityPlots/' + save)

#-----------------------------------------------------------------------------------------------------------------------
def rm_nan(x):
    x1 = x[np.logical_not(np.isnan(x))]
    x2 = x1[np.logical_not(np.isinf(x1))]
    x = x2
    return x

def rm_nan_2D(x, y):
    x1 = x[np.logical_not(np.isnan(x))]
    y1 = y[np.logical_not(np.isnan(x))]
    x2 = x1[np.logical_not(np.isnan(y1))]
    y2 = y1[np.logical_not(np.isnan(y1))]
    x3 = x2[np.logical_not(np.isinf(x2))]
    y3 = y2[np.logical_not(np.isinf(x2))]
    x4 = x3[np.logical_not(np.isinf(y3))]
    y4 = y3[np.logical_not(np.isinf(y3))]
    x = x4
    y = y4
    return x, y

# ----------------------------------------------------------------------------------------------------------------------
# Change detection method
# ----------------------------------------------------------------------------------------------------------------------
def change_detection(lat, lon, time, SIFnorm, WTD, clim_corr=True):
    i = np.zeros((len(lat), len(lon), len(time)))
    dSIFnorm = np.zeros((len(lat), len(lon), len(time)))
    dWTDz = np.zeros((len(lat), len(lon), len(time)))
    WTDz_avg = np.zeros((len(lat), len(lon), len(time)))  # Contains the average WTD between t1 and t2
    if clim_corr:
        SIF_clim = np.zeros((len(lat), len(lon), len(time)))
        WTD_clim = np.zeros((len(lat), len(lon), len(time)))
        DOY = np.zeros(len(time))
        for time1 in range(len(time)):
            DOY[time1] = time[time1].timetuple().tm_yday
        DOY = DOY.astype(int)
    for lon1 in range(len(lon)):
        for lat1 in range(len(lat)):
            if clim_corr:
                SIF_Ser = pd.Series(SIFnorm[lat1, lon1, :], index=time)
                WTD_Ser = pd.Series(WTD[lat1, lon1, :], index=time)
                SIF_Ser = SIF_Ser.resample('1D').bfill()
                WTD_Ser = WTD_Ser.resample('1D').bfill()
                SIF_clim1 = calc_climatology(SIF_Ser)
                WTD_clim1 = calc_climatology(WTD_Ser)
                SIF_clim8d = np.zeros(len(time))
                WTD_clim8d = np.zeros(len(time))
                for DOY1, time1 in zip(DOY, range(len(time))):
                    SIF_clim8d[time1] = SIF_clim1[DOY1]
                    WTD_clim8d[time1] = WTD_clim1[DOY1]
                    SIF_clim8d[time1] = np.where(np.isnan(SIFnorm[lat1, lon1, time1]), np.nan, SIF_clim8d[time1])
                    WTD_clim8d[time1] = np.where(np.isnan(WTD[lat1, lon1, time1]), np.nan, WTD_clim8d[time1])
                    SIF_clim[lat1, lon1, time1] = SIF_clim8d[time1]
                    WTD_clim[lat1, lon1, time1] = WTD_clim8d[time1]
            for time1 in range(len(time)):
                if not time1 + 1 == len(time):
                    SIFnorm1 = SIFnorm[lat1, lon1, time1]
                    SIFnorm2 = SIFnorm[lat1, lon1, time1 + 1]
                    # WTDz1 = WTD_zScore[time1]
                    # WTDz2 = WTD_zScore[time1 + 1]
                    WTDz1 = WTD[lat1, lon1, time1]
                    WTDz2 = WTD[lat1, lon1, time1 + 1]
                    if clim_corr:
                        SIFclim1 = SIF_clim[lat1, lon1, time1]
                        SIFclim2 = SIF_clim[lat1, lon1, time1+1]
                    if not np.isnan(SIFnorm1) and not np.isnan(SIFnorm2) and not np.isnan(WTDz1) and not np.isnan(WTDz2):
                        dSIFnorm[lat1, lon1, time1] = SIFnorm2 - SIFnorm1
                        dWTDz[lat1, lon1, time1] = WTDz2 - WTDz1
                        if clim_corr:
                            dSIFclim = SIFclim2 - SIFclim1
                            dSIFnorm[lat1, lon1, time1] = dSIFnorm[lat1, lon1, time1] - dSIFclim
                            i[lat1, lon1, time1] = (dSIFnorm[lat1, lon1, time1] / dWTDz[lat1, lon1, time1])
                        else:
                            i[lat1, lon1, time1] = (dSIFnorm[lat1, lon1, time1] / dWTDz[lat1, lon1, time1])
                        WTDz_avg[lat1, lon1, time1] = np.mean((WTDz1, WTDz2))
                    else:
                        dSIFnorm[lat1, lon1, time1] = np.nan
                        dWTDz[lat1, lon1, time1] = np.nan
                        i[lat1, lon1, time1] = np.nan
                        WTDz_avg[lat1, lon1, time1] = np.nan
    i[:, :, len(time)-1] = np.nan
    # i = np.where(i > np.nanquantile(i, 0.05), i, np.nan)
    # i = np.where(i < np.nanquantile(i, 0.95), i, np.nan)
    WTDz_avg[:, :, len(time)-1] = np.nan
    if clim_corr:
        return i, dSIFnorm, dWTDz, WTDz_avg, SIF_clim, WTD_clim
    else:
        return i, dSIFnorm, dWTDz, WTDz_avg
"""
# Remove outliers of i
i = np.where(i > np.nanquantile(i, 0.01), i, np.nan)
i = np.where(i < np.nanquantile(i, 0.99), i, np.nan)
"""

# ----------------------------------------------------------------------------------------------------------------------
# Anomalies approach
# ----------------------------------------------------------------------------------------------------------------------


def calc_anom(Ser, longterm=False, window_size=35):

    anom = calc_anom_pytesmo(Ser, climatology=calc_climatology(Ser, wraparound=True) if longterm else None, window_size=window_size)
    anom.name = Ser.name

    return anom

def anomalies(lat, lon, time, SIF, WTD):
    SIF_anom = np.zeros((len(lat), len(lon), len(time)))
    WTD_anom = np.zeros((len(lat), len(lon), len(time)))
    i = np.zeros((len(lat), len(lon), len(time)))
    DOY = np.zeros(len(time))
    for time1 in range(len(time)):
        DOY[time1] = time[time1].timetuple().tm_yday
    DOY = DOY.astype(int)
    for lat1 in range(len(lat)):
        for lon1 in range(len(lon)):
            if not np.isnan(SIF[lat1, lon1, :]).all() and not np.isnan(WTD[lat1, lon1, :]).all():
                SIF_Ser = pd.Series(SIF[lat1, lon1, :], index=time)
                SIF_Ser = SIF_Ser.resample('1D').bfill()
                WTD_Ser = pd.Series(WTD[lat1, lon1, :], index=time)
                WTD_Ser = WTD_Ser.resample('1D').bfill()
                SIF_anom1 = calc_anom(SIF_Ser, longterm=True)
                WTD_anom1 = calc_anom(WTD_Ser, longterm=True)
                SIF_anom8d = np.zeros(len(time))
                WTD_anom8d = np.zeros(len(time))
                SIF_anom8d = SIF_anom1.resample('8D').first()
                WTD_anom8d = WTD_anom1.resample('8D').first()
                SIF_anom8d = np.where(np.isnan(SIF[lat1, lon1, :]), np.nan, SIF_anom8d)
                WTD_anom8d = np.where(np.isnan(SIF[lat1, lon1, :]), np.nan, WTD_anom8d)
                SIF_anom[lat1, lon1, :] = SIF_anom8d
                WTD_anom[lat1, lon1, :] = WTD_anom8d
                i[lat1, lon1, :] = SIF_anom[lat1, lon1, :]/WTD_anom[lat1, lon1, :]
            else:
                SIF_anom[lat1, lon1, :] = np.nan
                WTD_anom[lat1, lon1, :] = np.nan
                i[lat1, lon1, :] = np.nan
        print(str(lat1))
    return i, SIF_anom, WTD_anom


# ----------------------------------------------------------------------------------------------------------------------
# Longterm anomalies
# ----------------------------------------------------------------------------------------------------------------------
def long_anom(lat, lon, time, SIF, WTD):
    SIF_lAnom = np.zeros((len(lat), len(lon), len(time)))
    WTD_lAnom = np.zeros((len(lat), len(lon), len(time)))
    for lat1 in range(len(lat)):
        for lon1 in range(len(lon)):
            if not np.isnan(SIF[lat1, lon1, :]).all() and not np.isnan(WTD[lat1, lon1, :]).all():
                # Convert to panda series
                SIF_ser = pd.Series(SIF[lat1, lon1, :], index=time)
                WTD_ser = pd.Series(WTD[lat1, lon1, :], index=time)
                SIF_ser = SIF_ser.resample('1D').bfill()
                WTD_ser = WTD_ser.resample('1D').bfill()
                SIF_longAnom = calc_anom(SIF_ser, longterm=True) - calc_anom(SIF_ser, longterm=False)
                WTD_longAnom = calc_anom(WTD_ser, longterm=True) - calc_anom(WTD_ser, longterm=False)
                SIF_longAnom8d = np.zeros(len(time))
                WTD_longAnom8d = np.zeros(len(time))
                SIF_longAnom8d = SIF_longAnom.resample('8D').first()
                WTD_longAnom8d = WTD_longAnom.resample('8D').first()
                SIF_longAnom8d = np.where(np.isnan(SIF[lat1, lon1, :]), np.nan, SIF_longAnom8d)
                WTD_longAnom8d = np.where(np.isnan(SIF[lat1, lon1, :]), np.nan, WTD_longAnom8d)
                SIF_lAnom[lat1, lon1, :] = SIF_longAnom8d
                WTD_lAnom[lat1, lon1, :] = WTD_longAnom8d
            else:
                SIF_lAnom[lat1, lon1, :] = np.nan
                WTD_lAnom[lat1, lon1, :] = np.nan
        print(str(round((lat1/(len(lat))*100)))+'%')
    return SIF_lAnom, WTD_lAnom

def seasonality(lat, lon, time, x):
    seas = np.zeros((len(lat), len(lon), len(time)))
    for lat1 in range(len(lat)):
        for lon1 in range(len(lon)):
            if not np.isnan(seas[lat1, lon1, :]).all():
                # Convert to panda series
                x_ser = pd.Series(x[lat1, lon1, :], index=time)
                x_ser = x_ser.resample('1D').bfill()
                x_seas = x_ser - calc_anom(x_ser, longterm=False)
                x_seas8d = np.zeros(len(time))
                x_seas8d = x_seas.resample('8D').first()
                x_seas8d = np.where(np.isnan(x[lat1, lon1, :]), np.nan, x_seas8d)
                seas[lat1, lon1, :] = x_seas8d
            else:
                seas[lat1, lon1, :] = np.nan
        print(str(round((lat1/(len(lat))*100)))+'%')
    return seas

def lt_climatology(lat, lon, time, x):
    lt_clim = np.zeros((len(lat), len(lon), len(time)))
    for lat1 in range(len(lat)):
        for lon1 in range(len(lon)):
            if not np.isnan(x[lat1, lon1, :]).all():
                # Convert to panda series
                x_ser = pd.Series(x[lat1, lon1, :], index=time)
                x_ser = x_ser.resample('1D').bfill()
                x_clim1 = x_ser - calc_anom(x_ser, longterm=True)
                x_clim8d = np.zeros(len(time))
                x_clim8d = x_clim1.resample('8D').first()
                x_clim8d = np.where(np.isnan(x[lat1, lon1, :]), np.nan, x_clim8d)
                lt_clim[lat1, lon1, :] = x_clim8d
            else:
                lt_clim[lat1, lon1, :] = np.nan
        print(str(round((lat1/(len(lat))*100)))+'%')
    return lt_clim


# ----------------------------------------------------------------------------------------------------------------------
# Shortterm anomalies
# ----------------------------------------------------------------------------------------------------------------------
def short_anom(lat, lon, time, SIF, WTD):
    SIF_sAnom = np.zeros((len(lat), len(lon), len(time)))
    WTD_sAnom = np.zeros((len(lat), len(lon), len(time)))
    for lat1 in range(len(lat)):
        for lon1 in range(len(lon)):
            if not np.isnan(SIF[lat1, lon1, :]).all() and not np.isnan(WTD[lat1, lon1, :]).all():
                # Convert to panda series
                SIF_ser = pd.Series(SIF[lat1, lon1, :], index=time)
                WTD_ser = pd.Series(WTD[lat1, lon1, :], index=time)
                SIF_ser = SIF_ser.resample('1D').bfill()
                WTD_ser = WTD_ser.resample('1D').bfill()
                SIF_shortAnom = calc_anom(SIF_ser, longterm=False)
                WTD_shortAnom = calc_anom(WTD_ser, longterm=False)
                SIF_shortAnom8d = np.zeros(len(time))
                WTD_shortAnom8d = np.zeros(len(time))
                SIF_shortAnom8d = SIF_shortAnom.resample('8D').first()
                WTD_shortAnom8d = WTD_shortAnom.resample('8D').first()
                SIF_shortAnom8d = np.where(np.isnan(SIF[lat1, lon1, :]), np.nan, SIF_shortAnom8d)
                WTD_shortAnom8d = np.where(np.isnan(SIF[lat1, lon1, :]), np.nan, WTD_shortAnom8d)
                SIF_sAnom[lat1, lon1, :] = SIF_shortAnom8d
                WTD_sAnom[lat1, lon1, :] = WTD_shortAnom8d
            else:
                SIF_sAnom[lat1, lon1, :] = np.nan
                WTD_sAnom[lat1, lon1, :] = np.nan
        print(str(round((lat1/(len(lat))*100)))+'%')
    return SIF_sAnom, WTD_sAnom


# ----------------------------------------------------------------------------------------------------------------------
# 1D Model calibration

# One pixel calibration of the Water Stress Model
# ----------------------------------------------------------------------------------------------------------------------
def cal_WaterStressModel(SIF_anom, WTD_anom, WTD, p, time):
    x1 = WTD_anom
    x2 = WTD
    y = SIF_anom
    if not np.isnan(x1).all() and not np.isnan(x2).all() and not np.isnan(y).all():
        df = pd.DataFrame(np.transpose(np.asarray([x1, x2, y])), index=time,
                          columns=['WTD_sAnom', 'WTD', 'SIF_sAnom'])
        df.dropna(inplace=True)
        n_corr = correct_n(df)
        # n_s = len(np.where(np.logical_not(np.isnan(y)))[0])
        x1a = x1[np.logical_not(np.isnan(x1))]
        x2a = x2[np.logical_not(np.isnan(x1))]
        ya = y[np.logical_not(np.isnan(x1))]
        x1b = x1a[np.logical_not(np.isnan(x2a))]
        x2b = x2a[np.logical_not(np.isnan(x2a))]
        yb = ya[np.logical_not(np.isnan(x2a))]
        x1 = x1b[np.logical_not(np.isnan(yb))]
        x2 = x2b[np.logical_not(np.isnan(yb))]
        y = yb[np.logical_not(np.isnan(yb))]
        Y, X = dmatrices('y ~ x1 + x1 : x2')
        model = sm.OLS(y, X[:, 1:3]).fit(nobs_corr=n_corr)
        # model = sm.OLS(y, X[:, 1:3]).fit()
        # Y_r, X_r = dmatrices('y ~ x1')
        # model_r = sm.OLS(y, X_r[:, 1:3]).fit()
        # fp_values_s[lat1, lon1] = CalcP_f_statistic(model, model_r, n_corr[lat1, lon1], 2)
        fp_values = model.f_pvalue
        p_values = model.pvalues
        if fp_values <= p:
            coef = model.params
            WTD_opt = -coef[0] / coef[1]
            Rsq = model.rsquared
        else:
            coef = model.params
            WTD_opt = -coef[0] / coef[1]
            Rsq = model.rsquared
            # p_values_s[lat1, lon1, :] = np.nan
            # n_corr[lat1, lon1] = np.nan
    else:
        coef = np.nan
        fp_values = np.nan
        WTD_opt= np.nan
        n_corr = np.nan
        n = np.nan
        rho = np.nan
        Rsq = np.nan
        p_values = np.nan
    return coef, WTD_opt, n_corr, fp_values, p_values, Rsq


# ----------------------------------------------------------------------------------------------------------------------
# Bootstrap uncertainty analysis
# ----------------------------------------------------------------------------------------------------------------------
def Bootstrap_uncertainty(SIF_Anom, WTD_Anom, WTD, nRuns, p, time):
    WTDopt_bootstrap = np.zeros((nRuns))
    for i in range(nRuns):
        time_rand = np.random.choice(SIF_Anom.shape[0], SIF_Anom.size, replace=True)
        SIFanom_rand = SIF_Anom[time_rand]
        WTDanom_rand = WTD_Anom[time_rand]
        WTD_rand = WTD[time_rand]
        WTD_opt = cal_WaterStressModel(SIFanom_rand, WTDanom_rand, WTD_rand, p, time)[1]
        WTDopt_bootstrap[i] = WTD_opt
    WTDopt_mean = np.nanmean(WTDopt_bootstrap)
    WTD_CI_5 = np.nanquantile(WTD_opt, 0.05)
    WTD_CI_95 = np.nanquantile(WTD_opt, 0.95)
    return WTDopt_mean, WTD_CI_5, WTD_CI_95


