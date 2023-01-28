# Import packages
from thesis_pub_tools import *
from scipy import optimize
from scipy import stats
import sklearn.preprocessing
import sklearn.linear_model as lm
import statsmodels.api as sm
import statsmodels.formula.api as smf
import datetime
from patsy import dmatrices
# from GitClone.validation_good_practice.ancillary.metrics import correct_n

# ----------------------------------------------------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------------------------------------------------
# Main settings
peat_min = 0.50
area = 'NH'
SIF_PAR = True
Time_2021 = True
SP_res = 0.2
p = 1
# ----------------------------
corr = '_c30_pt' + str(round(peat_min * 100)) + '_SPres' + str(round(SP_res*10)) + '_pValue' + str(int(round(p*100)))

# ----------------------------------------------------------------------------------------------------------------------
# Import data
# ----------------------------------------------------------------------------------------------------------------------
if area == 'Alaska':
    lat_lim = [61, 65]
    lon_lim = [-161, -142]
    int_lon = 5
    int_lat = 1
    figsize = (6, 2)
    title = 'Alaska'
    res = 'l'
elif area == 'Boreal_Plains':
    lat_lim = [57, 60]
    lon_lim = [-125, -112]
    int_lon = 5
    int_lat = 1
    figsize = (6, 2)
    title = 'Boreal Plains'
    res = 'l'
elif area == 'Hudson_Bay':
    lat_lim = [51, 55]
    lon_lim = [-97, -80]
    int_lon = 5
    int_lat = 1
    figsize = (6, 2)
    title = 'Hudson Bay'
    res = 'l'
elif area == 'Obi_River':
    lat_lim = [58, 62]
    lon_lim = [67, 78]
    int_lon = 2
    int_lat = 1
    figsize = (6, 3)
    title = 'Obi River'
    res = 'l'
elif area == 'NH':
    lat_lim = []
    lon_lim = []
    int_lon = 20
    int_lat = 5
    figsize = (10, 1.65)
    title = ''
    res = 'c'
if SIF_PAR:
    SIF_variable = 'sif_norm'
    corr = corr + '_SIFnorm'
else:
    SIF_variable = 'sif_dc'
    corr = corr + '_SIF'
if Time_2021 and SP_res == 0.05:
    lon, lat, time = data_analysis_dimension(SIF2021_file5, lat_lim=lat_lim, lon_lim=lon_lim)
    SIFnorm = data_analysis_variable(SIF2021_file5, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD = data_analysis_variable(WTD2021_file5, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    i_unit = ' (m$^{-3}$sr$^{-1}\mu$m$^{-1}$)'
    corr = corr + '_2021'
    nYears = 4
elif Time_2021 and SP_res == 0.2:
    lon, lat, time = data_analysis_dimension(SIF2021_file2, lat_lim=lat_lim, lon_lim=lon_lim)
    SIFnorm = data_analysis_variable(SIF2021_file2, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD = data_analysis_variable(WTD2021_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD_10y = data_analysis_variable(WTD10y_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    i_unit = ' (m$^{-3}$sr$^{-1}\mu$m$^{-1}$)'
    corr = corr + '_2021'
    nYears = 4
elif not Time_2021 and SP_res == 0.05:
    lon, lat, time = data_analysis_dimension(SIF_file5, lat_lim=lat_lim, lon_lim=lon_lim)
    SIFnorm = data_analysis_variable(SIF_file5, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD = data_analysis_variable(WTD_file5, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    PAR = data_analysis_variable(PAR_file5, 'PAR_avg', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    nYears = 3
else:
    lon, lat, time = data_analysis_dimension(SIF_file2, lat_lim=lat_lim, lon_lim=lon_lim)
    SIFnorm = data_analysis_variable(SIF_file2, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD = data_analysis_variable(WTD_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    PAR = data_analysis_variable(PAR_file2, 'PAR_avg', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD_10y = data_analysis_variable(WTD10y_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    nYears = 3


# ----------------------------------------------------------------------------------------------------------------------
# Short term anomalies
# ----------------------------------------------------------------------------------------------------------------------
SIF_sAnom, WTD_sAnom = short_anom(lat, lon, time, SIFnorm, WTD)
# SIF_seas = seasonality(lat, lon, time, SIFnorm)
# WTD_seas = seasonality(lat, lon, time, WTD)
# SIF_clim = lt_climatology(lat, lon, time, SIFnorm)
# WTD_clim = lt_climatology(lat, lon, time, WTD)


def func(WTD, WTD_anom, x1, x2):
    SIFnorm_anom = x1*WTD_anom+x2*WTD_anom*WTD
    return SIFnorm_anom


coef_s = np.zeros((len(lat), len(lon), 2))
p_values_s = np.zeros((len(lat), len(lon), 2))
fp_values_s = np.zeros((len(lat), len(lon)))
days_wet_s = np.zeros((len(lat), len(lon)))
days_dry_s = np.zeros((len(lat), len(lon)))
n_corr = np.zeros((len(lat), len(lon)))
slope_s = np.zeros((len(lat), len(lon), len(time)))
for lat1 in range(len(lat)):
    print('Calibration: '+str(round((lat1/len(lat))*100, 2))+'%')
    for lon1 in range(len(lon)):
        x1 = WTD_sAnom[lat1, lon1, :]
        x2 = WTD[lat1, lon1, :]
        y = SIF_sAnom[lat1, lon1, :]
        if not np.isnan(x1).all() and not np.isnan(x2).all() and not np.isnan(y).all():
            # df = pd.DataFrame(np.transpose(np.asarray([x1, x2, y])), index=time, columns=['WTD_sAnom', 'WTD', 'SIF_sAnom'])
            # n_corr[lat1, lon1] = correct_n(df)
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
            model = sm.OLS(y, X[:, 1:3]).fit()
            fp_values_s[lat1, lon1] = model.f_pvalue
            p_values_s[lat1, lon1, :] = model.pvalues
            if p_values_s[lat1, lon1, 1] <= p:
                coef_s[lat1, lon1, :] = model.params
                for time1 in range(len(time)):
                    slope_s[lat1, lon1, time1] = coef_s[lat1, lon1, 1] * WTD[lat1, lon1, time1] + coef_s[lat1, lon1, 0]
                if not np.isnan(SIFnorm[lat1, lon1, :]).all() and coef_s[lat1, lon1, 1] < 0:
                    days_wet_s[lat1, lon1] = (len(np.where(slope_s[lat1, lon1, :] < 0)[0]) / nYears) * 8
                    days_dry_s[lat1, lon1] = (len(np.where(slope_s[lat1, lon1, :] > 0)[0]) / nYears) * 8
                elif not np.isnan(SIFnorm[lat1, lon1, :]).all() and coef_s[lat1, lon1, 1] > 0:
                    days_wet_s[lat1, lon1] = -10
                    days_dry_s[lat1, lon1] = -10
                else:
                    days_wet_s[lat1, lon1] = np.nan
                    days_dry_s[lat1, lon1] = np.nan
            else:
                coef_s[lat1, lon1, :] = np.nan
                p_values_s[lat1, lon1, :] = np.nan
                days_wet_s[lat1, lon1] = np.nan
                days_dry_s[lat1, lon1] = np.nan
        else:
            coef_s[lat1, lon1, :] = np.nan
            p_values_s[lat1, lon1, :] = np.nan
            days_wet_s[lat1, lon1] = np.nan
            days_dry_s[lat1, lon1] = np.nan


# ----------------------------------------------------------------------------------------------------------------------
# Long term anomalies
# ----------------------------------------------------------------------------------------------------------------------
SIF_lAnom, WTD_lAnom = long_anom(lat, lon, time, SIFnorm, WTD)
# SIF_seas = seasonality(lat, lon, time, SIFnorm)
# WTD_seas = seasonality(lat, lon, time, WTD)
# SIF_clim = lt_climatology(lat, lon, time, SIFnorm)
# WTD_clim = lt_climatology(lat, lon, time, WTD)


def func(WTD, WTD_anom, x1, x2):
    x1 = np.repeat(x1[:, :, np.newaxis], len(WTD[0, 0, :]), axis=2)
    x2 = np.repeat(x2[:, :, np.newaxis], len(WTD[0, 0, :]), axis=2)
    SIFnorm_anom = x1*WTD_anom+x2*WTD_anom*WTD
    return SIFnorm_anom


coef_l = np.zeros((len(lat), len(lon), 2))
p_values_l = np.zeros((len(lat), len(lon), 2))
fp_values_l = np.zeros((len(lat), len(lon)))
days_wet_l = np.zeros((len(lat), len(lon)))
days_dry_l = np.zeros((len(lat), len(lon)))
slope_l = np.zeros((len(lat), len(lon), len(time)))
for lat1 in range(len(lat)):
    print('Calibration: '+str(round((lat1 / len(lat)) * 100, 2)) + '%')
    for lon1 in range(len(lon)):
        x1 = WTD_lAnom[lat1, lon1, :]
        x2 = WTD[lat1, lon1, :]
        y = SIF_lAnom[lat1, lon1, :]
        if not np.isnan(x1).all() and not np.isnan(x2).all() and not np.isnan(y).all():
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
            model = sm.OLS(y, X[:, 1:3]).fit()
            fp_values_l[lat1, lon1] = model.f_pvalue
            p_values_l[lat1, lon1, :] = model.pvalues
            if p_values_l[lat1, lon1, 1] <= p:
                coef_l[lat1, lon1, :] = model.params
                for time1 in range(len(time)):
                    slope_l[lat1, lon1, time1] = coef_l[lat1, lon1, 1] * WTD[lat1, lon1, time1] + coef_l[lat1, lon1, 0]
                if not np.isnan(SIFnorm[lat1, lon1, :]).all() and coef_l[lat1, lon1, 1] < 0:
                    days_wet_l[lat1, lon1] = (len(np.where(slope_l[lat1, lon1, :] < 0)[0]) / nYears) * 8
                    days_dry_l[lat1, lon1] = (len(np.where(slope_l[lat1, lon1, :] > 0)[0]) / nYears) * 8
                elif not np.isnan(SIFnorm[lat1, lon1, :]).all() and coef_l[lat1, lon1, 1] > 0:
                    days_dry_l[lat1, lon1] = -10
                    days_wet_l[lat1, lon1] = -10
                else:
                    days_dry_l[lat1, lon1] = np.nan
                    days_wet_l[lat1, lon1] = np.nan
            else:
                coef_l[lat1, lon1, :] = np.nan
                p_values_l[lat1, lon1, :] = np.nan
                days_dry_l[lat1, lon1] = np.nan
                days_wet_l[lat1, lon1] = np.nan
        else:
            coef_l[lat1, lon1, :] = np.nan
            p_values_l[lat1, lon1, :] = np.nan
            days_dry_l[lat1, lon1] = np.nan
            days_wet_l[lat1, lon1] = np.nan

# ----------------------------------------------------------------------------------------------------------------------
# Maps
# ----------------------------------------------------------------------------------------------------------------------
days_dry_sW = days_dry_s[:, 0:550]
days_dry_sE = days_dry_s[:, 774:]
days_dry_lW = days_dry_l[:, 0:550]
days_dry_lE = days_dry_l[:, 774:]

lonsW, latsW = np.meshgrid(lon[0:550], lat)
lonsE, latsE = np.meshgrid(lon[774:], lat)
cmap = 'plasma'
cmap = plt.cm.get_cmap('plasma')
cmap.set_under('darkgrey')
fig = plt.figure(figsize=(9.5, 4.6))
gs = fig.add_gridspec(2, 2, wspace=0.1)
plt_img1 = np.ma.masked_invalid(np.ma.squeeze(days_dry_sW))
plt_img2 = np.ma.masked_invalid(np.ma.squeeze(days_dry_sE))
plt_img3 = np.ma.masked_invalid(np.ma.squeeze(days_dry_lW))
plt_img4 = np.ma.masked_invalid(np.ma.squeeze(days_dry_lE))
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title('(a)', loc='left')
m1 = Basemap(projection='mill', llcrnrlat=latsW.min(), llcrnrlon=lonsW.min(), urcrnrlat=latsW.max(), urcrnrlon=lonsW.max(),
             resolution=res)
m1.drawcountries()
m1.drawcoastlines()
Basemap.drawlsmask(m1, land_color='whitesmoke', ocean_color="aliceblue")
m1.drawmeridians(np.arange(round(lonsW.min()), round(lonsW.max()), round(int_lon)))
m1.drawparallels(np.arange(round(latsW.min()), round(latsW.max()), round(int_lat)), labels=[1, 0, 0, 1])
im = m1.pcolormesh(lonsW, latsW, plt_img1, cmap=cmap, latlon=True)
im.set_clim(0, 122)
# cb = m1.colorbar(label='SIF\n'+SIF_unit, shrink=0.50)

ax2 = fig.add_subplot(gs[0, 1])
ax2.set_title('(a)', loc='left')
# ax2.set_title('x$_{1, avg}$: '+str(round(np.nanmean(coef_s[:, :, 0]), 2))+', x$_{2, avg}$: '+str(round(np.nanmean(coef_s[:, :, 1]), 2)), loc='right')
m2 = Basemap(projection='mill', llcrnrlat=latsE.min(), llcrnrlon=lonsE.min(), urcrnrlat=latsE.max(), urcrnrlon=lonsE.max(),
             resolution=res)
m2.drawcountries()
m2.drawcoastlines()
Basemap.drawlsmask(m2, land_color='whitesmoke', ocean_color="aliceblue")
m2.drawmeridians(np.arange(round(lonsE.min()), round(lonsE.max()), round(int_lon)))
m2.drawparallels(np.arange(round(latsE.min()), round(latsE.max()), round(int_lat)))
im = m2.pcolormesh(lonsE, latsE, plt_img2, cmap=cmap, latlon=True)
im.set_clim(0, 122)

ax3 = fig.add_subplot(gs[1, 0])
ax3.set_title('(b)', loc='left')
m3 = Basemap(projection='mill', llcrnrlat=latsW.min(), llcrnrlon=lonsW.min(), urcrnrlat=latsW.max(), urcrnrlon=lonsW.max(),
             resolution=res)
m3.drawcountries()
m3.drawcoastlines()
Basemap.drawlsmask(m3, land_color='whitesmoke', ocean_color="aliceblue")
m3.drawmeridians(np.arange(round(lonsW.min()), round(lonsW.max()), round(int_lon)), labels=[1, 0, 0, 1])
m3.drawparallels(np.arange(round(latsW.min()), round(latsW.max()), round(int_lat)), labels=[1, 0, 0, 1])
im = m3.pcolormesh(lonsW, latsW, plt_img3, cmap=cmap, latlon=True)
im.set_clim(0, 122)
# cb = m2.colorbar(label='WTD'+WTD_unit)

ax4 = fig.add_subplot(gs[1, 1])
ax4.set_title('(b)', loc='left')
# ax4.set_title('x$_{1, avg}$: '+str(round(np.nanmean(coef_l[:, :, 0]), 2))+', x$_{2, avg}$: '+str(round(np.nanmean(coef_l[:, :, 1]), 2)), loc='right')
m4 = Basemap(projection='mill', llcrnrlat=latsE.min(), llcrnrlon=lonsE.min(), urcrnrlat=latsE.max(), urcrnrlon=lonsE.max(),
             resolution=res)
m4.drawcountries()
m4.drawcoastlines()
Basemap.drawlsmask(m4, land_color='whitesmoke', ocean_color="aliceblue")
m4.drawmeridians(np.arange(round(lonsE.min()), round(lonsE.max()), round(int_lon)), labels=[1, 0, 0, 1])
m4.drawparallels(np.arange(round(latsE.min()), round(latsE.max()), round(int_lat)))
im = m4.pcolormesh(lonsE, latsE, plt_img4, cmap=cmap, latlon=True)
im.set_clim(0, 122)

fig.tight_layout()
cb_dry = plt.colorbar(ax=(ax1, ax2, ax3, ax4), label='% of drought regime days', location='bottom', shrink=0.2, pad=0.1, extend='min')
# cb_wet = plt.colorbar(ax=(ax1, ax2, ax3, ax4), label='% of waterlogging regime days', location='bottom', shrink=0.2, pad=0.1, extend='min')
cb_dry.set_ticks([-10, 0, 30, 61, 91, 122])
cb_dry.set_ticklabels(['x2>0', '0', '25', '50', '75', '100'])
# cb_wet.set_ticks([-10, 0, 30, 61, 91, 122])
# cb_wet.set_ticklabels(['', '', '', '', '', ''])

plt.savefig(save_maps+"GRL/maps_NumberDays2021_pix"+corr)
