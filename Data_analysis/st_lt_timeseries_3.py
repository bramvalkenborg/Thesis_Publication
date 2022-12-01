# Import packages
from thesis_pub_tools import *
from scipy import optimize
from scipy import stats
import sklearn.preprocessing
import sklearn.linear_model as lm
import statsmodels.api as sm
import datetime
from patsy import dmatrices
import matplotlib.colors as clr
from matplotlib.pyplot import cm

# ----------------------------------------------------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------------------------------------------------
# Main settings
peat_min = 0.50
area = 'NH'
SIF_PAR = True
Time_2021 = True
SP_res = 0.2
p = 0.05
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
    # lat_lim = [51, 53]
    # lon_lim = [-86, -84]
    lat_lim = [63, 64]
    lon_lim = [-117, -115]
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
    OW = data_analysis_variable(OW2021_file5, 'OW', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    i_unit = ' (m$^{-3}$sr$^{-1}\mu$m$^{-1}$)'
    corr = corr + '_2021'
    nYears = 4
elif Time_2021 and SP_res == 0.2:
    lon, lat, time = data_analysis_dimension(SIF2021_file2, lat_lim=lat_lim, lon_lim=lon_lim)
    SIFnorm = data_analysis_variable(SIF2021_file2, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD = data_analysis_variable(WTD2021_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD_10y = data_analysis_variable(WTD10y_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    OW = data_analysis_variable(OW2021_file2, 'OW', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res,
                                Time_2021=Time_2021)
    i_unit = ' (m$^{-3}$sr$^{-1}\mu$m$^{-1}$)'
    corr = corr + '_2021'
    nYears = 4
elif not Time_2021 and SP_res == 0.05:
    lon, lat, time = data_analysis_dimension(SIF_file5, lat_lim=lat_lim, lon_lim=lon_lim)
    SIFnorm = data_analysis_variable(SIF_file5, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD = data_analysis_variable(WTD_file5, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    PAR = data_analysis_variable(PAR_file5, 'PAR_avg', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    # OW = data_analysis_variable(OW2021_file5, 'OW', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res,
    #                            Time_2021=Time_2021)
    nYears = 3
else:
    lon, lat, time = data_analysis_dimension(SIF_file2, lat_lim=lat_lim, lon_lim=lon_lim)
    SIFnorm = data_analysis_variable(SIF_file2, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD = data_analysis_variable(WTD_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    PAR = data_analysis_variable(PAR_file2, 'PAR_avg', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD_10y = data_analysis_variable(WTD10y_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    # OW = data_analysis_variable(OW2021_file5, 'OW', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res,
    #                            Time_2021=Time_2021)
    nYears = 3

SIF_sAnom, WTD_sAnom = short_anom(lat, lon, time, SIFnorm, WTD)
SIF_lAnom, WTD_lAnom = long_anom(lat, lon, time, SIFnorm, WTD)
SIF_seas = seasonality(lat, lon, time, SIFnorm)
WTD_seas = seasonality(lat, lon, time, WTD)
SIF_clim = lt_climatology(lat, lon, time, SIFnorm)
WTD_clim = lt_climatology(lat, lon, time, WTD)


def lm(lat, lon, time, WTD_Anom, WTD, SIF_Anom, nYears):
    x1 = WTD_Anom[lat, lon, :]
    x2 = WTD[lat, lon, :]
    y = SIF_Anom[lat, lon, :]
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
        fp_values = model.f_pvalue
        p_values = model.pvalues
        coef = model.params
        if fp_values <= 0.05:
            for time1 in range(len(time)):
                slope = coef[1] * WTD[:, :, time1] + coef[0]
                Rsq_adj = model.rsquared_adj
                days_wet = (len(np.where(slope < 0)[0]) / nYears) * 8
                days_dry = (len(np.where(slope > 0)[0]) / nYears) * 8
                WTD_opt = -coef[0]/coef[1]
        else:
            coef = np.array((np.nan, np.nan))
            fp_values = np.nan
            slope = np.nan
            days_wet = np.nan
            days_dry = np.nan
            Rsq_adj = np.nan
            WTD_opt = np.nan
    else:
        coef = np.array((np.nan, np.nan))
        fp_values = np.nan
        slope = np.nan
        days_wet = np.nan
        days_dry = np.nan
        Rsq_adj = np.nan
        WTD_opt = np.nan
    return coef, slope, fp_values, days_wet, days_dry, Rsq_adj, WTD_opt


def func(WTDanom, WTD, a, b):
    return a * WTDanom + b * WTDanom * WTD


interaction_st = np.zeros((len(lat), len(lon)))
coef_st = np.zeros((len(lat), len(lon), 2))
WTD_opt_st = np.zeros((len(lat), len(lon)))
Rsq_st = np.zeros((len(lat), len(lon)))
for lat1 in range(len(lat)):
    for lon1 in range(len(lon)):
        x1 = WTD_sAnom[lat1, lon1, :]
        x2 = WTD[lat1, lon1, :]
        y = SIF_sAnom[lat1, lon1, :]
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
            fp_values = model.f_pvalue
            p_values = model.pvalues
            if fp_values <= p:
                coef_st[lat1, lon1, :] = model.params
                WTD_opt_st[lat1, lon1] = -coef_st[lat1, lon1, 0] / coef_st[lat1, lon1, 1]
                interaction_st[lat1, lon1] = True
                Rsq_st[lat1, lon1] = model.rsquared
            else:
                Y, X = dmatrices('y ~ x1')
                model = sm.OLS(y, X[:, 1:3]).fit()
                fp_values = model.f_pvalue
                interaction_st[lat1, lon1] = False
                Rsq_st[lat1, lon1] = np.nan
                if fp_values <= p:
                    coef_st[lat1, lon1, 0] = model.params
                    coef_st[lat1, lon1, 1] = np.nan
                    Rsq_st[lat1, lon1] = np.nan
                else:
                    coef_st[lat1, lon1, :] = np.nan
                    Rsq_st[lat1, lon1] = np.nan
        else:
            coef_st[lat1, lon1, :] = np.nan
            interaction_st[lat1, lon1] = False
            Rsq_st[lat1, lon1] = np.nan

slope_st = np.zeros((len(lat), len(lon), len(time)))
for time1 in range(len(time)):
    slope_st[:, :, time1] = coef_st[:, :, 1]*WTD[:, :, time1]+coef_st[:, :, 0]

interaction_lt = np.zeros((len(lat), len(lon)))
coef_lt = np.zeros((len(lat), len(lon), 2))
WTD_opt_lt = np.zeros((len(lat), len(lon)))
Rsq_lt = np.zeros((len(lat), len(lon)))
for lat1 in range(len(lat)):
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
            fp_values = model.f_pvalue
            p_values = model.pvalues
            if fp_values <= p:
                coef_lt[lat1, lon1, :] = model.params
                WTD_opt_lt[lat1, lon1] = -coef_lt[lat1, lon1, 0] / coef_lt[lat1, lon1, 1]
                interaction_lt[lat1, lon1] = True
                Rsq_lt[lat1, lon1] = model.rsquared
            else:
                Y, X = dmatrices('y ~ x1')
                model = sm.OLS(y, X[:, 1:3]).fit()
                fp_values = model.f_pvalue
                interaction_lt[lat1, lon1] = False
                Rsq_lt[lat1, lon1] = np.nan
                if fp_values <= p:
                    coef_lt[lat1, lon1, 0] = model.params
                    coef_lt[lat1, lon1, 1] = np.nan
                    Rsq_lt[lat1, lon1] = np.nan
                else:
                    coef_lt[lat1, lon1, :] = np.nan
                    Rsq_lt[lat1, lon1] = np.nan
        else:
            coef_lt[lat1, lon1, :] = np.nan
            interaction_lt[lat1, lon1] = False
            Rsq_lt[lat1, lon1] = np.nan

slope_lt = np.zeros((len(lat), len(lon), len(time)))
for time1 in range(len(time)):
    slope_lt[:, :, time1] = coef_lt[:, :, 1]*WTD[:, :, time1]+coef_lt[:, :, 0]

pixel = (interaction_lt==True) & (interaction_st==True) & (WTD_opt_st > np.nanmin(WTD, axis=2)) & (WTD_opt_st < np.nanmax(WTD, axis=2)) & (WTD_opt_lt > np.nanmin(WTD, axis=2)) & (WTD_opt_lt < np.nanmax(WTD, axis=2))
pixel_true = np.where(pixel)
Rsq_st = np.where(pixel, Rsq_st, np.nan)
Rsq_lt = np.where(pixel, Rsq_lt, np.nan)
"""
lat1 = 53.3
lon1 = -86.7
ilat = (np.abs(lat-lat1)).argmin()
ilon = (np.abs(lon-lon1)).argmin()
"""

"""
slope_st = lm(ilat, ilon, time, WTD_sAnom, WTD, SIF_sAnom, nYears)[1]
slope_lt = lm(ilat, ilon, time, WTD_lAnom, WTD, SIF_lAnom, nYears)[1]
WTD_opt_st = lm(ilat, ilon, time, WTD_sAnom, WTD, SIF_sAnom, nYears)[5]
WTD_opt_lt = lm(ilat, ilon, time, WTD_lAnom, WTD, SIF_lAnom, nYears)[5]
x1_st = lm(ilat, ilon, time, WTD_sAnom, WTD, SIF_sAnom, nYears)[0][0]
x2_st = lm(ilat, ilon, time, WTD_sAnom, WTD, SIF_sAnom, nYears)[0][1]
x1_lt = lm(ilat, ilon, time, WTD_lAnom, WTD, SIF_lAnom, nYears)[0][0]
x2_lt = lm(ilat, ilon, time, WTD_lAnom, WTD, SIF_lAnom, nYears)[0][1]
for ilat in pixel_true[0][0:1]:
    for ilon in pixel_true[1][0:1]:
"""
np.where(Rsq_st == np.nanmax(Rsq_st))
np.where(Rsq_lt == np.nanmax(Rsq_lt))

for p in np.random.randint(len(pixel_true[0]), size=200):
    # ilat = pixel_true[0][p]
    # ilon = pixel_true[1][p]
    ilat = 4
    ilon = 1
    slope_st_pixel = slope_st[ilat, ilon, :]
    slope_lt_pixel = slope_lt[ilat, ilon, :]
    WTD_opt_st_pixel = WTD_opt_st[ilat, ilon]
    WTD_opt_lt_pixel = WTD_opt_lt[ilat, ilon]
    x1_st = coef_st[ilat, ilon, 0]
    x2_st = coef_st[ilat, ilon, 1]
    x1_lt = coef_lt[ilat, ilon, 0]
    x2_lt = coef_lt[ilat, ilon, 1]

    SIF_sAnom_wet = np.where(WTD[ilat, ilon, :] > WTD_opt_st[ilat, ilon], SIF_sAnom[ilat, ilon, :], np.nan)
    SIF_sAnom_dry = np.where(WTD[ilat, ilon, :] < WTD_opt_st[ilat, ilon], SIF_sAnom[ilat, ilon, :], np.nan)
    SIF_lAnom_wet = np.where(WTD[ilat, ilon, :] > WTD_opt_lt[ilat, ilon], SIF_lAnom[ilat, ilon, :], np.nan)
    SIF_lAnom_dry = np.where(WTD[ilat, ilon, :] < WTD_opt_lt[ilat, ilon], SIF_lAnom[ilat, ilon, :], np.nan)
    WTD_sAnom_wet = np.where(WTD[ilat, ilon, :] > WTD_opt_st[ilat, ilon], WTD_sAnom[ilat, ilon, :], np.nan)
    WTD_sAnom_dry = np.where(WTD[ilat, ilon, :] < WTD_opt_st[ilat, ilon], WTD_sAnom[ilat, ilon, :], np.nan)
    WTD_lAnom_wet = np.where(WTD[ilat, ilon, :] > WTD_opt_lt[ilat, ilon], WTD_lAnom[ilat, ilon, :], np.nan)
    WTD_lAnom_dry = np.where(WTD[ilat, ilon, :] < WTD_opt_lt[ilat, ilon], WTD_lAnom[ilat, ilon, :], np.nan)
    fig = plt.figure(figsize=(10, 5))
    gs = fig.add_gridspec(4, 2, wspace=0.1, height_ratios=[1, 1, 0.7, 0.7])
    gs1 = gs[0].subgridspec(1, 4, wspace=0.05)
    gs2 = gs[1].subgridspec(1, 4, wspace=0.05)
    gs3 = gs[2].subgridspec(1, 4, wspace=0.05)
    gs4 = gs[3].subgridspec(1, 4, wspace=0.05)
    gs5 = gs[2:, 0].subgridspec(1, 2)
    gs6 = gs[2:, 1].subgridspec(1, 2)
    ax1a = fig.add_subplot(gs1[0, 0])
    ax1b = fig.add_subplot(gs1[0, 1], sharey=ax1a)
    ax1c = fig.add_subplot(gs1[0, 2], sharey=ax1a)
    ax1d = fig.add_subplot(gs1[0, 3], sharey=ax1a)
    ax2a = fig.add_subplot(gs2[0, 0])
    ax2b = fig.add_subplot(gs2[0, 1], sharey=ax2a)
    ax2c = fig.add_subplot(gs2[0, 2], sharey=ax2a)
    ax2d = fig.add_subplot(gs2[0, 3], sharey=ax2a)
    ax3a = fig.add_subplot(gs3[0, 0])
    ax3b = fig.add_subplot(gs3[0, 1], sharey=ax3a)
    ax3c = fig.add_subplot(gs3[0, 2], sharey=ax3a)
    ax3d = fig.add_subplot(gs3[0, 3], sharey=ax3a)
    ax4a = fig.add_subplot(gs4[0, 0])
    ax4b = fig.add_subplot(gs4[0, 1], sharey=ax4a)
    ax4c = fig.add_subplot(gs4[0, 2], sharey=ax4a)
    ax4d = fig.add_subplot(gs4[0, 3], sharey=ax4a)
    ax5a = fig.add_subplot(gs5[0, 1])
    ax5b = fig.add_subplot(gs5[0, 0], sharey=ax5a)
    ax6a = fig.add_subplot(gs6[0, 1])
    ax6b = fig.add_subplot(gs6[0, 0], sharey=ax6a)
    plt.setp(ax1a.get_xticklabels(), visible=False)
    plt.setp(ax1b.get_xticklabels(), visible=False)
    plt.setp(ax1c.get_xticklabels(), visible=False)
    plt.setp(ax1d.get_xticklabels(), visible=False)
    plt.setp(ax1b.get_yticklabels(), visible=False)
    plt.setp(ax1c.get_yticklabels(), visible=False)
    plt.setp(ax1d.get_yticklabels(), visible=False)
    plt.setp(ax2a.get_xticklabels(), visible=False)
    plt.setp(ax2b.get_xticklabels(), visible=False)
    plt.setp(ax2c.get_xticklabels(), visible=False)
    plt.setp(ax2d.get_xticklabels(), visible=False)
    plt.setp(ax2a.get_yticklabels(), visible=False)
    plt.setp(ax2b.get_yticklabels(), visible=False)
    plt.setp(ax2c.get_yticklabels(), visible=False)
    plt.setp(ax2d.get_yticklabels(), visible=False)
    plt.setp(ax3b.get_yticklabels(), visible=False)
    plt.setp(ax3c.get_yticklabels(), visible=False)
    plt.setp(ax3d.get_yticklabels(), visible=False)
    plt.setp(ax4a.get_yticklabels(), visible=False)
    plt.setp(ax4b.get_yticklabels(), visible=False)
    plt.setp(ax4c.get_yticklabels(), visible=False)
    plt.setp(ax4d.get_yticklabels(), visible=False)
    plt.setp(ax5a.get_yticklabels(), visible=False)
    plt.setp(ax6a.get_yticklabels(), visible=False)
    # plt.setp(ax5.get_yticklabels(), visible=False)
    # plt.setp(ax6.get_zticklabels(), visible=False)
    # plt.setp(ax6.get_yticklabels(), visible=False)
    ax1a.grid()
    ax1b.grid()
    ax1c.grid()
    ax1d.grid()
    ax2a.grid()
    ax2b.grid()
    ax2c.grid()
    ax2d.grid()
    ax3a.grid()
    ax3b.grid()
    ax3c.grid()
    ax3d.grid()
    ax4a.grid()
    ax4b.grid()
    ax4c.grid()
    ax4d.grid()
    ax5a.grid()
    ax5b.grid()
    ax6a.grid()
    ax6b.grid()

    ax1a.set_ylabel('SIF$_{n}$\n'+SIFy_unit, fontsize=10)
    ax3a.set_ylabel('WTD (m)', fontsize=10)
    # ax4b.set_xlabel('Time')
    # fig.text(0.72, 0.43, 'Time', ha='center')
    # fig.text(0.30, 0.43, 'Time', ha='center')
    ax1a.set_title('(a)', loc='left', fontsize=10)
    # ax1b.set_title('Short term', fontsize=12)
    fig.text(0.30, 0.92, 'Short-term', ha='center', fontsize=12)
    ax2a.set_title('(d)', loc='left', fontsize=10)
    # ax2b.set_title('Long term', fontsize=12)
    fig.text(0.72, 0.92, 'Long-term', ha='center', fontsize=12)
    ax3a.set_title('(b)', loc='left', fontsize=10)
    ax4a.set_title('(e)', loc='left', fontsize=10)

    ax1a.set_xlim(datetime.datetime(2018, 6, 1), datetime.datetime(2018, 10, 1))
    ax1b.set_xlim(datetime.datetime(2019, 6, 1), datetime.datetime(2019, 10, 1))
    ax1c.set_xlim(datetime.datetime(2020, 6, 1), datetime.datetime(2020, 10, 1))
    ax1d.set_xlim(datetime.datetime(2021, 6, 1), datetime.datetime(2021, 10, 1))
    ax2a.set_xlim(datetime.datetime(2018, 6, 1), datetime.datetime(2018, 10, 1))
    ax2b.set_xlim(datetime.datetime(2019, 6, 1), datetime.datetime(2019, 10, 1))
    ax2c.set_xlim(datetime.datetime(2020, 6, 1), datetime.datetime(2020, 10, 1))
    ax2d.set_xlim(datetime.datetime(2021, 6, 1), datetime.datetime(2021, 10, 1))
    ax3a.set_xlim(datetime.datetime(2018, 6, 1), datetime.datetime(2018, 10, 1))
    ax3b.set_xlim(datetime.datetime(2019, 6, 1), datetime.datetime(2019, 10, 1))
    ax3c.set_xlim(datetime.datetime(2020, 6, 1), datetime.datetime(2020, 10, 1))
    ax3d.set_xlim(datetime.datetime(2021, 6, 1), datetime.datetime(2021, 10, 1))
    ax4a.set_xlim(datetime.datetime(2018, 6, 1), datetime.datetime(2018, 10, 1))
    ax4b.set_xlim(datetime.datetime(2019, 6, 1), datetime.datetime(2019, 10, 1))
    ax4c.set_xlim(datetime.datetime(2020, 6, 1), datetime.datetime(2020, 10, 1))
    ax4d.set_xlim(datetime.datetime(2021, 6, 1), datetime.datetime(2021, 10, 1))
    ax1a.set_ylim(-0.3, 0.6)
    ax1b.set_ylim(-0.3, 0.6)
    ax1c.set_ylim(-0.3, 0.6)
    ax1d.set_ylim(-0.3, 0.6)
    ax2a.set_ylim(-0.3, 0.6)
    ax2b.set_ylim(-0.3, 0.6)
    ax2c.set_ylim(-0.3, 0.6)
    ax2d.set_ylim(-0.3, 0.6)
    ax3a.set_ylim(-0.3, 0)
    ax3b.set_ylim(-0.3, 0)
    ax3c.set_ylim(-0.3, 0)
    ax3d.set_ylim(-0.3, 0)
    ax4a.set_ylim(-0.3, 0)
    ax4b.set_ylim(-0.3, 0)
    ax4c.set_ylim(-0.3, 0)
    ax4d.set_ylim(-0.3, 0)
    line1, = ax1a.plot(time, SIFnorm[ilat, ilon, :], color='darkorange', label='Orig')
    line2, = ax1a.plot(time, SIF_seas[ilat, ilon, :], color='gray', label='Seas')
    ax1b.plot(time, SIFnorm[ilat, ilon, :], color='darkorange')
    ax1b.plot(time, SIF_seas[ilat, ilon, :], color='gray')
    ax1c.plot(time, SIFnorm[ilat, ilon, :], color='darkorange')
    ax1c.plot(time, SIF_seas[ilat, ilon, :], color='gray')
    ax1d.plot(time, SIFnorm[ilat, ilon, :], color='darkorange')
    ax1d.plot(time, SIF_seas[ilat, ilon, :], color='gray')

    ax2a.plot(time, SIF_seas[ilat, ilon, :], color='gray')
    line3, = ax2a.plot(time, SIF_clim[ilat, ilon, :], color='orangered', label='Clim')
    ax2b.plot(time, SIF_seas[ilat, ilon, :], color='gray')
    ax2b.plot(time, SIF_clim[ilat, ilon, :], color='orangered')
    ax2c.plot(time, SIF_seas[ilat, ilon, :], color='gray')
    ax2c.plot(time, SIF_clim[ilat, ilon, :], color='orangered')
    ax2d.plot(time, SIF_seas[ilat, ilon, :], color='gray')
    ax2d.plot(time, SIF_clim[ilat, ilon, :], color='orangered')

    line4, = ax3a.plot(time, WTD[ilat, ilon, :], color='tab:blue', label='0rig')
    line5, = ax3a.plot(time, WTD_seas[ilat, ilon, :], color='gray', label='Seas')
    ax3b.plot(time, WTD[ilat, ilon, :], color='tab:blue')
    ax3b.plot(time, WTD_seas[ilat, ilon, :], color='gray')
    ax3c.plot(time, WTD[ilat, ilon, :], color='tab:blue')
    ax3c.plot(time, WTD_seas[ilat, ilon, :], color='gray')
    ax3d.plot(time, WTD[ilat, ilon, :], color='tab:blue')
    ax3d.plot(time, WTD_seas[ilat, ilon, :], color='gray')

    ax4a.plot(time, WTD_seas[ilat, ilon, :], color='gray')
    line6, = ax4a.plot(time, WTD_clim[ilat, ilon, :], color='mediumslateblue', label='Clim')
    ax4b.plot(time, WTD_seas[ilat, ilon, :], color='gray')
    ax4b.plot(time, WTD_clim[ilat, ilon, :], color='mediumslateblue')
    ax4c.plot(time, WTD_seas[ilat, ilon, :], color='gray')
    ax4c.plot(time, WTD_clim[ilat, ilon, :], color='mediumslateblue')
    ax4d.plot(time, WTD_seas[ilat, ilon, :], color='gray')
    ax4d.plot(time, WTD_clim[ilat, ilon, :], color='mediumslateblue')
    first_legend = plt.figlegend(handles=[line1, line2, line3], bbox_to_anchor=(0.944, 0.993), fontsize='10')
    plt.gca().add_artist(first_legend)
    plt.legend(handles=[line4, line5, line6], loc='upper left', bbox_to_anchor=(2.25, 2.2), fontsize='10')

    ax1a.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax1a.xaxis.set_minor_locator(dates.MonthLocator())
    ax1b.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax1b.xaxis.set_minor_locator(dates.MonthLocator())
    ax1c.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax1c.xaxis.set_minor_locator(dates.MonthLocator())
    ax1d.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax1d.xaxis.set_minor_locator(dates.MonthLocator())
    ax2a.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax2a.xaxis.set_minor_locator(dates.MonthLocator())
    ax2b.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax2b.xaxis.set_minor_locator(dates.MonthLocator())
    ax2c.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax2c.xaxis.set_minor_locator(dates.MonthLocator())
    ax2d.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax2d.xaxis.set_minor_locator(dates.MonthLocator())
    ax3a.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax3a.xaxis.set_minor_locator(dates.MonthLocator())
    ax3b.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax3b.xaxis.set_minor_locator(dates.MonthLocator())
    ax3c.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax3c.xaxis.set_minor_locator(dates.MonthLocator())
    ax3d.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax3d.xaxis.set_minor_locator(dates.MonthLocator())
    ax4a.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax4a.xaxis.set_minor_locator(dates.MonthLocator())
    ax4b.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax4b.xaxis.set_minor_locator(dates.MonthLocator())
    ax4c.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax4c.xaxis.set_minor_locator(dates.MonthLocator())
    ax4d.xaxis.set_major_locator(dates.MonthLocator(bymonth=(7, 9)))
    ax4d.xaxis.set_minor_locator(dates.MonthLocator())

    ax3a.xaxis.set_major_formatter(dates.DateFormatter('%m/%y'))
    ax3b.xaxis.set_major_formatter(dates.DateFormatter('%m/%y'))
    ax3c.xaxis.set_major_formatter(dates.DateFormatter('%m/%y'))
    ax3d.xaxis.set_major_formatter(dates.DateFormatter('%m/%y'))
    ax4a.xaxis.set_major_formatter(dates.DateFormatter('%m/%y'))
    ax4b.xaxis.set_major_formatter(dates.DateFormatter('%m/%y'))
    ax4c.xaxis.set_major_formatter(dates.DateFormatter('%m/%y'))
    ax4d.xaxis.set_major_formatter(dates.DateFormatter('%m/%y'))

    ax1a.spines['right'].set_visible(False)
    ax1b.spines['left'].set_visible(False)
    ax1b.tick_params(left=False)
    ax1b.spines['right'].set_visible(False)
    ax1c.spines['left'].set_visible(False)
    ax1c.tick_params(left=False)
    ax1c.spines['right'].set_visible(False)
    ax1d.spines['left'].set_visible(False)
    ax1d.tick_params(left=False)
    ax2a.spines['right'].set_visible(False)
    ax2b.spines['left'].set_visible(False)
    ax2b.tick_params(left=False)
    ax2b.spines['right'].set_visible(False)
    ax2c.spines['left'].set_visible(False)
    ax2c.tick_params(left=False)
    ax2c.spines['right'].set_visible(False)
    ax2d.spines['left'].set_visible(False)
    ax2d.tick_params(left=False)
    ax3a.spines['right'].set_visible(False)
    ax3b.spines['left'].set_visible(False)
    ax3b.tick_params(left=False)
    ax3b.spines['right'].set_visible(False)
    ax3c.spines['left'].set_visible(False)
    ax3c.tick_params(left=False)
    ax3c.spines['right'].set_visible(False)
    ax3d.spines['left'].set_visible(False)
    ax3d.tick_params(left=False)
    ax4a.spines['right'].set_visible(False)
    ax4b.spines['left'].set_visible(False)
    ax4b.tick_params(left=False)
    ax4b.spines['right'].set_visible(False)
    ax4c.spines['left'].set_visible(False)
    ax4c.tick_params(left=False)
    ax4c.spines['right'].set_visible(False)
    ax4d.spines['left'].set_visible(False)
    ax4d.tick_params(left=False)

    d = 1  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=8,
                linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1a.plot([1, 1], [0, 1], transform=ax1a.transAxes, **kwargs)
    ax1b.plot([0, 0], [0, 1], transform=ax1b.transAxes, **kwargs)
    ax1b.plot([1, 1], [0, 1], transform=ax1b.transAxes, **kwargs)
    ax1c.plot([0, 0], [0, 1], transform=ax1c.transAxes, **kwargs)
    ax1c.plot([1, 1], [0, 1], transform=ax1c.transAxes, **kwargs)
    ax1d.plot([0, 0], [0, 1], transform=ax1d.transAxes, **kwargs)
    ax2a.plot([1, 1], [0, 1], transform=ax2a.transAxes, **kwargs)
    ax2b.plot([0, 0], [0, 1], transform=ax2b.transAxes, **kwargs)
    ax2b.plot([1, 1], [0, 1], transform=ax2b.transAxes, **kwargs)
    ax2c.plot([0, 0], [0, 1], transform=ax2c.transAxes, **kwargs)
    ax2c.plot([1, 1], [0, 1], transform=ax2c.transAxes, **kwargs)
    ax2d.plot([0, 0], [0, 1], transform=ax2d.transAxes, **kwargs)
    ax3a.plot([1, 1], [0, 1], transform=ax3a.transAxes, **kwargs)
    ax3b.plot([0, 0], [0, 1], transform=ax3b.transAxes, **kwargs)
    ax3b.plot([1, 1], [0, 1], transform=ax3b.transAxes, **kwargs)
    ax3c.plot([0, 0], [0, 1], transform=ax3c.transAxes, **kwargs)
    ax3c.plot([1, 1], [0, 1], transform=ax3c.transAxes, **kwargs)
    ax3d.plot([0, 0], [0, 1], transform=ax3d.transAxes, **kwargs)
    ax4a.plot([1, 1], [0, 1], transform=ax4a.transAxes, **kwargs)
    ax4b.plot([0, 0], [0, 1], transform=ax4b.transAxes, **kwargs)
    ax4b.plot([1, 1], [0, 1], transform=ax4b.transAxes, **kwargs)
    ax4c.plot([0, 0], [0, 1], transform=ax4c.transAxes, **kwargs)
    ax4c.plot([1, 1], [0, 1], transform=ax4c.transAxes, **kwargs)
    ax4d.plot([0, 0], [0, 1], transform=ax4d.transAxes, **kwargs)

    fig.tight_layout()

    ax5b.set_title('(c)', loc='left', fontsize=10)
    ax6b.set_title('(f)', loc='left', fontsize=10)
    WTDplot_stDry = np.linspace(np.nanmin(WTD[ilat, ilon, :]), WTD_opt_st[ilat, ilon], 5)
    WTDplot_stWet = np.linspace(WTD_opt_st[ilat, ilon], np.nanmax(WTD[ilat, ilon, :]), 5)
    WTDplot_ltDry = np.linspace(np.nanmin(WTD[ilat, ilon, :]), WTD_opt_lt[ilat, ilon], 5)
    WTDplot_ltWet = np.linspace(WTD_opt_lt[ilat, ilon], np.nanmax(WTD[ilat, ilon, :]), 5)
    # x_st = np.linspace(np.nanmin(WTD_sAnom[ilat, ilon, :]), np.nanmax(WTD_sAnom[ilat, ilon, :]), 100)
    # x_lt = np.linspace(np.nanmin(WTD_lAnom[ilat, ilon, :]), np.nanmax(WTD_lAnom[ilat, ilon, :]), 100)
    x_st = np.linspace(-0.064, 0.064, 100)
    x_lt = np.linspace(-0.099, 0.099, 100)
    colors = ['red', 'orangered', 'darksalmon', 'dimgrey', 'cornflowerblue', 'royalblue', 'blue']
    colors_dry = ['red', 'orangered', 'darksalmon', 'dimgrey']
    colors_wet = ['dimgrey', 'cornflowerblue', 'royalblue', 'blue']
    nodes = [0.0, 0.25, 0.40, 0.50, 0.60, 0.75, 1.0]
    nodes_dry = [0.0, 0.50, 0.90, 1.0]
    nodes_wet = [0.0, 0.10, 0.50, 1.0]
    cmap = LinearSegmentedColormap.from_list('mycmap', list(zip(nodes, colors)))
    cmap_dry = LinearSegmentedColormap.from_list('mycmap', list(zip(nodes_dry, colors_dry)))
    cmap_wet = LinearSegmentedColormap.from_list('mycmap', list(zip(nodes_wet, colors_wet)))
    color_dry = iter(cmap_dry(np.linspace(0, 1, 6)[1:]))
    color_wet = iter(cmap_wet(np.linspace(0, 1, 6)))
    for iWTD in WTDplot_stWet:
        c = next(color_wet)
        ax5a.plot(x_st, coef_st[ilat, ilon, 0]*x_st + coef_st[ilat, ilon, 1]*x_st*iWTD, alpha=0.75, c=c)
    color_dry = iter(cmap_dry(np.linspace(0, 1, 6)[1:]))
    color_wet = iter(cmap_wet(np.linspace(0, 1, 6)))
    for iWTD in WTDplot_stDry:
        c = next(color_dry)
        ax5b.plot(x_st, coef_st[ilat, ilon, 0]*x_st + coef_st[ilat, ilon, 1]*x_st*iWTD, alpha=0.75, c=c)
    color_dry = iter(cmap_dry(np.linspace(0, 1, 6)[1:]))
    color_wet = iter(cmap_wet(np.linspace(0, 1, 6)))
    for iWTD in WTDplot_ltWet:
        c = next(color_wet)
        ax6a.plot(x_lt, coef_lt[ilat, ilon, 0]*x_lt + coef_lt[ilat, ilon, 1]*x_lt*iWTD, alpha=0.75, c=c)
    color_dry = iter(cmap_dry(np.linspace(0, 1, 6)[1:]))
    color_wet = iter(cmap_wet(np.linspace(0, 1, 6)))
    for iWTD in WTDplot_ltDry:
        c = next(color_dry)
        ax6b.plot(x_lt, coef_lt[ilat, ilon, 0]*x_lt + coef_lt[ilat, ilon, 1]*x_lt*iWTD, alpha=0.75, c=c)
    ax5a.set_xlim(-0.07, 0.07)
    ax5b.set_xlim(-0.07, 0.07)
    ax5a.set_ylim(-0.42, 0.42)
    ax5b.set_ylim(-0.42, 0.42)
    # ax6.set_xlim(-0.4, 0)
    # ax6.set_ylim(-10, 10)
    ax5b.set_ylabel('SIF$_{n, anom}$\n'+SIFy_unit)
    fig.text(0.28, 0.03, 'WTD$_{anom}$ (m)', ha='center')
    fig.text(0.68, 0.03, 'WTD$_{anom}$ (m)', ha='center')
    ax6a.set_xticks((-0.07, 0, 0.07))
    ax6b.set_xticks((-0.07, 0, 0.07))
    norm_st = mpl.colors.TwoSlopeNorm(vmin=np.nanmin(WTD[ilat, ilon, :]), vcenter=WTD_opt_st[ilat, ilon],
                                      vmax=np.nanmax(WTD[ilat, ilon, :]))
    norm_lt = mpl.colors.TwoSlopeNorm(vmin=np.nanmin(WTD[ilat, ilon, :]), vcenter=WTD_opt_lt[ilat, ilon],
                                      vmax=np.nanmax(WTD[ilat, ilon, :]))
    ax5a.scatter(WTD_sAnom_wet, SIF_sAnom_wet, c=WTD[ilat, ilon, :], cmap=cmap, norm=norm_st, s=5)
    ax5b.scatter(WTD_sAnom_dry, SIF_sAnom_dry, c=WTD[ilat, ilon, :], cmap=cmap,
                 norm=norm_st, s=5)
    ax6a.scatter(WTD_lAnom_wet, SIF_lAnom_wet, c=WTD[ilat, ilon, :], cmap=cmap, norm=norm_lt, s=5)
    ax6b.scatter(WTD_lAnom_dry, SIF_lAnom_dry, c=WTD[ilat, ilon, :], cmap=cmap,
                 norm=norm_lt, s=5)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.1, hspace=0.8)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm_st, cmap=cmap), ax=(ax5a, ax5b), shrink=0.8)
    cb.set_ticks([-0.1, -0.2, -0.3])
    cb.set_ticklabels(['-0.1', '-0.2', '-0.3'])
    cb2 = fig.colorbar(plt.cm.ScalarMappable(norm=norm_lt, cmap=cmap), ax=(ax6a, ax6b), shrink=0.8, label='WTD (m)')
    cb2.set_ticks([-0.1, -0.2, -0.3])
    cb2.set_ticklabels(['-0.1', '-0.2', '-0.3'])

    plt.savefig(save_maps+"GRL/pixels/pixel_v3_Lat"+str(ilat)+"_Lon"+str(ilon)+corr)
    plt.close()
    print(str(p))

len(np.where(slope_st_pixel < 0)[0])
len(np.where(slope_st_pixel > 0)[0])
len(np.where(slope_lt_pixel < 0)[0])
len(np.where(slope_lt_pixel > 0)[0])