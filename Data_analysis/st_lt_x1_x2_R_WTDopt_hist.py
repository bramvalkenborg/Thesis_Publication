# Import packages
from thesis_pub_tools import *
from scipy import optimize
from scipy import stats
import sklearn.preprocessing
import sklearn.linear_model as lm
import statsmodels_adapted.api as sm
import statsmodels.formula.api as smf
import datetime
from patsy import dmatrices
from validation_good_practice.ancillary.metrics import correct_n
# import GitClone.statsmodels.statsmodels.regression.linear_model as sm

# ----------------------------------------------------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------------------------------------------------
# Main settings
peat_min = 0.50
area = 'test'
SIF_PAR = False
SIF_NIR = True
if SIF_NIR and SIF_PAR:
    print('ERROR: no normalization by both PAR and NIR possible')
    exit()
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
elif area == 'test':
    lat_lim = [51, 53]
    lon_lim = [-86, -84]
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
if not SIF_NIR:
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
else:
    if Time_2021 and SP_res == 0.05:
        lon, lat, time = data_analysis_dimension(SIF2021_file5, lat_lim=lat_lim, lon_lim=lon_lim)
        SIFnorm = data_analysis_variable(SIF2021_file5, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        NIR = data_analysis_variable(SIF2021_file5, 'NIR', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        WTD = data_analysis_variable(WTD2021_file5, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        i_unit = ' (m$^{-3}$sr$^{-1}\mu$m$^{-1}$)'
        SIFnorm = SIFnorm/NIR
        corr = corr + '_2021'
        nYears = 4
    elif Time_2021 and SP_res == 0.2:
        lon, lat, time = data_analysis_dimension(SIF2021_file2, lat_lim=lat_lim, lon_lim=lon_lim)
        SIFnorm = data_analysis_variable(SIF2021_file2, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        NIR = data_analysis_variable(SIF2021_file2, 'NIR', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim,
                                     SP_res=SP_res, Time_2021=Time_2021)
        WTD = data_analysis_variable(WTD2021_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        WTD_10y = data_analysis_variable(WTD10y_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        i_unit = ' (m$^{-3}$sr$^{-1}\mu$m$^{-1}$)'
        SIFnorm = SIFnorm / NIR
        corr = corr + '_2021'
        nYears = 4
    elif not Time_2021 and SP_res == 0.05:
        lon, lat, time = data_analysis_dimension(SIF_file5, lat_lim=lat_lim, lon_lim=lon_lim)
        SIFnorm = data_analysis_variable(SIF_file5, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        NIR = data_analysis_variable(SIF_file5, 'NIR', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim,
                                     SP_res=SP_res, Time_2021=Time_2021)
        WTD = data_analysis_variable(WTD_file5, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        PAR = data_analysis_variable(PAR_file5, 'PAR_avg', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        SIFnorm = SIFnorm / NIR
        nYears = 3
    else:
        lon, lat, time = data_analysis_dimension(SIF_file2, lat_lim=lat_lim, lon_lim=lon_lim)
        SIFnorm = data_analysis_variable(SIF_file2, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        NIR = data_analysis_variable(SIF_file2, 'NIR', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim,
                                     SP_res=SP_res, Time_2021=Time_2021)
        WTD = data_analysis_variable(WTD_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        PAR = data_analysis_variable(PAR_file2, 'PAR_avg', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        WTD_10y = data_analysis_variable(WTD10y_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
        SIFnorm = SIFnorm / NIR
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


def CalcP_f_statistic(model, model_restricted, n_corr, N_pred):
    ssr = model.ssr
    ssr_r = model_restricted.ssr
    dfn = N_pred
    dfd = n_corr - N_pred - 1
    f = (ssr_r-ssr)/ssr
    p_value = stats.f.cdf(f, dfn, dfd)
    return p_value


coef_s = np.zeros((len(lat), len(lon), 2))
p_values_s = np.zeros((len(lat), len(lon), 2))
fp_values_s = np.zeros((len(lat), len(lon)))
WTD_opt_s = np.zeros((len(lat), len(lon)))
n_corr_s = np.zeros((len(lat), len(lon)))
rho_s = np.zeros((len(lat), len(lon)))
n_s = np.zeros((len(lat), len(lon)))
Rsq_s = np.zeros((len(lat), len(lon)))
for lat1 in range(len(lat)):
    print('Calibration: ' + str(round((lat1 / len(lat)) * 100, 2)) + '%')
    for lon1 in range(len(lon)):
        x1 = WTD_sAnom[lat1, lon1, :]
        x2 = WTD[lat1, lon1, :]
        y = SIF_sAnom[lat1, lon1, :]
        if not np.isnan(x1).all() and not np.isnan(x2).all() and not np.isnan(y).all():
            df = pd.DataFrame(np.transpose(np.asarray([x1, x2, y])), index=time, columns=['WTD_sAnom', 'WTD', 'SIF_sAnom'])
            df.dropna(inplace=True)
            n_corr_s[lat1, lon1], rho_s[lat1, lon1] = correct_n(df)
            n_s[lat1, lon1] = len(np.where(np.logical_not(np.isnan(SIF_sAnom[lat1, lon1, :])))[0])
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
            model = sm.OLS(y, X[:, 1:3]).fit(nobs_corr=n_corr_s[lat1, lon1])
            # model = sm.OLS(y, X[:, 1:3]).fit()
            # Y_r, X_r = dmatrices('y ~ x1')
            # model_r = sm.OLS(y, X_r[:, 1:3]).fit()
            # fp_values_s[lat1, lon1] = CalcP_f_statistic(model, model_r, n_corr[lat1, lon1], 2)
            fp_values_s[lat1, lon1] = model.f_pvalue
            p_values_s[lat1, lon1, :] = model.pvalues
            if fp_values_s[lat1, lon1] <= p:
                coef_s[lat1, lon1, :] = model.params
                WTD_opt_s[lat1, lon1] = -coef_s[lat1, lon1, 0]/coef_s[lat1, lon1, 1]
                Rsq_s[lat1, lon1] = model.rsquared
            else:
                coef_s[lat1, lon1, :] = model.params
                WTD_opt_s[lat1, lon1] = -coef_s[lat1, lon1, 0]/coef_s[lat1, lon1, 1]
                Rsq_s[lat1, lon1] = model.rsquared
                # p_values_s[lat1, lon1, :] = np.nan
                # n_corr[lat1, lon1] = np.nan
        else:
            coef_s[lat1, lon1, :] = np.nan
            fp_values_s[lat1, lon1] = np.nan
            WTD_opt_s[lat1, lon1] = np.nan
            n_corr_s[lat1, lon1] = np.nan
            n_s[lat1, lon1] = np.nan
            rho_s[lat1, lon1] = np.nan
            Rsq_s[lat1, lon1] = np.nan


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
WTD_opt_l = np.zeros((len(lat), len(lon)))
n_corr_l = np.zeros((len(lat), len(lon)))
rho_l = np.zeros((len(lat), len(lon)))
n_l = np.zeros((len(lat), len(lon)))
Rsq_l = np.zeros((len(lat), len(lon)))
for lat1 in range(len(lat)):
    print('Calibration: ' + str(round((lat1 / len(lat)) * 100, 2)) + '%')
    for lon1 in range(len(lon)):
        x1 = WTD_lAnom[lat1, lon1, :]
        x2 = WTD[lat1, lon1, :]
        y = SIF_lAnom[lat1, lon1, :]
        if not np.isnan(x1).all() and not np.isnan(x2).all() and not np.isnan(y).all():
            df = pd.DataFrame(np.transpose(np.asarray([x1, x2, y])), index=time, columns=['WTD_lAnom', 'WTD', 'SIF_lAnom'])
            df.dropna(inplace=True)
            n_corr_l[lat1, lon1], rho_l[lat1, lon1] = correct_n(df)
            n_l[lat1, lon1] = len(np.where(np.logical_not(np.isnan(SIF_lAnom[lat1, lon1, :])))[0])
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
            model = sm.OLS(y, X[:, 1:3]).fit(nobs_corr=n_corr_l[lat1, lon1])
            # model = sm.OLS(y, X[:, 1:3]).fit()
            # Y_r, X_r = dmatrices('y ~ x1')
            # model_r = sm.OLS(y, X_r[:, 1:3]).fit()
            # fp_values_l[lat1, lon1] = CalcP_f_statistic(model, model_r, n_corr[lat1, lon1], 2)
            fp_values_l[lat1, lon1] = model.f_pvalue
            p_values_l[lat1, lon1, :] = model.pvalues
            if fp_values_l[lat1, lon1] <= p:
                coef_l[lat1, lon1, :] = model.params
                WTD_opt_l[lat1, lon1] = -coef_l[lat1, lon1, 0]/coef_l[lat1, lon1, 1]
                Rsq_l[lat1, lon1] = model.rsquared
            else:
                coef_l[lat1, lon1, :] = model.params
                WTD_opt_l[lat1, lon1] = -coef_l[lat1, lon1, 0]/coef_l[lat1, lon1, 1]
                Rsq_l[lat1, lon1] = model.rsquared
                # p_values_l[lat1, lon1, :] = np.nan
                # n_corr[lat1, lon1] = np.nan
        else:
            coef_l[lat1, lon1, :] = np.nan
            fp_values_l[lat1, lon1] = np.nan
            WTD_opt_l[lat1, lon1] = np.nan
            n_corr_l[lat1, lon1] = np.nan
            rho_l[lat1, lon1] = np.nan
            n_l[lat1, lon1] = np.nan
            Rsq_l[lat1, lon1] = np.nan
"""
plt.figure()
plt.scatter(n_corr_s, rho_s)
plt.xlabel('n$_{corr}$')
plt.ylabel('rho')
plt.title('Short-term')

plt.figure()
plt.hist(np.resize(rho_s, rho_s.size))
plt.xlim(-0.1, 0.9)
plt.xlabel('rho$_{short term}$')
plt.figure()
plt.hist(np.resize(rho_l, rho_l.size))
plt.xlim(-0.1, 0.9)
plt.xlabel('rho$_{long term}$')

plt.figure()
plt.scatter(n_corr_l, rho_l)
plt.xlabel('n$_{corr}$')
plt.ylabel('rho')
plt.title('Long-term')


n_corr_s = np.resize(n_corr_s, n_corr_s.size)
n_corr_l = np.resize(n_corr_l, n_corr_l.size)
n_s = np.resize(n_s, n_s.size)
n_l = np.resize(n_l, n_l.size)
plt.figure()
plt.hist(n_corr_s, histtype='step', label='Corrected')
plt.hist(n_s, histtype='step', label='Initial')
plt.title('Short-term')
plt.xlabel('Number of datapoints per pixel')
plt.legend()
plt.figure()
plt.hist(n_corr_l, histtype='step', label='Corrected')
plt.hist(n_l, histtype='step', label='Initial')
plt.title('Long-term')
plt.legend()
plt.xlabel('Number of datapoints per pixel')
"""

"""
plt.figure()
plt.hlines(2.71, 0, 150)
rho_mean = np.nanmean(np.nanmean(rho, axis=0), axis=0)
plt.plot(range(0, 150, 10), rho_mean)
tau_mean = np.nanmean(np.nanmean(tau, axis=0), axis=0)
plt.plot(range(0, 150, 10), tau_mean)
plt.xlabel('n_lags')
plt.ylabel('rho')
plt.grid()
"""
# ----------------------------------------------------------------------------------------------------------------------
# WTDopt and x2 graphs
# ----------------------------------------------------------------------------------------------------------------------
# x1 values: all and only-significant
x1_st_sig = rm_nan(np.where(fp_values_s <= 0.05, coef_s[:, :, 0], np.nan))
x1_st = rm_nan(coef_s[:, :, 0])
x1_st_hist = [x1_st, x1_st_sig]

x1_lt_sig = rm_nan(np.where(fp_values_l <= 0.05, coef_l[:, :, 0], np.nan))
x1_lt = rm_nan(coef_l[:, :, 0])
x1_lt_hist = [x1_lt, x1_lt_sig]

# x2 values
x2_st_sig = rm_nan(np.where(fp_values_s <= 0.05, coef_s[:, :, 1], np.nan))
x2_st = rm_nan(coef_s[:, :, 1])
x2_st_hist = [x2_st, x2_st_sig]

x2_lt_sig = rm_nan(np.where(fp_values_l <= 0.05, coef_l[:, :, 1], np.nan))
x2_lt = rm_nan(coef_l[:, :, 1])
x2_lt_hist = [x2_lt, x2_lt_sig]

# WTD_opt values
WTDopt_st_sig = rm_nan(np.where(fp_values_s <= 0.05, WTD_opt_s, np.nan))
WTDopt_st = rm_nan(WTD_opt_s)
WTDopt_st_hist = [WTDopt_st, WTDopt_st_sig]

WTDopt_lt_sig = rm_nan(np.where(fp_values_l <= 0.05, WTD_opt_l, np.nan))
WTDopt_lt = rm_nan(WTD_opt_l)
WTDopt_lt_hist = [WTDopt_lt, WTDopt_lt_sig]

# R values
Rsq_st_sig = rm_nan(np.where(fp_values_s <= 0.05, Rsq_s, np.nan))
Rsq_st = rm_nan(Rsq_s)
Rsq_st_hist = [Rsq_st, Rsq_st_sig]

Rsq_lt_sig = rm_nan(np.where(fp_values_l <= 0.05, Rsq_l, np.nan))
Rsq_lt = rm_nan(Rsq_l)
Rsq_lt_hist = [Rsq_lt, Rsq_lt_sig]

WTD_AvgMax = np.nanmedian(np.nanmax(WTD, axis=2))
WTD_AvgMin = np.nanmedian(np.nanmin(WTD, axis=2))

Hyp_true_s = (WTD_opt_s <= WTD_AvgMax) & (WTD_opt_s >= WTD_AvgMin) & (coef_s[:, :, 1] < 0)
Hyp_true_l = (WTD_opt_l <= WTD_AvgMax) & (WTD_opt_l >= WTD_AvgMin) & (coef_l[:, :, 1] < 0)

Hyp_true_s = (WTD_opt_s <= np.nanmax(WTD, axis=2)) & (WTD_opt_s >= np.nanmin(WTD, axis=2)) & (coef_s[:, :, 1] < 0)
Hyp_true_l = (WTD_opt_l <= np.nanmax(WTD, axis=2)) & (WTD_opt_l >= np.nanmin(WTD, axis=2)) & (coef_l[:, :, 1] < 0)

Hyp_true_s = (WTD_opt_s <= np.nanmax(WTD, axis=2)) & (WTD_opt_s >= np.nanmin(WTD, axis=2))
Hyp_true_l = (WTD_opt_l <= np.nanmax(WTD, axis=2)) & (WTD_opt_l >= np.nanmin(WTD, axis=2))

Hyp_true_s_sig = (WTD_opt_s <= WTD_AvgMax) & (WTD_opt_s >= WTD_AvgMin) & (coef_s[:, :, 1] < 0) & (fp_values_s <= 0.05)
Hyp_true_l_sig = (WTD_opt_l <= WTD_AvgMax) & (WTD_opt_l >= WTD_AvgMin) & (coef_l[:, :, 1] < 0) & (fp_values_l <= 0.05)


fig = plt.figure(figsize=(10, 4.2))
gs = fig.add_gridspec(2, 4)
ax1 = fig.add_subplot(gs[0, 0])
ax1b = ax1.twinx()
ax2 = fig.add_subplot(gs[0, 1])
ax2b = ax2.twinx()
ax3 = fig.add_subplot(gs[0, 2])
ax3b = ax3.twinx()
ax4 = fig.add_subplot(gs[0, 3])
ax4b = ax4.twinx()
ax5 = fig.add_subplot(gs[1, 0])
ax5b = ax5.twinx()
ax6 = fig.add_subplot(gs[1, 1])
ax6b = ax6.twinx()
ax7 = fig.add_subplot(gs[1, 2])
ax7b = ax7.twinx()
ax8 = fig.add_subplot(gs[1, 3])
ax8b = ax8.twinx()

ax1.set_title('(a)', loc='left')
ax2.set_title('(b)', loc='left')
ax3.set_title('(c)', loc='left')
ax4.set_title('(d)', loc='left')
ax5.set_title('(e)', loc='left')
ax6.set_title('(f)', loc='left')
ax7.set_title('(g)', loc='left')
ax8.set_title('(h)', loc='left')

ax1.tick_params(axis='y', labelcolor='tab:orange')
ax2.tick_params(axis='y', labelcolor='tab:orange')
ax3.tick_params(axis='y', labelcolor='tab:orange')
ax4.tick_params(axis='y', labelcolor='tab:orange')
ax5.tick_params(axis='y', labelcolor='tab:orange')
ax6.tick_params(axis='y', labelcolor='tab:orange')
ax7.tick_params(axis='y', labelcolor='tab:orange')
ax8.tick_params(axis='y', labelcolor='tab:orange')

ax1b.tick_params(axis='y', labelcolor='tab:blue')
ax2b.tick_params(axis='y', labelcolor='tab:blue')
ax3b.tick_params(axis='y', labelcolor='tab:blue')
ax4b.tick_params(axis='y', labelcolor='tab:blue')
ax5b.tick_params(axis='y', labelcolor='tab:blue')
ax6b.tick_params(axis='y', labelcolor='tab:blue')
ax7b.tick_params(axis='y', labelcolor='tab:blue')
ax8b.tick_params(axis='y', labelcolor='tab:blue')

# ax1.hist(x1_st_sig, bins=20, range=(-25, 25), histtype='step', color=['tab:orange'], label=['Only SGFNT'])
# ax2.hist(x2_st_sig, bins=20, range=(-100, 100), histtype='step', color=['tab:orange'])
ax1.hist(x1_st_sig, bins=20, range=(np.nanmin(x1_st_sig), np.nanmax(x1_st_sig)), histtype='step', color=['tab:orange'], label=['Only SGFNT'])
ax2.hist(x2_st_sig, bins=20, range=(np.nanmin(x2_st_sig), np.nanmax(x2_st_sig)), histtype='step', color=['tab:orange'])
ax3.hist(WTDopt_st_sig, bins=20, range=(-1.4, 0.8), histtype='step', color=['tab:orange'])
ax4.hist(Rsq_st_sig, bins=20, range=(0, 0.70), histtype='step', color=['tab:orange'])

# ax5.hist(x1_lt_sig, bins=20, range=(-8, 8), histtype='step', color=['tab:orange'])
# ax6.hist(x2_lt_sig, bins=20, range=(-35, 35), histtype='step', color=['tab:orange'])
ax5.hist(x1_lt_sig, bins=20, range=(np.nanmin(x1_lt_sig), np.nanmax(x1_lt_sig)), histtype='step', color=['tab:orange'])
ax6.hist(x2_lt_sig, bins=20, range=(-0.20, 0.20), histtype='step', color=['tab:orange'])
ax7.hist(WTDopt_lt_sig, bins=20, range=(-1.4, 0.8), histtype='step', color=['tab:orange'])
ax8.hist(Rsq_lt_sig, bins=20, range=(0, 0.70), histtype='step', color=['tab:orange'])

# ax1b.hist(x1_st, bins=20, range=(-25, 25), histtype='step', color=['tab:blue'], label=['All'])
# ax2b.hist(x2_st, bins=20, range=(-100, 100), histtype='step', color=['tab:blue'])
ax1b.hist(x1_st, bins=20, range=(np.nanmin(x1_st), np.nanmax(x1_st)), histtype='step', color=['tab:blue'], label=['All'])
ax2b.hist(x2_st, bins=20, range=(np.nanmin(x2_st), np.nanmax(x2_st)), histtype='step', color=['tab:blue'])
ax3b.hist(WTDopt_st, bins=20, range=(-1.4, 0.8), histtype='step', color=['tab:blue'])
ax4b.hist(Rsq_st, bins=20, range=(0, 0.70), histtype='step', color=['tab:blue'])

# ax5b.hist(x1_lt, bins=20, range=(-8, 8), histtype='step', color=['tab:blue'])
# ax6b.hist(x2_lt, bins=20, range=(-35, 35), histtype='step', color=['tab:blue'])
ax5b.hist(x1_lt, bins=20, range=(np.nanmin(x1_lt), np.nanmax(x1_lt)), histtype='step', color=['tab:blue'])
ax6b.hist(x2_lt, bins=20, range=(-0.20, 0.20), histtype='step', color=['tab:blue'])
ax7b.hist(WTDopt_lt, bins=20, range=(-1.4, 0.8), histtype='step', color=['tab:blue'])
ax8b.hist(Rsq_lt, bins=20, range=(0, 0.70), histtype='step', color=['tab:blue'])

ax5.set_xlabel(r'$\alpha$ (sr$^{-1}$mm$^{-2}$)')
ax6.set_xlabel(r'$\beta$ (sr$^{-1}$mm$^{-3}$)')
ax7.set_xlabel(r'$WTD_{opt}$ (m)')
ax8.set_xlabel(r'$R^2$ (-)')

ax3b.vlines(x=WTD_AvgMax, ymin=0, ymax=1900, linestyles='dashed')
ax3b.vlines(x=WTD_AvgMin, ymin=0, ymax=1900, linestyles='dashed')
ax3b.vlines(x=np.nanmedian(WTD_opt_s), ymin=0, ymax=1900, linestyles='dotted')
ax7b.vlines(x=WTD_AvgMax, ymin=0, ymax=1300, linestyles='dashed')
ax7b.vlines(x=WTD_AvgMin, ymin=0, ymax=1300, linestyles='dashed')
ax7b.vlines(x=np.nanmedian(WTD_opt_l), ymin=0, ymax=1300, linestyles='dotted')

# ax1.set_ylim(0, 1500)
# ax2.set_ylim(0, 1500)
# ax3.set_ylim(0, 1100)
# ax4.set_ylim(0, 1500)

fig.tight_layout()
fig.subplots_adjust(left=0.11, bottom=0.14, right=0.95)
fig.text(0.04, 0.5, 'Number of pixels', va='center', rotation='vertical')
# fig.text(0.25, 0.04, 'WTD$_{opt}$ (m)', va='center')
# fig.text(0.70, 0.04, r'$\beta$ (sr$^{-1}$mm$^{-3}$)', va='center')
fig.legend()
plt.savefig(save_maps+"GRL/Hist_WTDopt_x2"+corr)


WTD_opt_s_significant = np.where(fp_values_s <= 0.05, WTD_opt_s, np.nan)
WTD_opt_s_Nsignificant = np.where(fp_values_s > 0.05, WTD_opt_s, np.nan)
WTD_opt_l_significant = np.where(fp_values_l <= 0.05, WTD_opt_l, np.nan)
WTD_opt_l_Nsignificant = np.where(fp_values_l > 0.05, WTD_opt_l, np.nan)

x1_s_significant = np.where(fp_values_s <= 0.05, coef_s[:, :, 0], np.nan)
x1_s_Nsignificant = np.where(fp_values_s > 0.05, coef_s[:, :, 0], np.nan)
x1_l_significant = np.where(fp_values_l <= 0.05, coef_l[:, :, 0], np.nan)
x1_l_Nsignificant = np.where(fp_values_l > 0.05, coef_l[:, :, 0], np.nan)

x2_s_significant = np.where(fp_values_s <= 0.05, coef_s[:, :, 1], np.nan)
x2_s_Nsignificant = np.where(fp_values_s > 0.05, coef_s[:, :, 1], np.nan)
x2_l_significant = np.where(fp_values_l <= 0.05, coef_l[:, :, 1], np.nan)
x2_l_Nsignificant = np.where(fp_values_l > 0.05, coef_l[:, :, 1], np.nan)

print(len(np.where((WTD_opt_s_significant <= WTD_AvgMax) & (WTD_opt_s_significant >= WTD_AvgMin))[0]) / len(np.where(np.logical_not(np.isnan(WTD_opt_s_significant)))[0]))
print(len(np.where((WTD_opt_l_significant <= WTD_AvgMax) & (WTD_opt_l_significant >= WTD_AvgMin))[0]) / len(np.where(np.logical_not(np.isnan(WTD_opt_l_significant)))[0]))
print(len(np.where((WTD_opt_s_Nsignificant < WTD_AvgMax) & (WTD_opt_s_Nsignificant > WTD_AvgMin))[0]) / len(np.where(np.logical_not(np.isnan(WTD_opt_s_Nsignificant)))[0]))
print(len(np.where((WTD_opt_l_Nsignificant < WTD_AvgMax) & (WTD_opt_l_Nsignificant > WTD_AvgMin))[0]) / len(np.where(np.logical_not(np.isnan(WTD_opt_l_Nsignificant)))[0]))

print(len(np.where(x1_s_significant < 0)[0]) / len(np.where(np.logical_not(np.isnan(x1_s_significant)))[0]))
print(len(np.where(x1_l_significant < 0)[0]) / len(np.where(np.logical_not(np.isnan(x1_l_significant)))[0]))
print(len(np.where(x1_s_Nsignificant < 0)[0]) / len(np.where(np.logical_not(np.isnan(x1_s_Nsignificant)))[0]))
print(len(np.where(x1_l_Nsignificant < 0)[0]) / len(np.where(np.logical_not(np.isnan(x1_l_Nsignificant)))[0]))

print(len(np.where(x2_s_significant < 0)[0]) / len(np.where(np.logical_not(np.isnan(x2_s_significant)))[0]))
print(len(np.where(x2_l_significant < 0)[0]) / len(np.where(np.logical_not(np.isnan(x2_l_significant)))[0]))
print(len(np.where(x2_s_Nsignificant < 0)[0]) / len(np.where(np.logical_not(np.isnan(x2_s_Nsignificant)))[0]))
print(len(np.where(x2_l_Nsignificant < 0)[0]) / len(np.where(np.logical_not(np.isnan(x2_l_Nsignificant)))[0]))

print(len(np.where(fp_values_s <= 0.05)[0]))
print(len(np.where(fp_values_l <= 0.05)[0]))
