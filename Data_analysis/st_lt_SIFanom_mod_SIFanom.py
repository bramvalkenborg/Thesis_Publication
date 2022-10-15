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

# ----------------------------------------------------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------------------------------------------------
# Main settings
peat_min = 0.50
area = 'Hudson_Bay'
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
    # OW = data_analysis_variable(OW2021_file5, 'OW', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    i_unit = ' (m$^{-3}$sr$^{-1}\mu$m$^{-1}$)'
    corr = corr + '_2021'
    nYears = 4
elif Time_2021 and SP_res == 0.2:
    lon, lat, time = data_analysis_dimension(SIF2021_file2, lat_lim=lat_lim, lon_lim=lon_lim)
    SIFnorm = data_analysis_variable(SIF2021_file2, SIF_variable, peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD = data_analysis_variable(WTD2021_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    WTD_10y = data_analysis_variable(WTD10y_file2, 'WTD', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res, Time_2021=Time_2021)
    # OW = data_analysis_variable(OW2021_file2, 'OW', peat_min=peat_min, lat_lim=lat_lim, lon_lim=lon_lim, SP_res=SP_res,
    #                             Time_2021=Time_2021)
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


def lm(lat, lon, time, WTD_Anom, WTD, SIF_Anom, SIFnorm, nYears):
    x1 = WTD_Anom
    x2 = WTD
    y = SIF_Anom
    days_wet = np.zeros((len(lat), len(lon)))
    days_dry = np.zeros((len(lat), len(lon)))
    slope = np.zeros((len(lat), len(lon), len(time)))
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
        rsquared_adj = model.rsquared_adj
        if fp_values <= 0.05:
            for time1 in range(len(time)):
                slope[:, :, time1] = coef[1] * WTD[:, :, time1] + coef[0]
            for lat1 in range(len(lat)):
                for lon1 in range(len(lon)):
                    if not np.isnan(SIFnorm[lat1, lon1, :]).all():
                        days_wet[lat1, lon1] = (len(np.where(slope[lat1, lon1, :] < 0)[0]) / nYears) * 8
                        days_dry[lat1, lon1] = (len(np.where(slope[lat1, lon1, :] > 0)[0]) / nYears) * 8
                    else:
                        days_wet[lat1, lon1] = np.nan
                        days_dry[lat1, lon1] = np.nan
        else:
            coef = np.nan
            slope[:, :, :] = np.nan
            days_wet[:, :] = np.nan
            days_dry[:, :] = np.nan
    else:
        coef = np.nan
        slope[:, :, :] = np.nan
        days_wet[:, :] = np.nan
        days_dry[:, :] = np.nan
    return coef, slope, fp_values, rsquared_adj, days_wet, days_dry


def lme(lat, lon, time, WTD_Anom, WTD, SIF_Anom, SIFnorm, nYears):
    x1 = WTD_Anom
    x2 = WTD
    y = SIF_Anom
    days_wet = np.zeros((len(lat), len(lon)))
    days_dry = np.zeros((len(lat), len(lon)))
    slope = np.zeros((len(lat), len(lon), len(time)))
    pixel_id_2D = np.arange(1, days_wet.size + 1, dtype='int').reshape((len(lat), len(lon)))
    pixel_id = np.repeat(pixel_id_2D[:, :, np.newaxis], len(WTD[0, 0, :]), axis=2)
    WTD_opt = np.zeros((len(lat), len(lon)))
    re_coef = np.zeros((len(lat), len(lon)))
    if not np.isnan(x1).all() and not np.isnan(x2).all() and not np.isnan(y).all():
        x1a = x1[np.logical_not(np.isnan(x1))]
        x2a = x2[np.logical_not(np.isnan(x1))]
        ya = y[np.logical_not(np.isnan(x1))]
        pixel_id_a = pixel_id[np.logical_not(np.isnan(x1))]
        x1b = x1a[np.logical_not(np.isnan(x2a))]
        x2b = x2a[np.logical_not(np.isnan(x2a))]
        yb = ya[np.logical_not(np.isnan(x2a))]
        pixel_id_b = pixel_id_a[np.logical_not(np.isnan(x2a))]
        x1 = x1b[np.logical_not(np.isnan(yb))]
        x2 = x2b[np.logical_not(np.isnan(yb))]
        y = yb[np.logical_not(np.isnan(yb))]
        pixel_id_c = pixel_id_b[np.logical_not(np.isnan(yb))]
        data = pd.DataFrame(np.transpose(np.array([pixel_id_c, y, x1, x2])), columns=['id', 'SIFanom', 'WTDanom', 'WTD'])
        md = smf.mixedlm("SIFanom ~ WTDanom + WTDanom:WTD", data, re_formula="0+WTDanom:WTD", groups=data["id"])
        model = md.fit()
        # fp_values = model.f_pvalue
        p_values = model.pvalues
        coef = model.params
        re = model.random_effects
        if p_values[0] <= p and p_values[1] <= p and p_values[2] <= p:
            for ilat in range(len(lat)):
                print(str(round((ilat / len(lat)), 2) * 100) + '%')
                for ilon in range(len(lon)):
                    if not np.isnan(SIF_Anom[ilat, ilon, :]).all():
                        pix = pixel_id_2D[ilat, ilon]
                        re_coef[ilat, ilon] = coef[2]+re[pix].values
                        WTD_opt[ilat, ilon] = -coef[1]/(coef[2]+re[pix].values)
                    else:
                        WTD_opt[ilat, ilon] = np.nan
                        re_coef[ilat, ilon] = np.nan
            for time1 in range(len(time)):
                slope[:, :, time1] = re_coef[:, :] * WTD[:, :, time1] + coef[1]
            for ilat in range(len(lat)):
                print(str(round((ilat / len(lat)), 2) * 100) + '%')
                for ilon in range(len(lon)):
                    if not np.isnan(SIFnorm[ilat, ilon, :]).all():
                        days_wet[ilat, ilon] = (len(np.where(slope[ilat, ilon, :] < 0)[0]) / nYears) * 8
                        days_dry[ilat, ilon] = (len(np.where(slope[ilat, ilon, :] > 0)[0]) / nYears) * 8
                    else:
                        days_wet[ilat, ilon] = np.nan
                        days_dry[ilat, ilon] = np.nan
        else:
            coef = np.nan
            slope[:, :, :] = np.nan
            days_wet[:, :] = np.nan
            days_dry[:, :] = np.nan
            WTD_opt[:, :] = np.nan
            re_coef[:, :] = np.nan
    else:
        coef = np.nan
        slope[:, :, :] = np.nan
        days_wet[:, :] = np.nan
        days_dry[:, :] = np.nan
        WTD_opt[:, :] = np.nan
        re_coef[:, :] = np.nan
    return coef, slope, days_wet, days_dry, WTD_opt, re_coef


def func(WTD, WTD_anom, x1, x2):
    x2 = np.repeat(x2[:, :, np.newaxis], len(WTD[0, 0, :]), axis=2)
    SIFnorm_anom = x1*WTD_anom+x2*WTD_anom*WTD
    return SIFnorm_anom


# coef_l, slope_l, fp_values_l, rsquared_adj_l, days_wet_NH_l, days_dry_NH_l = lm(lat, lon, time, WTD_lAnom, WTD, SIF_lAnom, SIFnorm, nYears)
# coef_s, slope_s, fp_values_s, rsquared_adj_s, days_wet_NH_s, days_dry_NH_s = lm(lat, lon, time, WTD_lAnom, WTD, SIF_lAnom, SIFnorm, nYears)
coef_l, slope_l, days_wet_NH_l, days_dry_NH_l, WTD_opt_l, re_coef_l = lme(lat, lon, time, WTD_sAnom, WTD, SIF_sAnom, SIFnorm, nYears)
coef_s, slope_s, days_wet_NH_s, days_dry_NH_s, WTD_opt_s, re_coef_s = lme(lat, lon, time, WTD_sAnom, WTD, SIF_sAnom, SIFnorm, nYears)


SIFnorm_stAnom_mod = func(WTD, WTD_sAnom, coef_s[1], re_coef_s)
SIFnorm_ltAnom_mod = func(WTD, WTD_lAnom, coef_l[1], re_coef_l)

density_plots(SIF_lAnom, SIFnorm_ltAnom_mod, axis=[-1, 1, -1, 1], xlabel='SIF$_{ltAnom}$', ylabel='SIF$_{ltAnom, mod}$', title='NH based', grid=True)
density_plots(SIF_sAnom, SIFnorm_stAnom_mod, axis=[-1, 1, -1, 1], xlabel='SIF$_{stAnom}$', ylabel='SIF$_{stAnom, mod}$', title='NH based', grid=True)


SIFnorm_stAnom_mod_dry = np.where(WTD < np.nanquantile(WTD, 0.90), SIFnorm_stAnom_mod, np.nan)
SIFnorm_ltAnom_mod_dry = np.where(WTD < np.nanquantile(WTD, 0.90), SIFnorm_ltAnom_mod, np.nan)
SIFnorm_stAnom_dry = np.where(WTD < np.nanquantile(WTD, 0.90), SIF_sAnom, np.nan)
SIFnorm_ltAnom_dry = np.where(WTD < np.nanquantile(WTD, 0.90), SIF_lAnom, np.nan)

SIFnorm_stAnom_mod_wet = np.where(WTD > np.nanquantile(WTD, 0.10), SIFnorm_stAnom_mod, np.nan)
SIFnorm_ltAnom_mod_wet = np.where(WTD > np.nanquantile(WTD, 0.10), SIFnorm_ltAnom_mod, np.nan)
SIFnorm_stAnom_wet = np.where(WTD > np.nanquantile(WTD, 0.10), SIF_sAnom, np.nan)
SIFnorm_ltAnom_wet = np.where(WTD > np.nanquantile(WTD, 0.10), SIF_lAnom, np.nan)


density_plots(SIFnorm_stAnom_dry, SIFnorm_stAnom_mod_dry)

density_plots(SIFnorm_ltAnom_dry, SIFnorm_ltAnom_mod_dry)

density_plots(SIFnorm_stAnom_wet, SIFnorm_stAnom_mod_wet)

density_plots(SIFnorm_ltAnom_wet, SIFnorm_ltAnom_mod_wet)
