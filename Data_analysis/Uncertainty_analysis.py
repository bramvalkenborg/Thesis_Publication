# Import packages
from thesis_pub_tools import *
import datetime

# ----------------------------------------------------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------------------------------------------------
# Main settings
peat_min = 0.50
area = 'test'
SIF_PAR = True
Time_2021 = True
SP_res = 0.2
p = 1
output_file = '/data/leuven/336/vsc33653/OUTPUT_pub/uncertainty_WTDopt.txt'
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
    lat_lim = [63, 64]
    lon_lim = [-117, -115]
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


# Bootstrap analysis
WTDopt_s = np.zeros((len(lat), len(lon)))
WTD_CI_5_s = np.zeros((len(lat), len(lon)))
WTD_CI_95_s = np.zeros((len(lat), len(lon)))
for ilat in range(len(lat)):
    print(str(ilat))
    for ilon in range(len(lon)):
        if not np.isnan(WTD).all() and not np.isnan(SIF_sAnom).all() and not np.isnan(WTD_sAnom).all():
            WTDopt_mean_pix_s, WTD_CI_5_pix_s, WTD_CI_95_pix_s = Bootstrap_uncertainty(SIF_sAnom[ilat, ilon, :], WTD_sAnom[ilat, ilon, :], WTD[ilat, ilon, :], 10, p, time)
            WTDopt_s[ilat, ilon] = WTDopt_mean_pix_s
            WTD_CI_5_s[ilat, ilon] = WTD_CI_5_pix_s
            WTD_CI_95_s[ilat, ilon] = WTD_CI_95_pix_s
        else:
            WTDopt_s[ilat, ilon] = np.nan
            WTD_CI_5_s[ilat, ilon] = np.nan
            WTD_CI_95_s[ilat, ilon] = np.nan

WTD_CI_5_mean_s = np.nanmean(WTD_CI_5_s)
WTD_CI_95_mean_s = np.nanmean(WTD_CI_95_s)
WTDopt_mean_s = np.nanmean(WTDopt_s)

# Bootstrap analysis: Long term
WTDopt_l = np.zeros((len(lat), len(lon)))
WTD_CI_5_l = np.zeros((len(lat), len(lon)))
WTD_CI_95_l = np.zeros((len(lat), len(lon)))
for ilat in range(len(lat)):
    print('Calibration: ' + str(round((ilat / len(lat)) * 100, 2)) + '%')
    for ilon in range(len(lon)):
        if not np.isnan(WTD).all() and not np.isnan(SIF_lAnom).all() and not np.isnan(WTD_lAnom).all():
            WTDopt_mean_pix_l, WTD_CI_5_pix_l, WTD_CI_95_pix_l = Bootstrap_uncertainty(SIF_lAnom[ilat, ilon, :], WTD_lAnom[ilat, ilon, :], WTD[ilat, ilon, :], 10, p, time)
            WTDopt_l[ilat, ilon] = WTDopt_mean_pix_l
            WTD_CI_5_l[ilat, ilon] = WTD_CI_5_pix_l
            WTD_CI_95_l[ilat, ilon] = WTD_CI_95_pix_l
        else:
            WTDopt_l[ilat, ilon] = np.nan
            WTD_CI_5_l[ilat, ilon] = np.nan
            WTD_CI_95_l[ilat, ilon] = np.nan

WTD_CI_5_mean_l = np.nanmean(WTD_CI_5_l)
WTD_CI_95_mean_l = np.nanmean(WTD_CI_95_l)
WTDopt_mean_l = np.nanmean(WTDopt_l)

with open(output_file, 'a') as f:
    f.write('--------------------------------------------------------------------------------------\n')
    f.write('Uncertainty WTDopt                  Date:'+str(datetime.datetime.today())+'\n')
    f.write('--------------------------------------------------------------------------------------\n')
    f.wirte('Settings:\n')
    f.write('   peat mask: ' + str(peat_min) + '\n')
    f.write('   area: '+str(area)+'\n')
    f.write('   p: '+str(p)+'\n')
    f.write('   Normalized by PAR?: '+str(SIF_PAR)+'\n')
    f.write('   2021 included?: '+str(Time_2021)+'\n')
    f.write('   Spatial resolution: '+str(SP_res)+'\n')
    f.write('Short term analysis:\n')
    f.write('   Upper limit: '+str(WTD_CI_95_mean_s)+'\n')
    f.write('   WTDopt: ' + str(WTD_CI_95_mean_s) + '\n')
    f.write('   Lower limit: ' + str(WTD_CI_5_mean_s) + '\n\n')
    f.write('Long term analysis:\n')
    f.write('   Upper limit: ' + str(WTD_CI_95_mean_l) + '\n')
    f.write('   WTDopt: ' + str(WTD_CI_95_mean_l) + '\n')
    f.write('   Lower limit: ' + str(WTD_CI_5_mean_l) + '\n')
    f.write('-------------------------------------------------------------\n')
    f.write('------------------------------------------------------------------------------------\n\n')
