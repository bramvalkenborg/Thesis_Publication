from thesis_pub_tools import *
from validation_good_practice.ancillary.metrics import correct_n
from patsy import dmatrices
import statsmodels_adapted.api as sm

input_file = 'CA_MER_GPP_analysis.csv'
input_dir = '/data/leuven/336/vsc33653/Data/GPP_FluxTower/'
output_file = '/data/leuven/336/vsc33653/OUTPUT_pub/GPP_output.txt'
p = 1

# Load the data
Df = pd.read_csv(input_dir+input_file)
Df['Time'] = pd.to_datetime(Df['TIMESTAMP_END'])
del Df['TIMESTAMP_END']
Df = Df.set_index(['Time'])
WTD = pd.Series(Df['WTD'])
WTD = WTD/100
GPP = pd.Series(Df['GPP'])
PAR = pd.Series(Df['PAR'])

# Remove data out of the growing season (= Jun - Sep)
GPP = GPP.drop(GPP.index[GPP.index.month.isin([1, 2, 3, 4, 5, 10, 11, 12])])
PAR = PAR.drop(PAR.index[PAR.index.month.isin([1, 2, 3, 4, 5, 10, 11, 12])])
WTD = WTD.drop(WTD.index[WTD.index.month.isin([1, 2, 3, 4, 5, 10, 11, 12])])

# Average to different resolution, coarser resolution reduces the error
GPP_8D = GPP.resample('8D').mean()
WTD_8D = WTD.resample('8D').mean()
PAR_8D = PAR.resample('8D').mean()

# Normalization
GPPn_8D = GPP_8D/PAR_8D

# Refill the timeseries to a daily resolution, this is needed for the anomaly calculations
GPPn_1D = GPPn_8D.resample('1D').bfill()
WTD_1D = WTD_8D.resample('1D').bfill()

# Calculate the anomalies
WTD_shortAnom = calc_anom(WTD_1D, longterm=False)
GPPn_shortAnom = calc_anom(GPPn_1D, longterm=False)

# Change the resolution of the 1D anomalies to 8D anomalies
WTD_shortAnom8d = WTD_shortAnom.resample('8D').first()
GPPn_shortAnom8d = GPPn_shortAnom.resample('8D').first()


# Model short term anomalies
coef_s = np.zeros(2)
p_values_s = np.zeros(2)
# fp_values_s = np.zeros((len(lat), len(lon)))
# WTD_opt_s = np.zeros((len(lat), len(lon)))
# n_corr_s = np.zeros((len(lat), len(lon)))
# rho_s = np.zeros((len(lat), len(lon)))
# n_s = np.zeros((len(lat), len(lon)))
# Rsq_s = np.zeros((len(lat), len(lon)))
x1 = WTD_shortAnom8d
x2 = WTD_8D
y = GPPn_shortAnom8d
if not np.isnan(x1).all() and not np.isnan(x2).all() and not np.isnan(y).all():
    df = pd.DataFrame(np.transpose(np.asarray([x1, x2, y])), index=x2.index, columns=['WTD_sAnom', 'WTD', 'SIF_sAnom'])
    df.dropna(inplace=True)
    n_corr_s = correct_n(df)
    # n_s[lat1, lon1] = len(np.where(np.logical_not(np.isnan(SIF_sAnom[lat1, lon1, :])))[0])
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
    model = sm.OLS(y, X[:, 1:3]).fit(nobs_corr=n_corr_s)
    # model = sm.OLS(y, X[:, 1:3]).fit()
    # Y_r, X_r = dmatrices('y ~ x1')
    # model_r = sm.OLS(y, X_r[:, 1:3]).fit()
    # fp_values_s[lat1, lon1] = CalcP_f_statistic(model, model_r, n_corr[lat1, lon1], 2)
    fp_values_s = model.f_pvalue
    p_values_s[:] = model.pvalues
    if fp_values_s <= p:
        coef_s[:] = model.params
        WTD_opt_s = -coef_s[0]/coef_s[1]
        Rsq_s = model.rsquared
    else:
        coef_s[:] = model.params
        WTD_opt_s = -coef_s[0]/coef_s[1]
        Rsq_s = model.rsquared
        # p_values_s[lat1, lon1, :] = np.nan
        # n_corr[lat1, lon1] = np.nan
else:
    coef_s[:] = np.nan
    fp_values_s = np.nan
    WTD_opt_s = np.nan
    n_corr_s = np.nan
    n_s = np.nan
    rho_s = np.nan
    Rsq_s = np.nan

print('-------------------------------------------------------------')
print('Short term results')
print('-------------------------------------------------------------')
print('Alpha: '+str(coef_s[0]))
print('Beta: '+str(coef_s[1]))
print('WTDopt: '+str(WTD_opt_s))
print('fp_values_s: '+str(fp_values_s))
print('R2: '+str(Rsq_s))
print('-------------------------------------------------------------')

with open(output_file, 'a') as f:
    f.write('-------------------------------------------------------------\n')
    f.write('Short term results                           Date:'+str(datetime.datetime.today()))
    f.write('-------------------------------------------------------------\n')
    f.write('Alpha: '+str(coef_s[0])+'\n')
    f.write('Beta: '+str(coef_s[1])+'\n')
    f.write('WTDopt: '+str(WTD_opt_s)+'\n')
    f.write('fp_values: '+str(fp_values_s)+'\n')
    f.write('R2: '+str(Rsq_s)+'\n')
    f.write('-------------------------------------------------------------\n')
    f.write('Input file: '+input_dir+input_file+'\n')
    f.write('-------------------------------------------------------------\n')
    f.write('-------------------------------------------------------------\n\n')

GPPn_longAnom = calc_anom(GPPn_1D, longterm=True) - calc_anom(GPPn_1D, longterm=False)
WTD_longAnom = calc_anom(WTD_1D, longterm=True) - calc_anom(WTD_1D, longterm=False)
WTD_longAnom8d = WTD_longAnom.resample('8D').first()
GPPn_longAnom8d = GPPn_longAnom.resample('8D').first()

# Model long term anomalies
coef_l = np.zeros(2)
p_values_l = np.zeros(2)
x1 = WTD_longAnom8d
x2 = WTD_8D
y = GPPn_longAnom8d
x1 = WTD_longAnom8d
x2 = WTD_8D
y = GPPn_longAnom8d
# Model long term anomalies
coef_l = np.zeros(2)
p_values_l = np.zeros(2)
x1 = WTD_longAnom8d
x2 = WTD_8D
y = GPPn_longAnom8d
if not np.isnan(x1).all() and not np.isnan(x2).all() and not np.isnan(y).all():
    df = pd.DataFrame(np.transpose(np.asarray([x1, x2, y])), index=x2.index, columns=['WTD_lAnom', 'WTD', 'SIF_lAnom'])
    df.dropna(inplace=True)
    n_corr_l = correct_n(df)
    # n_s[lat1, lon1] = len(np.where(np.logical_not(np.isnan(SIF_sAnom[lat1, lon1, :])))[0])
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
    model = sm.OLS(y, X[:, 1:3]).fit(nobs_corr=n_corr_l)
    # model = sm.OLS(y, X[:, 1:3]).fit()
    # Y_r, X_r = dmatrices('y ~ x1')
    # model_r = sm.OLS(y, X_r[:, 1:3]).fit()
    # fp_values_s[lat1, lon1] = CalcP_f_statistic(model, model_r, n_corr[lat1, lon1], 2)
    fp_values_l = model.f_pvalue
    p_values_l[:] = model.pvalues
    if fp_values_l <= p:
        coef_l[:] = model.params
        WTD_opt_l = -coef_l[0]/coef_l[1]
        Rsq_l = model.rsquared
    else:
        coef_l[:] = model.params
        WTD_opt_l = -coef_l[0]/coef_l[1]
        Rsq_l = model.rsquared
else:
    coef_l[:] = np.nan
    fp_values_l = np.nan
    WTD_opt_l = np.nan
    n_corr_l = np.nan
    n_l = np.nan
    rho_l = np.nan
    Rsq_l = np.nan

print('-------------------------------------------------------------')
print('Long term results')
print('-------------------------------------------------------------')
print('Alpha: '+str(coef_l[0]))
print('Beta: '+str(coef_l[1]))
print('WTDopt: '+str(WTD_opt_l))
print('fp_values_s: '+str(fp_values_l))
print('R2: '+str(Rsq_l))
print('-------------------------------------------------------------')

with open(output_file, 'a') as f:
    f.write('-------------------------------------------------------------\n')
    f.write('Short term results                           Date:'+str(datetime.datetime.today()))
    f.write('-------------------------------------------------------------\n')
    f.write('Alpha: '+str(coef_s[0])+'\n')
    f.write('Beta: '+str(coef_s[1])+'\n')
    f.write('WTDopt: '+str(WTD_opt_s)+'\n')
    f.write('fp_values: '+str(fp_values_s)+'\n')
    f.write('R2: '+str(Rsq_s)+'\n')
    f.write('-------------------------------------------------------------\n')
    f.write('Input file: '+input_dir+input_file+'\n\n')
    f.write('-------------------------------------------------------------\n')
    f.write('-------------------------------------------------------------\n\n')

'''
# ----------------------------------------------------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------------------------------------------------
slope_lt = np.zeros(len(GPPn_longAnom8d))
slope_st = np.zeros(len(GPPn_longAnom8d))
for time1 in range(len(GPPn_longAnom8d)):
    slope_lt[time1] = coef_l[1]*WTD[time1]+coef_l[0]
    slope_st[time1] = coef_s[1] * WTD[time1] + coef_s[0]

slope_st_pixel = slope_st[:]
slope_lt_pixel = slope_lt[:]
WTD_opt_st_pixel = WTD_opt_s
WTD_opt_lt_pixel = WTD_opt_l
x1_st = coef_s[0]
x2_st = coef_s[1]
x1_lt = coef_l[0]
x2_lt = coef_l[1]
SIF_sAnom_wet = np.where(WTD > WTD_opt_s, SIF_sAnom[:], np.nan)
SIF_sAnom_dry = np.where(WTD < WTD_opt_s, SIF_sAnom[:], np.nan)
SIF_lAnom_wet = np.where(WTD > WTD_opt_l, SIF_lAnom[:], np.nan)
SIF_lAnom_dry = np.where(WTD < WTD_opt_l, SIF_lAnom[:], np.nan)
WTD_sAnom_wet = np.where(WTD > WTD_opt_s, WTD_sAnom[:], np.nan)
WTD_sAnom_dry = np.where(WTD < WTD_opt_s, WTD_sAnom[:], np.nan)
WTD_lAnom_wet = np.where(WTD > WTD_opt_l, WTD_lAnom[:], np.nan)
WTD_lAnom_dry = np.where(WTD < WTD_opt_l, WTD_lAnom[:], np.nan)
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
ax5a.set_title('(c)', loc='left', fontsize=10)
ax6a.set_title('(f)', loc='left', fontsize=10)
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
'''




