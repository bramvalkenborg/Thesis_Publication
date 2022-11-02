from thesis_pub_tools import *
from validation_good_practice.ancillary.metrics import correct_n
from patsy import dmatrices
import statsmodels_adapted.api as sm

write = False
description = 'Resolution 8D, GPP'

input_file = 'EFDC_L2_Flx_BEJal_2006_v04_30m.txt'
input_dir = '/data/leuven/336/vsc33653/Data/GPP_FluxTower/'
output_file = '/data/leuven/336/vsc33653/OUTPUT_pub/GPP_output.txt'
p = 1

# Load the data
Df = pd.read_csv(input_dir+input_file, na_values=-9999)
Df['Time'] = pd.to_datetime(Df['TIMESTAMP_END'], format='%Y%m%d%H%M')
del Df['TIMESTAMP_END']
del Df['TIMESTAMP_START']
del Df['DTime']

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
GPPn_8D = GPP_8D

# Refill the timeseries to a daily resolution, this is needed for the anomaly calculations
GPPn_1D = GPPn_8D.resample('1D').bfill()
WTD_1D = WTD_8D.resample('1D').bfill()
PAR_1D = PAR_8D.resample('1D').bfill()

# Calculate the anomalies
WTD_shortAnom = calc_anom(WTD_1D, longterm=False)
GPPn_shortAnom = calc_anom(GPPn_1D, longterm=False)
PAR_shortAnom = calc_anom(PAR_1D, longterm=False)

# Change the resolution of the 1D anomalies to 8D anomalies
WTD_shortAnom8d = WTD_shortAnom.resample('8D').first()
GPPn_shortAnom8d = GPPn_shortAnom.resample('8D').first()
PAR_shortAnom8d = PAR_shortAnom.resample('8D').first()

# Model short term anomalies
coef_s, WTD_opt_s, n_corr_s, fp_values_s, p_values_s, Rsq_s = cal_WaterStressModel(GPPn_shortAnom8d, WTD_shortAnom8d, WTD_8D, WTD_8D.index)


print('-------------------------------------------------------------')
print('Short term results')
print('-------------------------------------------------------------')
print('Alpha: '+str(coef_s[0]))
print('Beta: '+str(coef_s[1]))
print('WTDopt: '+str(WTD_opt_s))
print('fp_values_s: '+str(fp_values_s))
print('R2: '+str(Rsq_s))
print('-------------------------------------------------------------')

if write:
    with open(output_file, 'a') as f:
        f.write('-------------------------------------------------------------\n')
        f.write('Short term results                           Date:'+str(datetime.datetime.today())+'\n')
        f.write('-------------------------------------------------------------\n')
        f.write(description+'\n')
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
PAR_longAnom = calc_anom(PAR_1D, longterm=True) - calc_anom(PAR_1D, longterm=False)
WTD_longAnom8d = WTD_longAnom.resample('8D').first()
PAR_longAnom8d = PAR_longAnom.resample('8D').first()
GPPn_longAnom8d = GPPn_longAnom.resample('8D').first()

# Model long term anomalies
coef_l, WTD_opt_l, n_corr_l, fp_values_l, p_values_l, Rsq_l = cal_WaterStressModel(GPPn_longAnom8d, WTD_longAnom8d, WTD_8D, WTD_8D.index)


print('-------------------------------------------------------------')
print('Long term results')
print('-------------------------------------------------------------')
print('Alpha: '+str(coef_l[0]))
print('Beta: '+str(coef_l[1]))
print('WTDopt: '+str(WTD_opt_l))
print('fp_values_s: '+str(fp_values_l))
print('R2: '+str(Rsq_l))
print('-------------------------------------------------------------')

if write:
    with open(output_file, 'a') as f:
        f.write('-------------------------------------------------------------\n')
        f.write('Long term results                           Date:'+str(datetime.datetime.today())+'\n')
        f.write('-------------------------------------------------------------\n')
        f.write(description+'\n')
        f.write('-------------------------------------------------------------\n')
        f.write('Alpha: '+str(coef_l[0])+'\n')
        f.write('Beta: '+str(coef_l[1])+'\n')
        f.write('WTDopt: '+str(WTD_opt_l)+'\n')
        f.write('fp_values: '+str(fp_values_l)+'\n')
        f.write('R2: '+str(Rsq_l)+'\n')
        f.write('-------------------------------------------------------------\n')
        f.write('Input file: '+input_dir+input_file+'\n')
        f.write('-------------------------------------------------------------\n')
        f.write('-------------------------------------------------------------\n\n')


# Calculation of the seasonality
def seas(x_ser, time):
    x_ser = x_ser.resample('1D').bfill()
    x_seas = x_ser - calc_anom(x_ser, longterm=False)
    x_seas8d = np.zeros(len(time))
    x_seas8d = x_seas.resample('8D').first()
    return x_seas8d


def clim(x_ser, time):
    x_ser = x_ser.resample('1D').bfill()
    x_clim1 = x_ser - calc_anom(x_ser, longterm=True)
    x_clim8d = np.zeros(len(time))
    x_clim8d = x_clim1.resample('8D').first()
    return x_clim8d


WTD_seas_8D = seas(WTD_8D, time=WTD_8D.index)
WTD_clim_8D = clim(WTD_8D, time=WTD_8D.index)
GPPn_seas_8D = seas(GPPn_8D, time=GPPn_8D.index)
GPPn_clim_8D = clim(GPPn_8D, time=GPPn_8D.index)
GPP_seas_8D = seas(GPP_8D, time=GPPn_8D.index)
GPP_clim_8D = clim(GPP_8D, time=GPPn_8D.index)

# Timeseries
fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1)
GPP_8D.plot(ax=ax1, color='red')
GPP_seas_8D.plot(ax=ax1, color='grey')
GPP_clim_8D.plot(ax=ax1, color='grey', linestyle='--')
PAR_8D.plot(ax=ax2, color='green')
GPPn_8D.plot(ax=ax3, color='orange')
GPPn_seas_8D.plot(ax=ax3, color='grey')
GPPn_clim_8D.plot(ax=ax3, color='grey', linestyle='--')
WTD_8D.plot(ax=ax4, color='blue')
WTD_seas_8D.plot(ax=ax4, color='grey')
WTD_clim_8D.plot(ax=ax4, color='grey', linestyle='--')
plt.hlines(WTD_opt_s, color='grey', linestyles=':', xmin=WTD_8D.index.min(), xmax=WTD_8D.index.max())
plt.hlines(WTD_opt_l, color='slategray', linestyles='-.', xmin=WTD_8D.index.min(), xmax=WTD_8D.index.max())
ax1.set_ylabel('GPP')
ax2.set_ylabel('PAR')
ax3.set_ylabel('GPP/PAR')
ax4.set_ylabel('WTD')


plt.figure()
plt.scatter(GPP_8D, PAR_8D)
plt.xlabel('GPP')
plt.ylabel('PAR')

plt.figure()
plt.scatter(PAR, GPP, s=0.5, c=GPP.index.month, cmap='brg')
plt.xlabel('PAR')
plt.ylabel('GPP')
plt.colorbar()

plt.figure()
plt.scatter(PAR_shortAnom8d, GPPn_shortAnom8d)
plt.xlabel('PAR_shortAnom')
plt.ylabel('GPP_shortAnom')
plt.xlim((-600, 600))
plt.ylim((-4, 4))
plt.grid()

plt.figure()
plt.scatter(PAR_longAnom8d, GPPn_longAnom8d)
plt.xlabel('PAR_longAnom')
plt.ylabel('GPP_longAnom')
plt.grid()
plt.xlim((-600, 600))
plt.ylim((-4, 4))

scenario_1_s = ((WTD_8D.dropna() > WTD_opt_s) & (WTD_shortAnom8d.dropna() > 0) & (GPPn_shortAnom8d.dropna() < 0))
scenario_2_s = ((WTD_8D.dropna() > WTD_opt_s) & (WTD_shortAnom8d.dropna() < 0) & (GPPn_shortAnom8d.dropna() > 0))
scenario_3_s = ((WTD_8D.dropna() < WTD_opt_s) & (WTD_shortAnom8d.dropna() > 0) & (GPPn_shortAnom8d.dropna() > 0))
scenario_4_s = ((WTD_8D.dropna() < WTD_opt_s) & (WTD_shortAnom8d.dropna() < 0) & (GPPn_shortAnom8d.dropna() < 0))
hypothesis_s = scenario_1_s | scenario_2_s | scenario_3_s | scenario_4_s
count_s = hypothesis_s.value_counts()
percentage_True_s = count_s[True] / (count_s[True] + count_s[False]) * 100

scenario_1_l = ((WTD_8D.dropna() > WTD_opt_l) & (WTD_longAnom8d.dropna() > 0) & (GPPn_longAnom8d.dropna() < 0))
scenario_2_l = ((WTD_8D.dropna() > WTD_opt_l) & (WTD_longAnom8d.dropna() < 0) & (GPPn_longAnom8d.dropna() > 0))
scenario_3_l = ((WTD_8D.dropna() < WTD_opt_l) & (WTD_longAnom8d.dropna() > 0) & (GPPn_longAnom8d.dropna() > 0))
scenario_4_l = ((WTD_8D.dropna() < WTD_opt_l) & (WTD_longAnom8d.dropna() < 0) & (GPPn_longAnom8d.dropna() < 0))
hypothesis_l = scenario_1_l | scenario_2_l | scenario_3_l | scenario_4_l
count_l = hypothesis_l.value_counts()
percentage_True_l = count_l[True] / (count_l[True] + count_l[False]) * 100

GPPn_shortAnom8d_True = GPPn_shortAnom8d.where(hypothesis_s==True, np.nan)
GPPn_shortAnom8d_False = GPPn_shortAnom8d.where(hypothesis_s==False, np.nan)
WTD_shortAnom8d_True = WTD_shortAnom8d.where(hypothesis_s==True, np.nan)
WTD_shortAnom8d_False = WTD_shortAnom8d.where(hypothesis_s==False, np.nan)

GPPn_longAnom8d_True = GPPn_longAnom8d.where(hypothesis_l==True, np.nan)
GPPn_longAnom8d_False = GPPn_longAnom8d.where(hypothesis_l==False, np.nan)
WTD_longAnom8d_True = WTD_longAnom8d.where(hypothesis_l==True, np.nan)
WTD_longAnom8d_False = WTD_longAnom8d.where(hypothesis_l==False, np.nan)

fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, sharex=True)
GPPn_shortAnom8d_True.plot.bar(ax=ax1, xticks=[], color='red')
WTD_shortAnom8d_True.plot.bar(ax=ax2, xticks=[], color='blue')
GPPn_shortAnom8d_False.plot.bar(ax=ax1, xticks=[], color='grey')
WTD_shortAnom8d_False.plot.bar(ax=ax2, xticks=[], color='grey')

GPPn_longAnom8d_True.plot.bar(ax=ax3, xticks=[], color='red')
WTD_longAnom8d_True.plot.bar(ax=ax4, xticks=[], color='blue')
GPPn_longAnom8d_False.plot.bar(ax=ax3, xticks=[], color='grey')
WTD_longAnom8d_False.plot.bar(ax=ax4, xticks=[], color='grey')

ax1.set_ylabel('GPPn_short')
ax2.set_ylabel('WTD_short')
ax3.set_ylabel('GPPn_long')
ax4.set_ylabel('WTD_long')

fig, (ax1, ax2) = plt.subplots(1, 2)
norm_st = mpl.colors.TwoSlopeNorm(vmin=np.nanmin(WTD_8D), vcenter=WTD_opt_s,
                                  vmax=np.nanmax(WTD_8D))
WTD_shortAnom8d_w = WTD_shortAnom8d[WTD_8D <= WTD_opt_s]
WTD_shortAnom8d_d = WTD_shortAnom8d[WTD_8D > WTD_opt_s]
GPPn_shortAnom8d_w = GPPn_shortAnom8d[WTD_8D <= WTD_opt_s]
GPPn_shortAnom8d_d = GPPn_shortAnom8d[WTD_8D > WTD_opt_s]
WTD_8D_w = WTD_8D[WTD_8D <= WTD_opt_s]
WTD_8D_d = WTD_8D[WTD_8D > WTD_opt_s]
ax1.scatter(WTD_shortAnom8d_w, GPPn_shortAnom8d_w, c=WTD_8D_w, cmap='seismic', norm=norm_st)
ax2.scatter(WTD_shortAnom8d_d, GPPn_shortAnom8d_d, c=WTD_8D_d, cmap='seismic', norm=norm_st)
ax1.grid()
ax2.grid()
