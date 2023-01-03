import pandas as pd
import numpy as np
from thesis_pub_tools import calc_anom, cal_WaterStressModel
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from GPP_datasets import load_dataset
import copy

# ----------------------------------------------------------------------------------------------------------------------
# INPUT
# ----------------------------------------------------------------------------------------------------------------------
source_path = '/data/leuven/336/vsc33653/'
source_path = '/data/leuven/317/vsc31786/'
#source_path = '/Users/bramvalkenborg/Library/CloudStorage/OneDrive-KULeuven/Thesis/Publication/Resubmission/'

# Select the right input file, output directory and output file
#input_file = 'CA_MER_GPP_analysis.csv'
#input_file = 'US-Los_HH_2000-2022.csv'
#input_file = 'FI_LOM_GPP.csv'
input_file = 'mb_met_flux_data_1998_2018.txt'
#input_file = 'FLX_SE-Deg_FLUXNET2015_FULLSET_HH_2001-2020_beta-3.csv'
#input_file = 'SFZ_S_2010_2020_hh_Michl.csv'

# Option to write the data to a separate file
write = False
description = 'Resolution 8D, GPP/PAR'

input_dir = source_path + 'Data/GPP_FluxTower/'
input_dir = source_path + 'peatland_data/GPP_FluxTower/'

output_file = source_path + 'OUTPUT_pub/GPP_output.txt'
output_file = source_path + 'GPP/GPP_output.txt'
output_file = '/scratch/leuven/317/vsc31786/GPP_output.txt'
output_path_figure = '/staging/leuven/stg_00024/OUTPUT/michelb/FIG_tmp/GPP_analysis/'

# Set the input parameters
p = 1
# new cloud filter definition --> 0.5 means PAR needs to be 50% of max PAR in that month
cloud_filter = 0.50
description = description + ', cloud filter: ' + str(cloud_filter)
# ----------------------------------------------------------------------------------------------------------------------

# Load the column
WTD, PAR, GPP = load_dataset(input_dir, input_file)

# Remove data out of the growing season (= Jun - Sep)
GPP = GPP.drop(GPP.index[GPP.index.month.isin([1, 2, 3, 4, 5, 10, 11, 12])])
PAR = PAR.drop(PAR.index[PAR.index.month.isin([1, 2, 3, 4, 5, 10, 11, 12])])
WTD = WTD.drop(WTD.index[WTD.index.month.isin([1, 2, 3, 4, 5, 10, 11, 12])])

# Focus on noon, around SIF measurements
# revised: centered around hour of highest GPP value to catch possible inaccuracies due to time zone issues
hours = PAR.index.hour
df_temp = pd.DataFrame({'GPP':GPP,'hour':hours})
#df_temp.boxplot(column='GPP',by='hour')
# determine hour with max value of GPP
maxvalue = df_temp.groupby('hour').quantile(0.5).max()
GPPmax_hour = np.where(df_temp.groupby('hour').quantile(0.5)['GPP'].values==maxvalue.values[0])[0][0]
cond_hours = (hours >= GPPmax_hour-1) & (hours <= GPPmax_hour+1)

GPP = GPP[cond_hours]
PAR = PAR[cond_hours]
WTD = WTD[cond_hours]

# plot time series for checking
#plt.subplot(3,1,1)
#plt.plot(GPP,'.')
#plt.subplot(3,1,2)
#plt.plot(WTD,'.')
#plt.subplot(3,1,3)
#plt.plot(PAR,'.')
#plt.close()

# Cloud filter: remove the lowest PAR data of every month
PAR_cloud_filtered = copy.deepcopy(PAR)
GPP_cloud_filtered = copy.deepcopy(GPP)
WTD_cloud_filtered = copy.deepcopy(WTD)
for m in [6, 7, 8, 9]:
    months = PAR.index.month
    PAR_m = np.where(months == m, PAR, np.nan)
    PAR_max = np.nanmax(PAR_m)
    #PAR_cloud = PAR[PAR > np.nanquantile(PAR_m, cloud_filter)]
    #GPP = GPP[PAR > np.nanquantile(PAR_m, cloud_filter)]
    #WTD = WTD[PAR > np.nanquantile(PAR_m, cloud_filter)]
    #PAR = PAR_cloud
    PAR_cloud_filtered[(PAR < cloud_filter*PAR_max) & (months == m)] = np.nan
    GPP_cloud_filtered[(PAR < cloud_filter*PAR_max) & (months == m)] = np.nan
    WTD_cloud_filtered[(PAR < cloud_filter*PAR_max) & (months == m)] = np.nan

PAR = PAR_cloud_filtered
GPP = GPP_cloud_filtered
WTD = WTD_cloud_filtered

# catch outliers due to problems with PAR sensor
# Approach: Remove GPP/PAR values that deviate by more than 30% from the median of the current noon window
GPPn = GPP/PAR
df_temp = pd.DataFrame({'GPPn':GPPn})
GPPn_rolling = df_temp.rolling(10,min_periods=2).median()['GPPn']
rel_dif = np.abs(GPPn-GPPn_rolling)/(0.5*(GPPn+GPPn_rolling))
filt_outlier = rel_dif>0.3
PAR[filt_outlier] = np.nan
GPP[filt_outlier] = np.nan
WTD[filt_outlier] = np.nan
############################################

cond_nan = np.isnan(WTD) | np.isnan(GPP) | np.isnan(PAR)
WTD[cond_nan.values] = np.nan
GPP[cond_nan.values] = np.nan
PAR[cond_nan.values] = np.nan

# Average to different resolution, coarser resolution reduces the error
GPP_8D = GPP.resample('8D').mean()
count8D = GPP.resample('8D').count()
WTD_8D = WTD.resample('8D').mean()
PAR_8D = PAR.resample('8D').mean()
GPP_8D[count8D < 0.2*np.max(count8D)] = np.nan
WTD_8D[count8D < 0.2*np.max(count8D)] = np.nan
PAR_8D[count8D < 0.2*np.max(count8D)] = np.nan

# Normalization
GPPn_8D = GPP_8D/PAR_8D

# Refill the timeseries to a daily resolution, this is needed for the anomaly calculations
# GPPn_1D = GPPn_8D.resample('1D').bfill()
# WTD_1D = WTD_8D.resample('1D').bfill()
# PAR_1D = PAR_8D.resample('1D').bfill()
GPPn_1D = GPPn_8D.resample('1D').ffill(limit=7)     # set fill limit to avoid filling of long gaps
WTD_1D = WTD_8D.resample('1D').ffill(limit=7)
PAR_1D = PAR_8D.resample('1D').ffill(limit=7)

# Calculate the anomalies
WTD_shortAnom = calc_anom(WTD_1D, longterm=False)
GPPn_shortAnom = calc_anom(GPPn_1D, longterm=False)
PAR_shortAnom = calc_anom(PAR_1D, longterm=False)

# Change the resolution of the 1D anomalies to 8D anomalies
#WTD_shortAnom8d = WTD_shortAnom.resample('8D').first()
#GPPn_shortAnom8d = GPPn_shortAnom.resample('8D').first()
#PAR_shortAnom8d = PAR_shortAnom.resample('8D').first()
WTD_shortAnom8d = WTD_shortAnom.resample('8D').mean()
GPPn_shortAnom8d = GPPn_shortAnom.resample('8D').mean()
PAR_shortAnom8d = PAR_shortAnom.resample('8D').mean()

# Model short term anomalies
coef_s, WTD_opt_s, n_corr_s, fp_values_s, p_values_s, Rsq_s = cal_WaterStressModel(GPPn_shortAnom8d, WTD_shortAnom8d, WTD_8D, 1, WTD_8D.index)

print('-------------------------------------------------------------')
print('Short-term results')
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
WTD_longAnom8d = WTD_longAnom.resample('8D').mean()
PAR_longAnom8d = PAR_longAnom.resample('8D').mean()
GPPn_longAnom8d = GPPn_longAnom.resample('8D').mean()

# Model long term anomalies
coef_l, WTD_opt_l, n_corr_l, fp_values_l, p_values_l, Rsq_l = cal_WaterStressModel(GPPn_longAnom8d, WTD_longAnom8d, WTD_8D, 1, WTD_8D.index)

# GPPn anomalies vs WTD
dpi = 300
fig = plt.figure(figsize = (12,5))
ax1 = plt.subplot(1,2,1)
p1=plt.plot(WTD_8D,GPPn_shortAnom8d,'.')
GPPn_shortAnom8d_mod = (coef_s[0]+coef_s[1]*WTD_8D)*WTD_shortAnom8d
p2=plt.plot(WTD_8D,GPPn_shortAnom8d_mod,'r.')
plt.ticklabel_format(axis='y',style='sci',scilimits=(1,4))
plt.xlabel('Water level (m)')
plt.ylabel('GPP/PAR short-term anomaly (-)')
plt.legend(['obs','mod'])
plt.text(0.5, 0.01, 'a: %2.4f b: %2.4f WTDopt: %2.2f' % (coef_s[0],coef_s[1],WTD_opt_s), horizontalalignment='center',transform=ax1.transAxes)
ax2 = plt.subplot(1,2,2)
plt.plot(WTD_8D,GPPn_longAnom8d,'.')
GPPn_longAnom8d_mod = (coef_l[0]+coef_l[1]*WTD_8D)*WTD_longAnom8d
plt.plot(WTD_8D,GPPn_longAnom8d_mod,'r.')
plt.ticklabel_format(axis='y',style='sci',scilimits=(1,4))
plt.xlabel('Water level (m)')
plt.ylabel('GPP/PAR long-term anomaly (-)')
plt.legend(['obs','mod'])
plt.text(0.5, 0.01, 'a: %2.4f b: %2.4f WTDopt: %2.2f' % (coef_l[0],coef_l[1],WTD_opt_l), horizontalalignment='center',transform=ax2.transAxes)
plt.savefig(output_path_figure+input_file+'.png',dpi=dpi)
plt.close()


# GPPn anomalies vs WTD
dpi = 300
fig = plt.figure(figsize = (12,5))
ax1 = plt.subplot(1,2,1)
p1=plt.plot(WTD_shortAnom8d,GPPn_shortAnom8d,'.')
GPPn_shortAnom8d_mod = (coef_s[0]+coef_s[1]*WTD_8D)*WTD_shortAnom8d
p2=plt.plot(WTD_shortAnom8d,GPPn_shortAnom8d_mod,'r.')
plt.ticklabel_format(axis='y',style='sci',scilimits=(1,4))
plt.xlabel('Water level short-term anomaly (m)')
plt.ylabel('GPP/PAR short-term anomaly (-)')
plt.legend(['obs','mod'])
plt.text(0.5, 0.01, 'a: %2.4f b: %2.4f WTDopt: %2.2f' % (coef_s[0],coef_s[1],WTD_opt_s), horizontalalignment='center',transform=ax1.transAxes)
ax2 = plt.subplot(1,2,2)
plt.plot(WTD_longAnom8d,GPPn_longAnom8d,'.')
GPPn_longAnom8d_mod = (coef_l[0]+coef_l[1]*WTD_8D)*WTD_longAnom8d
plt.plot(WTD_longAnom8d,GPPn_longAnom8d_mod,'r.')
plt.ticklabel_format(axis='y',style='sci',scilimits=(1,4))
plt.xlabel('Water level short-term anomaly (m)')
plt.ylabel('GPP/PAR long-term anomaly (-)')
plt.legend(['obs','mod'])
plt.text(0.5, 0.01, 'a: %2.4f b: %2.4f WTDopt: %2.2f' % (coef_l[0],coef_l[1],WTD_opt_l), horizontalalignment='center',transform=ax2.transAxes)
plt.savefig(output_path_figure+input_file+'_anomalies.png',dpi=dpi)
plt.close()


print('-------------------------------------------------------------')
print('Long term results')
print('-------------------------------------------------------------')
print('Alpha: '+str(coef_l[0]))
print('Beta: '+str(coef_l[1]))
print('WTDopt: '+str(WTD_opt_l))
print('fp_values_s: '+str(fp_values_l))
print('R2: '+str(Rsq_l))
print('-------------------------------------------------------------')
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
fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(10, 7))
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
fig.tight_layout()
if write:
    plt.savefig(source_path + 'OUTPUT_pub/'+input_file[:-4]+'_Timeseries.png')

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
