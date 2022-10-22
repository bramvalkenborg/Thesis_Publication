from thesis_pub_tools import *
from validation_good_practice.ancillary.metrics import correct_n
from patsy import dmatrices
import statsmodels_adapted.api as sm

write = True
description = 'Resolution 8D'
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
WTD_longAnom8d = WTD_longAnom.resample('8D').first()
GPPn_longAnom8d = GPPn_longAnom.resample('8D').first()

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




