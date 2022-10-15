from thesis_pub_tools import *
from validation_good_practice.ancillary.metrics import correct_n
from patsy import dmatrices
import statsmodels_adapted.api as sm

input_file = 'CA_MER_GPP_analysis.csv'
input_dir = '/data/leuven/336/vsc33653/Data/GPP_FluxTower/'
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

GPP_2001 = GPP[GPP.index.year == 2001]
GPP_2004 = GPP[GPP.index.year == 2004]

GPP_2001 = GPP_2001.resample('H').mean()
GPP_2004 = GPP_2004.resample('H').mean()

mask_2001 = (GPP_2001.index.day_of_year > 252) & (GPP_2001.index.day_of_year <= 264)
GPP_2001_DOY = GPP_2001.loc[mask_2001]
mask_2004 = (GPP_2004.index.day_of_year > 252) & (GPP_2004.index.day_of_year <= 264)
GPP_2004_DOY = GPP_2004.loc[mask_2004]

plt.figure(figsize=(8, 4))
plt.plot(GPP_2001_DOY.index, GPP_2001_DOY)
GPP_2001_DOY_ticks = GPP_2001_DOY.resample('D').mean()
plt.xticks(ticks=np.array(GPP_2001_DOY_ticks.index), labels=np.array(GPP_2001_DOY_ticks.index.day_of_year))
plt.ylim((0, 14))
plt.grid()
plt.tight_layout()

plt.figure(figsize=(8, 4))
plt.plot(GPP_2004_DOY.index, GPP_2004_DOY)
GPP_2004_DOY_ticks = GPP_2004_DOY.resample('D').mean()
plt.xticks(ticks=np.array(GPP_2004_DOY_ticks.index), labels=np.array(GPP_2004_DOY_ticks.index.day_of_year))
plt.ylim((0, 14))
plt.grid()
plt.tight_layout()
