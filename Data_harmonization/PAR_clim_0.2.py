from SIF_tools.python.L2_tools import convert_time
import numpy as np
import datetime
from matplotlib import dates
from netCDF4 import Dataset
import cftime
from pytesmo.time_series.anomaly import calc_climatology
import pandas as pd

# Load all the data
SIF_nc = Dataset('/data/leuven/336/vsc33653/OUTPUT/Gridding/NH_SIF_2021_dN_omitnan_1day_0.2deg_c30.nc')
time_SIF = SIF_nc.variables['time'][:]
time_SIFgreg = convert_time(time_SIF, 'days since 1970-01-01', u"gregorian")
month_SIF = time_SIFgreg.astype('datetime64[M]').astype(int) % 12 + 1
lat_SIF = SIF_nc.variables['lat'][:]
lon_SIF = SIF_nc.variables['lon'][:]
# Load the peatmap, to save calculation time, only peatpixels are interpolated
Peat_nc = Dataset('/data/leuven/336/vsc33653/OUTPUT/Gridding/NH_PeatMap_0.2deg.nc')
peat = Peat_nc.variables['peatmap'][:, :, 0]
peat_mask = (peat != 0)

# Load all the data
PAR0a_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                   'CERES_0-2deg_int_t0_t125_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-20210731.nc')
PAR0b_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                   'CERES_0-2deg_int_t125_t225_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-20210731.nc')
PAR0c_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                   'CERES_0-2deg_int_t225_t325_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-20210731.nc')
PAR1_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t325_t425_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-20210731.nc')
PAR2_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t425_t525_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-20210731.nc')
PAR3_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t525_t625_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-20210731.nc')
PAR4_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t625_t725_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-20210731.nc')
PAR5_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t725_t825_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-20210731.nc')
PAR6_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t825_t930_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-20210731.nc')
PAR7_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t930_t1099_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-20210731.nc')
PAR8_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/CERES_0-2deg_int_t1099_t1340_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-20210731.nc', mode='r+')
PAR_int = np.zeros((len(lat_SIF), len(lon_SIF), len(time_SIF)))
PAR_int[:, :, 0:125] = PAR0a_nc.variables['PAR_int'][:, :, 0:125]
PAR_int[:, :, 125:225] = PAR0b_nc.variables['PAR_int'][:, :, 125:225]
PAR_int[:, :, 225:325] = PAR0c_nc.variables['PAR_int'][:, :, 225:325]
PAR_int[:, :, 325:425] = PAR1_nc.variables['PAR_int'][:, :, 325:425]
PAR_int[:, :, 425:525] = PAR2_nc.variables['PAR_int'][:, :, 425:525]
PAR_int[:, :, 525:625] = PAR3_nc.variables['PAR_int'][:, :, 525:625]
PAR_int[:, :, 625:725] = PAR4_nc.variables['PAR_int'][:, :, 625:725]
PAR_int[:, :, 725:825] = PAR5_nc.variables['PAR_int'][:, :, 725:825]
PAR_int[:, :, 825:930] = PAR6_nc.variables['PAR_int'][:, :, 825:930]
PAR_int[:, :, 930:1099] = PAR7_nc.variables['PAR_int'][:, :, 930:1099]
PAR_int[:, :, 1099:] = PAR8_nc.variables['PAR_int'][:, :, 1099:]
PAR_int = np.where(PAR_int == -999, np.nan, PAR_int)
PAR_int[:, :, 1341:len(time_SIF)] = np.nan

PAR_clim = np.zeros((len(lat_SIF), len(lon_SIF), len(time_SIF)))
for lat1 in range(len(lat_SIF)):
    for lon1 in range(len(lon_SIF)):
        if peat_mask[lat1, lon1]:
            PAR_pix = PAR_int[lat1, lon1, :]
            if not np.isnan(PAR_pix).all():
                PAR_ser = pd.Series(PAR_pix, index=time_SIFgreg)
                PAR_ser = PAR_ser.resample('8D').mean()
                PAR_ser = PAR_ser.resample('1D').bfill()
                PAR_clim_pix = calc_climatology(PAR_ser)
                # PAR_clim8d = PAR_clim_pix.resample('8D').first()
                # PAR_clim8d = np.where(np.isnan(PAR_pix), np.nan, PAR_clim8d)
                # PAR_clim[lat1, lon1, :] = PAR_clim_pix
                # PAR_int[lat1, lon1, 1341:len(time_SIF)] = np.array(PAR_clim_pix[212:273])
                PAR8_nc['PAR_int'][lat1, lon1, 1341:len(time_SIF)] = np.array(PAR_clim_pix[212:273])
                print('lon: '+str(lon1)+', lat: '+str(lat1))
            else:
                PAR_clim[lat1, lon1, 1341:len(time_SIF)] = np.nan
                PAR8_nc['PAR_int'][lat1, lon1, 1341:len(time_SIF)] = np.nan
                # print('lon: ' + str(lon1) + ', lat: ' + str(lat1))
        else:
            PAR_clim[lat1, lon1, 1341:len(time_SIF)] = np.nan
            PAR8_nc['PAR_int'][lat1, lon1, 1341:len(time_SIF)] = np.nan
            # print('lon: ' + str(lon1) + ', lat: ' + str(lat1))
PAR8_nc.close()
