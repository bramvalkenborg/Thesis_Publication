"""
------------------------------------------------------------------------------------------------------------------------
Normalises SIF measurements by 8-day averaged PARob (CERES). Ran after the interpolation step.
->  PARob data is split into different files for different time periods
------------------------------------------------------------------------------------------------------------------------
"""

from netCDF4 import Dataset
import numpy as np
from SIF_tools.python.L2_tools import convert_time
from pytesmo.time_series.anomaly import calc_climatology

# Load all the data
PAR0a_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                   'CERES_0-2deg_int_t0_t125_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
PAR0b_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                   'CERES_0-2deg_int_t125_t225_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
PAR0c_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                   'CERES_0-2deg_int_t225_t325_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
PAR1_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t325_t425_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
PAR2_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t425_t525_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
PAR3_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t525_t625_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
PAR4_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t625_t725_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
PAR5_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t725_t825_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
PAR6_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t825_t930_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
PAR7_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t930_t1099_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
PAR8_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t1099_t1340_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
PAR8_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/'
                  'CERES_0-2deg_int_t1340_t1401_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc')
# SIF_nc = Dataset('/data/leuven/336/vsc33653/OUTPUT/Gridding/NH_SIF_2021_dN_omitnan_8day_0.2deg_c30.nc', mode='r+')
SIF_nc = Dataset('/data/leuven/336/vsc33653/OUTPUT/Gridding/NH_SIF_2021_dN_omitnan_8day_0.2deg_c30.nc', mode='r+')
try:
    SIF_nc.createVariable('sif_norm', 'f4', ('lat', 'lon', 'time'), fill_value=-999)
except:
    print('Variable sif_norm already exists')

lat1 = SIF_nc.variables['lat'][:]
lon1 = SIF_nc.variables['lon'][:]
time_SIF = SIF_nc.variables['time'][:]
time_SIFgreg = convert_time(time_SIF, 'days since 1970-01-01', u"gregorian")
time_PAR = PAR2_nc.variables['time'][:]

# Create a new PAR netCDF file containing the 8-day averaged PAR
PAR_nc = Dataset('/data/leuven/336/vsc33653/OUTPUT/Gridding/NH_PAR_2021_8day_0.2deg.nc', mode='w', format='NETCDF4')
PAR_nc.createDimension('lat', len(lat1))
PAR_nc.createDimension('lon', len(lon1))
PAR_nc.createDimension('time', len(time_SIF))

PAR_nc.createVariable('lat', 'f4', ('lat',))
PAR_nc.createVariable('lon', 'f4', ('lon',))
PAR_nc.createVariable('time', 'f4', ('time',))
PAR_nc.createVariable('PAR_avg', 'f4', ('lat', 'lon', 'time'), fill_value=-999)

PAR_nc['lon'][:] = lon1
PAR_nc['lat'][:] = lat1
PAR_nc['time'][:] = time_SIF

# range(int(np.min(time_SIF)), int(np.max(time_SIF)), 8)
for i_tSIF in range(0, len(time_SIF)):
    SIF = SIF_nc.variables['sif_dc'][:, :, i_tSIF]
    print('Timestep: '+str(time_SIFgreg[i_tSIF]))
    i_timePAR = np.where(time_PAR == time_SIF[i_tSIF])[0][0]
    # Loading PAR data from the correct file
    if (i_timePAR >= 0) & (i_timePAR < 125):
        PAR = PAR0a_nc.variables['PAR_int'][:, :, i_timePAR:i_timePAR + 8]
    elif (i_timePAR >= 125) & (i_timePAR < 225):
        PAR = PAR0b_nc.variables['PAR_int'][:, :, i_timePAR:i_timePAR + 8]
    elif (i_timePAR >= 225) & (i_timePAR < 325):
        PAR = PAR0c_nc.variables['PAR_int'][:, :, i_timePAR:i_timePAR + 8]
    elif (i_timePAR >= 325) & (i_timePAR < 425):
        PAR = PAR1_nc.variables['PAR_int'][:, :, i_timePAR:i_timePAR + 8]
    elif (i_timePAR >= 425) & (i_timePAR < 525):
        PAR = PAR2_nc.variables['PAR_int'][:, :, i_timePAR:i_timePAR + 8]
    elif (i_timePAR >= 525) & (i_timePAR < 625):
        PAR = PAR3_nc.variables['PAR_int'][:, :, i_timePAR:i_timePAR + 8]
    elif (i_timePAR >= 625) & (i_timePAR < 725):
        PAR = PAR4_nc.variables['PAR_int'][:, :, i_timePAR:i_timePAR + 8]
    elif (i_timePAR >= 725) & (i_timePAR < 825):
        PAR = PAR5_nc.variables['PAR_int'][:, :, i_timePAR:i_timePAR + 8]
    elif (i_timePAR >= 825) & (i_timePAR < 930):
        PAR = PAR6_nc.variables['PAR_int'][:, :, i_timePAR:i_timePAR + 8]
    elif (i_timePAR >= 930) & (i_timePAR < 1099):
        PAR = PAR7_nc.variables['PAR_int'][:, :, i_timePAR:i_timePAR + 8]
    elif i_timePAR >= 1099:
        PAR = PAR8_nc.variables['PAR_int'][:, :, i_timePAR:i_timePAR + 8]
    # Averaging the data of 8-day period
    PAR_8day = np.nanmean(PAR, axis=2)
    PAR_nc['PAR_avg'][:, :, i_tSIF] = PAR_8day
    # Normalizing the 8-day SIF by the 8-day PAR
    if not np.isnan(SIF).all():
        cond = (PAR_8day < 0.1) | np.isnan(PAR_8day) | SIF.mask | (SIF.data == -999.0)
        SIF[cond] = np.nan
        PAR_8day[cond] = np.nan
        SIF_nc['sif_norm'][:, :, i_tSIF] = SIF/PAR_8day
        print('Normalized')
    else:
        SIF_nc['sif_norm'][:, :, i_tSIF] = np.nan

SIF_nc.close()
PAR_nc.close()
