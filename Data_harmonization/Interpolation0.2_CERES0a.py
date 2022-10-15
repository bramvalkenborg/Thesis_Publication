from SIF_tools.python.L2_tools import convert_time
import numpy as np
import datetime
from matplotlib import dates
from netCDF4 import Dataset
import cftime

# Load all the data
CERES_nc = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/CERES_grid0-2deg_SYN1deg-1H_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130_NH.nc')
SIF_nc = Dataset('/data/leuven/336/vsc33653/OUTPUT/Gridding/NH_SIF_2021_dN_omitnan_1day_0.2deg_c30.nc')
t_CERES = CERES_nc.variables['time'][:]
t_CERES = convert_time(t_CERES, 'days since 2000-03-01 00:00:00', u"gregorian")
time_SIF = SIF_nc.variables['time'][:]
time_SIFgreg = convert_time(time_SIF, 'days since 1970-01-01', u"gregorian")
month_SIF = time_SIFgreg.astype('datetime64[M]').astype(int) % 12 + 1
lat_SIF = SIF_nc.variables['lat'][:]
lon_SIF = SIF_nc.variables['lon'][:]
# Load the peatmap, to save calculation time, only peatpixels are interpolated
Peat_nc = Dataset('/data/leuven/336/vsc33653/OUTPUT/Gridding/NH_PeatMap_0.2deg.nc')
peat = Peat_nc.variables['peatmap'][:, :, 0]
peat_mask = (peat != 0)

#PAR_direct = CERES_nc.variables['adj_sfc_par_direct_all_1h'][:, :, :]
#PAR_diffuse = CERES_nc.variables['adj_sfc_par_diff_all_1h'][:, :, :]
# PAR = PAR_direct

nSIF = 1099
#t_SIF = t_SIF[:, :, 0:975]

# create an array that for each SIF value contains the 2 closest PAR measurements (in time)
# i_tCERES_lower = np.zeros((len(t_SIF), 1))
# i_tCERES_upper = np.zeros((len(t_SIF), 1))
# PAR_lower = np.zeros((len(t_SIF), PAR.shape[1], PAR.shape[2]))
# PAR_upper = np.zeros((len(t_SIF), PAR.shape[1], PAR.shape[2]))
# for i_tSIF in range(len(t_SIF)):
    # time = np.full(t_CERES.shape, t_SIF[i_tSIF])
    # t_diff = abs(t_CERES - time)
    # i_time = np.where(t_diff == t_diff.min())[0]
    # i_tCERES_lower[i_tSIF] = i_time[0]
    # i_tCERES_upper[i_tSIF] = i_time[1]
    # PAR_lower[i_tSIF, :, :] = PAR[i_time[0], :, :]
    # PAR_upper[i_tSIF, :, :] = PAR[i_time[1], :, :]
    # print('SIF time: '+str((i_tSIF/len(t_SIF))*100)+'%')

PARint = Dataset('/scratch/leuven/336/vsc33653/gridding/CERES/CERES_0-2deg_int_t0_t125_NH_SYN1deg_Terra-Aqua-MODIS_Ed4.1_Subset_20171101-2021130.nc',
                    mode='w', format='NETCDF4')
PARint.createDimension('time', len(time_SIF))
PARint.createDimension('lat', len(lat_SIF))
PARint.createDimension('lon', len(lon_SIF))
PARint.createVariable('time', 'f4', ('time',))
PARint.createVariable('lat', 'f4', ('lat',))
PARint.createVariable('lon', 'f4', ('lon',))
PARint['time'][:] = time_SIF
PARint['lat'][:] = lat_SIF
PARint['lon'][:] = lon_SIF
PARint.createVariable('PAR_int', 'f4', ('lat', 'lon', 'time'), fill_value=-999)
# PAR with exact time and not time_sif

cond_SIF = (month_SIF == 5) | (month_SIF == 6) | (month_SIF == 7) | (month_SIF == 8) | (month_SIF == 9) | (month_SIF == 10)
tCERES_x = [0, 0]
PAR_y = [0, 0]
for t in range(0, 125):
    # We only want an interpolation for the months of the growing season
    if cond_SIF[t]:
        # Loading the exact time of TROPOMI measurement
        print('Loading data: ' + str(t))
        t_SIF = convert_time(SIF_nc.variables['time_exact'][:, :, t], 'seconds since 1970-01-01 00:00:00', u"gregorian")
        month = t_SIF.astype('datetime64[M]').astype(int) % 12 + 1
        cond = (month == 5) | (month == 6) | (month == 7) | (month == 8) | (month == 9) | (month == 10)
        # if all of the pixels don't belong for growing season data, this timestep should be skipped
        if (cond == False).all():
            print('Timestep: ' + str(t))
            continue
        for lat in range(len(lat_SIF)):
            for lon in range(len(lon_SIF)):
                # Only interpolate for peatland pixels
                if peat_mask[lat, lon]:
                    # t_SIF = t_SIF[lat, lon]
                    # Only interpolate for the growing season
                    if cond[lat, lon]:
                        tSIF = t_SIF[lat, lon]
                        time = np.full(t_CERES.shape, tSIF)
                        t_diff = abs(t_CERES - time)
                        i_time = np.where(t_diff == t_diff.min())[0][0]
                        if tSIF > t_CERES[i_time]:
                            tCERES_x[0] = t_CERES[i_time]
                            tCERES_x[1] = t_CERES[i_time + 1]
                            PAR_y[0] = CERES_nc.variables['adj_sfc_par_direct_all_1h'][i_time, lat, lon] + CERES_nc.variables['adj_sfc_par_diff_all_1h'][i_time, lat, lon]
                            PAR_y[1] = CERES_nc.variables['adj_sfc_par_direct_all_1h'][i_time+1, lat, lon] + CERES_nc.variables['adj_sfc_par_diff_all_1h'][i_time+1, lat, lon]
                        if tSIF < t_CERES[i_time]:
                            tCERES_x[0] = t_CERES[i_time - 1]
                            tCERES_x[1] = t_CERES[i_time]
                            PAR_y[0] = CERES_nc.variables['adj_sfc_par_direct_all_1h'][i_time-1, lat, lon] + CERES_nc.variables['adj_sfc_par_diff_all_1h'][i_time-1, lat, lon]
                            PAR_y[1] = CERES_nc.variables['adj_sfc_par_direct_all_1h'][i_time, lat, lon] + CERES_nc.variables['adj_sfc_par_diff_all_1h'][i_time, lat, lon]
                        # tCERES_x[0] = int(tCERES_x[0].item().strftime("%Y%m%d%H%M%S"))
                        # tCERES_x[1] = int(tCERES_x[1].item().strftime("%Y%m%d%H%M%S"))
                        # PAR_int[t, lat, lon] = np.interp(int(tSIF.strftime("%Y%m%d%H%M%S")), tCERES_x, PAR_y)
                        PARint['PAR_int'][lat, lon, t] = PAR_y[0] + (tSIF - tCERES_x[0]).seconds / (tCERES_x[1] - tCERES_x[0]).seconds * (PAR_y[1] - PAR_y[0])
                        print('Timestep: ' + str(t) + ' | ' + 'lon: ' + str(lon) + ', lat: ' + str(
                            lat) + ' | Peat soil | Interpolation')
                        # PAR_y1 = np.transpose(np.asarray(PAR_y))
                        # tCERES_x1 = np.transpose(tCERES_x)
                        # PAR_int[lat, lon, t] = np.interp(tSIF_num, tCERES_x1, PAR_y1)
                    else:
                        PARint['PAR_int'][lat, lon, t] = np.nan
                        print('Timestep: ' + str(t) + ' | ' + 'lon: ' + str(lon) + ', lat: ' + str(
                           lat))
                    tCERES_x = [0, 0]
                    PAR_y = [0, 0]
                else:
                    PARint['PAR_int'][lat, lon, t] = np.nan
                print('Timestep: ' + str(t) + ' | ' + 'lon: ' + str(lon) + ', lat: ' + str(lat))
PARint.close()
