import netCDF4 as nc
import numpy as np
from SIF_tools.python.L2_tools import convert_time
from datetime import datetime


# 10 year data
WTD1 = nc.Dataset('/staging/leuven/stg_00024/OUTPUT/michelb/SMOS_mwRTM_EASEv2_M09_PEATCLSM_IcPmNO_DA/output/SMAP_EASEv2_M09/cat/ens_avg/daily_images.nc')
WTD = nc.Dataset('/staging/leuven/stg_00024/OUTPUT/michelb/SMOS_mwRTM_EASEv2_M09_PEATCLSM_IcPmNO_DA_v2/output/SMAP_EASEv2_M09/cat/ens_avg/daily_images.nc')
t_start = '2010-01-01'

lon1 = WTD.variables['lon'][:]
lat1 = WTD.variables['lat'][:]
time = WTD.variables['time'][:]
# Select the pixels where the max rzmc over time (axis=0) is larger than 0.65
peat_mask = np.nanmax(WTD.variables['rzmc'][:, :, :].data, axis=0) > 0.65

# Calculate the size of each cell in the lat and lon direction, outer pixels cannot be taken into account because for
# them the pixel size cannot be calculated. Therefore the outer pixels of lat1 and lon1 are removed.
# 1: upper left corner of a pixel
# 2: upper right corner of a pixel
# 3: lower right corner of a pixel
# 4: lower left corner of a pixel
# For the outer edges of the area (vectors, no 2D matrix)
d_lat = np.diff(lat1)
lat1_bnds1 = lat1[0:len(lat1)-1] - (d_lat/2)
lat1_bnds1 = np.insert(lat1_bnds1, len(lat1_bnds1), lat1[len(lat1)-1] - (d_lat[len(d_lat)-1]/2))
lat1_bnds2 = lat1[0:len(lat1)-1] - (d_lat/2)
lat1_bnds2 = np.insert(lat1_bnds2, len(lat1_bnds2), lat1[len(lat1)-1] - (d_lat[len(d_lat)-1]/2))
lat1_bnds3 = lat1[0:len(lat1)-1] + (d_lat/2)
lat1_bnds3 = np.insert(lat1_bnds3, len(lat1_bnds3), lat1[len(lat1)-1] + (d_lat[len(d_lat)-1]/2))
lat1_bnds4 = lat1[0:len(lat1)-1] + (d_lat/2)
lat1_bnds4 = np.insert(lat1_bnds4, len(lat1_bnds4), lat1[len(lat1)-1] + (d_lat[len(d_lat)-1]/2))

d_lon = np.diff(lon1)
lon1_bnds1 = lon1[0:len(lon1)-1] - (d_lon/2)
lon1_bnds1 = np.insert(lon1_bnds1, len(lon1_bnds1), lon1[len(lon1)-1] - (d_lon[len(d_lon)-1]/2))
lon1_bnds2 = lon1[0:len(lon1)-1] + (d_lon/2)
lon1_bnds2 = np.insert(lon1_bnds2, len(lon1_bnds2), lon1[len(lon1)-1] + (d_lon[len(d_lon)-1]/2))
lon1_bnds3 = lon1[0:len(lon1)-1] + (d_lon/2)
lon1_bnds3 = np.insert(lon1_bnds3, len(lon1_bnds3), lon1[len(lon1)-1] + (d_lon[len(d_lon)-1]/2))
lon1_bnds4 = lon1[0:len(lon1)-1] - (d_lon/2)
lon1_bnds4 = np.insert(lon1_bnds4, len(lon1_bnds4), lon1[len(lon1)-1] - (d_lon[len(d_lon)-1]/2))

lon = np.tile(lon1, len(lat1))
lon_bnds1 = np.tile(lon1_bnds1, len(lat1))
lon_bnds2 = np.tile(lon1_bnds2, len(lat1))
lon_bnds3 = np.tile(lon1_bnds3, len(lat1))
lon_bnds4 = np.tile(lon1_bnds4, len(lat1))

lat = np.repeat(lat1, len(lon1))
lat_bnds1 = np.repeat(lat1_bnds1, len(lon1))
lat_bnds2 = np.repeat(lat1_bnds2, len(lon1))
lat_bnds3 = np.repeat(lat1_bnds3, len(lon1))
lat_bnds4 = np.repeat(lat1_bnds4, len(lon1))

t_greg1 = convert_time(time, "hours since 2000-01-01 00:00", u"gregorian")
i_t_greg1 = np.where(t_greg1 == datetime.fromisoformat(t_start))

# Total amount of gridcells lat * lon, all pixels incl. non-peat (not over time)
ngridcells = np.shape(WTD.variables['zbar'][0, :, :])[0] * np.shape(WTD.variables['zbar'][0, :, :])[1]
# Reshape the peatmask to a 1D vector: ngridcells*1
peat_mask_1D = np.reshape(peat_mask, (ngridcells, 1))
for t in range(int(i_t_greg1[0]), len(time)):
    t_greg = t_greg1.tolist()[t].isoformat(' ')
    WTD_time = nc.Dataset('/scratch/leuven/336/vsc33653/gridding/WTD/WTD_'+t_greg[0:4]+'_'+t_greg[5:7]+'_'+t_greg[8:10]+'.nc',
                          mode='w', format='NETCDF4')
    #MB: leave the followin three lines in please, for
    #t_greg = str(t_greg1.tolist()[t].year)+"_"+str(t_greg1.tolist()[t].month).zfill(2)+"_"+str(t_greg1.tolist()[t].day).zfill(2)
    #WTD_time = nc.Dataset('/scratch/leuven/317/vsc31786/output/WTD/WTD_'+t_greg[0:4]+'_'+t_greg[5:7]+'_'+t_greg[8:10],
    #                      mode='w', format='NETCDF4')
    # WTD_time.createDimension('lon', len(lon))
    # WTD_time.createDimension('lat', len(lat))

    # Create the different dimensions and variables
    # WTD_time.createDimension('TIME', len(lon[peat_mask_1D[:, 0]]))
    WTD_time.createDimension('TIME', len(lon))
    WTD_time.createDimension('bnds', 4)

    WTD_time.createVariable('TIME', 'f4', ('TIME',))
    WTD_time.createVariable('lon', 'f4', ('TIME',))
    WTD_time.createVariable('lat', 'f4', ('TIME',))
    WTD_time.createVariable('zbar', 'f4', ('TIME',))
    WTD_time.createVariable('sfmc', 'f4', ('TIME',))
    WTD_time.createVariable('rzmc', 'f4', ('TIME',))
    #WTD_time.createVariable('tsurf', 'f4', ('TIME',))
    WTD_time.createVariable('tp1', 'f4', ('TIME',))
    #WTD_time.createVariable('lhflux', 'f4', ('TIME',))
    WTD_time.createVariable('bnds', 'f4', ('bnds',))
    WTD_time.createVariable('lat_bnds', 'f4', ('bnds', 'TIME'))
    WTD_time.createVariable('lon_bnds', 'f4', ('bnds', 'TIME'))

    # Complete the different dimensions and variables with values
    # WTD_time['bnds'][:] = [1, 2, 3, 4]
    # WTD_time['TIME'][:] = range(len(lon[peat_mask_1D[:, 0]]))
    # WTD_time['lon'][:] = lon[peat_mask_1D[:, 0]]
    # WTD_time['lat'][:] = lat[peat_mask_1D[:, 0]]
    # WTD_time['lat_bnds'][0, :] = lat_bnds1[peat_mask_1D[:, 0]]
    # WTD_time['lat_bnds'][1, :] = lat_bnds2[peat_mask_1D[:, 0]]
    # WTD_time['lat_bnds'][2, :] = lat_bnds3[peat_mask_1D[:, 0]]
    # WTD_time['lat_bnds'][3, :] = lat_bnds4[peat_mask_1D[:, 0]]
    # WTD_time['lon_bnds'][0, :] = lon_bnds1[peat_mask_1D[:, 0]]
    # WTD_time['lon_bnds'][1, :] = lon_bnds2[peat_mask_1D[:, 0]]
    # WTD_time['lon_bnds'][2, :] = lon_bnds3[peat_mask_1D[:, 0]]
    # WTD_time['lon_bnds'][3, :] = lon_bnds4[peat_mask_1D[:, 0]]

    WTD_time['bnds'][:] = [1, 2, 3, 4]
    WTD_time['TIME'][:] = range(len(lon))
    WTD_time['lon'][:] = lon
    WTD_time['lat'][:] = lat
    WTD_time['lat_bnds'][0, :] = lat_bnds1
    WTD_time['lat_bnds'][1, :] = lat_bnds2
    WTD_time['lat_bnds'][2, :] = lat_bnds3
    WTD_time['lat_bnds'][3, :] = lat_bnds4
    WTD_time['lon_bnds'][0, :] = lon_bnds1
    WTD_time['lon_bnds'][1, :] = lon_bnds2
    WTD_time['lon_bnds'][2, :] = lon_bnds3
    WTD_time['lon_bnds'][3, :] = lon_bnds4
    # Variables are reshaped to the same vector
    # In order to not average the wtd of a mineral soil with a peat soil, all mineral soil water table depth are set to
    # nan by using the peatmask.
    # Before in a while loop per lat, takes too much time to calculate

    # zbar = np.reshape(WTD.variables['zbar'][t, 0:100, 1000:1100], (ngridcells, ))[peat_mask_1D]
    # rzmc = np.reshape(WTD.variables['rzmc'][t, 0:100, 1000:1100], (ngridcells, ))[peat_mask_1D]
    # sfmc = np.reshape(WTD.variables['sfmc'][t, 0:100, 1000:1100], (ngridcells, ))[peat_mask_1D]
    # tp1 = np.reshape(WTD.variables['tp1'][t, 0:100, 1000:1100], (ngridcells, ))[peat_mask_1D]
    zbar = WTD.variables['zbar'][t, :, :]
    zbar = np.where(peat_mask, zbar, np.nan)
    zbar = np.reshape(zbar, (ngridcells, ))
    rzmc = WTD.variables['rzmc'][t, :, :]
    rzmc = np.where(peat_mask, rzmc, np.nan)
    rzmc = np.reshape(rzmc, (ngridcells,))
    sfmc = WTD.variables['sfmc'][t, :, :]
    sfmc = np.where(peat_mask, sfmc, np.nan)
    sfmc = np.reshape(sfmc, (ngridcells,))
    tp1 = WTD.variables['tp1'][t, :, :]
    tp1 = np.where(peat_mask, tp1, np.nan)
    tp1 = np.reshape(tp1, (ngridcells,))
    WTD_time['zbar'][:] = zbar[:]
    #WTD_time['srfexc'][:] = srfexc[:]
    #WTD_time['rzexc'][:] = rzexc[:]
    #WTD_time['catdef'][:] = catdef[:]
    #WTD_time['ar1'][:] = ar1[:]
    #WTD_time['ar2'][:] = ar2[:]
    WTD_time['sfmc'][:] = sfmc[:]
    WTD_time['rzmc'][:] = rzmc[:]
    #WTD_time['tsurf'][:] = tsurf[:]
    WTD_time['tp1'][:] = tp1[:]
    #WTD_time['tpN'][:] = tpN[:]
    #WTD_time['shflux'][:] = shflux[:]
    #WTD_time['lhflux'][:] = lhflux[:]
    #WTD_time['evap'][:] = evap[:]
    #WTD_time['runoff'][:] = runoff[:]
    print(t_greg[0:4] + '_' + t_greg[5:7] + '_' + t_greg[8:10])
    WTD_time.close()
