import netCDF4 as nc
import numpy as np

llcrnrlat = 44
urcrnrlat = 71
llcrnrlon = -171
urcrnrlon = 96

peat = nc.Dataset('/staging/leuven/stg_00024/OUTPUT/michelb/l_data/model_param/PEATMAP_HWSD/PEATMAP_HWSD_mask.nc4',
               mode='r')
lon1 = peat.variables['longitude'][:]
lat1 = peat.variables['latitude'][:]
peatmap1 = peat.variables['PEATMAP'][:, :]

i_lonmin = (np.abs(lon1-llcrnrlon)).argmin()
i_lonmax = (np.abs(lon1-urcrnrlon)).argmin()
i_latmin = (np.abs(lat1-llcrnrlat)).argmin()
i_latmax = (np.abs(lat1-urcrnrlat)).argmin()
lon1 = lon1[i_lonmin:i_lonmax]
lat1 = lat1[i_latmin:i_latmax]
peatmap1 = peatmap1[i_latmin:i_latmax, i_lonmin:i_lonmax]

# Calculate the size of each cell in the lat and lon direction, outer pixels cannot be taken into account because for
# them the pixel size cannot be calculated. Therefore the outer pixels of lat1 and lon1 are removed.
# 1: upper left corner of a pixel
# 2: upper right corner of a pixel
# 3: lower right corner of a pixel
# 4: lower left corner of a pixel
# CHECK THE SIGNS!!!!!!
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

"""
# In order to read it in Julia, vectors are used and not 2D arrays. Each element of the vector is a pixel.
# lon1 shows the different longitudes for a certain latitude, therefore lon1 should be repeated len(lat1) times.
lon_bnds1 = np.zeros((1, 1))
lon_bnds2 = np.zeros((1, 1))
lon_bnds3 = np.zeros((1, 1))
lon_bnds4 = np.zeros((1, 1))
lon = np.zeros((1, 1))
i = 0
while i <= (len(lat1)-1):
    lon = np.insert(lon, len(lon), lon1)
    lon_bnds1 = np.insert(lon_bnds1, len(lon_bnds1), lon1_bnds1)
    lon_bnds2 = np.insert(lon_bnds2, len(lon_bnds2), lon1_bnds2)
    lon_bnds3 = np.insert(lon_bnds3, len(lon_bnds3), lon1_bnds3)
    lon_bnds4 = np.insert(lon_bnds4, len(lon_bnds4), lon1_bnds4)
    print("longitude: " + str(i))
    i += 1
lon = np.delete(lon, 0)
lon_bnds1 = np.delete(lon_bnds1, 0)
lon_bnds2 = np.delete(lon_bnds2, 0)
lon_bnds3 = np.delete(lon_bnds3, 0)
lon_bnds4 = np.delete(lon_bnds4, 0)

lat_bnds1 = np.zeros((1, 1))
lat_bnds2 = np.zeros((1, 1))
lat_bnds3 = np.zeros((1, 1))
lat_bnds4 = np.zeros((1, 1))
lat = np.zeros((1, 1))
for i in range(len(lat1)):
    lat = np.insert(lat, len(lat), np.full(len(lon1), lat1[i]))
    lat_bnds1 = np.insert(lat_bnds1, len(lat_bnds1), np.full(len(lon1), lat1_bnds1[i]))
    lat_bnds2 = np.insert(lat_bnds2, len(lat_bnds2), np.full(len(lon1), lat1_bnds2[i]))
    lat_bnds3 = np.insert(lat_bnds3, len(lat_bnds3), np.full(len(lon1), lat1_bnds3[i]))
    lat_bnds4 = np.insert(lat_bnds4, len(lat_bnds4), np.full(len(lon1), lat1_bnds4[i]))
    print("latitude: " + str(i))
lat = np.delete(lat, 0)
lat_bnds1 = np.delete(lat_bnds1, 0)
lat_bnds2 = np.delete(lat_bnds2, 0)
lat_bnds3 = np.delete(lat_bnds3, 0)
lat_bnds4 = np.delete(lat_bnds4, 0)
"""

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

# Create netCDF file
# Not really a date, but used for Julia
peat_netCDF = nc.Dataset('/scratch/leuven/336/vsc33653/gridding/PeatMap/NH_PeatMap_2020_01_01.nc',
                         mode='w', format='NETCDF4')
# WTD_time.createDimension('lon', len(lon))
# WTD_time.createDimension('lat', len(lat))
peat_netCDF.createDimension('TIME', len(lon))
peat_netCDF.createDimension('bnds', 4)

peat_netCDF.createVariable('TIME', 'f4', ('TIME',))
peat_netCDF.createVariable('lon', 'f4', ('TIME',))
peat_netCDF.createVariable('lat', 'f4', ('TIME',))
peat_netCDF.createVariable('PEATMAP', 'f4', ('TIME',))
peat_netCDF.createVariable('bnds', 'f4', ('bnds',))
peat_netCDF.createVariable('lat_bnds', 'f4', ('bnds', 'TIME'))
peat_netCDF.createVariable('lon_bnds', 'f4', ('bnds', 'TIME'))

peat_netCDF['bnds'][:] = [1, 2, 3, 4]
peat_netCDF['TIME'][:] = range(len(lon))
peat_netCDF['lon'][:] = lon
peat_netCDF['lat'][:] = lat
peat_netCDF['lat_bnds'][0, :] = lat_bnds1
peat_netCDF['lat_bnds'][1, :] = lat_bnds2
peat_netCDF['lat_bnds'][2, :] = lat_bnds3
peat_netCDF['lat_bnds'][3, :] = lat_bnds4
peat_netCDF['lon_bnds'][0, :] = lon_bnds1
peat_netCDF['lon_bnds'][1, :] = lon_bnds2
peat_netCDF['lon_bnds'][2, :] = lon_bnds3
peat_netCDF['lon_bnds'][3, :] = lon_bnds4
PEATMAP = np.zeros((1, 1))
i = 0
while i <= (len(lat1)-1):
    PEATMAP = np.insert(PEATMAP, len(PEATMAP), peatmap1[i, :])
    i += 1
PEATMAP = np.delete(PEATMAP, 0)
peat_netCDF['PEATMAP'][:] = PEATMAP[:]
peat_netCDF.close()
