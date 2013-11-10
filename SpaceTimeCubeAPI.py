# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Space Time Cube API

# <codecell>

import numpy as NUM
import SSDataObject as SSDO
from netCDF4 import Dataset

# <headingcell level=2>

# Create the Space Time Cube

# <headingcell level=5>

# Give me the cell size, Get a Space Time Cube

# <codecell>

class CubeObject(object):
    def __init__(self, ssdo, inputTimeField, cell_size_xy = 1000, cell_size_t = 30):
        self.ssdo = ssdo
        self.dataset = ssdo.inputFC
        self.timeField = inputTimeField
        self.cellSize_xy = cell_size_xy
        self.cellSize_t = cell_size_t

        self.extent = ssdo.extent
        self.spatialRef = ssdo.spatialRefName
        self.info = ssdo.info

    def basicScale(self):
        return self.createCube(self.cellSize_xy)

    def detailScale(self):
        cellSize_xy =  round (self.cellSize_xy / 2)
        return self.createCube(cellSize_xy)

    def coarseScale(self):
        cellSize_xy = (self.cellSize_xy * 2)
        return self.createCube(cellSize_xy)

    def createCube(self, cellSize_xy):
        cellNumber_x = round((self.extent.XMax - self.extent.XMin) / cellSize_xy)
        cellNumber_y = round((self.extent.YMax - self.extent.YMin) / cellSize_xy)
        X = self.ssdo.xyCoords[:,0]
        Y = self.ssdo.xyCoords[:,1]
        time = self.ssdo.fields[self.timeField].data
        time = NUM.array([i for i in time], NUM.datetime64)
        startDateTime = NUM.datetime64('1970-01-01 00:00:00')
        T = time - startDateTime
        self.startTime = NUM.amin(T) + NUM.datetime64('1970-01-01 00:00:00')
        T = NUM.array([i.item().days for i in T], int)
        startT = NUM.amin(T)
        endT = NUM.amax(T)
        cellNumber_t = round((endT - startT) / self.cellSize_t) + 1
        X = (X - self.extent.XMin) / self.cellSize_xy
        Y = (self.extent.YMax - Y) / self.cellSize_xy
        T = (T - startT) / self.cellSize_t
        X = NUM.floor(X)
        Y = NUM.floor(Y)
        T = NUM.floor(T)
        CellIdList = (cellNumber_x * cellNumber_y * T) + (cellNumber_x * Y) + X
        BothEnds = NUM.array([0, (cellNumber_t * cellNumber_x * cellNumber_y -1)])
        CellIdList = NUM.concatenate((CellIdList, BothEnds), axis=0)
        CellIdList = NUM.array(CellIdList, dtype = 'int32')
        counts = NUM.bincount(CellIdList)
        counts[BothEnds[0]] = counts[BothEnds[0]] - 1
        counts[BothEnds[1]] = counts[BothEnds[1]] - 1
        return counts.reshape(cellNumber_t, cellNumber_x, cellNumber_y)

# <rawcell>

# Requirements:
# 1. point data
# 2. time field
# 3. output is NetCDF4 database

# <headingcell level=5>

# Space Time Cube is stored in NetCDF4 database

# <codecell>

outputFC = 'C:/Data/SPAC/crime.nc'
outputDataset = Dataset(outputFC, 'w', format = 'NETCDF4')

# <rawcell>

# Why NetCDF4?
#
# 1. support multidimensional
# 2. fast computation
# 3. easy to obtain data and extend the database
# 4. archive old data

# <codecell>

inputFC = 'C:/Data/SPAC/crime.shp'
inputTimeField = 'DATE1'
ssdo = SSDO.SSDataObject(inputFC)
ssdo.obtainData(ssdo.oidName, [inputTimeField], dateStr = True)
cube1 = CubeObject(ssdo, inputTimeField, 1000, 30)
basic = cube1.basicScale()

# <rawcell>

# Three scales in cube API
# 1. basic scale as user defined
# 2. detailed scale: half cell size
# 3. course scale: double cell size

# <rawcell>

# Specify NetCDF details
#
# Need to honor Climate and Forecast (CF) conventions
# 1. Create dimensions and global attributes (projections...)
# 2. Create variables to store the cube
#
# ex: Time variables

# <codecell>

# Create netCDF4 dimensions and global attributes
outputDataset.createDimension('time', basic.shape[0])
outputDataset.createDimension('x', basic.shape[1])
outputDataset.createDimension('y', basic.shape[2])


# Create netCDF4 variables and attributes
time = outputDataset.createVariable('time','f8',('time'))
lat = outputDataset.createVariable('lat','f4', ('y','x'))
lon = outputDataset.createVariable('lon','f4', ('y','x'))
crime = outputDataset.createVariable('crime','f4',('time','y','x',))

projection = outputDataset.createVariable('transverse_mercator','i4')
xCoordinate = outputDataset.createVariable('X_Coordina','i4',('y','x'))
yCoordinate = outputDataset.createVariable('Y_Coordina','i4',('y','x'))
x = outputDataset.createVariable('x','f4',('x'))
y = outputDataset.createVariable('y','f4', ('y'))
timeInterval = 30* 60* 60* 24

time.long_name = 'time'
time.standard_name = 'time'
time.units = 'seconds since'+ `cube1.startTime`
time.calendar = 'gregorian'
time._ChunkSize = cube1.cellSize_t
time._CoordinateAxisType = 'Time'
##
lat.long_name = 'Latitude'
lat.standard_name = 'Latitude'
lat.coordinates = 'x y'
lat.esri_pe_string = 'PROJCS[\"NAD_1983_StatePlane_Illinois_East_FIPS_1201\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137.0,298.257222101]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",300000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-88.33333333333333],PARAMETER[\"Scale_Factor\",0.999975],PARAMETER[\"Latitude_Of_Origin\",36.66666666666666],UNIT[\"Meter\",1.0]]'
lat.grid_mapping = 'transverse_mercator'
##
lon.long_name = 'Longitude'
lon.standard_name = 'Longitude'
lon.coordinates = 'x y'
lon.esri_pe_string = 'PROJCS[\"NAD_1983_StatePlane_Illinois_East_FIPS_1201\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137.0,298.257222101]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",300000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-88.33333333333333],PARAMETER[\"Scale_Factor\",0.999975],PARAMETER[\"Latitude_Of_Origin\",36.66666666666666],UNIT[\"Meter\",1.0]]'
lon.grid_mapping = 'transverse_mercator'
##
crime.long_name = 'crime'
crime.standard_name = 'crime'
crime.coordinates = 'time x y'
crime.missing_value = -9999.
crime._ChunkSize = cube1.cellSize_t, cube1.cellSize_xy, cube1.cellSize_xy

projection.grid_mapping_name = 'transverse_mercator'
projection.longitude_of_central_meridian = -88.33333333333333
projection.latitude_of_projection_origin = 36.66666666666666
projection.scale_factor_at_central_meridian = 0.999975
projection.false_easting = 300000.0
projection.false_northing = 0.0
projection._CoordinateTransformType = 'Projection'
projection._CoordinateAxisTypes = 'GeoX GeoY'

xCoordinate.long_name = 'X_Coordina'
xCoordinate.standard_name = 'X_Coordina'
xCoordinate.coordinates = 'x y'
xCoordinate.esri_pe_string = 'PROJCS[\"NAD_1983_StatePlane_Illinois_East_FIPS_1201\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137.0,298.257222101]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",300000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-88.33333333333333],PARAMETER[\"Scale_Factor\",0.999975],PARAMETER[\"Latitude_Of_Origin\",36.66666666666666],UNIT[\"Meter\",1.0]]'
xCoordinate = 'transverse_mercator'

yCoordinate.long_name = 'Y_Coordina'
yCoordinate.standard_name = 'Y_Coordina'
yCoordinate.coordinates = 'x y'
yCoordinate.esri_pe_string = 'PROJCS[\"NAD_1983_StatePlane_Illinois_East_FIPS_1201\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137.0,298.257222101]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",300000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-88.33333333333333],PARAMETER[\"Scale_Factor\",0.999975],PARAMETER[\"Latitude_Of_Origin\",36.66666666666666],UNIT[\"Meter\",1.0]]'
yCoordinate = 'transverse_mercator'

x.long_name = 'x coordinate of projection'
x.standard_name = 'projection_x_coordinate'
x.grid_mapping = 'transverse_mercator'
x._CoordinateAxisType = 'GeoX'

y.long_name = 'y coordinate of projection'
y.standard_name = 'projection_y_coordinate'
y.grid_mapping = 'transverse_mercator'
y._CoordinateAxisType = 'GeoY'

# <headingcell level=5>

# Cube is Almost Done!

# <codecell>

time[:] = NUM.arange(0, basic.shape[0] * timeInterval, timeInterval)
x[:] = NUM.arange(cube1.extent.XMin, cube1.extent.XMax, cube1.cellSize_xy)
y[:] = NUM.arange(cube1.extent.YMax, cube1.extent.YMin, -cube1.cellSize_xy)
crime[:] = basic

outputDataset.close()

# <headingcell level=2>

# Do the Stats!

# <rawcell>

# To answer those questions, we need some functions

# <codecell>

class netCube(object):
    def __init__(self, dataset, variable, time=[0], lat=[0], lon=[0]):
        self.dataset = dataset
        self.variable = variable

    # Return the whole dataset
##    def whole(self):
##        return NUM.array(self.dataset.variables[self.variable], dtype = 'f')
    def whole(self):
        return self.dataset.variables[self.variable]

    # Extract one or more time-slice(Grid), retun a 3D numpy array
    def timeSlice(self, time):
        try:
            start_time, end_time = time[0], time[-1] + 1
        except:
            start_time, end_time = time, time + 1

        return NUM.array(self.dataset.variables[self.variable], dtype = 'f')[start_time:end_time, :, :]

    # Extract one or more rows, retun a 3D numpy array
    def row(self, lat):
        try:
            start_lat, end_lat = lat[0], lat[-1] + 1
        except:
            start_lat, end_lat = lat, lat + 1

        return NUM.array(self.dataset.variables[self.variable], dtype = 'f')[:, start_lat:end_lat, :]

    # Extract one or more columns, retun a 3D numpy array
    def col(self, lon):
        try:
            start_lon, end_lon = lon[0], lon[-1] + 1
        except:
            start_lon, end_lon = lon, lon + 1

        return NUM.array(self.dataset.variables[self.variable], dtype = 'f')[:, :, start_lon:end_lon]

    # Extract one cell from a series of time, return a 3D numpy array
    def timeSeries(self, lat, lon):
        try:
            start_lat, end_lat = lat[0], lat[-1] + 1
        except:
            start_lat, end_lat = lat, lat + 1
        try:
            start_lon, end_lon = lon[0], lon[-1] + 1
        except:
            start_lon, end_lon = lon, lon + 1

        return NUM.array(self.dataset.variables[self.variable], dtype = 'f')[:, start_lat:end_lat, start_lon:end_lon]


    # Extract space-time cube, return a 3D numpy array
    def timeCube(self, time, lat, lon):
        try:
            start_lat, end_lat = lat[0], lat[-1] + 1
        except:
            start_lat, end_lat = lat, lat + 1
        try:
            start_lon, end_lon = lon[0], lon[-1] + 1
        except:
            start_lon, end_lon = lon, lon + 1

        return NUM.array(self.dataset.variables[self.variable], dtype = 'f')[time, start_lat:end_lat, start_lon:end_lon]

    def GridToId(self, time, lat, lon):
        t_dim, lat_dim, lon_dim = self.dataset.variables[self.variable].shape
        IDs = []
        try:
            num_of_time = len(time)
            t = time
        except:
            num_of_time = 1
            t = []
            t.append(time)
        try:
            num_of_row = len(lat)
            n = lat
        except:
            num_of_row = 1
            n =[]
            n.append(lat)
        try:
            num_of_col = len(lon)
            k = lon
        except:
            num_of_col = 1
            k =[]
            k.append(lon)

        for i1 in range(num_of_time):
            for i2 in range(num_of_row):
                for i3 in range(num_of_col):
                    IDs.append (t[i1] *lat_dim *lon_dim + n[i2] *lon_dim + k[i3])
        return IDs

    def IdToGrid(self, id):
        t_dim, lat_dim, lon_dim = self.dataset.variables[self.variable].shape
        try:
            if len(id) >= 1:
                IDs = NUM.array(id)
                time = (IDs / (lat_dim *lon_dim))
                lat = (IDs % (lat_dim *lon_dim)) / lon_dim
                lon = ((IDs % (lat_dim *lon_dim)) % lon_dim)
            return NUM.unique(time), NUM.unique(lat), NUM.unique(lon)
        except:
            time = int(id / (lat_dim *lon_dim))
            lat = int((id % (lat_dim *lon_dim)) / lon_dim)
            lon = int((id % (lat_dim *lon_dim)) % lon_dim)
            return time, lat, lon

    def neighborByCellInfo(self, time, lat, lon, t_lag = 1,s_lag = 1):
        t_dim, lat_dim, lon_dim = self.dataset.variables[self.variable].shape
        start_time = time - t_lag
        end_time = time + t_lag
        start_lat = lat - s_lag
        end_lat = lat + s_lag
        start_lon = lon - s_lag
        end_lon = lon + s_lag
        if start_time < 0:
            start_time = 0
        if end_time >= t_dim:
            end_time = (t_dim-  1)

        if start_lat < 0:
            start_lat = 0
        if end_lat >= lat_dim:
            end_lat = (lat_dim - 1)

        if start_lon < 0:
            start_lon = 0
        if end_lon >= lon_dim:
            end_lon = (lon_dim - 1)
        return range(start_time, end_time+1), range(start_lat, end_lat+1), range(start_lon, end_lon+1)

    def neighborByID(self, id, t_lag = 1, s_lag = 1):
        time, lat, lon = self.IdToGrid(id)
        t, row, col = self.neighborByCellInfo(time, lat, lon, t_lag, s_lag)
        return self.GridToId(t,row,col)

    def getNeighborValuesByID(self, id, t_lag = 1, s_lag = 1):
        time, lat, lon = self.IdToGrid(id)
        t, row, col = self.neighborByCellInfo(time, lat, lon, t_lag, s_lag)
        return self.timeCube(t, row, col)

    def getNeighborValuesByCellInfo(self, time, lat, lon, t_lag = 1, s_lag = 1):
        t, row, col = self.neighborByCellInfo(time, lat, lon, t_lag, s_lag)
        return self.timeCube(t, row, col)

    def spatialLagByID(self, id):
        neighbor = self.getNeighborValuesByID(id, 0, 1)
        return (NUM.sum(neighbor) - NUM.sum(self.getNeighborValuesByID(id, 0, 0)))/(neighbor.size)

    def spatialLagByCellInfo(self, time, lat, lon):
        return self.spatialLagByID(self.GridToId(time, lat, lon))

    def sumByTime(self, Array):
        not_null = Array > -9999.
        Array[not_null == False] = 0.
        result = Array.sum(1).sum(1)
        Array[not_null == False] = -9999.
        return result

    def meanByTime(self, Array):
        not_null = Array > -9999.
        temp = self.sumByTime(Array)
        return (temp / not_null.sum(1).sum(1))

    def stdByTime(self, Array):
        not_null = Array > -9999.
        Array[not_null == False] = NUM.repeat(self.meanByTime(Array),((not_null == False).sum(1).sum(1)))
        result = Array.std(1).std(1)
        Array[not_null == False] = -9999.
        return result

    def zByTime(self, Array):
        not_null = Array > -9999.
        temp = NUM.array(Array)
        temp[not_null] = (Array[not_null] - NUM.repeat(self.meanByTime(Array), not_null.sum(1).sum(1))) / NUM.repeat(self.stdByTime(Array), not_null.sum(1).sum(1))
        return temp

# <headingcell level=5>

# Read the Cube stored in NetCDF4 database

# <codecell>

inputFC = 'C:/Data/SPAC/crime.nc'
inputDataset = Dataset(inputFC, 'a')
inputVariable = "crime"
crime = inputDataset.variables['crime']
nc1 = netCube(inputDataset, inputVariable)

# <headingcell level=3>

# Functions in the cube

# <headingcell level=6>

# Extraction data

# <rawcell>

# I　want to...
#
# have an overview of the whole data
# know the spatial distribution over time
# know the change in certain locations
##
### <codecell>
##
print nc1.whole().type
##print 'Whole dataset shape:\n', nc1.whole().shape
##print '1 time slice:\n', nc1.timeSlice(3)
##print '2 rows over time shape:\n', nc1.row([1,2]).shape
##print '5 columns over time shape:\n', nc1.col([1,2,3,4,5]).shape
##print '4 cells over time:\n', nc1.timeSeries([33,34],[55,56])
##
### <headingcell level=5>
##
### Neighbor functions
##
### <rawcell>
##
### I want to..
###
### group the elements
### know the average of the surrounding area
### construct the weight by the neighbors
##
### <codecell>
##
##print 'cell (3,58,120) to ID:\n', nc1.GridToId(3,58,120)
##print 'Neighbors of Grid(3,58,120) are\n', nc1.neighborByCellInfo(3,58,120)
##print 'Neighbor values of Grid(3,58,120) is \n', nc1.getNeighborValuesByCellInfo(3,58,120)
##print 'Spatial lags of ID 117640 is \n', nc1.spatialLagByID(117640)
##print 'Spatial lags of Grid(3,58,120) is \n', nc1.spatialLagByCellInfo(3,58,120)
##
### <headingcell level=5>
##
### Simple Stats
##
### <rawcell>
##
### I want to..
###
### know how different is between two time periods
### scale it
##
### <codecell>
##
##print 'Sum of Time Slice 5 are \n', nc1.sumByTime(nc1.timeSlice(5))
##print 'Standard deviation of whole dataset are \n', nc1.stdByTime(nc1.whole())
##print 'Z scores of Time Slice 5 are \n', nc1.zByTime(nc1.timeSlice(5))
##
### <headingcell level=3>
##
### Hot Spot
##
### <rawcell>
##
### I want to...
###
### know where and when is the hot spots and cold spots
### find out the trajectories (does the crime moves?)
### surveillence the crime incidents
##
### <codecell>
##
##import hotspot_grid as HOTSPOT
##
##Hotspot = inputDataset.createVariable('Hotspot','f4',('time','y','x',))
### Variable 'Hot Spot'
##Hotspot.long_name = 'Hotspot'
##Hotspot.standard_name = 'Hotspot'
##Hotspot.coordinates = 'time x y'
##Hotspot.units = 'meters'
##Hotspot.missing_value = -9999.
##Hotspot._ChunkSize = cube1.cellSize_t, cube1.cellSize_xy, cube1.cellSize_xy
##
##numTime, numRows, numCols = nc1.whole().shape
##for time in xrange(numTime):
##    cellValues = nc1.timeSlice(time).flatten()
##    gi, pv = HOTSPOT.hotspotGrid(cellValues, numRows, numCols)
##    Hotspot[time,:,:] = gi.reshape(numRows, numCols)
##
##inputDataset.close()
##
### <rawcell>
##
##
### <codecell>


