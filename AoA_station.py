# Script to merge the station output files of the Age of Air experiments of august 2017.
# Call as python AoA_station.py startdate enddate
# with startdate and enddate in format YYYYMMDD, example: python AoA_station.py 19900101 19900401

import numpy as np
import netCDF4 as nc
import os,sys
from scipy.io import netcdf

startdate=sys.argv[1]
enddate=sys.argv[2]
outputfilename='/scratch/ms/nl/nklk/ECEARTH-RUNS/aoa1/'+startdate+'_'+enddate+'.nc'

station_fh = nc.Dataset('/scratch/ms/nl/nklk/ECEARTH-RUNS/aoa1/station_output_00_'+startdate+'00_'+enddate+'00.nc')
n_stations = station_fh.variables['n_stations'][:]
n_records = station_fh.variables['n_records'][:]
station_fh.close()

surf_grd_full = np.ones([len(n_records),len(n_stations)])*-999
surf_slp_full = np.ones([len(n_records),len(n_stations)])*-999
surf_int_full = np.ones([len(n_records),len(n_stations)])*-999
land_grd_full = np.ones([len(n_records),len(n_stations)])*-999
land_slp_full = np.ones([len(n_records),len(n_stations)])*-999
land_int_full = np.ones([len(n_records),len(n_stations)])*-999
nh_grd_full = np.ones([len(n_records),len(n_stations)])*-999
nh_slp_full = np.ones([len(n_records),len(n_stations)])*-999
nh_int_full = np.ones([len(n_records),len(n_stations)])*-999
sh_grd_full = np.ones([len(n_records),len(n_stations)])*-999
sh_slp_full = np.ones([len(n_records),len(n_stations)])*-999
sh_int_full = np.ones([len(n_records),len(n_stations)])*-999
strat_grd_full = np.ones([len(n_records),len(n_stations)])*-999
strat_slp_full = np.ones([len(n_records),len(n_stations)])*-999
strat_int_full = np.ones([len(n_records),len(n_stations)])*-999
e90_grd_full = np.ones([len(n_records),len(n_stations)])*-999
e90_slp_full = np.ones([len(n_records),len(n_stations)])*-999
e90_int_full = np.ones([len(n_records),len(n_stations)])*-999

for root, dirs, files in os.walk('/scratch/ms/nl/nklk/ECEARTH-RUNS/aoa1'):
    for file in sorted(files):
        if file.endswith(enddate+'00.nc'):
            print(file)
            station_fh = nc.Dataset(os.path.join('/scratch/ms/nl/nklk/ECEARTH-RUNS/aoa1/',file),mode='r')
            surf_grd = station_fh.variables['surf_grd'][:]
            n_stations = station_fh.variables['n_stations'][:]
            n_records = station_fh.variables['n_records'][:]
            surf_slp = station_fh.variables['surf_slp'][:]
            surf_int = station_fh.variables['surf_int'][:]
            land_grd = station_fh.variables['land_grd'][:]
            land_slp = station_fh.variables['land_slp'][:]
            land_int = station_fh.variables['land_int'][:]
            nh_grd = station_fh.variables['nh_grd'][:]
            nh_slp = station_fh.variables['nh_slp'][:]
            nh_int = station_fh.variables['nh_int'][:]
            sh_grd = station_fh.variables['sh_grd'][:]
            sh_slp = station_fh.variables['sh_slp'][:]
            sh_int = station_fh.variables['sh_int'][:]
            strat_grd = station_fh.variables['strat_grd'][:]
            strat_slp = station_fh.variables['strat_slp'][:]
            strat_int = station_fh.variables['strat_int'][:]
            e90_grd = station_fh.variables['e90_grd'][:]
            e90_slp = station_fh.variables['e90_slp'][:]
            e90_int = station_fh.variables['e90_int'][:]
            
            #Global attributes
            height_sample_above_surf = station_fh.height_sample_above_surf[:]
            height_surface = station_fh.height_surface[:]
            height_sample = station_fh.height_sample[:]
            station_ls = station_fh.station_ls[:]
            station_ident = station_fh.station_ident[:]
            station_names = station_fh.station_names[:]
            station_js = station_fh.station_js[:]
            station_is = station_fh.station_is[:]
            station_region = station_fh.station_region[:]
            station_index = station_fh.station_index[:]
            station_nonsurf = station_fh.station_nonsurf[:]
            station_height = station_fh.station_height[:]
            station_lon = station_fh.station_lon[:]
            station_lat = station_fh.station_lat[:]
	    station_fh.close()

	    for time in range(len(n_records)):
		    for i, j in enumerate(surf_grd[time]):
			    if j > 0:
				    #print(i)
                                    #print(station_index[i])
                                    #print(Temperature[time,i])
                                    surf_grd_full[time,i] = surf_grd[time,i]
                                    surf_slp_full[time,i] = surf_slp[time,i]
                                    surf_int_full[time,i] = surf_int[time,i]
                                    land_grd_full[time,i] = land_grd[time,i]
                                    land_slp_full[time,i] = land_slp[time,i]
                                    land_int_full[time,i] = land_int[time,i]
                                    nh_grd_full[time,i] = nh_grd[time,i]
                                    nh_slp_full[time,i] = nh_slp[time,i]
                                    nh_int_full[time,i] = nh_int[time,i]
                                    sh_grd_full[time,i] = sh_grd[time,i]
                                    sh_slp_full[time,i] = sh_slp[time,i]
                                    sh_int_full[time,i] = sh_int[time,i]
                                    strat_grd_full[time,i] = strat_grd[time,i]
                                    strat_slp_full[time,i] = strat_slp[time,i]
                                    strat_int_full[time,i] = strat_int[time,i]
                                    e90_grd_full[time,i] = e90_grd[time,i]
                                    e90_slp_full[time,i] = e90_slp[time,i]
                                    e90_int_full[time,i] = e90_int[time,i]
                                    #Temperature_full[time,i] = Temperature[time,i]
                                    #Pressure_full[time,i] = Pressure[time,i]
                                    #uwind_full[time,i] = uwind[time,i]
                                    #vwind_full[time,i] = vwind[time,i]
                                    #windspeed_full[time,i] = windspeed[time,i]
                                    #winddir_full[time,i] = winddir[time,i]
                                    #blh_full[time,i] = blh[time,i]

#### - Writing netCDF file with the full values - ###
if os.path.exists(outputfilename):
	os.remove(outputfilename)
f = netcdf.netcdf_file(outputfilename, 'w')
f.history = 'Created for a test'
f.createDimension('n_records', len(n_records))
f.createDimension('n_stations', len(n_stations))
surf_grd_nc = f.createVariable('surf_grd', 'float', ('n_records','n_stations'))
surf_grd_nc[:] = surf_grd_full
surf_grd_nc.units = 'station_fh.variables[surf_grd].Unit'
surf_slp_nc = f.createVariable('surf_slp', 'float', ('n_records','n_stations'))
surf_slp_nc[:] = surf_slp_full
surf_slp_nc.units = 'station_fh.variables[surf_slp].Unit'
surf_int_nc = f.createVariable('surf_int', 'float', ('n_records','n_stations'))
surf_int_nc[:] = surf_int_full
surf_int_nc.units = 'station_fh.variables[surf_int].Unit'
land_grd_nc = f.createVariable('land_grd', 'float', ('n_records','n_stations'))
land_grd_nc[:] = land_grd_full
land_grd_nc.units = 'station_fh.variables[land_grd].Unit'
land_slp_nc = f.createVariable('land_slp', 'float', ('n_records','n_stations'))
land_slp_nc[:] = land_slp_full
land_slp_nc.units = 'station_fh.variables[land_slp].Unit'
land_int_nc = f.createVariable('land_int', 'float', ('n_records','n_stations'))
land_int_nc[:] = land_int_full
land_int_nc.units = 'station_fh.variables[land_int].Unit'
nh_grd_nc = f.createVariable('nh_grd', 'float', ('n_records','n_stations'))
nh_grd_nc[:] = nh_grd_full
nh_grd_nc.units = 'station_fh.variables[nh_grd].Unit'
nh_slp_nc = f.createVariable('nh_slp', 'float', ('n_records','n_stations'))
nh_slp_nc[:] = nh_slp_full
nh_slp_nc.units = 'station_fh.variables[nh_slp].Unit'
nh_int_nc = f.createVariable('nh_int', 'float', ('n_records','n_stations'))
nh_int_nc[:] = nh_int_full
nh_int_nc.units = 'station_fh.variables[nh_int].Unit'
sh_grd_nc = f.createVariable('sh_grd', 'float', ('n_records','n_stations'))
sh_grd_nc[:] = sh_grd_full
sh_grd_nc.units = 'station_fh.variables[sh_grd].Unit'
sh_slp_nc = f.createVariable('sh_slp', 'float', ('n_records','n_stations'))
sh_slp_nc[:] = sh_slp_full
sh_slp_nc.units = 'station_fh.variables[sh_slp].Unit'
sh_int_nc = f.createVariable('sh_int', 'float', ('n_records','n_stations'))
sh_int_nc[:] = sh_int_full
sh_int_nc.units = 'station_fh.variables[sh_int].Unit'
strat_grd_nc = f.createVariable('strat_grd', 'float', ('n_records','n_stations'))
strat_grd_nc[:] = strat_grd_full
strat_grd_nc.units = 'station_fh.variables[strat_grd].Unit'
strat_slp_nc = f.createVariable('strat_slp', 'float', ('n_records','n_stations'))
strat_slp_nc[:] = strat_slp_full
strat_slp_nc.units = 'station_fh.variables[strat_slp].Unit'
strat_int_nc = f.createVariable('strat_int', 'float', ('n_records','n_stations'))
strat_int_nc[:] = strat_int_full
strat_int_nc.units = 'station_fh.variables[strat_int].Unit'
e90_grd_nc = f.createVariable('e90_grd', 'float', ('n_records','n_stations'))
e90_grd_nc[:] = e90_grd_full
e90_grd_nc.units = 'station_fh.variables[e90_grd].Unit'
e90_slp_nc = f.createVariable('e90_slp', 'float', ('n_records','n_stations'))
e90_slp_nc[:] = e90_slp_full
e90_slp_nc.units = 'station_fh.variables[e90_slp].Unit'
e90_int_nc = f.createVariable('e90_int', 'float', ('n_records','n_stations'))
e90_int_nc[:] = e90_int_full
e90_int_nc.units = 'station_fh.variables[e90_int].Unit'
# Temperature_nc = f.createVariable('Temperature', 'float', ('n_records','n_stations'))
# Temperature_nc[:] = Temperature_full
# Temperature_nc.units = 'station_fh.variables[Temperature].Unit'
# Pressure_nc = f.createVariable('Pressure', 'float', ('n_records','n_stations'))
# Pressure_nc[:] = Pressure_full
# Pressure_nc.units = 'station_fh.variables[Pressure].Unit'
# uwind_nc = f.createVariable('uwind', 'float', ('n_records','n_stations'))
# uwind_nc[:] = uwind_full
# uwind_nc.units = 'station_fh.variables[uwind].Unit'
# vwind_nc = f.createVariable('vwind', 'float', ('n_records','n_stations'))
# vwind_nc[:] = vwind_full
# vwind_nc.units = 'station_fh.variables[vwind].Unit'
# windspeed_nc = f.createVariable('windspeed', 'float', ('n_records','n_stations'))
# windspeed_nc[:] = windspeed_full
# windspeed_nc.units = 'station_fh.variables[windspeed].Unit'
# winddir_nc = f.createVariable('winddir', 'float', ('n_records','n_stations'))
# winddir_nc[:] = winddir_full
# winddir_nc.units = 'station_fh.variables[winddir].Unit'
# blh_nc = f.createVariable('blh', 'float', ('n_records','n_stations'))
# blh_nc[:] = blh_full
# blh_nc.units = 'station_fh.variables[blh].Unit'
f.close()

