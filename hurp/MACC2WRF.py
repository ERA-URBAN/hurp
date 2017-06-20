#!/usr/bin/env python

import numpy as np
import netCDF4 as nc
import csv
from hurp.RDConverter import * #RDWGSConverter
conv = RDWGSConverter()
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from sys import platform
from itertools import product
from os.path import exists
import cPickle
#from HURP_Toolbox import msg
import datetime

def msg(string,type=2):
    import sys

    if type == 1: print(string),
    if type == 2: print(string)
    
    sys.stdout.flush()
    return
  
def macc2wrf(MACCpath, wrfchemipath, wrfoutpath, Times, imon):
    FS               = 8 # FontSize
    domains          = ['d01','d02']  #['d04','d03','d02','d01']
    Nsnaps           = 13 
    zTimes           = ['00z','12z']
    tracers          = ['NOx','PM10','CO2']
    selectedtracers  = tracers[:2]
    dows             = ('Mon','Tue','Wed','Thu','Fri','Sat','Sun')
    moys             = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')


    emis_cat_names   = [['01: Energy industries                            '],   #  0        Stationary: 0,1,2,3,4,         10,11,12
                        ['02: Non-industrial combustion                    '],   #  1        Traffic:              5,6,7,8,9
                        ['34: Industry                                     '],   #  2
                        ['05: Fossil fuel production and distribution      '],   #  3 
                        ['06: Solvent and other product use                '],   #  4
                        ['71: Road transport, exhaust, gasoline            '],   #  5
                        ['72: Road transport, exhaust, diesel              '],   #  6
                        ['73: Road transport, exhaust, LPG and natural gas '],   #  7
                        ['74: Road transport, gasoline evaporation         '],   #  8
                        ['75: Road transport, tyre, brake and road wear    '],   #  9
                        ['08: Non-road transport                           '],   # 10
                        ['09: Waste                                        '],   # 11
                        ['10: Agriculture                                  '] ]  # 12     
                        

    #--- input
    if not 'progress' in globals(): progress = list()

    if not 'TimeProfilesLoaded' in progress:
        TPfilename    = '%s/time_profiles_2011.nc'%MACCpath
        TPfile        = nc.Dataset(TPfilename,'r')
        country_code  = TPfile.variables['country_code' ][:]
        country_name  = TPfile.variables['country_name' ][:]
        emis_cat      = TPfile.variables['emis_cat'     ][:]
        emis_cat_code = TPfile.variables['emis_cat_code'][:]
        emis_cat_name = TPfile.variables['emis_cat_name'][:]
        time          = TPfile.variables['time'         ][:]
        date          = TPfile.variables['date'         ][:]
        time_factors  = TPfile.variables['time_factors' ][:]
        TPfile.close()
        
        country_tmp        = country_code.tostring()
        country_tmp        = [country_tmp[i:i+3] for i in range(0, len(country_tmp), 3)]
        country_tmp        = list(enumerate(country_tmp))
        country_ids        = dict()
        for i in range(len(country_tmp)):    
            country_ids[country_tmp[i][1]] = country_tmp[i][0]
        ncountries         = len(country_ids.keys())

        MACCdow = list()
        for it in range(date.shape[0]):
            tt = datetime.datetime(date[it,0],date[it,1],date[it,2],date[it,3],date[it,4],date[it,5])
            MACCdow.append(tt.weekday()) # Mon = 0, ..., Sun = 6
        MACCdow = np.asarray(MACCdow)    
        progress.append('TimeProfilesLoaded')


    if not 'MACCloaded' in progress:
        #--- init
        msg('Loading gridded MACC files: ',type=1)
        E_macc                     = dict()
        for tracer in selectedtracers: 
            MACCncfilename         = '%s/MACC_%s.nc'%(MACCpath,tracer)
            msg('%s...'%MACCncfilename,type=1)
            MACCncfile             = nc.Dataset(MACCncfilename,'r')
            E_macc[tracer]         = MACCncfile.variables['E_%s'%tracer][:] # kg per area 
            if not 'lonmacc' in globals():
                lonmacc            = MACCncfile.variables['Longitude'][:]
                latmacc            = MACCncfile.variables['Latitude' ][:]
                country            = MACCncfile.variables['Country'  ][:]
                LONmacc,LATmacc    = np.meshgrid(lonmacc,latmacc)
                ny,nx              = LONmacc.shape
            MACCncfile.close()
        progress.append('MACCloaded')
    if not 'WRFOUTloaded' in progress:
      msg('Loading lat/lon from wrfout...',type=1)
      LATwrf = dict()
      LONwrf = dict()
      for domain in domains:
          #wrfoutfilename           = '%s/wrfout_%s_2015-06-30_00:00:00'%(wrfoutpath,domain)
          wrfoutfilename           = '%s/wrfinput_%s'%(wrfoutpath,domain)
          wrfoutfile               = nc.Dataset(wrfoutfilename,'r')       
          LATwrf[domain]           = wrfoutfile.variables['XLAT' ][0,:] # centerpoints
          LONwrf[domain]           = wrfoutfile.variables['XLONG'][0,:] # centerpoints
          wrfoutfile.close()
      progress.append('WRFOUTloaded')
      msg('done.')
                        
    if not 'MACCinterpolatedToWRF' in progress:
        Xwrf                        = dict()
        Ywrf                        = dict()
        E_wrf                       = dict()
        Country_wrf                 = dict()
        picklefilename              = '%s/E_MACC_wrf.pickle'%(wrfoutpath)
        if exists(picklefilename):
            msg('Loading pre-cooked %s'%picklefilename,type=1)
            pickledata              = cPickle.load(open(picklefilename,'r'))
            Xwrf,Ywrf,E_wrf,Country_wrf,domains,tracers = pickledata
            del(pickledata)
        else:
            for domain in domains:
                nlatwrf,nlonwrf     = LONwrf[domain].shape
                Country_wrf[domain] = griddata((LONmacc.flatten(),LATmacc.flatten()),country.flatten() , (LONwrf[domain],LATwrf[domain]),method='nearest')
                if domain == 'd01':
                    # Note: d04 has a higher resolution than NEI --? interpolate or nearest neighbor
                    #       d03 has the same resolution as NEI, but may be shifted in x or y direction --> nearest neighbour regridding
                    #       d02 has a similar resolution as MACC --> nearest neighbor interpolation with griddata
                    #       d01 has a coarser resolution than MACC  --> regrid cell by cell
                    #E_wrf[domain,tracer][isnap,:,:] = griddata((LONmacc.flatten(),LATmacc.flatten()),E_macc[tracer][isnap].flatten() , (LONwrf[domain],LATwrf[domain]),method='nearest')
                    msg('Interpolating %s..'%(domain),type=1)
                    for tracer in selectedtracers:
                        E_wrf[domain,tracer] = np.zeros((Nsnaps,nlatwrf,nlonwrf))
                    
                    Xbl                 = LONwrf[domain]*np.nan # boundary left
                    Xbr                 = LONwrf[domain]*np.nan # boundary right
                    Ix                  = np.arange(1,nlonwrf)
                    Xbl[:,Ix]           = (LONwrf[domain][:,Ix]+LONwrf[domain][:,Ix-1])/2.
                    Ix                  = np.arange(nlonwrf-1)
                    Xbr[:,Ix]           = (LONwrf[domain][:,Ix]+LONwrf[domain][:,Ix+1])/2.
                                        
                    Xbl[:, 0]           = LONwrf[domain][:, 0] - (Xbr[:, 0]-LONwrf[domain][:, 0])
                    Xbr[:,-1]           = LONwrf[domain][:,-1] + (LONwrf[domain][:,-1]-Xbl[:,-1])
                                        
                    Ybb                 = LATwrf[domain]*np.nan # boundary bottom
                    Ybt                 = LATwrf[domain]*np.nan # boundary top
                    Iy                  = np.arange(1,nlatwrf)
                    Ybb[Iy,:]           = (LATwrf[domain][Iy,:]+LATwrf[domain][Iy-1,:])/2.
                    Iy                  = np.arange(nlatwrf-1)
                    Ybt[Iy,:]           = (LATwrf[domain][Iy,:]+LATwrf[domain][Iy+1,:])/2.
                                        
                    Ybb[ 0,:]           = LATwrf[domain][ 0,:] - (Ybt[ 0,:]-LATwrf[domain][ 0,:])
                    Ybt[-1,:]           = LATwrf[domain][-1,:] + (LATwrf[domain][-1,:]-Ybb[-1,:])
                    
                    where               = np.where
                    t0                  = datetime.datetime.now()
                    for iy,ix in product(range(nlatwrf),range(nlonwrf)):
                        Iy,Ix           = where( (LONmacc > Xbl[iy,ix]) & (LONmacc <= Xbr[iy,ix]) & (LATmacc > Ybb[iy,ix]) & (LATmacc <= Ybt[iy,ix]) )
                        if len(Ix)>0:
                            for tracer in selectedtracers:
                                E_wrf[domain,tracer][:,iy,ix] = E_macc[tracer][:,Iy.min():Iy.max()+1,Ix.min():Ix.max()+1].mean(axis=(1,2))
                    
                        if ((np.remainder(iy,10) == 0) & (ix == nlonwrf-1)):
                            t1           = datetime.datetime.now()
                            dt           = t1-t0
                            msg('iy=%3d of %3d; %7.4f'%(iy,nlatwrf,dt.microseconds/1.e6)) # v0: 0.39 sec; v1 (Ilat outside loop): 0.007 sec
                            t0 = datetime.datetime.now()
                elif domain == 'd02':
                    msg('Interpolating %s..'%(domain),type=1)            
                    for tracer in selectedtracers:
                        E_wrf[domain,tracer] = np.zeros((Nsnaps,nlatwrf,nlonwrf))
                        for isnap in range(Nsnaps):
                            E_wrf[domain,tracer][isnap,:,:] = griddata((LONmacc.flatten(),LATmacc.flatten()),E_macc[tracer][isnap].flatten() , (LONwrf[domain],LATwrf[domain]),method='nearest')
                else: 
                    msg('domain not found: %s'%domain)
            # Save the data as a picklefile, so you don't have to the the interpolation each time.                    
            pickledata = (LONwrf,LATwrf,E_wrf,Country_wrf,domains,tracers)                
            cPickle.dump(pickledata,open(picklefilename,'wb'))   
            del(pickledata)        
        progress.append('MACCinterpolatedToWRF')
        msg('done.')
                        

                      
    #imon = 3 # Apr                    
    if not 'wrfchemi_written' in progress:
        for domain,zTime,dow,moy in product(domains,zTimes,dows,moys[imon-1:imon]):
            print('Writing wrfchemifile for %s %s %s %s'%(domain,zTime,moy,dow))
            zTimeOffset = 12 if zTime == '12z' else 0
        
            nlat,nlon                 = LATwrf[domain].shape 
            nz                        =  1

            TP                        = np.ones((24,Nsnaps,nlat,nlon))
            It                        = np.where( (date[:,1] == imon) & (MACCdow == dows.index(dow)) )
            It                        = It[0][:24]
            for ic in range(ncountries):
                I                     = np.where(Country_wrf[domain] == ic)
                Ilat                  = I[0]
                Ilon                  = I[1]
                TP_country            = np.tile(time_factors[It,:,ic,np.newaxis,np.newaxis],(1,1,nlat,nlon))
                TP[:,:,Ilat,Ilon]     = TP_country[:,:,Ilat,Ilon]
            
            wrfoutfilename            = '%s/wrfinput_%s'%(wrfoutpath,domain)
            wrfoutfile                = nc.Dataset(wrfoutfilename,'r')
            wrfchemifilename          = '%s/MACC_wrfchemi_%s_%s_%s_%s'%(wrfchemipath,zTime,domain,moy,dow)
            wrfchemifile              = nc.Dataset(wrfchemifilename, 'w',format = 'NETCDF3_CLASSIC' )

            wrfchemifile.createDimension('Time',        None) # will be 12
            wrfchemifile.createDimension('bottom_top' , nz)
            wrfchemifile.createDimension('south_north', nlat)
            wrfchemifile.createDimension('west_east'  , nlon)
            wrfchemifile.createDimension('DateStrLen' , 19)

        
            ncTimes                   = wrfchemifile.createVariable('Times'               ,'c',dimensions=('Time','DateStrLen'))
            ncXLAT                    = wrfchemifile.createVariable('XLAT'                ,'d',dimensions=('south_north','west_east'))
            ncXLONG                   = wrfchemifile.createVariable('XLONG'               ,'d',dimensions=('south_north','west_east'))
            ncE_NOx_stat              = wrfchemifile.createVariable('E_NOx_stat'          ,'d',dimensions=('Time','bottom_top','south_north','west_east'),fill_value = -1.e34)
            ncE_NOx_traf              = wrfchemifile.createVariable('E_NOx_traf'          ,'d',dimensions=('Time','bottom_top','south_north','west_east'),fill_value = -1.e34)
            ncE_PM10_stat             = wrfchemifile.createVariable('E_PM10_stat'         ,'d',dimensions=('Time','bottom_top','south_north','west_east'),fill_value = -1.e34)
            ncE_PM10_traf             = wrfchemifile.createVariable('E_PM10_traf'         ,'d',dimensions=('Time','bottom_top','south_north','west_east'),fill_value = -1.e34)
            ncTimes                   = Times[zTime]
            ncXLAT.FieldType          = 104.
            ncXLAT.MemoryOrder        = 'XY'
            ncXLAT.description        = 'LATITUDE, SOUTH IS NEGATIVE'
            ncXLAT.units              = 'degree_north'
            ncXLAT.stagger            = ' ' 
            ncXLAT.coordinates        = 'XLONG XLAT'        
            ncXLAT[:]                 = LATwrf[domain]
            ncXLONG.FieldType         = 104.
            ncXLONG.MemoryOrder       = 'XY'
            ncXLONG.description       = 'LONGITUDE, WEST IS NEGATIVE'
            ncXLONG.units             = 'degree_north'
            ncXLONG.stagger           = ' ' 
            ncXLONG.coordinates       = 'XLONG XLAT' 
            ncXLONG[:]                = LONwrf[domain]
            ncE_NOx_stat.FieldType    = 104.
            ncE_NOx_stat.MemoryOrder  = 'XYZ'
            ncE_NOx_stat.description  = 'NO2 emissions from stationary sources'
            ncE_NOx_stat.units        = 'mole/km2/hr'
            ncE_NOx_stat.stagger      = ' ' 
            ncE_NOx_stat.coordinates  = 'XLONG XLAT'   
            Isnap                     = np.asarray((0,1,2,3,4,10,11,12)) #         Stationary: 0,1,2,3,4,         10,11,12
            for it in range(12):     
                ncE_NOx_stat[it,0,:,:]= (TP[it+zTimeOffset,Isnap,:,:] * E_wrf[domain,'NOx'][Isnap,:,:]).sum(axis=0)
            ncE_NOx_traf.FieldType    = 104.
            ncE_NOx_traf.MemoryOrder  = 'XYZ'
            ncE_NOx_traf.description  = 'NO2 emissions from road traffic'
            ncE_NOx_traf.units        = 'mole/km2/hr'
            ncE_NOx_traf.stagger      = ' ' 
            ncE_NOx_traf.coordinates  = 'XLONG XLAT'   
            Isnap                     = np.asarray((5,6,7,8,9) )         # Traffic:              5,6,7,8,9
            for it in range(12):     
                ncE_NOx_traf[it,0,:,:] = (TP[it+zTimeOffset,Isnap,:,:]*E_wrf[domain,'NOx'][Isnap,:,:]).sum(axis=0)

            ncE_PM10_stat.FieldType   = 104.
            ncE_PM10_stat.MemoryOrder = 'XYZ'
            ncE_PM10_stat.description = 'PM10 emissions from stationary sources'
            ncE_PM10_stat.units       = 'kg/km2/hr'
            ncE_PM10_stat.stagger     = ' ' 
            ncE_PM10_stat.coordinates = 'XLONG XLAT'  
            Isnap                     = np.asarray((0,1,2,3,4,10,11,12)) #         Stationary: 0,1,2,3,4,         10,11,12
            for it in range(12):     
                ncE_PM10_stat[it,0,:,:]= (TP[it+zTimeOffset,Isnap,:,:] * E_wrf[domain,'PM10'][Isnap,:,:]).sum(axis=0)

            ncE_PM10_traf.FieldType   = 104.
            ncE_PM10_traf.MemoryOrder = 'XYZ'
            ncE_PM10_traf.description = 'PM10 emissions from road traffic'
            ncE_PM10_traf.units       = 'kg/km2/hr'
            ncE_PM10_traf.stagger     = ' ' 
            ncE_PM10_traf.coordinates = 'XLONG XLAT'  
            Isnap                     = np.asarray((5,6,7,8,9))          # Traffic:              5,6,7,8,9
            for it in range(12):     
                ncE_PM10_traf[it,0,:,:]   = (TP[it+zTimeOffset,Isnap,:,:]*E_wrf[domain,'PM10'][Isnap,:,:]).sum(axis=0)
                
            wrfchemifile.CEN_LAT      = wrfoutfile.getncattr('CEN_LAT')
            wrfchemifile.CEN_LON      = wrfoutfile.getncattr('CEN_LON')
            wrfchemifile.TRUELAT1     = wrfoutfile.getncattr('TRUELAT1')
            wrfchemifile.TRUELAT2     = wrfoutfile.getncattr('TRUELAT2') 
            wrfchemifile.MOAD_CEN_LAT = wrfoutfile.getncattr('MOAD_CEN_LAT')
            wrfchemifile.STAND_LON    = wrfoutfile.getncattr('STAND_LON')
            wrfchemifile.POLE_LAT     = wrfoutfile.getncattr('POLE_LAT') 
            wrfchemifile.POLE_LON     = wrfoutfile.getncattr('POLE_LON')
            wrfchemifile.GMT          = wrfoutfile.getncattr('GMT')  
            wrfchemifile.JULYR        = wrfoutfile.getncattr('JULYR')
            wrfchemifile.JULDAY       = wrfoutfile.getncattr('JULDAY')
            wrfchemifile.MAP_PROJ     = wrfoutfile.getncattr('MAP_PROJ')
            wrfchemifile.MMINLU       = wrfoutfile.getncattr('MMINLU')
            wrfchemifile.NUM_LAND_CAT = wrfoutfile.getncattr('NUM_LAND_CAT')
            wrfchemifile.ISWATER      = wrfoutfile.getncattr('ISWATER')
            wrfchemifile.ISLAKE       = wrfoutfile.getncattr('ISLAKE')   
            wrfchemifile.ISICE        = wrfoutfile.getncattr('ISICE')    
            wrfchemifile.ISURBAN      = wrfoutfile.getncattr('ISURBAN')     
            wrfchemifile.ISOILWATER   = wrfoutfile.getncattr('ISOILWATER')
            
            wrfchemifile.close()
        progress.append('wrfchemi_written')

if __name__=="__main__":
    if platform == 'win32':
        MACCpath     = 'U:/AMS_Stimulus_2016/Data_share/Emissions/TNO_MACC'
        wrfchemipath = 'U:/AMS_Stimulus_2016/Data_share/Emissions/wrfchemi'
        wrfoutpath   = 'S:/data/HURP/Adam_wrfout'
    elif (platform == 'darwin' and exists('/Users/molen050/mnt/promise')): # CapeGrim
        MACCpath     = '/Users/molen050/mnt/promise/WRF/michiel/HURP/Data_share/Emissions/TNO_MACC'
        wrfchemipath = '/Users/molen050/mnt/promise/WRF/michiel/HURP/Data_share/Emissions/wrfchemi'
        wrfoutpath   = '/Users/molen050/mnt/promise/WRF/michiel/HURP/Data_share/Domains'
    elif (platform == 'linux2' ): # Cartesius
        print('Setting env for Cartesius.')
        MACCpath     = '/projects/0/aams/wrfv3/Data_share/Emissions/TNO_MACC'
        wrfchemipath = '/projects/0/aams/wrfv3/Data_share/Emissions/wrfchemi'
        wrfoutpath   = '/projects/0/aams/wrfv3/OutputMichiel/HURP_20150630-20150704'
    else:  # maunaloa
        print('Setting path for Maunaloa')
        MACCpath     = '/Storage/WRF/michiel/HURP/Data_share/Emissions/TNO_MACC'
        wrfchemipath = '/Storage/WRF/michiel/HURP/Data_share/Emissions/wrfchemi'
        wrfoutpath   = '/Storage/WRF/michiel/HURP/Data_share/Domains'
        
    Times            = dict()
    Times['00z']     = [['2014-01-01_00:00:00'],['2014-01-01_01:00:00'],['2014-01-01_02:00:00'],['2014-01-01_03:00:00'],
                        ['2014-01-01_04:00:00'],['2014-01-01_05:00:00'],['2014-01-01_06:00:00'],['2014-01-01_07:00:00'],
                        ['2014-01-01_08:00:00'],['2014-01-01_09:00:00'],['2014-01-01_10:00:00'],['2014-01-01_11:00:00']]

    Times['12z']     = [['2014-01-01_12:00:00'],['2014-01-01_13:00:00'],['2014-01-01_14:00:00'],['2014-01-01_15:00:00'],
                        ['2014-01-01_16:00:00'],['2014-01-01_17:00:00'],['2014-01-01_18:00:00'],['2014-01-01_19:00:00'],
                        ['2014-01-01_20:00:00'],['2014-01-01_21:00:00'],['2014-01-01_22:00:00'],['2014-01-01_23:00:00']]
    macc2wrf(MACCpath, wrfchemipath, wrfoutpath, Times)
