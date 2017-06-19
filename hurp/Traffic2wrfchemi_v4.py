#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import os
from sys import platform

def traffic2wrfchemi(datapath, wrfchemipath, month, day):

  # SNAP                      1      2      3      4      5      6      7      8      9     10     11     12     13     14
  TP_moy = dict() # Month of Year
  TP_dow = dict() # Day of Week
  TP_hod = dict() # Hour of Day 
  TP_moy['Jan'] = np.array([1.000, 1.060, 0.950, 1.060, 1.000, 1.200, 0.880, 0.450, 1.000, 1.060, 1.200, 1.000, 0.880, 0.880])
  TP_moy['Feb'] = np.array([1.000, 1.047, 0.960, 1.047, 1.000, 1.150, 0.920, 1.300, 1.000, 1.047, 1.200, 1.000, 0.920, 0.920])
  TP_moy['Mar'] = np.array([1.000, 1.035, 1.020, 1.035, 1.000, 1.050, 0.980, 2.350, 1.000, 1.035, 1.200, 1.000, 0.980, 0.980])
  TP_moy['Apr'] = np.array([1.000, 1.010, 1.000, 1.010, 1.000, 1.000, 1.030, 1.700, 1.000, 1.010, 0.800, 1.000, 1.030, 1.030])
  TP_moy['May'] = np.array([1.000, 0.985, 1.010, 0.985, 1.000, 0.900, 1.050, 0.850, 1.000, 0.985, 0.800, 1.000, 1.050, 1.050])
  TP_moy['Jun'] = np.array([1.000, 0.960, 1.030, 0.960, 1.000, 0.850, 1.060, 0.850, 1.000, 0.960, 0.800, 1.000, 1.060, 1.060])
  TP_moy['Jul'] = np.array([1.000, 0.965, 1.030, 0.965, 1.000, 0.800, 1.010, 0.850, 1.000, 0.965, 0.800, 1.000, 1.010, 1.010])
  TP_moy['Aug'] = np.array([1.000, 0.935, 1.010, 0.935, 1.000, 0.875, 1.020, 1.000, 1.000, 0.935, 0.800, 1.000, 1.020, 1.020])
  TP_moy['Sep'] = np.array([1.000, 0.995, 1.040, 0.995, 1.000, 0.950, 1.060, 1.100, 1.000, 0.995, 0.800, 1.000, 1.060, 1.060])
  TP_moy['Oct'] = np.array([1.000, 1.010, 1.030, 1.010, 1.000, 1.000, 1.050, 0.650, 1.000, 1.010, 1.200, 1.000, 1.050, 1.050])
  TP_moy['Nov'] = np.array([1.000, 1.023, 1.010, 1.023, 1.000, 1.075, 1.010, 0.450, 1.000, 1.023, 1.200, 1.000, 1.010, 1.010])
  TP_moy['Dec'] = np.array([1.000, 0.975, 0.910, 0.975, 1.000, 1.150, 0.930, 0.450, 1.000, 0.975, 1.200, 1.000, 0.930, 0.930])

  TP_dow['Mon'] = np.array([1.000, 1.050, 1.200, 1.050, 1.000, 1.060, 1.020, 1.000, 1.000, 1.050, 1.000, 1.000, 1.000, 1.020])
  TP_dow['Tue'] = np.array([1.000, 1.050, 1.200, 1.050, 1.000, 1.060, 1.060, 1.000, 1.000, 1.050, 1.000, 1.000, 1.000, 1.060])
  TP_dow['Wed'] = np.array([1.000, 1.050, 1.200, 1.050, 1.000, 1.060, 1.080, 1.000, 1.000, 1.050, 1.000, 1.000, 1.000, 1.080])
  TP_dow['Thu'] = np.array([1.000, 1.050, 1.200, 1.050, 1.000, 1.060, 1.100, 1.000, 1.000, 1.050, 1.000, 1.000, 1.000, 1.100])
  TP_dow['Fri'] = np.array([1.000, 1.050, 1.200, 1.050, 1.000, 1.060, 1.140, 1.000, 1.000, 1.050, 1.000, 1.000, 1.000, 1.140])
  TP_dow['Sat'] = np.array([1.000, 0.910, 0.500, 0.910, 1.000, 0.850, 0.810, 1.000, 1.000, 0.910, 1.000, 1.000, 1.000, 0.810])
  TP_dow['Sun'] = np.array([1.000, 0.840, 0.500, 0.840, 1.000, 0.850, 0.790, 1.000, 1.000, 0.840, 1.000, 1.000, 1.000, 0.790])

  TP_hod['00']  = np.array([1.000, 0.875, 0.500, 0.875, 1.000, 0.790, 0.190, 0.600, 1.000, 0.875, 1.000, 1.000, 1.000, 0.190])
  TP_hod['01']  = np.array([1.000, 0.875, 0.350, 0.875, 1.000, 0.720, 0.090, 0.600, 1.000, 0.875, 1.000, 1.000, 1.000, 0.090])
  TP_hod['02']  = np.array([1.000, 0.890, 0.200, 0.890, 1.000, 0.720, 0.060, 0.600, 1.000, 0.890, 1.000, 1.000, 1.000, 0.060])
  TP_hod['03']  = np.array([1.000, 0.910, 0.100, 0.910, 1.000, 0.710, 0.050, 0.600, 1.000, 0.910, 1.000, 1.000, 1.000, 0.050])
  TP_hod['04']  = np.array([1.000, 0.940, 0.100, 0.940, 1.000, 0.740, 0.090, 0.600, 1.000, 0.940, 1.000, 1.000, 1.000, 0.090])
  TP_hod['05']  = np.array([1.000, 0.975, 0.200, 0.975, 1.000, 0.800, 0.220, 0.650, 1.000, 0.975, 1.000, 1.000, 1.000, 0.220])
  TP_hod['06']  = np.array([1.000, 1.010, 0.750, 1.010, 1.000, 0.920, 0.860, 0.750, 1.000, 1.010, 1.000, 1.000, 1.000, 0.860])
  TP_hod['07']  = np.array([1.000, 1.045, 1.250, 1.045, 1.000, 1.080, 1.840, 0.900, 1.000, 1.045, 1.000, 1.000, 1.000, 1.840])
  TP_hod['08']  = np.array([1.000, 1.080, 1.400, 1.080, 1.000, 1.190, 1.860, 1.100, 1.000, 1.080, 1.000, 1.000, 1.000, 1.860])
  TP_hod['09']  = np.array([1.000, 1.110, 1.500, 1.110, 1.000, 1.220, 1.410, 1.350, 1.000, 1.110, 1.000, 1.000, 1.000, 1.410])
  TP_hod['10']  = np.array([1.000, 1.140, 1.500, 1.140, 1.000, 1.210, 1.240, 1.450, 1.000, 1.140, 1.000, 1.000, 1.000, 1.240])
  TP_hod['11']  = np.array([1.000, 1.150, 1.500, 1.150, 1.000, 1.210, 1.200, 1.600, 1.000, 1.150, 1.000, 1.000, 1.000, 1.200])
  TP_hod['12']  = np.array([1.000, 1.110, 1.500, 1.110, 1.000, 1.170, 1.320, 1.650, 1.000, 1.110, 1.000, 1.000, 1.000, 1.320])
  TP_hod['13']  = np.array([1.000, 1.120, 1.500, 1.120, 1.000, 1.150, 1.440, 1.750, 1.000, 1.120, 1.000, 1.000, 1.000, 1.440])
  TP_hod['14']  = np.array([1.000, 1.125, 1.500, 1.125, 1.000, 1.140, 1.450, 1.700, 1.000, 1.125, 1.000, 1.000, 1.000, 1.450])
  TP_hod['15']  = np.array([1.000, 1.080, 1.500, 1.080, 1.000, 1.130, 1.590, 1.550, 1.000, 1.080, 1.000, 1.000, 1.000, 1.590])
  TP_hod['16']  = np.array([1.000, 1.040, 1.500, 1.040, 1.000, 1.100, 2.030, 1.350, 1.000, 1.040, 1.000, 1.000, 1.000, 2.030])
  TP_hod['17']  = np.array([1.000, 1.005, 1.400, 1.005, 1.000, 1.070, 2.080, 1.100, 1.000, 1.005, 1.000, 1.000, 1.000, 2.080])
  TP_hod['18']  = np.array([1.000, 0.975, 1.250, 0.975, 1.000, 1.040, 1.510, 0.900, 1.000, 0.975, 1.000, 1.000, 1.000, 1.510])
  TP_hod['19']  = np.array([1.000, 0.950, 1.100, 0.950, 1.000, 1.020, 1.060, 0.750, 1.000, 0.950, 1.000, 1.000, 1.000, 1.060])
  TP_hod['20']  = np.array([1.000, 0.925, 1.000, 0.925, 1.000, 1.020, 0.740, 0.650, 1.000, 0.925, 1.000, 1.000, 1.000, 0.740])
  TP_hod['21']  = np.array([1.000, 0.905, 0.900, 0.905, 1.000, 1.010, 0.620, 0.600, 1.000, 0.905, 1.000, 1.000, 1.000, 0.620])
  TP_hod['22']  = np.array([1.000, 0.890, 0.800, 0.890, 1.000, 0.960, 0.610, 0.600, 1.000, 0.890, 1.000, 1.000, 1.000, 0.610])
  TP_hod['23']  = np.array([1.000, 0.875, 0.700, 0.875, 1.000, 0.880, 0.440, 0.600, 1.000, 0.875, 1.000, 1.000, 1.000, 0.440])

  # define months/days of week
  moys             = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
  dows             = ('Mon','Tue','Wed','Thu','Fri','Sat','Sun')

  infilename1                  = '%s/resulting_intensities_onGrid_binnenstad.csv'%datapath
  infilename2                  = '%s/resulting_intensities_onGrid_snelweg.csv'%datapath
  infilename3                  = '%s/NEI_wrfchemi_00z_d04_%s_%s'%(wrfchemipath,moys[month-1], dows[day])

  if not 'progress' in globals(): progress = list()
  if not 'dataloaded' in progress:
      import pdb; pdb.set_trace()
      #--- load grid from existing wrfchemi file
      print('Loading grid data from %s'%infilename3)
      ncfile                   = nc.Dataset(infilename3,'r')
      xlon                     = ncfile.variables['XLONG'][:]
      xlat                     = ncfile.variables['XLAT'] [:]
      ncfile.close()           
                              
      # --- initialise traffic intensity matrices                           
      I_city_L                 = xlat * 0 # veh/24h
      I_city_M                 = xlat * 0 # veh/24h
      I_city_H                 = xlat * 0 # veh/24h
      I_ring_L                 = xlat * 0 # veh/24h
      I_ring_M                 = xlat * 0 # veh/24h
      I_ring_H                 = xlat * 0 # veh/24h
      
      #--- load traffic intensities for city    
      print('Loading traffic intensity data from %s'%infilename1)
      infile                   = open(infilename1,'r')
      iline = -1               
      for line in infile:      
          iline                = iline+1
          if iline > 1:
              #,lat,lon,gridId,STRAATNAAM,KM_LV24_x,KM_MV24_x,KM_ZV24_x
              words            = line.split(',')
              if not words[5]: words[5] = 0
              if not words[6]: words[6] = 0
              if(not words[7] or (words[7] == '\r\n')): words[7] = 0
              order            = np.int  (words[0])
              lat              = np.float(words[1])
              lon              = np.float(words[2])
              gridID           = np.int  (words[3])
              NAME             =          words[4] 
              LVDU_x           = np.float(words[5]) # light  veh km's driven in grid cell in 24h
              MVDU_x           = np.float(words[6]) # middle veh km's driven in grid cell in 24h
              ZVDU_x           = np.float(words[7]) # heavy veh km's driven in grid cell in 24h
              
              distance         = np.square(xlon-lon) + np.square(xlat-lat)
              iy,ix            = np.where( distance == distance.min() )
            # if distance[iy,ix] > 1e-10: 
            #     print('(lat,lon) = (%12.10f,%12.10f),(xlat,xlon) = (%12.10f,%12.10f)'%(lat,lon,xlat[iy,ix],xlon[iy,ix]))
            #     print('(dlat,dlon)=(%f,%f)'%(xlat[iy,ix]-lat,xlon[iy,ix]-lon))
              
              if I_city_L[iy,ix] > 0: print('warning: I_city_L[%3d,%3d] is not empty!'%(iy,ix))
              if I_city_M[iy,ix] > 0: print('warning: I_city_M[%3d,%3d] is not empty!'%(iy,ix))
              if I_city_H[iy,ix] > 0: print('warning: I_city_H[%3d,%3d] is not empty!'%(iy,ix))
              I_city_L[iy,ix]    = I_city_L[iy,ix] + LVDU_x
              I_city_M[iy,ix]    = I_city_M[iy,ix] + MVDU_x
              I_city_H[iy,ix]    = I_city_H[iy,ix] + ZVDU_x
            
      infile.close()
      words = ''
      
      #--- load traffic intensities for ring
      print('Loading traffic intensity data from %s'%infilename2)
      infile                   = open(infilename2,'r')
      iline = -1               
      for line in infile:      
          iline                = iline+1
          if iline > 1:
              #city: ,lat,lon,gridId,NAME,      LVDU_x,MVDU_x,ZVDU_x,Shape_Length_x
              #ring: ,lat,lon,gridId,STRAATNAAM,LVDU_x,MVDU_x,ZVDU_x
              words            = line.split(',')
              print words[7]  
              words[7]         = words[7].rstrip()
              if   not words[5]      : words[5] = 0
              if   not words[6]      : words[6] = 0
              if  (not words[7] or (words[7] == '\r\n')): words[7] = 0

              order            = np.int  (words[0])
              lat              = np.float(words[1])
              lon              = np.float(words[2])
              gridID           = np.int  (words[3])
              NAME             =          words[4] 
              LVDU_x           = np.float(words[5]) # light  veh km's driven in grid cell in 24h
              MVDU_x           = np.float(words[6]) # middle veh km's driven in grid cell in 24h
              ZVDU_x           = np.float(words[7]) # heavy veh km's driven in grid cell in 24h
                          
              distance         = np.square(xlon-lon) + np.square(xlat-lat)
              iy,ix            = np.where( distance == distance.min() )
            # if distance[iy,ix] > 1e-10: 
            #     print('(lat,lon) = (%12.10f,%12.10f),(xlat,xlon) = (%12.10f,%12.10f)'%(lat,lon,xlat[iy,ix],xlon[iy,ix]))
            #     print('(dlat,dlon)=(%f,%f)'%(xlat[iy,ix]-lat,xlon[iy,ix]-lon))
              
              if I_ring_L[iy,ix] > 0: print('warning: I_ring_L[%3d,%3d] is not empty!'%(iy,ix))
              if I_ring_M[iy,ix] > 0: print('warning: I_ring_M[%3d,%3d] is not empty!'%(iy,ix))
              if I_ring_H[iy,ix] > 0: print('warning: I_ring_H[%3d,%3d] is not empty!'%(iy,ix))
              I_ring_L[iy,ix]    = I_ring_L[iy,ix] + LVDU_x
              I_ring_M[iy,ix]    = I_ring_M[iy,ix] + MVDU_x
              I_ring_H[iy,ix]    = I_ring_H[iy,ix] + ZVDU_x
      infile.close()

      progress.append('dataloaded')

  #--- Emissions
  g2mole_NOx     = 1./46.             # M_NO2 = 46 g/mole
  g2kg           = 0.001
  day2hr         = 1/24.
  A              = 0.1 * 0.1          # km2
  f_city_L_NOx   = 0.35  * g2mole_NOx # g/km/veh NOx in NO2-equivalenten;     Stad-normaal               2017; licht       verkeer
  f_city_M_NOx   = 5.63  * g2mole_NOx # g/km/veh NOx in NO2-equivalenten;     Stad-normaal               2017; middelzwaar verkeer
  f_city_H_NOx   = 6.78  * g2mole_NOx # g/km/veh NOx in NO2-equivalenten;     Stad-normaal               2017; zwaar       verkeer
  f_city_L_PM10  = 0.036 * g2kg       # g/km/veh PM10 (verbranding+slijtage); Stad-normaal               2017; licht       verkeer 
  f_city_M_PM10  = 0.175 * g2kg       # g/km/veh PM10 (verbranding+slijtage); Stad-normaal               2017; middelzwaar verkeer 
  f_city_H_PM10  = 0.183 * g2kg       # g/km/veh PM10 (verbranding+slijtage); Stad-normaal               2017; zwaar       verkeer 
  f_ring_L_NOx   = 0.27  * g2mole_NOx # g/km/veh NOx in NO2-equivalenten;     Snelweg-vrije doorstroming 2017; licht       verkeer
  f_ring_M_NOx   = 2.42  * g2mole_NOx # g/km/veh NOx in NO2-equivalenten;     Snelweg-vrije doorstroming 2017; middelzwaar verkeer
  f_ring_H_NOx   = 2.41  * g2mole_NOx # g/km/veh NOx in NO2-equivalenten;     Snelweg-vrije doorstroming 2017; zwaar       verkeer
  f_ring_L_PM10  = 0.024 * g2kg       # g/km/veh PM10 (verbranding+slijtage); Snelweg-vrije doorstroming 2017; licht       verkeer 
  f_ring_M_PM10  = 0.097 * g2kg       # g/km/veh PM10 (verbranding+slijtage); Snelweg-vrije doorstroming 2017; middelzwaar verkeer 
  f_ring_H_PM10  = 0.090 * g2kg       # g/km/veh PM10 (verbranding+slijtage); Snelweg-vrije doorstroming 2017; zwaar       verkeer 


  E_NOx_traf     = day2hr/A * (I_city_L * f_city_L_NOx  + I_city_M * f_city_M_NOx  + I_city_H * f_city_H_NOx  + \
                              I_ring_L * f_ring_L_NOx  + I_ring_M * f_ring_M_NOx  + I_ring_H * f_ring_H_NOx)      # veh-km/day * days2hr * g/km/veh * g2mole / km2 = mole/km2/hr
                  
  E_PM10_traf    = day2hr/A * (I_city_L * f_city_L_PM10 + I_city_M * f_city_M_PM10 + I_city_H * f_city_H_PM10 + \
                              I_ring_L * f_ring_L_PM10 + I_ring_M * f_ring_M_PM10 + I_ring_H * f_ring_H_PM10)     # veh-km/day * days2hr * g/km/veh * g2kg   / km2 =   kg/km2/hr
      

  if not 'wrfchemi_written' in progress:
      ncinfilenames = [os.path.join(wrfchemipath,filename) for filename in os.listdir(wrfchemipath) if (filename.startswith('NEI_wrfchemi') and 'd04' in filename)]
      for ncinfilename in ncinfilenames:
          zTime                     = ncinfilename[-15:-12]
          domain                    = 'd04'
          moy                       = ncinfilename[ -7: -4]
          dow                       = ncinfilename[ -3:   ]
          zTimeOffset               = 12 if zTime == '12z' else 0 
          
          ncoutfilename             = '%s/NEI_Traffic_wrfchemi_%s_%s_%s_%s'%(wrfchemipath,zTime,domain,moy,dow)
                                    
          ncinfile                  = nc.Dataset(ncinfilename,'r')
          Times                     = ncinfile.variables['Times'      ][:]
          xlat                      = ncinfile.variables['XLAT'       ][:]
          xlon                      = ncinfile.variables['XLONG'      ][:]
          E_NOx_stat                = ncinfile.variables['E_NOx_stat' ][:]
          E_PM10_stat               = ncinfile.variables['E_PM10_stat'][:]

          
          ncoutfile                 = nc.Dataset(ncoutfilename,'w')
          nz                        =  1
          nlat,nlon                 = xlat.shape
      
          ncoutfile.createDimension('Time',        None) # will be 12
          ncoutfile.createDimension('bottom_top' , nz)
          ncoutfile.createDimension('south_north', nlat)
          ncoutfile.createDimension('west_east'  , nlon)
          ncoutfile.createDimension('DateStrLen' , 19)
      
          ncTimes                   = ncoutfile.createVariable('Times'               ,'c',dimensions=('Time','DateStrLen'))
          ncXLAT                    = ncoutfile.createVariable('XLAT'                ,'d',dimensions=('south_north','west_east'))
          ncXLONG                   = ncoutfile.createVariable('XLONG'               ,'d',dimensions=('south_north','west_east'))
          ncE_NOx_stat              = ncoutfile.createVariable('E_NOx_stat'          ,'d',dimensions=('Time','bottom_top','south_north','west_east'),fill_value = -1.e34)
          ncE_NOx_traf              = ncoutfile.createVariable('E_NOx_traf'          ,'d',dimensions=('Time','bottom_top','south_north','west_east'),fill_value = -1.e34)
          ncE_PM10_stat             = ncoutfile.createVariable('E_PM10_stat'         ,'d',dimensions=('Time','bottom_top','south_north','west_east'),fill_value = -1.e34)
          ncE_PM10_traf             = ncoutfile.createVariable('E_PM10_traf'         ,'d',dimensions=('Time','bottom_top','south_north','west_east'),fill_value = -1.e34)
      
          ncTimes[:]                = Times
                                      
          ncXLAT.FieldType          = 104.
          ncXLAT.MemoryOrder        = 'XY'
          ncXLAT.description        = 'LATITUDE, SOUTH IS NEGATIVE'
          ncXLAT.units              = 'degree_north'
          ncXLAT.stagger            = ' ' 
          ncXLAT.coordinates        = 'XLONG XLAT'        
          ncXLAT[:]                 = xlat
          
          ncXLONG.FieldType         = 104.
          ncXLONG.MemoryOrder       = 'XY'
          ncXLONG.description       = 'LONGGITUDE, WEST IS NEGATIVE'
          ncXLONG.units             = 'degree_north'
          ncXLONG.stagger           = ' ' 
          ncXLONG.coordinates       = 'XLONG XLAT' 
          ncXLONG[:]                = xlon
      
          ncE_NOx_stat.FieldType    = 104.
          ncE_NOx_stat.MemoryOrder  = 'XYZ'
          ncE_NOx_stat.description  = 'NOx emissions from stationary sources'
          ncE_NOx_stat.units        = 'mole/km2/hr'
          ncE_NOx_stat.stagger      = ' ' 
          ncE_NOx_stat.coordinates  = 'XLONG XLAT'   
          ncE_NOx_stat[:]           = E_NOx_stat
      
          ncE_NOx_traf.FieldType    = 104.
          ncE_NOx_traf.MemoryOrder  = 'XYZ'
          ncE_NOx_traf.description  = 'NOx emissions from road traffic'
          ncE_NOx_traf.units        = 'mole/km2/hr'
          ncE_NOx_traf.stagger      = ' ' 
          ncE_NOx_traf.coordinates  = 'XLONG XLAT'   
          Isnap                     = 13
          for it in range(12):     
              TP                    = TP_hod['%02d'%(it+zTimeOffset)][Isnap] * TP_dow[dow][Isnap] * TP_moy[moy][Isnap]
              ncE_NOx_traf[it,0,:,:] = TP*E_NOx_traf
      #      TP2                   = np.tile(TP[...,np.newaxis,np.newaxis],(1,nlat,nlon))            
      
          ncE_PM10_stat.FieldType   = 104.
          ncE_PM10_stat.MemoryOrder = 'XYZ'
          ncE_PM10_stat.description = 'PM10 emissions from stationary sources'
          ncE_PM10_stat.units       = 'kg/km2/hr'
          ncE_PM10_stat.stagger     = ' ' 
          ncE_PM10_stat.coordinates = 'XLONG XLAT'  
          ncE_PM10_stat[:]          = E_PM10_stat
      
          ncE_PM10_traf.FieldType   = 104.
          ncE_PM10_traf.MemoryOrder = 'XYZ'
          ncE_PM10_traf.description = 'PM10 emissions from road traffic'
          ncE_PM10_traf.units       = 'kg/km2/hr'
          ncE_PM10_traf.stagger     = ' ' 
          ncE_PM10_traf.coordinates = 'XLONG XLAT'  
          Isnap                     = 13
          for it in range(12):     
              TP                    = TP_hod['%02d'%(it+zTimeOffset)][Isnap] * TP_dow[dow][Isnap] * TP_moy[moy][Isnap]
      #       TP2                   = np.tile(TP[...,np.newaxis,np.newaxis],(1,nlat,nlon))
              ncE_PM10_traf[it,0,:,:]   = TP*E_PM10_traf
              
          ncoutfile.CEN_LAT      = ncinfile.getncattr('CEN_LAT')
          ncoutfile.CEN_LON      = ncinfile.getncattr('CEN_LON')
          ncoutfile.TRUELAT1     = ncinfile.getncattr('TRUELAT1')
          ncoutfile.TRUELAT2     = ncinfile.getncattr('TRUELAT2') 
          ncoutfile.MOAD_CEN_LAT = ncinfile.getncattr('MOAD_CEN_LAT')
          ncoutfile.STAND_LON    = ncinfile.getncattr('STAND_LON')
          ncoutfile.POLE_LAT     = ncinfile.getncattr('POLE_LAT') 
          ncoutfile.POLE_LON     = ncinfile.getncattr('POLE_LON')
          ncoutfile.GMT          = ncinfile.getncattr('GMT')  
          ncoutfile.JULYR        = ncinfile.getncattr('JULYR')
          ncoutfile.JULDAY       = ncinfile.getncattr('JULDAY')
          ncoutfile.MAP_PROJ     = ncinfile.getncattr('MAP_PROJ')
          ncoutfile.MMINLU       = ncinfile.getncattr('MMINLU')
          ncoutfile.NUM_LAND_CAT = ncinfile.getncattr('NUM_LAND_CAT')
          ncoutfile.ISWATER      = ncinfile.getncattr('ISWATER')
          ncoutfile.ISLAKE       = ncinfile.getncattr('ISLAKE')   
          ncoutfile.ISICE        = ncinfile.getncattr('ISICE')    
          ncoutfile.ISURBAN      = ncinfile.getncattr('ISURBAN')     
          ncoutfile.ISOILWATER   = ncinfile.getncattr('ISOILWATER')
          
          ncinfile.close()
          ncoutfile.close()
          progress.append('wrfchemi_written') 

  if False:        
      f = plt.figure(1)
      f.clf()
      vmax = 2000.
      ax = f.add_subplot(231)
      ax.pcolor(xlon,xlat,I_city_L,vmax=vmax)
      ax.set_title('I city light')
      
      ax = f.add_subplot(232)
      ax.pcolor(xlon,xlat,I_city_M,vmax=vmax)
      ax.set_title('I city middle')
      
      ax = f.add_subplot(233)
      ax.pcolor(xlon,xlat,I_city_H,vmax=vmax)
      ax.set_title('I city heavy')
      
      ax = f.add_subplot(234)
      ax.pcolor(xlon,xlat,I_ring_L,vmax=vmax)
      ax.set_title('I ring light')
      
      ax = f.add_subplot(235)
      ax.pcolor(xlon,xlat,I_ring_M,vmax=vmax)
      ax.set_title('I ring middle')
      
      ax = f.add_subplot(236)
      ax.pcolor(xlon,xlat,I_ring_H,vmax=vmax)
      ax.set_title('I ring heavy')
      
      f = plt.figure(2)
      f.clf()
      
      ax = f.add_subplot(121)
      ax.pcolor(xlon,xlat,E_NOx)
      ax.set_title('E_NOx')
      
      ax = f.add_subplot(122)
      ax.pcolor(xlon,xlat,E_PM10)
      ax.set_title('E_PM10')
      
      f = plt.figure(3)
      f.clf()
      ax = f.add_subplot(111)
      ax.pcolor(xlon,xlat,L_city)
      
if __name__=="__main__":
    month = 'Mar'

    if   platform == 'win32':
        trafficpath                     = 'u:\AMS_Stimulus_2016\Data_share\Verkeersintensiteiten per grid ring en binnenstad'
        wrfchemipath                 = 'u:\AMS_Stimulus_2016\Data_share\Emissions\wrfchemi'
    elif platform == 'linux2':
        trafficpath                     = '/projects/0/aams/wrfv3/Data_share/Verkeersintensiteiten per grid ring en binnenstad'
        wrfchemipath                 = '/projects/0/aams/wrfv3/Data_share/Emissions/wrfchemi'

    traffic2wrfchemi(trafficpath, wrfchemipath, month)
