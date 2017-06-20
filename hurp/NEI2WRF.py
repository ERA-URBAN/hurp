import numpy as np
import netCDF4 as nc
import csv
from hurp.RDConverter import * #RDWGSConverter
conv = RDWGSConverter()
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from sys import platform
#from mpl_toolkits.basemap import interp as mapinterp
from itertools import product
from os.path import exists
import cPickle

def msg(string,type=2):
    import sys

    if type == 1: print(string),
    if type == 2: print(string)
    
    sys.stdout.flush()
    return

def nei2wrf(NEIpath, wrfchemipath, wrfoutpath, Times, month):
    '''
    The purpose of this script is to write WRF-Chem emission files for HURP.
    There are three input sources:
    - NEI - Netherlands Emission Inventory: a list of area sources of CO2, NOx and PM10, with RD coordinates
            This data needs to be gridded.
    - wrfout - contains the information of the domains in HURP, one wrfout file for each domain. The coordinates are in WGS84 lat/lon, and spaced at fixed distance (km) intervals.

    Four wrfchemi output files will be written, one for each domain, containing the gridded emission data projected on the WRF grids.
    '''
    #--- settings
    FS               = 8 # FontSize
    domains          = ['d04','d03','d02','d01']
    kg2mol           = dict()
    kg2mol['CO2']    = 1000./44. # 1 kmol = 44 kg       | 1/44 kmol = 1 kg | 1 kg = 1000/44 mol
    kg2mol['NOx']    = 1000./30.8 # 1 kmol = 44 kg       | 1/44 kmol = 1 kg | 1 kg = 1000/44 mol
    kg2mol['PM10']   = 1.
    yr2hr            = 1./(365.*24.)
    Nsnaps           = 14 # 13 SNAPS + 1 for road traffic
    zTimes           = ['00z','12z']
    tracers          = ['NOx','PM10','CO2']
    selectedtracers  = tracers[:2]
    dows             = ('Sun', 'Mon','Tue','Wed','Thu','Fri','Sat')
    moys             = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

                                    # in NER output file                                 # In Ingrid Super's file
    SNAPnames                    = [ [ '01: Afvalverwijdering                       '  ],   #SNAP9='waste treatment and disposal';
                                    [ '02: Bouw                                    '  ],   # 
                                    [ '03: Chemische Industrie                     '  ],   #SNAP6='solvents and other product use';
                                    [ '04: Consumenten                             '  ],   # 
                                    [ '05: Drinkwatervoorziening                   '  ],   # 
                                    [ '06: Energiesector                           '  ],   #SNAP1='combustion in energy and transformation industries';
                                    [ '07: Handel, Diensten en Overheid (HDO)      '  ],   # SNAP2='non-industrial combustion plants';
                                    [ '08: Landbouw                                '  ],   #SNAP10='agriculture'
                                    [ '09: Natuur                                  '  ],   # 
                                    [ '10: Overige industrie                       '  ],   #SNAP3='combustion in manufacturing industry'; SNAP4='production processes'
                                    [ '11: Raffinaderijen                          '  ],   #SNAP5='extraction and distribution of fossil fuels and geothermal energy';
                                    [ '12: Riolering en waterzuiveringsinstallaties'  ],   #
                                    [ '13: Verkeer en vervoer - Except road traffic'  ],   #SNAP8='other mobile sources and machinery'
                                    [ '14: Verkeer en vervoer - Wegverkeer         '  ] ]  #SNAP7='road transport';
                            
                        
    # 01: Afvalverwerking                          = SNAP  9 'waste treatment and disposal';
    # 02: Bouw                                     = SNAP  3 'combustion in manufacturing industry'; SNAP4='production processes' : SNAP 3 & 4 combined                                                                                                                                     
    # 03: Chemische industrie                      = SNAP  6 'solvents and other product use';
    # 04: Consumenten                              = SNAP  3 'combustion in manufacturing industry'; SNAP4='production processes' : SNAP 3 & 4 combined                                                                                                                                   
    # 05: Drinkwatervoorziening                    = Constant genomen
    # 06: Energiesector                            = SNAP  1 'combustion in energy and transformation industries';
    # 07: Handel, Diensten en Overheid (HDO)       = SNAP 7 'road transport'
    # 08: Landbouw                                 = SNAP 10 'agriculture'             
    # 09: Natuur                                   = Constant genomen                                                                                                                      
    # 10: Overige industrie                        = SNAP3='combustion in manufacturing industry'; SNAP4='production processes' : SNAP 3 & 4 combined
    # 11: Raffinaderijen                           = SNAP  5 'extraction and distribution of fossil fuels and geothermal energy';                  
    # 12: Riolering en waterzuiveringsinstallaties = Constant genomen
    # 13: Verkeer en vervoer - Except road traffic = SNAP  8 'other mobile sources and machinery'
    # 14: Verkeer en vervoer - Wegverkeer          = SNAP  7 'road transport'        
    TP_moy = dict() # Month of Year
    TP_dow = dict() # Day of Week
    TP_hod = dict() # Hour of Day

    # SNAP                      1      2      3      4      5      6      7      8      9     10     11     12     13     14
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
                        
    #--- input
    if not 'progress' in globals(): progress = list()
    if not 'NEIloaded' in progress:
        #--- init
        msg('Loading gridded NEI files: ',type=1)
        E_nei             = dict()
        for tracer in selectedtracers: 
            NEIncfilename = '%s/NEI_%s.nc'%(NEIpath,tracer)
            msg('%s...'%NEIncfilename,type=1)
            NEIncfile     = nc.Dataset(NEIncfilename,'r')
            E_nei[tracer] = NEIncfile.variables['E_%s'%tracer][:] # kg per area 
            if not 'xnei' in globals():
                xnei      = NEIncfile.variables['Longitude'][:]
                ynei      = NEIncfile.variables['Latitude' ][:]
                Xnei,Ynei = np.meshgrid(xnei,ynei)
                ny,nx     = Xnei.shape
            NEIncfile.close()
        progress.append('NEIloaded')
        
    if not 'WRFOUTloaded' in progress:
        msg('Loading lat/lon from wrfout...',type=1)
        xlat = dict()
        xlon = dict()
        for domain in domains:
            wrfoutfilename = '%s/wrfinput_%s'%(wrfoutpath,domain)
            wrfoutfile     = nc.Dataset(wrfoutfilename,'r')       
            xlat[domain]   = wrfoutfile.variables['XLAT' ][0,:] # centerpoints
            xlon[domain]   = wrfoutfile.variables['XLONG'][0,:] # centerpoints
            wrfoutfile.close()
        progress.append('WRFOUTloaded')
        msg('done.')
        
        
    if not 'NEIinterpolatedToWRF' in progress:
        Xwrf             = dict()
        Ywrf             = dict()
        E_wrf            = dict()
        for tracer in selectedtracers:
            picklefilename = '%s/E_wrf_%s.pickle'%(wrfoutpath,tracer)
            if exists(picklefilename):
                msg('Loading pre-cooked %s'%picklefilename,type=1)
                pickledata = cPickle.load(open(picklefilename,'r'))
                Xwrf,Ywrf,E_wrf_tmp,domains,tracers = pickledata
                for domain in domains:
                    E_wrf[domain,tracer] = E_wrf_tmp[domain,tracer]
                del(E_wrf_tmp)
            else:

                # Note: d04 has a higher resolution than NEI --? interpolate or nearest neighbor
                #       d03 has the same resolution as NEI, but may be shifted in x or y direction --> nearest neighbour regridding
                #       d02 has a coarser resolution than NEI (but still fits in the NEI coverage?) --> mass conservative interpolation
                #       d01 has a coarser resolution than NEI and also extends beyond the NEI domain --> Use MACC (how about MACC resolution?)
                msg('Interpolating...',type=1)
                for domain in domains:
                    nlat,nlon    = xlat[domain].shape
                    if (tracer == selectedtracers[0]):
                        msg('Converting (lon,lat) in (x,y) for %s'%domain,type=1) 
                        Xwrf[domain] = np.empty((nlat,nlon))
                        Ywrf[domain] = np.empty((nlat,nlon))
                        for ilat,ilon in product(range(nlat), range(nlon)):
                            x,y   = conv.fromWgsToRd([xlat[domain][ilat,ilon],xlon[domain][ilat,ilon]])
            #               msg xlon[domain][ilat,ilon], xlat[domain][ilat,ilon], x, y
                            Xwrf[domain][ilat,ilon] = x
                            Ywrf[domain][ilat,ilon] = y
                        
                    msg('Interpolating %s...'%tracer,type=1)
                    E_wrf[domain,tracer] = np.zeros((Nsnaps,nlat,nlon))
                    for isnap in range(Nsnaps):
                        msg('snap %02d...'%isnap,type=1) 
                        E_wrf[domain,tracer][isnap,:,:] = griddata((Xnei.flatten(),Ynei.flatten()),E_nei[tracer][isnap].flatten() , (Xwrf[domain],Ywrf[domain]),method='nearest')
                    # of
                    # sst_smooth = basemap.interp(sst_coarse, lons_sub[0,:], lats_sub[:,0], *np.meshgrid(lons, lats), order=1)
                    #E_NOx_wrf [domain] = mapinterp(E_NOx_nei [Isnaps].sum(axis=0),xnei,ynei, Xwrf[domain],Ywrf[domain],order=0) # order 0:Nearest
                    #E_PM10_wrf[domain] = mapinterp(E_PM10_nei[Isnaps].sum(axis=0),xnei,ynei, Xwrf[domain],Ywrf[domain],order=0) # order 1: bilinear

                E_wrf_tmp = {}
                for domain in domains:
                    E_wrf_tmp[domain,tracer] = E_wrf[domain,tracer]
                    pickledata     = (Xwrf,Ywrf,E_wrf_tmp,domains,tracers)
                    cPickle.dump(pickledata,open(picklefilename,'wb'))
                del(E_wrf_tmp)
            
        progress.append('NEIinterpolatedToWRF')
        msg('done.')

    moy = moys[month-1]
    if not 'wrfchemi_written' in progress:
        for domain,zTime,dow in product(domains,zTimes,dows):
            print('Writing wrfchemifile for %s %s %s %s'%(domain,zTime,moy,dow))
            #its = range(12) if zTime == '00z' else range(12,24)
            zTimeOffset = 12 if zTime == '12z' else 0
        
            
            wrfoutfilename            = '%s/wrfinput_%s'%(wrfoutpath,domain)
            wrfoutfile                = nc.Dataset(wrfoutfilename,'r')
            wrfchemifilename          = '%s/NEI_wrfchemi_%s_%s_%s_%s'%(wrfchemipath,zTime,domain,moy,dow)
            wrfchemifile              = nc.Dataset(wrfchemifilename, 'w',format = 'NETCDF3_CLASSIC' )
        
            nz                        =  1
            nlat,nlon                 = xlat[domain].shape

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
            ncXLAT[:]                 = xlat[domain]
            ncXLONG.FieldType         = 104.
            ncXLONG.MemoryOrder       = 'XY'
            ncXLONG.description       = 'LONGGITUDE, WEST IS NEGATIVE'
            ncXLONG.units             = 'degree_north'
            ncXLONG.stagger           = ' ' 
            ncXLONG.coordinates       = 'XLONG XLAT'
            ncXLONG[:]                = xlon[domain]

            ncE_NOx_stat.FieldType    = 104.
            ncE_NOx_stat.MemoryOrder  = 'XYZ'
            ncE_NOx_stat.description  = 'NOx emissions from stationary sources'
            ncE_NOx_stat.units        = 'mole/km2/hr'
            ncE_NOx_stat.stagger      = ' ' 
            ncE_NOx_stat.coordinates  = 'XLONG XLAT'   
            Isnap                     = range(13)
            for it in range(12):     
                TP                    = TP_hod['%02d'%(it+zTimeOffset)][Isnap] * TP_dow[dow][Isnap] * TP_moy[moy][Isnap]
                TP2                   = np.tile(TP[...,np.newaxis,np.newaxis],(1,nlat,nlon))
                ncE_NOx_stat[it,0,:,:]= (TP2 * E_wrf[domain,'NOx'][Isnap,:,:]).sum(axis=0)
    #           ncE_NOx_stat[it,:]    = (TP2 * E_wrf[domain,'NOx'][Isnap]).sum(axis=0)     
        
            ncE_NOx_traf.FieldType    = 104.
            ncE_NOx_traf.MemoryOrder  = 'XYZ'
            ncE_NOx_traf.description  = 'NOx emissions from road traffic'
            ncE_NOx_traf.units        = 'mole/km2/hr'
            ncE_NOx_traf.stagger      = ' ' 
            ncE_NOx_traf.coordinates  = 'XLONG XLAT'   
            Isnap                     = 13
            for it in range(12):     
                TP                    = TP_hod['%02d'%(it+zTimeOffset)][Isnap] * TP_dow[dow][Isnap] * TP_moy[moy][Isnap]
                ncE_NOx_traf[it,0,:,:] = TP*E_wrf[domain,'NOx'][Isnap,:,:]
    #           TP2                   = np.tile(TP[...,np.newaxis,np.newaxis],(1,nlat,nlon))            

            ncE_PM10_stat.FieldType   = 104.
            ncE_PM10_stat.MemoryOrder = 'XYZ'
            ncE_PM10_stat.description = 'PM10 emissions from stationary sources'
            ncE_PM10_stat.units       = 'kg/km2/hr'
            ncE_PM10_stat.stagger     = ' ' 
            ncE_PM10_stat.coordinates = 'XLONG XLAT'  
            Isnap                     = range(13)
            for it in range(12):     
                TP                    = TP_hod['%02d'%(it+zTimeOffset)][Isnap] * TP_dow[dow][Isnap] * TP_moy[moy][Isnap]
                TP2                   = np.tile(TP[...,np.newaxis,np.newaxis],(1,nlat,nlon))
                ncE_PM10_stat[it,0,:,:]= (TP2 * E_wrf[domain,'PM10'][Isnap,:,:]).sum(axis=0)

            ncE_PM10_traf.FieldType   = 104.
            ncE_PM10_traf.MemoryOrder = 'XYZ'
            ncE_PM10_traf.description = 'PM10 emissions from road traffic'
            ncE_PM10_traf.units       = 'kg/km2/hr'
            ncE_PM10_traf.stagger     = ' ' 
            ncE_PM10_traf.coordinates = 'XLONG XLAT'  
            Isnap                     = 13
            for it in range(12):     
                TP                    = TP_hod['%02d'%(it+zTimeOffset)][Isnap] * TP_dow[dow][Isnap] * TP_moy[moy][Isnap]
    #           TP2                   = np.tile(TP[...,np.newaxis,np.newaxis],(1,nlat,nlon))
                ncE_PM10_traf[it,0,:,:]   = TP*E_wrf[domain,'PM10'][Isnap,:,:]
                
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

    #--- output
    if not 'dofig' in globals(): dofig = [None, False, True, False]
    dofig = 4*(False,)
    #E_CO2_nei = np.ma.masked_where(E_CO2_nei == 0, E_CO2_nei)
    #E_NOx_nei = np.ma.masked_where(E_NOx_nei == 0, E_NOx_nei)
    if dofig[1]:
        tracer = 'NOx'
        f = plt.figure(1)   
        f.clf()
        ax = f.add_subplot(111)
        ax.scatter(data[tracer][:,0],data[tracer][:,1],s=20,c=np.log10(data[tracer][:,5]),marker='.',edgecolors=None)   
        msg('fig 1 done.')
    plt.show()
        
    if dofig[2]:
        Isnap = range(14)
        tracer = 'NOx'
        f = plt.figure(2)
        f.clf()
        ax = f.add_subplot(111)
        Z = E_nei[tracer][Isnap,:].sum(axis=0); Z = np.ma.masked_where(Z == 0, Z)
        h=ax.pcolor(Xnei,Ynei,np.log10(Z))
        for domain in domains[:3]:
            xwrf = Xwrf[domain]
            ywrf = Ywrf[domain]
            ax.plot([xwrf.min(), xwrf.max(), xwrf.max(),xwrf.min(),xwrf.min()],[ywrf.min(),ywrf.min(),ywrf.max(),ywrf.max(),ywrf.min()],'k-')
            cb = f.colorbar(h,orientation='vertical')   
            pow = list()
            for YTL in cb.ax.get_ymajorticklabels(): pow.append(np.double(YTL.get_text().replace(u'\N{MINUS SIGN}', '-')))
            cb.ax.set_yticklabels(np.power(10,pow),fontsize=FS)

        msg('fig 2 done.')
        dofig[2] = False

    if dofig[3]:
        msg('Doing figure 3...',type=1)
        f = plt.figure(3)
        f.clf()
        for domain in domains[::-1]:
            tracer = 'NOx'
        
            ax= f.add_subplot(111)
            Z = E_wrf_stat[domain,tracer]+E_wrf_traf[domain,tracer]; Z = np.ma.masked_where(Z == 0, Z)
            h=ax.pcolor(Xwrf[domain],Ywrf[domain],np.log10(Z),vmin=0,vmax=np.log10(2000))
            xwrf = Xwrf[domain]
            ywrf = Ywrf[domain]
            ax.plot([xwrf.min(), xwrf.max(), xwrf.max(),xwrf.min(),xwrf.min()],[ywrf.min(),ywrf.min(),ywrf.max(),ywrf.max(),ywrf.min()],'k-')

        cb = f.colorbar(h,orientation='vertical')   
        pow = list()
        for YTL in cb.ax.get_ymajorticklabels(): pow.append(np.double(YTL.get_text().replace(u'\N{MINUS SIGN}', '-')))
        cb.ax.set_yticklabels(np.power(10,pow),fontsize=FS)
        
        msg('Done figure 3.')


if __name__=="__main__":
    if platform == 'win32':
        NEIpath      = 'U:/AMS_Stimulus_2016/Data_share/Emissions/NEI'
        wrfchemipath = 'U:/AMS_Stimulus_2016/Data_share/Emissions/wrfchemi'
        wrfoutpath   = 'S:/data/HURP/Adam_wrfout'
    elif (platform == 'darwin' and exists('/Users/molen050/mnt/promise')): # CapeGrim
        NEIpath      = '/Users/molen050/mnt/promise/WRF/michiel/HURP/Data_share/Emissions/NEI'
        wrfchemipath = '/Users/molen050/mnt/promise/WRF/michiel/HURP/Data_share/Emissions/wrfchemi'
        wrfoutpath   = '/Users/molen050/mnt/promise/WRF/michiel/HURP/Data_share/Domains'
    elif (platform == 'linux2' ): # Cartesius
        print('Setting env for Cartesius.')
        NEIpath      = '/projects/0/aams/wrfv3/Data_share/Emissions/NEI'
        wrfchemipath = '/projects/0/aams/wrfv3/Data_share/Emissions/wrfchemi'
        wrfoutpath   = '/projects/0/aams/wrfv3/OutputMichiel/HURP_20150630-20150704'    
    else:  # maunaloa
        NEIpath      = '/Storage/WRF/michiel/HURP/Data_share/Emissions/NEI'
        wrfchemipath = '/Storage/WRF/michiel/HURP/Data_share/Emissions/wrfchemi'
        wrfoutpath   = '/Storage/WRF/michiel/HURP/Data_share/Domains'
        
    Times            = dict()
    Times['00z']     = [['2014-01-01_00:00:00'],['2014-01-01_01:00:00'],['2014-01-01_02:00:00'],['2014-01-01_03:00:00'],
                        ['2014-01-01_04:00:00'],['2014-01-01_05:00:00'],['2014-01-01_06:00:00'],['2014-01-01_07:00:00'],
                        ['2014-01-01_08:00:00'],['2014-01-01_09:00:00'],['2014-01-01_10:00:00'],['2014-01-01_11:00:00']]
                        
    Times['12z']     = [['2014-01-01_12:00:00'],['2014-01-01_13:00:00'],['2014-01-01_14:00:00'],['2014-01-01_15:00:00'],
                        ['2014-01-01_16:00:00'],['2014-01-01_17:00:00'],['2014-01-01_18:00:00'],['2014-01-01_19:00:00'],
                        ['2014-01-01_20:00:00'],['2014-01-01_21:00:00'],['2014-01-01_22:00:00'],['2014-01-01_23:00:00']]
    nei2wrf(NEIpath, wrfchemipath, wrfoutpath, Times)
