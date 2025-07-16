#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

-----------------------------------------------------------------
"""

#------------------------------------------------------------------------------
#                   NOTES 
# RTTOV IS RUN OUTSIDE. 
# CURRENTLY THE CODE BELOW IS HARDCODED FOR AMSRE? BUT SHOULD BE FOR MHS AT 20.00
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# to ensure that pyrttov is importable

import gc
import numpy as np
import matplotlib.pyplot as plt
import Tools2Plot as T2P
from Tools2RunRTTOV import read_wrf, run_IFS_rttov14, filter_pixels, find_pixels, filter_pixels_monotonic, run_IFS_rttov14_version2
from config import config_folders
#import Plots4Analysis as P4A
from netCDF4 import Dataset
import wrf
import xarray as xr
import sys
import seaborn as sns 
import matplotlib as mpl


plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  

#------------------
def check_main(instrument, HHtime, mp_version, server):
        
    plotpath, folders = config_folders(server)

    # Select server and folder locations
    #--------------------------------------------------------------------------
    if 'yakaira' in server:
        paths_rttov  = '/home/vito.galligani/Work/RTTOV/rttov14.0_beta'
        upfolder     = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/'
        general_path = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RTTOVout/'
        prov         = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')
        mhs_noaa19_dir  = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/PMW/GPM/RS/V07/PMW/1C-MHS-NOAA19/2018/11/10/'

    elif 'cnrm' in server:
        upfolder    = '/home/galliganiv/'
        sys.path.insert(1,'/home/galliganiv/ACMS_hail/src')
        processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'+mp_physics

    #--------------------------------------------------------------------------
    # server dependent files    
    from package_functions import pressure2height
    import package_functions as funs
    
    if mp_version == 6:
        mp_physics = 'WRF-WSM6'
        ncfolder  = upfolder+'WRFOUT/WSM6_domain3_NoahMP/'
        ncfile    = ncfolder+'wrfout_d02_2018-11-10_'+HHtime+':00'

    elif mp_version == 6.1:
        mp_physics = 'WRFWSM6_YSU12hr'
        ncfolder  = upfolder+'WRFOUT_1domain/WSM6_1domaintest_12hrs_YSU/'    
        ncfile    = ncfolder+'wrfout_d01_2018-11-10_'+HHtime+':00'  
        
    
    # Select instrument configurations
    #--------------------------------------------------------------------------
    if 'MHS' in instrument:
        nchan = 5
    elif 'AMSR' in instrument:
        nchan = 14


    A = read_wrf(ncfile)
    return A

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
timetime = '20:00'
A = check_main('MHS', timetime , 6.1, 'yakaira')

server = 'yakaira'

if 'yakaira' in server: 
    upfolder     = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/'
    sys.path.insert(1,'/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/src')
    processedFolder = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RTTOVout/Processed/'+'WRF-WSM6'
    mhs_noaa19_dir  = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/PMW/GPM/RS/V07/PMW/1C-MHS-NOAA19/2018/11/10/'
    plotpath, folders = config_folders(server)    
    prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')
    d_cs    = xr.open_dataset(processedFolder+'/'+'output_tb_'+'MHS'+'rttov_processed_clearsky.nc')
    cmaps = T2P.GMI_colormap()     


ii = 230
ij = 550
ji = 100
jj = 600

fig, axes = plt.subplots(nrows=1, ncols=2, constrained_layout=True,figsize=[16,8])
axes[0].pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['MHs_domain_obs'].data[0,:,:], cmap=cmaps['turbo_r'], shading='auto')
axes[1].pcolormesh(A['XLONG'][ii:ij,ji:jj], A['XLAT'][ii:ij,ji:jj], np.sum(A['swc'],0)[ii:ij,ji:jj], cmap='viridis', shading='auto')
axes[0].plot(prov[:,0],prov[:,1],color='r');
axes[1].plot(prov[:,0],prov[:,1],color='r');
for index in range(2):
    axes[index].set_xlim([-66,-58])
    axes[index].set_ylim([-36,-31])   
plt.suptitle('Time:  '+ timetime)    
fig.savefig(plotpath+'/domain1/MHS_realobservations_andselecteddomain4rttov.png', dpi=300, transparent=False, bbox_inches='tight')   

ii = 390
ij = 550
ji = 100
jj = 419

fig, axes = plt.subplots(nrows=1, ncols=2, constrained_layout=True,figsize=[16,8])
axes[0].pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['MHs_domain_obs'].data[0,:,:], cmap=cmaps['turbo_r'], shading='auto')
axes[1].pcolormesh(A['XLONG'][ii:ij,ji:jj], A['XLAT'][ii:ij,ji:jj], np.sum(A['swc'],0)[ii:ij,ji:jj], cmap='viridis', shading='auto')
axes[0].plot(prov[:,0],prov[:,1],color='r');
axes[1].plot(prov[:,0],prov[:,1],color='r');
for index in range(2):
    axes[index].set_xlim([-66,-58])
    axes[index].set_ylim([-36,-31])   
plt.suptitle('Time:  '+ timetime)    
fig.savefig(plotpath+'/domain1/MHS_realobservations_andselecteddomain4rttov_zoom.png', dpi=300, transparent=False, bbox_inches='tight')   






