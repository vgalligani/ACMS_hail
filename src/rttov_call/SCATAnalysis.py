#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 13:36:42 2025

@author: galliganiv
"""

import gc
import numpy as np
import matplotlib.pyplot as plt
import Tools2Plot as T2P
from Tools2RunRTTOV import read_wrf, run_IFS_rttov14, filter_pixels, find_pixels, filter_pixels_monotonic
from config import config_folders
#import Plots4Analysis as P4A
from netCDF4 import Dataset
import wrf
import xarray as xr
import sys
import seaborn as sns 
import matplotlib as mpl
import pandas as pd
import netCDF4 as nc


sys.path.insert(1,'/home/galliganiv/ACMS_hail/src')

plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def SingleCol(fdat):
    """
    -------------------------------------------------------------
    Loads a single column ascii file. Note here that these files
    are stored as follows:

    Do ichan = 1, nchannels
      Do ilayer = 1, nlevels
      Enddo
    Enddo
    -------------------------------------------------------------
    OUT   dset    Dataset
    IN    fdat    Filename
    -------------------------------------------------------------
    """
    with open(fdat,"r") as f :
        lines   = f.readlines()
        dset    = pd.DataFrame(l.strip().split() for l in lines)
    dset        = np.array(dset)
    dset        = dset[:,0].astype(float)
    f.close()
    return dset


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def identify_index(server, d_asExp1, dWRF, zoom, index_i, index_j): 
    
    cmaps = T2P.GMI_colormap() 
    cmap = plt.cm.viridis  # Colormap
    cmap.set_bad(color='gray')  # Color for NaN values

    # Some basic info for colormaps  
    if 'yakaira' in server:
        prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')    
        fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'

    elif 'cnrm' in server:
        prov = np.genfromtxt("/home/galliganiv/ACMS_hail/src/provincias.txt", delimiter='')    
        fn = '/home/galliganiv/ACMS_hail/src/etopo1_bedrock.nc'   
        
    # Identify profiles of interest
    tb  = d_asExp1['rttov_as']
    lon = d_asExp1['wrf_lon']
    lat = d_asExp1['wrf_lat']
    
    # Topografia
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
    
    # Plot figure of an example of the fooprint search 
    fig, axes = plt.subplots(nrows=1, ncols=3, constrained_layout=True,figsize=[8*3,8])
    pcm = axes[0].pcolormesh(lon, lat, tb[0,:,:],  cmap=cmaps['turbo_r'], vmin=200, vmax=300)
    cbar = plt.colorbar(pcm, ax=axes[0], shrink=1)
    axes[0].plot(lon[index_i, index_j], lat[index_i, index_j],  marker='s', color='m', markersize=5)
    axes[0].plot(prov[:,0],prov[:,1],color='w'); 
    if zoom == 0: 
        axes[0].set_xlim([-67,-62])   
        axes[0].set_ylim([-34.5,-30])   
    else:
        axes[0].set_xlim([-64.5,-62])   
        axes[0].set_ylim([-34,-32])     
    
    
    #Agrego el cuadrado con el que grafica matplotlib? 
    axes[0].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
    axes[0].set_xlabel('Longitude')
    axes[0].set_ylabel('Latitude')
    axes[0].set_title('Simulated BT 89GHz')
    
    
    # Agrego integrated swc
    gcm = axes[1].pcolormesh(lon, lat, dWRF['WRF_intqs'],  cmap=cmaps['turbo_r'], vmin=0, vmax=10)
    cbar = plt.colorbar(gcm, ax=axes[1], shrink=1)
    axes[1].plot(lon[index_i, index_j], lat[index_i, index_j],  marker='s', color='m', markersize=5)
    axes[1].plot(prov[:,0],prov[:,1],color='w');
    if zoom == 0: 
        axes[1].set_xlim([-67,-62])   
        axes[1].set_ylim([-34.5,-30])   
    else:
        axes[1].set_xlim([-64.5,-62])   
        axes[1].set_ylim([-34,-32])           
        
    #Agrego el cuadrado con el que grafica matplotlib? 
    axes[1].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
    axes[1].set_xlabel('Longitude')
    axes[1].set_ylabel('Latitude')
    axes[1].set_title('integrated snow content')
    
    # Agrego integrated gwc
    gcm = axes[2].pcolormesh(lon, lat, dWRF['WRF_intqg'],  cmap=cmaps['turbo_r'], vmin=0, vmax=10)
    cbar = plt.colorbar(gcm, ax=axes[2], shrink=1)
    axes[2].plot(lon[index_i, index_j], lat[index_i, index_j], marker='s', color='m', markersize=5)
    axes[2].plot(prov[:,0],prov[:,1],color='w');
    if zoom == 0: 
        axes[2].set_xlim([-67,-62])   
        axes[2].set_ylim([-34.5,-30])   
    else:
        axes[2].set_xlim([-64.5,-62])   
        axes[2].set_ylim([-34,-32])           
        
    #Agrego el cuadrado con el que grafica matplotlib? 
    axes[2].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
    axes[2].set_xlabel('Longitude')
    axes[2].set_ylabel('Latitude')
    axes[2].set_title('integrated grau content')
    
    
    return

#-----------------------------------------------------------------------------
mp_version = 6
instrument = 'MHS'

if mp_version == 6:
    mp_physics = 'WRF-WSM6'
        
outfile         = 'output_tb_'+instrument
processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'+mp_physics
d_asExp1        = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_exp1.nc')
WRFvars         = Dataset(processedFolder+'/'+'wrfdata_processed.nc')

#- Run for these files rttov on its own and save its bulk scat profiles for
# rain-only, snow-only, grau-only 
identify_index('cnrm', d_asExp1, WRFvars, zoom=0, index_i=120, index_j=280)
identify_index('cnrm', d_asExp1, WRFvars, zoom=0, index_i=120, index_j=300)
identify_index('cnrm', d_asExp1, WRFvars, zoom=0, index_i=135, index_j=310)


# # Get dimensions
# nchan, ix, iy = d_asExp1['rttov_as'].shape

# scatfolders = '/home/vito.galligani/rttov/rttov_test/run_hail/mietables/'
# exp         = 'mhs_eqmassWSM6_rsg_s10g2'


# do i need to add i y j a esto? lee outtb_scatt to get dimensions

# hydro_ssa = np.reshape( SingleCol(scatfolders+exp+'optp_hydro_ssa.txt'),(nchan,nlev_rttov-1,5))
# hydro_asy = np.reshape( SingleCol(scatfolders+exp+'optp_hydro_asy.txt'),(nchan,nlev_rttov-1,5))
# hydro_ext = np.reshape( SingleCol(scatfolders+exp+'optp_hydro_ext.txt'),(nchan,nlev_rttov-1,5))






