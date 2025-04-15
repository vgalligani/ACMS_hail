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
from Tools2RunRTTOV import read_wrf, run_IFS_rttov14_singleProf, SelectProf
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
from package_functions import pressure2height


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
def reorder_bssp_mhs(var):
    
    nchan = 5
    nlev  = 45
    
    var1    = np.zeros( (nchan, 5, nlev-1) ); var1[:]=np.nan  
    counter = 0
    for ichan in range(nchan):
        for ilay in range(nlev-1): 
            for itype in range(5):
                var1[ichan,itype,ilay] = var[counter]
                counter = counter+1
                
    return var1
    

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def get_SingleProf_rttov(mp_version, server, HHtime, instrument, i, j, rsgExp): 

    nlev = 45
    
    if 'MHS' in instrument: 
        nchan = 5
        
    if mp_version == 6:
        mp_physics = 'WRF-WSM6'

    plotpath, folders = config_folders(server)

    outfolder   = folders['read_out_dir']+mp_physics+'/'+mp_physics+'_20181110_'+HHtime+'_'+instrument
    outfolder   = outfolder+'_atlas_satzen__i'+str(i)+'_j'+str(j)+'__eqmassWSM6_rsg_'+rsgExp+'/'
       
    outfile_tb  = np.genfromtxt(outfolder+'output_as_tb_'+instrument+'_eqmassWSM6_rsg_'+rsgExp)
    outfile_asy = np.genfromtxt(outfolder+'optp_hydro_asy.txt')
    outfile_ext = np.genfromtxt(outfolder+'optp_hydro_ext.txt')
    outfile_ssa = np.genfromtxt(outfolder+'optp_hydro_ssa.txt')
    
    outfile_asy   = reorder_bssp_mhs(outfile_asy)
    outfile_ext   = reorder_bssp_mhs(outfile_ext)
    outfile_ssa   = reorder_bssp_mhs(outfile_ssa)

    scat = {'asy':outfile_asy,'ext':outfile_ext,'ssa':outfile_ssa}


    return outfile_tb, scat

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
        
    plotpath, folders = config_folders(server)

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
    fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[8*3,8])
    pcm = axes[0].pcolormesh(lon, lat, tb[1,:,:],  cmap=cmaps['turbo_r'], vmin=200, vmax=300)
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
    axes[0].set_title('Simulated BT 157GHz')
    
    
    # Agrego integrated swc
    gcm = axes[1].pcolormesh(lon, lat, dWRF['WRF_intqs'],  cmap=cmaps['turbo_r'], vmin=0, vmax=20)
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
    gcm = axes[2].pcolormesh(lon, lat, dWRF['WRF_intqg'],  cmap=cmaps['turbo_r'], vmin=0, vmax=40)
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
    
    # Agrego integrated rwc
    gcm = axes[3].pcolormesh(lon, lat, dWRF['WRF_intqr'],  cmap=cmaps['turbo_r'], vmin=0, vmax=20)
    cbar = plt.colorbar(gcm, ax=axes[3], shrink=1)
    axes[3].plot(lon[index_i, index_j], lat[index_i, index_j], marker='s', color='m', markersize=5)
    axes[3].plot(prov[:,0],prov[:,1],color='w');
    if zoom == 0: 
        axes[3].set_xlim([-67,-62])   
        axes[3].set_ylim([-34.5,-30])   
    else:
        axes[3].set_xlim([-64.5,-62])   
        axes[3].set_ylim([-34,-32])           
        
    #Agrego el cuadrado con el que grafica matplotlib? 
    axes[3].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
    axes[3].set_xlabel('Longitude')
    axes[3].set_ylabel('Latitude')
    axes[3].set_title('integrated rain content')    
    
    # Agrego perfiles de qx
    axes[4].plot(dWRF['WRF_qr'][:,index_i, index_j], dWRF['wrf_zalt'][:,index_i, index_j]/1e3, linestyle='-', color='k', linewidth=1.2, label='rain' )
    axes[4].plot(dWRF['WRF_qs'][:,index_i, index_j], dWRF['wrf_zalt'][:,index_i, index_j]/1e3, linestyle='-', color='darkblue',linewidth=1.2, label='snow' )
    axes[4].plot(dWRF['WRF_qg'][:,index_i, index_j], dWRF['wrf_zalt'][:,index_i, index_j]/1e3, linestyle='-', color='darkred',linewidth=1.2, label='grau' )
    axes[4].plot(dWRF['WRF_qi'][:,index_i, index_j], dWRF['wrf_zalt'][:,index_i, index_j]/1e3, linestyle='-', color='darkgreen',linewidth=1.2, label='ice' )
    axes[4].set_ylim([0, 20])
    axes[4].set_xlabel('WRF qx')
    axes[4].set_ylabel('Altitude [km]')
    axes[4].set_title('qx Profile')
    axes[4].legend()
    
    fig.savefig(plotpath+'/RTTOV/'+'Map4_i_j'+str(index_i)+'_'+str(index_j)+'.png', dpi=300,transparent=False)   
    
    
    return

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def bulkPlot(server, experiment, d_asExp1, dWRF, scat, zoom, index_i, index_j): 
    
    cmaps = T2P.GMI_colormap() 
    cmap = plt.cm.viridis  # Colormap
    cmap.set_bad(color='gray')  # Color for NaN values

    # Some basic info for colormaps  
    if 'yakaira' in server:
        prov     = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')    
        fn       = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
        upfolder = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/'
        
    elif 'cnrm' in server:
        prov     = np.genfromtxt("/home/galliganiv/ACMS_hail/src/provincias.txt", delimiter='')    
        fn       = '/home/galliganiv/ACMS_hail/src/etopo1_bedrock.nc'   
        upfolder = '/home/galliganiv/'       
        
    plotpath, folders = config_folders(server)

    # Identify profiles of interest
    tb  = d_asExp1['rttov_as']
    lon = d_asExp1['wrf_lon']
    lat = d_asExp1['wrf_lat']
    
    # Get profile
    
    ncfolder   = upfolder+'WRFOUT/WSM6_domain3_NoahMP/'
    ncfile   = ncfolder+'wrfout_d02_2018-11-10_'+'20:30'+':00'
    A        = read_wrf(ncfile)
    Ai       = SelectProf(A,index_i,index_j)
    
    # Topografia
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
    
    # Plot figure of an example of the fooprint search 
    fig, axes = plt.subplots(nrows=3, ncols=5, constrained_layout=True,figsize=[8*3,8*3])
    pcm = axes[0,0].pcolormesh(lon, lat, tb[1,:,:],  cmap=cmaps['turbo_r'], vmin=200, vmax=300)
    cbar = plt.colorbar(pcm, ax=axes[0,0], shrink=1)
    axes[0,0].plot(lon[index_i, index_j], lat[index_i, index_j],  marker='s', color='m', markersize=5)
    axes[0,0].plot(prov[:,0],prov[:,1],color='w'); 
    if zoom == 0: 
        axes[0,0].set_xlim([-67,-62])   
        axes[0,0].set_ylim([-34.5,-30])   
    else:
        axes[0,0].set_xlim([-64.5,-62])   
        axes[0,0].set_ylim([-34,-32])     
    #Agrego el cuadrado con el que grafica matplotlib? 
    axes[0,0].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
    axes[0,0].set_xlabel('Longitude')
    axes[0,0].set_ylabel('Latitude')
    axes[0,0].set_title('Simulated BT 157GHz')
        
    
    # Agrego perfiles de qx
    #breakpoint()Ai['pressure'][::-1]
    wrfzalt = pressure2height(Ai['pressure'][:], Ai['T'][:])/1e3
    
    axes[0,1].plot(dWRF['WRF_qr'][:,index_i, index_j], wrfzalt, linestyle='-', color='k', linewidth=1.2, label='rain' )
    axes[0,1].plot(dWRF['WRF_qs'][:,index_i, index_j], wrfzalt, linestyle='-', color='darkblue',linewidth=1.2, label='snow' )
    axes[0,1].plot(dWRF['WRF_qg'][:,index_i, index_j], wrfzalt, linestyle='-', color='darkred',linewidth=1.2, label='grau' )
    axes[0,1].plot(dWRF['WRF_qi'][:,index_i, index_j], wrfzalt, linestyle='-', color='darkgreen',linewidth=1.2, label='ice' )
    axes[0,1].set_ylim([0, 20])
    axes[0,1].set_xlabel('WRF qx')
    axes[0,1].set_ylabel('Altitude [km]')
    axes[0,1].set_title('qx Profile')
    axes[0,1].legend()
    
    
    # Agrego perfiles de qx
    counter = 0
    ichans  = [0,1,4]
    ichans_title = ['89GHz', '157GHz', '190GHz']
    for ichan in ichans:
        axes[counter,2].plot(scat['asy'][ichan,0,::-1], wrfzalt , linestyle='-', color='k', linewidth=1.2, label='rain' )
        axes[counter,2].plot(scat['asy'][ichan,1,::-1], wrfzalt , linestyle='-', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,2].plot(scat['asy'][ichan,2,::-1], wrfzalt , linestyle='-', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,2].plot(scat['asy'][ichan,4,::-1], wrfzalt , linestyle='-', color='darkgreen',linewidth=1.2, label='ice' )
        axes[counter,2].set_ylim([0, 20])
        axes[counter,2].set_xlabel('rttov asy')
        axes[counter,2].set_ylabel('Altitude [km]')
        axes[counter,2].set_title('Bulk asy Profile (' + ichans_title[counter] +')')
        axes[counter,2].legend()    
    
        # Agrego perfiles de qx
        axes[counter,3].plot(scat['ext'][ichan,0,::-1], wrfzalt , linestyle='-', color='k', linewidth=1.2, label='rain' )
        axes[counter,3].plot(scat['ext'][ichan,1,::-1], wrfzalt , linestyle='-', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,3].plot(scat['ext'][ichan,2,::-1], wrfzalt , linestyle='-', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,3].plot(scat['ext'][ichan,4,::-1], wrfzalt , linestyle='-', color='darkgreen',linewidth=1.2, label='ice' )
        axes[counter,3].set_ylim([0, 20])
        axes[counter,3].set_xlabel('rttov ext')
        axes[counter,3].set_ylabel('Altitude [km]')
        axes[counter,3].set_title('Bulk ext Profile (' + ichans_title[counter] +')')
        axes[counter,3].legend()    
    
        # Agrego perfiles de qx
        axes[counter,4].plot(scat['ssa'][ichan,0,::-1], wrfzalt , linestyle='-', color='k', linewidth=1.2, label='rain' )
        axes[counter,4].plot(scat['ssa'][ichan,1,::-1], wrfzalt , linestyle='-', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,4].plot(scat['ssa'][ichan,2,::-1], wrfzalt , linestyle='-', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,4].plot(scat['ssa'][ichan,4,::-1], wrfzalt , linestyle='-', color='darkgreen',linewidth=1.2, label='ice' )
        axes[counter,4].set_ylim([0, 20])
        axes[counter,4].set_xlabel('rttov ssa')
        axes[counter,4].set_ylabel('Altitude [km]')
        axes[counter,4].set_title('Bulk ssa Profile (' + ichans_title[counter] +')')
        axes[counter,4].legend()   
        
        counter=counter+1

    # Remove empty xlabels:
    fig.delaxes(axes[1,0])
    fig.delaxes(axes[2,0])
    fig.delaxes(axes[1,1])
    fig.delaxes(axes[2,1])
    
    fig.text(0.1, 0.62, '(i,j)='+str(index_i)+', '+str(index_j), fontsize=14, 
             ha='center', va='center', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round, pad=0.5') )
    
    fig.savefig(plotpath+'/RTTOV/'+'BulkScattering_equalmass_i_j'+str(index_i)+'_'+str(index_j)+'_'+experiment+'.png', dpi=300,transparent=False)   
        
    
    return

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def bulkPlotcompare(server, experiments, d_asExp1, d_asExp2, d_asExp3, dWRF, scat1, scat2, scat3, zoom, index_i, index_j): 
    
        
    cmaps = T2P.GMI_colormap() 
    cmap = plt.cm.viridis  # Colormap
    cmap.set_bad(color='gray')  # Color for NaN values

    # Some basic info for colormaps  
    if 'yakaira' in server:
        prov     = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')    
        fn       = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
        upfolder = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/'
        
    elif 'cnrm' in server:
        prov     = np.genfromtxt("/home/galliganiv/ACMS_hail/src/provincias.txt", delimiter='')    
        fn       = '/home/galliganiv/ACMS_hail/src/etopo1_bedrock.nc'   
        upfolder = '/home/galliganiv/'       
        
    plotpath, folders = config_folders(server)

    # Identify profiles of interest
    tb1  = d_asExp1['rttov_as']
    tb2  = d_asExp2['rttov_as']
    tb3  = d_asExp3['rttov_as']

    lon = d_asExp1['wrf_lon']
    lat = d_asExp1['wrf_lat']
    
    # Get profile
    
    ncfolder   = upfolder+'WRFOUT/WSM6_domain3_NoahMP/'
    ncfile   = ncfolder+'wrfout_d02_2018-11-10_'+'20:30'+':00'
    A        = read_wrf(ncfile)
    Ai       = SelectProf(A,index_i,index_j)
    
    # Topografia
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
    
    # Plot figure of an example of the fooprint search 
    fig, axes = plt.subplots(nrows=3, ncols=5, constrained_layout=True,figsize=[8*3,8*3])
    pcm = axes[0,0].pcolormesh(lon, lat, tb1[1,:,:],  cmap=cmaps['turbo_r'], vmin=200, vmax=300)
    cbar = plt.colorbar(pcm, ax=axes[0,0], shrink=1)
    axes[0,0].plot(lon[index_i, index_j], lat[index_i, index_j],  marker='s', color='m', markersize=5)
    axes[0,0].plot(prov[:,0],prov[:,1],color='w'); 
    if zoom == 0: 
        axes[0,0].set_xlim([-67,-62])   
        axes[0,0].set_ylim([-34.5,-30])   
    else:
        axes[0,0].set_xlim([-64.5,-62])   
        axes[0,0].set_ylim([-34,-32])     
    #Agrego el cuadrado con el que grafica matplotlib? 
    axes[0,0].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
    axes[0,0].set_xlabel('Longitude')
    axes[0,0].set_ylabel('Latitude')
    axes[0,0].set_title('Simulated BT 157GHz ('+experiments[0]+')')
        
    pcm = axes[1,0].pcolormesh(lon, lat, tb2[1,:,:],  cmap=cmaps['turbo_r'], vmin=200, vmax=300)
    cbar = plt.colorbar(pcm, ax=axes[1,0], shrink=1)
    axes[1,0].plot(lon[index_i, index_j], lat[index_i, index_j],  marker='s', color='m', markersize=5)
    axes[1,0].plot(prov[:,0],prov[:,1],color='w'); 
    if zoom == 0: 
        axes[1,0].set_xlim([-67,-62])   
        axes[1,0].set_ylim([-34.5,-30])   
    else:
        axes[1,0].set_xlim([-64.5,-62])   
        axes[1,0].set_ylim([-34,-32])     
    #Agrego el cuadrado con el que grafica matplotlib? 
    axes[1,0].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
    axes[1,0].set_xlabel('Longitude')
    axes[1,0].set_ylabel('Latitude')
    axes[1,0].set_title('Simulated BT 157GHz ('+experiments[1]+')')
    
    pcm = axes[2,0].pcolormesh(lon, lat, tb3[1,:,:],  cmap=cmaps['turbo_r'], vmin=200, vmax=300)
    cbar = plt.colorbar(pcm, ax=axes[2,0], shrink=1)
    axes[2,0].plot(lon[index_i, index_j], lat[index_i, index_j],  marker='s', color='m', markersize=5)
    axes[2,0].plot(prov[:,0],prov[:,1],color='w'); 
    if zoom == 0: 
        axes[2,0].set_xlim([-67,-62])   
        axes[2,0].set_ylim([-34.5,-30])   
    else:
        axes[2,0].set_xlim([-64.5,-62])   
        axes[2,0].set_ylim([-34,-32])     
    #Agrego el cuadrado con el que grafica matplotlib? 
    axes[2,0].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
    axes[2,0].set_xlabel('Longitude')
    axes[2,0].set_ylabel('Latitude')
    axes[2,0].set_title('Simulated BT 157GHz ('+experiments[2]+')')
    
    # Agrego perfiles de qx
    #breakpoint()Ai['pressure'][::-1]
    wrfzalt = pressure2height(Ai['pressure'][:], Ai['T'][:])/1e3
    
    axes[0,1].plot(dWRF['WRF_qr'][:,index_i, index_j], wrfzalt, linestyle='-', color='k', linewidth=1.2, label='rain' )
    axes[0,1].plot(dWRF['WRF_qs'][:,index_i, index_j], wrfzalt, linestyle='-', color='darkblue',linewidth=1.2, label='snow' )
    axes[0,1].plot(dWRF['WRF_qg'][:,index_i, index_j], wrfzalt, linestyle='-', color='darkred',linewidth=1.2, label='grau' )
    axes[0,1].plot(dWRF['WRF_qi'][:,index_i, index_j], wrfzalt, linestyle='-', color='darkgreen',linewidth=1.2, label='ice' )
    axes[0,1].set_ylim([0, 20])
    axes[0,1].set_xlabel('WRF qx')
    axes[0,1].set_ylabel('Altitude [km]')
    axes[0,1].set_title('qx Profile')
    axes[0,1].legend()
    
    
    # Agrego perfiles de qx
    counter = 0
    ichans  = [0,1,4]
    ichans_title = ['89GHz', '157GHz', '190GHz']
    for ichan in ichans:
        axes[counter,2].plot(scat1['asy'][ichan,0,::-1], wrfzalt , linestyle='-', color='k', linewidth=1.2, label='rain' )
        axes[counter,2].plot(scat1['asy'][ichan,1,::-1], wrfzalt , linestyle='-', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,2].plot(scat1['asy'][ichan,2,::-1], wrfzalt , linestyle='-', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,2].plot(scat1['asy'][ichan,4,::-1], wrfzalt , linestyle='-', color='darkgreen',linewidth=1.2, label='ice' )
        axes[counter,2].legend()    
        #
        axes[counter,2].plot(scat2['asy'][ichan,0,::-1], wrfzalt , linestyle='--', color='k', linewidth=1.2, label='rain' )
        axes[counter,2].plot(scat2['asy'][ichan,1,::-1], wrfzalt , linestyle='--', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,2].plot(scat2['asy'][ichan,2,::-1], wrfzalt , linestyle='--', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,2].plot(scat2['asy'][ichan,4,::-1], wrfzalt , linestyle='--', color='darkgreen',linewidth=1.2, label='ice' )        
        #
        axes[counter,2].plot(scat3['asy'][ichan,0,::-1], wrfzalt , linestyle=':', color='k', linewidth=1.2, label='rain' )
        axes[counter,2].plot(scat3['asy'][ichan,1,::-1], wrfzalt , linestyle=':', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,2].plot(scat3['asy'][ichan,2,::-1], wrfzalt , linestyle=':', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,2].plot(scat3['asy'][ichan,4,::-1], wrfzalt , linestyle=':', color='darkred',linewidth=1.2, label='grau' )
        #
        axes[counter,2].set_ylim([0, 20])
        axes[counter,2].set_xlabel('rttov asy')
        axes[counter,2].set_ylabel('Altitude [km]')
        axes[counter,2].set_title('Bulk asy Profile (' + ichans_title[counter] +')')
    
        # Agrego perfiles de qx
        axes[counter,3].plot(scat1['ext'][ichan,0,::-1], wrfzalt , linestyle='-', color='k', linewidth=1.2, label='rain' )
        axes[counter,3].plot(scat1['ext'][ichan,1,::-1], wrfzalt , linestyle='-', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,3].plot(scat1['ext'][ichan,2,::-1], wrfzalt , linestyle='-', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,3].plot(scat1['ext'][ichan,4,::-1], wrfzalt , linestyle='-', color='darkgreen',linewidth=1.2, label='ice' )
        #
        axes[counter,3].plot(scat2['ext'][ichan,0,::-1], wrfzalt , linestyle='--', color='k', linewidth=1.2, label='rain' )
        axes[counter,3].plot(scat2['ext'][ichan,1,::-1], wrfzalt , linestyle='--', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,3].plot(scat2['ext'][ichan,2,::-1], wrfzalt , linestyle='--', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,3].plot(scat2['ext'][ichan,4,::-1], wrfzalt , linestyle='--', color='darkgreen',linewidth=1.2, label='ice' )        
        #
        axes[counter,3].plot(scat3['ext'][ichan,0,::-1], wrfzalt , linestyle=':', color='k', linewidth=1.2, label='rain' )
        axes[counter,3].plot(scat3['ext'][ichan,1,::-1], wrfzalt , linestyle=':', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,3].plot(scat3['ext'][ichan,2,::-1], wrfzalt , linestyle=':', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,3].plot(scat3['ext'][ichan,4,::-1], wrfzalt , linestyle=':', color='darkgreen',linewidth=1.2, label='ice' )        
        #
        axes[counter,3].set_ylim([0, 20])
        axes[counter,3].set_xlabel('rttov ext')
        axes[counter,3].set_ylabel('Altitude [km]')
        axes[counter,3].set_title('Bulk ext Profile (' + ichans_title[counter] +')')
    
        # Agrego perfiles de qx
        
        axes[counter,4].plot(np.nan, np.nan, linestyle='-', color='gray', linewidth=1.2, label=experiments[0] )
        axes[counter,4].plot(np.nan, np.nan, linestyle='--', color='gray', linewidth=1.2, label=experiments[1] )
        axes[counter,4].plot(np.nan, np.nan, linestyle=':', color='gray', linewidth=1.2, label=experiments[2] )
        axes[counter,4].legend()    

        axes[counter,4].plot(scat1['ssa'][ichan,0,::-1], wrfzalt , linestyle='-', color='k', linewidth=1.2, label='rain' )
        axes[counter,4].plot(scat1['ssa'][ichan,1,::-1], wrfzalt , linestyle='-', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,4].plot(scat1['ssa'][ichan,2,::-1], wrfzalt , linestyle='-', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,4].plot(scat1['ssa'][ichan,4,::-1], wrfzalt , linestyle='-', color='darkgreen',linewidth=1.2, label='ice' )
        #
        axes[counter,4].plot(scat2['ssa'][ichan,0,::-1], wrfzalt , linestyle='--', color='k', linewidth=1.2, label='rain' )
        axes[counter,4].plot(scat2['ssa'][ichan,1,::-1], wrfzalt , linestyle='--', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,4].plot(scat2['ssa'][ichan,2,::-1], wrfzalt , linestyle='--', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,4].plot(scat2['ssa'][ichan,4,::-1], wrfzalt , linestyle='--', color='darkgreen',linewidth=1.2, label='ice' )        
        #
        axes[counter,4].plot(scat3['ssa'][ichan,0,::-1], wrfzalt , linestyle=':', color='k', linewidth=1.2, label='rain' )
        axes[counter,4].plot(scat3['ssa'][ichan,1,::-1], wrfzalt , linestyle=':', color='darkblue',linewidth=1.2, label='snow' )
        axes[counter,4].plot(scat3['ssa'][ichan,2,::-1], wrfzalt , linestyle=':', color='darkred',linewidth=1.2, label='grau' )
        axes[counter,4].plot(scat3['ssa'][ichan,4,::-1], wrfzalt , linestyle=':', color='darkgreen',linewidth=1.2, label='ice' )        
        #
        axes[counter,4].set_ylim([0, 20])
        axes[counter,4].set_xlabel('rttov ssa')
        axes[counter,4].set_ylabel('Altitude [km]')
        axes[counter,4].set_title('Bulk ssa Profile (' + ichans_title[counter] +')')
        
        counter=counter+1

    # Remove empty xlabels:
    #fig.delaxes(axes[1,0])
    #fig.delaxes(axes[2,0])
    fig.delaxes(axes[1,1])
    fig.delaxes(axes[2,1])
    
    fig.text(0.3, 0.62, '(i,j)='+str(index_i)+', '+str(index_j), fontsize=14, 
             ha='center', va='center', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round, pad=0.5') )
    
    fig.savefig(plotpath+'/RTTOV/'+'Compare_BulkScattering_equalmass_i_j'+str(index_i)+'_'+str(index_j)+'_'+experiment+'.png', dpi=300,transparent=False)   
        
    
    return

#-----------------------------------------------------------------------------
def run_main(mp_version, instrument, experiment):
    
    if mp_version == 6:
        mp_physics = 'WRF-WSM6'
            
    outfile         = 'output_tb_'+instrument
    processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'+mp_physics

    d_asExp1        = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_'+experiment+'.nc')
    WRFvars         = Dataset(processedFolder+'/'+'wrfdata_processed.nc')
    
    #- Run for these files rttov on its own and save its bulk scat profiles for
    # rain-only, snow-only, grau-only 
    do_plot = 0 
    if do_plot == 1: 
        identify_index('cnrm', d_asExp1, WRFvars, zoom=1, index_i=120, index_j=280)
        identify_index('cnrm', d_asExp1, WRFvars, zoom=1, index_i=120, index_j=300)
        identify_index('cnrm', d_asExp1, WRFvars, zoom=1, index_i=135, index_j=310)
        identify_index('cnrm', d_asExp1, WRFvars, zoom=1, index_i=148, index_j=320)
    
    # # Make profiles for the i, j single profile and recompile rttov to save bulk properties
    # i=148, j=320 # max int_g 
    # i=120, j=300 # moderate all
    # i=135, j=310 # snow == grau aprox. 
    makeProfs = 0
    if makeProfs == 1:
        ipnd = 1 
        run_IFS_rttov14_singleProf(ipnd, mp_version, '20:30', 'MHS', 'cnrm', 148, 320)
        run_IFS_rttov14_singleProf(ipnd, mp_version, '20:30', 'MHS', 'cnrm', 120, 300)
        run_IFS_rttov14_singleProf(ipnd, mp_version, '20:30', 'MHS', 'cnrm', 135, 310)
    
    # Read scat-outputs and plot
    #-------------------------------------------------------------------------------------------------
    outfile_tb, bulkscat = get_SingleProf_rttov(mp_version, 'cnrm', '20:30', 'MHS', 148, 320, 's14g2') 
    bulkPlot('cnrm', 'rsg_'+experiment, d_asExp1, WRFvars, bulkscat, 1, index_i=148, index_j=320) # for 89GHz 
       
    outfile_tb, bulkscat = get_SingleProf_rttov(mp_version, 'cnrm', '20:30', 'MHS', 120, 300, 's14g2') 
    bulkPlot('cnrm', 'rsg_'+experiment, d_asExp1, WRFvars, bulkscat, 1, index_i=120, index_j=300) # for 89GHz 
           
    outfile_tb, bulkscat = get_SingleProf_rttov(mp_version, 'cnrm', '20:30', 'MHS', 135, 310, 's14g2') 
    bulkPlot('cnrm', 'rsg_'+experiment, d_asExp1, WRFvars, bulkscat, 1, index_i=135, index_j=310) # for 89GHz  
    
    
    return

#-----------------------------------------------------------------------------
def run_main_compare(mp_version, instrument, experiments):
    
    if mp_version == 6:
        mp_physics = 'WRF-WSM6'
            
    outfile         = 'output_tb_'+instrument
    processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'+mp_physics

    d_asExp1        = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_'+experiments[0]+'.nc')
    d_asExp2        = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_'+experiments[1]+'.nc')
    d_asExp3        = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_'+experiments[2]+'.nc')

    WRFvars         = Dataset(processedFolder+'/'+'wrfdata_processed.nc')
           
    # Read scat-outputs and plot
    #-------------------------------------------------------------------------------------------------
    outfile_tb1, bulkscat1 = get_SingleProf_rttov(mp_version, 'cnrm', '20:30', 'MHS', 148, 320, experiments[0]) 
    outfile_tb2, bulkscat2 = get_SingleProf_rttov(mp_version, 'cnrm', '20:30', 'MHS', 148, 320, experiments[1]) 
    outfile_tb3, bulkscat3 = get_SingleProf_rttov(mp_version, 'cnrm', '20:30', 'MHS', 148, 320, experiments[2]) 

    bulkPlotcompare('cnrm', experiments, d_asExp1, d_asExp2, d_asExp3, 
                    WRFvars, bulkscat1, bulkscat2, bulkscat3, 1, index_i=148, index_j=320) 
      
    return

#------------------------------------------------------------------------------
# Bulk scattering properties: for equal mass experiments 
#------------------------------------------------------------------------------
mp_version = 6
instrument = 'MHS'
run_main(mp_version, instrument, 's9g2')
run_main(mp_version, instrument, 's10g2')
run_main(mp_version, instrument, 's14g2')
run_main_compare(mp_version, instrument, ['s9g2','s10g2','s14g2'])

 