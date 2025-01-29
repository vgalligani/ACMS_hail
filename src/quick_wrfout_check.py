#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:24:02 2024

Initial plots to analyse WRF qx fields and radar observations (45 dBZ contours)


help: https://medium.com/@alexgo1/seamless-remote-ipython-kernel-debugging-with-spyder-ide-d34a0a44baaf

@author: vito.galligani
"""

# Importar librerias necesarias 
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import integrate
import wrf
import numpy as np
import warnings
import pyart 
import wradlib as wrl
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import glob
import os 
from PIL import Image
import datetime
import pandas
import matplotlib
import matplotlib.dates as mdates
from matplotlib.cm import get_cmap
import warnings
from matplotlib import ticker

import vito_functions 
from vito_functions import mixr2massconc
from vito_functions import hhmm_to_minutes, minutes_to_hhmm
import vito_plots as VP
from vito_plots  import pyplot_rings
from scipy.optimize import fsolve
import vito_plots


plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  


#------------------------------------------------------------------------------
def treat_RMA1(rfile):
    
    # radar data
    radar       = pyart.io.read(rfile) 
    start_index = radar.sweep_start_ray_index['data'][0] 
    end_index   = radar.sweep_end_ray_index['data'][0]	
    
    # RMA1
    reflectivity_name = 'TH'
    rlats        = radar.gate_latitude['data'][start_index:end_index] 
    rlons        = radar.gate_longitude['data'][start_index:end_index] 
    ZH           = radar.fields[reflectivity_name]['data'][start_index:end_index] 
    RHOHV        = radar.fields['RHOHV']['data'][start_index:end_index] 
    
    # Att. correction
    #iso0 = np.ma.mean(radar.fields['height']['data'][np.where(np.abs(radar.fields['sounding_temperature']['data']) < 0.1)])
    #radar.fields['height_over_iso0'] = deepcopy(radar.fields['height'])
    #radar.fields['height_over_iso0']['data'] -= iso0   
    #spec_at, pia_dict, Z_correc_ZPHI, spec_diff_at, pida_dict, cor_zdr = pyart.correct.calculate_attenuation_zphi(
    #        radar, temp_field='sounding_temperature', iso0_field='height_over_iso0', 
    #        refl_field='TH', phidp_field='PHIDP',
    #        a_coef=0.06, beta=0.8, temp_ref='height_over_iso0')
    #radar.add_field_like(reflectivity_name, 'corrTH', Z_correc_ZPHI['data'], replace_existing=True)   
        
    # Filtrar cosas de RHOHV
    mask = (ZH < 20) & (RHOHV < 0.8)
    ZH[mask] = np.nan
    
    # Filtro a mano para este caso un sector del PPI
    ignore_r = 190
    ignore_g = 200
    
    rlon = rlons[:ignore_r,ignore_g:]
    rlat = rlats[:ignore_r,ignore_g:]
    rZH  = ZH[:ignore_r,ignore_g:]
    
    return rlon, rlat, rZH

#------------------------------------------------------------------------------
def plot_qxs_P3(title, itime, ncfile, rfile, savedir, radar_name): 
    
    if 'RMA1' in radar_name:
        rlons, rlats, ZH = treat_RMA1(rfile)
        print('Note ZH here is uncorrected for attenuation')
        
        # Info RMA1
        radarLAT = -31.441389
        radarLON = -64.191944
        
    elif 'DOW6' in radar_name:
        reflectivity_name = 'DBZHC'
        
        radar       = pyart.io.read(rfile) 
        start_index = radar.sweep_start_ray_index['data'][0] 
        end_index   = radar.sweep_end_ray_index['data'][0]	
        
        rlats        = radar.gate_latitude['data'][start_index:end_index] 
        rlons        = radar.gate_longitude['data'][start_index:end_index] 
        ZH           = radar.fields[reflectivity_name]['data'][start_index:end_index] 

        # Info DOW6
        radarLAT = radar.latitude['data'][0]
        radarLON = radar.longitude['data'][0]

    elif 'DOW7' in radar_name:
        reflectivity_name = 'DBZHC'
        
        radar       = pyart.io.read(rfile) 
        start_index = radar.sweep_start_ray_index['data'][0] 
        end_index   = radar.sweep_end_ray_index['data'][0]	
        
        rlats        = radar.gate_latitude['data'][start_index:end_index] 
        rlons        = radar.gate_longitude['data'][start_index:end_index] 
        ZH           = radar.fields[reflectivity_name]['data'][start_index:end_index] 
        RHOHV        = radar.fields['RHOHV']['data'][start_index:end_index] 
        
        # Filtrar cosas de RHOHV
        mask = (ZH < 20) & (RHOHV < 0.8)
        ZH[mask] = np.nan

        # Info DOW7
        radarLAT = radar.latitude['data'][0]
        radarLON = radar.longitude['data'][0]
        
    elif 'CSAPR2' in radar_name:
        reflectivity_name = 'attenuation_corrected_reflectivity_h'   
        
        radar       = pyart.io.read(rfile) 
        start_index = radar.sweep_start_ray_index['data'][0] 
        end_index   = radar.sweep_end_ray_index['data'][0]	
        
        rlats        = radar.gate_latitude['data'][start_index:end_index] 
        rlons        = radar.gate_longitude['data'][start_index:end_index] 
        ZH           = radar.fields[reflectivity_name]['data'][start_index:end_index] 
        RHOHV        = radar.fields['copol_correlation_coeff']['data'][start_index:end_index] 

        # Filtrar cosas de RHOHV
        mask = (ZH < 20) & (RHOHV < 0.8)
        ZH[mask] = np.nan
        
        # Info CSAPR2
        radarLAT = radar.latitude['data'][0]
        radarLON = radar.longitude['data'][0]
        
    temp     = wrf.g_temp.get_tk(ncfile)
    pressure = wrf.g_pressure.get_pressure(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL
    Re       = 6.3781e6
    z_level  = Re*geopo_p/(Re-geopo_p)
    lat      =  wrf.getvar( ncfile,"lat") 
    lon      =  wrf.getvar( ncfile,"lon")

    qr = mixr2massconc( np.squeeze(ncfile.variables["QRAIN"][itime,:,:,:]  ), pressure, temp )        
    qi = mixr2massconc( np.squeeze(ncfile.variables["QICE"][itime,:,:,:]   ), pressure, temp )        
    qc = mixr2massconc( np.squeeze(ncfile.variables["QCLOUD"][itime,:,:,:] ), pressure, temp )       

    qr_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qr)) , z_level, axis=0)
    qi_int = integrate.trapz(np.ma.array(qi, mask=np.isnan(qi)) , z_level, axis=0)
    qc_int = integrate.trapz(np.ma.array(qc, mask=np.isnan(qc)) , z_level, axis=0)
    
    qr_int[qr_int<0.0001] = np.nan
    qi_int[qi_int<0.0001] = np.nan
    qc_int[qc_int<0.0001] = np.nan

    fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[12,2])
    pcm0 = axes[0].pcolormesh(lon, lat, qr_int, vmin=0,  vmax=12)
    axes[0].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)
    pcm3 = axes[1].pcolormesh(lon, lat, qc_int, vmin=0,  vmax=12)
    axes[1].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)
    pcm4 = axes[2].pcolormesh(lon, lat, qi_int, vmin=0,  vmax=12)
    axes[2].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)
    axes[1].plot(np.nan, np.nan, 'darkred',  linewidth=1.5, label='45dBZ cont.')
    axes[1].legend(fontsize=9, loc='upper left')


    cbar = plt.colorbar(pcm4, ax=axes[2], shrink=1, label='qx [kg/m$^2$]')
    cbar.cmap.set_under('white')
    
    for i in range(3):
        axes[i].grid()
        axes[i].set_xlim([-65.5,-62])
        axes[i].set_ylim([-33.5,-31.3])

    axes[0].set_title('P3 q_rain')                                                  
    axes[1].set_title('P3 q_cld')  
    axes[2].set_title('P3 q_ice(tot)') 

    plt.suptitle(title)
    
    for i in [120, 240]: 
        [lat_radius, lon_radius] = pyplot_rings(radarLAT,radarLON,i)
        for iax in range(3): 
            axes[iax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)  
            
    fig.delaxes(axes[3])
    fig.delaxes(axes[4])
    
    fig.savefig(savedir+title+'.png', dpi=300,transparent=False)
    #plt.close()
    
    #fig, axes = plt.subplots(nrows=1, ncols=2, constrained_layout=True,figsize=[12,2])
    #axes[0].contour(rlons, rlats, ZH, [45], colors='darkred', linewidths=1.5)
    plt.show() 
    
    
    return

#------------------------------------------------------------------------------
def plot_qxs(title, itime, ncfile, rfile, savedir, radar_name): 
    
    if 'RMA1' in radar_name:
        rlons, rlats, ZH = treat_RMA1(rfile)
        print('Note ZH here is uncorrected for attenuation')
        
        # Info RMA1
        radarLAT = -31.441389
        radarLON = -64.191944
        
    elif 'DOW6' in radar_name:
        reflectivity_name = 'DBZHC'
        
        radar       = pyart.io.read(rfile) 
        start_index = radar.sweep_start_ray_index['data'][0] 
        end_index   = radar.sweep_end_ray_index['data'][0]	
        
        rlats        = radar.gate_latitude['data'][start_index:end_index] 
        rlons        = radar.gate_longitude['data'][start_index:end_index] 
        ZH           = radar.fields[reflectivity_name]['data'][start_index:end_index] 

        # Info DOW6
        radarLAT = radar.latitude['data'][0]
        radarLON = radar.longitude['data'][0]

    elif 'DOW7' in radar_name:
        reflectivity_name = 'DBZHC'
        
        radar       = pyart.io.read(rfile) 
        start_index = radar.sweep_start_ray_index['data'][0] 
        end_index   = radar.sweep_end_ray_index['data'][0]	
        
        rlats        = radar.gate_latitude['data'][start_index:end_index] 
        rlons        = radar.gate_longitude['data'][start_index:end_index] 
        ZH           = radar.fields[reflectivity_name]['data'][start_index:end_index] 
        RHOHV        = radar.fields['RHOHV']['data'][start_index:end_index] 
        
        # Filtrar cosas de RHOHV
        mask = (ZH < 20) & (RHOHV < 0.8)
        ZH[mask] = np.nan

        # Info DOW7
        radarLAT = radar.latitude['data'][0]
        radarLON = radar.longitude['data'][0]
        
    elif 'CSAPR2' in radar_name:
        reflectivity_name = 'attenuation_corrected_reflectivity_h'   
        
        radar       = pyart.io.read(rfile) 
        start_index = radar.sweep_start_ray_index['data'][0] 
        end_index   = radar.sweep_end_ray_index['data'][0]	
        
        rlats        = radar.gate_latitude['data'][start_index:end_index] 
        rlons        = radar.gate_longitude['data'][start_index:end_index] 
        ZH           = radar.fields[reflectivity_name]['data'][start_index:end_index] 
        RHOHV        = radar.fields['copol_correlation_coeff']['data'][start_index:end_index] 

        # Filtrar cosas de RHOHV
        mask = (ZH < 20) & (RHOHV < 0.8)
        ZH[mask] = np.nan
        
        # Info CSAPR2
        radarLAT = radar.latitude['data'][0]
        radarLON = radar.longitude['data'][0]
        
    temp     = wrf.g_temp.get_tk(ncfile)
    pressure = wrf.g_pressure.get_pressure(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL
    Re       = 6.3781e6
    z_level  = Re*geopo_p/(Re-geopo_p)
    lat      =  wrf.getvar( ncfile,"lat") 
    lon      =  wrf.getvar( ncfile,"lon")

    qr = mixr2massconc( np.squeeze(ncfile.variables["QRAIN"][itime,:,:,:]  ), pressure, temp )        
    qs = mixr2massconc( np.squeeze(ncfile.variables["QSNOW"][itime,:,:,:]  ), pressure, temp )        
    qi = mixr2massconc( np.squeeze(ncfile.variables["QICE"][itime,:,:,:]   ), pressure, temp )        
    qc = mixr2massconc( np.squeeze(ncfile.variables["QCLOUD"][itime,:,:,:] ), pressure, temp )       
    qg = mixr2massconc( np.squeeze(ncfile.variables["QGRAUP"][itime,:,:,:] ), pressure, temp ) 

    qr_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qr)) , z_level, axis=0)
    qs_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qs)) , z_level, axis=0)
    qg_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qg)) , z_level, axis=0)
    qi_int = integrate.trapz(np.ma.array(qi, mask=np.isnan(qi)) , z_level, axis=0)
    qc_int = integrate.trapz(np.ma.array(qc, mask=np.isnan(qc)) , z_level, axis=0)
    
    qr_int[qr_int<0.0001] = np.nan
    qs_int[qs_int<0.0001] = np.nan
    qg_int[qg_int<0.0001] = np.nan
    qi_int[qi_int<0.0001] = np.nan
    qc_int[qc_int<0.0001] = np.nan

    fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[12,2])
    pcm0 = axes[0].pcolormesh(lon, lat, qr_int, vmin=0,  vmax=12)
    axes[0].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)
    pcm1 = axes[1].pcolormesh(lon, lat, qs_int, vmin=0,  vmax=12)
    axes[1].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)
    pcm2 = axes[2].pcolormesh(lon, lat, qg_int, vmin=0,  vmax=12)
    axes[2].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)
    pcm3 = axes[3].pcolormesh(lon, lat, qc_int, vmin=0,  vmax=12)
    axes[3].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)
    pcm4 = axes[4].pcolormesh(lon, lat, qi_int, vmin=0,  vmax=1)
    axes[4].plot(np.nan, np.nan, 'darkred',  linewidth=1.5, label='45dBZ cont.')
    axes[4].legend(fontsize=9, loc='lower left')
    axes[4].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)

    cbar = plt.colorbar(pcm3, ax=axes[3], shrink=1)
    cbar.cmap.set_under('white')

    cbar = plt.colorbar(pcm4, ax=axes[4], shrink=1, label='qx [kg/m$^2$]')
    cbar.cmap.set_under('white')
    
    for i in range(5):
        axes[i].grid()
        axes[i].set_xlim([-65.5,-62])
        axes[i].set_ylim([-33.5,-31.3])

    axes[0].set_title('WSM6 q_rain')                                                  
    axes[1].set_title('WSM6 q_snow')                                                        
    axes[2].set_title('WSM6 q_grau')  
    axes[3].set_title('WSM6 q_cld')  
    axes[4].set_title('WSM6 q_ice') 

    plt.suptitle(title)
    
    for i in [120, 240]: 
        [lat_radius, lon_radius] = pyplot_rings(radarLAT,radarLON,i)
        for iax in range(5): 
            axes[iax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)  
    
    fig.savefig(savedir+'/'+title+'.png', dpi=300,transparent=False)
    #plt.close()
    
    #fig, axes = plt.subplots(nrows=1, ncols=2, constrained_layout=True,figsize=[12,2])
    #axes[0].contour(rlons, rlats, ZH, [45], colors='darkred', linewidths=1.5)
    plt.show() 
    
    
    return

#------------------------------------------------------------------------------
def plot_qxs_qitot(title, itime, ncfile, rfile, savedir, radar_name): 
    
    if 'RMA1' in radar_name:
        rlons, rlats, ZH = treat_RMA1(rfile)
        print('Note ZH here is uncorrected for attenuation')
        
        # Info RMA1
        radarLAT = -31.441389
        radarLON = -64.191944
        
    elif 'DOW6' in radar_name:
        reflectivity_name = 'DBZHC'
        
        radar       = pyart.io.read(rfile) 
        start_index = radar.sweep_start_ray_index['data'][0] 
        end_index   = radar.sweep_end_ray_index['data'][0]	
        
        rlats        = radar.gate_latitude['data'][start_index:end_index] 
        rlons        = radar.gate_longitude['data'][start_index:end_index] 
        ZH           = radar.fields[reflectivity_name]['data'][start_index:end_index] 

        # Info DOW6
        radarLAT = radar.latitude['data'][0]
        radarLON = radar.longitude['data'][0]

    elif 'DOW7' in radar_name:
        reflectivity_name = 'DBZHC'
        
        radar       = pyart.io.read(rfile) 
        start_index = radar.sweep_start_ray_index['data'][0] 
        end_index   = radar.sweep_end_ray_index['data'][0]	
        
        rlats        = radar.gate_latitude['data'][start_index:end_index] 
        rlons        = radar.gate_longitude['data'][start_index:end_index] 
        ZH           = radar.fields[reflectivity_name]['data'][start_index:end_index] 
        RHOHV        = radar.fields['RHOHV']['data'][start_index:end_index] 
        
        # Filtrar cosas de RHOHV
        mask = (ZH < 20) & (RHOHV < 0.8)
        ZH[mask] = np.nan

        # Info DOW7
        radarLAT = radar.latitude['data'][0]
        radarLON = radar.longitude['data'][0]
        
    elif 'CSAPR2' in radar_name:
        reflectivity_name = 'attenuation_corrected_reflectivity_h'   
        
        radar       = pyart.io.read(rfile) 
        start_index = radar.sweep_start_ray_index['data'][0] 
        end_index   = radar.sweep_end_ray_index['data'][0]	
        
        rlats        = radar.gate_latitude['data'][start_index:end_index] 
        rlons        = radar.gate_longitude['data'][start_index:end_index] 
        ZH           = radar.fields[reflectivity_name]['data'][start_index:end_index] 
        RHOHV        = radar.fields['copol_correlation_coeff']['data'][start_index:end_index] 

        # Filtrar cosas de RHOHV
        mask = (ZH < 20) & (RHOHV < 0.8)
        ZH[mask] = np.nan
        
        # Info CSAPR2
        radarLAT = radar.latitude['data'][0]
        radarLON = radar.longitude['data'][0]
        
    temp     = wrf.g_temp.get_tk(ncfile)
    pressure = wrf.g_pressure.get_pressure(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL
    Re       = 6.3781e6
    z_level  = Re*geopo_p/(Re-geopo_p)
    lat      =  wrf.getvar( ncfile,"lat") 
    lon      =  wrf.getvar( ncfile,"lon")

    qr = mixr2massconc( np.squeeze(ncfile.variables["QRAIN"][itime,:,:,:]  ), pressure, temp )        
    qs = mixr2massconc( np.squeeze(ncfile.variables["QSNOW"][itime,:,:,:]  ), pressure, temp )        
    qi = mixr2massconc( np.squeeze(ncfile.variables["QICE"][itime,:,:,:]   ), pressure, temp )        
    qc = mixr2massconc( np.squeeze(ncfile.variables["QCLOUD"][itime,:,:,:] ), pressure, temp )       
    qg = mixr2massconc( np.squeeze(ncfile.variables["QGRAUP"][itime,:,:,:] ), pressure, temp ) 

    qr_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qr)) , z_level, axis=0)
    qs_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qs)) , z_level, axis=0)
    qg_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qg)) , z_level, axis=0)
    qi_int = integrate.trapz(np.ma.array(qi, mask=np.isnan(qi)) , z_level, axis=0)
    qc_int = integrate.trapz(np.ma.array(qc, mask=np.isnan(qc)) , z_level, axis=0)
    
    qr_int[qr_int<0.0001] = np.nan
    qs_int[qs_int<0.0001] = np.nan
    qg_int[qg_int<0.0001] = np.nan
    qi_int[qi_int<0.0001] = np.nan
    qc_int[qc_int<0.0001] = np.nan

    fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[12,2])
    pcm0 = axes[0].pcolormesh(lon, lat, qr_int, vmin=0,  vmax=12)
    axes[0].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)
    pcm3 = axes[1].pcolormesh(lon, lat, qc_int, vmin=0,  vmax=12)
    axes[1].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)    
    pcm1 = axes[2].pcolormesh(lon, lat, (qs_int+qg_int+qi_int), vmin=0,  vmax=12)
    axes[2].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)
    #pcm2 = axes[2].pcolormesh(lon, lat, qg_int, vmin=0,  vmax=12)
    #axes[2].contour(rlons, rlats, ZH, [45], colors=(['darkred']), linewidths=1.5)
    #pcm4 = axes[4].pcolormesh(lon, lat, qi_int, vmin=0,  vmax=1)
    axes[1].plot(np.nan, np.nan, 'darkred',  linewidth=1.5, label='45dBZ cont.')
    axes[1].legend(fontsize=9, loc='upper left')

    cbar = plt.colorbar(pcm1, ax=axes[2], shrink=1)
    cbar.cmap.set_under('white')

    #cbar = plt.colorbar(pcm4, ax=axes[4], shrink=1, label='qx [kg/m$^2$]')
    #cbar.cmap.set_under('white')
    
    for i in range(3):
        axes[i].grid()
        axes[i].set_xlim([-65.5,-62])
        axes[i].set_ylim([-33.5,-31.3])

    axes[0].set_title('WSM6 q_rain')                                                  
    axes[1].set_title('WSM6 q_cld')                                                        
    axes[2].set_title('WSM6 q_ice(tot)')  
    #axes[3].set_title('WSM6 q_cld')  
    #axes[4].set_title('WSM6 q_ice') 

    plt.suptitle(title)
    fig.delaxes(axes[3])
    fig.delaxes(axes[4])
    
    for i in [120, 240]: 
        [lat_radius, lon_radius] = pyplot_rings(radarLAT,radarLON,i)
        for iax in range(3): 
            axes[iax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)  
    
    fig.savefig(savedir+'/'+title+'qi_tot.png', dpi=300,transparent=False)
    plt.close()
    
    #fig, axes = plt.subplots(nrows=1, ncols=2, constrained_layout=True,figsize=[12,2])
    #axes[0].contour(rlons, rlats, ZH, [45], colors='darkred', linewidths=1.5)
    plt.show() 
    
    
    return
#------------------------------------------------------------------------------
def plot_qxs_onlyWRF_P3(title, itime, ncfile, savedir, DOW6file, DOW7file, CSAPR2file):
    
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    
    radar         = pyart.io.read(DOW6file) 
    radarLAT_DOW6 = radar.latitude['data'][0]
    radarLON_DOW6 = radar.longitude['data'][0]
    DOW6_range = round(np.nanmax(radar.range['data'])/1000)

    radar         = pyart.io.read(DOW7file) 
    radarLAT_DOW7 = radar.latitude['data'][0]
    radarLON_DOW7 = radar.longitude['data'][0]
    DOW7_range = round(np.nanmax(radar.range['data'])/1000)
    
    radar         = pyart.io.read(CSAPR2file)     
    radarLAT_CSPR2 = radar.latitude['data'][0]
    radarLON_CSPR2 = radar.longitude['data'][0]
    CSAPR2_range = round(np.nanmax(radar.range['data'])/1000)
        
    temp     = wrf.g_temp.get_tk(ncfile)
    pressure = wrf.g_pressure.get_pressure(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL
    Re       = 6.3781e6
    z_level  = Re*geopo_p/(Re-geopo_p)
    lat      =  wrf.getvar( ncfile,"lat") 
    lon      =  wrf.getvar( ncfile,"lon")

    qr = mixr2massconc( np.squeeze(ncfile.variables["QRAIN"][itime,:,:,:]  ), pressure, temp )        
    qi = mixr2massconc( np.squeeze(ncfile.variables["QICE"][itime,:,:,:]   ), pressure, temp )        
    qc = mixr2massconc( np.squeeze(ncfile.variables["QCLOUD"][itime,:,:,:] ), pressure, temp )       

    qr_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qr)) , z_level, axis=0)
    qi_int = integrate.trapz(np.ma.array(qi, mask=np.isnan(qi)) , z_level, axis=0)
    qc_int = integrate.trapz(np.ma.array(qc, mask=np.isnan(qc)) , z_level, axis=0)
    
    qr_int[qr_int<0.0001] = np.nan
    qi_int[qi_int<0.0001] = np.nan
    qc_int[qc_int<0.0001] = np.nan

    fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[12,2])
    pcm0 = axes[0].pcolormesh(lon, lat, qr_int, vmin=0,  vmax=12)
    pcm3 = axes[1].pcolormesh(lon, lat, qc_int, vmin=0,  vmax=12)
    pcm4 = axes[2].pcolormesh(lon, lat, qi_int, vmin=0,  vmax=12)

    cbar = plt.colorbar(pcm4, ax=axes[2], shrink=1, label='qx [kg/m$^2$]')
    cbar.cmap.set_under('white')
    
    for i in range(5):
        axes[i].grid()
        axes[i].set_xlim([-65.5,-62])
        axes[i].set_ylim([-33.5,-31.3])

    axes[0].set_title('P3 q_rain')                                                  
    axes[1].set_title('P3 q_cld')  
    axes[2].set_title('P3 q_ice(tot)') 

    fig.delaxes(axes[4])
    
    plt.suptitle(title)
        
    # RMA1 
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)
    for iax in range(3): 
        axes[iax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)  
    # DOW6
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_DOW6,radarLON_DOW6,DOW6_range)
    for iax in range(3): 
        axes[iax].plot(lon_radius, lat_radius, 'darkred', linewidth=0.8)  
    # DOW7
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_DOW7,radarLON_DOW7,DOW7_range)
    for iax in range(3): 
        axes[iax].plot(lon_radius, lat_radius, 'magenta', linewidth=0.8)  
    # CSPR2
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_CSPR2,radarLON_CSPR2,CSAPR2_range)
    for iax in range(3): 
        axes[iax].plot(lon_radius, lat_radius, 'darkgreen', linewidth=0.8)  
    
    axes[3].plot(np.nan, np.nan, 'k',  linewidth=1.5, label='RMA1 (120km)')
    axes[3].plot(np.nan, np.nan, 'darkred',  linewidth=1.5, label='DOW6 ('+str(DOW6_range)+'km)')
    axes[3].plot(np.nan, np.nan, 'magenta',  linewidth=1.5, label='DOW7 ('+str(DOW7_range)+'km)')
    axes[3].plot(np.nan, np.nan, 'darkgreen',  linewidth=1.5, label='CSAPR-2 ('+str(CSAPR2_range)+'km)')
    axes[3].legend(fontsize=9, loc='upper left')
    
        
    fig.savefig(savedir+title+'.png', dpi=300,transparent=False)
    #plt.close()
    plt.show() 
    
    return

#------------------------------------------------------------------------------
def plot_qxs_onlyWRF(title, itime, ncfile, savedir, DOW6file, DOW7file, CSAPR2file):
    
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    
    radar         = pyart.io.read(DOW6file) 
    radarLAT_DOW6 = radar.latitude['data'][0]
    radarLON_DOW6 = radar.longitude['data'][0]
    DOW6_range = round(np.nanmax(radar.range['data'])/1000)

    radar         = pyart.io.read(DOW7file) 
    radarLAT_DOW7 = radar.latitude['data'][0]
    radarLON_DOW7 = radar.longitude['data'][0]
    DOW7_range = round(np.nanmax(radar.range['data'])/1000)
    
    radar         = pyart.io.read(CSAPR2file)     
    radarLAT_CSPR2 = radar.latitude['data'][0]
    radarLON_CSPR2 = radar.longitude['data'][0]
    CSAPR2_range = round(np.nanmax(radar.range['data'])/1000)
        
    temp     = wrf.g_temp.get_tk(ncfile)
    pressure = wrf.g_pressure.get_pressure(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL
    Re       = 6.3781e6
    z_level  = Re*geopo_p/(Re-geopo_p)
    lat      =  wrf.getvar( ncfile,"lat") 
    lon      =  wrf.getvar( ncfile,"lon")

    qr = mixr2massconc( np.squeeze(ncfile.variables["QRAIN"][itime,:,:,:]  ), pressure, temp )        
    qs = mixr2massconc( np.squeeze(ncfile.variables["QSNOW"][itime,:,:,:]  ), pressure, temp )        
    qi = mixr2massconc( np.squeeze(ncfile.variables["QICE"][itime,:,:,:]   ), pressure, temp )        
    qc = mixr2massconc( np.squeeze(ncfile.variables["QCLOUD"][itime,:,:,:] ), pressure, temp )       
    qg = mixr2massconc( np.squeeze(ncfile.variables["QGRAUP"][itime,:,:,:] ), pressure, temp ) 

    qr_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qr)) , z_level, axis=0)
    qs_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qs)) , z_level, axis=0)
    qg_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qg)) , z_level, axis=0)
    qi_int = integrate.trapz(np.ma.array(qi, mask=np.isnan(qi)) , z_level, axis=0)
    qc_int = integrate.trapz(np.ma.array(qc, mask=np.isnan(qc)) , z_level, axis=0)
    
    qr_int[qr_int<0.0001] = np.nan
    qs_int[qs_int<0.0001] = np.nan
    qg_int[qg_int<0.0001] = np.nan
    qi_int[qi_int<0.0001] = np.nan
    qc_int[qc_int<0.0001] = np.nan

    fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[12,2])
    pcm0 = axes[0].pcolormesh(lon, lat, qr_int, vmin=0,  vmax=12)
    pcm1 = axes[1].pcolormesh(lon, lat, qs_int, vmin=0,  vmax=12)
    pcm2 = axes[2].pcolormesh(lon, lat, qg_int, vmin=0,  vmax=12)
    pcm3 = axes[3].pcolormesh(lon, lat, qc_int, vmin=0,  vmax=12)
    pcm4 = axes[4].pcolormesh(lon, lat, qi_int, vmin=0,  vmax=1)

    cbar = plt.colorbar(pcm3, ax=axes[3], shrink=1)
    cbar.cmap.set_under('white')

    cbar = plt.colorbar(pcm4, ax=axes[4], shrink=1, label='qx [kg/m$^2$]')
    cbar.cmap.set_under('white')
    
    for i in range(5):
        axes[i].grid()
        axes[i].set_xlim([-65.5,-62])
        axes[i].set_ylim([-33.5,-31.3])

    axes[0].set_title('WSM6 q_rain')                                                  
    axes[1].set_title('WSM6 q_snow')                                                        
    axes[2].set_title('WSM6 q_grau')  
    axes[3].set_title('WSM6 q_cld')  
    axes[4].set_title('WSM6 q_ice') 

    plt.suptitle(title)
        
    # RMA1 
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)
    for iax in range(5): 
        axes[iax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)  
    # DOW6
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_DOW6,radarLON_DOW6,DOW6_range)
    for iax in range(5): 
        axes[iax].plot(lon_radius, lat_radius, 'darkred', linewidth=0.8)  
    # DOW7
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_DOW7,radarLON_DOW7,DOW7_range)
    for iax in range(5): 
        axes[iax].plot(lon_radius, lat_radius, 'magenta', linewidth=0.8)  
    # CSPR2
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_CSPR2,radarLON_CSPR2,CSAPR2_range)
    for iax in range(5): 
        axes[iax].plot(lon_radius, lat_radius, 'darkgreen', linewidth=0.8)  
    
    axes[4].plot(np.nan, np.nan, 'k',  linewidth=1.5, label='RMA1 (120km)')
    axes[4].plot(np.nan, np.nan, 'darkred',  linewidth=1.5, label='DOW6 ('+str(DOW6_range)+'km)')
    axes[4].plot(np.nan, np.nan, 'magenta',  linewidth=1.5, label='DOW7 ('+str(DOW7_range)+'km)')
    axes[4].plot(np.nan, np.nan, 'darkgreen',  linewidth=1.5, label='CSAPR-2 ('+str(CSAPR2_range)+'km)')
    axes[4].legend(fontsize=9, loc='lower left')
    
        
    fig.savefig(savedir+'/'+title+'.png', dpi=300,transparent=False)
    #plt.close()
    
    #fig, axes = plt.subplots(nrows=1, ncols=2, constrained_layout=True,figsize=[12,2])
    #axes[0].contour(rlons, rlats, ZH, [45], colors='darkred', linewidths=1.5)
    plt.show() 
    
    return

#------------------------------------------------------------------------------
def general_plot(radar1, radar2, radar3, radar4, savedir, timetit):
    
    # define 
    cmap = VP.colormaps('ref')
    vmax = 60
    vmin = 0
    
    fig = plt.figure(figsize=(8,8)) 
    # radar1 (rma1)
    radar       = pyart.io.read(radar1) 
    start_index = radar.sweep_start_ray_index['data'][0] 
    end_index   = radar.sweep_end_ray_index['data'][0]	
    rlats        = radar.gate_latitude['data'][start_index:end_index] 
    rlons        = radar.gate_longitude['data'][start_index:end_index] 
    ZH           = radar.fields['TH']['data'][start_index:end_index] 
    pcm1=plt.pcolormesh(rlons,rlats, ZH, cmap=cmap, vmax=vmax, vmin=vmin)
    plt.xlim([-65.5,-62])
    plt.ylim([-33.5,-31.3])   
    plt.title('RMA1 at '+timetit)
    plt.colorbar(pcm1,label='(dBZ)')
    fig.savefig(savedir+'general_radar_deployment_RMA1'+timetit+'.png', dpi=300, transparent=False)
    plt.show()
    
    fig = plt.figure(figsize=(8,8)) 
    # radar2 (DOW6)
    radar       = pyart.io.read(radar2) 
    start_index = radar.sweep_start_ray_index['data'][0] 
    end_index   = radar.sweep_end_ray_index['data'][0]	
    rlats        = radar.gate_latitude['data'][start_index:end_index] 
    rlons        = radar.gate_longitude['data'][start_index:end_index] 
    ZH           = radar.fields['DBZHC']['data'][start_index:end_index] 
    pcm1=plt.pcolormesh(rlons,rlats, ZH, cmap=cmap, vmax=vmax, vmin=vmin)
    plt.xlim([-65.5,-62])
    plt.ylim([-33.5,-31.3])   
    plt.title('DOW6 at '+timetit)
    plt.colorbar(pcm1,label='(dBZ)')
    fig.savefig(savedir+'general_radar_deployment_DOW6'+timetit+'.png', dpi=300,transparent=False)
    plt.show()
    
    fig = plt.figure(figsize=(8,8)) 
    # radar3 (DOW7)
    radar       = pyart.io.read(radar3) 
    start_index = radar.sweep_start_ray_index['data'][0] 
    end_index   = radar.sweep_end_ray_index['data'][0]	
    rlats        = radar.gate_latitude['data'][start_index:end_index] 
    rlons        = radar.gate_longitude['data'][start_index:end_index] 
    ZH           = radar.fields['DBZHC']['data'][start_index:end_index] 
    pcm1=plt.pcolormesh(rlons,rlats,ZH, cmap=cmap, vmax=vmax, vmin=vmin)
    plt.xlim([-65.5,-62])
    plt.ylim([-33.5,-31.3])   
    plt.title('DOW7 at '+timetit)
    plt.colorbar(pcm1,label='(dBZ)')
    fig.savefig(savedir+'general_radar_deployment_DOW7'+timetit+'.png', dpi=300,transparent=False)
    plt.show()
    
    fig = plt.figure(figsize=(8,8)) 
    # radar4 (CSPR2)
    radar       = pyart.io.read(radar4) 
    start_index = radar.sweep_start_ray_index['data'][0] 
    end_index   = radar.sweep_end_ray_index['data'][0]	
    rlats        = radar.gate_latitude['data'][start_index:end_index] 
    rlons        = radar.gate_longitude['data'][start_index:end_index] 
    ZH           = radar.fields['attenuation_corrected_reflectivity_h' ]['data'][start_index:end_index] 
    pcm1=plt.pcolormesh(rlons,rlats,ZH, cmap=cmap, vmax=vmax, vmin=vmin)
    plt.xlim([-65.5,-62])
    plt.ylim([-33.5,-31.3])   
    plt.title('CSPR2 at '+timetit)    
    plt.colorbar(pcm1,label='(dBZ)')
    fig.savefig(savedir+'general_radar_deployment_CSAPR2'+timetit+'.png', dpi=300,transparent=False)
    plt.show()
    
    return



#------------------------------------------------------------------------------
def make_radar_gif(radar_folder, wildfiles, time1infile, all_times_minutes, radar_name, output_gif_path, tempsavedir): 

    file_pattern = os.path.join(radar_folder, wildfiles)
    file_list    = sorted(glob.glob(file_pattern))
    times        = [filename[time1infile:time1infile+4] for filename in file_list]

    # Step 1: Round each time to the nearest 10 minutes
    rounded_minutes = [int(round(hhmm_to_minutes(t) / 10) * 10) for t in times]

    # Step 2: Create an array with NaN values for all intervals
    time_array = np.full(len(all_times_minutes), np.nan, dtype=object)
    filename_array = np.full(len(all_times_minutes), np.nan, dtype=object)

    # Step 3: Fill in times (and filenames) where data is available
    for i, rounded_time in enumerate(rounded_minutes):
        index = np.where(all_times_minutes == rounded_time)[0]
        if index.size > 0:
            time_array[index[0]] = minutes_to_hhmm(rounded_time)
            filename_array[index[0]] = file_list[i]

    # Convert all_times_minutes back to HHMM format for easy reading
    all_times_hhmm = [minutes_to_hhmm(m) for m in all_times_minutes]
    
    if 'RMA1' in radar_name:
        TH_name = 'TH'
    elif 'DOW6' in radar_name:
        TH_name = 'DBZHC'
    elif 'DOW7' in radar_name:
        TH_name = 'DBZHC'
    elif 'CSAPR2' in radar_name:
        TH_name = 'attenuation_corrected_reflectivity_h'   
        
    images = []
    #for ii, filename in enumerate(file_list):

    for i, ii in enumerate(time_array):
        
        if np.isnan(float(time_array[i])) == True:     # np.nan and there is no observation for that time. i.e. empty figure
            
            fig = plt.figure(figsize=(8,8)) 
            plt.xlim([-65.5,-62])
            plt.ylim([-33.5,-31.3])   
            plt.title(radar_name+' at '+all_times_hhmm[i])
    
            # Save the figure as a temporary image file
            temp_filename = f"{all_times_hhmm[i]}.png"
            plt.savefig(tempsavedir+radar_name+temp_filename) 
            images.append(Image.open(tempsavedir+radar_name+temp_filename))
            plt.close(fig)  # Close the figure to free memory 
        
        else: 
        
            fig = plt.figure(figsize=(8,8)) 
    
            radar       = pyart.io.read(filename_array[i]) 
            start_index = radar.sweep_start_ray_index['data'][0] 
            end_index   = radar.sweep_end_ray_index['data'][0]	
            rlats        = radar.gate_latitude['data'][start_index:end_index] 
            rlons        = radar.gate_longitude['data'][start_index:end_index] 
            ZH           = radar.fields[TH_name]['data'][start_index:end_index] 
            
            # Info radar
            radarLAT = radar.latitude['data'][0]
            radarLON = radar.longitude['data'][0]
            
            # PLOT
            pcm1 = plt.pcolormesh(rlons,rlats, ZH, cmap=VP.colormaps('ref'), vmax=60, vmin=0)    
            plt.xlim([-65.5,-62])
            plt.ylim([-33.5,-31.3])   
            plt.title(radar_name+' at '+all_times_hhmm[i])
    
            # Save the figure as a temporary image file
            temp_filename = f"{all_times_hhmm[i]}.png"
            plt.savefig(tempsavedir+radar_name+temp_filename) 
            images.append(Image.open(tempsavedir+radar_name+temp_filename))
            plt.close(fig)  # Close the figure to free memory

    # Step 3: Save as a GIF
    output_gif_path = output_gif_path+radar_name+'.gif'
    images[0].save(output_gif_path, save_all=True, append_images=images[1:], duration=1000, loop=0)
    print(f"GIF saved at {output_gif_path}")
    
    return 

#------------------------------------------------------------------------------
# explore variables in wrf netcdf file 
def expand_var_dict(var_dict, group):
    for var_key, _ in group.variables.items():
      if group.name not in var_dict.keys():
        var_dict[group.name] = list()
      var_dict[group.name].append(var_key)
    
    for _, sub_group in group.groups.items():
      expand_var_dict(var_dict, sub_group)
      
      
#------------------------------------------------------------------------------
# explore variables and descriptions
def print_wrf_variables(wrf_ncfile):
    
    all_vars = []
    all_descriptions = []
    
    
    for var_key, _ in wrf_ncfile.variables.items():
        if var_key != 'Times':
            print(var_key)
            all_vars.append(var_key)
            all_descriptions.append(wrf_ncfile[var_key].description)
            
    return all_vars, all_descriptions

#------------------------------------------------------------------------------
def return_accums(ncfile):
 
    RAINC =   np.squeeze(ncfile.variables["RAINC"][0,:,:])      
    RAINSH =  np.squeeze(ncfile.variables["RAINSH"][0,:,:])    
    RAINNC =  np.squeeze(ncfile.variables["RAINNC"][0,:,:])      
    SNOWNC =  np.squeeze(ncfile.variables["SNOWNC"][0,:,:])       
    GRAUPELNC = np.squeeze(ncfile.variables["GRAUPELNC"][0,:,:])
    HAILNC = np.squeeze(ncfile.variables["HAILNC"][0,:,:]) 
    CLDFRA = np.squeeze(ncfile.variables["CLDFRA"][0,:,:,:])
    
    lat      =  wrf.getvar( ncfile,"lat") 
    lon      =  wrf.getvar( ncfile,"lon")
    
    return RAINC, RAINSH, RAINNC, SNOWNC, GRAUPELNC, HAILNC, CLDFRA, lon, lat

#------------------------------------------------------------------------------
def plot_accums(wrf_file_wsm6, wrf_file_p3, title, savedir, DOW6file, DOW7file, CSAPR2file):
    
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    
    radar         = pyart.io.read(DOW6file) 
    radarLAT_DOW6 = radar.latitude['data'][0]
    radarLON_DOW6 = radar.longitude['data'][0]
    DOW6_range = round(np.nanmax(radar.range['data'])/1000)

    radar         = pyart.io.read(DOW7file) 
    radarLAT_DOW7 = radar.latitude['data'][0]
    radarLON_DOW7 = radar.longitude['data'][0]
    DOW7_range = round(np.nanmax(radar.range['data'])/1000)
    
    radar         = pyart.io.read(CSAPR2file)     
    radarLAT_CSPR2 = radar.latitude['data'][0]
    radarLON_CSPR2 = radar.longitude['data'][0]
    CSAPR2_range = round(np.nanmax(radar.range['data'])/1000)
    
    [RAINC, RAINSH, RAINNC, SNOWNC, GRAUPELNC, HAILNC, CLDFRA, lon, lat] = return_accums(wrf_file_wsm6)
    [p3RAINC, p3RAINSH, p3RAINNC, p3SNOWNC, p3GRAUPELNC, p3HAILNC, p3CLDFRA, p3lon, p3lat] = return_accums(wrf_file_p3)
    

    fig, axes = plt.subplots(nrows=1, ncols=3, constrained_layout=True,figsize=[12,6])
    pcm0 = axes[0].pcolormesh(lon, lat, RAINNC, vmin=0,  vmax=150)    
    pcm1 = axes[1].pcolormesh(p3lon, p3lat, p3RAINNC, vmin=0,  vmax=150)
    
    axes[0].set_title('RAINNC WSM6')  
    axes[1].set_title('RAINNC P3 3MOM LF')  
    axes[2].set_title('RAINNC P3 3MOM no LF')  


    cbar = plt.colorbar(pcm1, ax=axes[1], shrink=1,  label='Rain accum [mm]')
    cbar.cmap.set_under('white')

    #cbar = plt.colorbar(pcm2, ax=axes[2], shrink=1,  label='Rain accum [mm]')
    #cbar.cmap.set_under('white')
    
    for i in range(3):
        axes[i].grid()
        axes[i].set_xlim([-65.5,-62])
        axes[i].set_ylim([-33.5,-31.3])
        
    plt.suptitle(title)

    # radar rings    
    # RMA1 
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)
    for iax in range(3): 
        axes[iax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)  
    # DOW6
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_DOW6,radarLON_DOW6,DOW6_range)
    for iax in range(3): 
        axes[iax].plot(lon_radius, lat_radius, 'darkred', linewidth=0.8)  
    # DOW7
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_DOW7,radarLON_DOW7,DOW7_range)
    for iax in range(3): 
        axes[iax].plot(lon_radius, lat_radius, 'magenta', linewidth=0.8)  
    # CSPR2
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_CSPR2,radarLON_CSPR2,CSAPR2_range)
    for iax in range(3): 
        axes[iax].plot(lon_radius, lat_radius, 'darkgreen', linewidth=0.8)  
    
    axes[2].plot(np.nan, np.nan, 'k',  linewidth=1.5, label='RMA1 (120km)')
    axes[2].plot(np.nan, np.nan, 'darkred',  linewidth=1.5, label='DOW6 ('+str(DOW6_range)+'km)')
    axes[2].plot(np.nan, np.nan, 'magenta',  linewidth=1.5, label='DOW7 ('+str(DOW7_range)+'km)')
    axes[2].plot(np.nan, np.nan, 'darkgreen',  linewidth=1.5, label='CSAPR-2 ('+str(CSAPR2_range)+'km)')
    axes[2].legend(fontsize=9, loc='lower left')
    
    plt.show()    
    fig.savefig(savedir+'comparePrecipAccum'+title+'.png', dpi=300,transparent=False)
    #plt.close()

    return

#------------------------------------------------------------------------------
def return_qxs_WSM6(ncfile, domainlons, domainlats, z_interp):

    # Need to interpolate first to common z_level grid 
    z_interp =  [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 
                 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20]    
    
    temp     = wrf.g_temp.get_tk(ncfile)
    pressure = wrf.g_pressure.get_pressure(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL
    Re       = 6.3781e6
    z_level  = Re*geopo_p/(Re-geopo_p)
    lats      =  wrf.getvar( ncfile,"lat") 
    lons      =  wrf.getvar( ncfile,"lon")

    qr = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QRAIN"][0,:,:,:]  ),  "ght_msl", z_interp)          
    qs = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QSNOW"][0,:,:,:]  ),  "ght_msl", z_interp)          
    qi = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QICE"][0,:,:,:]   ),  "ght_msl", z_interp)          
    qc = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QCLOUD"][0,:,:,:] ),  "ght_msl", z_interp)         
    qg = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QGRAUP"][0,:,:,:] ),  "ght_msl", z_interp)   
    
    
    # Create a mask for the domain
    lat_mask = (lats >= domainlats[0]) & (lats <= domainlats[1])
    lon_mask = (lons >= domainlons[0]) & (lons <= domainlons[1])
    mask = lat_mask & lon_mask 
    
    # Apply the mask to the variable q
    qr_masked = np.where(mask, qr, np.nan)  # Replace values outside the domain with NaN
    qi_masked = np.where(mask, qi, np.nan)  # Replace values outside the domain with NaN
    qc_masked = np.where(mask, qc, np.nan)  # Replace values outside the domain with NaN
    qs_masked = np.where(mask, qs, np.nan)  # Replace values outside the domain with NaN
    qg_masked = np.where(mask, qg, np.nan)  # Replace values outside the domain with NaN
    z_level_masked = np.where(mask, z_level, np.nan)    
    
    return qr_masked, qs_masked, qi_masked, qc_masked, qg_masked, qg_masked+qi_masked+qs_masked, z_level_masked

#------------------------------------------------------------------------------
def return_qxs_p3(ncfile, domainlons, domainlats, z_interp):
    
    temp     = wrf.g_temp.get_tk(ncfile)
    pressure = wrf.g_pressure.get_pressure(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL
    Re       = 6.3781e6
    z_level  = Re*geopo_p/(Re-geopo_p)
    lats      =  wrf.getvar( ncfile,"lat") 
    lons      =  wrf.getvar( ncfile,"lon")

    qr = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QRAIN"][0,:,:,:]), "ght_msl", z_interp)   
    qi = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QICE"][0,:,:,:]),  "ght_msl", z_interp)    
    qc = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QCLOUD"][0,:,:,:]), "ght_msl", z_interp)   
    
    # Create a mask for the domain
    lat_mask = (lats >= domainlats[0]) & (lats <= domainlats[1])
    lon_mask = (lons >= domainlons[0]) & (lons <= domainlons[1])
    
    mask = lat_mask & lon_mask
    
    # Apply the mask to the variable q
    qr_masked = np.where(mask, qr, np.nan)  # Replace values outside the domain with NaN
    qi_masked = np.where(mask, qi, np.nan)  # Replace values outside the domain with NaN
    qc_masked = np.where(mask, qc, np.nan)  # Replace values outside the domain with NaN
    z_level_masked = np.where(mask, z_level, np.nan)
    
    return qr_masked, qi_masked, qc_masked, z_level_masked

#------------------------------------------------------------------------------
def plot_domain_qxs(wrf_file_wsm6, wrf_file_p3, title, savedir):
    
    # Domain    
    domainlons = [-65.5,-62]
    domainlats = [-33.5,-31.3] 
    
    # Need to interpolate first to common z_level grid 
    z_interp =  [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 
                 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20] 
    
    [qr, qs, qi, qc, qg, qitot, z_lev] = return_qxs_WSM6(wrf_file_wsm6, domainlons, domainlats, z_interp)
    [qr_p3, qi_p3, qc_p3, z_lev_p3] = return_qxs_p3(wrf_file_p3, domainlons, domainlats, z_interp)
    
    
    # Calculate the average inside the domain, ignoring NaN values
    #average_q = np.nanmean(q_masked)
   
    fig, axes = plt.subplots(nrows=1, ncols=3, constrained_layout=True,figsize=[12,6])
    axes[0].plot(np.nanmean(qr*1000, axis = (1, 2)), z_interp, color='darkred', linewidth=1.2, label='qr')    
    axes[0].plot(np.nanmean(qs*1000, axis = (1, 2)), z_interp, color='darkblue', linewidth=1.2, label='qs')    
    axes[0].plot(np.nanmean(qg*1000, axis = (1, 2)), z_interp, color='darkgreen', linewidth=1.2, label='qg')    
    axes[0].plot(np.nanmean(qi*1000, axis = (1, 2)), z_interp, color='indigo', linewidth=1.2, label='qi')    
    axes[0].plot(np.nanmean(qc*1000, axis = (1, 2)), z_interp, color='dodgerblue', linewidth=1.2, label='qc')    
    axes[0].plot(np.nanmean(qitot*1000, axis = (1, 2)), z_interp, color='k', linewidth=1.2, label='qg+qs+qi')  
    axes[0].set_xlabel('Mixing ratio (g/kg)')
    axes[0].set_ylabel('Height (km)')
    axes[0].set_title('WSM6')
    axes[0].legend()
    axes[0].grid(True)
    axes[0].set_xlim([0, 0.4])
    axes[0].set_ylim([0,20])

    axes[1].plot(np.nanmean(qr_p3*1000, axis = (1, 2)), z_interp, color='darkred', linewidth=1.2, label='qr')    
    axes[1].plot(np.nanmean(qi_p3*1000, axis = (1, 2)), z_interp, color='k', linewidth=1.2, label='qi (tot)')    
    axes[1].plot(np.nanmean(qc_p3*1000, axis = (1, 2)), z_interp, color='dodgerblue', linewidth=1.2, label='qc')    
    axes[1].set_xlabel('Mixing ratio (g/kg)')
    axes[1].set_ylabel('Height (km)')
    axes[1].set_title('P3 3MOM LF')
    axes[1].legend()
    axes[1].set_xlim([0, 0.35])
    axes[1].grid(True)
    axes[1].set_ylim([0,20])
                
    plt.suptitle('Domain average vertical profiles at' +title)    
    plt.show()    
    fig.savefig(savedir+'compareDomain_qxs'+title+'.png', dpi=300,transparent=False)

    return

#------------------------------------------------------------------------------
def get_domainAccum(folderWRFOUT, domainlons, domainlats, timeNr):
    
    # Domain-average evolution of precip rate
    file_list = sorted(glob.glob(folderWRFOUT))
    times          = [filename[timeNr:timeNr+16] for filename in file_list]

    accumPrecip_average = [] 

    for filename in file_list:
        wrf_file = Dataset(filename,'r')
        RAINNC =  np.squeeze(wrf_file.variables["RAINNC"][0,:,:])   
        lats   =  wrf.getvar( wrf_file,"lat") 
        lons   =  wrf.getvar( wrf_file,"lon")
        # Create a mask for the domain
        lat_mask = (lats >= domainlats[0]) & (lats <= domainlats[1])
        lon_mask = (lons >= domainlons[0]) & (lons <= domainlons[1])    
        mask = lat_mask & lon_mask
        # Apply the mask to the variable q
        RAINNC_masked = np.where(mask, RAINNC, np.nan)  # Replace values outside the domain with NaN    
        accumPrecip_average.append( np.nanmean(RAINNC_masked) )  
    
    return accumPrecip_average, times
    



#------------------------------------------------------------------------------
def find_0( arr, target ):
    
    # Find the indices and values closest to 273 in each column
    diff =  np.abs(arr.data - target)
    closest_indices = np.nanargmin(diff, axis=0)
    ### closest_values = arr.data[closest_indices, np.arange(arr.shape[1])]
    #closest_values = height_arr.data[closest_indices, np.arange(arr.shape[1])]

    return closest_indices

def find_00( arr, target, arrZ):
    
    x_dim = arr.data.shape[1]
    y_dim = arr.data.shape[2]
    
    output_array_T = np.empty((x_dim, y_dim))
    output_array_Z = np.empty((x_dim, y_dim))
    
    # Iterate over each (x, y) position to find the closest value along the z-axis
    for i in range(x_dim):
        for j in range(y_dim):
            # Find the index in z-dimension closest to the target value
            closest_index = np.argmin(np.abs(arr.data[:, i, j] - target))
            # Store the value at this index
            output_array_T[i, j] = arr.data[closest_index, i, j]
            output_array_Z[i, j] = arrZ.data[closest_index, i, j]

    return output_array_T, output_array_Z

#------------------------------------------------------------------------------
def plot_common_transect(ncfile1, ncfile2, title, savedir):

    
    # Need to interpolate first to common z_level grid 
    z_interp = 1e3*np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 
                 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20])
    
    Re       = 6.3781e6
    # ncfile1: WSM6
    temp_WSM6     = wrf.g_temp.get_tk(ncfile1)
    geopo_p_WSM6  = wrf.g_geoht.get_height(ncfile1) # geopotential height as Mean Sea Level (MSL)
    z_WSM6  = Re*geopo_p_WSM6/(Re-geopo_p_WSM6)

    # ncfile2: P3 3MOM LF
    temp_P3_3mom_LF_1     = wrf.g_temp.get_tk(ncfile2)
    geopo_p_P3_3mom_LF_1  = wrf.g_geoht.get_height(ncfile2) # geopotential height as Mean Sea Level (MSL)
    z_P3_3mom_LF  = Re*geopo_p_P3_3mom_LF_1/(Re-geopo_p_P3_3mom_LF_1)

    #z_WSM6       = wrf.getvar(ncfile1, "z")
    #z_P3_3mom_LF = wrf.getvar(ncfile2, "z")
    
    # Create the start point and end point for the cross section for WSM6
    start_point_WSM6_1 = wrf.CoordPair(lat=y_wsm6[0], lon=x_wsm6[0])
    end_point_WSM6_1 = wrf.CoordPair(lat=y_wsm6[1], lon=x_wsm6[1])
    start_point_WSM6_2 = wrf.CoordPair(lat=y_wsm6[2], lon=x_wsm6[2])
    end_point_WSM6_2 = wrf.CoordPair(lat=y_wsm6[3], lon=x_wsm6[3])

    # Create the start point and end point for the cross section for P3_3MOM_LF
    start_point_P3_3mom_LF_1 = wrf.CoordPair(lat=y_P3_3mom_LF[0], lon=x_P3_3mom_LF[0])
    end_point_P3_3mom_LF_1 = wrf.CoordPair(lat=y_P3_3mom_LF[1], lon=x_P3_3mom_LF[1])
    start_point_P3_3mom_LF_2 = wrf.CoordPair(lat=y_P3_3mom_LF[2], lon=x_P3_3mom_LF[2])
    end_point_P3_3mom_LF_2 = wrf.CoordPair(lat=y_P3_3mom_LF[3], lon=x_P3_3mom_LF[3])

    [qr_WSM6, qs_WSM6, qi_WSM6, qc_WSM6, qg_WSM6, 
     qr_int_WSM6, qs_int_WSM6, qi_int_WSM6, qc_int_WSM6, qg_int_WSM6] = get_q_ints6(ncfile1)
    qi_sum_WSM6 = qi_WSM6 + qi_WSM6 + qg_WSM6
    
    #
    [qr_P3_3mom_LF, qi_P3_3mom_LF, qc_P3_3mom_LF, 
     qr_int_P3_3mom_LF, qi_int_P3_3mom_LF, qc_int_P3_3mom_LF] = get_q_ints3(ncfile2)

    
    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section.
    qitot_cross_WSM6_1 = wrf.vertcross(qi_sum_WSM6, z_WSM6, levels=z_interp, wrfin=ncfile1, start_point=start_point_WSM6_1,
                       end_point=end_point_WSM6_1, latlon=True, meta=True)
    qitot_cross_WSM6_2 = wrf.vertcross(qi_sum_WSM6, z_WSM6,  levels=z_interp, wrfin=ncfile1, start_point=start_point_WSM6_2,
                       end_point=end_point_WSM6_2, latlon=True, meta=True)
    #
    qitot_cross_P3_3mom_LF_1 = wrf.vertcross(qi_P3_3mom_LF, z_P3_3mom_LF, levels=z_interp, wrfin=ncfile2, start_point=start_point_P3_3mom_LF_1,
                       end_point=end_point_P3_3mom_LF_1, latlon=True, meta=True)
    qitot_cross_P3_3mom_LF_2 = wrf.vertcross(qi_P3_3mom_LF, z_P3_3mom_LF,  levels=z_interp, wrfin=ncfile2, start_point=start_point_P3_3mom_LF_2,
                       end_point=end_point_P3_3mom_LF_2, latlon=True, meta=True)

    #
    fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True,figsize=[16,8])
    wsm6_contours_1 = axes[0,0].contourf( np.arange(0, qitot_cross_WSM6_1.shape[-1], 1), 
                                          wrf.to_np(qitot_cross_WSM6_1.coords["vertical"]), 
                                          wrf.to_np(qitot_cross_WSM6_1), cmap=get_cmap("viridis"), vmin=0.0001, vmax=0.005)
    
    wsm6_contours_2 = axes[1,0].contourf( np.arange(0, qitot_cross_WSM6_2.shape[-1], 1), 
                                          wrf.to_np(qitot_cross_WSM6_2.coords["vertical"]), 
                                          wrf.to_np(qitot_cross_WSM6_2), cmap=get_cmap("viridis"), vmin=0.0002, vmax=0.005)

    P3_3mom_LF_contours_1 = axes[0,1].contourf( np.arange(0, qitot_cross_P3_3mom_LF_1.shape[-1], 1), 
                                                wrf.to_np(qitot_cross_P3_3mom_LF_1.coords["vertical"]), 
                                                wrf.to_np(qitot_cross_P3_3mom_LF_1), cmap=get_cmap("viridis"), vmin=0.0002, vmax=0.005)
    
    P3_3mom_LF_contours_2 = axes[1,1].contourf( np.arange(0, qitot_cross_P3_3mom_LF_2.shape[-1], 1), 
                                                wrf.to_np(qitot_cross_P3_3mom_LF_2.coords["vertical"]),  
                                                wrf.to_np(qitot_cross_P3_3mom_LF_2), cmap=get_cmap("viridis"), vmin=0.0002, vmax=0.005)

    axes[0,0].set_title('WSM6 transect Nr. 1')
    axes[1,0].set_title('WSM6 transect Nr. 2')
    axes[0,1].set_title('P3 3mom LF transect Nr. 1')
    axes[1,1].set_title('P3 3mom LF transect Nr. 2')
    
    cbar =plt.colorbar(wsm6_contours_1, ax=axes[0,0])
    cbar.cmap.set_under('white')
    cbar =plt.colorbar(wsm6_contours_2, ax=axes[1,0])
    cbar.cmap.set_under('white')
    cbar =plt.colorbar(P3_3mom_LF_contours_1, ax=axes[0,1])
    cbar.cmap.set_under('white')
    cbar =plt.colorbar(P3_3mom_LF_contours_2, ax=axes[1,1])
    cbar.cmap.set_under('white')
    
    #--------------------------------------------------------------------------     # now add 0C level    
    [output_array_T, output_array_Z] = find_00(temp_WSM6, 273, z_WSM6)
    HGT = wrf.getvar(ncfile1, "HGT", timeidx=-1)
    HGT.data = output_array_Z
    isoline_hgt = wrf.interpline(HGT, wrfin=ncfile1, start_point=start_point_WSM6_1, end_point=end_point_WSM6_1)
    xs = np.arange(0, qitot_cross_WSM6_1.shape[-1], 1)
    axes[0,0].plot(xs, wrf.to_np(isoline_hgt), '-r', linewidth=1.2)

    isoline_hgt = wrf.interpline(HGT, wrfin=ncfile1, start_point=start_point_WSM6_2, end_point=end_point_WSM6_2)
    xs = np.arange(0, qitot_cross_WSM6_2.shape[-1], 1)
    axes[1,0].plot(xs, wrf.to_np(isoline_hgt), '-r', linewidth=1.2)

    [output_array_T, output_array_Z] = find_00(temp_P3_3mom_LF_1, 273, z_P3_3mom_LF)
    HGT = wrf.getvar(ncfile2, "HGT", timeidx=-1)
    HGT.data = output_array_Z
    isoline_hgt = wrf.interpline(HGT, wrfin=ncfile2, start_point=start_point_P3_3mom_LF_1, end_point=end_point_P3_3mom_LF_1)
    xs = np.arange(0, qitot_cross_P3_3mom_LF_1.shape[-1], 1)
    axes[0,1].plot(xs, wrf.to_np(isoline_hgt), '-r', linewidth=1.2)

    isoline_hgt = wrf.interpline(HGT, wrfin=ncfile2, start_point=start_point_P3_3mom_LF_2, end_point=end_point_P3_3mom_LF_2)
    xs = np.arange(0, qitot_cross_P3_3mom_LF_2.shape[-1], 1)
    axes[1,1].plot(xs, wrf.to_np(isoline_hgt), '-r', linewidth=1.2)
    
    #--------------------------------------------------------------------------
    # Set the y-ticks to be height.
    #------ esta la info en la meta data dde qitot_cross_WSM6_1, etc
    #vert_vals = wrf.to_np(qitot_cross_WSM6_1.coords["vertical"])
    #v_ticks = np.arange(vert_vals.shape[0])
    #axes[0,0].set_yticks(v_ticks[::3])
    #axes[0,0].set_yticklabels(vert_vals[::3], fontsize=12)
    axes[0,0].set_ylabel("Height (m)", fontsize=12)

    # Set the y-ticks to be height.
    #vert_vals = wrf.to_np(qitot_cross_WSM6_2.coords["vertical"])
    #v_ticks = np.arange(vert_vals.shape[0])
    #axes[1,0].set_yticks(v_ticks[::3])
    #axes[1,0].set_yticklabels(vert_vals[::3], fontsize=12)
    axes[1,0].set_ylabel("Height (m)", fontsize=12)
    
    # Set the y-ticks to be height.
    #vert_vals = wrf.to_np(qitot_cross_P3_3mom_LF_1.coords["vertical"])
    #v_ticks = np.arange(vert_vals.shape[0])
    #axes[0,1].set_yticks(v_ticks[::3])
    #axes[0,1].set_yticklabels(vert_vals[::3], fontsize=12)
    axes[0,1].set_ylabel("Height (m)", fontsize=12)
        
    # Set the y-ticks to be height.
    #vert_vals = wrf.to_np(qitot_cross_P3_3mom_LF_2.coords["vertical"])
    #v_ticks = np.arange(vert_vals.shape[0])
    #axes[1,1].set_yticks(v_ticks[::3])
    #axes[1,1].set_yticklabels(vert_vals[::3], fontsize=12)
    axes[1,1].set_ylabel("Height (m)", fontsize=12)    
    
    #--------------------------------------------------------------------------
    # Set the X-ticks to be LATITUDE.
    coord_pairs = wrf.to_np(qitot_cross_WSM6_1.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [round(pair.lat,2) for pair in wrf.to_np(coord_pairs)]
    axes[0,0].set_xticks(x_ticks[::10])
    axes[0,0].set_xticklabels(x_labels[::10],fontsize=12)                                 
    axes[0,0].set_xlabel("Latitude", fontsize=12)    

    # Set the X-ticks to be LONGITUDE 
    coord_pairs = wrf.to_np(qitot_cross_WSM6_2.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [round(pair.lon,2) for pair in wrf.to_np(coord_pairs)]
    axes[1,0].set_xticks(x_ticks[::20])
    axes[1,0].set_xticklabels(x_labels[::20],fontsize=12)   
    axes[1,1].set_xlabel("Longitude", fontsize=12)    

    # Set the X-ticks to be LATITUDE.
    coord_pairs = wrf.to_np(qitot_cross_P3_3mom_LF_1.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [round(pair.lat,2) for pair in wrf.to_np(coord_pairs)]
    axes[0,1].set_xticks(x_ticks[::10])
    axes[0,1].set_xticklabels(x_labels[::10],fontsize=12)                                 
    axes[0,1].set_xlabel("Latitude", fontsize=12)    

    # Set the X-ticks to be LONGITUDE 
    coord_pairs = wrf.to_np(qitot_cross_P3_3mom_LF_2.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [round(pair.lon,2) for pair in wrf.to_np(coord_pairs)]
    axes[1,1].set_xticks(x_ticks[::20])
    axes[1,1].set_xticklabels(x_labels[::20],fontsize=12)   
    axes[1,1].set_xlabel("Longitude", fontsize=12)  
        
    #------
    plt.suptitle( 'q_i (tot) Transects at ' + title + 'UTC')
    fig.savefig(savedir+title+'qi_tot_compare_contourmap_Transect.png', dpi=300,transparent=False)
    plt.show() 
    
    return 

#------------------------------------------------------------------------------
def plot_common_transect_P3micro(ncfile, x, y, title, savedir, fig, axes, ncol):
    
    # Need to interpolate first to common z_level grid 
    z_interp = 1e3*np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 
                 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20])
    
    Re       = 6.3781e6
    
    temp     = wrf.g_temp.get_tk(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL)
    z        = Re*geopo_p/(Re-geopo_p)

    # Create the start point and end point for the cross section for P3_3MOM_LF
    start_point = wrf.CoordPair(lat=y[0], lon=x[0])
    end_point   = wrf.CoordPair(lat=y[1], lon=x[1])

    [qr, qi, qc, qr_int, qi_int, qc_int, qir, qib, qil] = get_q_ints3(ncfile)
    
    QNCLOUD = wrf.getvar(ncfile, "QNCLOUD")     # 'cloud water Number concentration' ( kg-1)
    QNICE   = wrf.getvar(ncfile, "QNICE")       # 'Ice Number concentration' ( kg-1)
    QNRAIN  = wrf.getvar(ncfile, "QNRAIN")      # 'Rain Number concentration' ( kg-1)

    QZI  = wrf.getvar(ncfile, "QZI")            # Sq root of sixth moment ice mixing ratio times number mixing ratio'  ??? 

    # Compute the vertical cross-section interpolation.  Also, include the lat/lon points along the cross-section.
    qr_cross = wrf.vertcross(qr, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True)
    qc_cross = wrf.vertcross(qc, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True)
    qi_cross = wrf.vertcross(qi, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True)
    qir_cross = wrf.vertcross(qir, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True)
    qib_cross = wrf.vertcross(qib, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True)
    qil_cross = wrf.vertcross(qil, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True)
    #
    ncloud_cross = wrf.vertcross(QNCLOUD, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    nice_cross   = wrf.vertcross(QNICE, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    nrain_cross  = wrf.vertcross(QNRAIN, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    
    contour_x = np.arange(0, qi_cross.shape[-1], 1)
    contour_y = wrf.to_np(qi_cross.coords["vertical"])

    vmin = np.floor(np.log10(0.00002))-1; vmax = np.ceil(np.log10(0.01))+1; 
    lev_exp_qrain = np.logspace( vmin, vmax, num=int(vmax-vmin+1))
    vmin = np.floor(np.log10(0.000002))-1; vmax = np.ceil(np.log10(0.01))+1; 
    lev_exp_qcloud = np.logspace(vmin, vmax, num=int(vmax-vmin+1))
    vmin = np.floor(np.log10(0.00002))-1; vmax = np.ceil(np.log10(0.005))+1; 
    lev_exp_qice = np.logspace(vmin, vmax, num=int(vmax-vmin+1))
    vmin = np.floor(np.log10(0.000002))-1; vmax = np.ceil(np.log10(0.001))+1;
    lev_exp_qil = np.logspace(vmin, vmax, num=int(vmax-vmin+1))
        
    lev_exp_rain  = np.logspace(-3, 5, num=5-(-3)+1)
    lev_exp_cloud = np.logspace(-3, 10, num=10-(-3)+1)
    lev_exp_ice   = np.logspace(-3, 10, num=10-(-3)+1)

    contours_map_0 = axes[0,ncol].contourf(contour_x, contour_y, wrf.to_np(qr_cross), lev_exp_qrain, locator=ticker.LogLocator(), cmap=get_cmap("viridis")) #, vmin=0.0002, vmax=0.01)
    contours_map_1 = axes[1,ncol].contourf(contour_x, contour_y, wrf.to_np(qc_cross), lev_exp_qcloud, locator=ticker.LogLocator(),  cmap=get_cmap("viridis")) #, vmin=0.00002, vmax=0.05)
    contours_map_2 = axes[2,ncol].contourf(contour_x, contour_y, wrf.to_np(qi_cross), lev_exp_qice,  locator=ticker.LogLocator(), cmap=get_cmap("viridis")) #, vmin=0.0002, vmax=0.005)
    contours_map_3 = axes[3,ncol].contourf(contour_x, contour_y, wrf.to_np(qir_cross), lev_exp_qice, locator=ticker.LogLocator(), cmap=get_cmap("viridis")) #, vmin=0.0002, vmax=0.005)
    contours_map_4 = axes[4,ncol].contourf(contour_x, contour_y, wrf.to_np(qib_cross), lev_exp_qil, locator=ticker.LogLocator(), cmap=get_cmap("viridis"))
    contours_map_5 = axes[5,ncol].contourf(contour_x, contour_y, wrf.to_np(qil_cross), lev_exp_qil, locator=ticker.LogLocator(), cmap=get_cmap("viridis"))
    contours_map_6 = axes[6,ncol].contourf(contour_x, contour_y, wrf.to_np(nrain_cross),  lev_exp_rain, locator=ticker.LogLocator(), cmap=get_cmap("viridis"))
    contours_map_7 = axes[7,ncol].contourf(contour_x, contour_y, wrf.to_np(ncloud_cross), lev_exp_cloud, locator=ticker.LogLocator(), cmap=get_cmap("viridis"))
    contours_map_8 = axes[8,ncol].contourf(contour_x, contour_y, wrf.to_np(nice_cross),  lev_exp_ice, locator=ticker.LogLocator(), cmap=get_cmap("viridis"))
    
    axes[0,ncol].set_title('qr (kg/m3)')
    axes[1,ncol].set_title('qc (kg/m3)')
    axes[2,ncol].set_title('qi (kg/m3)')
    axes[3,ncol].set_title('qir (kg/m3)')
    axes[4,ncol].set_title('qib (kg/m3)')
    axes[5,ncol].set_title('qil (kg/m3)')
    axes[6,ncol].set_title('qnrain (1/kg)')
    axes[7,ncol].set_title('qncloud (1/kg)')
    axes[8,ncol].set_title('qnice (1/kg)')

    cbar =plt.colorbar(contours_map_0, ax=axes[0,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_1, ax=axes[1,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_2, ax=axes[2,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_3, ax=axes[3,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_4, ax=axes[4,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_5, ax=axes[5,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_6, ax=axes[6,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_7, ax=axes[7,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_8, ax=axes[8,ncol]); cbar.cmap.set_under('white')

    #--------------------------------------------------------------------------     # now add 0C level    
    [output_array_T, output_array_Z] = find_00(temp, 273, z)
    HGT = wrf.getvar(ncfile, "HGT", timeidx=-1)
    HGT.data = output_array_Z
    isoline_hgt = wrf.interpline(HGT, wrfin=ncfile, start_point=start_point, end_point=end_point)
    xs = np.arange(0, qi_cross.shape[-1], 1)
    for ii in range(9):
        axes[ii, ncol].plot(xs, wrf.to_np(isoline_hgt), '-r', linewidth=1.2)
        axes[ii, ncol].set_ylabel("Height (m)", fontsize=12)

    #--------------------------------------------------------------------------
    # Set the X-ticks to be LATITUDE.
    coord_pairs = wrf.to_np(qi_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [round(pair.lat,2) for pair in wrf.to_np(coord_pairs)]
    for ii in range(9):
        axes[ii,ncol].set_xticks(x_ticks[::10])
        axes[ii,ncol].set_xticklabels(x_labels[::10],fontsize=12)                                 
        axes[ii,ncol].set_xlabel("Latitude", fontsize=12)    

    # Set the X-ticks to be LONGITUDE 
    coord_pairs = wrf.to_np(qi_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [round(pair.lon,2) for pair in wrf.to_np(coord_pairs)]
    for ii in range(9):
        axes[ii,ncol].set_xticks(x_ticks[::20])
        axes[ii,ncol].set_xticklabels(x_labels[::20],fontsize=12)   
        axes[ii,ncol].set_xlabel("Longitude", fontsize=12)  
        
    return 

#------------------------------------------------------------------------------
def mu_solveequation(mu, right_eq):
  """Defines the equation to solve."""
  return ( (6+mu)*(5+mu)*(4+mu) / ((3+mu)*(2+mu)*(1+mu)) ) - right_eq



 
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def plot_common_transect_level_MAP(title, savedir):

    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    
    #---------------------------------------------------------
    # WSM6 files
    folderWSM6    = os.path.join('/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_WSM6_v4.5.2/', 'wrfout_d02*00:00')
    file_list_WSM6    = sorted(glob.glob(folderWSM6))
    timeNr_WSM6   = 92
    times          = [filename[timeNr_WSM6:timeNr_WSM6+16] for filename in file_list_WSM6]

    #---------------------------------------------------------
    # P3 files
    folderP3_1    = os.path.join('/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_P3_3MOM_LF_v4.5.2/', 'wrfout_d02*00:00')
    file_listP3_1     = sorted(glob.glob(folderP3_1))
    timeNr_P3_1   = 98
    times_p3          = [filename[timeNr_P3_1:timeNr_P3_1+16] for filename in file_listP3_1]
    
    #---------------------------------------------------------
    Re       = 6.3781e6
    
    ilevels_plot = np.array([0,1,2,3,4,5,6,7,8,9,10,15,20,25,30])

    for item, filename in enumerate(file_list_WSM6):
        if int(times[item][11:13]) >16:
            print(filename)
            print(item)
        
            itime = times[item]
        
            wrf_file = Dataset(filename,'r')
            #
            temp     = wrf.g_temp.get_tk(wrf_file)
            pressure = wrf.g_pressure.get_pressure(wrf_file)
            geopo_p  = wrf.g_geoht.get_height(wrf_file) # geopotential height as Mean Sea Level (MSL
            z_level  = Re*geopo_p/(Re-geopo_p)
            lat      =  wrf.getvar( wrf_file,"lat") 
            lon      =  wrf.getvar( wrf_file,"lon")
    
            qr = mixr2massconc( np.squeeze(wrf_file.variables["QRAIN"][0,:,:,:]  ), pressure, temp )        
            qi = mixr2massconc( np.squeeze(wrf_file.variables["QICE"][0,:,:,:]   ), pressure, temp )        
            qc = mixr2massconc( np.squeeze(wrf_file.variables["QCLOUD"][0,:,:,:] ), pressure, temp )       
            qs = mixr2massconc( np.squeeze(wrf_file.variables["QSNOW"][0,:,:,:] ), pressure, temp )       
            qg = mixr2massconc( np.squeeze(wrf_file.variables["QGRAUP"][0,:,:,:] ), pressure, temp )       
            
            qr = qr.where(qr >= 1e-5, np.nan)
            qs = qs.where(qs >= 1e-5, np.nan)
            qg = qg.where(qg >= 1e-5, np.nan)
            qi = qi.where(qi >= 1e-5, np.nan)
            
            #qr[qr<1E-5] = np.nan
            #qs[qs<1E-5] = np.nan
            #qg[qg<1E-5] = np.nan
            #qi[qi<1E-5] = np.nan   
            #qc[qc<1E-5] = np.nan
    
            # Similarly get the P3_1 info
            wrf_file_P3 = Dataset(file_listP3_1[item],'r')
            print(file_listP3_1[item])
           
            temp_p3     = wrf.g_temp.get_tk(wrf_file_P3)
            pressure_p3 = wrf.g_pressure.get_pressure(wrf_file_P3)
            geopo_p_p3  = wrf.g_geoht.get_height(wrf_file_P3) # geopotential height as Mean Sea Level (MSL
            z_level_p3  = Re*geopo_p_p3/(Re-geopo_p_p3)
            lat_p3      =  wrf.getvar( wrf_file_P3,"lat") 
            lon_p3      =  wrf.getvar( wrf_file_P3,"lon")
    
            qr_p3 = mixr2massconc( np.squeeze(wrf_file_P3.variables["QRAIN"][0,:,:,:]  ), pressure_p3, temp_p3 )        
            qi_p3 = mixr2massconc( np.squeeze(wrf_file_P3.variables["QICE"][0,:,:,:]   ), pressure_p3, temp_p3 )        
            qc_p3 = mixr2massconc( np.squeeze(wrf_file_P3.variables["QCLOUD"][0,:,:,:] ), pressure_p3, temp_p3 )       

            qr_p3 = qr_p3.where(qr_p3 >= 1e-5, np.nan)
            qi_p3 = qi_p3.where(qi_p3 >= 1e-5, np.nan)
            
            for iz in ilevels_plot:
                    
                fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True,figsize=[12,12])
                pcm0 = axes[0,0].pcolormesh(lon, lat, qr[iz,:,:], vmin=1e-5, vmax=0.001)
                pcm1 = axes[1,0].pcolormesh(lon, lat, qi[iz,:,:] + qs[iz,:,:] + qg[iz,:,:], vmin=1e-5, vmax=0.001)       
                
                cbar = plt.colorbar(pcm0, ax=axes[0,0], shrink=1, label='qr [kg/m$^2$]')
                cbar = plt.colorbar(pcm1, ax=axes[1,0], shrink=1, label='qi+qs+qg [kg/m$^2$]')

                cbar.cmap.set_under('white')
                
                for i in range(2):
                    axes[i,0].grid()
                    axes[i,0].set_xlim([-65.5,-62])
                    axes[i,0].set_ylim([-33.5,-31.3])
                    
                    axes[0,0].set_title('WSM6 q_rain')                                                  
                    axes[1,0].set_title('WSM6 q_i(tot)')  

                pcm2 = axes[0,1].pcolormesh(lon_p3, lat_p3, qr_p3[iz,:,:], vmin=1e-5, vmax=0.001)
                pcm3 = axes[1,1].pcolormesh(lon_p3, lat_p3, qi_p3[iz,:,:],vmin=1e-5, vmax=0.001)        
                
                cbar = plt.colorbar(pcm2, ax=axes[0,1], shrink=1, label='qr [kg/m$^2$]')
                cbar = plt.colorbar(pcm3, ax=axes[1,1], shrink=1, label='qi [kg/m$^2$]')

                cbar.cmap.set_under('white')
                
                for i in range(2):
                    axes[i,1].grid()
                    axes[i,1].set_xlim([-65.5,-62])
                    axes[i,1].set_ylim([-33.5,-31.3])
                    
                    axes[0,1].set_title('P3 3MOM LF q_rain')                                                  
                    axes[1,1].set_title('P3 3MOM LF q_i')  

                #plt.suptitle( 'Level qx (z='+ iz +') ' + title + 'UTC')
                plt.suptitle(title+ ' (level:'+ str(iz)+' '+itime +'UTC)')

                # RMA1 
                [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)
                for iax in range(2): 
                    axes[iax,0].plot(lon_radius, lat_radius, 'k', linewidth=0.8, label='RMA1 (120km)')
                    axes[iax,1].plot(lon_radius, lat_radius, 'k', linewidth=0.8, label='RMA1 (120km)')

                axes[0,0].legend(fontsize=10, loc='lower left')
                 
                # #
                # xx = np.vstack([x_wsm6[[0,2]],x_wsm6[[1,3]]])
                # yy = np.vstack([y_wsm6[[0,2]],y_wsm6[[1,3]]])
                # axes[0].plot(xx, yy, '-or' , linewidth=1.2)
                # #
                # xx = np.vstack([x_P3_3mom_LF[[0,2]],x_P3_3mom_LF[[1,3]]])
                # yy = np.vstack([y_P3_3mom_LF[[0,2]],y_P3_3mom_LF[[1,3]]])
                # axes[1].plot(xx, yy, '-or' , linewidth=1.2)
                #
                fig.savefig(savedir+title+'qi_tot_compare_contourmap.png', dpi=300,transparent=False)
                fig.savefig(savedir+title+itime+'_iz'+str(iz)+'.png', dpi=300,transparent=False)
                plt.close()

    return

#------------------------------------------------------------------------------
def plot_layer_maps(folderWRFOUT, timeNr, title, savedir):

    # Domain-average evolution of precip rate
    file_list = sorted(glob.glob(folderWRFOUT))
    times          = [filename[timeNr:timeNr+16] for filename in file_list]
    
    Re       = 6.3781e6
    
    ilevels_plot = np.array([0,1,2,3,4,5,6,7,8,9,10,15,20,25,30])

    for item, filename in enumerate(file_list):
        if int(times[item][11:13]) >16:
            print(filename)
        
            itime = times[item]
        
            wrf_file = Dataset(filename,'r')
            #
            temp     = wrf.g_temp.get_tk(wrf_file)
            pressure = wrf.g_pressure.get_pressure(wrf_file)
            geopo_p  = wrf.g_geoht.get_height(wrf_file) # geopotential height as Mean Sea Level (MSL
            z_level  = Re*geopo_p/(Re-geopo_p)
            lat      =  wrf.getvar( wrf_file,"lat") 
            lon      =  wrf.getvar( wrf_file,"lon")
    
            qr = mixr2massconc( np.squeeze(wrf_file.variables["QRAIN"][0,:,:,:]  ), pressure, temp )        
            qi = mixr2massconc( np.squeeze(wrf_file.variables["QICE"][0,:,:,:]   ), pressure, temp )        
            qc = mixr2massconc( np.squeeze(wrf_file.variables["QCLOUD"][0,:,:,:] ), pressure, temp )       
            qs = mixr2massconc( np.squeeze(wrf_file.variables["QSNOW"][0,:,:,:] ), pressure, temp )       
            qg = mixr2massconc( np.squeeze(wrf_file.variables["QGRAUP"][0,:,:,:] ), pressure, temp )       
    
                # get max value for the common variables:
            vminnn= 0.001
            vmaxxx= np.nanmax((qr, qi + qs + qg))
                
            for iz in ilevels_plot:
                    
                fig, axes = plt.subplots(nrows=1, ncols=6, constrained_layout=True,figsize=[12,2])
                pcm0 = axes[0].pcolormesh(lon, lat, qr[iz,:,:])
                pcm1 = axes[1].pcolormesh(lon, lat, qc[iz,:,:])
                pcm2 = axes[2].pcolormesh(lon, lat, qi[iz,:,:] + qs[iz,:,:] + qg[iz,:,:])        
                pcm3 = axes[3].pcolormesh(lon, lat, qi[iz,:,:])
                pcm4 = axes[4].pcolormesh(lon, lat, qs[iz,:,:])
                pcm5 = axes[5].pcolormesh(lon, lat, qg[iz,:,:])   
                
                cbar = plt.colorbar(pcm0, ax=axes[0], shrink=1, label='qx [kg/m$^2$]')
                cbar = plt.colorbar(pcm1, ax=axes[1], shrink=1, label='qx [kg/m$^2$]')
                cbar = plt.colorbar(pcm2, ax=axes[2], shrink=1, label='qx [kg/m$^2$]')
                cbar = plt.colorbar(pcm3, ax=axes[3], shrink=1, label='qx [kg/m$^2$]')
                cbar = plt.colorbar(pcm4, ax=axes[4], shrink=1, label='qx [kg/m$^2$]')
                cbar = plt.colorbar(pcm5, ax=axes[5], shrink=1, label='qx [kg/m$^2$]')

                cbar.cmap.set_under('white')
                
                for i in range(6):
                    axes[i].grid()
                    axes[i].set_xlim([-65.5,-62])
                    axes[i].set_ylim([-33.5,-31.3])
                    
                    axes[0].set_title('q_rain')                                                  
                    axes[1].set_title('q_cld')  
                    axes[2].set_title('q_i(tot)') 
                    axes[3].set_title('q_i') 
                    axes[4].set_title('q_s') 
                    axes[5].set_title('q_g') 
    
                plt.suptitle(title+ ' (level:'+ str(iz)+' '+itime +'UTC)')
                
                fig.savefig(savedir+title+itime+'_iz'+str(iz)+'.png', dpi=300,transparent=False)
                plt.close() 
            
        
    return

#------------------------------------------------------------------------------
def plot_WSM6_all(strfile_WSM6, rma1_str, dow6_str, dow7_str, csapr2_str, savedir):
    
    #==========================    WSM6    =======================================
    #=============================================================================  
    # Make WRF qx plots + 45dbz contours for all radars
    for ii, istr in enumerate(strfile_WSM6): 
        
        wrf_file = Dataset(istr,'r')
        itime = istr[-8:-3] 
        plot_qxs(title='WSM6 ('+itime+') + RMA1', itime=0, ncfile=wrf_file, rfile=rma1_str[ii], savedir=savedir, radar_name='RMA1')
        plot_qxs(title='WSM6 ('+itime+') + DOW6', itime=0, ncfile=wrf_file, rfile=dow6_str[ii], savedir=savedir, radar_name='DOW6')
        plot_qxs(title='WSM6 ('+itime+') + DOW7', itime=0, ncfile=wrf_file, rfile=dow7_str[ii], savedir=savedir, radar_name='DOW7')
        plot_qxs(title='WSM6 ('+itime+') + CSAPR2', itime=0, ncfile=wrf_file, rfile=csapr2_str[ii], savedir=savedir, radar_name='CSAPR2')
        
        plot_qxs_qitot(title='WSM6 ('+itime+') + RMA1', itime=0, ncfile=wrf_file, rfile=rma1_str[ii], savedir=savedir, radar_name='RMA1')
        plot_qxs_qitot(title='WSM6 ('+itime+') + DOW6', itime=0, ncfile=wrf_file, rfile=dow6_str[ii], savedir=savedir, radar_name='DOW6')
        plot_qxs_qitot(title='WSM6 ('+itime+') + DOW7', itime=0, ncfile=wrf_file, rfile=dow7_str[ii], savedir=savedir, radar_name='DOW7')
        plot_qxs_qitot(title='WSM6 ('+itime+') + CSAPR2', itime=0, ncfile=wrf_file, rfile=csapr2_str[ii], savedir=savedir, radar_name='CSAPR2')


    # Make WRF plots with actual radar range per hour and no Zh contours
    istr = strfile_WSM6[4]
    wrf_file = Dataset(istr,'r')
    itime = istr[-8:-3] 
    plot_qxs_onlyWRF(title='WSM6 ('+itime+')', itime=0, ncfile=wrf_file, savedir=savedir, 
                          DOW6file=dow6_str[4], DOW7file=dow7_str[4], CSAPR2file=csapr2_str[4])
    
    istr = strfile_WSM6[1]
    wrf_file = Dataset(istr,'r')
    itime = istr[-8:-3] 
    plot_qxs_onlyWRF(title='WSM6 ('+itime+')', itime=0, ncfile=wrf_file, savedir=savedir, 
                          DOW6file=dow6_str[1], DOW7file=dow7_str[1], CSAPR2file=csapr2_str[1])
    
    return


def plot_radar_obs(rma1_str, dow6_str, dow7_str, csapr2_str, savedir, 
             cswr_dir, csapr2_dir, output_gif_path, rma1_dir, tempsavedir, strfile_WSM6):
    

    #========================= OBSERVATIONS ======================================
    #=============================================================================  
    # Make a plot of all radar data for itime=1 (17:40) y itime=4 (2020)
    general_plot(rma1_str[1], dow6_str[1], dow7_str[1], csapr2_str[1], savedir+'/', '1740UTC')
    general_plot(rma1_str[4], dow6_str[4], dow7_str[4], csapr2_str[4], savedir+'/', '2020UTC')

    #Make gifs for qrain, q_snow, q_grau, Zh(RMA1), Zh(DOW7), ZH(DOW6), ZH(CSAPR2) on common timeline 17:00 - 20:30
    # Generate every 10-minute interval in minutes within the range (minutes from midnight)
    start_time = 17*60  #17:00
    end_time = 20*60+30 #20:30
    all_times_minutes = np.arange(start_time, end_time + 10, 10)

    # DOW6
    make_radar_gif(cswr_dir, 'cfrad.20181110*DOW6*s02*.nc', 86, all_times_minutes, 'DOW6', output_gif_path+'/', tempsavedir)
    # DOW7
    make_radar_gif(cswr_dir, 'cfrad.20181110*DOW7*s02*.nc', 86, all_times_minutes, 'DOW7', output_gif_path+'/', tempsavedir)
    # CSAPR2 
    make_radar_gif(csapr2_dir, 'corcsapr2cfrppiqcM1.b1.20181110.*.nc', 101, all_times_minutes, 'CSAPR2', output_gif_path+'/', tempsavedir) 
    # RMA1
    make_radar_gif(rma1_dir, 'cfrad.20181110*01.nc', 82, all_times_minutes, 'RMA1', output_gif_path+'/', tempsavedir) 

    
    
    return

#------------------------------------------------------------------------------
# RUN MAIN
def main():
    
    server = 'yakaira'
    
    if 'yakaira' in server: 
        # Leer wrfout y variables de interes: control-ndg
        #folder_WSM6  = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_WSM6_v4.5.2_all/'
        # re-run with correct namelist.input 
        folder_WSM6 = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFOUT/WSM6/'
        folder_P3_3MOM_LF = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_P3_3MOM_LF_v4.5.2/' 
        folder_P3_3MOM_noLF = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_P3_3MOM_noLF_v4.5.2/' 
        savedir = '/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/Plots/WRF_WSM6/'
        # Within this folder, define the name of a sub-folder according to date
        savedir       = os.path.join(savedir,datetime.datetime.now().strftime('%d%m%Y'))
        # If the latter sub-folder does not exist, create it.
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        tempsavedir = '/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/Plots/tmp/'
        output_gif_path = '/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/Plots/gifs/'

    
        savedirP3 = '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/WRF_P3/'
        csapr2_dir = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/CSAPR2/'
        rma1_dir = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RMA1/'
        cswr_dir = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/CSWRdata/'
        save_dir_compare = '/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/Plots/'
        save_dir_compare       = os.path.join(save_dir_compare,datetime.datetime.now().strftime('%d%m%Y'))
        # If the latter sub-folder does not exist, create it.
        if not os.path.exists(save_dir_compare):
            os.makedirs(save_dir_compare)        
        if not os.path.exists(save_dir_compare+'/WRFcompare'):
            os.makedirs(save_dir_compare+'/WRFcompare')        

        
        
    else:
        # Leer wrfout y variables de interes: control-ndg
        folder_WSM6  = '/home/vito.galligani/disk1/derek/Work/HAILCASE_10112018_datos/WRFout_WSM6_v4.5.2/'
        savedir = '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/WRF_WSM6/'
        csapr2_dir = '/home/vito.galligani/disk1/derek/Work/HAILCASE_10112018_datos/CSAPR2/'
        rma1_dir = '/home/vito.galligani/disk1/derek/Work/HAILCASE_10112018_datos/RMA1/'
        tempsavedir = '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/tmp/'
    
    # wrf files
    strfile_WSM6 = [folder_WSM6+'wrfout_d02_2018-11-10_17:00:00', 
               folder_WSM6+'wrfout_d02_2018-11-10_17:40:00',
               folder_WSM6+'wrfout_d02_2018-11-10_19:40:00',
               folder_WSM6+'wrfout_d02_2018-11-10_20:00:00',
               folder_WSM6+'wrfout_d02_2018-11-10_20:20:00']
    
    strfile_P3_3MOM_LF = [folder_P3_3MOM_LF+'wrfout_d02_2018-11-10_17:00:00', 
               folder_P3_3MOM_LF+'wrfout_d02_2018-11-10_17:40:00',
               folder_P3_3MOM_LF+'wrfout_d02_2018-11-10_19:40:00',
               folder_P3_3MOM_LF+'wrfout_d02_2018-11-10_20:00:00',
               folder_P3_3MOM_LF+'wrfout_d02_2018-11-10_20:20:00']
    
    # Datos radar [CSAPR-2]
    csapr2_str = [csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.170003.nc',
                  csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.174503.nc',
                  csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.194503.nc',
                  csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.200003.nc',
                  csapr2_dir+'corcsapr2cfrppiqcM1.b1.20181110.201746.nc']
    
    # Datos de radar [RMA1]
    rma1_str_02 = [rma1_dir+'cfrad.20181110_170419.0000_to_20181110_170542.0000_RMA1_0301_02.nc',
                  rma1_dir+'cfrad.20181110_174533.0000_to_20181110_174657.0000_RMA1_0301_02.nc',
                  rma1_dir+'cfrad.20181110_194101.0000_to_20181110_194221.0000_RMA1_0301_02.nc',
                  rma1_dir+'cfrad.20181110_200546.0000_to_20181110_200709.0000_RMA1_0301_02.nc',
                  rma1_dir+'cfrad.20181110_202213.0000_to_20181110_202332.0000_RMA1_0301_02.nc']
    
    rma1_str = [rma1_dir+'cfrad.20181110_170542.0000_to_20181110_171223.0000_RMA1_0301_01.nc',
                  rma1_dir+'cfrad.20181110_173842.0000_to_20181110_174523.0000_RMA1_0301_01.nc',
                  rma1_dir+'cfrad.20181110_194224.0000_to_20181110_194913.0000_RMA1_0301_01.nc',
                  rma1_dir+'cfrad.20181110_200709.0000_to_20181110_201358.0000_RMA1_0301_01.nc',
                  rma1_dir+'cfrad.20181110_202339.0000_to_20181110_203020.0000_RMA1_0301_01.nc']
    
    # Datos de radar [DOW6]
    dow6_str = [cswr_dir+'cfrad.20181110_170011_DOW6high_v215_s01_el0.49_SUR.nc',
                  cswr_dir+'cfrad.20181110_174011_DOW6high_v227_s01_el0.49_SUR.nc',
                  cswr_dir+'cfrad.20181110_194512_DOW6high_v260_s02_el0.63_SUR.nc',    
                  cswr_dir+'cfrad.20181110_200014_DOW6high_v266_s02_el0.92_SUR.nc',
                  cswr_dir+'cfrad.20181110_202015_DOW6high_v274_s02_el0.90_SUR.nc']                                                         # FALTA
    
    # Datos de radar [DOW7]
    dow7_str = [cswr_dir+'cfrad.20181110_170011_DOW7high_v216_s01_el0.95_SUR.nc',
                cswr_dir+'cfrad.20181110_172011_DOW7high_v222_s01_el0.97_SUR.nc',
                cswr_dir+'cfrad.20181110_194054_DOW7high_v256_s01_el0.70_SUR.nc',  #LE LLUEVE ENCIMA! 
                cswr_dir+'cfrad.20181110_200023_DOW7high_v264_s03_el1.14_SUR.nc',  #LE LLUEVE ENCIMA! 
                cswr_dir+'cfrad.20181110_202012_DOW7high_v271_s02_el0.99_SUR.nc']  #JUSTO PASÓ
    

    folder_WSM6 = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_WSM6_v4.5.2_domain2/'
    elev = 3 # RMA1 elev 3 == 1.78 for .01.nc
    VP.plot_MAXRADAR_WRF_evolution(folder_WSM6, 111, 'test', save_dir_compare+'/', '', '', elev)
    VP.plot_MAXRADAR_WRF_evolution_radar(folder_WSM6, 111, 'test', save_dir_compare+'/', '', '', elev)

    breakpoint()

    


    # Plot subplots of WRFOUT qx with radar rings frorm all radars at different selected times:
    #==================
    plot_WSM6_all(strfile_WSM6, rma1_str, dow6_str, dow7_str, csapr2_str, savedir)
        
    # Plot radar observations (and gifs)
    #==================
    plot_radar_obs(rma1_str, dow6_str, dow7_str, csapr2_str, savedir, 
                 cswr_dir, csapr2_dir, output_gif_path, rma1_dir, tempsavedir, strfile_WSM6)
            
    # re-do for new P3
    # ========================== Compare wsm6 y p3 ================================
    for ii, istr in enumerate(strfile_P3_3MOM_LF): 
        wrf_file_wsm6 = Dataset(strfile_WSM6[ii],'r')
        wrf_file_p3 = Dataset(strfile_P3_3MOM_LF[ii],'r')
        itime = istr[-8:-3] 
        # Plot maps of precip accumulation
        plot_accums(wrf_file_wsm6, wrf_file_p3, itime, save_dir_compare+'/WRFcompare/', dow6_str[1], dow7_str[1], csapr2_str[1])
        # Plot domain average vertical profiles of qx. 
        plot_domain_qxs(wrf_file_wsm6, wrf_file_p3, itime, save_dir_compare+'/WRFcompare/')

    # ========================== Compare wsm6 y p3 PRECIP ACCUM ===================
    # Domain-average evolution of precip rate
    domainlons = [-65.5,-62]
    domainlats = [-33.5,-31.3] 
    #folder_WSM6    = os.path.join('/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_WSM6_v4.5.2/', 'wrfout_d02*')
    folder_WSM6    = os.path.join('/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFOUT/WSM6/', 'wrfout_d02*')
    [accumPrecip_average_wsm6, times_wsm6] = get_domainAccum(folder_WSM6, domainlons, domainlats,85) 

    # re-do for new P3
    folder_p3    = os.path.join('/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_P3_3MOM_LF_v4.5.2/', 'wrfout_d02*')
    [accumPrecip_average_p3, times_p3] = get_domainAccum(folder_p3, domainlons, domainlats, 98) 
    
    # ========================== Compare wsm6 y p3 common transect ?  ===================
    # look first at 20:20 Utc
    wrf_file_WSM6       = Dataset(strfile_WSM6[4],'r')
    wrf_file_p3_3mom_lf = Dataset(strfile_P3_3MOM_LF[4],'r')
    itime = strfile_WSM6[4][-8:-3] 
    #
    # check transects (plot two transects x= x11, x12, x21, x22)
    x_wsm6=np.array([-63.45, -63.45, -63.75, -62.0])
    y_wsm6=np.array([-32.0, -33.0, -32.00, -33.35])
    x_P3_3mom_LF=np.array([-63.6, -63.6, -63.8, -62.0])
    y_P3_3mom_LF=np.array([-32.0, -33.0, -32.00, -33.5])   
    
    fig = plt.figure(figsize=(8,8)) 
    formatter = matplotlib.dates.DateFormatter('%H:%M')
    ax = fig.add_subplot(111)
    # #- WSM6
    x_time = [datetime.datetime.strptime(d,'%Y-%m-%d_%H:%M') for d in times_wsm6]
    plt.plot(x_time, accumPrecip_average_wsm6, 'darkred', linewidth=1.2, label='WSM6')
    # #- P3
    x_time_p3 = [datetime.datetime.strptime(d,'%Y-%m-%d_%H:%M') for d in times_p3]
    plt.plot(x_time_p3, accumPrecip_average_p3, 'darkblue', linewidth=1.2, label='P3 3MOM LF')
    # #
    plt.legend()
    plt.xlabel('Time (UTC)')
    plt.ylabel('Accumulated Precip (mm)')
    ax.xaxis.set_major_formatter(formatter)
    plt.setp(ax.get_xticklabels(), rotation=15)
    plt.show()
    fig.savefig(save_dir_compare+'compareDomain_PrecipAccum.png', dpi=300,transparent=False)

    VP.plot_common_transect_MAP(wrf_file_WSM6, wrf_file_p3_3mom_lf, itime, dow6_str[4], dow7_str[4], csapr2_str[4], save_dir_compare+'/WRFcompare/', 'QICE', x_wsm6, y_wsm6, x_P3_3mom_LF, y_P3_3mom_LF)
    VP.plot_common_transect_MAP(wrf_file_WSM6, wrf_file_p3_3mom_lf, itime, dow6_str[4], dow7_str[4], csapr2_str[4], save_dir_compare+'/WRFcompare/', 'MAXRADAR_WRF', x_wsm6, y_wsm6, x_P3_3mom_LF, y_P3_3mom_LF)
    plot_common_transect(wrf_file_WSM6, wrf_file_p3_3mom_lf, itime, save_dir_compare+'/WRFcompare/')

    # Plot all microphysics
    fig, axes = plt.subplots(nrows=9, ncols=2, constrained_layout=True,figsize=[16,30])
    #----  the first transect
    x    = [-63.6, -63.6]
    y    = [-32.0, -33.0]
    ncol = 0
    plot_common_transect_P3micro(wrf_file_p3_3mom_lf, x, y, itime, save_dir_compare, fig, axes, ncol)
    # add the second transect 
    x    = [-63.8, -62.0]
    y    = [-32.00, -33.5]
    ncol = 1
    plot_common_transect_P3micro(wrf_file_p3_3mom_lf, x, y, itime, save_dir_compare, fig, axes, ncol)
    fig.savefig(save_dir_compare+'P3micro_compare_contourmap_Transect.png', dpi=300,transparent=False)
    plt.show() 

    fig, axes = plt.subplots(nrows=8, ncols=2, constrained_layout=True,figsize=[16,30])
    # #----  the first transect
    x    = [-63.6, -63.6]
    y    = [-32.0, -33.0]
    ncol = 0
    VP.plot_common_transect_P3micro_ICE_derived(wrf_file_p3_3mom_lf, x, y, itime, save_dir_compare, fig, axes, ncol)
    # # add the second transect 
    x    = [-63.8, -62.0]
    y    = [-32.00, -33.5]
    ncol = 1
    VP.plot_common_transect_P3micro_ICE_derived(wrf_file_p3_3mom_lf, x, y, itime, save_dir_compare, fig, axes, ncol)
    fig.savefig(save_dir_compare+'P3micro_compare_contourmap_Transect_ICE_CALC.png', dpi=300,transparent=False)
    plt.show() 

    #- WSM6
    fig, axes = plt.subplots(nrows=6, ncols=2, constrained_layout=True,figsize=[16,30])
    #----  the first transect
    x    = [-63.6, -63.6]
    y    = [-32.0, -33.0]
    ncol = 0
    VP.plot_common_transect_WSM6micro_ICE_derived(wrf_file_WSM6, x, y, itime, save_dir_compare, fig, axes, ncol)
    # add the second transect 
    x    = [-63.8, -62.0]
    y    = [-32.00, -33.5]
    ncol = 1
    VP.plot_common_transect_WSM6micro_ICE_derived(wrf_file_WSM6, x, y, itime, save_dir_compare, fig, axes, ncol)
    fig.savefig(save_dir_compare+'WSM6micro_compare_contourmap_Transect_ICE_CALC.png', dpi=300,transparent=False)
    plt.show() 

    warnings.filterwarnings("ignore")
    folder_WSM6    = os.path.join('/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_WSM6_v4.5.2/', 'wrfout_d02*00:00')
    plot_common_transect_level_MAP('Layer_qx', save_dir_compare+'level_plots/' ) 

    
    return

main()







#------------------------------------------------------------------------------ FOR WSM6
# [all_vars, all_descriptions] = print_wrf_variables(wrf_file)
#
# 'RAINC' ===> 'ACCUMULATED TOTAL CUMULUS PRECIPITATION', [mm]
# 'RAINSH' ===>  'ACCUMULATED SHALLOW CUMULUS PRECIPITATION', [mm]
# 'RAINNC' ===>   'ACCUMULATED TOTAL GRID SCALE PRECIPITATION', [mm]
# 'SNOWNC' ===>   'ACCUMULATED TOTAL GRID SCALE SNOW AND ICE', [mm]
# 'GRAUPELNC' ===>   'ACCUMULATED TOTAL GRID SCALE GRAUPEL', [mm]
# 'HAILNC' ===>   'ACCUMULATED TOTAL GRID SCALE HAIL', [mm]
# 'CLDFRA' ===>   'CLOUD FRACTION',
#
#------------------------------------------------------------------------------ FOR P3
#[all_vars, all_descriptions] = print_wrf_variables(wrf_file)
#
# 'QVAPOR' ===> 'Water vapor mixing ratio'
# 'QCLOUD' ===> 'Cloud water mixing ratio',
# 'QRAIN' ===> 'Rain water mixing ratio',
# 'QICE' ===> 'Ice mixing ratio',
# 'QNCLOUD' ===> 'cloud water Number concentration',
# 'QNICE' ===> 'Ice Number concentration',
# 'QNRAIN' ===> 'Rain Number concentration',
# 'QIR' ===> 'Rime ice mass-1 mixing ratio'
# 'QIB' ===> 'Rime ice volume-1 mixing ratio'
# 'QZI' ===> 'Sq root of sixth moment ice mixing ratio times number mixing ratio'
# 'QIL' ===>  'Liquid water mixing ratio accumulated on ice',
#  
# 'RAINC' ===> 'ACCUMULATED TOTAL CUMULUS PRECIPITATION', [mm]
# 'RAINSH' ===>  'ACCUMULATED SHALLOW CUMULUS PRECIPITATION', [mm]
# 'RAINNC' ===>   'ACCUMULATED TOTAL GRID SCALE PRECIPITATION', [mm]
# 'SNOWNC' ===>   'ACCUMULATED TOTAL GRID SCALE SNOW AND ICE', [mm]
# 'GRAUPELNC' ===>   'ACCUMULATED TOTAL GRID SCALE GRAUPEL', [mm]
# 'HAILNC' ===>   'ACCUMULATED TOTAL GRID SCALE HAIL', [mm]
# 'CLDFRA' ===> 'CLOUD FRACTION', 
# 'REFL_10CM' ===>    'Radar reflectivity (lamda = 10 cm)',

#
# NOTES: RAINC, RAINSH, SNOWNC, GRAUPELNC, HAILNC: 0
# CLDFRA IS POR NIVEL
# 
#------------------------------------------------------------------------------ 



# ========================== P3_3MOM_LF =======================================
# =============================================================================  
# Make WRF qx plots + 45dbz contours for all radars
# for ii, istr in enumerate(strfile_P3_3MOM_LF): 
#      wrf_file = Dataset(istr,'r')
#      itime = istr[-8:-3] 
#      plot_qxs_P3(title='P3 3MOM LF ('+itime+') + RMA1', itime=0, ncfile=wrf_file, rfile=rma1_str[ii], savedir=savedirP3, radar_name='RMA1')
#      plot_qxs_P3(title='P3 3MOM LF ('+itime+') + DOW6', itime=0, ncfile=wrf_file, rfile=dow6_str[ii], savedir=savedirP3, radar_name='DOW6')
#      plot_qxs_P3(title='P3 3MOM LF ('+itime+') + DOW7', itime=0, ncfile=wrf_file, rfile=dow7_str[ii], savedir=savedirP3, radar_name='DOW7')
#      plot_qxs_P3(title='P3 3MOM LF ('+itime+') + CSAPR2', itime=0, ncfile=wrf_file, rfile=csapr2_str[ii], savedir=savedirP3, radar_name='CSAPR2')

# # Make WRF plots with actual radar range per hour and no Zh contours
# istr = strfile_P3_3MOM_LF[4]
# wrf_file = Dataset(istr,'r')
# itime = istr[-8:-3] 
# plot_qxs_onlyWRF_P3(title='P3 3MOM LF ('+itime+')', itime=0, ncfile=wrf_file, savedir=savedirP3, 
#                       DOW6file=dow6_str[4], DOW7file=dow7_str[4], CSAPR2file=csapr2_str[4])
 
# istr = strfile_P3_3MOM_LF[1]
# wrf_file = Dataset(istr,'r')
# itime = istr[-8:-3] 
# plot_qxs_onlyWRF_P3(title='P3 3MOM LF ('+itime+')', itime=0, ncfile=wrf_file, savedir=savedirP3, 
#                       DOW6file=dow6_str[1], DOW7file=dow7_str[1], CSAPR2file=csapr2_str[1])
# =============================================================================














