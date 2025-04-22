#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------
@purpose : Run rttov for the 2018-11-10 hail case (WRF WSM6 and P3)
@author  : Vito Galligani
@email   : vito.galligani@gmail.com
@condaenv: conda activate 
-----------------------------------------------------------------
@main    : Makes the necessary RTTOV v14 input profiles for all-sky simulations 
          
           remember use for reference: https://github.com/lkugler/RTTOV-WRF/blob/master/rttov_wrf.py
           and also to open remote session:
               nohup python -m spyder_kernels.console -matplotlib='inline' -f=./remotemachine.json &4
               sshfs vito.galligani@yakaira.cima.fcen.uba.ar:/home/vito.galligani/Work remote_project
               
              
@PROFILES: Some profiles are ignored based on pressure (not monotonic and outside rttov valid limits)

@TODO: not yet consistent with P3-scheme
-----------------------------------------------------------------
"""


#------------------------------------------------------------------------------
#                   NOTES 
# 1) RTTOV IS RUN OUTSIDE. CURRENTLY THE CODE BELOW IS HARDCODED FOR AMSRE? BUT SHOULD BE FOR MHS AT 20.00
# 2) pensar en agregar checks de invalid profiles:
#      Input water vapour profile exceeds allowed minimum:
# For example:
# 2025/03/14  13:17:29  Limit               =  0.1000E-10
# 2025/03/14  13:17:29  Upper level p (hPa) =    116.9933
# 2025/03/14  13:17:29  Value               =  0.1642E-28
#
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# to ensure that pyrttov is importable

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
import h5py

plt.matplotlib.rc('font', family='serif', size = 18)
plt.rcParams['xtick.labelsize']=18
plt.rcParams['ytick.labelsize']=18  


#----------------------------------------------------------------------------------------------------
# footprint operator     
def plot_footprint(server, d_cs): 
    
    lons = d_cs['wrf_lon'].data
    lats = d_cs['wrf_lat'].data
    tb0  = d_cs['rttov_as'].data
    
    plotpath, folders = config_folders(server)

    # Some basic info for colormaps  
    if 'yakaira' in server:
        prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')    
        fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
        mhs_noaa19_dir  = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/PMW/GPM/RS/V07/PMW/1C-MHS-NOAA19/2018/11/10/'

    elif 'cnrm' in server:
        prov = np.genfromtxt("/home/galliganiv/ACMS_hail/src/provincias.txt", delimiter='')    
        fn = '/home/galliganiv/ACMS_hail/src/etopo1_bedrock.nc'   
        mhs_noaa19_dir  = '/home/galliganiv/Work/HAILCASE_10112018/PMW/GPM/RS/V07/PMW/1C-MHS-NOAA19/2018/11/10/'

    
    mhs_noaa19_file = mhs_noaa19_dir+'1C.NOAA19.MHS.XCAL2021-V.20181110-S201428-E215628.050287.V07A.HDF5'
    extent = [ -70, -50, -40, -20]

    #-----
    Nrlimitswath   = 30
    with h5py.File(mhs_noaa19_file, "r") as f:
        Tc = f['/S1/Tc'][:,:,:]  # shape: [channel, y, x]
        lat = f['/S1/Latitude'][:,:]  # shape: [y, x]
        lon = f['/S1/Longitude'][:,:]  # shape: [y, x]
        
    # Crop by extent
    mask = (
        (lon >= extent[0]) & (lon <= extent[1]) &
        (lat >= extent[2]) & (lat <= extent[3])
        )
    rows, cols = np.where(mask)
    # If any points are found, crop arrays to the bounding box
    if rows.size > 0 and cols.size > 0:
        row_min, row_max = rows.min(), rows.max()
        col_min, col_max = cols.min(), cols.max()
        
        # Crop all arrays accordingly
        Tc_crop = Tc[row_min:row_max+1, :, :]
        lat_crop = lat[row_min:row_max+1, :]
        lon_crop = lon[row_min:row_max+1, :]
    else:
        raise ValueError("No data points found within the specified extent.")
    
    lat_cropt = lat_crop.T
    lon_cropt = lon_crop.T
    Tc_cropt  = Tc_crop.T
    # Create labeled DataArray
    da_subset = xr.DataArray(
        Tc_cropt,
        dims=("nchan","x", "y"),
        coords={
            "lat": (("x", "y"), lat_cropt),
            "lon": (("x", "y"), lon_cropt),
            },
        name="Tc",
        attrs={
            "units": "K",
            "description": "Brightness temperature at 89 GHz",
            "source": "NOAA MHS"
            })    
    
    #-----------------------------------------------------------------------------
    domain_mhs_obs = da_subset.data
    domain_mhs_obs = domain_mhs_obs[:,:Nrlimitswath,:].copy()

    pointInterpTB = T2P.get_pointdataInterp(lats, lons, tb0, da_subset, Nrlimitswath) 
    T2P.get_footprintVals_poster(lats, lons, tb0, da_subset, Nrlimitswath, domain_mhs_obs, plotpath, server, 'as_test') 

    return

    
#--------------------------------------------
def main():
    
    server     = 'cnrm'
    instrument = 'MHS'
    HHtime     = '20:30'
    mp_version = 6
    mp_physics = 'WRF-WSM6'

    processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'+mp_physics
    outfile   = 'output_tb_'+instrument

    expname     = 'rttov_processed_allsky_eqMass_rsg_s10g3.nc'
    d_liuliu    = xr.open_dataset(processedFolder+'/'+outfile+expname)
    plot_footprint(server, d_liuliu)


if __name__ == "__main__":
    main()

# =============================================================================
#         #----- SIMPLE PLOTS of qxs for model cs 
#         # Plot qxs maps 
#         T2P.plot_test_qxints(dwrf, ds, lats, lons)
#         # Plot Tbs maps 
#         data_footprint = np.ma.masked_greater_equal(dwrf['MHSfootprintmean_intTot'], 1)
#         data_regrid    = np.ma.masked_greater_equal(dwrf['MHS_intTot'], 1)
#         cmaps = T2P.GMI_colormap() 
#         prov  = np.genfromtxt("/home/galliganiv/ACMS_hail/src/provincias.txt", delimiter='')      
#         VAR   = np.ma.array(ds['rttov_cs_footprintmean'].data[0,:,:], mask=data_footprint.mask) 
#         fig = plt.figure(figsize=[8,8])
#         plt.pcolormesh(ds['MHS_lon'], ds['MHS_lat'], VAR, cmap=cmaps['turbo_r'], 
#                                  shading='auto') 
#         plt.xlim([-65.5,-62])   #[-68,-62]); 
#         plt.ylim([-34.5,-31])   #[-36,-31])
#         plt.plot(prov[:,0],prov[:,1],color='w')
#         # Diff in clear sky map        
#         VARdiff = np.ma.array( ds['MHs_domain_obs'].data[0,:,:]-VAR, mask=data_footprint.mask) 
#         fig = plt.figure(figsize=[8,8])
#         plt.pcolormesh(ds['MHS_lon'], ds['MHS_lat'], VARdiff, cmap=cmaps['turbo_r'], 
#                                  shading='auto') 
#         plt.xlim([-65.5,-62])   #[-68,-62]); 
#         plt.ylim([-34.5,-31])   #[-36,-31])
#         plt.plot(prov[:,0],prov[:,1],color='w'); plt.colorbar()
# =============================================================================





