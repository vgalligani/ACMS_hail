#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------

-----------------------------------------------------------------
"""

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
import os
import contextlib
import warnings
import logging

plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.ERROR)  # or CRITICAL

@contextlib.contextmanager
def suppress_stdout():
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout



# =============================================================================
def main_basic():

    server = 'yakaira'
    main_Process_cs_andWRF('MHS', '20:30', 6, server)
    print('Finished clear sky and wrf grids')    
    
    return

# =============================================================================
def main_Process_cs_andWRF(instrument, HHtime, mp_version, server): 

    plotpath, folders = config_folders(server)

    if (mp_version == 6):
      mp_physics = 'WRF-WSM6'
    elif mp_version == 6.1:
        mp_physics = 'WRFWSM6_YSU12hr'
                
    # Select server and folder locations
    #--------------------------------------------------------------------------
    if 'yakaira' in server: 
        upfolder     = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/'
        sys.path.insert(1,'/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/src')
        processedFolder = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RTTOVout/Processed/'+mp_physics
      
    elif 'cnrm' in server:
        upfolder    = '/home/galliganiv/'   
        sys.path.insert(1,'/home/galliganiv/ACMS_hail/src')
        processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'+mp_physics

    # Select instrument configurations
    #--------------------------------------------------------------------------
    if 'MHS' in instrument: 
        nchan = 5
    elif 'AMSR' in instrument:
        nchan = 14
        
    #--------------------------------------------------------------------------
    # server dependent files    
    from package_functions import pressure2height
    import package_functions as funs
    
    #--------------------------------------------
    # Filter monotonic pressure and h2o > 0.1E-10
    skipProfs = filter_pixels_monotonic(mp_version, HHtime, server)

    if mp_version == 6:

        mp_physics = 'WRF-WSM6'
        ncfile     = upfolder+'WRFOUT/WSM6_domain3_NoahMP/wrfout_d02_2018-11-10_'+HHtime+':00'
        ncdata     = Dataset(ncfile,'r') 
        [qr, qs, qi, qc, qg, qr_int, qs_int, qi_int, qc_int, qg_int] = funs.get_q_ints6(ncdata)
        int_titles = ['qr','qc','qi','qs','qg']

    elif mp_version == 6.1:
        
        mp_physics = 'WRFWSM6_YSU12hr'
        ncfolder  = upfolder+'WRFOUT_1domain/WSM6_1domaintest_12hrs_YSU/'    
        ncfile    = ncfolder+'wrfout_d01_2018-11-10_'+HHtime+':00'  
        
        ncdata     = Dataset(ncfile,'r') 
        [qr, qs, qi, qc, qg, qr_int, qs_int, qi_int, qc_int, qg_int] = funs.get_q_ints6(ncdata)
        int_titles = ['qr','qc','qi','qs','qg']
        
        # RTTOVout folders and file. 
        outfolder         = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RTTOVout/WRFWSM6_YSU12hr/rttov_sieron/'
        WSM6_file         = np.genfromtxt(outfolder+'output_cs_tb_'+instrument)
   
    
    skipProfs = set(skipProfs)  # Make lookup faster

    ij_list = []
    counter = 0
    for i in range(390, 550):
        for j in range(100, 410):
            counter += 1
            if counter not in skipProfs:
                ij_list.append((i, j))  # This profile was used
             
                
    A    = read_wrf(ncfile)        
    tb_csrttov = np.full((5, A['XLONG'].shape[0], A['XLONG'].shape[1]   ), np.nan)  # or whatever shape you need
    for idx, (i, j) in enumerate(ij_list):
        tb_csrttov[:, i, j] = WSM6_file[idx, :]  # transpose from [49499, 5] to [5, i, j]
    
                            
    lats = A['XLAT']
    lons = A['XLONG']
    q    = A['h20']
    qinttot = np.nansum( [qi_int.data, qc_int.data, qs_int.data,  
                        qr_int.data, qg_int.data], axis=0 ) 

    # para make zalt efficient:
    # Extract your subdomain
    p_sub = A['pressure'][:, 390:550, 100:410]  # shape (nlev, 160, 310)
    T_sub = A['T'][:, 390:550, 100:410]         # shape (nlev, 160, 310)

    # Reshape to apply along each column (n_profiles = 160*310)
    nlev, ni, nj = p_sub.shape
    n_profiles = ni * nj

    p_flat = p_sub.reshape(nlev, n_profiles)
    T_flat = T_sub.reshape(nlev, n_profiles)

    # Now apply pressure2height to each profile (axis 0 is vertical)
    zalt_flat = np.stack([
        pressure2height(p_flat[:, k], T_flat[:, k])
        for k in range(n_profiles)
    ], axis=1)  # shape (nlev, n_profiles)

    # Reshape back to (nlev, ni, nj)
    zalt_sub = zalt_flat.reshape(nlev, ni, nj)
    
    # Initialize full array with NaNs and insder zalt_sub
    zalt = np.full((nlev, A['XLONG'].shape[0], A['XLONG'].shape[1]   ), np.nan)
    zalt[:, 390:550, 100:410] = zalt_sub[:,:,:]


    #----- SIMPLE PLOTS: make plots with integrated qx forzen y rain, y todos los canales
    # y stats. comparar diferencia clear sky con obs.# improves with atlas? 

    ds = T2P.MHS_cs_sims(lons, lats, q, tb_csrttov, plotpath, server)
    ds.to_netcdf(processedFolder+'/'+'rttov_processed_clearsky.nc', 'w')

    # Ok. so far we have ds and das with cs and as Tbs (average, std, gaussian, min and max)
    # and i would like the same for the WRF outputs. 
    #
    # I need to conduct the same type of 'averaging' over the different qints. 
    #-------------------------------------------------------------------------
    # 1) Point q_int on lat/lon grid
    #-------------------------------------------------------------------------
    # 2) Interpolated to MHS lat/lon cut domain
    #-------------------------------------------------------------------------
    
    WRFinterp_qtot_int  = T2P.regrid2D_toMHSgrid(lats, lons, ds, qinttot.reshape(-1))
    WRFinterp_landMask  = T2P.regrid2D_toMHSgrid(lats, lons, ds, A['LANDMASK'].reshape(-1))
    WRFinterp_qr_int    = T2P.regrid2D_toMHSgrid(lats, lons, ds, qr_int.data.reshape(-1))
    WRFinterp_qs_int    = T2P.regrid2D_toMHSgrid(lats, lons, ds, qs_int.data.reshape(-1))
    WRFinterp_qg_int    = T2P.regrid2D_toMHSgrid(lats, lons, ds, qg_int.data.reshape(-1))
    WRFinterp_qi_int    = T2P.regrid2D_toMHSgrid(lats, lons, ds, qi_int.data.reshape(-1))
    WRFinterp_qc_int    = T2P.regrid2D_toMHSgrid(lats, lons, ds, qc_int.data.reshape(-1))
    
    # 3) Average over footprint 
    #-------------------------------------------------------------------------
    WRF_intqtot_fpt_mean, WRF_intqtot_fpt_min, WRF_intqtot_fpt_max, WRF_intqtot_fpt_std  = T2P.average_over_footprint(lats, lons, ds, qinttot)
    WRF_intqr_fpt_mean, WRF_intqr_fpt_min, WRF_intqr_fpt_max, WRF_intqr_fpt_std          = T2P.average_over_footprint(lats, lons, ds, qr_int.data)
    WRF_intqs_fpt_mean,WRF_intqs_fpt_min, WRF_intqs_fpt_max,  WRF_intqs_fpt_std    = T2P.average_over_footprint(lats, lons, ds, qs_int.data)
    WRF_intqg_fpt_mean, WRF_intqg_fpt_min, WRF_intqg_fpt_max, WRF_intqg_fpt_std    = T2P.average_over_footprint(lats, lons, ds, qg_int.data)
    WRF_intqi_fpt_mean, WRF_intqi_fpt_min, WRF_intqi_fpt_max, WRF_intqi_fpt_std    = T2P.average_over_footprint(lats, lons, ds, qi_int.data)
    WRF_intqc_fpt_mean, WRF_intqc_fpt_min, WRF_intqc_fpt_max, WRF_intqc_fpt_std    = T2P.average_over_footprint(lats, lons, ds, qc_int.data)

    # plt.pcolormesh(ds['MHS_lon'], ds['MHS_lat'], WRF_intqtot_fpt_mean); plt.colorbar()
    # 3) GAUSSIAN ?
    #-------------------------------------------------------------------------                                              
    WRF_intqtot_gaussian = T2P.overGaussian(lats, lons, ds, qinttot)
    WRF_intqr_gaussian   = T2P.overGaussian(lats, lons, ds, qr_int.data)
    WRF_intqs_gaussian   = T2P.overGaussian(lats, lons, ds, qs_int.data)
    WRF_intqg_gaussian   = T2P.overGaussian(lats, lons, ds, qg_int.data)
    WRF_intqi_gaussian   = T2P.overGaussian(lats, lons, ds, qi_int.data)
    WRF_intqc_gaussian   = T2P.overGaussian(lats, lons, ds, qc_int.data)


    # plt.pcolormesh(ds['MHS_lon'], ds['MHS_lat'], var); plt.colorbar()
    # Organizar aux data en un dataframe:
    dwrf =  xr.Dataset({
         "wrf_lat":    (["lat","lon"], lats), 
         "wrf_lon":    (["lat","lon"], lons),
         "wrf_zalt":   (["var","lat","lon"], zalt),
         "WRF_qr":     (["var","lat","lon"],  qr.data),
         "WRF_qs":     (["var","lat","lon"],  qs.data),   
         "WRF_qg":     (["var","lat","lon"],  qg.data),   
         "WRF_qi":     (["var","lat","lon"],  qi.data),   
         "WRF_qc":     (["var","lat","lon"],  qc.data),   
         "WRF_intqr":  (["lat","lon"],  qr_int.data),
         "WRF_intqs":  (["lat","lon"],  qs_int.data),
         "WRF_intqg":  (["lat","lon"],  qg_int.data),
         "WRF_intqc":  (["lat","lon"],  qc_int.data),
         "WRF_intqi":  (["lat","lon"],  qi_int.data),
         "WRF_intTot":  (["lat","lon"], qinttot),                 
         "WRF_LandMask":   (["lat","lon"],  A['LANDMASK']),
         #
         "MHSinterp_intqr":  (["lat2","lon2"],  WRFinterp_qr_int),
         "MHSinterp_intqs":  (["lat2","lon2"],  WRFinterp_qs_int),
         "MHSinterp_intqg":  (["lat2","lon2"],  WRFinterp_qg_int),
         "MHSinterp_intqc":  (["lat2","lon2"],  WRFinterp_qc_int),
         "MHSinterp_intqi":  (["lat2","lon2"],  WRFinterp_qi_int),
         "MHSinterp_intTot":  (["lat2","lon2"],  WRFinterp_qtot_int),
         "MHSinterp_LandMask":  (["lat2","lon2"],  WRFinterp_landMask),
         #
         "MHSfootprintmean_intqr":  (["lat2","lon2"],  WRF_intqr_fpt_mean),
         "MHSfootprintmean_intqs":  (["lat2","lon2"],  WRF_intqs_fpt_mean),
         "MHSfootprintmean_intqg":  (["lat2","lon2"],  WRF_intqg_fpt_mean),
         "MHSfootprintmean_intqi":  (["lat2","lon2"],  WRF_intqi_fpt_mean),
         "MHSfootprintmean_intqc":  (["lat2","lon2"],  WRF_intqc_fpt_mean),
         "MHSfootprintmean_intTot":  (["lat2","lon2"], WRF_intqtot_fpt_mean), 
         #
         "MHSfootprintmin_intqr":  (["lat2","lon2"],  WRF_intqr_fpt_min),
         "MHSfootprintmin_intqs":  (["lat2","lon2"],  WRF_intqs_fpt_min),
         "MHSfootprintmin_intqg":  (["lat2","lon2"],  WRF_intqg_fpt_min),
         "MHSfootprintmin_intqi":  (["lat2","lon2"],  WRF_intqi_fpt_min),
         "MHSfootprintmin_intqc":  (["lat2","lon2"],  WRF_intqc_fpt_min),
         "MHSfootprintmin_intTot":  (["lat2","lon2"], WRF_intqtot_fpt_min),              
         #
         "MHSfootprintmax_intqr":  (["lat2","lon2"],  WRF_intqr_fpt_max),
         "MHSfootprintmax_intqs":  (["lat2","lon2"],  WRF_intqs_fpt_max),
         "MHSfootprintmax_intqg":  (["lat2","lon2"],  WRF_intqg_fpt_max),
         "MHSfootprintmax_intqi":  (["lat2","lon2"],  WRF_intqi_fpt_max),
         "MHSfootprintmax_intqc":  (["lat2","lon2"],  WRF_intqc_fpt_max),
         "MHSfootprintmax_intTot":  (["lat2","lon2"], WRF_intqtot_fpt_max), 
         #             
         "MHSfootprintstd_intqr":  (["lat2","lon2"],  WRF_intqr_fpt_std),
         "MHSfootprintstd_intqs":  (["lat2","lon2"],  WRF_intqs_fpt_std),
         "MHSfootprintstd_intqg":  (["lat2","lon2"],  WRF_intqg_fpt_std),
         "MHSfootprintstd_intqi":  (["lat2","lon2"],  WRF_intqi_fpt_std),
         "MHSfootprintstd_intqc":  (["lat2","lon2"],  WRF_intqc_fpt_std),
         "MHSfootprintstd_intTot":  (["lat2","lon2"], WRF_intqtot_fpt_std),       
         #
         "MHSGaussian_intqr":  (["lat2","lon2"],  WRF_intqr_gaussian),
         "MHSGaussian_intqs":  (["lat2","lon2"],  WRF_intqs_gaussian),
         "MHSGaussian_intqg":  (["lat2","lon2"],  WRF_intqg_gaussian),
         "MHSGaussian_intqi":  (["lat2","lon2"],  WRF_intqi_gaussian),
         "MHSGaussian_intqc":  (["lat2","lon2"],  WRF_intqc_gaussian),
         "MHSGaussian_intTot":  (["lat2","lon2"],  WRF_intqtot_gaussian) 
         })
    #dwrf.to_netcdf(processedFolder+'/'+'wrfdata_processed.nc', 'w')
    
        
    return

# =============================================================================


# =============================================================================
def main_Process(instrument, HHtime, mp_version, server, allskyfile, expname, rttov_approach): 

    plotpath, folders = config_folders(server)
    
    # Select server and folder locations
    #--------------------------------------------------------------------------
    if 'yakaira' in server: 
        upfolder     = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/'
        sys.path.insert(1,'/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/src')
        processedFolder = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RTTOVout/Processed/'
      
    elif 'cnrm' in server:
        upfolder    = '/home/galliganiv/'   
        sys.path.insert(1,'/home/galliganiv/ACMS_hail/src')
        processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'
        
    #--------------------------------------------------------------------------
    # server dependent files    
    from package_functions import pressure2height
    import package_functions as funs
        
    #--------------------------------------------------------------------------
    if mp_version == 6:

        mp_physics = 'WRF-WSM6'
        ncfile     = upfolder+'WRFOUT/WSM6_domain3_NoahMP/wrfout_d02_2018-11-10_'+HHtime+':00'
        ncdata     = Dataset(ncfile,'r') 
        [qr, qs, qi, qc, qg, qr_int, qs_int, qi_int, qc_int, qg_int] = funs.get_q_ints6(ncdata)
        int_titles = ['qr','qc','qi','qs','qg']
        skipProfs = filter_pixels_monotonic(mp_version, '20:30', server)   

    elif mp_version == 6.1:
        
        mp_physics = 'WRFWSM6_YSU12hr'
        ncfolder  = upfolder+'WRFOUT_1domain/WSM6_1domaintest_12hrs_YSU/'    
        ncfile    = ncfolder+'wrfout_d01_2018-11-10_'+HHtime+':00'  
        
        ncdata     = Dataset(ncfile,'r') 
        [qr, qs, qi, qc, qg, qr_int, qs_int, qi_int, qc_int, qg_int] = funs.get_q_ints6(ncdata)
        int_titles = ['qr','qc','qi','qs','qg']
        
        # RTTOVout folders and file. 
        outfolder = f'/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RTTOVout/WRFWSM6_YSU12hr/rttov_{rttov_approach}/tests/'
        skipProfs = filter_pixels_monotonic(mp_version, '20:00', server)   
    
    #--------------------------------------------------------------------------
    processedFolder = processedFolder+mp_physics
    
    #--------------------------------------------------------------------------
    A           = read_wrf(ncfile)
    allskyfile='output_tb_allsky__'+allskyfile
    print('Running for: '+ allskyfile)
    WSM6_allsky = np.genfromtxt(outfolder+allskyfile) 
                                    
    # Select instrument configurations
    #--------------------------------------------------------------------------
    if 'MHS' in instrument: 
        nchan = 5
    elif 'AMSR' in instrument:
        nchan = 14
           
    #------        
    ij_list = []
    counter = 0
    for i in range(390, 550):
        for j in range(100, 410):
            counter += 1
            if counter not in skipProfs:
                ij_list.append((i, j))  # This profile was used
             
    tb_asrttov = np.full((5, A['XLONG'].shape[0], A['XLONG'].shape[1]), np.nan)  # or whatever shape you need
    for idx, (i, j) in enumerate(ij_list):
        tb_asrttov[:, i, j] = WSM6_allsky[idx, :]  # transpose from [49499, 5] to [5, i, j]
        
    lats = A['XLAT']
    lons = A['XLONG']
        
    # Pre-process like this if MHS                
    if 'MHS' in instrument: 
        das1 = T2P.MHS_as_sims(lons, lats, tb_asrttov[:,:,:], plotpath, server, expname)
        das1.to_netcdf(processedFolder+'/'+f'rttov_processed_allsky_{rttov_approach}_'+expname+'.nc', 'w')
        das1.close()
        gc.collect()
        del tb_asrttov

    return





    
# =============================================================================
# THIS SUBSECTION IS TO MAKE_PROFS:
# nohup python main_process_efficient.py > output.log 2>&1 &
if __name__ == "__main__":
    
# =============================================================================
#     if len(sys.argv) < 3:
#         print("Usage w. grausp: python your_script.py <exp> <isnow> ")
#         sys.exit(1)
#     
#     exp     = sys.argv[1]    
#     isnow   = int(sys.argv[2])     # or float() if needed
#     igrau   = int(sys.argv[3])     # or float() if needed
#     
# 
#     # Only suppress output from noisy libraries or processes
#     with suppress_stdout():
#         # # rttov iwcrescaling: 
#         # allskyfile = f'sliu{isnow}_gliu{igrau}'+exp            
#         # expname    = 'allsky_liusnow'+str(isnow)+'_liugrau'+str(igrau)+exp
#         # main_Process('MHS', '20:30', 6.1, 'yakaira', allskyfile, expname, 'iwcrescaling' )
# 
#         # rttov eqmass: 
#         allskyfile = f'sliu{isnow}_gliu{igrau}'+exp            
#         expname    = 'allsky_liusnow'+str(isnow)+'_liugrau'+str(igrau)+exp
#         main_Process('MHS', '20:30', 6.1, 'yakaira', allskyfile, expname, 'eqmass')
#         print('Finished allsky for isnow: ', str(isnow) )  
#         
#         # rttov sieron
#         #allskyfile = f'sliu{isnow}_gliu{igrau}dummy_2snowprofile_cirsgs_6nhydro_'+exp
#         #expname    = 'allsky_liusnow'+str(isnow)+'liugrau'+str(igrau)+exp        
#         #main_Process('MHS', '20:30', 6.1, 'yakaira', allskyfile, expnamem 'sieron' )
# =============================================================================
    
    # Process clear sky and WRF variables
    #main_Process_cs_andWRF('MHS', '20:00', 6.1, 'yakaira')
    
    exp = ''
    isnow1 = 8
    isnow2 = 9
    
    for igrau in [5]: #5
        allskyfile = f'sliu{isnow1}_gliu{igrau}_aggregates{isnow2}'+'2snowprofile_temp2'+exp            
        expname    = 'allsky__liusnow'+str(isnow1)+'gliu'+str(igrau)+'_aggregates'+str(isnow2)+'_2snowprofile_temp2'
        main_Process('MHS', '20:30', 6.1, 'yakaira', allskyfile, expname, 'sieron')
        plt.close()
        print('Finished allsky for igrau: '+ str(igrau))  
        
    
    
    
        
        

    
    
    
    
    
    
    
    
