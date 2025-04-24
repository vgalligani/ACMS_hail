#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------
@purpose : Run rttov for the 2018-11-10 hail case (WRF WSM6 and P3)
@author  : Vito Galligani
@email   : vito.galligani@gmail.com
@condaenv: conda activate 
-----------------------------------------------------------------
@main    : 1) Makes the necessary RTTOV v14 input profiles for all-sky simulations 
           2) Creates .nc files with valid i,j grid for Clear sky and hydro experiments
             [note that for liu experiments I use bash script to run in parallel all .nc creations]

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
# RTTOV IS RUN OUTSIDE. 
# CURRENTLY THE CODE BELOW IS HARDCODED FOR AMSRE? BUT SHOULD BE FOR MHS AT 20.00
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


plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def eqmass_exp(folders, outfoldereq, mp_physics, HHtime, instrument, experiment, nchan): #where exp='__eqMass_WSM6_rsg'
	
	# main output folder
    main_folder = folders['read_out_dir']+mp_physics+'/'
    subfolder = mp_physics+'_20181110_'+HHtime+'_'+instrument+'_atlas_satzen_input_'
        
    counter = 0
    for i in range(11):
        for  j in range(11):
            Expname   = experiment+'_sliu'+str(i)+'_gliu'+str(j)
            file_folder = main_folder + subfolder+Expname+'/'
            file_       = 'output_as_tb_'+instrument+Expname
            if (counter == 0):
                init = np.genfromtxt(file_folder+file_)
                tb = np.zeros(( 11, 11, init.shape[0], init.shape[1] )); tb[:]=np.nan
            tb[i,j,:,:] =  np.genfromtxt(file_folder+file_)
            counter=counter+1

        return tb       
 
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def eqmass_exp_one(folders, outfoldereq, mp_physics, HHtime, instrument, experiment, nchan, i, j): #where exp='__eqMass_WSM6_rsg'
	
	# main output folder
    main_folder = folders['read_out_dir']+mp_physics+'/'
    subfolder = mp_physics+'_20181110_'+HHtime+'_'+instrument+'_atlas_satzen_input_'
        
    Expname   = experiment+'_sliu'+str(i)+'_gliu'+str(j)
    file_folder = main_folder + subfolder+Expname+'/'
    file_       = 'output_as_tb_'+instrument+Expname
    tb =  np.genfromtxt(file_folder+file_)

    return tb       

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def eqmass_exp_one_iwc(folders, outfoldereq, mp_physics, HHtime, instrument, experiment, nchan, i, j, iwc_exp): #where exp='__eqMass_WSM6_rsg'

        # main output folder
    main_folder = folders['read_out_dir']+mp_physics+'/'
    subfolder = mp_physics+'_20181110_'+HHtime+'_'+instrument+'_atlas_satzen_input_'

    Expname   = experiment+'_sliu'+str(i)+'_gliu'+str(j)+'__half'+iwc_exp
    file_folder = main_folder + subfolder+Expname+'/'
    file_       = 'output_as_tb_'+instrument+experiment+'_sliu'+str(i)+'_gliu'+str(j)
    tb =  np.genfromtxt(file_folder+file_)

    return tb

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def eqmass_exp_grausp(folders, outfoldereq, mp_physics, HHtime, instrument, experiment, nchan, i): #where exp='__eqMass_WSM6_rsg'

    # main output folder
    main_folder = folders['read_out_dir']+mp_physics+'/'
    subfolder = mp_physics+'_20181110_'+HHtime+'_'+instrument+'_atlas_satzen_input_'

    Expname   = experiment+'_sliu'+str(i)+'grausp'
    file_folder = main_folder + subfolder+Expname+'/'
    file_       = 'output_as_tb_'+instrument+Expname
    tb =  np.genfromtxt(file_folder+file_)

    return tb


    
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main_makeProfs(instrument, HHtime, mp_version, server): 

    plotpath, folders = config_folders(server)
    
    # Select server and folder locations
    #--------------------------------------------------------------------------
    if 'yakaira' in server: 
        upfolder     = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/'
        sys.path.insert(1,'/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/src')
                
    elif 'cnrm' in server:
        upfolder    = '/home/galliganiv/'   
        sys.path.insert(1,'/home/galliganiv/ACMS_hail/src')
            
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
    skipProfs_monotonic = filter_pixels_monotonic(mp_version, HHtime, server)

    # Make profiles and run rttov on terminal 
    ipnd=1
    run_IFS_rttov14(ipnd, mp_version, HHtime, instrument, skipProfs_monotonic, server)

    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main_Process_cs_andWRF(instrument, HHtime, mp_version, server): 

    plotpath, folders = config_folders(server)

    if (mp_version == 6):
      mp_physics = 'WRF-WSM6'
    
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

    
    # Select server and folder locations
    flag_name = 'rttov14_'+instrument+'_'+mp_physics+'_2018-11-10_'+HHtime

    # RTTOVout folders
    outfolder         = mp_physics+'_20181110_'+HHtime+'_'+instrument +'_atlas_satzen_input/'
    WSM6_file         = np.genfromtxt(folders['read_out_dir']+mp_physics+'/'+outfolder+'output_cs_tb_'+instrument)
   
    #------        
    # Load all profiles
    A    = read_wrf(ncfile)
    toti = A['XLONG'].shape[0]
    totj = A['XLONG'].shape[1]       
    lats         = np.zeros(A['XLONG'].shape); lats[:]=np.nan
    lons         = np.zeros(A['XLONG'].shape); lons[:]=np.nan
    q            = np.zeros(A['h20'].shape);      q[:]=np.nan
    qinttot      = np.zeros(A['XLONG'].shape);   qinttot[:]=np.nan  
    rows, cols   = A['XLONG'].shape  # Original dimensions
    zalt         = np.zeros(A['pressure'].shape);      zalt[:]=np.nan

    #---------------------------------------------------------------------
    tb_csrttov   = np.zeros( (nchan,rows,cols) );   tb_csrttov[:]=np.nan 
    tb1          = np.zeros( (nchan,rows,cols) );   tb1[:]=np.nan 
     
    #--------
    # Read scatt rttov outputs
    #--- allskytesting (delfat test with official hydro_table. no psd nonsistency and
    # overlap_param == cloud_overlap_2col_weighted and per_hydro_frac == false)
    outfile_as_test   = 'output_as_tb_'+instrument+'test'        
    outfoldetest      = mp_physics+'_20181110_'+HHtime+'_'+instrument +'_atlas_satzen_input_test/'
    WSM6_file_as_test = np.genfromtxt(folders['read_out_dir']+mp_physics+'/'+outfoldetest+outfile_as_test)
    tb_asrttov_test   = np.zeros( (nchan,rows,cols) );   tb_asrttov_test[:]=np.nan 

    # Save some aux WRF data to dataframe
    qinttot = np.nansum( [qi_int.data, qc_int.data, qs_int.data,  
                        qr_int.data, qg_int.data], axis=0 )
            
    counter = 0
    rttov_counter = 0
    for i in range(A['XLONG'].shape[0]): 
        for j in range(A['XLONG'].shape[1]): 

            counter = counter+1
            
            lats[i,j] = A['XLAT'].data[i,j]
            lons[i,j] = A['XLONG'].data[i,j]
                
            if counter in skipProfs:
                print('Skipping this profile '+str(counter))
                q[:,i,j]    = np.nan
                qi_int.data[i,j] = np.nan
                qs_int.data[i,j] = np.nan
                qg_int.data[i,j] = np.nan
                qr_int.data[i,j] = np.nan
                qc_int.data[i,j] = np.nan
                qr.data[:,i,j]   = np.nan
                qs.data[:,i,j]   = np.nan
                qi.data[:,i,j]   = np.nan
                qc.data[:,i,j]   = np.nan
                qg.data[:,i,j]   = np.nan
                qinttot[i,j]     = np.nan
                zalt[:,i,j]      = np.nan
                
                #tb1[:,i,j] = np.nan
                tb_csrttov[:,i,j]       = np.nan
                tb_asrttov_test[:,i,j]  = np.nan
                #tb_asrttov_eqMass_rsg[:,:,:,i,j] = np.nan
                #tb_asrttov_rsg[:,:,:,i,j] = np.nan
                
            else:
                rttov_counter     = rttov_counter+1
                q[:,i,j]          = A['h20'].data[:,i,j]
                zalt[:,i,j]       = pressure2height(A['pressure'][:,i,j], A['T'][:,i,j])
                
                tb_asrttov_test[:,i,j] = WSM6_file_as_test[rttov_counter-1,:]
                tb_csrttov[:,i,j]      = WSM6_file[rttov_counter-1,:]
                
            
    #----- SIMPLE PLOTS: make plots with integrated qx forzen y rain, y todos los canales
    # y stats. comparar diferencia clear sky con obs.# improves with atlas? 
            
    if 'MHS' in instrument: 
        
        outfile           = 'output_tb_'+instrument

        # Plot simulation (WRF resolution) and Observations (real resolution)
        ds = T2P.MHS_cs_sims(lons, lats, q, tb_csrttov, plotpath, server)
        ds.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_clearsky.nc', 'w')

        das = T2P.MHS_as_sims(lons, lats, tb_asrttov_test, plotpath, server, '_test')
        das.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_test.nc', 'w')
        
    if 'AMSR' in instrument: 
        T2P.plot_simple_AMSR2_TEST(lons, lats, q, tb_csrttov, tb1, plotpath)
        T2P.plot_simple_AMSR2_comparison(lons, lats, q, tb_csrttov, tb1, plotpath, 'amsr-2')

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
    dwrf.to_netcdf(processedFolder+'/'+'wrfdata_processed.nc', 'w')
    
        
    return

        
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main_Process_Expliu(instrument, HHtime, mp_version, server, skipProfs, eqMass_do, isnow, igrau): 

    if (mp_version == 6):
      mp_physics = 'WRF-WSM6'

    plotpath, folders = config_folders(server)
    
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
    
    if mp_version == 6:
        mp_physics = 'WRF-WSM6'
        ncfile     = upfolder+'WRFOUT/WSM6_domain3_NoahMP/wrfout_d02_2018-11-10_'+HHtime+':00'
        ncdata     = Dataset(ncfile,'r') 
        [qr, qs, qi, qc, qg, qr_int, qs_int, qi_int, qc_int, qg_int] = funs.get_q_ints6(ncdata)
        int_titles = ['qr','qc','qi','qs','qg']

    flag_name = 'rttov14_'+instrument+'_'+mp_physics+'_2018-11-10_'+HHtime

    # RTTOVout folders
    outfoldereq       = mp_physics+'_20181110_'+HHtime+'_'+instrument +'_atlas_satzen_input__eqmassWSM6_'
   
    #------        
    # Load all profiles
    A    = read_wrf(ncfile)
    toti = A['XLONG'].shape[0]
    totj = A['XLONG'].shape[1]       
    lats         = np.zeros(A['XLONG'].shape); lats[:]=np.nan
    lons         = np.zeros(A['XLONG'].shape); lons[:]=np.nan
    rows, cols = A['XLONG'].shape
    shape_ = A['XLONG'].shape
 
    if (eqMass_do == 1): 
        	#--- eqmassWSM6_rsg_s10g2: equal mass PSD consistency with WRF and 
        # read for all liu - liu combinations w/ snow and grau
        # default cloud overlap settings as above (+renormalization)
        #exp_asrttov_eqmass_rsgliu  = eqmass_exp(folders, outfoldereq, mp_physics, HHtime, instrument, '_eqMass_WSM6_rsg', nchan)
        exp_asrttov_eqmass_rsgliu = eqmass_exp_one(folders, outfoldereq, mp_physics, HHtime, instrument, '_eqMass_WSM6_rsg', nchan, isnow, igrau)
    elif (eqMass_do == 0): 
        	#--- WSM6_rsg_s10g2: Dmax PSD consistency with WRF and 
        # read for all liu - liu combinations w/ snow and grau
        # default cloud overlap settings as above (re-normalization)
        #exp_asrttov_rsgliu  = eqmass_exp(folders, outfoldereq, mp_physics,HHtime, instrument, '_WSM6_rsg', nchan)
        exp_asrttov_rsgliu = eqmass_exp_one(folders, outfoldereq, mp_physics, HHtime, instrument, '_WSM6_rsg', nchan, isnow, igrau)
    elif (eqMass_do == 3):
        exp_asrttov_rsgliu_gsp = eqmass_exp_grausp(folders, outfoldereq, mp_physics, HHtime, instrument, '_WSM6_rsg', nchan, isnow)

    elif (eqMass_do == 4):
        exp_asrttov_eqmass_rsgliu_gsp = eqmass_exp_grausp(folders, outfoldereq, mp_physics, HHtime, instrument, '_eqMass_WSM6_rsg_', nchan, isnow)


    outfile = 'output_tb_'+instrument
    tb_asrttov_eqMass_rsg  = np.zeros( (nchan,rows,cols) );   tb_asrttov_eqMass_rsg[:]=np.nan 
    tb_asrttov_rsg         = np.zeros( (nchan,rows,cols) );   tb_asrttov_rsg[:]=np.nan 
    tb_asrttov_rsgliu_gsp  = np.zeros( (nchan,rows,cols) );   tb_asrttov_rsgliu_gsp[:]=np.nan
    tb_asrttov_eqMass_rsgliu_gsp  = np.zeros( (nchan,rows,cols) );   tb_asrttov_eqMass_rsgliu_gsp[:]=np.nan

    counter = 0
    rttov_counter = 0
    
    for idx, (i,j) in enumerate(np.ndindex(shape_)):
        lats[i,j] = A['XLAT'].data[i,j]
        lons[i,j] = A['XLONG'].data[i,j]
        
        if idx in skipProfs:
            tb_asrttov_eqMass_rsg[:,i,j] = np.nan
            tb_asrttov_rsg[:,i,j] = np.nan
        
        else:
            if (eqMass_do==1):
                tb_asrttov_eqMass_rsg[:,i,j] = exp_asrttov_eqmass_rsgliu[rttov_counter-1,:]
            elif (eqMass_do==2): 
                tb_asrttov_rsg[:,i,j] = exp_asrttov_rsgliu[rttov_counter-1,:]
            elif (eqMass_do==3):	    
                tb_asrttov_rsgliu_gsp[:,i,j] = exp_asrttov_rsgliu_gsp[rttov_counter-1,:]
            elif (eqMass_do==4):
                tb_asrttov_eqMass_rsgliu_gsp[:,i,j] = exp_asrttov_eqmass_rsgliu_gsp[rttov_counter-1,:]

            rttov_counter=rttov_counter+1
                    
    # Pre-process like this if MHS                
    if 'MHS' in instrument: 
        

        if (eqMass_do == 0): 
            das1 = T2P.MHS_as_sims(lons, lats, tb_asrttov_rsg[:,:,:], plotpath, server, '_rsg_s'+str(isnow)+'g'+str(igrau))
            das1.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s'+str(isnow)+'g'+str(igrau)+'.nc', 'w')
            das1.close()
            gc.collect()
            del tb_asrttov_rsg
            
        elif (eqMass_do == 1): 
            das1 = T2P.MHS_as_sims(lons, lats, tb_asrttov_eqMass_rsg[:,:,:], plotpath, server, '_rsg_s'+str(isnow)+'g'+str(igrau))
            das1.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_eqMass_rsg_s'+str(isnow)+'g'+str(igrau)+'.nc', 'w')
            das1.close()
            gc.collect()
            del tb_asrttov_eqMass_rsg
        elif (eqMass_do == 3):
            das1 = T2P.MHS_as_sims(lons, lats, tb_asrttov_rsgliu_gsp[:,:,:], plotpath, server, '_rsg_s'+str(isnow)+'grausp')
            das1.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s'+str(isnow)+'grausp.nc', 'w')
            das1.close()
            gc.collect()
            del tb_asrttov_rsgliu_gsp

        elif (eqMass_do == 4):
            das1 = T2P.MHS_as_sims(lons, lats, tb_asrttov_eqMass_rsgliu_gsp[:,:,:], plotpath, server, '_rsg_s'+str(isnow)+'grausp')
            das1.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_eqMass_rsg_s'+str(isnow)+'grausp.nc', 'w')
            das1.close()
            gc.collect()
            del tb_asrttov_eqMass_rsgliu_gsp


    return

def main_Process_Expliu_iwc(instrument, HHtime, mp_version, server, skipProfs, eqMass_do, isnow, igrau, iwc_name):

    if (mp_version == 6):
      mp_physics = 'WRF-WSM6'

    plotpath, folders = config_folders(server)

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

    if mp_version == 6:
        mp_physics = 'WRF-WSM6'
        ncfile     = upfolder+'WRFOUT/WSM6_domain3_NoahMP/wrfout_d02_2018-11-10_'+HHtime+':00'
        ncdata     = Dataset(ncfile,'r')
        [qr, qs, qi, qc, qg, qr_int, qs_int, qi_int, qc_int, qg_int] = funs.get_q_ints6(ncdata)
        int_titles = ['qr','qc','qi','qs','qg']

    flag_name = 'rttov14_'+instrument+'_'+mp_physics+'_2018-11-10_'+HHtime

    # RTTOVout folders
    outfoldereq       = mp_physics+'_20181110_'+HHtime+'_'+instrument +'_atlas_satzen_input__eqmassWSM6_'

    #------        
    # Load all profiles
    A    = read_wrf(ncfile)
    toti = A['XLONG'].shape[0]
    totj = A['XLONG'].shape[1]
    lats         = np.zeros(A['XLONG'].shape); lats[:]=np.nan
    lons         = np.zeros(A['XLONG'].shape); lons[:]=np.nan
    rows, cols = A['XLONG'].shape
    shape_ = A['XLONG'].shape

    if (eqMass_do == 1):
        exp_asrttov_eqmass_rsgliu = eqmass_exp_one_iwc(folders, outfoldereq, mp_physics, HHtime, instrument, '_eqMass_WSM6_rsg', nchan, isnow, igrau, iwc_name)
    elif (eqMass_do == 0):
        exp_asrttov_rsgliu = eqmass_exp_one_iwc(folders, outfoldereq, mp_physics, HHtime, instrument, '_WSM6_rsg', nchan, isnow, igrau, iwc_name)

    outfile           = 'output_tb_'+instrument
    tb_asrttov_eqMass_rsg  = np.zeros( (nchan,rows,cols) );   tb_asrttov_eqMass_rsg[:]=np.nan
    tb_asrttov_rsg         = np.zeros( (nchan,rows,cols) );   tb_asrttov_rsg[:]=np.nan
    counter = 0
    rttov_counter = 0

    for idx, (i,j) in enumerate(np.ndindex(shape_)):
        lats[i,j] = A['XLAT'].data[i,j]
        lons[i,j] = A['XLONG'].data[i,j]

        if idx in skipProfs:
            tb_asrttov_eqMass_rsg[:,i,j] = np.nan
            tb_asrttov_rsg[:,i,j] = np.nan

        else:
            if (eqMass_do==1):
                tb_asrttov_eqMass_rsg[:,i,j] = exp_asrttov_eqmass_rsgliu[rttov_counter-1,:]
            elif (eqMass_do==0):
                tb_asrttov_rsg[:,i,j] = exp_asrttov_rsgliu[rttov_counter-1,:]
            rttov_counter=rttov_counter+1


    # Pre-process like this if MHS                
    if 'MHS' in instrument:
        if (eqMass_do == 0):
            das1 = T2P.MHS_as_sims(lons, lats, tb_asrttov_rsg[:,:,:], plotpath, server, '_rsg_s'+str(isnow)+'g'+str(igrau)+iwc_name)
            das1.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s'+str(isnow)+'g'+str(igrau)+iwc_name+'.nc', 'w')
            das1.close()
            gc.collect()
            del tb_asrttov_rsg

        elif (eqMass_do == 1):
            das1 = T2P.MHS_as_sims(lons, lats, tb_asrttov_eqMass_rsg[:,:,:], plotpath, server, '_rsg_s'+str(isnow)+'g'+str(igrau)+iwc_name)
            das1.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_eqMass_rsg_s'+str(isnow)+'g'+str(igrau)+iwc_name+'.nc', 'w')
            das1.close()
            gc.collect()

    return


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#main_makeProfs('MHS', '20:30', 6, 'cnrm')
#main_Process_cs_andWRF('MHS', '20:30', 6, 'yakaira')

#--------------------------------------------
def main(isnow, igrau, eqMass_do):
    
    server = 'yakaira'
    skipProfs = filter_pixels_monotonic(6, '20:30', server)   
    main_Process_Expliu('MHS', '20:30', 6, server, skipProfs, eqMass_do=eqMass_do, isnow=isnow, igrau=igrau)
    print('Finished running for isnow: '+str(isnow)+' and igrau: '+str(igrau))

def main_basic():
    server = 'yakaira'
    main_Process_cs_andWRF('MHS', '20:30', 6, server)
    print('Finished clear sky and wrf grids')
    
def main_sp(isnow,eqMass):
    server = 'yakaira'
    skipProfs = filter_pixels_monotonic(6, '20:30', server)
    main_Process_Expliu('MHS', '20:30', 6, server, skipProfs, eqMass_do=eqMass, isnow=isnow, igrau=0)
    print('Finished running for isnow: '+str(isnow))

def main_halfiwc(isnow,igrau,eqMass,iwcname):
    server = 'yakaira'
    skipProfs = filter_pixels_monotonic(6, '20:30', server)
    main_Process_Expliu_iwc('MHS', '20:30', 6, server, skipProfs, eqMass_do=eqMass, isnow=isnow, igrau=igrau, iwc_name=iwcname)
    print('Finished running for isnow: '+str(isnow))

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python your_script.py <isnow> <igrau> <?>")
        sys.exit(1)
    
    isnow = int(sys.argv[1])   # or float() if needed
    igrau = int(sys.argv[2])   # or float() if needed
    ieqMass =  int(sys.argv[3]) 
    iwcname =  str(sys.argv[4]) 
    main_basic()
    #main(isnow, igrau, ieqMass)
    #main_sp(isnow,ieqMass)
    #main_halfiwc(isnow, igrau, ieqMass, iwcname)
        
