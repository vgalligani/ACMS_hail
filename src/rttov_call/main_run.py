#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------
@purpose : Run rttov (or ARTS?) for the 2018-11-10 hail case (WRF WSM6 and P3)
@author  : Vito Galligani
@email   : vito.galligani@gmail.com
@condaenv: conda activate 
-----------------------------------------------------------------
@main    : Runs clear-sky RTTOV v13.4 or 14? 
           using run_example_fwd_VITO.sh
          
           remember use for reference: https://github.com/lkugler/RTTOV-WRF/blob/master/rttov_wrf.py
           and also to open remote session:
               nohup python -m spyder_kernels.console -matplotlib='inline' -f=./remotemachine.json &4
               sshfs vito.galligani@yakaira.cima.fcen.uba.ar:/home/vito.galligani/Work remote_project
               
               
@PROFILES: HARD-CODED some profiles that have invalid values. 
profile number =     2324

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

#sys.path.insert(1,'/home/galliganiv/ACMS_hail/src')
sys.path.insert(1,'/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/src')

from package_functions import pressure2height

plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  


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
            nter=counter+1

        return tb       
 
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main(makeProfs, instrument, HHtime, mp_version, server): 

    import package_functions as funs

    plotpath, folders = config_folders(server)
    
    # Select server and folder locations
    #--------------------------------------------------------------------------
    if 'yakaira' in server: 
        upfolder     = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/'
                
    elif 'cnrm' in server:
        upfolder    = '/home/galliganiv/'   
            
    # Select instrument configurations
    #--------------------------------------------------------------------------
    if 'MHS' in instrument: 
        nchan = 5
    elif 'AMSR' in instrument:
        nchan = 14
    
    #--------------------------------------------
    # Filter some profiles 
    # find_pixels(568.3049)
    #--------------------------------------------
    # Filter monotonic pressure and h2o > 0.1E-10
    skipProfs_monotonic = filter_pixels_monotonic(mp_version, HHtime, server)
    #skipProfs = filter_pixels(mp_version, HHtime, server)
    #for iprofflag in skipProfs_monotonic:    
    #        skipProfs.append(iprofflag)
    
    # I added this for yakaira 16/04
    #skipProfs, counts = np.unique(skipProfs_monotonic,return_counts=True)
    skipProfs = skipProfs_monotonic 
    #--------------------------------------------
    # Write profiles or plot quick results
    #--------------------------------------------
    if makeProfs == 1: 
        # Make profiles and run rttov on terminal 
        ipnd=1
        run_IFS_rttov14(ipnd, mp_version, HHtime, instrument, skipProfs, server)

        
    elif makeProfs == 0:
        
        print('Preprocessing and map plots')
       
        if mp_version == 6:
            mp_physics = 'WRF-WSM6'
            ncfile     = upfolder+'WRFOUT/WSM6_domain3_NoahMP/wrfout_d02_2018-11-10_'+HHtime+':00'
            ncdata     = Dataset(ncfile,'r') 
            [qr, qs, qi, qc, qg, qr_int, qs_int, qi_int, qc_int, qg_int] = funs.get_q_ints6(ncdata)
            int_titles = ['qr','qc','qi','qs','qg']

        flag_name = 'rttov14_'+instrument+'_'+mp_physics+'_2018-11-10_'+HHtime
    
       
        # Processed data:
        if 'yakaira' in server: 
            processedFolder = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RTTOVout/Processed/'+mp_physics
 
        elif 'cnrm' in server:
            processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'+mp_physics
                
    
        # RTTOVout folders
        outfolder         = mp_physics+'_20181110_'+HHtime+'_'+instrument +'_atlas_satzen_input/'
        outfoldereq       = mp_physics+'_20181110_'+HHtime+'_'+instrument +'_atlas_satzen_input__eqmassWSM6_'

        # Read clear-sky file file
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
        tb_asrttov_eqMass_rsg   = np.zeros( (11,11,nchan,rows,cols) );   tb_asrttov_eqMass_rsg[:]=np.nan 
        tb_asrttov_rsg          = np.zeros( (11,11,nchan,rows,cols) );   tb_asrttov_rsg[:]=np.nan 

	
        #--- eqmassWSM6_rsg_s10g2: equal mass PSD consistency with WRF and 
        # read for all liu - liu combinations w/ snow and grau
        # default cloud overlap settings as above (+renormalization)
        exp_asrttov_eqmass_rsgliu  = eqmass_exp(folders, outfoldereq, mp_physics, 
		HHtime, instrument, '_eqMass_WSM6_rsg', nchan)
        gc.collect()
	
        #--- WSM6_rsg_s10g2: Dmax PSD consistency with WRF and 
        # read for all liu - liu combinations w/ snow and grau
        # default cloud overlap settings as above (re-normalization)
        exp_asrttov_rsgliu  = eqmass_exp(folders, outfoldereq, mp_physics,
                HHtime, instrument, '_WSM6_rsg', nchan)
        gc.collect()


        # Save some aux WRF data to dataframe
        #domainlons = [-65.5,-62]
        #domainlats = [-33.5,-31.3] 
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
                    tb_asrttov_eqMass_rsg[:,:,:,i,j] = np.nan
                    tb_asrttov_rsg[:,:,:,i,j] = np.nan
                    
                else:
                    rttov_counter     = rttov_counter+1
                    q[:,i,j]          = A['h20'].data[:,i,j]
                    zalt[:,i,j]       = pressure2height(A['pressure'][:,i,j], A['T'][:,i,j])
                    
                    tb_asrttov_test[:,i,j] = WSM6_file_as_test[rttov_counter-1,:]
                    tb_csrttov[:,i,j]      = WSM6_file[rttov_counter-1,:]
                    
                    for ilius in range(11):
                        for iliug in range(11):
                            tb_asrttov_eqMass_rsg[ilius,iliug,:,i,j] = exp_asrttov_eqmass_rsgliu[ilius, iliug, rttov_counter-1,:]
                            tb_asrttov_rsg[ilius,iliug,:,i,j]        = exp_asrttov_rsgliu[ilius, iliug, rttov_counter-1,:]

                    #tb_asrttov_rsg[0,:,i,j]  = WSM6_file_as_rsg_s1g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[1,:,i,j]  = WSM6_file_as_rsg_s2g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[2,:,i,j]  = WSM6_file_as_rsg_s3g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[3,:,i,j]  = WSM6_file_as_rsg_s4g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[4,:,i,j]  = WSM6_file_as_rsg_s5g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[5,:,i,j]  = WSM6_file_as_rsg_s6g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[6,:,i,j]  = WSM6_file_as_rsg_s7g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[7,:,i,j]  = WSM6_file_as_rsg_s8g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[8,:,i,j]  = WSM6_file_as_rsg_s9g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[9,:,i,j]  = WSM6_file_as_rsg_s10g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[10,:,i,j] = WSM6_file_as_rsg_s11g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[11,:,i,j] = WSM6_file_as_rsg_s12g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[12,:,i,j] = WSM6_file_as_rsg_s13g2[rttov_counter-1,:]
                    #tb_asrttov_rsg[13,:,i,j] = WSM6_file_as_rsg_s14g2[rttov_counter-1,:]

                    #tb1[:,i,j] = WSM6_atlas_file[rttov_counter-1,:]
        
        #----- SIMPLE PLOTS: make plots with integrated qx forzen y rain, y todos los canales
        # y stats. comparar diferencia clear sky con obs.# improves with atlas? 
        
        
        if 'MHS' in instrument: 
            
            outfile           = 'output_tb_'+instrument

            # Plot simulation (WRF resolution) and Observations (real resolution)
            ds = T2P.MHS_cs_sims(lons, lats, q, tb_csrttov, plotpath, server)
            ds.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_clearsky.nc', 'w')

            das = T2P.MHS_as_sims(lons, lats, tb_asrttov_test, plotpath, server, '_test')
            das.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_test.nc', 'w')

            #for issp in range(14):
            #    das1 = T2P.MHS_as_sims(lons, lats, tb_asrttov_rsg[issp,:,:,:], plotpath, server, '_rsg_s'+str(issp+1)+'g2')
            #    das1.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s'+str(issp+1)+'g2.nc', 'w')
            #    das1.close()
            for ilius in range(11):
            	for iliug in range(11):
                    das1 = T2P.MHS_as_sims(lons, lats, tb_asrttov_rsg[ilius, iliug,:,:,:], plotpath, server, '_rsg_s'+str(ilius)+'g'+str(iliug))
                    das1.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s'+str(ilius)+'g'+str(iliug)+'.nc', 'w')
                    das1.close()

                    das1 = T2P.MHS_as_sims(lons, lats, tb_asrttov_eqMass_rsg[ilius, iliug,:,:,:], plotpath, server, '_rsg_s'+str(ilius)+'g'+str(iliug))
                    das1.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_eqMass_rsg_s'+str(ilius)+'g'+str(iliug)+'.nc', 'w')
                    das1.close()

            # das2 = T2P.MHS_as_sims(lons, lats, tb_asrttov_rsg_s11g2, plotpath, server, '_rsg_s11g2')
            # das2.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s11g2.nc', 'w')

            # das3 = T2P.MHS_as_sims(lons, lats, tb_asrttov_rsg_s3g2, plotpath, server, '_rsg_s3g2')
            # das3.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s3g2.nc', 'w')

            # das4 = T2P.MHS_as_sims(lons, lats, tb_asrttov_rsg_s7g2, plotpath, server, '_rsg_s7g2')
            # das4.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s7g2.nc', 'w')

            # das5 = T2P.MHS_as_sims(lons, lats, tb_asrttov_rsg_s8g2, plotpath, server, '_rsg_s8g2')
            # das5.to_netcdf(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s8g2.nc', 'w')

            #T2P.plot_simple_MHS_comparison(lons, lats, q, tb0, tb1, plotpath, 'mhs',server)
            
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
        WRF_intqr_gaussian   = T2P.overGaussian(lats, lons, ds, qr_int.dataeqMass_)
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
        

    elif makeProfs == 2:
        
        print('Statistics and plots')
        
        if mp_version == 6:
            mp_physics = 'WRF-WSM6'
            
        # 1st Read nc files 
        processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'+mp_physics
        
        #------------------------------------------------------------------------------------------
        outfile   = 'output_tb_'+instrument
        WRFvars   = Dataset(processedFolder+'/'+'wrfdata_processed.nc')
        # clearsky
        d_cs      = Dataset(processedFolder+'/'+outfile+'rttov_processed_clearsky.nc')

        #------------------------------------------------------------------------------------------    
        # allskytesting (delfat test with official hydro_table. no psd nonsistency and
        # overlap_param == cloud_overlap_2col_weighted and per_hydro_frac == false)
        d_asTest  = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_test.nc')

        #------------------------------------------------------------------------------------------
        # allskytesting exp1        
        # eqmassWSM6_rsg_s10g2: equal mass PSD consistency with WRF and 
        # default cloud overlap settings as above
        d_asExp3   = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s3g2.nc')
        d_asExp7   = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s7g2.nc')
        d_asExp8   = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s8g2.nc')
        d_asExp10  = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s10g2.nc')
        d_asExp11  = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s11g2.nc')
        
        breakpoint()
        plt.pcolormesh(d_asExp3['wrf_lon'],d_asExp3['wrf_lat'], d_asExp3['rttov_as']) 
        
        
        #- Clear-sky pixels: For simulations I use cloudmask and for observations? 
        cloudmask_1 = np.ma.masked_greater_equal(WRFvars['MHSfootprintmean_intTot'], 1) # footprint mean
        cloudmask_2 = np.ma.masked_greater_equal(WRFvars['WRF_intTot'], 1)              # point model grid
        cloudmask_3 = np.ma.masked_greater_equal(WRFvars['MHSinterp_intTot'], 1)        # obs interp grid 
        cloudmask_4 = np.ma.masked_greater_equal(WRFvars['MHSGaussian_intTot'], 1)      # Gaussian interp grid 
        cloudobsmask =  np.ma.masked_greater_equal( d_cs['MHs_domain_obs'][3,:,:]-d_cs['MHs_domain_obs'][4,:,:], 0.2)     # very rought estimate: obs interp grid 
        # I also test a more restrictive condition to MHS observed: d_cs['MHs_domain_obs'][1,:,:]-280)<0
        cloudobsmask2 =  np.ma.masked_less_equal( d_cs['MHs_domain_obs'][1,:,:]-280, 0)
        
        cloudmasks = dict({'footprintmean':cloudmask_1, 'WRFgrid':cloudmask_2, 'MHSinterp':cloudmask_3, 
                         'footprintGaussian':cloudmask_4, 'roughOBS':cloudobsmask2})

        #- Histograms for clear sky (using 1) rough cloud masks and 2) all pixels
        filename = 'Stats_clearsky_allpixels_initHist_tests'
        T2P.simple_histograms(d_cs, 'rttov clear-sky', filename, plotpath)
        #
        filename = 'Stats_clearsky_cloudmasks_initHist_tests'
        T2P.simple_histograms_cloudmasks(d_cs, cloudmasks, 'rttov clear-sky (cloudmask)', filename, plotpath)

        #- Plot stats (mean, min, max, std for each simulation)
        EXPsTitles = ['obs', 'footprint mean', 'wrf grid', 'interp. nearest', 'Gaussian antenna']
        colorExp   = ['darkblue','darkred','red','magenta','darkgreen']
        cldmasks   = [cloudobsmask2, cloudmask_1, cloudmask_2, cloudmask_3, cloudmask_4]
        nchan      = 5 
        filename   = 'Stats_clearsky_cloudmasks_initstats'
        T2P.plot1ststats_cs(EXPsTitles, d_cs, cldmasks, nchan, colorExp, plotpath, filename)

        #----------------------------------------------------------------------
        #Look at min/max/std of clear-sky simulations with Cloudy Flags            
        EXPs        = ['MHs_domain_obs', 'rttov_cs_footprintmean', 'rttov_cs', 'rttov_cs_pointInterpNearest', 'rttov_cs_Gaussianantennasigma_']
        titles_axes = ['89.0', '157.0', '183.311$\pm$1', '183.311$\pm$3', '190.311']
        
        dsmean = np.zeros((nchan,nchan))
        dsstd  = np.zeros((nchan,nchan))
        dsmin  = np.zeros((nchan,nchan))
        dsmax  = np.zeros((nchan,nchan))
        ii = 0
        for EXP in EXPs:    
            dsmean[ii], dsstd[ii], dsmin[ii], dsmax[ii] = T2P.calcStats(d_cs[EXP], nchan, cldmasks[ii]) 
            ii=ii+1
        #----------------------------------------------------------------------
        name  = 'ClearSky_stats' 
        Stats = open(processedFolder+'/'+str(name)+'.txt','w')
        print('-----------------------------------', file = Stats)
        print('--- MHS obs w/ strong CloudFlag ---', file = Stats)
        print(' mean: ', np.round(dsmean[0],3), file = Stats)
        print(' std: ' , np.round(dsstd[0],3), file = Stats)
        print(' min: ' , np.round(dsmin[0],3), file = Stats)
        print(' max: ' , np.round(dsmin[0],3), file = Stats)
        print('--- footprintmean, wrfgrid, pointInterp, Gaussian ---', file = Stats)
        print(' mean',  np.round(dsmean[1],3), np.round(dsmean[2],3), np.round(dsmean[3],3), np.round(dsmean[4],3), file = Stats)
        print(' std',   np.round(dsstd[1],3), np.round(dsstd[2],3), np.round(dsstd[3],3), np.round(dsstd[4],3), file = Stats)
        print(' min: ', np.round(dsmin[1],3), np.round(dsmin[2],3), np.round(dsmin[3],3), np.round(dsmin[4],3), file = Stats)
        print(' max: ', np.round(dsmax[1],3), np.round(dsmax[2],3), np.round(dsmax[3],3), np.round(dsmax[4],3), file = Stats)
        print('-----------------------------------', file = Stats)
        Stats.close()
        #----------------------------------------------------------------------
        name  = 'ClearSky_diffstats' 
        Stats = open(processedFolder+'/'+str(name)+'.txt','w')
        print('-----------------------------------', file = Stats)
        print('--- MHS obs w/ strong CloudFlag diff with Gaussian ---', file = Stats)
        for i in range(nchan):
            obs = np.ma.array( d_cs['MHs_domain_obs'][i,:,:].data, mask=cloudobsmask2.mask) 
            sim = np.ma.array( d_cs['rttov_cs_Gaussianantennasigma_'][i,:,:].data, mask=cloudmask_4.mask) 
            print(r'CH ('+titles_axes[i]+'): '+str(np.nanmean(obs-sim)), file = Stats)
            #--------------------------------------------------- 
            # levels = np.arange(-10,10,1)
            # cmap = plt.get_cmap('viridis', len(levels)-1)
            # norm = mpl.colors.BoundaryNorm(boundaries=levels,ncolors=cmap.N)
            # plt.pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], 
            #                (obs-sim),
            #                cmap=cmap, norm=norm); 
            #--------------------------------------------------- 
        print('-----------------------------------', file = Stats)
        Stats.close()
        
        #----------------------------------------------------------------------
        # Std of the mean footprint size clear sky. 
        T2P.plot_MHSstd(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['rttov_cs_footprintstd'], 'FootprintSTD_rttovclearsky', plotpath, 'cnrm')
        T2P.simple_footprint_Std_histograms(d_cs, 'clearsky Footprint std', 'FootprintSTD_rttovclearskyHIST', plotpath) 


        #----------------------------------------------------------------------
        #----------------------------------------------------------------------
        # All sky histogram checks look at Plots4Interc.py in harddrive
        T2P.plot_MHS(d_asTest['MHS_lon'], d_asTest['MHS_lat'], d_asTest['rttov_as_Gaussianantennasigma_'] ,'test', plotpath, 'cnrm')
        filename = 'Stats_allsky_initHist_tests'
        T2P.simple_histograms_as(d_asTest, 'rttov all-sky (test)', filename, plotpath)

        #----------------------------------------------------------------------
        # Std of the mean footprint size. 
        T2P.plot_MHSstd(d_asTest['MHS_lon'], d_asTest['MHS_lat'], d_asTest['rttov_as_footprintstd'], 'FootprintSTD_rttovallskytest', plotpath, 'cnrm')
        T2P.simple_footprint_allsky_Std_histograms(d_asTest, 'allsky Footprint std', 'FootprintSTD_rttovallskyHIST', plotpath) 

        #- Define CloudSignal dT
        dT_rttov_footprintmean = np.zeros((d_cs['rttov_cs_footprintmean'].shape)); dT_rttov_footprintmean[:]=np.nan
        dT_rttov_pointInterp   = np.zeros((d_cs['rttov_cs_footprintmean'].shape)); dT_rttov_pointInterp[:]=np.nan
        dT_rttov_WRFgrid       = np.zeros((d_cs['rttov_cs'].shape)); dT_rttov_WRFgrid[:]=np.nan
        dT_rttov_Gaussian      = np.zeros((d_cs['rttov_cs_footprintmean'].shape)); dT_rttov_Gaussian[:]=np.nan
        for i in range(nchan):
            dT_rttov_footprintmean[i,:,:] = d_cs['rttov_cs_footprintmean'][i,:,:] - d_asTest['rttov_as_footprintmean'][i,:,:]
            dT_rttov_pointInterp[i,:,:]   = d_cs['rttov_cs_pointInterpNearest'][i,:,:] - d_asTest['rttov_as_pointInterpNearest'][i,:,:]
            dT_rttov_Gaussian[i,:,:]      = d_cs['rttov_cs_Gaussianantennasigma_'][i,:,:] - d_asTest['rttov_as_Gaussianantennasigma_'][i,:,:]
            dT_rttov_WRFgrid[i,:,:]       = d_cs['rttov_cs'][i,:,:] - d_asTest['rttov_as'][i,:,:]
    
        #----------------------------------------------------------------------
        EXPs        = ['MHs_domain_obs', 'rttov_as_footprintmean', 'rttov_as', 'rttov_as_pointInterpNearest', 'rttov_as_Gaussianantennasigma_']
        titles_axes = ['89.0', '157.0', '183.311$\pm$1', '183.311$\pm$3', '190.311']
        
        dsmean = np.zeros((nchan,nchan))
        dsstd  = np.zeros((nchan,nchan))
        dsmin  = np.zeros((nchan,nchan))
        dsmax  = np.zeros((nchan,nchan))
        ii = 0
        for EXP in EXPs:    
            dsmean[ii], dsstd[ii], dsmin[ii], dsmax[ii] = T2P.calcStats_nomask(d_asTest[EXP], nchan) 
            ii=ii+1

        dTmean = np.zeros((nchan,nchan))
        dTstd  = np.zeros((nchan,nchan))
        dTmin  = np.zeros((nchan,nchan))
        dTmax  = np.zeros((nchan,nchan))
        dTmean[1], dTstd[1], dTmin[1], dTmax[1] = T2P.calcStats_nomask(dT_rttov_footprintmean, nchan) 
        dTmean[2], dTstd[2], dTmin[2], dTmax[2] = T2P.calcStats_nomask(dT_rttov_WRFgrid , nchan) 
        dTmean[3], dTstd[3], dTmin[3], dTmax[3] = T2P.calcStats_nomask(dT_rttov_pointInterp, nchan) 
        dTmean[4], dTstd[4], dTmin[4], dTmax[4] = T2P.calcStats_nomask(dT_rttov_Gaussian , nchan) 
            
        #----------------------------------------------------------------------
        # aca lo que importa es el max de los diferentes metodos son muy diferentes!
        #----------------------------------------------------------------------
        name  = 'AllSky_stats' 
        Stats = open(processedFolder+'/'+str(name)+'.txt','w')
        print('--------- BT -----------', file = Stats)
        print('--- footprintmean, wrfgrid, pointInterp, Gaussian ---', file = Stats)
        print(' mean',  np.round(dsmean[1],3), np.round(dsmean[2],3), np.round(dsmean[3],3), np.round(dsmean[4],3), file = Stats)
        print(' std',   np.round(dsstd[1],3), np.round(dsstd[2],3), np.round(dsstd[3],3), np.round(dsstd[4],3), file = Stats)
        print(' min: ', np.round(dsmin[1],3), np.round(dsmin[2],3), np.round(dsmin[3],3), np.round(dsmin[4],3), file = Stats)
        print(' max: ', np.round(dsmax[1],3), np.round(dsmax[2],3), np.round(dsmax[3],3), np.round(dsmax[4],3), file = Stats)
        print('-----------------------------------', file = Stats)
        print('-------- BTclear-BTallsky ---------', file = Stats)
        print('--- footprintmean, wrfgrid, pointInterp, Gaussian ---', file = Stats)
        print(' mean',  np.round(dTmean[1],3), np.round(dTmean[2],3), np.round(dTmean[3],3), np.round(dTmean[4],3), file = Stats)
        print(' std',   np.round(dTstd[1],3), np.round(dTstd[2],3), np.round(dTstd[3],3), np.round(dTstd[4],3), file = Stats)
        print(' min: ', np.round(dTmin[1],3), np.round(dTmin[2],3), np.round(dTmin[3],3), np.round(dTmin[4],3), file = Stats)
        print(' max: ', np.round(dTmax[1],3), np.round(dTmax[2],3), np.round(dTmax[3],3), np.round(dTmax[4],3), file = Stats)
        Stats.close()


        x = d_asTest['rttov_as_footprintmean']
        y = d_asTest['rttov_as_Gaussianantennasigma_'] 
        z = d_asTest['rttov_as_pointInterpNearest']
        #----------------------------------------------------------------------
        fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[30,7])  
        for i in range(5):
            axes[i].scatter(x[i,:,:].ravel(), y[i,:,:].ravel(), color='darkblue', alpha=0.8, label='Gaussian antenna')
            axes[i].scatter(x[i,:,:].ravel(), z[i,:,:].ravel(), color='darkred', alpha=0.3, label='Point Interp.')
            axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
            axes[i].grid(True)
        axes[0].set_xlabel('BT rttov footprint mean')
        axes[0].set_ylabel('BT rttov')
        axes[0].legend()
        fig.savefig(plotpath+'/RTTOV/'+'scatterplot_rttovac_test.png', dpi=300,transparent=False)     
        
        x   = d_asTest['rttov_as_footprintmean']
        y   = d_asTest['rttov_as_Gaussianantennasigma_'] 
        z   = d_asTest['rttov_as_pointInterpNearest']
        obs = d_asTest['MHs_domain_obs']
        #----------------------------------------------------------------------
        fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[30,7])  
        for i in range(5):
            axes[i].scatter(obs[i,:,:].ravel(), x[i,:,:].ravel(), color='darkblue', alpha=0.8, label='Footprint mean')
            axes[i].scatter(obs[i,:,:].ravel(), y[i,:,:].ravel(), color='darkgreen', alpha=0.3, label='Gaussian antenna')
            axes[i].scatter(obs[i,:,:].ravel(), z[i,:,:].ravel(), color='darkred', alpha=0.3, label='Point Interp.')
            axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
            axes[i].grid(True)
        axes[0].set_xlabel('BT obs')
        axes[0].set_ylabel('BT rttov')
        axes[0].legend()
        fig.savefig(plotpath+'/RTTOV/'+'scatterplotwOBS_rttovac_test.png', dpi=300,transparent=False)     


        #----------------------------------------------------------------------
        fig, axes = plt.subplots(nrows=4, ncols=5, constrained_layout=True,figsize=[32,14])  
        for i in range(5):
            axes[0,i].scatter(dT_rttov_footprintmean[i,:,:].ravel(), d_asTest['rttov_as_footprintmin'][i,:,:].ravel(), color='darkblue', alpha=0.8)
            #axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
            axes[0,i].grid(True)
        axes[0,0].set_ylabel('min BT in footprint')
        #
        for i in range(5):
            var = d_asTest['rttov_as_footprintmean'][i,:,:].ravel()  
            axes[1,i].scatter(dT_rttov_footprintmean[i,:,:].ravel(), ( var - d_asTest['rttov_as_footprintmin'][i,:,:].ravel()), color='darkred', alpha=0.8)
            #axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
            axes[1,i].grid(True)
        axes[1,0].set_ylabel('Diff in Footprint - min in footprint')        
        #
        for i in range(5):
            axes[2,i].scatter(dT_rttov_footprintmean[i,:,:].ravel(), d_asTest['rttov_as_footprintstd'][i,:,:].ravel(), color='red', alpha=0.8)
            #axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
            axes[2,i].grid(True)
        axes[2,0].set_ylabel('Footprint std')      
        #
        for i in range(5):
            var = d_asTest['rttov_as_Gaussianantennasigma_'][i,:,:].ravel()  
            axes[3,i].scatter(dT_rttov_footprintmean[i,:,:].ravel(), (var - d_asTest['rttov_as_footprintmean'][i,:,:].ravel()), color='blue', alpha=0.8)
            #axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
            axes[3,i].grid(True)
        axes[3,0].set_ylabel('diff Gaussian-Footprint')   
        axes[3,0].set_xlabel('rttov dT(Footsprintmean) (clear-allsky) ')
        #
        fig.savefig(plotpath+'/RTTOV/'+'scatterplot_FootprintDt_VS_footprintmintest.png', dpi=300,transparent=False)     
        #----------------------------------------------------------------------

        # mask = (dT_rttov_footprintmean[0,:,:] == 0) & ~np.isnan(dT_rttov_footprintmean[0,:,:])
        # indices_ = np.argwhere(mask)
        # dT_rttov_footprintmean[0,11,27]
        # CHECK THE ANTENNA PATTERN OF [11,27]
        
        
    

        












                   
        # obs = d_asTest['MHs_domain_obs'][2,:,:].data
        # sim = d_asTest['rttov_as_Gaussianantennasigma_'][2,:,:].data
        # levels = np.arange(-10,10,1)
        # cmap = plt.get_cmap('viridis', len(levels)-1)
        # norm = mpl.colors.BoundaryNorm(boundaries=levels,ncolors=cmap.N)
        # plt.pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], 
        #                (obs-sim),
        #                cmap=cmap, norm=norm);  
        

        # meanrms = np.zeros((3,5)); meanrms[:]=np.nan
        # meansims = np.zeros((3,5)); meansims[:]=np.nan
        # meanobs = [] 
        # rmsobs = [] 
        # fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[30,7])    
        # for i in range(5):
        #     varsim1 = np.ma.array( ds['rttov_cs_footprintmean'].data[i,:,:], mask=cloudmask_1.mask) 
        #     varsim2 = np.ma.array( ds['rttov_cs'].data[i,:,:],               mask=cloudmask_2.mask) 
        #     varsim3 = np.ma.array( ds['rttov_cs_pointInterp'].data[i,:,:],   mask=cloudmask_3.mask)   # update to rttov_cs_pointInterpNearest
        #     varobs  = np.ma.array( ds['MHs_domain_obs'].data[i,:,:],         mask=cloudobsmask.mask) 
        #     sns.histplot( data=varobs.flatten(),  color='darkblue', fill=False, ax=axes[i], element='step', stat="density", bins=np.arange(200,320,1), label='MHS (BT89<270)')
        #     sns.histplot( data=varsim1.flatten(), color='darkred', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (footprint mean)')
        #     sns.histplot( data=varsim2.flatten(), color='red',     fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (point model grid)')
        #     sns.histplot( data=varsim3.flatten(), color='magenta', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (interp obs grid)')
        #     axes[i].set_title(titles_axes[i]+' GHz')
        #     axes[i].set_xlim([200,310])
        #     meansims[0,i] = np.round( np.nanmean( varsim1.flatten() ) ,2) 
        #     meansims[1,i] = np.round( np.nanmean( varsim2.flatten() ) ,2) 
        #     meansims[2,i] = np.round( np.nanmean( varsim3.flatten() ) ,2) 
        #     meanobs.append( np.round( np.nanmean( varobs.flatten() ), 2))

        #     meanrms[0,i] = np.round( np.nanstd( varsim1.flatten() ) ,2) 
        #     meanrms[1,i] = np.round( np.nanstd( varsim2.flatten() ) ,2) 
        #     meanrms[2,i] = np.round( np.nanstd( varsim3.flatten() ) ,2) 
        #     rmsobs.append( np.round( np.nanstd( varobs.flatten() ), 2))
            
        #     string = 'Mean (std) \n'
        #     string = string + 'rttov_cs footprint av.: ' + str( meansims[0,i] ) + 'K ('+ str( meanrms[0,i] ) + ') \n'
        #     string = string + 'rttov_cs (model grid): ' + str( meansims[1,i] ) + 'K ('+ str( meanrms[1,i] ) + ') \n'
        #     string = string + 'rttov_cs (nearest interp): ' + str( meansims[2,i] ) + 'K ('+ str( meanrms[2,i] ) + ') \n'
        #     string = string + 'obs: ' + str( np.nanmean( varobs.flatten() )) + 'K ('+ str( rmsobs[i] ) + ') \n'
            
        #     axes[i].text(.02, .98, string, ha='left', va='top', bbox=dict(facecolor='gray', 
        #                     edgecolor='black', boxstyle='round') ,  transform=axes[i].transAxes)
    
        # axes[0].set_xlabel('Brightness Temperature (K)')
        # plt.suptitle('Sims w/ CloudFlag (qtot_int < 1)')
        # fig.savefig(plotpath+'/RTTOV/'+'Stats_clearsky_qtotthresholds_init_tests_cloud1.png', dpi=300,transparent=False)    
        
        # #-----------------------------------------------------------------------------
        # #----- SIMPLE PLOTS for stats (cloudmask)        
        # cloudmask_1 = np.ma.masked_greater_equal(dwrf['MHSfootprintmean_intTot'], 0.1) # footprint mean
        # cloudmask_2 = np.ma.masked_greater_equal(dwrf['WRF_intTot'], 0.1) # point model grid
        # cloudmask_3 = np.ma.masked_greater_equal(dwrf['MHS_intTot'], 0.1) # obs interp grid 
        # cloudobsmask =  np.ma.masked_less_equal(ds['MHs_domain_obs'][0,:,:], 280) # obs interp grid 

        # meanrms = np.zeros((3,5)); meanrms[:]=np.nan
        # meansims = np.zeros((3,5)); meansims[:]=np.nan
        # meanobs = [] 
        # rmsobs = [] 
        # fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[30,7])    
        # for i in range(5):
        #     varsim1 = np.ma.array( ds['rttov_cs_footprintmean'].data[i,:,:], mask=cloudmask_1.mask) 
        #     varsim2 = np.ma.array( ds['rttov_cs'].data[i,:,:],               mask=cloudmask_2.mask) 
        #     varsim3 = np.ma.array( ds['rttov_cs_pointInterp'].data[i,:,:],   mask=cloudmask_3.mask)   # update to rttov_cs_pointInterpNearest
        #     varobs  = np.ma.array( ds['MHs_domain_obs'].data[i,:,:],         mask=cloudobsmask.mask) 
        #     sns.histplot( data=varobs.flatten(),  color='darkblue', fill=False, ax=axes[i], element='step', stat="density", bins=np.arange(200,320,1), label='MHS (BT89<270)')
        #     sns.histplot( data=varsim1.flatten(), color='darkred', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (footprint mean)')
        #     sns.histplot( data=varsim2.flatten(), color='red',     fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (point model grid)')
        #     sns.histplot( data=varsim3.flatten(), color='magenta', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (interp obs grid)')
        #     axes[i].set_title(titles_axes[i]+' GHz')
        #     axes[i].set_xlim([200,310])
        #     meansims[0,i] = np.round( np.nanmean( varsim1.flatten() ) ,2) 
        #     meansims[1,i] = np.round( np.nanmean( varsim2.flatten() ) ,2) 
        #     meansims[2,i] = np.round( np.nanmean( varsim3.flatten() ) ,2) 
        #     meanobs.append( np.round( np.nanmean( varobs.flatten() ), 2))

        #     meanrms[0,i] = np.round( np.nanstd( varsim1.flatten() ) ,2) 
        #     meanrms[1,i] = np.round( np.nanstd( varsim2.flatten() ) ,2) 
        #     meanrms[2,i] = np.round( np.nanstd( varsim3.flatten() ) ,2) 
        #     rmsobs.append( np.round( np.nanstd( varobs.flatten() ), 2))
            
        #     string = 'Mean (std) \n'
        #     string = string + 'rttov_cs footprint av.: ' + str( meansims[0,i] ) + 'K ('+ str( meanrms[0,i] ) + ') \n'
        #     string = string + 'rttov_cs (model grid): ' + str( meansims[1,i] ) + 'K ('+ str( meanrms[1,i] ) + ') \n'
        #     string = string + 'rttov_cs (nearest interp): ' + str( meansims[2,i] ) + 'K ('+ str( meanrms[2,i] ) + ') \n'
        #     string = string + 'obs: ' + str( np.nanmean( varobs.flatten() )) + 'K ('+ str( rmsobs[i] ) + ') \n'
            
        #     axes[i].text(.02, .98, string, ha='left', va='top', bbox=dict(facecolor='gray', 
        #                     edgecolor='black', boxstyle='round') ,  transform=axes[i].transAxes)
    
        # axes[0].set_xlabel('Brightness Temperature (K)')
        # plt.suptitle('Sims w/ CloudFlag (qtot_int < 0.1)')
        # fig.savefig(plotpath+'/RTTOV/'+'Stats_clearsky_qtotthresholds_init_tests_cloud01.png', dpi=300,transparent=False) 



    return


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def make_plots(instrument, HHtime, mp_version, server): 

    makeProfs = 0
    main(makeProfs, instrument, HHtime, mp_version, server)  

    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def make_profs(instrument, HHtime, mp_version, server):

    #main(EXP=EXP2, makeProfs=1, instrument='AMSR2', HHtime='20', mp_version=mp_version)  
    makeProfs = 1
    main(makeProfs, instrument, HHtime, mp_version, server)  
    
    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def make_stats(instrument, HHtime, mp_version, server):

    #main(EXP=EXP2, makeProfs=1, instrument='AMSR2', HHtime='20', mp_version=mp_version)  
    makeProfs = 2
    main(makeProfs, instrument, HHtime, mp_version, server)  
    
    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#make_profs('MHS', '20:30', 6, 'cnrm')
make_plots('MHS', '20:30', 6, 'yakaira')
#make_stats('MHS', '20:30', 6, 'cnrm')










            
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





