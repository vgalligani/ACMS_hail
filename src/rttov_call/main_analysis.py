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



plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def cut_extent_2snowsample(isnow_ranges, d_cs, tb_rttov):
    
    # use extent to limit 
    extent = [ -64, -50, -40, -20]
    var_cut    = []
    for isnow in range(isnow_ranges):
        var_cut1,lon_cut, lat_cut = T2P.get_maskedDomain(extent, d_cs,  tb_rttov[isnow,0,:,:,:])
        var_cut.append(var_cut1)
      
    var_cut = np.array(var_cut)            

    return var_cut, lon_cut, lat_cut

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def cut_extent_liuliu_gaus(d_cs, tb_rttov):
    
    # use extent to limit 
    extent = [ -64, -50, -40, -20]
    var_cut    = []
    for isnow in range(11):
        rowi1 = []
        for igrau in range(11):
            var_cut1, lon_cut, lat_cut = T2P.get_maskedDomain_gaussian(extent, d_cs,  tb_rttov[isnow,igrau,:,:,:])
            rowi1.append(var_cut1)
        var_cut.append(rowi1)
    
    var_cut = np.array(var_cut)         
    
    return var_cut, lon_cut, lat_cut
        
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def cut_extent_liuliu(d_cs, tb_rttov):
    
    # use extent to limit 
    extent = [ -64, -50, -40, -20]
    var_cut    = []
    for isnow in range(11):
        rowi1 = []
        for igrau in range(11):
            var_cut1, lon_cut, lat_cut = T2P.get_maskedDomain(extent, d_cs,  tb_rttov[isnow,igrau,:,:,:])
            rowi1.append(var_cut1)
        var_cut.append(rowi1)
    
    var_cut = np.array(var_cut)         
    
    return var_cut, lon_cut, lat_cut

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_experiment_arrays_liuliu(processedFolder, experiment, instrument):
    
    tb_as_liuliu = [] 
    tb_as_liuliu_gaus = []
    for i in range(11):
        rowi = []
        rowig = [] 
        for  j in range(11):
            expname     = experiment+str(i)+'g'+str(j)+'.nc'
            d_liuliu    = xr.open_dataset(processedFolder+'/'+'output_tb_'+instrument+expname)
            var         = d_liuliu['rttov_as'].values
            var_gaus    = d_liuliu['rttov_as_Gaussianantennasigma_'].values
            rowi.append(var)
            rowig.append(var_gaus)
        tb_as_liuliu_gaus.append(rowig)
        tb_as_liuliu.append(rowi)
    tb_as_liuliu = np.array(tb_as_liuliu)
    tb_as_liuliu_gaus = np.array(tb_as_liuliu_gaus)

    return tb_as_liuliu, tb_as_liuliu_gaus

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_experiment_array_grausp(experiment, instrument, processedFolder, lonlon): 
    
    tb_as_eqMass_liu_gsf       = np.zeros((11,1,5,lonlon.shape[0],lonlon.shape[1])); 
    tb_as_eqMass_liu_gsf[:]    = np.nan
    tb_as_eqMass_liu_gsf_gaus  = []

    for i in range(11):
        expname     = experiment+str(i)+'grausp'+'.nc'
        d_liugrausp = xr.open_dataset(processedFolder+'/'+'output_tb_'+instrument+expname)
        var         = d_liugrausp['rttov_as'].values
        var_gaus    = d_liugrausp['rttov_as_Gaussianantennasigma_'].values
        tb_as_eqMass_liu_gsf[i,0,:,:,:] = (var)
        tb_as_eqMass_liu_gsf_gaus.append( var_gaus )
    tb_as_eqMass_liu_gsf_gaus      = np.array(tb_as_eqMass_liu_gsf_gaus)          
    
    return tb_as_eqMass_liu_gsf, tb_as_eqMass_liu_gsf_gaus

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_iwc_experiment(experiment_1, processedFolder, instrument, lonlon):
    
    tb_liuliu_snowhalfiwc     = np.zeros((2,1,5,lonlon.shape[0],lonlon.shape[1])); tb_liuliu_snowhalfiwc[:]=np.nan
    tb_liuliu_grauhalfiwc     = np.zeros((2,1,5,lonlon.shape[0],lonlon.shape[1])); tb_liuliu_grauhalfiwc[:]=np.nan
    tb_liuliu_grausnowhalfiwc = np.zeros((2,1,5,lonlon.shape[0],lonlon.shape[1])); tb_liuliu_grausnowhalfiwc[:]=np.nan
    tb_liuliu_rainhalfiwc = np.zeros((2,1,5,lonlon.shape[0],lonlon.shape[1])); tb_liuliu_rainhalfiwc[:]=np.nan

    tb_liuliu_noiwc     = np.zeros((2,1,5,lonlon.shape[0],lonlon.shape[1])); tb_liuliu_noiwc[:]=np.nan

    liuopts = [3,9]
    counter = 0
    for i in liuopts:
        expname = experiment_1+str(i)+'g'+str(i)+'grau_iwc'+'.nc'  
        d_1     = xr.open_dataset(processedFolder+'/'+'output_tb_'+instrument+expname)
        var     = d_1['rttov_as'].values
        tb_liuliu_grauhalfiwc[counter,0,:,:,:] = var
        
        expname = experiment_1+str(i)+'g'+str(i)+'snow_iwc'+'.nc'  
        d_1     = xr.open_dataset(processedFolder+'/'+'output_tb_'+instrument+expname)
        var     = d_1['rttov_as'].values
        tb_liuliu_snowhalfiwc[counter,0,:,:,:] = var                 

        expname = experiment_1+str(i)+'g'+str(i)+'rain_iwc'+'.nc'  
        d_1     = xr.open_dataset(processedFolder+'/'+'output_tb_'+instrument+expname)
        var     = d_1['rttov_as'].values
        tb_liuliu_rainhalfiwc[counter,0,:,:,:] = var          
        
        expname = experiment_1+str(i)+'g'+str(i)+'grausnow_iwc'+'.nc'  
        d_1     = xr.open_dataset(processedFolder+'/'+'output_tb_'+instrument+expname)
        var     = d_1['rttov_as'].values
        tb_liuliu_grausnowhalfiwc[counter,0,:,:,:] = var  

        expname = experiment_1+str(i)+'g'+str(i)+'noiwc'+'.nc'  
        d_1     = xr.open_dataset(processedFolder+'/'+'output_tb_'+instrument+expname)
        var     = d_1['rttov_as'].values
        tb_liuliu_noiwc[counter,0,:,:,:] = var  
        counter = counter+1
    
    return tb_liuliu_grauhalfiwc, tb_liuliu_snowhalfiwc, tb_liuliu_grausnowhalfiwc, tb_liuliu_noiwc, tb_liuliu_rainhalfiwc

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_singleiwc_experiment(experiment_1, processedFolder, instrument, lonlon):
    
    tb_liuliu_onlysnow  = np.zeros((2,1,5,lonlon.shape[0],lonlon.shape[1])); tb_liuliu_onlysnow[:]=np.nan
    tb_liuliu_onlygrau  = np.zeros((2,1,5,lonlon.shape[0],lonlon.shape[1])); tb_liuliu_onlygrau[:]=np.nan
    tb_liuliu_onlyrain  = np.zeros((1,1,5,lonlon.shape[0],lonlon.shape[1])); tb_liuliu_onlyrain[:]=np.nan
    tb_liuliu_onlyice   = np.zeros((1,1,5,lonlon.shape[0],lonlon.shape[1])); tb_liuliu_onlyice[:]=np.nan

    liuopts = [3,9]
    counter = 0
    for i in liuopts:
        expname = experiment_1+str(i)+'g'+str(i)+'onlygrau'+'.nc'  
        d_1     = xr.open_dataset(processedFolder+'/'+'output_tb_'+instrument+expname)
        var     = d_1['rttov_as'].values
        tb_liuliu_onlygrau[counter,0,:,:,:] = var
        
        expname = experiment_1+str(i)+'g'+str(i)+'onlysnow'+'.nc'  
        d_1     = xr.open_dataset(processedFolder+'/'+'output_tb_'+instrument+expname)
        var     = d_1['rttov_as'].values
        tb_liuliu_onlysnow[counter,0,:,:,:] = var      

        counter = counter+1
    
    i = 9
    expname = experiment_1+str(i)+'g'+str(i)+'onlyrain'+'.nc'  
    d_1     = xr.open_dataset(processedFolder+'/'+'output_tb_'+instrument+expname)
    var     = d_1['rttov_as'].values
    tb_liuliu_onlyrain[0,0,:,:,:] = var          
        
    expname = experiment_1+str(i)+'g'+str(i)+'onlyice'+'.nc'  
    d_1     = xr.open_dataset(processedFolder+'/'+'output_tb_'+instrument+expname)
    var     = d_1['rttov_as'].values
    tb_liuliu_onlyice[0,0,:,:,:] = var  

    
    return tb_liuliu_onlysnow, tb_liuliu_onlygrau, tb_liuliu_onlyrain, tb_liuliu_onlyice

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main_Process_exp(instrument, HHtime, mp_version, server): 

    # Some internal defitions:
    do_map_plot = 0

    plotpath, folders = config_folders(server)
    
    if mp_version == 6:
        mp_physics = 'WRF-WSM6'
        nchan=5
        
    # Select server and folder locations
    #--------------------------------------------------------------------------
    if 'yakaira' in server: 
        upfolder     = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/'
        sys.path.insert(1,'/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/src')
        processedFolder = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RTTOVout/Processed/'+mp_physics
        mhs_noaa19_dir  = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/PMW/GPM/RS/V07/PMW/1C-MHS-NOAA19/2018/11/10/'
                
    elif 'cnrm' in server:
        upfolder    = '/home/galliganiv/'   
        sys.path.insert(1,'/home/galliganiv/ACMS_hail/src')
        processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'+mp_physics
        mhs_noaa19_dir  = '/home/galliganiv/Work/HAILCASE_10112018/PMW/GPM/RS/V07/PMW/1C-MHS-NOAA19/2018/11/10/'
        
        
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

    print('Statistics and plots')
    extent = [ -64, -50, -40, -20]
            
    # 1st Read nc files 
    #------------------------------------------------------------------------------------------
    # WRF data 
    WRFvars   = xr.open_dataset(processedFolder+'/'+'wrfdata_processed.nc')
    lonlon    = WRFvars['wrf_lon'].data
    latlat    = WRFvars['wrf_lat'].data
    nx,ny=np.shape(lonlon)

    #------------------------------------------------------------------------------------------    
    # rttov clearsky simulations
    outfile   = 'output_tb_'+instrument
    d_cs      = xr.open_dataset(processedFolder+'/'+outfile+'rttov_processed_clearsky.nc')

    #------------------------------------------------------------------------------------------    
    # rttov allskytesting (delfat test with official hydro_table. no psd nonsistency and
    # overlap_param == cloud_overlap_2col_weighted and per_hydro_frac == false)
    d_asTest  = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_test.nc')

    #------------------------------------------------------------------------------------------
    # WRF-WSM6 consistent experiment (LiuLiu combinations)
    tb_as_liuliu, tb_as_liuliu_gaus = get_experiment_arrays_liuliu(processedFolder,  'rttov_processed_allsky_rsg_s', instrument)
    print('Read WRF-WSM6 liuliu')
    
    #------------------------------------------------------------------------------------------
    # eq. Mass WRF-WSM6 consistent experiment 
    tb_as_eqMass_liuliu, tb_as_eqMass_liuliu_gaus = get_experiment_arrays_liuliu(processedFolder,  'rttov_processed_allsky_eqMass_rsg_s', instrument)     
    print('Read WRF-eqMass WSM6 liuliu')
    #------------------------------------------------------------------------------
    if (do_map_plot == 1):
        title = 'rttov_nativegrid_eqMassWSM6_rsg_allisnow_wigrau_'
        for igrau in range(11):
            T2P.make_all_liuliu_maps(lonlon, latlat, tb_as_eqMass_liuliu[:,igrau,:,:], server, igrau, plotpath, title)     

        title = 'rttov_nativegrid_WSM6_rsg_allisnow_wigrau_'
        T2P.make_obsdummy_liuliu_maps(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['MHs_domain_obs'].data, server, plotpath) 
        for igrau in range(11):
            T2P.make_all_liuliu_maps(lonlon, latlat, tb_as_liuliu[:,igrau,:,:], server, igrau, plotpath, title) 
            
    #------------------------------------------------------------------------------------------
    # eq. Plot all the SSP options for snow and grau = soft sphere for EqMass and WSM6 options        
    # READ NETCDFs output_tb_MHSrttov_processed_allsky_eqMass_rsg_s9g9grau_iwc.nc
    tb_as_eqMass_liu_gsf, tb_as_eqMass_liu_gsf_gaus = get_experiment_array_grausp('rttov_processed_allsky_eqMass_rsg_s', instrument, processedFolder, lonlon)
    tb_as_liu_gsf, tb_as_liu_gsf_gaus = get_experiment_array_grausp('rttov_processed_allsky_rsg_s', instrument, processedFolder, lonlon)
    print('Read WRF-eqMass and WSM6 liu grau soft shere')

    if (do_map_plot == 1):
        title = 'rttov_nativegrid_eqMassWSM6_rsg_allisnow_wgrausp'
        fig, axes = plt.subplots(nrows=11, ncols=4, constrained_layout=True,figsize=[10,15])
        for isnow in range(11):
            T2P.plot_MHS_row_WRFgrid(lonlon, latlat, tb_as_eqMass_liu_gsf[isnow,0,:,:], server, axes, isnow, 
                                 'EqMass_WSM6_rsg_s'+str(isnow)+'graupel_softsphere')
        fig.savefig(plotpath+'/RTTOV/'+title+'.png', dpi=300,transparent=False)   
        plt.close()
        
        title = 'FIX_rttov_nativegrid_WSM6_rsg_allisnow_wgrausp'
        fig, axes = plt.subplots(nrows=11, ncols=4, constrained_layout=True,figsize=[10,15])
        for isnow in range(11):
            pcm = T2P.plot_MHS_row_WRFgrid(lonlon, latlat, tb_as_liu_gsf[isnow,0,:,:], server, axes, isnow, 
                                 'WSM6_rsg_s'+str(isnow)+'graupel_softsphere')
        cbar = fig.colorbar(pcm, ax=axes,orientation='vertical',fraction=0.2, pad=0.04)
        fig.savefig(plotpath+'/RTTOV/'+title+'.png', dpi=300,transparent=False)   
        plt.close()    

    #------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------
    # eq. Get experiments of half_snow_iwc and half_grau_iwc for ssp=9 and ssp=3 for snow and grau
    # for eqMass and WSM6
    tb_as_eqMass_liuliu_grauhalfiwc, tb_as_eqMass_liuliu_snowhalfiwc, tb_as_eqMass_liuliu_grausnowhalfiwc, tb_as_eqMass_liuliu_noiwc, tb_as_eqMass_liuliu_rainhalfiwc = get_iwc_experiment( 'rttov_processed_allsky_eqMass_rsg_s', processedFolder, instrument, lonlon)
    tb_as_liuliu_grauhalfiwc, tb_as_liuliu_snowhalfiwc, tb_as_liuliu_grausnowhalfiwc, tb_as_liuliu_noiwc, tb_as_liuliu_rainhalfiwc =  get_iwc_experiment( 'rttov_processed_allsky_rsg_s', processedFolder, instrument, lonlon)
    print('Finish reading iwc sensitivity checks - with noiwc')

    tb_eqMass_liuliu_onlysnow, tb_eqMass_liuliu_onlygrau, tb_eqMass_liuliu_onlyrain, tb_eqMass_liuliu_onlyice = get_singleiwc_experiment('rttov_processed_allsky_eqMass_rsg_s', processedFolder, instrument, lonlon)
    tb_WSM6_liuliu_onlysnow, tb_WSM6_liuliu_onlygrau, tb_WSM6_liuliu_onlyrain, tb_WSM6_liuliu_onlyice = get_singleiwc_experiment('rttov_processed_allsky_rsg_s', processedFolder, instrument, lonlon)

    do_map_plot = 1
    if (do_map_plot == 1):  
        # Plot IWC sensitivity maps     (only snow,rain,ice, grau)    
        T2P.make_maps_iwc_onlyexps(lonlon, latlat, tb_as_eqMass_liuliu, tb_eqMass_liuliu_onlysnow, 
                                   tb_eqMass_liuliu_onlygrau, tb_eqMass_liuliu_onlyrain, tb_eqMass_liuliu_onlyice , plotpath, 'eqMassWSM6', server)

        T2P.make_maps_iwc_onlyexps(lonlon, latlat, tb_as_liuliu, tb_WSM6_liuliu_onlysnow, tb_WSM6_liuliu_onlygrau, 
                                   tb_WSM6_liuliu_onlyrain, tb_WSM6_liuliu_onlyice, plotpath, 'WSM6', server)

        # Plot IWC sensitivity difference maps         
        print('Read WRF-eqMass and WSM6 liu iwc half testing')
        T2P.make_DIFF_maps_iwc_onlyexps(lonlon, latlat, tb_as_eqMass_liuliu, tb_eqMass_liuliu_onlysnow, 
                                        tb_eqMass_liuliu_onlygrau, tb_eqMass_liuliu_onlyrain, tb_eqMass_liuliu_onlyice , plotpath, 'eqMassWSM6', server)
        T2P.make_DIFF_maps_iwc_onlyexps(lonlon, latlat, tb_as_liuliu, tb_WSM6_liuliu_onlysnow, tb_WSM6_liuliu_onlygrau, 
                                        tb_WSM6_liuliu_onlyrain, tb_WSM6_liuliu_onlyice, plotpath, 'WSM6', server)
        

    do_map_plot = 0
    if (do_map_plot == 1):  
        # Plot IWC sensitivity maps         
        T2P.make_maps_iwc_exps(lonlon, latlat, tb_as_eqMass_liuliu,tb_as_eqMass_liuliu_snowhalfiwc,tb_as_eqMass_liuliu_grauhalfiwc, 
                       tb_as_eqMass_liuliu_grausnowhalfiwc, tb_as_eqMass_liuliu_noiwc, tb_as_eqMass_liuliu_rainhalfiwc, plotpath, 'eqMassWSM6', server)

        T2P.make_maps_iwc_exps(lonlon, latlat, tb_as_liuliu, tb_as_liuliu_snowhalfiwc, tb_as_liuliu_grauhalfiwc, 
                       tb_as_liuliu_grausnowhalfiwc, tb_as_liuliu_noiwc, tb_as_liuliu_rainhalfiwc, plotpath, 'WSM6', server)

        # Plot IWC sensitivity difference maps         
        print('Read WRF-eqMass and WSM6 liu iwc half testing')
        T2P.make_DIFF_maps_iwc_exps(lonlon, latlat, tb_as_eqMass_liuliu,tb_as_eqMass_liuliu_snowhalfiwc,tb_as_eqMass_liuliu_grauhalfiwc, 
                       tb_as_eqMass_liuliu_grausnowhalfiwc, tb_as_eqMass_liuliu_noiwc, tb_as_eqMass_liuliu_rainhalfiwc, plotpath, 'eqMassWSM6', server)

        T2P.make_DIFF_maps_iwc_exps(lonlon, latlat, tb_as_liuliu, tb_as_liuliu_snowhalfiwc, tb_as_liuliu_grauhalfiwc, 
                       tb_as_liuliu_grausnowhalfiwc, tb_as_liuliu_noiwc, tb_as_liuliu_rainhalfiwc, plotpath, 'WSM6', server)
        
                
    #------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------
    # Get with a smaller domain! 
    var_cut_halfgrau, lon_wrf_cut, lat_wrf_cut = cut_extent_2snowsample(2, d_cs, tb_as_liuliu_grauhalfiwc)
    var_cut_halfsnow, _, _         = cut_extent_2snowsample(2, d_cs, tb_as_liuliu_snowhalfiwc)
    var_cut_eqMass_halfsnow, _, _  = cut_extent_2snowsample(2, d_cs, tb_as_eqMass_liuliu_snowhalfiwc)
    var_cut_eqMass_halfgrau, _, _  = cut_extent_2snowsample(2, d_cs, tb_as_eqMass_liuliu_grauhalfiwc)

    #------------------------------------------------------------------------------------------            
    #------------------------------------------------------------------------------------------
    # Get cloud masks in case i want to keep only clouds 
    WRF_intTot_cut,_,_ = T2P.get_maskedDomain2d(extent, d_cs, WRFvars['WRF_intTot'].data) 
    cloudmask_WRF  = np.ma.masked_less_equal(WRF_intTot_cut, 0.1) # point model gri
    cloudmask_obs1 = np.ma.masked_less_equal( d_cs['MHs_domain_obs'][3,:,:]-d_cs['MHs_domain_obs'][4,:,:], 0.2)     # very rought estimate: obs interp grid 
    # I also test a more restrictive condition to MHS observed: d_cs['MHs_domain_obs'][1,:,:]-280)<0
    cloudmask_obs2 = np.ma.masked_greater_equal( d_cs['MHs_domain_obs'][1,:,:]-280, 0)
    
    #------------------------------------------------------------------------------------------            
    # Also get a smaller domian of simulations 
    #------------------------------------------------------------------------------------------
    var_cut_liu, _, _    = cut_extent_liuliu(d_cs, tb_as_liuliu)
    var_cut_eqMliu, _ , _  = cut_extent_liuliu(d_cs, tb_as_eqMass_liuliu)    

    #------------------------------------------------------------------------------------------
    var_cut_liu_gaus, _, _    = cut_extent_liuliu_gaus(d_cs, tb_as_liuliu_gaus)
    var_cut_eqMliu_gaus, _, _ = cut_extent_liuliu_gaus(d_cs, tb_as_eqMass_liuliu_gaus)    

    #------------------------------------------------------------------------------------------
    var_cut_eqgrausp, _, _ = cut_extent_2snowsample(11, d_cs, tb_as_eqMass_liu_gsf)
    var_cut_grausp, _, _   = cut_extent_2snowsample(11, d_cs, tb_as_liu_gsf)

    #------------------------------------------------------------------------------------------
    var_cut_grauhalf, _, _       = cut_extent_2snowsample(2, d_cs, tb_as_liuliu_grauhalfiwc)
    var_cut_snowhalf, _, _       = cut_extent_2snowsample(2, d_cs, tb_as_liuliu_snowhalfiwc)
    var_cut_eqMasswsm6_snowhalf, _, _ = cut_extent_2snowsample(2, d_cs, tb_as_eqMass_liuliu_snowhalfiwc)
    var_cut_eqMasswsm6_grauhalf, _, _ = cut_extent_2snowsample(2, d_cs, tb_as_eqMass_liuliu_grauhalfiwc)

    if (do_map_plot == 1):  
        T2P.make_obsdummy_liuliu_maps(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['MHs_domain_obs'].data, server, plotpath) 
        
    do_map_plot = 1
    if (do_map_plot == 1):
        fig, axes = plt.subplots(nrows=7, ncols=4, constrained_layout=True,figsize=[15,20])
        T2P.plot_MHS_row_WRFgrid_largeobs(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['MHs_domain_obs'].data, server, axes, 0, 'OBS')        
        #T2P.plot_MHS_row_WRFgrid_large(lon.T, lat.T, Tc.T, server, axes, 0, 'OBS')        
        T2P.plot_MHS_row_WRFgrid_large(lonlon, latlat, tb_as_eqMass_liuliu[3,3,:,:,:], server, axes, 1, 'EqMass_WSM6 (s3g3)')        
        T2P.plot_MHS_row_WRFgrid_large(lonlon, latlat, tb_as_eqMass_liuliu[9,9,:,:,:], server, axes, 2, 'EqMass WSM6 (s9g9)')       
        T2P.plot_MHS_row_WRFgrid_large(lonlon, latlat, tb_as_eqMass_liu_gsf[9,0,:,:,:], server, axes, 3, 'EqMass_WSM6 (s9gsp)')       
        T2P.plot_MHS_row_WRFgrid_large(lonlon, latlat, tb_as_liuliu[3,3,:,:,:], server, axes, 4, 'WSM6 (s3g3)')        
        T2P.plot_MHS_row_WRFgrid_large(lonlon, latlat, tb_as_liuliu[9,9,:,:,:], server, axes, 5, 'WSM6 (s9g9)')       
        T2P.plot_MHS_row_WRFgrid_large(lonlon, latlat, tb_as_liu_gsf[9,0,:,:,:], server, axes, 6, 'WSM6 (s9gsp)')       
        fig.savefig(plotpath+'/RTTOV/'+'sample_experiment_maps_zoomout.png', dpi=300,transparent=False)   
        plt.close()

        fig, axes = plt.subplots(nrows=7, ncols=4, constrained_layout=True,figsize=[15,20])
        T2P.plot_MHS_row_WRFgrid_largeobs(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['MHs_domain_obs'].data, server, axes, 0, 'OBS')        
        #T2P.plot_MHS_row_WRFgrid_large(lon.T, lat.T, Tc.T, server, axes, 0, 'OBS')        
        T2P.plot_MHS_row_WRFgrid_large(d_cs['MHS_lon'], d_cs['MHS_lat'], tb_as_eqMass_liuliu_gaus[3,3,:,:,:], server, axes, 1, 'EqMass_WSM6 (s3g3)')        
        T2P.plot_MHS_row_WRFgrid_large(d_cs['MHS_lon'], d_cs['MHS_lat'], tb_as_eqMass_liuliu_gaus[9,9,:,:,:], server, axes, 2, 'EqMass WSM6 (s9g9)')       
        T2P.plot_MHS_row_WRFgrid_large(d_cs['MHS_lon'], d_cs['MHS_lat'], tb_as_eqMass_liu_gsf_gaus[9,:,:,:], server, axes, 3, 'EqMass_WSM6 (s9gsp)')       
        T2P.plot_MHS_row_WRFgrid_large(d_cs['MHS_lon'], d_cs['MHS_lat'], tb_as_liuliu_gaus[3,3,:,:,:], server, axes, 4, 'WSM6 (s3g3)')        
        T2P.plot_MHS_row_WRFgrid_large(d_cs['MHS_lon'], d_cs['MHS_lat'], tb_as_liuliu_gaus[9,9,:,:,:], server, axes, 5, 'WSM6 (s9g9)')       
        T2P.plot_MHS_row_WRFgrid_large(d_cs['MHS_lon'], d_cs['MHS_lat'], tb_as_liu_gsf_gaus[9,:,:,:], server, axes, 6, 'WSM6 (s9gsp)')       
        fig.savefig(plotpath+'/RTTOV/'+'gaussian_sample_experiment_maps_zoomout.png', dpi=300,transparent=False)   
        plt.close()

    #------------------------------------------------------------------------
    # Gaussian wrf
    WRF_gausintTot_cut, lon_gauscrop, lat_gauscrop = T2P.get_maskedDomain2d_gaussian(extent, d_cs, WRFvars['MHSGaussian_intTot'].data) 
    cloudmask_gausWRF  = np.ma.masked_less_equal(WRF_gausintTot_cut, 0.1) # point model gri

    #------------------------------------------------------------------------
    #----- checking cloud masks and domain! 
    check_masksandcuts = 0
    if check_masksandcuts == 1: 
        
        from Tools2Plot import GMI_colormap
        cmaps = GMI_colormap() 
        cmap = plt.cm.viridis  # Colormap
        cmap.set_bad(color='gray')  # Color for NaN values
        fig = plt.figure(figsize=[8,8])
        plt.pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], tb_as_liuliu_gaus[9,9,1,:,:], cmap=cmaps['turbo_r'], shading='auto', vmin=50, vmax=300)
        fig = plt.figure(figsize=[8,8])
        plt.pcolormesh(lonlon, latlat, tb_as_liuliu[9,9,1,:,:], cmap=cmaps['turbo_r'], shading='auto', vmin=50, vmax=300)
    
        fig = plt.figure(figsize=[8,8])
        plt.pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], tb_as_liuliu_gaus[9,9,1,:,:], cmap=cmaps['turbo_r'], shading='auto', vmin=50, vmax=300)
    
        fig = plt.figure(figsize=[8,8])
        plt.pcolormesh(lon_wrf_cut, lat_wrf_cut, var_cut_liu[9,9,1,:,:], cmap=cmaps['turbo_r'], shading='auto', vmin=50, vmax=300)
    
        fig = plt.figure(figsize=[8,8])
        plt.pcolormesh(lon_gauscrop, lat_gauscrop, var_cut_liu_gaus[9,9,1,:,:], cmap=cmaps['turbo_r'], shading='auto', vmin=50, vmax=300)
            
        #--- cut also the gaussian versions! 
        onlycloud = np.ma.array( var_cut_liu[9,9,1,:,:], mask=cloudmask_WRF.mask) 
        fig = plt.figure(figsize=[8,8])
        plt.pcolormesh(lon_wrf_cut, lat_wrf_cut, onlycloud, cmap=cmaps['turbo_r'], shading='auto', vmin=50, vmax=300)
    
        onlycloud_gauss = np.ma.array( var_cut_liu_gaus[9,9,1,:,:], mask=cloudmask_gausWRF.mask) 
        fig = plt.figure(figsize=[8,8])
        plt.pcolormesh(lon_gauscrop, lat_gauscrop, onlycloud_gauss, cmap=cmaps['turbo_r'], shading='auto', vmin=50, vmax=300)
    
        #--- cut also the gaussian versions! 
        # check obs maks 
        fig = plt.figure(figsize=[8,8])
        plt.pcolormesh(d_cs['MHS_lon'][:20,14:35], d_cs['MHS_lat'][:20,14:35],  d_cs['MHs_domain_obs'][1,:20,14:35], cmap=cmaps['turbo_r'], shading='auto', vmin=50, vmax=300)

        onlycloud_mhs = np.ma.array( d_cs['MHs_domain_obs'].data[0,:20,14:35], mask=cloudmask_obs2.mask[:20,14:35])
        fig = plt.figure(figsize=[8,8])
        plt.pcolormesh(d_cs['MHS_lon'][:20,14:35], d_cs['MHS_lat'][:20,14:35], onlycloud_mhs, cmap=cmaps['turbo_r'], shading='auto', vmin=50, vmax=300)

    
    breakpoint() 
    
    
    #- histogram plots
    # Simple log histograms
    do_map_plot = 1
    if (do_map_plot == 1):
        
        maskdo = 0  # all pixels (not just cloudy)        
        T2P.make_hists_liu_Perisnow_grausp_simple(d_cs, cloudmask_obs2, var_cut_eqMliu, var_cut_eqgrausp, cloudmask_WRF, plotpath, 'WSM6_eqMass', maskdo) 
        T2P.make_hists_liu_Perisnow_grausp_simple(d_cs, cloudmask_obs2, var_cut_liu, var_cut_grausp, cloudmask_WRF, plotpath, 'WSM6', maskdo) 

        T2P.make_hists_liu_withgrausp_simple(d_cs, cloudmask_obs2, var_cut_liu, var_cut_grausp , cloudmask_WRF, plotpath, 'WSM6', maskdo) 
        T2P.make_hists_liu_withgrausp_simple(d_cs, cloudmask_obs2, var_cut_eqMliu, var_cut_eqgrausp, cloudmask_WRF, plotpath, 'WSM6_eqMass', maskdo) 
    

        T2P.make_hists_liu_simple(d_cs, cloudmask_obs2, var_cut_liu, cloudmask_WRF, plotpath, 'WSM6', maskdo)
        T2P.make_hists_liu_simple(d_cs, cloudmask_obs2, var_cut_eqMliu, cloudmask_WRF, plotpath, 'WSM6_eqMass', maskdo)

        T2P.make_hists_liu_Perisnow_simple(d_cs, cloudmask_obs2, var_cut_liu, cloudmask_WRF, plotpath, 'WSM6', maskdo)
        T2P.make_hists_liu_Perisnow_simple(d_cs, cloudmask_obs2, var_cut_eqMliu, cloudmask_WRF, plotpath, 'WSM6_eqMass', maskdo)    


    do_map_plot = 1
    # Simple log histograms
    if (do_map_plot == 1):
        #------------------------------------------------------------------------------------------
        # Plot the distributions for the half_snow and half_grau experiments
        base_colors = sns.color_palette('Paired')         
        base_colors[10] = base_colors[10+1]
        all_shades = []

        #------------------------------------------------------
        #---- All in one figure 
        fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True,figsize=[8,8]) 
        axes = axes.flatten()  
        axes[0].plot([],[], linewidth=1.2, linestyle='-', color= base_colors[3], label='Liu 3')
        axes[0].plot([],[], linewidth=1.2, linestyle='-', color= base_colors[9], label='Liu 9')
        axes[0].legend()
        axes[1].plot([],[], linewidth=1.2, linestyle='-', color='gray', label='WRF')
        axes[1].plot([],[], linewidth=1.2, linestyle='-.', color='gray', label='half snow')
        axes[1].plot([],[], linewidth=1.2, linestyle='-.', color='gray', label='half grau')
        axes[1].legend()
        T2P.add_hists_2fig_simple( var_cut_eqMliu[3,3,:,:,:], cloudmask_WRF, base_colors[3], axes, 'WRF (liu=3)', '-', maskdo)    
        T2P.add_hists_2fig_simple( var_cut_eqMass_halfsnow[0,:,:,:], cloudmask_WRF, base_colors[3], axes, 'half_snow (liu=3)', '--', maskdo)    
        T2P.add_hists_2fig_simple( var_cut_eqMass_halfgrau[0,:,:,:], cloudmask_WRF, base_colors[3], axes, 'half_grau (liu=3)', '-.', maskdo)    
        T2P.add_hists_2fig_simple( var_cut_eqMliu[9,9,:,:,:], cloudmask_WRF, base_colors[9], axes, 'WRF (liu=9)', '-', maskdo)    
        T2P.add_hists_2fig_simple( var_cut_eqMass_halfsnow[1,:,:,:], cloudmask_WRF, base_colors[9], axes, 'half_snow (liu=9)', '--', maskdo)    
        T2P.add_hists_2fig_simple( var_cut_eqMass_halfgrau[1,:,:,:], cloudmask_WRF, base_colors[9], axes, 'half_grau (liu=9)', '-.', maskdo)  
        for index, i in enumerate([0,1,3,4]):
            varobs2 = d_cs['MHs_domain_obs'].data[i,:,:]
            varobs2 = varobs2.flatten()  
            counts, bin_edges = np.histogram(varobs2, bins=np.arange(10,300,10),density=True)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            axes[index].semilogy(bin_centers, counts, linewidth=1.2,linestyle='-', color='k', label='Obs.')
        plt.suptitle('eqMass_WSM6')
        if maskdo == 1:
            fig.savefig(plotpath+f'/RTTOV/loghist_cloudy_rttov_distribution_eqMassWSM6_halfiwcs.png', dpi=300,transparent=False)   
        else:
            fig.savefig(plotpath+f'/RTTOV/loghist_all_rttov_distribution_eqMassWSM6_halfiwcs.png', dpi=300,transparent=False)   
            
        #plt.close()
        
        fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True,figsize=[8,8]) 
        axes = axes.flatten()  
        axes[0].plot([],[], linewidth=1.2, linestyle='-', color= base_colors[3], label='Liu 3')
        axes[0].plot([],[], linewidth=1.2, linestyle='-', color= base_colors[9], label='Liu 9')
        axes[0].legend()
        axes[1].plot([],[], linewidth=1.2, linestyle='-', color='gray', label='WRF')
        axes[1].plot([],[], linewidth=1.2, linestyle='-.', color='gray', label='half snow')
        axes[1].plot([],[], linewidth=1.2, linestyle='-.', color='gray', label='half grau')
        axes[1].legend()
        T2P.add_hists_2fig_simple( var_cut_liu[3,3,:,:,:], cloudmask_WRF, base_colors[3], axes, 'WRF (liu=3)', '-', maskdo)    
        T2P.add_hists_2fig_simple( var_cut_halfsnow[0,:,:,:], cloudmask_WRF, base_colors[3], axes, 'half_snow (liu=3)', '--', maskdo)    
        T2P.add_hists_2fig_simple( var_cut_halfgrau[0,:,:,:], cloudmask_WRF, base_colors[3], axes, 'half_grau (liu=3)', '-.', maskdo)    
        T2P.add_hists_2fig_simple( var_cut_liu[9,9,:,:,:], cloudmask_WRF, base_colors[9], axes, 'WRF (liu=9)', '-', maskdo)    
        T2P.add_hists_2fig_simple( var_cut_halfsnow[1,:,:,:], cloudmask_WRF, base_colors[9], axes, 'half_snow (liu=9)', '--', maskdo)    
        T2P.add_hists_2fig_simple( var_cut_halfgrau[1,:,:,:], cloudmask_WRF, base_colors[9], axes, 'half_grau (liu=9)', '-.', maskdo)  
        for index, i in enumerate([0,1,3,4]):
            varobs2 = d_cs['MHs_domain_obs'].data[i,:,:]
            varobs2 = varobs2.flatten()        
            counts, bin_edges = np.histogram(varobs2, bins=np.arange(10,300,10),density=True)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            axes[index].semilogy(bin_centers, counts, linewidth=1.2,linestyle='-', color='k', label='Obs.')
        plt.suptitle('WSM6')
        if maskdo == 1:
            fig.savefig(plotpath+f'/RTTOV/loghist_cloudy_rttov_distribution_WSM6_halfiwcs.png', dpi=300,transparent=False)   
        else:
            fig.savefig(plotpath+f'/RTTOV/loghist_all_rttov_distribution_WSM6_halfiwcs.png', dpi=300,transparent=False)   
            
        #plt.close()        
    #------------------------------------------------------------------------------------------
    # Plot histograms of all posible ways for liu9sp 
    
    
    do_map_plot = 0
    # Hisotgrams with kde seaborn
    if (do_map_plot == 1):
        #----------- AT WRF GRID 
        # WSM6
        T2P.make_hists_liu(d_cs, cloudmask_obs2, var_cut_liu, cloudmask_WRF, plotpath, 'WSM6')
        T2P.make_hists_liu_Perisnow(d_cs, cloudmask_obs2, var_cut_liu, cloudmask_WRF, plotpath, 'WSM6')
    
        # eqMass WSM6
        T2P.make_hists_liu(d_cs, cloudmask_obs2, var_cut_eqMliu, cloudmask_WRF, plotpath, 'WSM6_eqMass')
        T2P.make_hists_liu_Perisnow(d_cs, cloudmask_obs2, var_cut_eqMliu, cloudmask_WRF, plotpath, 'WSM6_eqMass')    

        # Adding soft spheres 
        T2P.make_hists_liu_withgrausp(d_cs, cloudmask_obs2, var_cut_liu, var_cut_grausp , cloudmask_WRF, plotpath, 'WSM6') 
        T2P.make_hists_liu_withgrausp(d_cs, cloudmask_obs2, var_cut_eqMliu, var_cut_eqgrausp, cloudmask_WRF, plotpath, 'WSM6_eqMass') 
    
        T2P.make_hists_liu_Perisnow_grausp(d_cs, cloudmask_obs2, var_cut_eqMliu, var_cut_eqgrausp, cloudmask_WRF, plotpath, 'WSM6_eqMass') 
        T2P.make_hists_liu_Perisnow_grausp(d_cs, cloudmask_obs2, var_cut_liu, var_cut_grausp, cloudmask_WRF, plotpath, 'WSM6') 

        #------------------------------------------------------------------------------------------
        # Plot the distributions for the half_snow and half_grau experiments
        base_colors = sns.color_palette('Paired')         
        base_colors[10] = base_colors[10+1]
        all_shades = []

        #------------------------------------------------------
        #---- All in one figure 
        fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True,figsize=[8,8]) 
        axes = axes.flatten()  
        T2P.add_hists_2fig( var_cut_eqMliu[3,3,:,:,:], cloudmask_WRF, base_colors[3], axes, 'WRF (liu=3)', '-')    
        T2P.add_hists_2fig( var_cut_eqMass_halfsnow[0,:,:,:], cloudmask_WRF, base_colors[3], axes, 'half_snow (liu=3)', '--')    
        T2P.add_hists_2fig( var_cut_eqMass_halfgrau[0,:,:,:], cloudmask_WRF, base_colors[3], axes, 'half_grau (liu=3)', '-.')    
        T2P.add_hists_2fig( var_cut_eqMliu[9,9,:,:,:], cloudmask_WRF, base_colors[9], axes, 'WRF (liu=9)', '-')    
        T2P.add_hists_2fig( var_cut_eqMass_halfsnow[1,:,:,:], cloudmask_WRF, base_colors[9], axes, 'half_snow (liu=9)', '--')    
        T2P.add_hists_2fig( var_cut_eqMass_halfgrau[1,:,:,:], cloudmask_WRF, base_colors[9], axes, 'half_grau (liu=9)', '-.')  
        for index, i in enumerate([0,1,3,4]):
            varobs2 = np.ma.array( d_cs['MHs_domain_obs'].data[i,:20,14:35], mask=cloudmask_obs2.mask[:20,14:35]).flatten()
            varobs2 = varobs2.flatten()        
            sns.kdeplot( data = varobs2, color='k', ax=axes[index], label='Obs', bw_adjust=0.2) 
        plt.suptitle('eqMass_WSM6')
        fig.savefig(plotpath+f'/RTTOV/cloudy_rttov_distribution_eqMassWSM6_halfiwcs.png', dpi=300,transparent=False)   
        plt.close()
        
        fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True,figsize=[8,8]) 
        axes = axes.flatten()  
        T2P.add_hists_2fig( var_cut_liu[3,3,:,:,:], cloudmask_WRF, base_colors[3], axes, 'WRF (liu=3)', '-')    
        T2P.add_hists_2fig( var_cut_halfsnow[0,:,:,:], cloudmask_WRF, base_colors[3], axes, 'half_snow (liu=3)', '--')    
        T2P.add_hists_2fig( var_cut_halfgrau[0,:,:,:], cloudmask_WRF, base_colors[3], axes, 'half_grau (liu=3)', '-.')    
        T2P.add_hists_2fig( var_cut_liu[9,9,:,:,:], cloudmask_WRF, base_colors[9], axes, 'WRF (liu=9)', '-')    
        T2P.add_hists_2fig( var_cut_halfsnow[1,:,:,:], cloudmask_WRF, base_colors[9], axes, 'half_snow (liu=9)', '--')    
        T2P.add_hists_2fig( var_cut_halfgrau[1,:,:,:], cloudmask_WRF, base_colors[9], axes, 'half_grau (liu=9)', '-.')  
        for index, i in enumerate([0,1,3,4]):
            varobs2 = np.ma.array( d_cs['MHs_domain_obs'].data[i,:20,14:35], mask=cloudmask_obs2.mask[:20,14:35]).flatten()
            varobs2 = varobs2.flatten()        
            sns.kdeplot( data = varobs2, color='k', ax=axes[index], label='Obs', bw_adjust=0.2) 
        plt.suptitle('WSM6')
        fig.savefig(plotpath+f'/RTTOV/cloudy_rttov_distribution_WSM6_halfiwcs.png', dpi=300,transparent=False)   
        plt.close()        


    #------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------
    # Calculate HDI index at the two differente resolutions
    cloudmask_WRFgaus = np.ma.masked_less_equal(WRFvars['MHSGaussian_intTot'], 1)      # Gaussian interp grid 
    HDI_wrfresolution_eqMass, HDI_wrfgaussian_eqMass = T2P.calc_stats( nchan, d_cs, var_cut_eqMliu, 
                        tb_as_eqMass_liuliu_gaus, cloudmask_obs2, cloudmask_WRF, cloudmask_WRFgaus)
    
    HDI_wrfresolution_WSM6, HDI_wrfgaussian_WSM6 = T2P.calc_stats( nchan, d_cs, var_cut_liu, 
                        tb_as_liuliu_gaus, cloudmask_obs2, cloudmask_WRF, cloudmask_WRFgaus)
    if (do_map_plot == 1):
        T2P.make_plots_hdi(nchan, plotpath, HDI_wrfresolution_eqMass, HDI_wrfgaussian_eqMass,  'eqMass_WSM6')
        T2P.make_plots_hdi(nchan, plotpath, HDI_wrfresolution_WSM6, HDI_wrfgaussian_WSM6,  'WSM6')
        



    # and make plots
    #make_plots_hdi(nchan, plotpath, HDI_wrfresolution, HDI_wrfgaussian)
        


    # Calculate HDI index at the two differente resolutions maybe choose the least diferente and 



    
    # difference between gauss and not gauss? 
    
    


    #------------------------------------------------------------------------------------------

    #kwargs      = dict({'element':'step', 'fill':False, 'stat':"density", 
    #                    'bins':np.arange(50,320,10), 'alpha':0.8}) #'log_scale':(False,True)}) # maybe stat=density
    # kde_kws = dict({'bw_adjust': 0.1, 'common_norm':False }) 

    # comapre with rttov_as_Gaussianantennasigma_
    #import Tools2Plot as T2P
    # represent in terms of deltaTb for gaussian, mean, etc. and also how much it changes over the most scattering
    # and least scattering experiments the deltaTB. only in term of deltaTB. 
    # also do soft spehre! 
    


    return
        


        
        
        
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
    
    
    

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#main_Process_exp('MHS', '20:30', 6, 'cnrm')



    
    
    
    
#--------------------------------------------
def main():
    
    server     = 'yakaira'
    instrument = 'MHS'
    HHtime     = '20:30'
    mp_version = 6
    main_Process_exp('MHS', '20:30', 6, server)


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





