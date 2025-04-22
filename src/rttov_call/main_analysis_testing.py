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


plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main_Process_exp(instrument, HHtime, mp_version, server): 

    # Some internal defitions:
    do_map_plot = 0
    grid = 'wrf_nativ_grid'

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

    print('Statistics and plots')
    
    if mp_version == 6:
        mp_physics = 'WRF-WSM6'
        nchan=5
        
    # 1st Read nc files 
    processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'+mp_physics
    
    #------------------------------------------------------------------------------------------
    # WRF data 
    outfile   = 'output_tb_'+instrument
    WRFvars   = xr.open_dataset(processedFolder+'/'+'wrfdata_processed.nc')
    lonlon    = WRFvars['wrf_lon'].data
    latlat    = WRFvars['wrf_lat'].data
    nx,ny=np.shape(lonlon)

    #------------------------------------------------------------------------------------------    
    # rttov clearsky simulations
    d_cs      = xr.open_dataset(processedFolder+'/'+outfile+'rttov_processed_clearsky.nc')

    #------------------------------------------------------------------------------------------    
    # rttov allskytesting (delfat test with official hydro_table. no psd nonsistency and
    # overlap_param == cloud_overlap_2col_weighted and per_hydro_frac == false)
    d_asTest  = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_test.nc')

    #------------------------------------------------------------------------------------------
    # MAPS on WRF native grid 
    if 'wrf_nativ_grid' in grid:
        # WRF-WSM6 consistent experiment (Repeat with eqMass and all liu-liu combinations)
        # READ NETCDFs
        #------------------------------------------------------------------------------------------
        # eq. Mass WRF-WSM6 consistent experiment 
        # READ NETCDFs
        tb_as_eqMass_liuliu = [] 
        tb_as_eqMass_liuliu_gaus = [] 
        for i in range(11):
            rowi = []
            rowigaus = []
            for  j in range(11):
                expname     = 'rttov_processed_allsky_eqMass_rsg_s'+str(i)+'g'+str(j)+'.nc'
                d_liuliu    = xr.open_dataset(processedFolder+'/'+outfile+expname)
                var         = d_liuliu['rttov_as'].values
                var_gaus    = d_liuliu['rttov_as_Gaussianantennasigma_'].values
                rowi.append(var)
                rowigaus.append(var_gaus)
            tb_as_eqMass_liuliu.append(rowi)
            tb_as_eqMass_liuliu_gaus.append(rowigaus)

        tb_as_eqMass_liuliu      = np.array(tb_as_eqMass_liuliu)            
        tb_as_eqMass_liuliu_gaus = np.array(tb_as_eqMass_liuliu_gaus)            
    
        #------------------------------------------------------------------------------
        if (do_map_plot == 1):
            for igrau in range(11):
                title = 'FIX_rttov_nativegrid_eqMassWSM6_rsg_allisnow_wigrau_'
                T2P.make_all_liuliu_maps(lonlon, latlat, tb_as_eqMass_liuliu[:,igrau,:,:], 'cnrm', igrau, plotpath, title)     

    #------------------------------------------------------------------------------------------
    # All outputs saved in netcdf belong to the same extent data: [ -70, -50, -40, -20]
    # but what about this other;    
    extent = [ -64, -50, -40, -20]
    var_cut_liu    = []
    var_cut_eqMliu = []
    
    for isnow in range(11):
        rowi1 = []
        rowi2 = []
        
        for igrau in range(11):
            var_cut2 = T2P.get_maskedDomain(extent, d_cs,  tb_as_eqMass_liuliu[isnow,igrau,:,:,:])
            rowi2.append(var_cut2)    
        var_cut_eqMliu.append(rowi2)

    var_cut_eqMliu = np.array(var_cut_eqMliu)   

    #------------------------------------------------------------------------------------------
    # I would also like to mask clouds to focus on histograms of scatt! 
    WRF_intTot_cut    = T2P.get_maskedDomain2d(extent, d_cs, WRFvars['WRF_intTot'].data) 
    cloudmask_WRF     = np.ma.masked_less_equal(WRF_intTot_cut, 0.1) # point model gri
    cloudmask_WRFgaus = np.ma.masked_less_equal(WRFvars['MHSGaussian_intTot'], 1)      # Gaussian interp grid 
    cloudmask_obs1    = np.ma.masked_less_equal( d_cs['MHs_domain_obs'][3,:,:]-d_cs['MHs_domain_obs'][4,:,:], 0.2)     # very rought estimate: obs interp grid 
    # I also test a more restrictive condition to MHS observed: d_cs['MHs_domain_obs'][1,:,:]-280)<0
    cloudmask_obs2 = np.ma.masked_greater_equal( d_cs['MHs_domain_obs'][1,:,:]-280, 0)
    
    # eqMass WSM6
    #T2P.make_hists_liu(d_cs, cloudmask_obs2, var_cut_eqMliu, cloudmask_WRF, plotpath, 'WSM6_eqMass_fix')
    #T2P.make_hists_liu_Perisnow(d_cs, cloudmask_obs2, var_cut_eqMliu, cloudmask_WRF, plotpath, 'WSM6_eqMass_fix')    
    
    breakpoint() 
    # Calculate HDI index at the two differente resolutions
    HDI_wrfresolution, HDI_wrfgaussian = T2P.calc_stats( nchan, d_cs, var_cut_eqMliu, 
                        tb_as_eqMass_liuliu_gaus, cloudmask_obs2, cloudmask_WRF, cloudmask_WRFgaus)

    # and make plots
    make_plots_hdi(nchan, plotpath, HDI_wrfresolution, HDI_wrfgaussian)
    


    # Calculate HDI index at the two differente resolutions




    #--
    
    # Histoplot params
    #------------------------------------------------------
    ichan_title = ['89.0', '157.0', '183.311$\pm$3', '190.311']
    chan_indx   = [0,1,3,4]
    base_colors = sns.color_palette('Paired')         
    all_shades = []
    for base in base_colors:
        shades = sns.light_palette(base, n_colors=11, input='rgb') #, reverse=True)
        all_shades.append(shades)
    #------------------------------------------------------
    #---- All in one figure 
    fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True,figsize=[8,8]) 
    axes = axes.flatten()
    for index, i in enumerate(chan_indx):
        varobs2 = np.ma.array( d_cs['MHs_domain_obs'].data[i,:,:], mask=cloudmask_obs2.mask).flatten()
        varobs2 = varobs2.flatten()     
        for igrau in range(11): 
            varliu  = np.ma.array( var_cut_eqMliu[0, igrau,i,:,:], mask=cloudmask_WRF.mask) 
            varliu  = varliu.flatten() 
            sns.kdeplot( data = varliu, color=all_shades[0][igrau], ax=axes[index], bw_adjust=0.1)#**kde_kws) 
            varliu_g  = np.ma.array( tb_as_eqMass_liuliu_gaus[0, igrau,i,:,:], mask=cloudmask_WRFgaus.mask) 
            varliu_g  = varliu_g.flatten() 
            sns.kdeplot( data = varliu_g, color='gray', ax=axes[index], bw_adjust=0.2)#**kde_kws) 

        axes[index].set_title(ichan_title[index]+' GHz')
        axes[index].set_xlim([50,310])
        sns.kdeplot( data = varobs2, color='k', ax=axes[index], label='Obs', bw_adjust=0.2) #**kde_kws) 
    
    #------------------------------------------------------    
    # check the histograms for 89GGHz
    data_obs = np.ma.array( d_cs['MHs_domain_obs'].data[0,:,:], mask=cloudmask_obs2.mask).flatten()
    hist_obs, bin_edges = np.histogram(data_obs[:], bins=np.arange(10,320,5), density=True)

    data_sim = np.ma.array( var_cut_eqMliu[0,0,0,:,:], mask=cloudmask_WRF.mask) 
    data_sim = data_sim.flatten() 
    hist_sim, _  = np.histogram(data_sim[:], bins=np.arange(10,320,5), density=True)
                
    data_sim_gaus = np.ma.array( tb_as_eqMass_liuliu_gaus[0,0,0,:,:], mask=cloudmask_WRFgaus.mask) 
    data_sim_gaus = data_sim_gaus.flatten() 
    hist_sim_gaus, _  = np.histogram(data_sim_gaus[:], bins=np.arange(10,320,5), density=True)
    
    fig = plt.figure(figsize=[8,8])
    plt.plot(hist_obs, 'k-', label='obs')
    plt.plot(hist_sim, color='darkblue', label='wrf grid')
    plt.plot(hist_sim_gaus, color='darkred', label='gaussian')
    
                
                
    #------------------------------------------------------    
    plot as a functoipn of freqwuencyu deltaTB due to cloud-clear                      
                
    T2P.plot_MHS_colorbarlims_poster(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['MHs_domain_obs'].data,'testing_poster_obs',
                              50, 300, plotpath, 'cnrm','testing_poster_obs')
    T2P.plot_MHS_colorbarlims_poster(d_cs['MHS_lon'], d_cs['MHS_lat'], tb_as_eqMass_liuliu_gaus[10,8,:,:,:],'testing_poster_gaussian',
                              50, 300, plotpath, 'cnrm','testing_poster_gaussian')

    
    T2P.plot_MHS(d_cs['wrf_lon'], d_cs['wrf_lat'], d_cs['MHs_domain_obs'].data, 50, 300,
                 'testing_poster', plotpath, 'cnrm', 'testing_poster')


d(lonlon, latlat, tbtb, title, vminn, vmaxx, plotpath, server, exp):
    # fig, axes = plt.subplots(nrows=1, ncols=3, constrained_layout=True,figsize=[15,7])
    # pcm=axes[0].pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['MHs_domain_obs'].data[0,:,:]);
    # plt.colorbar(pcm, ax=axes[0]); axes[0].set_title('OBS')

    

    # comapre with rttov_as_Gaussianantennasigma_
    #import Tools2Plot as T2P
    # represent in terms of deltaTb for gaussian, mean, etc. and also how much it changes over the most scattering
    # and least scattering experiments the deltaTB. only in term of deltaTB. 
    # also do soft spehre! 
    


    return



    # 2) I don't compute stats based on FG as they are too different.  
    # Like skewness. Here I need to chose one way to keep the same resolution! choose gaussian? 
    # and replicate above HDI with the gausssian too !
    FG = np.zeros(nchan, 11, 11); FG[:]=np.nan
    for i in range(nchan):
        for isnow in range(11):
            for igrau in range(11):
                FG_ = d_cs['MHs_domain_obs'].data[i,:,:] - data_sim_rescaled[isnow,igrau,i,:,:]






    # fig, axes = plt.subplots(nrows=1, ncols=3, constrained_layout=True,figsize=[15,7])
    # pcm=axes[0].pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['MHs_domain_obs'].data[0,:,:]);
    # plt.colorbar(pcm, ax=axes[0]); axes[0].set_title('OBS')
    # pcm=axes[1].pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], tb_as_eqMass_liuliu_gaus[0,0,0,:,:]); 
    # plt.colorbar(pcm, ax=axes[1]); axes[1].set_title('sim')
    # pcm=axes[2].pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['MHs_domain_obs'].data[0,:,:] - tb_as_eqMass_liuliu_gaus[0,0,0,:,:]); 
    # plt.colorbar(pcm, ax=axes[2]); axes[2].set_title('diff')
    # for i in range(3):
    #     axes[i].set_ylim([-35,-30])
    #     axes[i].set_xlim([-70,-60])
        
        
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
    
    server     = 'cnrm'
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





