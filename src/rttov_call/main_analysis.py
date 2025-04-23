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
    
    if mp_version == 6:
        mp_physics = 'WRF-WSM6'
        nchan=5
        
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

    print('Statistics and plots')
            
    # 1st Read nc files 
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
        tb_as_liuliu = [] 
        tb_as_liuliu_gaus = []
        for i in range(11):
            rowi = []
            rowig = [] 
            for  j in range(11):
                expname     = 'rttov_processed_allsky_rsg_s'+str(i)+'g'+str(j)+'.nc'
                d_liuliu    = xr.open_dataset(processedFolder+'/'+outfile+expname)
                var         = d_liuliu['rttov_as'].values
                var_gaus    = d_liuliu['rttov_as_Gaussianantennasigma_'].values
                rowi.append(var)
                rowig.append(var_gaus)
            tb_as_liuliu_gaus.append(rowig)
            tb_as_liuliu.append(rowi)
        tb_as_liuliu = np.array(tb_as_liuliu)
        tb_as_liuliu_gaus = np.array(tb_as_liuliu_gaus)
         
        #if (do_map_plot == 1):
        T2P.make_obsdummy_liuliu_maps(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['MHs_domain_obs'].data, 
                                      server, plotpath) 
        #------------------------------------------------------------------------------
        if (do_map_plot == 1):
            for igrau in range(11):
                T2P.make_all_liuliu_maps(lonlon, latlat, tb_as_liuliu[:,igrau,:,:], server, igrau, plotpath) 
                
        #------------------------------------------------------------------------------------------
        # eq. Mass WRF-WSM6 consistent experiment 
        # READ NETCDFs
        tb_as_eqMass_liuliu      = [] 
        tb_as_eqMass_liuliu_gaus = []
        for i in range(11):
            rowi = []
            rowig = []
            for  j in range(11):
                expname     = 'rttov_processed_allsky_eqMass_rsg_s'+str(i)+'g'+str(j)+'.nc'
                d_liuliu    = xr.open_dataset(processedFolder+'/'+outfile+expname)
                var         = d_liuliu['rttov_as'].values
                var_gaus    = d_liuliu['rttov_as_Gaussianantennasigma_'].values
                rowi.append(var)
                rowig.append(var_gaus)
                
            tb_as_eqMass_liuliu.append(rowi)
            tb_as_eqMass_liuliu_gaus.append(rowig)
        tb_as_eqMass_liuliu      = np.array(tb_as_eqMass_liuliu)            
        tb_as_eqMass_liuliu_gaus = np.array(tb_as_eqMass_liuliu_gaus)            
   
        #------------------------------------------------------------------------------
        if (do_map_plot == 1):
            for igrau in range(11):
                title = 'FIX_rttov_nativegrid_eqMassWSM6_rsg_allisnow_wigrau_'
                T2P.make_all_liuliu_maps(lonlon, latlat, tb_as_eqMass_liuliu[:,igrau,:,:], server, igrau, plotpath, title)     

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
            var_cut1 = T2P.get_maskedDomain(extent, d_cs,  tb_as_liuliu[isnow,igrau,:,:,:])
            var_cut2 = T2P.get_maskedDomain(extent, d_cs,  tb_as_eqMass_liuliu[isnow,igrau,:,:,:])
            rowi1.append(var_cut1)
            rowi2.append(var_cut2)    
        var_cut_liu.append(rowi1)
        var_cut_eqMliu.append(rowi2)
    var_cut_liu = np.array(var_cut_liu)            
    var_cut_eqMliu = np.array(var_cut_eqMliu)            

    #------------------------------------------------------------------------------------------
    # I would also like to mask clouds to focus on histograms of scatt! 
    WRF_intTot_cut = T2P.get_maskedDomain2d(extent, d_cs, WRFvars['WRF_intTot'].data) 
    cloudmask_WRF  = np.ma.masked_less_equal(WRF_intTot_cut, 0.1) # point model gri
    cloudmask_obs1 = np.ma.masked_less_equal( d_cs['MHs_domain_obs'][3,:,:]-d_cs['MHs_domain_obs'][4,:,:], 0.2)     # very rought estimate: obs interp grid 
    # I also test a more restrictive condition to MHS observed: d_cs['MHs_domain_obs'][1,:,:]-280)<0
    cloudmask_obs2 = np.ma.masked_greater_equal( d_cs['MHs_domain_obs'][1,:,:]-280, 0)
           
    #- histogram plots
    if (do_map_plot == 1):
        # WSM6
        T2P.make_hists_liu(d_cs, cloudmask_obs2, var_cut_liu, cloudmask_WRF, plotpath, 'WSM6')
        T2P.make_hists_liu_Perisnow(d_cs, cloudmask_obs2, var_cut_liu, cloudmask_WRF, plotpath, 'WSM6')
    
        # eqMass WSM6
        T2P.make_hists_liu(d_cs, cloudmask_obs2, var_cut_eqMliu, cloudmask_WRF, plotpath, 'WSM6_eqMass')
        T2P.make_hists_liu_Perisnow(d_cs, cloudmask_obs2, var_cut_eqMliu, cloudmask_WRF, plotpath, 'WSM6_eqMass')    
    
    # Calculate HDI index at the two differente resolutions
    cloudmask_WRFgaus = np.ma.masked_less_equal(WRFvars['MHSGaussian_intTot'], 1)      # Gaussian interp grid 
    HDI_wrfresolution_eqMass, HDI_wrfgaussian_eqMass = T2P.calc_stats( nchan, d_cs, var_cut_eqMliu, 
                        tb_as_eqMass_liuliu_gaus, cloudmask_obs2, cloudmask_WRF, cloudmask_WRFgaus)
    
    HDI_wrfresolution_WSM6, HDI_wrfgaussian_WSM6 = T2P.calc_stats( nchan, d_cs, var_cut_liu, 
                        tb_as_liuliu_gaus, cloudmask_obs2, cloudmask_WRF, cloudmask_WRFgaus)
    if (do_map_plot == 1):
        T2P.make_plots_hdi(nchan, plotpath, HDI_wrfresolution_eqMass, HDI_wrfgaussian_eqMass,  'eqMass_WSM6')
        T2P.make_plots_hdi(nchan, plotpath, HDI_wrfresolution_WSM6, HDI_wrfgaussian_WSM6,  'WSM6')
    
    breakpoint() 



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





