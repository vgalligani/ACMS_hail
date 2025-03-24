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
sys.path.insert(1,'/home/galliganiv/ACMS_hail/src')

plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  

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
            
    skipProfs = skipProfs_monotonic
    #--------------------------------------------
    # Write profiles or plot quick results
    #--------------------------------------------
    if makeProfs == 1: 
        # Make profiles and run rttov on terminal 
        ipnd=1
        run_IFS_rttov14(ipnd, mp_version, HHtime, instrument, skipProfs, server)

        
    else:
       
        if mp_version == 6:
            mp_physics = 'WRF-WSM6'
            ncfile     = upfolder+'WRFOUT/WSM6_domain3_NoahMP/wrfout_d02_2018-11-10_'+HHtime+':00'
            ncdata     = Dataset(ncfile,'r') 
            [qr, qs, qi, qc, qg, qr_int, qs_int, qi_int, qc_int, qg_int] = funs.get_q_ints6(ncdata)
            int_titles = ['qr','qc','qi','qs','qg']

        flag_name = 'rttov14_'+instrument+'_'+mp_physics+'_2018-11-10_'+HHtime
    
        # Processed data:
        processedFolder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/'+mp_physics
    
        # Read file
        outfolder = mp_physics+'_20181110_'+HHtime+'_'+instrument +'_atlas_satzen_input/'
        outfile   = 'output_cs_tb_'+instrument
        WSM6_file = np.genfromtxt(folders['read_out_dir']+mp_physics+'/'+outfolder+outfile)
        
        # Load all profiles
        A    = read_wrf(ncfile)
        toti = A['XLONG'].shape[0]
        totj = A['XLONG'].shape[1]
        
        lats = np.zeros(A['XLONG'].shape); lats[:]=np.nan
        lons = np.zeros(A['XLONG'].shape); lons[:]=np.nan
        q    = np.zeros(A['h20'].shape);      q[:]=np.nan
        qinttot = np.zeros(A['XLONG'].shape);   qinttot[:]=np.nan  
        rows, cols = A['XLONG'].shape  # Original dimensions
        tb_csrttov   = np.zeros( (nchan,rows,cols) );   tb_csrttov[:]=np.nan 
        tb1   = np.zeros( (nchan,rows,cols) );   tb1[:]=np.nan 
    
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
                    
                    #tb1[:,i,j] = np.nan
                    tb_csrttov[:,i,j]  = np.nan
                    
                else:
                    rttov_counter     = rttov_counter+1
                    q[:,i,j]          = A['h20'].data[:,i,j]
                    tb_csrttov[:,i,j] = WSM6_file[rttov_counter-1,:]
                    #tb1[:,i,j] = WSM6_atlas_file[rttov_counter-1,:]
        
        #----- SIMPLE PLOTS: make plots with integrated qx forzen y rain, y todos los canales
        # y stats. comparar diferencia clear sky con obs.# improves with atlas? 
        if 'MHS' in instrument: 
            
            # Plot simulation (WRF resolution) and Observations (real resolution)
            ds = T2P.MHS_cs_sims(lons, lats, q, tb_csrttov, plotpath, server)
            ds.to_netcdf(processedFolder+'/'+outfile+'rttov_processed.nc', 'w')
            
            #T2P.plot_simple_MHS_comparison(lons, lats, q, tb0, tb1, plotpath, 'mhs',server)
            
        if 'AMSR' in instrument: 
            T2P.plot_simple_AMSR2_TEST(lons, lats, q, tb_csrttov, tb1, plotpath)
            T2P.plot_simple_AMSR2_comparison(lons, lats, q, tb_csrttov, tb1, plotpath, 'amsr-2')

        # I need to conduct the same type of 'averaging' over the different qints. 
        #-------------------------------------------------------------------------
        # 1) Point q_int on lat/lon grid
        #-------------------------------------------------------------------------
        # 2) Interpolated to MHS lat/lon cut domain
        #-------------------------------------------------------------------------
        MHSinterp_qtot_int  = T2P.regrid2D_toMHSgrid(lats, lons, ds, qinttot.reshape(-1))
        MHSinterp_landMask  = T2P.regrid2D_toMHSgrid(lats, lons, ds, A['LANDMASK'].reshape(-1))
        MHSinterp_qr_int    = T2P.regrid2D_toMHSgrid(lats, lons, ds, qr_int.data.reshape(-1))
        MHSinterp_qs_int    = T2P.regrid2D_toMHSgrid(lats, lons, ds, qs_int.data.reshape(-1))
        MHSinterp_qg_int    = T2P.regrid2D_toMHSgrid(lats, lons, ds, qg_int.data.reshape(-1))
        MHSinterp_qi_int    = T2P.regrid2D_toMHSgrid(lats, lons, ds, qi_int.data.reshape(-1))
        MHSinterp_qc_int    = T2P.regrid2D_toMHSgrid(lats, lons, ds, qc_int.data.reshape(-1))
        # 3) Average over footprint 
        #-------------------------------------------------------------------------
        MHS_qtot_int_footprint_mean  = T2P.average_over_footprint(lats, lons, ds, qinttot)
        MHS_qr_int_footprint_mean    = T2P.average_over_footprint(lats, lons, ds, qr_int.data)
        MHS_qs_int_footprint_mean    = T2P.average_over_footprint(lats, lons, ds, qs_int.data)
        MHS_qg_int_footprint_mean    = T2P.average_over_footprint(lats, lons, ds, qg_int.data)
        MHS_qi_int_footprint_mean    = T2P.average_over_footprint(lats, lons, ds, qi_int.data)
        MHS_qc_int_footprint_mean    = T2P.average_over_footprint(lats, lons, ds, qc_int.data)
                
        # Organizar aux data en un dataframe:
        dwrf =  xr.Dataset({
             "wrf_lat":    (["lat","lon"], lats), 
             "wrf_lon":    (["lat","lon"], lons),
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
             "MHS_intqr":  (["lat2","lon2"],  MHSinterp_qr_int),
             "MHS_intqs":  (["lat2","lon2"],  MHSinterp_qs_int),
             "MHS_intqg":  (["lat2","lon2"],  MHSinterp_qg_int),
             "MHS_intqc":  (["lat2","lon2"],  MHSinterp_qc_int),
             "MHS_intqi":  (["lat2","lon2"],  MHSinterp_qi_int),
             "MHS_intTot":  (["lat2","lon2"],  MHSinterp_qtot_int),
             "MHS_LandMask":  (["lat2","lon2"],  MHSinterp_landMask),
             "MHSfootprintmean_intqr":  (["lat2","lon2"],  MHS_qr_int_footprint_mean),
             "MHSfootprintmean_intqs":  (["lat2","lon2"],  MHS_qs_int_footprint_mean),
             "MHSfootprintmean_intqg":  (["lat2","lon2"],  MHS_qg_int_footprint_mean),
             "MHSfootprintmean_intqi":  (["lat2","lon2"],  MHS_qi_int_footprint_mean),
             "MHSfootprintmean_intqc":  (["lat2","lon2"],  MHS_qc_int_footprint_mean),
             "MHSfootprintmean_intTot":  (["lat2","lon2"], MHS_qtot_int_footprint_mean)
             })
        dwrf.to_netcdf(processedFolder+'/'+'wrfdata_processed.nc', 'w')
        
        
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
        import seaborn as sns 
        #kwargs = dict(hist_kws={'alpha':0.6}, kde_kws={'linewidth':2})
        
        #----- SIMPLE PLOTS for stats (ALL)
        titles_axes = ['89.0', '157.0', '183.311$\pm$3', '183.311$\pm$1', '190.311']
        fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[30,7])    
        for i in range(5):
            sns.histplot( data=ds['MHs_domain_obs'].data[i,:,:].flatten(), 
                         kde="true", color='darkblue', ax=axes[i],  label='MHS', stat="density")
            sns.histplot( data=ds['rttov_cs_footprintmean'].data[i,:,:].flatten(), 
                         kde="true", color='darkred', ax=axes[i],  label='rttov_cs (footprint mean)', stat="density")            
            sns.histplot( data=ds['rttov_cs'].data[i,:,:].flatten(), 
                         kde="true", color='red', ax=axes[i],  label='rttov_cs (point model grid)', stat="density")     
            sns.histplot( data=ds['rttov_cs_pointInterp'].data[i,:,:].flatten(), 
                         kde="true", color='magenta', ax=axes[i],  label='rttov_cs (interp obs grid)', stat="density")     
            axes[i].set_title(titles_axes[i]+' GHz')
            axes[i].set_xlim([200,320])
        axes[0].legend()
        axes[0].set_xlabel('Brightness Temperature (K)')
        plt.suptitle('All pixels')
        fig.savefig(plotpath+'/RTTOV/'+'Stats_clearsky_allpixels_init_tests.png', dpi=300,transparent=False)           

        #-----------------------------------------------------------------------------
        #----- SIMPLE PLOTS for stats (cloudmask)        
        cloudmask_1 = np.ma.masked_greater_equal(dwrf['MHSfootprintmean_intTot'], 1) # footprint mean
        cloudmask_2 = np.ma.masked_greater_equal(dwrf['WRF_intTot'], 1) # point model grid
        cloudmask_3 = np.ma.masked_greater_equal(dwrf['MHS_intTot'], 1) # obs interp grid 
        cloudobsmask =  np.ma.masked_less_equal(ds['MHs_domain_obs'][0,:,:], 270) # obs interp grid 

        meanrms = np.zeros((3,5)); meanrms[:]=np.nan
        meansims = np.zeros((3,5)); meansims[:]=np.nan
        meanobs = [] 
        rmsobs = [] 
        fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[30,7])    
        for i in range(5):
            varsim1 = np.ma.array( ds['rttov_cs_footprintmean'].data[i,:,:], mask=cloudmask_1.mask) 
            varsim2 = np.ma.array( ds['rttov_cs'].data[i,:,:],               mask=cloudmask_2.mask) 
            varsim3 = np.ma.array( ds['rttov_cs_pointInterp'].data[i,:,:],   mask=cloudmask_3.mask)   # update to rttov_cs_pointInterpNearest
            varobs  = np.ma.array( ds['MHs_domain_obs'].data[i,:,:],         mask=cloudobsmask.mask) 
            sns.histplot( data=varobs.flatten(),  color='darkblue', fill=False, ax=axes[i], element='step', stat="density", bins=np.arange(200,320,1), label='MHS (BT89<270)')
            sns.histplot( data=varsim1.flatten(), color='darkred', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (footprint mean)')
            sns.histplot( data=varsim2.flatten(), color='red',     fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (point model grid)')
            sns.histplot( data=varsim3.flatten(), color='magenta', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (interp obs grid)')
            axes[i].set_title(titles_axes[i]+' GHz')
            axes[i].set_xlim([200,310])
            meansims[0,i] = np.round( np.nanmean( varsim1.flatten() ) ,2) 
            meansims[1,i] = np.round( np.nanmean( varsim2.flatten() ) ,2) 
            meansims[2,i] = np.round( np.nanmean( varsim3.flatten() ) ,2) 
            meanobs.append( np.round( np.nanmean( varobs.flatten() ), 2))

            meanrms[0,i] = np.round( np.nanstd( varsim1.flatten() ) ,2) 
            meanrms[1,i] = np.round( np.nanstd( varsim2.flatten() ) ,2) 
            meanrms[2,i] = np.round( np.nanstd( varsim3.flatten() ) ,2) 
            rmsobs.append( np.round( np.nanstd( varobs.flatten() ), 2))
            
            string = 'Mean (std) \n'
            string = string + 'rttov_cs footprint av.: ' + str( meansims[0,i] ) + 'K ('+ str( meanrms[0,i] ) + ') \n'
            string = string + 'rttov_cs (model grid): ' + str( meansims[1,i] ) + 'K ('+ str( meanrms[1,i] ) + ') \n'
            string = string + 'rttov_cs (nearest interp): ' + str( meansims[2,i] ) + 'K ('+ str( meanrms[2,i] ) + ') \n'
            string = string + 'obs: ' + str( np.nanmean( varobs.flatten() )) + 'K ('+ str( rmsobs[i] ) + ') \n'
            
            axes[i].text(.02, .98, string, ha='left', va='top', bbox=dict(facecolor='gray', 
                            edgecolor='black', boxstyle='round') ,  transform=axes[i].transAxes)
    
        axes[0].set_xlabel('Brightness Temperature (K)')
        plt.suptitle('Sims w/ CloudFlag (qtot_int < 1)')
        fig.savefig(plotpath+'/RTTOV/'+'Stats_clearsky_qtotthresholds_init_tests_cloud1.png', dpi=300,transparent=False)    
        
        #-----------------------------------------------------------------------------
        #----- SIMPLE PLOTS for stats (cloudmask)        
        cloudmask_1 = np.ma.masked_greater_equal(dwrf['MHSfootprintmean_intTot'], 0.1) # footprint mean
        cloudmask_2 = np.ma.masked_greater_equal(dwrf['WRF_intTot'], 0.1) # point model grid
        cloudmask_3 = np.ma.masked_greater_equal(dwrf['MHS_intTot'], 0.1) # obs interp grid 
        cloudobsmask =  np.ma.masked_less_equal(ds['MHs_domain_obs'][0,:,:], 280) # obs interp grid 

        meanrms = np.zeros((3,5)); meanrms[:]=np.nan
        meansims = np.zeros((3,5)); meansims[:]=np.nan
        meanobs = [] 
        rmsobs = [] 
        fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[30,7])    
        for i in range(5):
            varsim1 = np.ma.array( ds['rttov_cs_footprintmean'].data[i,:,:], mask=cloudmask_1.mask) 
            varsim2 = np.ma.array( ds['rttov_cs'].data[i,:,:],               mask=cloudmask_2.mask) 
            varsim3 = np.ma.array( ds['rttov_cs_pointInterp'].data[i,:,:],   mask=cloudmask_3.mask)   # update to rttov_cs_pointInterpNearest
            varobs  = np.ma.array( ds['MHs_domain_obs'].data[i,:,:],         mask=cloudobsmask.mask) 
            sns.histplot( data=varobs.flatten(),  color='darkblue', fill=False, ax=axes[i], element='step', stat="density", bins=np.arange(200,320,1), label='MHS (BT89<270)')
            sns.histplot( data=varsim1.flatten(), color='darkred', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (footprint mean)')
            sns.histplot( data=varsim2.flatten(), color='red',     fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (point model grid)')
            sns.histplot( data=varsim3.flatten(), color='magenta', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (interp obs grid)')
            axes[i].set_title(titles_axes[i]+' GHz')
            axes[i].set_xlim([200,310])
            meansims[0,i] = np.round( np.nanmean( varsim1.flatten() ) ,2) 
            meansims[1,i] = np.round( np.nanmean( varsim2.flatten() ) ,2) 
            meansims[2,i] = np.round( np.nanmean( varsim3.flatten() ) ,2) 
            meanobs.append( np.round( np.nanmean( varobs.flatten() ), 2))

            meanrms[0,i] = np.round( np.nanstd( varsim1.flatten() ) ,2) 
            meanrms[1,i] = np.round( np.nanstd( varsim2.flatten() ) ,2) 
            meanrms[2,i] = np.round( np.nanstd( varsim3.flatten() ) ,2) 
            rmsobs.append( np.round( np.nanstd( varobs.flatten() ), 2))
            
            string = 'Mean (std) \n'
            string = string + 'rttov_cs footprint av.: ' + str( meansims[0,i] ) + 'K ('+ str( meanrms[0,i] ) + ') \n'
            string = string + 'rttov_cs (model grid): ' + str( meansims[1,i] ) + 'K ('+ str( meanrms[1,i] ) + ') \n'
            string = string + 'rttov_cs (nearest interp): ' + str( meansims[2,i] ) + 'K ('+ str( meanrms[2,i] ) + ') \n'
            string = string + 'obs: ' + str( np.nanmean( varobs.flatten() )) + 'K ('+ str( rmsobs[i] ) + ') \n'
            
            axes[i].text(.02, .98, string, ha='left', va='top', bbox=dict(facecolor='gray', 
                            edgecolor='black', boxstyle='round') ,  transform=axes[i].transAxes)
    
        axes[0].set_xlabel('Brightness Temperature (K)')
        plt.suptitle('Sims w/ CloudFlag (qtot_int < 0.1)')
        fig.savefig(plotpath+'/RTTOV/'+'Stats_clearsky_qtotthresholds_init_tests_cloud01.png', dpi=300,transparent=False) 



        # ojo que griddata es neasrest neighbour: aplicar gaussuab! 
        # limpiar codigo
        

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
#make_profs('MHS', '20:30', 6, 'cnrm')
make_plots('MHS', '20:30', 6, 'cnrm')



            






