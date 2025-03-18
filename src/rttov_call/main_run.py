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


import numpy as np
import matplotlib.pyplot as plt
import Tools2Plot as T2P
from Tools2RunRTTOV import read_wrf, run_IFS_rttov14, filter_pixels, find_pixels, filter_pixels_monotonic
from config import config_folders



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main(makeProfs, instrument, HHtime, mp_version, server): 

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

        flag_name = 'rttov14_'+instrument+'_'+mp_physics+'_2018-11-10_'+HHtime
    
        # Read file
        outfolder = mp_physics+'_20181110_'+HHtime+'_'+instrument +'_atlas_satzen_input/'
        outfile   = 'output_cs_tb_'+instrument
        WSM6_file = np.genfromtxt(folders['read_out_dir']+mp_physics+'/'+outfolder+outfile)
        
        #skips = [2325, 8947, 10251, 10247] #, 2325, 2761, 2796, 3647] #13/03/2025: skip 10247?         
        #for iprofflag in skips:    
        #    skipProfs.append(iprofflag)    
        
        # Load all profiles
        A    = read_wrf(ncfile)
        toti = A['XLONG'].shape[0]
        totj = A['XLONG'].shape[1]
        
        lats = np.zeros(A['XLONG'].shape); lats[:]=np.nan
        lons = np.zeros(A['XLONG'].shape); lons[:]=np.nan
        q    = np.zeros(A['h20'].shape);   q[:]=np.nan
    
        rows, cols = A['XLONG'].shape  # Original dimensions
        tb0   = np.zeros( (nchan,rows,cols) );   tb0[:]=np.nan 
        tb1   = np.zeros( (nchan,rows,cols) );   tb1[:]=np.nan 
        
        
        counter = 0
        rttov_counter = 0
        for i in range(A['XLONG'].shape[0]): 
            for j in range(A['XLONG'].shape[1]): 
    
                counter = counter+1
                
                lats[i,j] = A['XLAT'].data[i,j]
                lons[i,j] = A['XLONG'].data[i,j]
                    
                if counter in skipProfs:
                    print('Skipping this profile '+str(counter))
                    q[:,i,j]  = np.nan
                    tb0[:,i,j] = np.nan
                    #tb1[:,i,j] = np.nan
                    
                else:
                    rttov_counter = rttov_counter+1
                    q[:,i,j]  = A['h20'].data[:,i,j]
                    tb0[:,i,j] = WSM6_file[rttov_counter-1,:]
                    #tb1[:,i,j] = WSM6_atlas_file[rttov_counter-1,:]
        
        breakpoint()
        #----- SIMPLE PLOTS: make plots with integrated qx forzen y rain, y todos los canales
        # y stats. comparar diferencia clear sky con obs.# improves with atlas? 
        if 'MHS' in instrument: 
            
            T2P.plot_simple_MHS(lons, lats, q, tb0, plotpath,server)
            #T2P.plot_simple_MHS_comparison(lons, lats, q, tb0, tb1, plotpath, 'mhs',server)
            
        if 'AMSR' in instrument: 
            T2P.plot_simple_AMSR2_TEST(lons, lats, q, tb0, tb1, plotpath)
            T2P.plot_simple_AMSR2_comparison(lons, lats, q, tb0, tb1, plotpath, 'amsr-2')
            
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



            






