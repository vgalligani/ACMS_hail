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
def main(makeProfs, instrument, HHtime, mp_version): 

    plotpath, folders = config_folders()
    
    #--------------------------------------------
    # Filter some profiles 
    # find_pixels(568.3049)
    #--------------------------------------------
    skipProfs_monotonic = filter_pixels_monotonic()
    skipProfs = filter_pixels()
    for iprofflag in skipProfs_monotonic:    
            skipProfs.append(iprofflag)
    
    #--------------------------------------------
    # Write profiles or plot quick results
    #--------------------------------------------
    if makeProfs == 1: 
        # Make profiles and run rttov on terminal 
        run_IFS_rttov14(ipnd=1, mp_version=mp_version, HHtime=HHtime, sensor_config=instrument, Flagged_Profs=skipProfs)
    
    else:
       
        if mp_version == 6:
            mp_physics = 'WRF-WSM6'
        
        flag_name = 'rttov14_'+mp_physics+'_2018-11-10_'+HHtime
        flag2name = flag_name+'_'+instrument+'satzen_OBSinterp'
        
        
        #atm_em1.0rttov14_WRF-WSM6_2018-11-10_20:30_MHSsatzen_OBSinterp.dat
        
        # Read file of initial reference of fixed emissivity == 1 for satzen=45
        WSM6_file         = np.genfromtxt(folders['read_out_dir']+'WRF_WSM6_20181110_20_AMSR2/output_tb_45_fixedEmissivity')
        
        # Read experiment of interest 
        EXP='check'
        #WSM6_atlas_file   = np.genfromtxt(folders['read_out_dir']+EXP)
    
        skips = [2325, 8947] #, 2325, 2761, 2796, 3647]
        for iprofflag in skips:    
            skipProfs.append(iprofflag)
            
        # WRFOUT 
        upfolder   = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/'
        mp_physics = 'WRF-WSM6'
        ncfolder   = upfolder+'update'#'WRFout_WSM6_v4.5.2/'
        ncfile     = ncfolder+'wrfout_d02_2018-11-10_'+HHtime+':00'
            
        # Load all profiles
        A    = read_wrf(ncfile)
        toti = A['XLONG'].shape[0]
        totj = A['XLONG'].shape[1]
        
        lats = np.zeros(A['XLONG'].shape); lats[:]=np.nan
        lons = np.zeros(A['XLONG'].shape); lons[:]=np.nan
        q    = np.zeros(A['h20'].shape);   q[:]=np.nan
    
        rows, cols = A['XLONG'].shape  # Original dimensions
        tb0   = np.zeros( (14,rows,cols) );   tb0[:]=np.nan 
        tb1   = np.zeros( (14,rows,cols) );   tb1[:]=np.nan 
        
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
                    tb1[:,i,j] = np.nan
                    
                else:
                    rttov_counter = rttov_counter+1
                    q[:,i,j]  = A['h20'].data[:,i,j]
                    #tb0[:,i,j] = WSM6_file[rttov_counter-1,:]
                    #tb1[:,i,j] = WSM6_atlas_file[rttov_counter-1,:]
        
        #----- SIMPLE PLOTS
        T2P.plot_simple_AMSR2_TEST(lons, lats, q, tb0, tb1, plotpath)
        T2P.plot_simple_AMSR2_comparison(lons, lats, q, tb0, tb1, plotpath, 'mhs')
            
    return


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def make_plots(makeProfs, mp_version): 

    # HARD CODEDE AT 20:30 AS A MHS FILE IS USED FOR ZENITHANGLE! 
    main(makeProfs=makeProfs, instrument='MHS', HHtime='20:30', mp_version=mp_version)  
    #main(EXP=EXP, makeProfs=0, instrument='AMSR2', HHtime='20:30')  

    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def make_profs(EXP,mp_version):

    main(EXP=EXP2, makeProfs=1, instrument='AMSR2', HHtime='20', mp_version=mp_version)  
    
    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
EXP0 = ''        # Surface emissivity == 1 (initial test a nadie ya le importa)
EXP1 = 'WRF_WSM6_20181110_20_AMSR2'+ '_atlas' +'/output_tb_45'     # TELSEM ATLAS and fixed satze rttov run
EXP2 = 'WRF_WSM6_20181110_20_AMSR2'+ '_atlas_satzen_input' +'/output_tb_AMSR2'  # TELSEM ATLAS and satzen input AMSR2 rttov run

# from terminal run run_example_fwd_WRF_HAIL_atlas_satzeninput.sh before switching to makeProfs=0 
make_plots(makeProfs=1, mp_version=6)



# no entiendo porque los resultados son dfierentes?  es porque le puse 55 en vez de 45? CHECK 


            






