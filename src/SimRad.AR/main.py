#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: vito.galligani

wrf2radar_simple has no regionmask

TO DO RADIAL VELOCITY
        # u, v, w] = w

"""

from SimRadAR import wrf2radar_simple 
import matplotlib.pyplot as plt
from SimRadAR import functions, readWRF
import numpy as np 
import config_folders
import os
from scipy import integrate
from SimRadAR import functions
import pyart 
import Plots4SimRad as P4SimRad

plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  

                    
#------------------------------------------------------------------------------
run_once = 0
#------------------------------------------------------------------------------
EXP   = 'WSM6_domain3_NoahMP'
mp    = 6
time  = '20:00'
rfile = 'cfrad.20181110_200709.0000_to_20181110_201358.0000_RMA1_0301_01.nc'

#------------------------------------------------------------------------------
folders=config_folders.config_folders('yakaira')
wrfoutfile = os.path.join(folders[EXP], 'wrfout_d02_2018-11-10_'+time+':00')

# 1) Plot WRF qxs
P4SimRad.plot_WRF_qxs(EXP, time, mp)                                                                                                                                                                                                             
                                                                                                                                                                                   
# 2) READ RADAR DATA and plot (inside plot PHIDP correct)
radar_folder     = folders['rma1_dir']
filename = os.path.join(radar_folder, rfile)
radar    = pyart.io.read(filename) 
radar = P4SimRad.plot_radar_grid(radar, EXP, time)

# 3) Run RadarSim FOR ALL FREQUENCIES:                
if run_once == 1:                                                                                  
    [lon, lat, qr, qs, qg, z_level, Zh_r_central, Zh_s_central, Zh_g_central, Zh_t_central, 
        Zdr_t_central, Zdr_r_central, Zdr_s_central, Zdr_g_central, 
        KDP_t_central, elev_center, azimuth_center, z_FOR_theta_radar, Zh_r_WRFGRID, 
        Zh_s_WRFGRID, Zh_g_WRFGRID, Zh_tot_WRFGRID] = wrf2radar_simple.main_wrf2radar(radar_site=3, 
                    mp=mp, ncfile=wrfoutfile, ftable=folders['LUT_WSM6'])
                           
# Plot the median of Zhh. vs. ZDR, y Zhh vs. KDP for C-band, S-band and observations
bandID = 1          # default plots for C-band
P4SimRad.plot_medianZhhvsPol(Zh_t_central, Zdr_t_central, KDP_t_central, bandID, radar, EXP, mp, time)

# Grids RadarSim outputs (A TODO ESTO CHECK COLMAX?)
P4SimRad.plot_WRF_grid_all(radar, EXP, time, mp, bandID, Zh_r_WRFGRID, Zh_s_WRFGRID, Zh_g_WRFGRID, Zh_tot_WRFGRID)
P4SimRad.plot_radar_grid_all(lon, lat, radar, EXP, mp, Zh_r_central, Zh_s_central, Zh_g_central, Zh_t_central, bandID, time)
P4SimRad.plot_CvsS(lon, lat, radar, Zh_r_central, Zh_s_central, Zh_g_central,  Zh_t_central, EXP, bandID, time)
P4SimRad.run_and_plot_stats(radar, Zh_r_central, Zh_s_central, Zh_g_central, Zh_t_central, EXP, time)
P4SimRad.plot_WRF_grid_all_ZH1km(radar, EXP, time, mp, bandID, Zh_r_WRFGRID, Zh_s_WRFGRID, Zh_g_WRFGRID, Zh_tot_WRFGRID)



# compare with the WRF simulated DBZ y la versión simple de Juan mejorar banda brillante? 
# se queda medio corto con los valores simulados de reflectividad 

        










     
    



