#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 16:31:51 2025

@author: galliganiv
"""

mietable_folder = '/home/galliganiv/Work/HAILCASE_10112018/RTTOVinout/Processed/mietables/'
mietable = 'rttov_hydrotable_noaa_mhs_eqmassWSM6_rsg_s10g2.dat'

# define iwc
def get_lwc(i_iwc):
    #wc_slope = 0.01 # for version 13.0 which is the one I used 
    #wc_offset = 301 # for version 13.0 which is the one I used 
    # Try my own values for n_lwc == 10
    wc_slope  = 0.1 
    wc_offset = 2.01
    
    return 1e-3 * 10**(wc_slope * (i_iwc - wc_offset))

# define temperature
def get_temp(i_temp):
    # for snow or grau (n_type) =  177.0    
    temp_offset = 177
    return temp_offset + i_temp

#------------------------------------------------
n_temp = 4   # n_temp = 4 leads to 177, 178, 179, 180. quite high value 177! for ice actually 
n_lwc = 4    # n_lwc  = 4 leads to 9.7e-7, 1e-6, 1.02 e-6, 1.05e-6. 
             #        = 401 leads to values between 9.77e-7 and 0.0098

#for ifreq loop 
#for i_band loop
#for i_mtpe loop

#id like water contents between 0.001 and 0 gm-3 
# what is rrttov doping? what is 9.7e-7 to 0.0098? 

for i_temp in range(n_temp):
    temp = get_temp(i_temp) 
    
    # aca llamo a scattering_one_temp donde
    for i_lwc in range(n_lwc):     
        get_lwc(i_lwc)
        # donde llamo scpectra y scattering_one_wc
        print(r'temperature: {temp}' )
    
