#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: vito.galligani
"""

# Importar librerias necesarias 
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import integrate
import wrf
import numpy as np
import warnings
import pyart 
import wradlib as wrl
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import glob
import os 
from PIL import Image
import datetime
import pandas
import matplotlib
import matplotlib.dates as mdates
from matplotlib.cm import get_cmap
import warnings
from matplotlib import ticker

import vito_functions 
import vito_plots as VP
from scipy.optimize import fsolve




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# RUN MAIN

server = 'yakaira'

if 'yakaira' in server: 
    # Leer wrfout y variables de interes: control-ndg
    folder_WSM6  = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_WSM6_v4.5.2_all/'
    folder_P3_3MOM_LF = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_P3_3MOM_LF_v4.5.2/' 
    folder_P3_3MOM_noLF = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_P3_3MOM_noLF_v4.5.2/' 
    savedir = '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/WRF_WSM6/'
    savedirP3 = '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/WRF_P3/'
    csapr2_dir = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/CSAPR2/'
    rma1_dir = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RMA1/'
    cswr_dir = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/CSWRdata/'
    tempsavedir = '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/tmp/'
    output_gif_path = '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/gifs/'
    save_dir_compare = '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/'
    
else:
    # Leer wrfout y variables de interes: control-ndg
    folder_WSM6  = '/home/vito.galligani/disk1/derek/Work/HAILCASE_10112018_datos/WRFout_WSM6_v4.5.2/'
    savedir = '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/WRF_WSM6/'
    csapr2_dir = '/home/vito.galligani/disk1/derek/Work/HAILCASE_10112018_datos/CSAPR2/'
    rma1_dir = '/home/vito.galligani/disk1/derek/Work/HAILCASE_10112018_datos/RMA1/'
    tempsavedir = '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/tmp/'

# wrf files
strfile_WSM6 = [folder_WSM6+'wrfout_d02_2018-11-10_17:00:00', 
           folder_WSM6+'wrfout_d02_2018-11-10_17:40:00',
           folder_WSM6+'wrfout_d02_2018-11-10_19:40:00',
           folder_WSM6+'wrfout_d02_2018-11-10_20:00:00',
           folder_WSM6+'wrfout_d02_2018-11-10_20:20:00']

strfile_P3_3MOM_LF = [folder_P3_3MOM_LF+'wrfout_d02_2018-11-10_17:00:00', 
           folder_P3_3MOM_LF+'wrfout_d02_2018-11-10_17:40:00',
           folder_P3_3MOM_LF+'wrfout_d02_2018-11-10_19:40:00',
           folder_P3_3MOM_LF+'wrfout_d02_2018-11-10_20:00:00',
           folder_P3_3MOM_LF+'wrfout_d02_2018-11-10_20:20:00']

# Datos radar [CSAPR-2]
csapr2_str = [csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.170003.nc',
              csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.174503.nc',
              csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.194503.nc',
              csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.200003.nc',
              csapr2_dir+'corcsapr2cfrppiqcM1.b1.20181110.201746.nc']

# Datos de radar [RMA1]
rma1_str_02 = [rma1_dir+'cfrad.20181110_170419.0000_to_20181110_170542.0000_RMA1_0301_02.nc',
              rma1_dir+'cfrad.20181110_174533.0000_to_20181110_174657.0000_RMA1_0301_02.nc',
              rma1_dir+'cfrad.20181110_194101.0000_to_20181110_194221.0000_RMA1_0301_02.nc',
              rma1_dir+'cfrad.20181110_200546.0000_to_20181110_200709.0000_RMA1_0301_02.nc',
              rma1_dir+'cfrad.20181110_202213.0000_to_20181110_202332.0000_RMA1_0301_02.nc']

rma1_str = [rma1_dir+'cfrad.20181110_170542.0000_to_20181110_171223.0000_RMA1_0301_01.nc',
              rma1_dir+'cfrad.20181110_173842.0000_to_20181110_174523.0000_RMA1_0301_01.nc',
              rma1_dir+'cfrad.20181110_194224.0000_to_20181110_194913.0000_RMA1_0301_01.nc',
              rma1_dir+'cfrad.20181110_200709.0000_to_20181110_201358.0000_RMA1_0301_01.nc',
              rma1_dir+'cfrad.20181110_202339.0000_to_20181110_203020.0000_RMA1_0301_01.nc']

# Datos de radar [DOW6]
dow6_str = [cswr_dir+'cfrad.20181110_170011_DOW6high_v215_s01_el0.49_SUR.nc',
              cswr_dir+'cfrad.20181110_174011_DOW6high_v227_s01_el0.49_SUR.nc',
              cswr_dir+'cfrad.20181110_194512_DOW6high_v260_s02_el0.63_SUR.nc',    
              cswr_dir+'cfrad.20181110_200014_DOW6high_v266_s02_el0.92_SUR.nc',
              cswr_dir+'cfrad.20181110_202015_DOW6high_v274_s02_el0.90_SUR.nc']                                                         # FALTA

# Datos de radar [DOW7]
dow7_str = [cswr_dir+'cfrad.20181110_170011_DOW7high_v216_s01_el0.95_SUR.nc',
            cswr_dir+'cfrad.20181110_172011_DOW7high_v222_s01_el0.97_SUR.nc',
            cswr_dir+'cfrad.20181110_194054_DOW7high_v256_s01_el0.70_SUR.nc',  #LE LLUEVE ENCIMA! 
            cswr_dir+'cfrad.20181110_200023_DOW7high_v264_s03_el1.14_SUR.nc',  #LE LLUEVE ENCIMA! 
            cswr_dir+'cfrad.20181110_202012_DOW7high_v271_s02_el0.99_SUR.nc']  #JUSTO PASÓ



# ========================== Compare wsm6 y p3 common transect ?  ===================
# look first at 20:20 Utc
wrf_file_WSM6       = Dataset(strfile_WSM6[4],'r')
wrf_file_p3_3mom_lf = Dataset(strfile_P3_3MOM_LF[4],'r')
itime = strfile_WSM6[4][-8:-3] 
#
# check transects (plot two transects x= x11, x12, x21, x22)
x_wsm6=np.array([-63.45, -63.45, -63.75, -62.0])
y_wsm6=np.array([-32.0, -33.0, -32.00, -33.35])
x_P3_3mom_LF=np.array([-63.6, -63.6, -63.8, -62.0])
y_P3_3mom_LF=np.array([-32.0, -33.0, -32.00, -33.5])   
    
#VP.plot_common_transect_MAP(wrf_file_WSM6, wrf_file_p3_3mom_lf, itime, dow6_str[4], 
#                         dow7_str[4], csapr2_str[4], save_dir_compare, 'MAXRADAR_WRF', x_wsm6, y_wsm6, x_P3_3mom_LF, y_P3_3mom_LF)

# ver evolucion de max Zh WRF for both WSM6 and P3: desde las 19:30 a 2040 ?
#VP.plot_MAXRADAR_WRF_evolution(folder_WSM6, 107, 'WSM6', save_dir_compare, folder_P3_3MOM_LF, 'P3_3MOM_LF')
#VP.plot_MAXRADAR_WRF_evolution(folder_P3_3MOM_LF, 109, 'P3_3MOM_LF', save_dir_compare, folder_WSM6, 'WSM6')

# ver evolucion de los wrfouts
#VP.make_wrfout_gif(folder_WSM6, 107, 'WSM6', output_gif_path) 
#VP.make_wrfout_gif(folder_P3_3MOM_LF, 109, 'P3_3MOM_LF', output_gif_path) 


# RHI PLOT FOR DOW6, CSPR2 AND RMA1 AT TIME OF SUPERCELL (Aca plot 4 degree elevation?)
dow6 = 'cfrad.20181110_201015_DOW6high_v270_s02_el0.90_SUR.nc'
#
radar       = pyart.io.read(cswr_dir+dow6) 
start_index = radar.sweep_start_ray_index['data'][0] 
end_index   = radar.sweep_end_ray_index['data'][0]	
rlats        = radar.gate_latitude['data'][start_index:end_index] 
rlons        = radar.gate_longitude['data'][start_index:end_index] 
ZH           = radar.fields['DBZHC']['data'][start_index:end_index] 
VEL           = radar.fields['VEL']['data'][start_index:end_index] 
#
radarLAT = radar.latitude['data'][0]
radarLON = radar.longitude['data'][0]


fig = plt.figure(figsize=(8,8)) 
pcm1 = plt.pcolormesh(rlons,rlats,ZH, cmap=VP.colormaps('ref'), vmax=60, vmin=0)    
plt.xlim([-64.4,-64])
plt.ylim([-31.7,-32.1])   
plt.title('DOW6 ZH at 2010 s02 elev 0.90')
plt.colorbar()
plt.show()


fig = plt.figure(figsize=(8,8)) 
pcm1 = plt.pcolormesh(rlons,rlats,VEL, cmap=pyart.graph.cm.BuDRd18, vmin=-20, vmax=20)    
plt.xlim([-64.4,-64])
plt.ylim([-31.7,-32.1])   
plt.title('DOW6 VEL at 2010 s02 elev 0.90')
plt.colorbar()
plt.show()




RMA1 = 'cfrad.20181110_201015_DOW6high_v270_s02_el0.90_SUR.nc'




