#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 15:20:56 2025

@author: galliganiv
"""

import os 
import datetime
import glob
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import wrf

import Plots4Analysis as P4A
import re 
import package_functions as funs
import pyart 
import config_folders
from Plots4Analysis  import pyplot_rings
from matplotlib.patches import Polygon
import netCDF4 as nc
import matplotlib.colors as mcolors
import package_functions
from matplotlib.cm import get_cmap
import pandas as pd 

plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  



def plot_sondeos():
    
    # import sys
    # #reload(sys)
    # #sys.setdefaultencoding('utf-8')
    # sys.path.insert(0,'/media/hernymet/datos/home/Dropbox/Hernan/drylines/programas_utiles/')
    # #from netCDF4 import num2date #, date2num     # Librerías netCDF para abrir los datos
    # #import netCDF4
    # import numpy as np   # Numerical python
    # import datetime
    # from grafica_sondeo_sharppy import grafica_sondeo_sharppy
    # from siphon.simplewebservice.wyoming import WyomingUpperAir
    # from metpy.units import pandas_dataframe_to_unit_arrays 
    
    # dt = datetime.datetime(1993, 4, 14, 0)
    # #est = 'SAZR'
    # #est = 'SAZN'
    # #est = 'SACO'
    # #est = 'SARE'
    # est = 'SAEZ'
    
    # # Read remote sounding data based on time (dt) and station
    # df = WyomingUpperAir.request_data(dt, est)
    # # Create dictionary of united arrays
    # data = pandas_dataframe_to_unit_arrays(df)
    
    # path_salida = '/media/hernymet/datos/home/Dropbox/Hernan/sondeos_relampago/figuras/'
    
        
    # tiempo_sondeo = dt
    # # Extraemos la latitud y longitud del primer dato (ubicación en superficie)
    # lat_sondeo = df['latitude'][0]
    # lon_sondeo = df['longitude'][0]
    
    # # Abrimos los datos y aplicamos el control de calidad
    # temp_sondeo = df['temperature'].values       # °C
    # td_sondeo = df['dewpoint'].values       # °C
    # pres_sondeo = df['pressure'].values      # hPa   ESTÁN TODOS GRILLADOS CON LA MISMA DISTRIBUCIÓN
    # z_sondeo = df['height'].values
    # #z_sondeo = np.flipud(z_sondeo.astype(float))
    
    # u_sondeo = df['u_wind'].values/1.94   # m/s    
    # v_sondeo = df['v_wind'].values/1.94   # m/s 
    
    # u_sondeo[pres_sondeo<100] = np.nan
    # v_sondeo[pres_sondeo<100] = np.nan
    
    # lat = lat_sondeo
    # lon = lon_sondeo
    # name = est + '_' + tiempo_sondeo.strftime("%Y-%m-%d_%H:%M_UTC")
    # grafica_sondeo_sharppy(pres_sondeo, z_sondeo,temp_sondeo,td_sondeo, u_sondeo, v_sondeo, lat, lon,  name, path_salida, skip=1, wyoming = True)

    return


def lista_datos():
    
    # Research Soundings: 
    # CSU, UI1, UI2, CSWR, CACTI (Villa Dolores)
    # Operational Upper-Air Sounding: Cordoba, Mendoza, Neuquen, 
    #                   Santa Maria, Santa Rosa, Santo Domingo CL, 
    #                   Villa Maria del Rio Seco
    stations = [
        {"Sondeo": "Cordoba", "Latitude":-31.298, "Longitude":-64.212}, 
        #{"Sondeo": "Santa Rosa", "Latitude":-36.591, "Longitude":-64.279}, 
        #{"Sondeo": "Villa Maria del Rio Seco", "Latitude":-29.906, "Longitude":-63.726},
        {"Sondeo": "CSU", "Latitude":-31.817, "Longitude":-64.285}, 
        {"Sondeo": "UI1", "Latitude":-31.541, "Longitude":-64.604}, 
        {"Sondeo": "UI2", "Latitude":-31.997, "Longitude":-64.223}, 
        {"Sondeo": "Scout1", "Latitude":-31.729, "Longitude":-63.845},      # OJO TRANSECTA
        {"Sondeo": "Scout2", "Latitude":-31.877, "Longitude":-64.321},       # OJO TRANSECTA
        {"Sondeo": "Villa Dolores", "Latitude":-31.950, "Longitude":-65.150}
        ]
    
    df = pd.DataFrame(stations)

    return df 




#------------------------------------------------------------------------------
def plot_domain_SurfObs(server):
    
   folders   = config_folders.config_folders(server)
   save_dir_compare = folders['save_dir_compare']+'/SurfaceObs/'
   if not os.path.exists(save_dir_compare):
        os.makedirs(save_dir_compare)       
    
   if 'yakaira' in server:
       prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')    
       fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
   elif 'cnrm' in server: 
       prov = np.genfromtxt("provincias.txt", delimiter='')    
       fn = 'etopo1_bedrock.nc'        
         
   df = lista_datos() 
       
   ds = nc.Dataset(fn)
   topo_lat = ds.variables['lat'][:]
   topo_lon = ds.variables['lon'][:]
   topo_dat = ds.variables['Band1'][:]/1e3
   
   lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
   
   radarLAT_RMA1 = -31.441389
   radarLON_RMA1 = -64.191944
   [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)   
   [lat_radius2, lon_radius2] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,220)   
   
   # De referencia plot contour of observed simulation 
   TH_name       = 'TH'
   file          = 'cfrad.20181110_202339.0000_to_20181110_203020.0000_RMA1_0301_01.nc'
   radar         = pyart.io.read(os.path.join(folders['rma1_dir'], file)) 

   # Configure a gatefilter to filter out copolar correlation coefficient values > 0.9
   gatefilter = pyart.filters.GateFilter(radar)
   gatefilter.exclude_transition()
   gatefilter.exclude_equal('RHOHV', 0.9)

   elev=3
   ZHelev18 = radar.get_field(elev, TH_name, copy=False)
   [lats_, lons_, _] = radar.get_gate_lat_lon_alt(elev, reset_gate_coords=False, filter_transitions=False)

   fig, ax = plt.subplots(figsize=(8,8)) 
   pcm = ax.pcolormesh(lons_, lats_, ZHelev18, cmap=P4A.colormaps('ref'), vmin=5,  vmax=65)
   cbar = plt.colorbar(pcm, ax=ax, shrink=1, label='Zh RMA1 elev 3')
   cbar.cmap.set_under('white')
   cbar.cmap.set_under('white')


   ax.grid()       
   ax.plot(prov[:,0],prov[:,1],color='k'); 
   ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
   ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
       
   ax.set_xlim([-65.2,-63.2]); 
   ax.set_ylim([-32.5,-31])
               
   #else:
   #ax.set_xlim([-65.5,-62]); 
   #ax.set_ylim([-35,-31])
   
   #ax.set_xlim([-68,-61]); 
   #ax.set_ylim([-35,-29])
   
   #ax.set_xlim([-70,-55]); 
   #ax.set_ylim([-41,-25])
      
   # agrego contorno de 500 y 1000m
   ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
   
   slatitudes  = df['Latitude'].tolist()
   slongitudes = df['Longitude'].tolist()
   radiosondelabel = df['Sondeo'].tolist()
   
   for i in range(len(slatitudes)):
       ax.scatter(slongitudes[i], slatitudes[i], marker='s', s=100, color='k', label=radiosondelabel[i])   
       ax.text(slongitudes[i], slatitudes[i], radiosondelabel[i], ha='left', va='bottom', fontsize=10, color='blue')
    
   # Shrink axis by 10%
   #box = ax.get_position()
   #ax.set_position([box.x0, box.y0+box.height*0.1, box.width, box.height*0.9])
   #ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3)
   
   ax.set_title('IOP surface and radiosondes obs. (Zh 20:30)')
   plt.show()
   fig.savefig(save_dir_compare+'IOP_observation_domain.png', dpi=300,transparent=False,bbox_inches='tight')
   #plt.close()
   
   return




#------------------------------------------------------------------------------
plot_domain_SurfObs('cnrm')
    