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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import Hodograph, SkewT
from metpy.units import units
from scipy.spatial import Delaunay
import metpy

plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  

#------------------------------------------------------------------------------
def makesndfigure(ncfile, figuretitle, server, lat, lon):

    from metpy.units import units
    
    xloc, yloc = wrf.to_np(wrf.ll_to_xy(ncfile, lat, lon))

    
    #------------------------------------------
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,220)   
    
    # De referencia plot contour of observed simulation 
    folders       = config_folders.config_folders(server)
    TH_name       = 'TH'
    file          = 'cfrad.20181110_202339.0000_to_20181110_203020.0000_RMA1_0301_01.nc'
    radar         = pyart.io.read(os.path.join(folders['rma1_dir'], file)) 
    
    # Configure a gatefilter to filter out copolar correlation coefficient values > 0.9
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_equal('RHOHV', 0.9)    

    if 'yakaira' in server:
        prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')    
        fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
    elif 'cnrm' in server:
        prov = np.genfromtxt("provincias.txt", delimiter='')    
        fn = 'etopo1_bedrock.nc' 
        
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
    #------------------------------------------

    p  = wrf.getvar(ncfile, "p", units="hPa")
    z  = wrf.getvar(ncfile, "z", units="m")
    t  = wrf.getvar(ncfile, "tc")
    td = wrf.getvar(ncfile, "td", units="degC")
    u  = wrf.getvar(ncfile, "uvmet", units="kt")[0,:]
    v  = wrf.getvar(ncfile, "uvmet", units="kt")[1,:]
    
    p = p.sel(south_north=yloc, west_east=xloc)
    z = z.sel(south_north=yloc, west_east=xloc)
    t = t.sel(south_north=yloc, west_east=xloc)
    td = td.sel(south_north=yloc, west_east=xloc)
    u = u.sel(south_north=yloc, west_east=xloc)
    v = v.sel(south_north=yloc, west_east=xloc)
    
    
    # Create a figure 
    fig = plt.figure(figsize=(12, 9), dpi=300.)

    # Ascribe the figure instance to a metpy SkewT object.
    skew = SkewT(fig)

    # Set the sounding's axis parameters: in this case, -60 to +40 C
    # for the skewed x axis and 1000-100 hPa for the logarithmic y axis.
    # We also specify to place tick marks and labels on both axes.
    skew.ax.set_xlim(-60.,40.)
    skew.ax.set_ylim(1000.,100.)
    skew.ax.tick_params(axis='both',labelsize=14)
    
    # Plot our fields: temperature, then dewpoint, then winds.
    # After doing so, we add a standard set of lines to the
    # plot: dry adiabats, moist adiabats (or pseudoadiabats),
    # and mixing-ratio lines. Finally, we label the axes and
    # give the plot a title.
    skew.plot(p,t,'r',linewidth=4)
    skew.plot(p,td,'g',linewidth=4)
    skew.plot_barbs(p,u,v)
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()
    skew.ax.set_xlabel("Temperature (°C)", fontsize=16, weight='bold')
    skew.ax.set_ylabel("Pressure (hPa)", fontsize=16, weight='bold')
    skew.ax.set_title(figuretitle, fontsize=16, weight='bold')
    
    # agregar algunas variables:
    # # CALCULAMOS ALGUNOS ÍNDICES ASOCIADOS AL SONDEO
    muparcel= metpy.calc.most_unstable_parcel(p,t,td)
    cape_cin = metpy.calc.cape_cin(p,t,td,mpcalc.parcel_profile(p, t[0], td[0]).data.to('degC'))
    mucape_cin = metpy.calc.most_unstable_cape_cin(p,t,td)
    srh = metpy.calc.storm_relative_helicity(z,u,v,depth=3000 * units.m )
    pwat = metpy.calc.precipitable_water(p,td)
    shear_06 = mpcalc.bulk_shear(p,u,v,z,depth=6000 * units.meter)
    shear_03 = mpcalc.bulk_shear(p,u,v,z,depth=3000 * units.meter)
    
    indices = {'CAPE': ['%.0f' % cape_cin[0].magnitude, 'J/kg'],\
    'CIN': ['%.0f' % cape_cin[1].magnitude, 'J/kg'],\
    'MUCAPE': ['%.0f' % mucape_cin[0].magnitude, 'J/kg'],\
    'MUCIN': ['%.0f' % mucape_cin[1].magnitude, 'J/kg'],\
    '0-3 km Shear': ['%.0f' % (np.sqrt(shear_03[0].magnitude**2 + shear_03[1].magnitude**2)*1.94), 'kts'],\
    '0-6 km Shear': ['%.0f' % (np.sqrt(shear_06[0].magnitude**2 + shear_06[1].magnitude**2)*1.94), 'kts'],\
    '0-3 km SRH': ['%.0f' % srh[2].magnitude, 'm2/s2'],\
    'PW': ['%.1f' % pwat.magnitude, 'mm']}        

        
    string = ''
    for key in (indices.keys()):
        #print(key)
        string = string + key + ': ' + str(indices[key][0]) + ' ' + indices[key][1] + '\n'
    skew.ax.text(.01, .99, string, ha='left', va='top', bbox=dict(facecolor='gray', 
                            edgecolor='black', boxstyle='round'),  transform=plt.gca().transAxes)
        
    
    # Create extra subplot for map? 
    # Create an inset axes object that is 40% width and height of the
    # figure and put it in the upper right hand corner.
    h = inset_axes(skew.ax, '40%', '40%', loc=1, borderpad=3) #,

    ZHelev18 = radar.get_field(3, TH_name, copy=False)
    ZHelev18[np.where(ZHelev18<4)]=np.nan
    [lats_, lons_, _] = radar.get_gate_lat_lon_alt(3, reset_gate_coords=False, filter_transitions=False)
    h.pcolormesh(lons_, lats_, ZHelev18, cmap=P4A.colormaps('ref'), vmin=5,  vmax=65)
    h.plot(prov[:,0],prov[:,1],color='k'); 
    h.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    h.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)        
    h.set_xlim([-65.2,-63.2]); 
    h.set_ylim([-32.5,-31])
    h.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
    h.scatter(lon, lat, s=100,color='k', marker='x') 
    
    
    
    
    
    # Show the image.
    # plt.show()
    
    return fig

#------------------------------------------------------------------------------
def generate_random_latlon(center_lat, center_lon, num_points=2, max_offset_km=30):
    """
    Generate random latitude/longitude points around a central point within a given distance.
    """
    lat_offset = (np.random.uniform(-1, 1, num_points) * max_offset_km) / 111  # Approx. 1 degree ~ 111 km
    lon_offset = (np.random.uniform(-1, 1, num_points) * max_offset_km) / (111 * np.cos(np.radians(center_lat)))
   
    random_lats = center_lat + lat_offset
    random_lons = center_lon + lon_offset
   
    return list(zip(random_lats, random_lons))

#------------------------------------------------------------------------------
def is_point_within_domain(pointlat, pointlon, lat_grid, lon_grid):
    points = np.column_stack( (lat_grid.data.ravel(), lon_grid.data.ravel()) )
    hull   = Delaunay(points)    
    
    return hull.find_simplex([pointlat, pointlon]) >= 0 # check if inside

#------------------------------------------------------------------------------
def plot_sondeos(EXP, random_points, domain, title, server, sndLoc):

    folders   = config_folders.config_folders(server)
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']

    if 'yakaira' in server:
        prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')    
        fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
    elif 'cnrm' in server:
        prov = np.genfromtxt("provincias.txt", delimiter='')    
        fn = 'etopo1_bedrock.nc'  
        
    prefix   = 'wrfout_'+domain+'_2018-12-13_'+title
    filename = os.path.join(WRFfolder, 'wrfout_'+domain+'_2018-11-10_'+title+':00')
    ncfile   = Dataset(filename,'r') 
    wrflat      = wrf.getvar( ncfile,"lat") 
    wrflon      = wrf.getvar( ncfile,"lon")
    
    # Loop a bit randomly around the snd geo location
    latloniteraton = 0 

    for i in range(len(random_points)):
        
        lat = round(random_points[i][0],2)
        lon = round(random_points[i][1],2)

        # check if lat/lon within domain:
        if (is_point_within_domain(lat, lon, wrflat, wrflon) ) == 1:     

            figuretitle = EXP+'('+domain+')'+' sounding at (' + str(lat) + ', ' + str(lon) + ') '+title+' UTC'
            fig         = makesndfigure(ncfile, figuretitle, server, lat, lon)
            figurename  = EXP+domain+'_'+title+'_'+sndLoc+'_'+str(latloniteraton)
            fig.savefig(save_dir_compare+'/'+EXP+'/Sondeos/'+figurename+'.png', dpi=300,transparent=False,bbox_inches='tight')
            latloniteraton = latloniteraton + 1
            plt.close(fig)
    return

#------------------------------------------------------------------------------
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
        {"Sondeo": "CSU", "Latitude":-32.534, "Longitude":-63.257}, 
        {"Sondeo": "UI1", "Latitude":-31.877, "Longitude":-64.567}, 
        {"Sondeo": "UI2", "Latitude":-31.475, "Longitude":-64.431}, 
        {"Sondeo": "Scout1", "Latitude":-33.031, "Longitude":-63.528},      # OJO TRANSECTA
        {"Sondeo": "Scout2", "Latitude":-31.618, "Longitude":-64.876},       # OJO TRANSECTA
        {"Sondeo": "Scout3", "Latitude":-32.012, "Longitude":-64.204},       # OJO TRANSECTA
        {"Sondeo": "Villa Yacanto", "Latitude":-32.130, "Longitude":-64.730},       # OJO TRANSECTA
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
    fig, ax = plt.subplots(figsize=(8,8))      
    ax.plot(prov[:,0],prov[:,1],color='k'); 
    ax.set_xlim([-65.2,-63.2]); 
    ax.set_ylim([-32.5,-31])
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
    plt.close()
   
    return df

#------------------------------------------------------------------------------
def plot_many_wrf_sondes(dfsondeos):
    
    slatitudes      = dfsondeos['Latitude'].tolist()
    slongitudes     = dfsondeos['Longitude'].tolist()
    radiosondelabel = dfsondeos['Sondeo'].tolist()
      
    EXP = ['WSM6', 'P3mp54']

    # Loop over launched radiosonde point observations
    for i in range(len(slatitudes)):
    
        # Get 20 random lat/lon points within 30km of the radiosonde launch 
        random_points  = generate_random_latlon(slatitudes[i], slongitudes[i])

        # Loop over experiments of interest and between 18:00 and 22:00
        for EXPi in EXP:
            for h in range(18, 22):
                for m in range(0, 60, 30): 
                       print(EXPi+'(Time: '+f"{h}:{m:02d}"+') for Sondeo: '+ radiosondelabel[i])
                       plot_sondeos(EXPi, random_points, 'd02', f"{h}:{m:02d}", 'cnrm_1312', radiosondelabel[i].replace(" ", "") )

    return

#------------------------------------------------------------------------------
def plot_sondeos_compare(random_points, title, server, sndLoc):

    folders   = config_folders.config_folders(server)
    save_dir_compare = folders['save_dir_compare']

    if 'yakaira' in server:
        prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')    
        fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
    elif 'cnrm' in server:
        prov = np.genfromtxt("provincias.txt", delimiter='')    
        fn = 'etopo1_bedrock.nc'  


    # Para agregar subplot en la esquina con datos del radar
    #------------------------------------------
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,220)   
    
    # topo    
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
    #------------------------------------------
    
    EXP    = ['WSM6', 'P3mp54']
    domain = ['d02','d02']
    colorsr = ['red','darkred','magenta']
    colorsg = ['blue', 'darkblue', 'cyan']
    
    # Loop a bit randomly around the snd geo location
    latloniteraton = 0 

    for i in range(len(random_points)):
        
        lat = round(random_points[i][0],2)
        lon = round(random_points[i][1],2)
    
        figuretitle = 'Sounding at (' + str(lat) + ', ' + str(lon) + ') '+title+' UTC'

        # New figure to compare all models under evaluation 
        fig = plt.figure(figsize=(12, 9), dpi=300.)
        skew = SkewT(fig)
        skew.ax.set_xlim(-60.,40.)
        skew.ax.set_ylim(1000.,100.)
        skew.ax.tick_params(axis='both',labelsize=14)   

        skew.plot(np.nan,np.nan,'red',linewidth=2, label='WSM6')    
        skew.plot(np.nan,np.nan,'darkred',linewidth=2, label='P3')
        
        skew.plot(np.nan,np.nan,'blue',linewidth=2, label='WSM6')    
        skew.plot(np.nan,np.nan,'darkblue',linewidth=2, label='P3')
        
        plt.legend(loc='upper left')
    
        iie = 0
        pwat = []
        for EXPi in EXP:
            
            WRFfolder = folders[EXPi]        
            filename  = os.path.join(WRFfolder, 'wrfout_'+domain[iie]+'_2018-12-13_'+title+':00')
            ncfile    = Dataset(filename,'r') 
            
            wrflat      = wrf.getvar( ncfile,"lat") 
            wrflon      = wrf.getvar( ncfile,"lon")

            if (is_point_within_domain(lat, lon, wrflat, wrflon) ) == 1:     
                p, t, td, u, v = get_data_sondeo(ncfile, lat, lon)
                pwatt =  metpy.calc.precipitable_water(p,td)  # p in hPa y td y C. check! 
                pwat.append( pwatt.magnitude )


                skew.plot(p,t,colorsr[iie],linewidth=2)        # reds
                skew.plot(p,td,colorsg[iie],linewidth=2)       # greens
                #skew.plot_barbs(p,u,v)
                skew.plot_dry_adiabats()
                skew.plot_moist_adiabats()
                skew.plot_mixing_lines()
            
                del p, t, td, u, v
            else:
                pwat.append(np.nan)
            iie=iie+1

        indices = {'PW(WSM6)': ['%.1f' % pwat[0], 'mm'],\
                   'PW(P3)': ['%.1f' % pwat[1], 'mm']}
        string = ''
        for key in (indices.keys()):
            #print(key)
            string = string + key + ': ' + str(indices[key][0]) + ' ' + indices[key][1] + '\n'
        skew.ax.text(.01, .01, string, ha='left', va='top', bbox=dict(facecolor='gray', 
                            edgecolor='black', boxstyle='round'),  transform=plt.gca().transAxes)
        
        # Continue making plot    
        skew.ax.set_xlabel("Temperature (°C)", fontsize=16, weight='bold')
        skew.ax.set_ylabel("Pressure (hPa)", fontsize=16, weight='bold')
        skew.ax.set_title(figuretitle, fontsize=16, weight='bold')
            
        # Create extra subplot for map? 
        # Create an inset axes object that is 40% width and height of the
        # figure and put it in the upper right hand corner.
        h = inset_axes(skew.ax, '40%', '40%', loc=1) #,

        h.plot(prov[:,0],prov[:,1],color='k'); 
        h.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
        h.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)        
        h.set_xlim([-65.2,-63.2]); 
        h.set_ylim([-32.5,-31])
        h.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
        h.scatter(lon, lat, s=100,color='k', marker='x') 
        
        figurename = 'ModelComparison_'+title+'_'+sndLoc+'_'+str(latloniteraton)
        fig.savefig(save_dir_compare+'/Comparison/Sondeos/'+sndLoc+'/'+figurename+'.png', dpi=300,transparent=False,bbox_inches='tight')
        latloniteraton = latloniteraton + 1
        plt.close(fig)
            
    
    return fig 

#-----------------------------------------------------------------------------
def get_data_sondeo(ncfile, lat, lon):

    xloc, yloc = wrf.to_np(wrf.ll_to_xy(ncfile, lat, lon))

    p  = wrf.getvar(ncfile, "p", units="hPa")
    z  = wrf.getvar(ncfile, "z", units="m")
    t  = wrf.getvar(ncfile, "tc")
    td = wrf.getvar(ncfile, "td", units="degC")
    u  = wrf.getvar(ncfile, "uvmet", units="kt")[0,:]
    v  = wrf.getvar(ncfile, "uvmet", units="kt")[1,:]
    
    p = p.sel(south_north=yloc, west_east=xloc)
    z = z.sel(south_north=yloc, west_east=xloc)
    t = t.sel(south_north=yloc, west_east=xloc)
    td = td.sel(south_north=yloc, west_east=xloc)
    u = u.sel(south_north=yloc, west_east=xloc)
    v = v.sel(south_north=yloc, west_east=xloc)
    
    return p,t,td,u,v

#------------------------------------------------------------------------------
def plot_compare_wrf_sondes(df, server):
    
    slatitudes      = dfsondeos['Latitude'].tolist()
    slongitudes     = dfsondeos['Longitude'].tolist()
    radiosondelabel = dfsondeos['Sondeo'].tolist()
      
    # Loop over launched radiosonde point observations
    for i in range(len(slatitudes)):
    
        # Get 20 random lat/lon points within 30km of the radiosonde launch 
        random_points  = generate_random_latlon(slatitudes[i], slongitudes[i])
        
        for h in range(17, 23):
            for m in range(0, 60, 30): 
                print('(Time: '+f"{h}:{m:02d}"+') for Sondeo: '+ radiosondelabel[i])
                plot_sondeos_compare(random_points, f"{h}:{m:02d}", server, radiosondelabel[i].replace(" ", ""))
        for m in range(0, 60, 30): 
            plot_sondeos_compare(random_points, f"23:{m:02d}", server, radiosondelabel[i].replace(" ", ""))

    return


    
#------------------------------------------------------------------------------
dfsondeos = plot_domain_SurfObs('cnrm_1312')
#plot_many_wrf_sondes(dfsondeos)
plot_compare_wrf_sondes(dfsondeos, 'cnrm_1312')

