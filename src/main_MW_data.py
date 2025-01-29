#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 11:09:57 2024

@author: vito.galligani

files:
    base_dir = "/home/vito.galligani/disk1/derek/Work/HAILCASE_10112018_datos/PMW/"  # where to download all GPM data
    1C-AMSR2-GCOM1 at 1739 
    file = base_dir+'1C-AMSR2-GCOM1'+'/2018/11/10/'+'1C.GCOMW1.AMSR2.XCAL2016-V.20181110-S172223-E190115.034487.V07A.HDF5'
    1C-MHS-NOAA19 at 2033 
    file = base_dir+'1C-MHS-NOAA19'+'/2018/11/10/'+'1C.NOAA19.MHS.XCAL2021-V.20181110-S201428-E215628.050287.V07A.HDF5'
"""

import datetime
import gpm
import fsspec
import numpy as np
import ximage  # noqa
import xarray as xr
from xarray.backends.api import open_datatree
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from gpm.utils.geospatial import (
    get_country_extent,
    get_geographic_extent_around_point,
    get_circle_coordinates_around_point,
)
from gpm.gv import volume_matching, compare_maps, calibration_summary
from pathlib import Path
from gpm.utils.geospatial import get_country_extent
from gpm.utils.geospatial import get_continent_extent
import cartopy.feature as cfeature

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def create_config_file(): 
    
    username_pps = "vito.galligani@cima.fcen.uba.ar"  # likely your mail, all in lowercase
    password_pps = "Maitena10"  # likely your mail, all in lowercase
    username_earthdata = "vito.galligani@cima.fcen.uba.ar"
    password_earthdata = "PierMaitena1*"
    base_dir = "/home/vito.galligani/disk1/derek/Work/HAILCASE_10112018_datos/PMW"  # where to download all GPM data
    gpm.define_configs(
        username_pps=username_pps,
        password_pps=password_pps,
        username_earthdata=username_earthdata,
        password_earthdata=password_earthdata,
        base_dir=base_dir)
    
    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def pyplot_rings(lat_radar,lon_radar,radius):

    import numpy as np

    R=12742./2.
    m=2.*np.pi*R/360.
    alfa=np.arange(-np.pi,np.pi,0.0001)

    nazim  = 360.0
    nbins  = 480.0
    binres = 0.5

    lat_radius = lat_radar + (radius/m)*np.sin(alfa)
    lon_radius = lon_radar + ((radius/m)*np.cos(alfa)/np.cos(lat_radius*np.pi/180))

    return lat_radius, lon_radius

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def download_data(start_time, end_time): 
    
    # Download all data in the time period of interest from all PMW Level 1C products
    products = gpm.available_products(product_categories="PMW", product_levels="1C")
    
    for i in products:
        
        if i not in ('1C-AMSRE-AQUA','1C-AMSUB-NOAA15','1C-AMSUB-NOAA16',
                     '1C-AMSUB-NOAA17','1C-ATMS-NOAA21', '1C-MHS-METOPC', 
                     '1C-MHS-NOAA18', '1C-SSMI-F08', '1C-SSMI-F10', 
                     '1C-SSMI-F11', '1C-SSMI-F13', '1C-SSMI-F14', '1C-SSMI-F15',
                     '1C-SSMIS-F19', '1C-TMI'): 
            # Download the data
            gpm.download(
                product=i,
                product_type="RS",
                version=7,
                start_time=start_time,
                end_time=end_time,
                storage="GES_DISC",
                force_download=False,
                verbose=True,
                progress_bar=True,
                check_integrity=False)
        else:
            print('Not downloading '+str(i))
            
    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Specify the time period you are interested in
start_time = datetime.datetime.strptime("2018-11-10 17:00:00", "%Y-%m-%d %H:%M:%S")
end_time = datetime.datetime.strptime("2018-11-11 00:00:00", "%Y-%m-%d %H:%M:%S")

# Download data once 
do_download = 0
if do_download == 1: 
    download_data(start_time, end_time)
    
PMWdir = '/home/vito.galligani/disk1/derek/Work/HAILCASE_10112018_datos/PMW/GPM/RS/V07/PMW/'
L1Cdatasets = ['1C-AMSR2-GCOMW1', '1C-ATMS-NPP', '1C-GMI-R', '1C-MHS-METOPB', 
            '1C-SSMIS-F17', '1C-ATMS-NOAA20', '1C-GMI', '1C-MHS-METOPA', '1C-MHS-NOAA19',
            '1C-SSMIS-F16', '1C-SSMIS-F18'] 

# Add limit provincias 
states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none',
        edgecolor='black')

# Info RMA1
radarLAT = -31.441389
radarLON = -64.191944

# Save figures here: 
fig_dir = '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/PMW_obs/'
    
for product in L1Cdatasets:

    print('Running plots for: '+ product)
    # Load the 1C dataset (get only Level 1C Tb)
    # - If scan_mode is not specified, it automatically load one!
    scanmodes = gpm.available_scan_modes(product=product, version=7) 
    
    if product in ('1C-MHS-METOPB','1C-SAPHIR-MT1', '1C-MHS-NOAA19', '1C-MHS-METOPA'):
        scan_mode = scanmodes[0] 
    elif product in ('1C-GMI-R', '1C-GMI'):
        scan_mode = scanmodes[1]         
    else:
        scan_mode = scanmodes[-2] 
    
    ds = gpm.open_dataset(product=product, product_type="RS", version=7,  
            start_time=start_time, end_time=end_time, variables='Tc', scan_mode=scan_mode)        
    
    # All orbits and crop by country
    #title = ds.gpm.title(add_timestep=True)
    #extent = get_country_extent("Argentina")
    #da = ds['Tc'].isel(pmw_frequency=0)
    #p = da.gpm.plot_map(vmin=150, vmax=300)
    #_ = p.axes.set_extent(extent)
    #_ = p.axes.set_title(label=title)

    # Crop by continent L1C Tc
    da = ds['Tc'].isel(pmw_frequency=0)
    extent = get_continent_extent("South America")
    extent_plot = (-70, -50, -40, -20)  # (xmin, xmax, ymin, ymax)
    list_isel_dict = da.gpm.get_crop_slices_by_extent(extent)
    # - Plot the swath crossing the country
    for isel_dict in list_isel_dict:
        da_subset = da.isel(isel_dict)
        slice_title = da_subset.gpm.title(add_timestep=True)
        p = da_subset.gpm.plot_map(vmin=150, vmax=300)
        p.axes.set_extent(extent_plot)
        p.axes.set_title(label=slice_title + ' ' + ds.coords['pmw_frequency'].data[0] )
        p.axes.add_feature(states_provinces,linewidth=0.4)
        [lat_radius, lon_radius] = pyplot_rings(radarLAT,radarLON,240)
        p.axes.plot(lon_radius, lat_radius, 'k', linewidth=0.8)  
    
        figfile = (slice_title + ' ' + ds.coords['pmw_frequency'].data[0]).replace(' ', '') 
        plt.savefig(fig_dir+figfile+'.png', dpi=300,transparent=False)    
        
        
        





    
    

    
    
    


    
        