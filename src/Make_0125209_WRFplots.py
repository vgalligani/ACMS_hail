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
import config_folders_final
from Plots4Analysis  import pyplot_rings
from matplotlib.patches import Polygon
import netCDF4 as nc
import matplotlib.colors as mcolors
import package_functions
from matplotlib.cm import get_cmap
from PIL import Image

import warnings
warnings.filterwarnings("ignore")

plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  

lonrange = [-67.5,-62]
latrange = [-35,-31]
prov     = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_ZH1km_WRF(EXP, title, domain, date, servidor, folders):

    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']

    if 'yakaira' in servidor:
        prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')
        fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
    
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
    
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,220)   
    
    prefix   = f'wrfout_{domain}_{date}_{title}'
    filename = os.path.join(WRFfolder, prefix+':00')
    
    ncfile       = Dataset(filename,'r')        
    z            = wrf.getvar( ncfile,"z") 
    zh           = wrf.getvar(ncfile, "REFL_10CM")
    REFL_10CM    = wrf.interplevel(zh, z, 3000)
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")

       
    fig, ax = plt.subplots(figsize=(8,8)) 
    pcm = ax.pcolormesh(lon, lat,  REFL_10CM, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
    cbar = plt.colorbar(pcm, ax=ax, shrink=1)
    cbar.cmap.set_under('white')
        
    ax.set_xlim(lonrange) 
    ax.set_ylim(latrange)
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
    # agrego contorno de 500 y 1000m y provincias 
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
        
    ax.grid()
    ax.set_title('Zh interp. at 3km at '+title)
    
    fig.savefig(save_dir_compare+'/'+EXP+f'/WRF_ZH3km_{domain}_{title}.png', dpi=300,
                transparent=False, bbox_inches='tight')
    plt.close()
    
    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def single_WRF_files(EXP, domain, title, save_dir_compare): 
    
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']
    
    fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
    
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,220)   
    
    start_time_1    = 1200
    
    all_files_ = sorted(glob.glob(os.path.join(WRFfolder, 'wrfout_d02_2019-01-25*')))

    prefix = f'wrfout_{domain}_2019-01-25_'
    start_index = all_files_[0].find(prefix) + len(prefix)
            
    # filter and keep only the hh:00 or hh:30 files
    all_files = [f for f in all_files_ if f[-5:-3] in ("00", "30")]
        
    # Filter filenames based on HH:MM > '12:00' and < 22:00 HARDCODED
    filtered_files = [
        file for file in all_files
        if start_time_1-10 < int(os.path.basename(file).split('_')[3].split(':')[0]) * 100 +
           int(os.path.basename(file).split('_')[3].split(':')[1]) < 2200 ]
        

    counter=0
    rain_prev=0
    for filename in filtered_files:
        
        times = filename[start_index:start_index+5]   

        ncfile       = Dataset(filename,'r')        
        rain         = wrf.getvar(ncfile, "RAINNC")
        if counter >0:
            rain = rain-rain_prev
        rain_prev    = rain.copy()
        lat          = wrf.getvar( ncfile,"lat") 
        lon          = wrf.getvar( ncfile,"lon")

        fig, ax = plt.subplots(figsize=(8,8)) 
        pcm = ax.pcolormesh(lon, lat,  rain, cmap=P4A.colormaps('ref'))
        cbar = plt.colorbar(pcm, ax=ax, shrink=1)


        ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
        ax.plot(prov[:,0],prov[:,1],color='k'); 
        ax.set_xlim(lonrange); 
        ax.set_ylim(latrange)
        ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
        ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)

        ax.set_title(title+' '+str(times))
        fig.savefig(save_dir_compare+'/WRF_RAINC_'+title+'_time'+str(times)+'.png', dpi=300,transparent=False, bbox_inches='tight')
        plt.close()
        counter=counter+1
        
        
    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
folders=config_folders_final.config_folders('yakaira')
    
EXP    = '2501_WSM6'
date   = '2019-01-25' 
domain = 'd02'

for h in range(16, 23):
    for m in range(0, 60, 30):
        plot_ZH1km_WRF(EXP, f"{h}:{m:02d}", domain, date, 'yakaira', folders)
plot_ZH1km_WRF(EXP, "23:00", domain, date, 'yakaira', folders)
plot_ZH1km_WRF(EXP, "23:30", domain, date, 'yakaira', folders)
plot_ZH1km_WRF(EXP, "00:00", domain, date, 'yakaira', folders)
plot_ZH1km_WRF(EXP, "00:30", domain, date, 'yakaira', folders)
plot_ZH1km_WRF(EXP, "01:00", domain, date, 'yakaira', folders)
plot_ZH1km_WRF(EXP, "01:30", domain, date, 'yakaira', folders)

single_WRF_files('d02', '2501_WSM6', folders)
                

# Find all matching PNG files and sort them (optional but recommended)
save_dir_compare = folders[EXP]['save_dir_compare']    
image_files = sorted(glob.glob(save_dir_compare+'/'+EXP+"/WRF_ZH3km_*.png"))

# Open all images and convert them to RGB (PDF doesn't support RGBA)
images = [Image.open(f).convert("RGB") for f in image_files]

# Save them as a single PDF
if images:
    images[0].save("output.pdf", save_all=True, append_images=images[1:])
else:
    print("No images found matching 'Zh*.png'")
