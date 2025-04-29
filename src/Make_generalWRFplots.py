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

import warnings
warnings.filterwarnings("ignore")

plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  


#------------------------------------------------------------------------------
def plot_radar_singletime(rfiles, time, elev, latrange, lonrange, colmax, folders):

    prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')

    # RMA1 
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    TH_name       = 'TH'
    time1infile  = 82
    
    #folders          = config_folders.config_folders('yakaira')
    radar_folder     = folders['rma1_dir']
    save_dir_compare = folders['save_dir_compare']
   
    filename = os.path.join(radar_folder, rfiles)
    radar    = pyart.io.read(filename) 

    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,220) 
    
    # Configure a gatefilter to filter out copolar correlation coefficient values > 0.9
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_equal('RHOHV', 0.9)
    
    ZHelev18 = radar.get_field(elev, TH_name, copy=False)
    [lats_, lons_, _] = radar.get_gate_lat_lon_alt(elev, reset_gate_coords=False, filter_transitions=False)

    fig, ax = plt.subplots(figsize=(8,8)) 
    if colmax == 1:
        gatefilter = pyart.filters.GateFilter(radar)
        gatefilter.exclude_transition()
        gatefilter.exclude_equal('RHOHV', 0.9)
        lons_,lats_,ZHcolmax = package_functions.get_colmax(radar, TH_name, gatefilter)
        pcm = ax.pcolormesh(lons_, lats_, ZHcolmax, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
        ax.set_title('RMA1 zoom obs COLMAX at ' + time)
        title = time+'_colmax'    
    else:
        pcm = ax.pcolormesh(lons_, lats_, ZHelev18, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
        ax.set_title('RMA1 zoom obs elev '+str(elev)+ ' at ' + time)
        title = time+'_elev_'+str(elev)
    
    cbar = plt.colorbar(pcm, ax=ax, shrink=1, label='Zh RMA1 elev 3')
    cbar.cmap.set_under('white')
    ax.grid()       
    ax.plot(prov[:,0],prov[:,1],color='k'); 
    ax.set_xlim(lonrange); 
    ax.set_ylim(latrange)
    

    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
    
    plt.show()
    fig.savefig(save_dir_compare+'/OBS/RMA1'+'/ZH_RMA1_obs_'+title+'.png', dpi=300,transparent=False, bbox_inches='tight')
    
    return


#------------------------------------------------------------------------------
def plot_radar_cspr2_singletime(latrange, lonrange):

    prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')

    folders          = config_folders_final.config_folders('yakaira')
    radar_folder     = folders['csapr2_dir']
    save_dir_compare = folders['save_dir_compare']
   
    file_list    = sorted(glob.glob(radar_folder+'*b1*.nc'))
    
    # RMA1 
    radarLAT = -32.12641
    radarLON = -64.72837
    TH_name  = 'attenuation_corrected_reflectivity_h'
    elev     = 1 #(==1.4996338)
    
    [lat_radius, lon_radius] = pyplot_rings(radarLAT,radarLON,30)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT,radarLON,60) 
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT,radarLON,100) 
    
    prefix = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/CSAPR2/corcsapr2cfrppiqcM1.b1.20181110'

    for filename in file_list:
        radar    = pyart.io.read(filename) 
        ZHelev18 = radar.get_field(elev, TH_name, copy=False)
        [lats_, lons_, _] = radar.get_gate_lat_lon_alt(elev, reset_gate_coords=False, filter_transitions=False)

        start_index = len(prefix)
        time        = filename[start_index+1:start_index+5]

        fig, ax = plt.subplots(figsize=(8,8)) 
        pcm = ax.pcolormesh(lons_, lats_, ZHelev18, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
        cbar = plt.colorbar(pcm, ax=ax, shrink=1, label=r'CSAPR2 Zh Att. corr. Zh (1.5$^o$)')
        cbar.cmap.set_under('white')
        ax.grid()       
        ax.plot(prov[:,0],prov[:,1],color='k'); 
        ax.set_xlim(lonrange); 
        ax.set_ylim(latrange)
        ax.set_title('CSAPR2 zoom obs elev '+str(elev)+ 'at ' + time)
        
        ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
        ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
        plt.show()
        fig.savefig(save_dir_compare+'/OBS'+'/CSAPR2/ZH_CSAPR2_obs_'+time+'.png', dpi=300,transparent=False, bbox_inches='tight')
        plt.close()
 
    return
#------------------------------------------------------------------------------
def discrete_cmapp(N, base_cmap=None):

    base  = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0,1,N))
    cmap_name = base.name+str(N)


    return base.from_list(cmap_name, color_list, N)

#------------------------------------------------------------------------------
def plot_VELradar_cspr2_singletime(latrange, lonrange):


    prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')

    folders          = config_folders_final.config_folders('yakaira')
    radar_folder     = folders['csapr2_dir']
    save_dir_compare = folders['save_dir_compare']
   
    file_list    = sorted(glob.glob(radar_folder+'*b1*.nc'))
    
    # RMA1 
    radarLAT = -32.12641
    radarLON = -64.72837
    TH_name  = 'attenuation_corrected_reflectivity_h'
    elev     = 1 #(==1.4996338)
    
    [lat_radius, lon_radius] = pyplot_rings(radarLAT,radarLON,30)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT,radarLON,60) 
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT,radarLON,100) 
    
    prefix    = 'corcsapr2cfrppiqcM1.b1.20181110.'
    folder    = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/CSAPR2/'
       
    #------ 1745
    file_1745 = 'corcsapr2cfrppiqcM1.b1.20181110.174503.nc' 
    radar     = pyart.io.read(folder+file_1745) 
    VEL       = radar.get_field(0, 'mean_doppler_velocity', copy=False)
    [lats_, lons_, _] = radar.get_gate_lat_lon_alt(0, reset_gate_coords=False, filter_transitions=False)
    start_index = len(prefix)
    time        = file_1745[start_index:start_index+5]
    
    fig, ax = plt.subplots(figsize=(8,8)) 
    pcm = ax.pcolormesh(lons_, lats_, VEL, cmap=discrete_cmapp(12,'seismic'), vmin=-12,  vmax=12)
    cbar = plt.colorbar(pcm, ax=ax, shrink=1, label=r'CSAPR2 vel (0.5$^o$)', ticks=np.arange(-12,12.01,2))
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    ax.grid()       
    ax.plot(prov[:,0],prov[:,1],color='k'); 
    ax.set_xlim(lonrange); 
    ax.set_ylim(latrange)
    ax.set_title('CSAPR2 zoom obs elev '+str(elev)+ 'at ' +time)
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
    
    plt.show()
    fig.savefig(save_dir_compare+'OBS'+'/CSAPR2/ZH_CSAPR2_obs_'+time+'.png', dpi=300,transparent=False, bbox_inches='tight')
    plt.close()

    #------ 1730
    file_1730 = 'corcsapr2cfrppiqcM1.b1.20181110.173003.nc' 
    radar     = pyart.io.read(folder+file_1730) 
    VEL       = radar.get_field(0, 'mean_doppler_velocity', copy=False)
    [lats_, lons_, _] = radar.get_gate_lat_lon_alt(0, reset_gate_coords=False, filter_transitions=False)
    start_index = len(prefix)
    time        = file_1730[start_index:start_index+5]
    fig, ax = plt.subplots(figsize=(8,8)) 
    pcm = ax.pcolormesh(lons_, lats_, VEL, cmap=discrete_cmapp(12,'seismic'), vmin=-12,  vmax=12)
    cbar = plt.colorbar(pcm, ax=ax, shrink=1, label=r'CSAPR2 vel (0.5$^o$)', ticks=np.arange(-12,12.01,2))
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    ax.grid()       
    ax.plot(prov[:,0],prov[:,1],color='k'); 
    ax.set_xlim(lonrange); 
    ax.set_ylim(latrange)
    ax.set_title('CSAPR2 zoom obs elev '+str(elev)+ 'at ' + time)
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
    
    plt.show()
    fig.savefig(save_dir_compare+'OBS'+'/CSAPR2/ZH_CSAPR2_obs_'+time+'.png', dpi=300,transparent=False, bbox_inches='tight')
    plt.close()


    
    #------
 
    return


#------------------------------------------------------------------------------
def plot_ZH1km_WRF_wWRF_shelicity(EXP, title, folders, domain):
    
    import matplotlib

    #folders   = config_folders.config_folders('yakaira')
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']

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
    
    
    prefix   = 'wrfout_'+domain+'_2018-11-10_'+title
    filename = os.path.join(WRFfolder, 'wrfout_'+domain+'_2018-11-10_'+title+':00')

    ncfile       = Dataset(filename,'r')        
    #------ READ WRF variables of interest 
    z            = wrf.getvar( ncfile,"z") 
    zh           = wrf.getvar(ncfile, "REFL_10CM")
    COLMAXREFL   = np.nanmax(zh, axis=0)
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")
    
    # ADEMAS SUMAR VIENTOS
    #----------    
    uvmet10 = wrf.getvar(ncfile, "uvmet10",units="kt")   # Default is ‘m s-1’.
    wspd_wdir10    = wrf.getvar(ncfile, "wspd_wdir10") 
    
    # Storm helicity
    #----------
    shelicity = -1.*(wrf.getvar(ncfile, 'helicity',top=1000.0))

    fig, ax = plt.subplots(figsize=(8,8)) 
    ax.plot(np.nan, np.nan, 'r', linewidth=1.5, label='COLMAX 35dBz')
    ax.plot(np.nan, np.nan, 'gray', linewidth=2.5, label='Topography')
    plt.legend()
    
    helicity_contour = ax.contourf(wrf.to_np(lon), wrf.to_np(lat), wrf.to_np(shelicity), np.arange(50,750,50), cmap=matplotlib.cm.get_cmap("gist_stern_r"))
    plt.colorbar(helicity_contour, ax=ax, orientation="horizontal", pad=.05, ticks=np.arange(50,750,50))
    ax.contour(lon, lat, COLMAXREFL, levels=[35], colors=['r'], linewidths=1.5)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
        
    if 'Maite' in title: 
        ax.set_xlim([-65,-63.9]); 
        ax.set_ylim([-33,-31.5])
                
    else:
        ax.set_xlim([-65.5,-62]); 
        ax.set_ylim([-35,-31])
        
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
    # agrego contorno de 500 y 1000m
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
        
    # wind barbs    
    resobarb = 15
    ax.barbs(wrf.to_np(lon[::resobarb,::resobarb]), wrf.to_np(lat[::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[0,::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[1,::resobarb,::resobarb]), length=6)
    
    
    #cbar.cmap.set_under('white')
    ax.grid()
    ax.set_title('0-1km AGL storm relative helicity (m2 s-2) at '+title)
    
    ax.text(x=-65, y=-34.9, s='120 and 220 km radar rings')
    #plt.show()
    fig.savefig(save_dir_compare+'/'+EXP+'/helicity/WRF_ZH1km_general_evolution_'+title+'stromhelicity.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()



    # Updraft Helicity
    #----------
    uhelicity = -1.*wrf.getvar(ncfile, 'updraft_helicity')

    fig, ax = plt.subplots(figsize=(8,8)) 
    ax.plot(np.nan, np.nan, 'r', linewidth=1.5, label='COLMAX 35dBz')
    ax.plot(np.nan, np.nan, 'gray', linewidth=2.5, label='Topography')
    plt.legend()
    
        
    levels = np.arange(0,20,2) 
    
    helicity_contour = ax.contourf(wrf.to_np(lon), wrf.to_np(lat), wrf.to_np(uhelicity), 10, cmap=matplotlib.cm.get_cmap("gist_stern_r"), levels=levels, extend='both', alpha=0.8)
    plt.colorbar(helicity_contour, ax=ax, orientation="horizontal", pad=.05, ticks=levels)
    
    ax.contour(lon, lat, COLMAXREFL, levels=[35], colors=['r'], linewidths=1.5)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
        
    if 'Maite' in title: 
        ax.set_xlim([-65,-63.9]); 
        ax.set_ylim([-33,-31.5])
                
    else:
        ax.set_xlim([-65.5,-62]); 
        ax.set_ylim([-35,-31])
        
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
    # agrego contorno de 500 y 1000m
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
        
    # wind barbs    
    resobarb = 15
    ax.barbs(wrf.to_np(lon[::resobarb,::resobarb]), wrf.to_np(lat[::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[0,::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[1,::resobarb,::resobarb]), length=6)
    
    
    #cbar.cmap.set_under('white')
    ax.grid()
    ax.set_title('Updraft helicity (m2 s-2) at '+title)
    
    ax.text(x=-65, y=-34.9, s='120 and 220 km radar rings')
    #plt.show()
    fig.savefig(save_dir_compare+'/'+EXP+'/helicity/WRF_ZH1km_general_evolution_'+title+'updrafthelicity.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()


    return
    
#------------------------------------------------------------------------------
def plot_ZH1km_WRF_wWRF_uhelicity_only(EXP, intplev, title, folders, domain):
    
    import matplotlib

    #folders   = config_folders.config_folders('yakaira')
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']

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
    
    prefix   = 'wrfout_'+domain+'_2018-11-10_'+title
    filename = os.path.join(WRFfolder, prefix+':00')
    ncfile       = Dataset(filename,'r')        
    
    #------ READ WRF variables of interest ------------------------------------
    z            = wrf.getvar( ncfile,"z") 
    zh           = wrf.getvar(ncfile, "REFL_10CM")
    COLMAXREFL   = np.nanmax(zh, axis=0)
    zh_1km      = wrf.interplevel(zh, z, 1000)
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")
    uvmet10 = wrf.getvar(ncfile, "uvmet10",units="kt")   # Default is ‘m s-1’.
    wspd_wdir10    = wrf.getvar(ncfile, "wspd_wdir10") 
    uhelicity = -1.*wrf.getvar(ncfile, 'updraft_helicity')

    #--------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8,8)) 
    
    ax.plot(np.nan, np.nan, 'r', linewidth=1.5, label='COLMAX 45dBz')
    ax.plot(np.nan, np.nan, 'gray', linewidth=2.5, label='Topography')

    if 'Maite' in title: 
        ax.set_xlim([-65,-63.9]); 
        ax.set_ylim([-33,-31.5])
    else:
        #ax.set_xlim([-65.5,-63]); 
        #ax.set_ylim([-32.5,-31.4])
        ax.set_xlim([-65.5,-62]); 
        ax.set_ylim([-32.5,-31])
        
    plt.legend(loc='lower right')
     
    levels = np.arange(0,20,2) 
    # helicity_contour = ax.contourf(wrf.to_np(lon), wrf.to_np(lat), wrf.to_np(uhelicity), 10, cmap=matplotlib.cm.get_cmap("gist_stern_r"), levels=levels, extend='both', alpha=0.8)
    # plt.colorbar(helicity_contour, ax=ax, orientation="horizontal", pad=.05, ticks=levels)
    
    helicity_contour = ax.pcolormesh(wrf.to_np(lon), wrf.to_np(lat), wrf.to_np(uhelicity), cmap=matplotlib.cm.get_cmap("gist_stern_r"), vmin=-10, vmax=80)
    plt.colorbar(helicity_contour, ax=ax)
    
    ax.contour(lon, lat, COLMAXREFL, levels=[45], colors=['r'], linewidths=1.5)
    ax.plot(prov[:,0],prov[:,1],color='k');     

        
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
    # agrego contorno de 500 y 1000m
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
        
    # wind barbs    
    resobarb = 15
    ax.barbs(wrf.to_np(lon[::resobarb,::resobarb]), wrf.to_np(lat[::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[0,::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[1,::resobarb,::resobarb]), length=6)
    
    
    ax.grid()
    ax.set_title('Updraft helicity (m2 s-2) at '+title)    
    plt.show()

    fig.savefig(save_dir_compare+'/'+EXP+'/helicity/WRF_ZH1km_general_evolution_'+title+'updrafthelicity_zoom.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()


    return#------------------------------------------------------------------------------
    
    
#------------------------------------------------------------------------------
def plot_WRF_intqx(EXP, title, mp, folders, domain):
    
    import matplotlib

    #folders   = config_folders.config_folders('yakaira')
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']

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
    
    prefix   = 'wrfout_'+domain+'_2018-11-10_'+title
    filename = os.path.join(WRFfolder, 'wrfout_'+domain+'_2018-11-10_'+title+':00')

    ncfile       = Dataset(filename,'r')        
    #------ READ WRF variables of interest 
    z            = wrf.getvar( ncfile,"z") 
    #zh           = wrf.getvar(ncfile, "REFL_10CM")
    #COLMAXREFL   = np.nanmax(zh, axis=0)
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")
    
    #-----
    # Need to interpolate first to common z_level grid 
    z_interp = 1e3*np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 
                 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20])
    Re       = 6.3781e6
    
    temp     = wrf.g_temp.get_tk(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL)
    z        = Re*geopo_p/(Re-geopo_p)

    if mp == 6: 
        [qr, qs, qi, qc, qg, qr_int, qs_int, qi_int, qc_int, qg_int] = funs.get_q_ints6(ncfile)
        int_titles = ['qr','qc','qi','qs','qg']
        subplots_tot = 5
    elif mp == 54:
        [qr_P3_3mom_LF, qi_P3_3mom_LF, qc_P3_3mom_LF, 
        qr_int, qi_int, qc_int, qir, qib, qil] = funs.get_q_ints3(ncfile)
        int_titles = ['qr','qc','qi']
        subplots_tot = 3 

            
    vminn = 0
    vmaxx = 12
    cmap  = get_cmap("viridis")
    vminn2 = 0
    vmaxx2 = 4        
    
    
    fig, ax = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=(16,(16/5))) 

    for i in range(subplots_tot):
        ax[i].plot(np.nan, np.nan, 'gray', linewidth=2.5, label='Topography')
        ax[i].grid()

        if i == 0:
            ax[i].legend()
    
    pcm0 = ax[0].pcolormesh(lon, lat, qr_int, cmap=cmap, vmin=vminn,  vmax=vmaxx)
    pcm1 = ax[1].pcolormesh(lon, lat, qc_int, cmap=cmap, vmin=vminn2,  vmax=vmaxx2)
    
    if mp == 54:
        pcm2 = ax[2].pcolormesh(lon, lat, qi_int, cmap=cmap, vmin=vminn2,  vmax=20)
    elif mp == 6:
        pcm2 = ax[2].pcolormesh(lon, lat, qi_int, cmap=cmap, vmin=vminn2,  vmax=vmaxx2)

        
    if mp == 6:    
        pcm3 = ax[3].pcolormesh(lon, lat, qs_int, cmap=cmap, vmin=vminn,  vmax=vmaxx)
        pcm4 = ax[4].pcolormesh(lon, lat, qg_int, cmap=cmap, vmin=vminn,  vmax=vmaxx)
        
    for i in range(subplots_tot): 
        ax[i].plot(prov[:,0],prov[:,1],color='k'); 
        
        if 'Maite' in title: 
            ax[i].set_xlim([-65,-63.9]); 
            ax[i].set_ylim([-33,-31.5])
                    
        else:
            ax[i].set_xlim([-67.5,-62]); 
            ax[i].set_ylim([-35,-31])
            #ax[i].set_xlim([-65.5,-62]); 
            #ax[i].set_ylim([-35,-31])
            
        ax[i].plot(lon_radius, lat_radius, 'k', linewidth=0.8)
        ax[i].plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
        # agrego contorno de 500 y 1000m
        ax[i].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
        
    for i in range(subplots_tot): 
        cbar = plt.colorbar(pcm0, ax=ax[i], shrink=1)
        cbar.cmap.set_under('white')
        
    for i in range(subplots_tot):
        ax[i].set_title('Integrated' +int_titles[i]+' '+title)
    
    if mp == 54:
        fig.delaxes(ax[3])
        fig.delaxes(ax[4])
        
    
    plt.show()
    
    fig.savefig(save_dir_compare+'/'+EXP+'/WRF_totqx_'+title, dpi=300,transparent=False,bbox_inches='tight')
    plt.close()    
    
    
    
    
    return


    
    
#------------------------------------------------------------------------------
def plot_ZH1km_WRF_wWRFwind(EXP, title, domain, intplev, folders):
    
    import matplotlib

    #folders   = config_folders.config_folders('yakaira')
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']


    print(save_dir_compare)

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
    
    prefix   = 'wrfout_'+domain+'_2018-11-10_'+title
    filename = os.path.join(WRFfolder, 'wrfout_'+domain+'_2018-11-10_'+title+':00')

    ncfile       = Dataset(filename,'r')        
    #------ READ WRF variables of interest 
    z            = wrf.getvar( ncfile,"z") 
    zh           = wrf.getvar(ncfile, "REFL_10CM")
    COLMAXREFL   = np.nanmax(zh, axis=0)
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")
    
    # ADEMAS SUMAR VIENTOS
    #----------    
    uvmet10 = wrf.getvar(ncfile, "uvmet10",units="kt")   # Default is ‘m s-1’.
    wspd_wdir10    = wrf.getvar(ncfile, "wspd_wdir10") 
    
    # MOIST CONVERGENCE
    #----------
    uvmet = wrf.getvar(ncfile, 'wspd_wdir', units="kt")
    wspd  = uvmet.isel(wspd_wdir=0)
    eth   = wrf.getvar(ncfile, 'eth')
    q     = wrf.getvar(ncfile, 'QVAPOR')
    dx=uvmet.XLAT[1,0]-uvmet.XLAT[0,0]
    dx=dx.values*111100

    #    MFC_advect = -( u*(dq/dx)+v*(dq/dy) )       ; advection term
    #    MFC_conv    = -q*( (du/dx)+  (dv/dy) )      ; con(div)-vergence
    #    MFC = MFC_advect + MFC_conv

    p = wrf.getvar(ncfile, 'p', units="hPa")
    u = wrf.getvar(ncfile, 'ua')
    v = wrf.getvar(ncfile, 'va')
    ufld = wrf.interplevel(u, p, intplev)
    vfld = wrf.interplevel(v, p, intplev)
    qfld = wrf.interplevel(q, p, intplev)
    
    grad_q_x,grad_q_y = np.gradient(qfld.values)
    grad_u_x,grad_u_y = np.gradient(ufld.values)
    grad_v_x,grad_v_y = np.gradient(vfld.values)
    
    MFC_advect=-1.* (ufld.values*grad_q_x/dx)+(vfld.values*grad_q_y/dx)
    MFC_conv=-1.*qfld.values*((grad_u_x/dx)+(grad_v_y/dx))
    
    #-------- for plots 
    cffield = qfld
    cffield.values= -1.*wrf.smooth2d(86400.*1000.*(MFC_advect + MFC_conv),30)

    cffield.attrs['description']=str(intplev)+' hPa moisture convergence and theta-e'
    cffield.attrs['units']='g kg-1 dy-1; K'
    ethinterp = wrf.smooth2d(wrf.interplevel(eth, p, intplev),30)
    #lfield= None
    lfield3 = None 
    ufld = ufld*1.94
    vfld = vfld*1.94

    fig, ax = plt.subplots(figsize=(8,8)) 
    ax.plot(np.nan, np.nan, 'r', linewidth=1.5, label='COLMAX 35dBz')
    ax.plot(np.nan, np.nan, 'gray', linewidth=2.5, label='Topography')
    plt.legend()
    
    wspd_contours = ax.contourf(wrf.to_np(lon), wrf.to_np(lat), wrf.to_np(wspd_wdir10[0,:,:]), np.linspace(4.0, 24.0, 6),
                             cmap=matplotlib.cm.get_cmap("YlOrBr"))
    plt.colorbar(wspd_contours, ax=ax, orientation="horizontal", pad=.05, ticks=np.arange(0,36,4))

    ax.contour(lon, lat, COLMAXREFL, levels=[35], colors=['r'], linewidths=1.5)
    
    ax.plot(prov[:,0],prov[:,1],color='k'); 
        
    if 'Maite' in title: 
        ax.set_xlim([-65,-63.9]); 
        ax.set_ylim([-33,-31.5])
                
    else:
        ax.set_xlim([-65.5,-62]); 
        ax.set_ylim([-35,-31])
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
    # agrego contorno de 500 y 1000m
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
        
    # wind barbs    
    resobarb = 15
    ax.barbs(wrf.to_np(lon[::resobarb,::resobarb]), wrf.to_np(lat[::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[0,::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[1,::resobarb,::resobarb]), length=6)
    
    
    #cbar.cmap.set_under('white')
    ax.grid()
    ax.set_title('10m wind (m/s) at '+title)
    
    ax.text(x=-65, y=-34.9, s='120 and 220 km radar rings')
    #plt.show()
    fig.savefig(save_dir_compare+'/'+EXP+'/winds/WRF_ZH1km_'+domain+'_general_evolution_'+title+'10mwind.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()
    
    #----------------------------------------------------
    # and moisture convernce ... 
    
    fig, ax = plt.subplots(figsize=(8,8)) 
    
    ax.plot(np.nan, np.nan, 'k', linewidth=2, label='COLMAX 35dBz')
    ax.plot(np.nan, np.nan, 'gray', linewidth=2.5, label='Topography')
    plt.legend()
    
    #---
    ax.contour(wrf.to_np(lon), wrf.to_np(lat), wrf.to_np(ethinterp), 10, colors="navy", levels=np.arange(250,400,5),linewidths=1.2)

    resobarb = 15
    ax.barbs(wrf.to_np(lon[::resobarb,::resobarb]), wrf.to_np(lat[::resobarb,::resobarb]), 
                 wrf.to_np(ufld[::resobarb,::resobarb]),
                 wrf.to_np(vfld[::resobarb,::resobarb]), length=5, linewidth=0.75, zorder=10)

    colors1 = plt.cm.YlOrRd_r(np.linspace(0, 0.8, 120))    
    colors2 = plt.cm.PiYG(np.linspace(0.8, 0.7, 60))
    colors = np.vstack((colors1,colors2))#, colors2, colors1))
    cmap_conv= mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    cmap_conv.set_over('white')
    cmap_conv.set_under('red')

    moistlevels = np.arange(-600,20,100) #np.arange(-500,-20,100)  #np.arange(-220,0,40) 
    # moistlevels = np.arange(-180,20,40)
    
    fcontours = ax.contourf(wrf.to_np(lon), wrf.to_np(lat), wrf.to_np(cffield), 10, cmap=cmap_conv, levels=moistlevels, extend='both',
               alpha=0.8)
    plt.colorbar(fcontours, ax=ax, orientation="horizontal", pad=.05,  extend='both', ticks=moistlevels)

    #---- 
    ax.plot(prov[:,0],prov[:,1],color='k'); 
    ax.contour(lon, lat, COLMAXREFL, levels=[35], colors=['k'], linewidths=2)
        
    if 'Maite' in title: 
        ax.set_xlim([-65,-63.9]); 
        ax.set_ylim([-33,-31.5])
                
    else:
        ax.set_xlim([-65.5,-62]); 
        ax.set_ylim([-35,-31])
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
    # agrego contorno de 500 y 1000m
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
         
    ax.grid()
    ax.set_title(str(intplev)+' hPa moisture convergence and theta-e at '+title)
    
    ax.text(x=-65, y=-34.9, s='120 and 220 km radar rings')
    #plt.show()
    fig.savefig(save_dir_compare+'/'+EXP+'/convergence/WRF_ZH1km_'+domain+'_general_evolution_'+title+'WRFvars_moistureconvergence'+str(intplev)+'_simpler.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()

    
    
    
    return#------------------------------------------------------------------------------
    
    
#------------------------------------------------------------------------------
def plot_ZH1km_WRF(EXP, title, domain, servidor, folders):

    #folders   = config_folders.config_folders(servidor)
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']

    if 'yakaira' in servidor:
        prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')
        fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'

    elif 'cnrm' in servidor:
        prov = np.genfromtxt("provincias.txt", delimiter='')
        fn = 'etopo1_bedrock.nc'
    
    
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
    
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,220)   
    
    #prefix   = 'wrfout_'+domain+'_2019-01-25_'+title
    prefix   = 'wrfout_'+domain+'_2018-11-10_'+title
    #prefix   = 'wrfout_'+domain+'_2018-11-04_'+title

    filename = os.path.join(WRFfolder, prefix+':00')

    
    ncfile       = Dataset(filename,'r')        
    z            = wrf.getvar( ncfile,"z") 
    zh           = wrf.getvar(ncfile, "REFL_10CM")
    REFL_10CM    = wrf.interplevel(zh, z, 3000)
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")

    # ADEMAS SUMAR VIENTOS
    #----------    
    uvmet10 = wrf.getvar(ncfile, "uvmet10")   # Default is ‘m s-1’.
        
    fig, ax = plt.subplots(figsize=(8,8)) 
    pcm = ax.pcolormesh(lon, lat,  REFL_10CM, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
    cbar = plt.colorbar(pcm, ax=ax, shrink=1)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
        
    if 'Maite' in title: 
        ax.set_xlim([-65,-63.9]); 
        ax.set_ylim([-33,-31.5])
                
    else:
        ax.set_xlim([-67.5,-62]); 
        ax.set_ylim([-35,-31])
        
        #ax.set_xlim([-65.5,-62]); 
        #ax.set_ylim([-32.5,-31])

        #ax.set_xlim([-65.5,-62]); orig. 
        #ax.set_ylim([-35,-31])
        
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
    # agrego contorno de 500 y 1000m
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
        
    # wind barbs    
    resobarb = 20
    ax.barbs(wrf.to_np(lon[::resobarb,::resobarb]), wrf.to_np(lat[::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[0,::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[1,::resobarb,::resobarb]), length=6)
    

    cbar.cmap.set_under('white')
    ax.grid()
    ax.set_title('Zh interp. at 3km at '+title)
    
    #ax.text(x=-65, y=-34.9, s='120 and 220 km radar rings')
    #plt.show()
    fig.savefig(save_dir_compare+'/'+EXP+'/WRF_ZH3km_'+domain+'_ZOOM_general_evolution_'+title+'.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()
    
    return



    
#------------------------------------------------------------------------------
def plot_ZH1km_WRFdate(EXP, title, domain, servidor, folders, date):

    #folders   = config_folders.config_folders(servidor)
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']

    if 'yakaira' in servidor:
        prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')
        fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'

    elif 'cnrm' in servidor:
        prov = np.genfromtxt("provincias.txt", delimiter='')
        fn = 'etopo1_bedrock.nc'
    
    
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
    
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,220)   
    
    #prefix   = 'wrfout_'+domain+'_2019-01-25_'+title
    prefix   = 'wrfout_'+domain+'_'+date+'_'+title
    #prefix   = 'wrfout_'+domain+'_2018-11-04_'+title

    filename = os.path.join(WRFfolder, prefix+':00')

    
    ncfile       = Dataset(filename,'r')        
    z            = wrf.getvar( ncfile,"z") 
    zh           = wrf.getvar(ncfile, "REFL_10CM")
    REFL_10CM    = wrf.interplevel(zh, z, 3000)
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")

    # ADEMAS SUMAR VIENTOS
    #----------    
    uvmet10 = wrf.getvar(ncfile, "uvmet10")   # Default is ‘m s-1’.
        
    fig, ax = plt.subplots(figsize=(8,8)) 
    pcm = ax.pcolormesh(lon, lat,  REFL_10CM, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
    cbar = plt.colorbar(pcm, ax=ax, shrink=1)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
        
    if 'Maite' in title: 
        ax.set_xlim([-65,-63.9]); 
        ax.set_ylim([-33,-31.5])
                
    else:
        ax.set_xlim([-67.5,-62]); 
        ax.set_ylim([-35,-31])
        
        #ax.set_xlim([-65.5,-62]); 
        #ax.set_ylim([-32.5,-31])

        #ax.set_xlim([-65.5,-62]); orig. 
        #ax.set_ylim([-35,-31])
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
    # agrego contorno de 500 y 1000m
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)
        
    # wind barbs    
    resobarb = 20
    ax.barbs(wrf.to_np(lon[::resobarb,::resobarb]), wrf.to_np(lat[::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[0,::resobarb,::resobarb]), 
             wrf.to_np(uvmet10[1,::resobarb,::resobarb]), length=6)
    

    cbar.cmap.set_under('white')
    ax.grid()
    ax.set_title('Zh interp. at 3km at '+title)
    
    #ax.text(x=-65, y=-34.9, s='120 and 220 km radar rings')
    #plt.show()
    fig.savefig(save_dir_compare+'/'+EXP+'/WRF_ZH3km_'+domain+'_ZOOM_general_evolution_'+title+'.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()
    
    return


    
    
#------------------------------------------------------------------------------
def plot_ZH1_COLMAX_WRF(EXP, title, folders):

    folders   = config_folders_final.config_folders('yakaira')
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']

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
    
    prefix   = 'wrfout_d02_2018-11-10_'+title
    filename = os.path.join(WRFfolder, 'wrfout_d02_2018-11-10_'+title+':00')

    
    ncfile       = Dataset(filename,'r')        
    z            = wrf.getvar( ncfile,"z") 
    zh           = wrf.getvar(ncfile, "REFL_10CM")
    REFL_10CM    = np.nanmax(zh, axis=0)
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")

    fig, ax = plt.subplots(figsize=(8,8)) 
    pcm = ax.pcolormesh(lon, lat,  REFL_10CM, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
    cbar = plt.colorbar(pcm, ax=ax, shrink=1)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
        
    if 'Maite' in title: 
        ax.set_xlim([-65,-63.9]); 
        ax.set_ylim([-33,-31.5])
                
    else:
        ax.set_xlim([-65.5,-62]); 
        ax.set_ylim([-35,-31])
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
    # agrego contorno de 500 y 1000m
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['k','k'], linewidths=2)
        
    cbar.cmap.set_under('white')
    ax.grid()
    ax.set_title('Zh COLMAX at '+title)
    
    ax.text(x=-65, y=-34.9, s='120 and 220 km radar rings')
    plt.show()
    fig.savefig(save_dir_compare+'/'+EXP+'/WRF_ZHcolmax_general_evolution_'+title+'.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()
     
    return#------------------------------------------------------------------------------
    
    
#------------------------------------------------------------------------------
def plot_general_WRF_evolution(WRFfolder, title, save_dir_compare, other_WRFfolder, OTHER_NAME):

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
    
    start_time_1    = 1200
    
    all_files_ = sorted(glob.glob(os.path.join(WRFfolder, 'wrfout_d02_2018-11-10*')))

    prefix = 'wrfout_d02_2018-11-10_'
    start_index = all_files_[0].find(prefix) + len(prefix)
            
    # filter and keep only the hh:00 or hh:30 files
    all_files = [f for f in all_files_ if f[-5:-3] in ("00", "30")]
        
    # Filter filenames based on HH:MM > '12:00' and < 22:00 HARDCODED
    filtered_files = [
        file for file in all_files
        if start_time_1-10 < int(os.path.basename(file).split('_')[3].split(':')[0]) * 100 +
           int(os.path.basename(file).split('_')[3].split(':')[1]) < 2200 ]
        
    fig, axes = plt.subplots(nrows=4, ncols=5, constrained_layout=True,figsize=[20,20])
    for ax, filename in zip(axes.flat, filtered_files):
        times = filename[start_index:start_index+5]   

        ncfile       = Dataset(filename,'r')        
        REFL_10CM    = np.nanmax(wrf.getvar(ncfile, "REFL_10CM"), axis=0)
        lat          = wrf.getvar( ncfile,"lat") 
        lon          = wrf.getvar( ncfile,"lon")

        pcm = ax.pcolormesh(lon, lat,  REFL_10CM, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
        cbar = plt.colorbar(pcm, ax=ax, shrink=1)
        ax.plot(prov[:,0],prov[:,1],color='k'); 


        #ax.set_xlim([-68,-62]); 
        #ax.set_ylim([-36,-30])
        
        if 'Maite' in title: 
            ax.set_xlim([-65,-63.9]); 
            ax.set_ylim([-33,-31.5])
                
        else:
            ax.set_xlim([-65.5,-62]); 
            ax.set_ylim([-35,-31])
            ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
            ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
        # agrego contorno de 500 y 1000m
        ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['k','k'], linewidths=2)
        
        cbar.cmap.set_under('white')
        ax.grid()
        ax.set_title(title+' '+times)
        plt.suptitle('COLMAX')
    
# =============================================================================
#         #read the other model output to add contours
#         THEFILE = os.path.basename(filename)
#         ncfile_other       = Dataset(other_WRFfolder+THEFILE,'r')
#         print(other_WRFfolder+THEFILE)
#         REFL_10CM = np.nanmax(wrf.getvar(ncfile_other, "REFL_10CM"), axis=0)
#         lat       =  wrf.getvar( ncfile_other,"lat") 
#         lon       =  wrf.getvar( ncfile_other,"lon")
#         cont = ax.contour(lon, lat, REFL_10CM, levels=[50], colors=['darkblue'], linewidths=1.2)
#         # ALSO ADD THE ACTUAL OBSERVATION CONTOUR? 
#        cont2 = ax.contour(lons_, lats_, colmax, levels=[50], colors=['magenta'], linewidths=2)
        
# =============================================================================
    ax.text(x=-65, y=-34.8, s='120 and 220 km radar rings')
        
    plt.show()
    if 'Maite' in title: 
        titleforfolder = title.replace('_Maitezoom', '')
        fig.savefig(save_dir_compare+'/'+titleforfolder+'/WRF_ZHcolmax_general_evolution_'+title+'.png', dpi=300,transparent=False)
        plt.close()
    else: 
        fig.savefig(save_dir_compare+'/'+title+'/WRF_ZHcolmax_general_evolution_'+title+'.png', dpi=300,transparent=False)
        plt.close()

    fig, axes = plt.subplots(nrows=4, ncols=5, constrained_layout=True,figsize=[20,20])
    for ax, filename in zip(axes.flat, filtered_files):
        times = filename[start_index:start_index+5]   

        ncfile       = Dataset(filename,'r')        
        #REFL_10CM    = np.nanmax(wrf.getvar(ncfile, "REFL_10CM"), axis=0)
        z            = wrf.getvar( ncfile,"z") 
        zh           = wrf.getvar(ncfile, "REFL_10CM")
        REFL_10CM    = wrf.interplevel(zh, z, 1000)
        lat          = wrf.getvar( ncfile,"lat") 
        lon          = wrf.getvar( ncfile,"lon")

        pcm = ax.pcolormesh(lon, lat,  REFL_10CM, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
        cbar = plt.colorbar(pcm, ax=ax, shrink=1)
        ax.plot(prov[:,0],prov[:,1],color='k'); 


        #ax.set_xlim([-68,-62]); 
        #ax.set_ylim([-36,-30])
        
        if 'Maite' in title: 
            ax.set_xlim([-65,-63.9]); 
            ax.set_ylim([-33,-31.5])
                
        else:
            ax.set_xlim([-65.5,-62]); 
            ax.set_ylim([-35,-31])
            ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
            ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
        # agrego contorno de 500 y 1000m
        ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['k','k'], linewidths=2)
        
        cbar.cmap.set_under('white')
        ax.grid()
        ax.set_title(title+' '+times)      
        plt.suptitle('Zh interp. at 1km')

# =============================================================================
#         #read the other model output to add contours
#         THEFILE = os.path.basename(filename)
#         ncfile_other       = Dataset(other_WRFfolder+THEFILE,'r')
#         print(other_WRFfolder+THEFILE)
#         REFL_10CM = np.nanmax(wrf.getvar(ncfile_other, "REFL_10CM"), axis=0)
#         lat       =  wrf.getvar( ncfile_other,"lat") 
#         lon       =  wrf.getvar( ncfile_other,"lon")
#         cont = ax.contour(lon, lat, REFL_10CM, levels=[50], colors=['darkblue'], linewidths=1.2)
#         # ALSO ADD THE ACTUAL OBSERVATION CONTOUR? 
#        cont2 = ax.contour(lons_, lats_, colmax, levels=[50], colors=['magenta'], linewidths=2)
        
# =============================================================================
    ax.text(x=-65, y=-34.8, s='120 and 220 km radar rings')
        
    plt.show()
    
    if 'Maite' in title: 
        titleforfolder = title.replace('_Maitezoom', '')
        fig.savefig(save_dir_compare+'/'+titleforfolder+'/WRF_Zh1km_general_evolution_'+title+'.png', dpi=300,transparent=False)
        plt.close()
    else:
        fig.savefig(save_dir_compare+'/'+title+'/WRF_Zh1km_general_evolution_'+title+'.png', dpi=300,transparent=False)
        plt.close()
    
    return#------------------------------------------------------------------------------
    
    
#------------------------------------------------------------------------------
def plot_general_radar_evolution(radar_folder, title, save_dir_compare, elev):

    prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')

    print(title)
    # RMA1 
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    TH_name       = 'TH'
    time1infile  = 82
    
    wildfiles    = 'cfrad.20181110*01.nc'  # try also with .02.nc maybe?     
    file_pattern = os.path.join(radar_folder, wildfiles)
    file_list    = sorted(glob.glob(file_pattern))
    timesr       = [filename[time1infile:time1infile+4] for filename in file_list]    
    
    # Step 1: Round each time to the nearest 10 minutes
    rounded_minutes = [int(round(funs.hhmm_to_minutes(t) / 10) * 10) for t in timesr]

    # Step 2: Create an array with NaN values for all intervals
    timereso   = 30
    start_time = 17*60 # remplazar por 12
    end_time   = 21*60 #(21*60)+30
        
    all_times_minutes = np.arange(start_time, end_time + timereso, timereso) 
    time_array        = np.full(len(all_times_minutes), np.nan, dtype=object)
    filename_array    = np.full(len(all_times_minutes), np.nan, dtype=object)
    
    # Step 3: Fill in times (and filenames) where data is available
    for i, rounded_time in enumerate(rounded_minutes):
        index = np.where(all_times_minutes == rounded_time)[0]
        if index.size > 0:
            time_array[index[0]] = funs.minutes_to_hhmm(rounded_time)
            filename_array[index[0]] = file_list[i]

    # Convert all_times_minutes back to HHMM format for easy reading
    all_times_hhmm = [funs.minutes_to_hhmm(m) for m in all_times_minutes]    
                
    fig, axes = plt.subplots(nrows=4, ncols=5, constrained_layout=True,figsize=[20,20])
    for ax, filename in zip(axes.flat, filename_array):
        times = filename[time1infile:time1infile+4]

        radar       = pyart.io.read(filename) 
        
        # Configure a gatefilter to filter out copolar correlation coefficient values > 0.9
        gatefilter = pyart.filters.GateFilter(radar)
        gatefilter.exclude_transition()
        gatefilter.exclude_equal('RHOHV', 0.9)
        
        ZHelev18 = radar.get_field(elev, TH_name, copy=False)
        [lats_, lons_, _] = radar.get_gate_lat_lon_alt(elev, reset_gate_coords=False, filter_transitions=False)
        pcm = ax.pcolormesh(lons_, lats_, ZHelev18, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
        cbar = plt.colorbar(pcm, ax=ax, shrink=1, label='Zh RMA1 elev 3')
        cbar.cmap.set_under('white')
        ax.grid()       

        ax.plot(prov[:,0],prov[:,1],color='k'); 
        
        #ax.set_xlim([-68,-62]); 
        #ax.set_ylim([-36,-30])

        ax.set_xlim([-65.5,-62]); 
        ax.set_ylim([-35,-31])
        

        ax.set_title(title+' '+times)          
        
    plt.show()
    fig.savefig(save_dir_compare+'OBS/'+title+'/radarOBS_general_evolution_'+title+'.png', dpi=300,transparent=False)
    
    return

def plot_domain(server, exp):
    
    
    folders=config_folders_final.config_folders(server)
    save_dir = folders['save_dir_compare']


    prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')

    folders=config_folders_final.config_folders('yakaira')
    
    ncfile_d01 = Dataset( folders[exp]+'wrfout_d01_2018-11-10_17:00:00','r')     
    ncfile_d02 = Dataset( folders[exp]+'wrfout_d02_2018-11-10_16:00:00','r')    
    
    latd01          = wrf.getvar( ncfile_d01,"lat") 
    lond01          = wrf.getvar( ncfile_d01,"lon")

    latd02          = wrf.getvar( ncfile_d02,"lat") 
    lond02          = wrf.getvar( ncfile_d02,"lon")
    
    LANDMASK    = wrf.getvar(ncfile_d01, "LANDMASK")
        
    lonmin = lond02.min()
    latmin = latd02.min()
    lonmax = lond02.max()
    latmax = latd02.max()
    vertices = [(lonmin,latmin),(lonmax,latmin),(lonmax,latmax),(lonmin,latmax)]
    polygon = Polygon(vertices, closed=True, edgecolor='blue', facecolor='lightblue', linewidth=2)
    
    fig, ax = plt.subplots(figsize=(8,8)) 
    ax.pcolormesh(lond01, latd01, LANDMASK);
    ax.add_patch(polygon)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
    ax.set_ylim([-43,-20])
    ax.set_title(exp)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    plt.show()
    fig.savefig(save_dir+'/domain_'+exp+'WRF.png', dpi=300,transparent=False)
    plt.close()
 
    return


#------------------------------------------------------------------------------
def single_WRF_files(WRFfolder, title, save_dir_compare): 

    
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
    
    start_time_1    = 1200
    
    all_files_ = sorted(glob.glob(os.path.join(WRFfolder, 'wrfout_d02_2018-11-10*')))

    prefix = 'wrfout_d02_2018-11-10_'
    start_index = all_files_[0].find(prefix) + len(prefix)
            
    # filter and keep only the hh:00 or hh:30 files
    all_files = [f for f in all_files_ if f[-5:-3] in ("00", "30")]
        
    # Filter filenames based on HH:MM > '12:00' and < 22:00 HARDCODED
    filtered_files = [
        file for file in all_files
        if start_time_1-10 < int(os.path.basename(file).split('_')[3].split(':')[0]) * 100 +
           int(os.path.basename(file).split('_')[3].split(':')[1]) < 2200 ]
        
    # for filename in filtered_files:
        
    #     times = filename[start_index:start_index+5]   

    #     ncfile       = Dataset(filename,'r')        
    #     #REFL_10CM    = np.nanmax(wrf.getvar(ncfile, "REFL_10CM"), axis=0)
    #     z            = wrf.getvar( ncfile,"z") 
    #     zh           = wrf.getvar(ncfile, "REFL_10CM")
    #     REFL_10CM    = wrf.interplevel(zh, z, 1000)
    #     lat          = wrf.getvar( ncfile,"lat") 
    #     lon          = wrf.getvar( ncfile,"lon")

    #     fig, ax = plt.subplots(figsize=(8,8)) 
    #     pcm = ax.pcolormesh(lon, lat,  REFL_10CM, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
    #     cbar = plt.colorbar(pcm, ax=ax, shrink=1)
    #     ax.plot(prov[:,0],prov[:,1],color='k'); 
    
    #     if 'Maite' in title: 
    #         ax.set_xlim([-65,-63.9]); 
    #         ax.set_ylim([-33,-31.5])
    #     else:
    #         ax.set_xlim([-65.5,-62]); 
    #         ax.set_ylim([-35,-31])
    #         ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    #         ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
            
    #     ax.set_title(title+' '+str(times))
    #     fig.savefig(save_dir_compare+'/WRF_Zh1km_'+title+'_time'+str(times)+'.png', dpi=300,transparent=False)
    #     plt.close()
        
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
        ax.plot(prov[:,0],prov[:,1],color='k'); 
    
        if 'Maite' in title: 
            ax.set_xlim([-65,-63.9]); 
            ax.set_ylim([-33,-31.5])
        else:
            ax.set_xlim([-65.5,-62]); 
            ax.set_ylim([-35,-31])
            ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
            ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
            
        ax.set_title(title+' '+str(times))
        fig.savefig(save_dir_compare+'/WRF_RAINC_'+title+'_time'+str(times)+'.png', dpi=300,transparent=False, bbox_inches='tight')
        plt.close()
        counter=counter+1
        
        
    return


#------------------------------------------------------------------------------
def plot_common_transect_level_MAP(EXP_WSM6, time):

    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944

    folders=config_folders_final.config_folders('yakaira')
    save_dir_compare=folders['save_dir_compare']    
    #---------------------------------------------------------
    # WSM6 files
    WSM6_file = os.path.join(folders[EXP_WSM6], 'wrfout_d02_2018-11-10_'+time+':00')
    
    #---------------------------------------------------------
    # P3 files
    #folderP3_1    = os.path.join('/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_P3_3MOM_LF_v4.5.2/', 'wrfout_d02*00:00')
    #file_listP3_1     = sorted(glob.glob(folderP3_1))
    #timeNr_P3_1   = 98
    #times_p3          = [filename[timeNr_P3_1:timeNr_P3_1+16] for filename in file_listP3_1]
    
    #---------------------------------------------------------
    Re       = 6.3781e6
    ilevels_plot = np.array([0,1,2,3,4,5,6,7,8,9,10,15,20,25,30])

    #---------------------------------------------------------
    wrf_file = Dataset(WSM6_file,'r')
    temp     = wrf.g_temp.get_tk(wrf_file)
    pressure = wrf.g_pressure.get_pressure(wrf_file)
    geopo_p  = wrf.g_geoht.get_height(wrf_file) # geopotential height as Mean Sea Level (MSL
    z_level  = Re*geopo_p/(Re-geopo_p)
    lat      =  wrf.getvar( wrf_file,"lat") 
    lon      =  wrf.getvar( wrf_file,"lon")
    
    qr = funs.mixr2massconc( np.squeeze(wrf_file.variables["QRAIN"][0,:,:,:]  ), pressure, temp )        
    qi = funs.mixr2massconc( np.squeeze(wrf_file.variables["QICE"][0,:,:,:]   ), pressure, temp )        
    qc = funs.mixr2massconc( np.squeeze(wrf_file.variables["QCLOUD"][0,:,:,:] ), pressure, temp )       
    qs = funs.mixr2massconc( np.squeeze(wrf_file.variables["QSNOW"][0,:,:,:] ), pressure, temp )       
    qg = funs.mixr2massconc( np.squeeze(wrf_file.variables["QGRAUP"][0,:,:,:] ), pressure, temp )       
                
    # Similarly get the P3_1 info
    # wrf_file_P3 = Dataset(file_listP3_1[item],'r')
    # print(file_listP3_1[item])
          
    # temp_p3     = wrf.g_temp.get_tk(wrf_file_P3)
    # pressure_p3 = wrf.g_pressure.get_pressure(wrf_file_P3)
    # geopo_p_p3  = wrf.g_geoht.get_height(wrf_file_P3) # geopotential height as Mean Sea Level (MSL
    # z_level_p3  = Re*geopo_p_p3/(Re-geopo_p_p3)
    # lat_p3      =  wrf.getvar( wrf_file_P3,"lat") 
    # lon_p3      =  wrf.getvar( wrf_file_P3,"lon")
    
    # qr_p3 = mixr2massconc( np.squeeze(wrf_file_P3.variables["QRAIN"][0,:,:,:]  ), pressure_p3, temp_p3 )        
    # qi_p3 = mixr2massconc( np.squeeze(wrf_file_P3.variables["QICE"][0,:,:,:]   ), pressure_p3, temp_p3 )        
    # qc_p3 = mixr2massconc( np.squeeze(wrf_file_P3.variables["QCLOUD"][0,:,:,:] ), pressure_p3, temp_p3 )       

            
    for iz in ilevels_plot:
            
        fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True,figsize=[12,12])
        pcm0 = axes[0,0].pcolormesh(lon, lat, qr[iz,:,:], vmin=1e-6, vmax=0.004)
        pcm1 = axes[1,0].pcolormesh(lon, lat, qi[iz,:,:] + qs[iz,:,:] + qg[iz,:,:], vmin=1e-6, vmax=0.004)       
        
        cbar = plt.colorbar(pcm0, ax=axes[0,0], shrink=1, label='qr [kg/m$^2$]')
        cbar.cmap.set_under('white')
        cbar = plt.colorbar(pcm1, ax=axes[1,0], shrink=1, label='qi+qs+qg [kg/m$^2$]')
        cbar.cmap.set_under('white')
        
        for i in range(2):
            axes[i,0].grid()
            axes[i,0].set_xlim([-65.5,-62])
            axes[i,0].set_ylim([-33.5,-31.3])
            
            axes[0,0].set_title('WSM6 q_rain')                                                  
            axes[1,0].set_title('WSM6 q_i(tot)')  

        #pcm2 = axes[0,1].pcolormesh(lon_p3, lat_p3, qr_p3[iz,:,:], vmin=1e-5, vmax=0.001)
        #pcm3 = axes[1,1].pcolormesh(lon_p3, lat_p3, qi_p3[iz,:,:],vmin=1e-5, vmax=0.001)        
        #
        #cbar = plt.colorbar(pcm2, ax=axes[0,1], shrink=1, label='qr [kg/m$^2$]')
        #cbar.cmap.set_under('white')
        #cbar = plt.colorbar(pcm3, ax=axes[1,1], shrink=1, label='qi [kg/m$^2$]')
        #cbar.cmap.set_under('white')
        
        # for i in range(2):
        #     axes[i,1].grid()
        #     axes[i,1].set_xlim([-65.5,-62])
        #     axes[i,1].set_ylim([-33.5,-31.3])
            
        #     axes[0,1].set_title('P3 3MOM LF q_rain')                                                  
        #     axes[1,1].set_title('P3 3MOM LF q_i')  

        #plt.suptitle( 'Level qx (z='+ iz +') ' + title + 'UTC')
        plt.suptitle('Layer qx '+ ' (level:'+ str(iz)+' '+time +'UTC)')

        # RMA1 
        [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)
        for iax in range(2): 
            axes[iax,0].plot(lon_radius, lat_radius, 'k', linewidth=0.8, label='RMA1 (120km)')
            axes[iax,1].plot(lon_radius, lat_radius, 'k', linewidth=0.8, label='RMA1 (120km)')

        axes[0,0].legend(fontsize=10, loc='lower left')
         
        # #
        #xx = np.vstack([x_wsm6[[0,2]],x_wsm6[[1,3]]])
        #yy = np.vstack([y_wsm6[[0,2]],y_wsm6[[1,3]]])
        #axes[0].plot(xx, yy, '-or' , linewidth=1.2)
        # #
        # xx = np.vstack([x_P3_3mom_LF[[0,2]],x_P3_3mom_LF[[1,3]]])
        # yy = np.vstack([y_P3_3mom_LF[[0,2]],y_P3_3mom_LF[[1,3]]])
        # axes[1].plot(xx, yy, '-or' , linewidth=1.2)
        #
        plt.show()
        fig.savefig(save_dir_compare+'/Layer_qx_leveliz'+str(iz)+'_time'+str(time)+'.png', dpi=300,transparent=False, bbox_inches='tight')


    return



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_WRF_diffvar_only(EXP1, EXP2, intplev, title):
    
    import matplotlib

    folders   = config_folders_final.config_folders('yakaira')
    WRFfolder1 = folders[EXP1]
    WRFfolder2 = folders[EXP2]
    save_dir_compare = folders['save_dir_compare']

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
    
    # Create the start point and end point for the cross section
    lon_points = [-63.9, -63.9]
    lat_points = [-34.5, -31.5]


    prefix   = 'wrfout_d02_2018-11-10_'+title

    #------ READ WRF variables of interest ------------------------------------
    filename1 = os.path.join(WRFfolder1, 'wrfout_d02_2018-11-10_'+title+':00')
    ncfile       = Dataset(filename1,'r')        
    
    z            = wrf.getvar( ncfile,"z") 
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")
    p            = wrf.getvar(ncfile, "pressure")
    theta        = wrf.getvar(ncfile, "theta", units="K")   # potential temperature 
    t2           = wrf.getvar(ncfile, "T2", timeidx=-1)     # 2-m temperature
    eth          = wrf.g_temp.get_eth(ncfile)               # equivalent potential temperature theta_e tita-e
    wspd         = wrf.getvar(ncfile, "wspd_wdir", units="kts")[0,:]
    eth_850      = wrf.interplevel(eth, p, intplev)
    theta_850    = wrf.interplevel(theta, p,intplev )
    ua           = wrf.getvar(ncfile, "ua", units="kt")
    va           = wrf.getvar(ncfile, "va", units="kt")
    u850         = wrf.interplevel(ua, p, intplev)
    v850         = wrf.interplevel(va, p, intplev)

    
    #------ READ WRF variables of interest ------------------------------------
    del ncfile
    filename2 = os.path.join(WRFfolder2, 'wrfout_d02_2018-11-10_'+title+':00')
    ncfile       = Dataset(filename2,'r')        
    
    z1            = wrf.getvar( ncfile,"z") 
    lat1          = wrf.getvar( ncfile,"lat") 
    lon1          = wrf.getvar( ncfile,"lon")
    p1            = wrf.getvar(ncfile, "pressure")
    theta1        = wrf.getvar(ncfile, "theta", units="K")   # potential temperature 
    t21           = wrf.getvar(ncfile, "T2", timeidx=-1)     # 2-m temperature
    eth1          = wrf.g_temp.get_eth(ncfile)               # equivalent potential temperature theta_e tita-e
    wspd1         = wrf.getvar(ncfile, "wspd_wdir", units="kts")[0,:]
    eth_8501      = wrf.interplevel(eth1, p1, intplev)
    theta_8501    = wrf.interplevel(theta1, p1, intplev)
    ua1           = wrf.getvar(ncfile, "ua", units="kt")
    va1           = wrf.getvar(ncfile, "va", units="kt")
    u8501         = wrf.interplevel(ua1, p1, intplev)
    v8501         = wrf.interplevel(va1, p1, intplev)
    
    #--------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8,8)) 
    pcm            = ax.contourf(lon, lat,  eth_850-eth_8501,  cmap=matplotlib.cm.get_cmap("viridis_r"), levels=np.arange(-20,20,5), vmin=-20, vmax=20)
    cbar = plt.colorbar(pcm, ax=ax, shrink=1, label='eth ('+str(intplev)+' hPa)')
    ax.plot(prov[:,0],prov[:,1],color='k'); 
    cbar.cmap.set_under('white')
    ax.set_xlim([-65.5,-62]); 
    ax.set_ylim([-34.5,-31])

    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
    
    # agrego vertical cross section selection:
    ax.plot(lon_points, lat_points, '-r', linewidth=1.2)
    # agrego contorno de 500 y 1000m
    
    # wind barbs at 850hPa 
    resobarb = 20
    ax.barbs(wrf.to_np(lon[::resobarb,::resobarb]), wrf.to_np(lat[::resobarb,::resobarb]), 
             wrf.to_np(u850[::resobarb,::resobarb])- wrf.to_np(u8501[::resobarb,::resobarb]), 
             wrf.to_np(v850[::resobarb,::resobarb])- wrf.to_np(v8501[::resobarb,::resobarb]), length=6)

    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)    
    ax.set_title(EXP1+'-'+EXP2+ ' ('+title+')')
    fig.savefig(save_dir_compare+'/'+'vertical_crossSection/diff_eth'+str(intplev)+'hPa_map'+title+'.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()


    return#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_WRF_hovmoller_thetae(EXP):
    
    import matplotlib
    import xarray as xr

    folders   = config_folders_final.config_folders('yakaira')
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']
    
    # Define bounds
    lat_min, lat_max = -34.5, -31
    lon_min, lon_max = -64.5, -64.3
    
    prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter=''); 

    theta_e_850_list = []
    time_list    = []
    for h in range(16, 23):
        for m in range(0, 60, 30):    
            prefix   = 'wrfout_d02_2018-11-10_'+f"{h}:{m:02d}"+':00'
            filename = os.path.join(WRFfolder, prefix) 
            ncfile   = Dataset(filename,'r')       
            p         = wrf.getvar(ncfile, "pressure")
            eth       = wrf.g_temp.get_eth(ncfile)   
            eth_850   = wrf.interplevel(eth, p, 850) 
        
            # recortar slice of interest for hovmoller:
            lat          = wrf.getvar( ncfile,"lat") 
            lon          = wrf.getvar( ncfile,"lon")
                
            #--------------------------
            # Convert time to a numeric format (if needed, extract from filename)
            time_str = f"{h}:{m:02d}"
            time_list.append(time_str)  
    
            #--------------------------
            # data: 
            theta_e_850_list.append(wrf.to_np(eth_850))              
        
        
    # Convert lists to numpy arrays
    theta_e_850_arr = np.array(theta_e_850_list)  # Shape: (time, lat, lon)
    data_xr = xr.DataArray(theta_e_850_arr, dims=["time", "lat", "lon"],                           
                           coords={"time": time_list, "lat": lat.data[:,0], "lon": lon.data[0,:]}) 
    hovmoller = data_xr.sel(lon =-64.3, lat = np.arange(-34.5,-31, 0.01), method='nearest' )



    #--------------------------------------------------------------------------
    lat_     = np.arange(-34.5,-31, 0.01) 
    times    = np.arange(0,14,1)
    timetime = np.zeros( (len(times), lat_.shape[0]) )
    latlat = np.zeros( (len(times), lat_.shape[0]) )
    for i in range(len(times)):
        latlat[i,:]=lat_
        timetime[i,:]=times[i]
    
    #--------------------------------------------------------------------------
    # Plot Hovmöller diagram (Time vs Longitude)
    fig=plt.figure(figsize=(10,8))
    plt.contourf(timetime, latlat, hovmoller, cmap="jet", levels=20)
    plt.colorbar(label="Equivalent Potential Temperature 850 hPa (K)")
    plt.ylabel("Latitude")
    plt.xlabel("Time (UTC)")
    plt.title(EXP+ " Hovmöller Diagram of 850 hPa Equivalent Potential Temperature")
    plt.xticks(ticks=np.arange(len(time_list)), labels=time_list, rotation=45)
    plt.show()
    fig.savefig(save_dir_compare+'/'+EXP+'/LonHovmoller_.png', dpi=300,transparent=False,bbox_inches='tight')


    return#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_WRF_var_only(EXP, intplev, title, folders, domain):    
    import matplotlib

    #folders   = config_folders.config_folders('yakaira')
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']

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
    
    prefix   = 'wrfout_'+domain+'_2018-11-10_'+title
    filename = os.path.join(WRFfolder, 'wrfout_'+domain+'_2018-11-10_'+title+':00')
    ncfile       = Dataset(filename,'r')        
    
    #------ READ WRF variables of interest ------------------------------------
    z            = wrf.getvar( ncfile,"z") 
    zh           = wrf.getvar(ncfile, "REFL_10CM")
    COLMAXREFL   = np.nanmax(zh, axis=0)
    zh_1km       = wrf.interplevel(zh, z, 1000)
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")
    uvmet10      = wrf.getvar(ncfile, "uvmet10",units="kt")   # Default is ‘m s-1’.
    wspd_wdir10  = wrf.getvar(ncfile, "wspd_wdir10") 
    #pertT        = wrf.getvar(ncfile, "T")                    # perturbation potential temperature  % (theta-t0) (K) 
    
    # I want vertical cross section equivalent potential temperature and potential
    # temperature perturbations, cross sectional winds. 
    
    # Create the start point and end point for the cross section
    lon_points = [-63.9, -63.9]
    lat_points = [-34.5, -31.5]

    #--------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8,8)) 
    pcm = ax.pcolormesh(lon, lat,  COLMAXREFL, cmap=P4A.colormaps('ref'), vmin=0,  vmax=60)
    cbar = plt.colorbar(pcm, ax=ax, shrink=1)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
    cbar.cmap.set_under('white')
    ax.set_xlim([-67.5,-62]); 
    #ax.set_xlim([-65.5,-62]); 
    ax.set_ylim([-34.5,-31])
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
    # agrego vertical cross section selection:
    ax.plot(lon_points, lat_points, '-r', linewidth=1.2)
    # agrego contorno de 500 y 1000m
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['k','k'], linewidths=2)    
    ax.set_title('Zh interp. at 1km at '+title)
    fig.savefig(save_dir_compare+'/'+EXP+'/vertical_crossSection/TransectaMap_'+title+'.png', dpi=300,transparent=False,bbox_inches='tight')
    #plt.close()
    #--------------------------------------------------------------------------
            
    p         = wrf.getvar(ncfile, "pressure")
    theta     = wrf.getvar(ncfile, "theta", units="K")   # potential temperature 
    t2        = wrf.getvar(ncfile, "T2", timeidx=-1)     # 2-m temperature
    eth       = wrf.g_temp.get_eth(ncfile)               # equivalent potential temperature theta_e tita-e
    wspd      = wrf.getvar(ncfile, "wspd_wdir", units="kts")[0,:]
    eth_850   = wrf.interplevel(eth, p, intplev)
    theta_850 = wrf.interplevel(theta, p,intplev)
    ua = wrf.getvar(ncfile, "ua", units="kt")
    va = wrf.getvar(ncfile, "va", units="kt")
    u850 = wrf.interplevel(ua, p, intplev)
    v850 = wrf.interplevel(va, p,intplev)
  
    #--------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8,8)) 
    #pcm            = ax.contourf(lon, lat,  eth_850,  cmap=matplotlib.cm.get_cmap("viridis_r"), levels=np.arange(290,390,10), vmin=280, vmax=380)
    pcm            = ax.contourf(lon, lat,  eth_850,  cmap=matplotlib.cm.get_cmap("viridis"), levels=np.arange(320,365,5), vmin=315, vmax=360)


    ax.contour(lon, lat, COLMAXREFL, levels=[30], colors=['r'], linewidths=1.5)

    cbar = plt.colorbar(pcm, ax=ax, shrink=1)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
    cbar.cmap.set_under('white')

    ax.set_xlim([-67.5,-62]); 
    #ax.set_xlim([-65.5,-62]); 
    ax.set_ylim([-34.5,-31])

    # add theta contours and labels
    #theta_contours = ax.contour(lon, lat, theta_850, levels=np.arange(250,400,10), colors=['k','k'], linewidths=2)    
    #labels = plt.clabel(theta_contours, inline=1, fontsize=12, fmt="%i")
    #for text in labels:
    #    text.set_clip_on(True)
    # Add a white background to each label
    #for text in labels:
    #    text.set_bbox(dict(facecolor='white', edgecolor='none', pad=2))    
    
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
    
    # agrego vertical cross section selection:
    ax.plot(lon_points, lat_points, '-r', linewidth=1.2)
    # agrego contorno de 500 y 1000m
    
    # wind barbs at 850hPa 
    resobarb = 20
    ax.barbs(wrf.to_np(lon[::resobarb,::resobarb]), wrf.to_np(lat[::resobarb,::resobarb]), 
             wrf.to_np(u850[::resobarb,::resobarb]), 
             wrf.to_np(v850[::resobarb,::resobarb]), length=6)

    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)    
    ax.set_title(EXP+ ' eth ('+str(intplev)+'hPa) '+title)
    fig.savefig(save_dir_compare+'/'+EXP+'/vertical_crossSection/eth'+str(intplev)+'hPa_map'+title+'.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()


    
    
    
    # Compute the vertical cross-section interpolation.  Also, include the
    # # lat/lon points along the cross-section.
    start_point = wrf.CoordPair(lat=lat_points[0], lon=lon_points[0])
    end_point   = wrf.CoordPair(lat=lat_points[1], lon=lon_points[1])    
    
    theta_cross = wrf.vertcross(theta, z, wrfin=ncfile, start_point=start_point,
                        end_point=end_point, latlon=True, meta=True)
    eth_cross   = wrf.vertcross(eth, z, wrfin=ncfile, start_point=start_point,
                        end_point=end_point, latlon=True, meta=True)
    coord_pairs = wrf.to_np(eth_cross.coords["xy_loc"])

    #--------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8,8)) 
    contours = ax.contourf(np.arange(coord_pairs.shape[0]), wrf.to_np(eth_cross["vertical"]),
                            wrf.to_np(eth_cross), cmap=matplotlib.cm.get_cmap("jet"), levels=np.arange(300,400,10))
    plt.colorbar(contours, ax=ax)
        
    theta_contours = ax.contour(np.arange(coord_pairs.shape[0]), wrf.to_np(theta_cross["vertical"]),
                            wrf.to_np(theta_cross), np.arange(270.,400.,10.), colors='black', linestyles='dashed')
    plt.clabel(theta_contours, inline=1, fontsize=12, fmt="%i")

    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in wrf.to_np(coord_pairs)]
    ax.set_xticks(x_ticks[::10])
    ax.set_xticklabels(x_labels[::10], rotation=45, fontsize=10)
    
    # Set the x-axis and  y-axis labels
    ax.set_xlabel("Latitude, Longitude", fontsize=12)
    ax.set_ylabel("Height (m)", fontsize=12)
    ax.set_ylim([0, 10000])
    
    ax.set_title("Equivalent Potential Temp. " + EXP + " ("+title+")" )
    ax.grid()
    fig.savefig(save_dir_compare+'/'+EXP+'/vertical_crossSection/eth_cross_section_'+title+'.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()


    return#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_WRF_TyTd_only(EXP, title):
    
    import matplotlib

    folders   = config_folders_final.config_folders('yakaira')
    WRFfolder = folders[EXP]
    save_dir_compare = folders['save_dir_compare']

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
    
    prefix   = 'wrfout_d02_2018-11-10_'+title
    filename = os.path.join(WRFfolder, 'wrfout_d02_2018-11-10_'+title+':00')
    ncfile       = Dataset(filename,'r')        
    
    #------ READ WRF variables of interest ------------------------------------
    z            = wrf.getvar( ncfile,"z") 
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")
    uvmet10      = wrf.getvar(ncfile, "uvmet10",units="kt")   # Default is ‘m s-1’.
    wspd_wdir10  = wrf.getvar(ncfile, "wspd_wdir10") 
    p         = wrf.getvar(ncfile, "pressure")
    theta     = wrf.getvar(ncfile, "theta", units="K")   # potential temperature 
    t2        = wrf.getvar(ncfile, "T2", timeidx=-1)     # 2-m temperature
    td_2m     = wrf.getvar( ncfile,"td2")


    
    # I want vertical cross section equivalent potential temperature and potential
    # temperature perturbations, cross sectional winds. 
    
    # Create the start point and end point for the cross section
    lon_points = [-63.9, -63.9]
    lat_points = [-34.5, -31.5]

    #--------------------------------------------------------------------------
    #-------------------------------------------------------------------------
    
    colors1 = plt.cm.BrBG(np.linspace(0., 0.4, 120))
    colors2 = plt.cm.YlOrRd(np.linspace(0.0, 0.55, 70))
    colors3 = plt.cm.Greens(np.linspace(0.1, 0.7, 80))
    colors4 = plt.cm.Blues(np.linspace(0.2, 0.9, 80))
    
    colors = np.vstack((colors1,colors2,colors3,colors4))#, colors2, colors1))
    cmap_td= mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    cmap_td.set_over((0,0,0.69))
    cmap_td.set_under((0.2,0.1,0))

    
    fig, ax = plt.subplots(figsize=(8,8)) 
    pcm            = ax.contourf(lon, lat,  td_2m,  cmap=cmap_td, levels=np.arange(-12,30,2), vmin=-12, vmax=28)
    cbar = plt.colorbar(pcm, ax=ax, shrink=1)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
    cbar.cmap.set_under('white')
    ax.set_xlim([-65.5,-62]); 
    ax.set_ylim([-34.5,-31])
    
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
    
    # agrego vertical cross section selection:
    ax.plot(lon_points, lat_points, '-r', linewidth=1.2)

    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)    
    ax.set_title(EXP+ ' Td 2m '+title)
    fig.savefig(save_dir_compare+'/'+EXP+'/td2m_map'+title+'.png', dpi=300,transparent=False,bbox_inches='tight')


    #-------------------------------------------------------------------------
    import colormaps_new as cmapss
    colors1=np.array(cmapss._HomeyerRainbow_data[0:80])
    colors2=np.array(cmapss._HomeyerRainbow_data[100:230])
    colors3=np.array(cmapss._viridis_data[0:50])
    colors = np.vstack((colors3,colors1,colors2))
    cmap_temp = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    
    fig, ax = plt.subplots(figsize=(8,8)) 
    pcm            = ax.contourf(lon, lat,  t2-273,  cmap=cmap_temp, levels=np.arange(0,42,2), vmin=0, vmax=40)
    cbar = plt.colorbar(pcm, ax=ax, shrink=1)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
    cbar.cmap.set_under('white')
    ax.set_xlim([-65.5,-62]); 
    ax.set_ylim([-34.5,-31])
    
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
    
    # agrego vertical cross section selection:
    ax.plot(lon_points, lat_points, '-r', linewidth=1.2)

    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)    
    ax.set_title(EXP+ ' T 2m '+title)
    fig.savefig(save_dir_compare+'/'+EXP+'/t2m_map'+title+'.png', dpi=300,transparent=False,bbox_inches='tight')
        
    # # Compute the vertical cross-section interpolation.  Also, include the
    # # lat/lon points along the cross-section.
    # start_point = wrf.CoordPair(lat=lat_points[0], lon=lon_points[0])
    # end_point   = wrf.CoordPair(lat=lat_points[1], lon=lon_points[1])    
    
    # theta_cross = wrf.vertcross(theta, z, wrfin=ncfile, start_point=start_point,
    #                    end_point=end_point, latlon=True, meta=True)
    # eth_cross   = wrf.vertcross(eth, z, wrfin=ncfile, start_point=start_point,
    #                    end_point=end_point, latlon=True, meta=True)
    # coord_pairs = wrf.to_np(eth_cross.coords["xy_loc"])

    # #--------------------------------------------------------------------------
    # fig, ax = plt.subplots(figsize=(8,8)) 
    # contours = ax.contourf(np.arange(coord_pairs.shape[0]), wrf.to_np(eth_cross["vertical"]),
    #                        wrf.to_np(eth_cross), cmap=matplotlib.cm.get_cmap("jet"), levels=np.arange(300,400,10))
    # plt.colorbar(contours, ax=ax)
        
    # theta_contours = ax.contour(np.arange(coord_pairs.shape[0]), wrf.to_np(theta_cross["vertical"]),
    #                        wrf.to_np(theta_cross), np.arange(270.,400.,10.), colors='black', linestyles='dashed')
    # plt.clabel(theta_contours, inline=1, fontsize=12, fmt="%i")

    # x_ticks = np.arange(coord_pairs.shape[0])
    # x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in wrf.to_np(coord_pairs)]
    # ax.set_xticks(x_ticks[::10])
    # ax.set_xticklabels(x_labels[::10], rotation=45, fontsize=10)
    
    # # Set the x-axis and  y-axis labels
    # ax.set_xlabel("Latitude, Longitude", fontsize=12)
    # ax.set_ylabel("Height (m)", fontsize=12)
    # ax.set_ylim([0, 10000])
    
    # ax.set_title("Equivalent Potential Temp. " + EXP + " ("+title+")" )
    # ax.grid()
    # fig.savefig(save_dir_compare+'/'+EXP+'/vertical_crossSection/eth_cross_section_'+title+'.png', dpi=300,transparent=False,bbox_inches='tight')
    # #plt.close()


    return

#------------------------------------------------------------------------------
def run_transects(): 

    # Agregar transectas para el EXP='WSM6_domain3_NoahMP' y su correspondiente p3
    #-------------------------------------------------------------
    x_wsm6 = np.array([-64.3, -63.5])
    y_wsm6 = np.array([-31.75, -32.65])
    all_WSM6 = sorted(glob.glob(os.path.join(folders['WSM6_domain3_NoahMP'], 'wrfout_d02_2018-11-10*')))    
    prefix = 'wrfout_d02_2018-11-10_'
    start_index = all_WSM6[0].find(prefix) + len(prefix)
    # =============================================================================
    # TODO ESTO TENGO QUE ADAPTAR Y VOVLER A CORRER CON LAS SIGUIENTES VARIABLES UPDATED:
    x_P3_3mom_LF = np.array([-63.6, -63.6, -63.8, -62.0])
    y_P3_3mom_LF = np.array([-32.0, -33.0, -32.00, -33.5])  
    EXP_P3 = 'P3_'    
    # =============================================================================

    # One plot per layer of qx for WSM6 and P3
    plot_common_transect_level_MAP('WSM6_domain3_NoahMP', '20:00')
    
    # MAX_DBZ with transects to analyse 
    P4A.plot_common_transect_MAP('WSM6_domain3_NoahMP', EXP_P3, '20:00', folders['dow6_str'][4], 
                folders['dow7_str'][4], folders['csapr2_str'][4], 'MAXRADAR_WRF', x_wsm6, y_wsm6, x_P3_3mom_LF, y_P3_3mom_LF)
    
    # DBZ interp 1km
    P4A.plot_common_transect_MAP('WSM6_domain3_NoahMP', EXP_P3, '20:00', folders['dow6_str'][4], 
                folders['dow7_str'][4], folders['csapr2_str'][4], 'interp1km', x_wsm6, y_wsm6, x_P3_3mom_LF, y_P3_3mom_LF)

    # Now Plot the transects
    P4A.plot_common_transect_WSM6micro_ICE_derived('WSM6_domain3_NoahMP', x_wsm6, y_wsm6, '20:00')
    #P4A.plot_common_transect_P3

    # Domain average qxs
    for ii, istr in enumerate(all_WSM6): 
        wrf_file_wsm6 = Dataset(all_WSM6[ii],'r')
        #wrf_file_p3 = Dataset(strfile_P3_3MOM_LF[ii],'r')
        itime = all_WSM6[ii][start_index:start_index+5]
        # Plot maps of precip accumulation
        #plot_accums(wrf_file_wsm6, wrf_file_p3, itime, save_dir_compare+'/WRFcompare/', dow6_str[1], dow7_str[1], csapr2_str[1])
        # Plot domain average vertical profiles of qx. 
        P4A.plot_domain_qxs(wrf_file_wsm6, wrf_file_wsm6, itime, folders['save_dir_compare']+'/WSM6_domain3_NoahMP/')    
    
    for h in range(18, 22):
        for m in range(0, 60, 30):
            plot_WRF_intqx('P3_3MOM_LF_domain3_NoahMP', f"{h}:{m:02d}", 54)
            plot_WRF_intqx('WSM6_domain3_NoahMP', f"{h}:{m:02d}", 6)
            
        
    return

#------------------------------------------------------------------------------
# RUN MAIN
def main(exp, folders): 
    
    # Plot a general domain subplot 4x5 map of WRF_DBZ every 30 min
    plot_general_WRF_evolution(WRFfolder=folders[exp], title=exp, save_dir_compare=folders['save_dir_compare'], 
                               other_WRFfolder='', OTHER_NAME='')

    plot_general_WRF_evolution(WRFfolder=folders[exp], title=exp+'_Maitezoom', save_dir_compare=folders['save_dir_compare'], 
                               other_WRFfolder='', OTHER_NAME='')
    
    
    # Make a foler for each experiment and save figure for each time
    newdir = os.path.join(folders['save_dir_compare'], exp) 
    if not os.path.exists(newdir):
            os.makedirs(newdir)
    single_WRF_files(WRFfolder=folders[exp], title=exp+'_Maitezoom', save_dir_compare=newdir)
    single_WRF_files(WRFfolder=folders[exp], title=exp, save_dir_compare=newdir)
    
    return

#------------------------------------------------------------------------------
def run_all(EXPs): 
    
    # borrar los config foldeders de adentro y pasarlos afurea! 
    folders=config_folders_final.config_folders('yakaira')

    plot_domain('yakaira', 'WSM6_domain2')
    plot_domain('yakaira', 'WSM6_domain3')
    plot_domain('yakaira', 'WSM6_domain3_NoahMP')
    plot_domain('yakaira', 'WSM6_domain4_NoahMP')
    plot_domain('yakaira', 'P3_3MOM_LF_domain3_NoahMP_highres')
    #plot_domain('yakaira', 'P3_3MOM_LF_domain3_NoahMP_lowres')
    #plot_domain('yakaira', 'P3_3MOM_LF_domain5_NoahMP')

    #main('WSM6_domain2', folders)
    #main('WSM6_domain3', folders)
    #main('WSM6_domain3_NoahMP', folders)
    #main('WSM6_domain4_NoahMP', folders)
    #main('P3_3MOM_LF_domain3_NoahMP', folders)
    #main('THOM_domain3_NoahMP', folders)
        
    # EXPS que analizo finalmente: 
    #EXPs = ['WSM6_domain3', 'WSM6_domain3_NoahMP', 'P3_3MOM_LF_domain3_NoahMP', 
    #        'THOM_domain3_NoahMP', 'WDM6_domain3_NoahMP', 'WSM6_domain3_YSU_noNoahMP']
    # , 'P3_3MOM_LF_domain_noNoah']
    #
    #EXPs = ['P3mp54_domain3_YSU_noNoahMP']    
    for EXP in EXPs:
        for h in range(12, 23):
            for m in range(0, 60, 30):
                plot_ZH1km_WRF_wWRFwind(EXP, f"{h}:{m:02d}", 'd02', 850, folders)
                plot_ZH1km_WRF_wWRFwind(EXP, f"{h}:{m:02d}", 'd02', 950, folders)
                plot_ZH1km_WRF_wWRF_shelicity(EXP, f"{h}:{m:02d}", folders, 'd02')
                plot_ZH1km_WRF(EXP, f"{h}:{m:02d}", 'd02', 'yakaira', folders)   #plot_ZH1km_WRF(EXP, title, domain, servidor):
                plot_ZH1km_WRF_wWRF_uhelicity_only(EXP, 850, f"{h}:{m:02d}", folders, 'd02')
                plot_ZH1km_WRF_wWRF_uhelicity_only(EXP, 950, f"{h}:{m:02d}", folders, 'd02')
                plot_WRF_var_only(EXP, 850, f"{h}:{m:02d}", folders, 'd02')
                plot_WRF_var_only(EXP, 950, f"{h}:{m:02d}", folders, 'd02')
                if 'WSM6' in EXP:	 
                    plot_WRF_intqx(EXP, f"{h}:{m:02d}",6, folders, 'd02')#,'yakaira')
                elif 'P3' in EXP:
                    plot_WRF_intqx(EXP, f"{h}:{m:02d}",54, folders, 'd02')#,'yakaira')
        #for m in range(0, 60, 30):            
        #    plot_ZH1km_WRF_wWRFwind(EXP, f"23:{m:02d}", 'd02', folders)
        #    plot_ZH1km_WRF_wWRF_shelicity(EXP, f"23:{m:02d}", folders, 'd02')
        #    plot_ZH1km_WRF(EXP, f"23:{m:02d}", 'd02', 'yakaira', folders)   #plot_ZH1km_WRF(EXP, title, domain, servidor):
        #    plot_ZH1km_WRF_wWRF_uhelicity_only(EXP, f"23:{m:02d}", folders, 'd02')
        #    plot_WRF_var_only(EXP, f"23:{m:02d}", folders, 'd02')
        #    if 'WSM6' in EXP:	 
        #        plot_WRF_intqx(EXP, f"23:{m:02d}",6, folders, 'd02')#,'yakaira')
        #    elif 'P3' in EXP:
        #        plot_WRF_intqx(EXP, f"23:{m:02d}",54, folders, 'd02')#,'yakaira')
        #            

    #EXPs = ['P3_3MOM_LF_domain3_NoahMP_highres']
    #for EXP in EXPs:
    #    for h in range(21, 23):
    #        for m in range(0, 60, 30):
    #            plot_ZH1km_WRF_wWRFwind(EXP, f"{h}:{m:02d}", 'd02', folders)
    #            plot_ZH1km_WRF_wWRF_shelicity(EXP, f"{h}:{m:02d}", folders, 'd02')
    #            plot_ZH1km_WRF(EXP, f"{h}:{m:02d}", 'd02', 'yakaira', folders)   #plot_ZH1km_WRF(EXP, title, domain, servidor):
    #            plot_ZH1km_WRF_wWRF_uhelicity_only(EXP, f"{h}:{m:02d}", folders, 'd02')
    #            plot_WRF_var_only(EXP, f"{h}:{m:02d}", folders, 'd02')
    #             if 'WSM6' in EXP:	 
    #                plot_WRF_intqx(EXP, f"{h}:{m:02d}",6, folders, 'd02')#,'yakaira')
    #            elif 'P3' in EXP:
    #                plot_WRF_intqx(EXP, f"{h}:{m:02d}",54, folders, 'd02')#,'yakaira')

    return


#------------------------------------------------------------------------------
def run_all1(): 

    folders=config_folders_final.config_folders('yakaira')
    
    EXPs = ['initcond_fromwrf_domain3_WSM6_d01P3_54']
    for EXP in EXPs:
        for h in range(15, 22):
            for m in range(0, 60, 30):
                plot_ZH1km_WRF_wWRFwind(EXP, f"{h}:{m:02d}", 'd01', folders)
                plot_ZH1km_WRF_wWRF_shelicity(EXP, f"{h}:{m:02d}", folders, 'd01')
                plot_ZH1km_WRF(EXP, f"{h}:{m:02d}", 'd01', 'yakaira', folders)   #plot_ZH1km_WRF(EXP, title, domain, servidor):
                plot_ZH1km_WRF_wWRF_uhelicity_only(EXP, f"{h}:{m:02d}", folders, 'd01')
                plot_WRF_var_only(EXP, f"{h}:{m:02d}", folders, 'd01')
                plot_WRF_intqx(EXP, f"{h}:{m:02d}",54, folders, 'd01')#,'yakaira')


    #for h in range(18, 22):
    #    for m in range(0, 60, 30):        
    #        plot_WRF_diffvar_only('WSM6_domain3_NoahMP','P3_3MOM_LF_domain3_NoahMP', f"{h}:{m:02d}")
    #        plot_WRF_var_only('WSM6_domain3_NoahMP', f"{h}:{m:02d}")
    #        plot_WRF_var_only('P3_3MOM_LF_domain3_NoahMP', f"{h}:{m:02d}")
    #
    #for h in range(18, 22):
    #    for m in range(0, 60, 30):    
    #        plot_WRF_TyTd_only('WSM6_domain3_NoahMP', f"{h}:{m:02d}")
    #        plot_WRF_TyTd_only('P3_3MOM_LF_domain3_NoahMP', f"{h}:{m:02d}")
    #        
    #plot_WRF_hovmoller_thetae('WSM6_domain3_NoahMP')
    #plot_WRF_hovmoller_thetae('P3_3MOM_LF_domain3_NoahMP')
    #plot_WRF_hovmoller_thetae('THOM_domain3_NoahMP')
   
    return

#------------------------------------------------------------------------------
def run_1213case():

    folders=config_folders_final.config_folders('yakaira')
    
    EXP =  '1312_WSM6check' # '2501_WSM6check' 
    date =  '2018-12-13' # '2019-01-25'

    for h in range(18, 23):
        for m in range(0, 60, 30):
            plot_ZH1km_WRFdate(EXP, f"{h}:{m:02d}", 'd02', 'yakaira', folders, date)   #plot_ZH1km_WRF(EXP, title, domain, servidor):
    plot_ZH1km_WRFdate(EXP, "23:00", 'd02', 'yakaira', folders, date)   #plot_ZH1km_WRF(EXP, title, domain, servidor):
    plot_ZH1km_WRFdate(EXP, "23:30", 'd02', 'yakaira', folders, date)   #plot_ZH1km_WRF(EXP, title, domain, servidor):
    plot_ZH1km_WRFdate(EXP, "00:00", 'd02', 'yakaira', folders, '2018-12-14')   #plot_ZH1km_WRF(EXP, title, domain, servidor):
    plot_ZH1km_WRFdate(EXP, "00:30", 'd02', 'yakaira', folders, '2018-12-14')   #plot_ZH1km_WRF(EXP, title, domain, servidor):

                
    return

#------------------------------------------------------------------------------
def run_obs_radar():

    folders=config_folders_final.config_folders('yakaira')
 
    # Evolution of RMA1 to understund evolution of supercell  every 30 min 
    plot_general_radar_evolution(radar_folder=folders['rma1_dir'], title='RMA1', save_dir_compare=folders['save_dir_compare'], elev=3)

    #PLOT CSARP2
    plot_radar_cspr2_singletime([-33,-31.5], [-65.5,-63.5])

    #PLOT ALL RMA1:
    radar_folder     = folders['rma1_dir']
    file_list    = sorted(glob.glob(radar_folder+'*01.nc'))
    prefix = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RMA1/cfrad.20181110_'
    start_index = len(prefix)
    for filename in file_list:
        time        = filename[start_index:start_index+4]
        plot_radar_singletime(filename, time, 3, [-35,-31], [-65.5,-62], colmax=0, folders=folders)
        plot_radar_singletime(filename, time, 3, [-35,-31], [-65.5,-62], colmax=1, folders=folders)
    
 
    return

#-----------------------------------------------------------------
def rhop3():
    
    data = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFOUT/P3_3MOM_LF_domain3_NoahMP/'
    file = 'wrfout_d02_2018-11-10_20:30:00'
    
    ncfile       = Dataset(data+file,'r')        
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")
    rho          = wrf.getvar( ncfile,"rho_ice")
    z            = wrf.getvar( ncfile,"z") 
    zh           = wrf.getvar(ncfile, "REFL_10CM")
    REFL_10CM    = wrf.interplevel(zh, z, 3000)    
    
    arr_masked = np.where(rho == 0, np.nan, rho)
    
    #--------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12,8)) 
    pcm = ax.pcolormesh(lon[:,320], z[:,:,320]/1e3,  arr_masked[:,:,320], cmap='viridis', vmin=0,  vmax=900)
    cbar = plt.colorbar(pcm, ax=ax, shrink=1)
    cbar.cmap.set_under('white')    
    ax.set_title('rho_ice WRF-P3 model output')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Height (km)')
    ax.grid(True)
    ax.set_ylim([0, 20])
    
    fig.savefig('/home/vito.galligani/rho_p3_poster.png', dpi=300,transparent=False, bbox_inches='tight')


    return

#----------------------------------------------------------------
def run_cnrm():
    
    folders=config_folders.config_folders('cnrm')

    #EXP = 'WSM6_domain3_NoahMP'; domain = 'd02'
    EXP = 'initcond_fromwrf_domain3_WSM6_d01P3_54_test2'; domain = 'd01'    
    for h in range(18, 23):
        for m in range(0, 60, 30): 
            plot_ZH1km_WRF(EXP, f"{h}:{m:02d}",domain,'cnrm')
            
    return

#-----------------------------------------------------------------
run_all(['WSM6_domain3','WSM6_domain3_plevels']) 
#run_obs_radar()
#plot_VELradar_cspr2_singletime([-33.4,-31], [-65.5,-63.5])
