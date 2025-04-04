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
import matplotlib

plt.matplotlib.rc('font', family='serif', size = 12)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12  



#------------------------------------------------------------------------------
def figure_stad_map(fig, ax, domain, time, params, subtitle, server, WRFfolder):
    
    
    #------ READ WRF variables of interest ------------------------------------
    prefix   = 'wrfout_'+domain+'_2018-12-13_'+time+':00'
    filename = os.path.join(WRFfolder, prefix)
    ncfile       = Dataset(filename,'r')        
    z            = wrf.getvar( ncfile,"z") 
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")
    uvmet10      = wrf.getvar( ncfile, "uvmet10",units="kt")   # Default is ‘m s-1’.
    wspd_wdir10  = wrf.getvar( ncfile, "wspd_wdir10") 
    p            = wrf.getvar( ncfile, "pressure")
    theta        = wrf.getvar( ncfile, "theta", units="K")   # potential temperature 
    t2           = wrf.getvar( ncfile, "T2", timeidx=-1)     # 2-m temperature
    td_2m        = wrf.getvar( ncfile, "td2")
    ter          = wrf.getvar( ncfile, 'ter', units="m")
    u            = wrf.getvar( ncfile, 'ua', units="kt")
    v            = wrf.getvar( ncfile, 'va', units="kt")
    q            = wrf.getvar( ncfile, 'QVAPOR')
    # other calcls
    ter_3d    = np.tile(ter.values,[41,1,1])
    z.values  = z.values-ter.values
    u2        = wrf.interplevel(u, z, 6000)
    v2        = wrf.interplevel(v, z, 6000)
    ufld      = uvmet10.isel(u_v=0)
    vfld      = uvmet10.isel(u_v=1)
    u2.values = u2.values-ufld.values
    v2.values = v2.values-vfld.values
    
    dbz        = wrf.getvar(ncfile, 'dbz')        
    REFL_10CM  = wrf.interplevel(dbz, z, 1000)
    REFL_10CM.values = np.ma.masked_less(REFL_10CM.values,10.)
        
    sounding_parameters = wrf.smooth2d(wrf.getvar( ncfile, 'cape_2d'),8)


    #--------------------------------------------------------------------------
    if ('td_2m' in params['VAR']):
        VAR=td_2m
    
    elif ('t2' in params['VAR']):
        VAR=t2-273
    
    elif ('CAPE' in params['VAR']):
        VAR=sounding_parameters.isel(mcape_mcin_lcl_lfc=0)
    
    elif ('Updraft_Helicity'in params['VAR']):
        VAR = -1.*wrf.getvar(ncfile, 'updraft_helicity')
    
    elif ('0-1km_Storm_Relative_Helicity' in params['VAR']):
        VAR = (wrf.getvar( ncfile, 'helicity', top=1000.0))*-1

    elif ('0-3km_Storm_Relative_Helicity' in params['VAR']):
        VAR = (wrf.getvar( ncfile, 'helicity', top=3000.0))*-1
        
        
    elif ('1km_Radar_Reflectivity' in params['VAR']):
        dbz        = wrf.getvar(ncfile, 'dbz')        
        z.values   = z.values-ter.values
        VAR        = wrf.interplevel(dbz, z, 1000)
        VAR.values = np.ma.masked_less(VAR.values,5.)
    
    elif ('potentialTemp850' in params['VAR']):
        eth = wrf.g_temp.get_eth(ncfile)
        VAR = wrf.interplevel(eth, p, 850)
    
    elif ('convergence' in params['VAR']):
        qfld = wrf.interplevel(q, p, 850)
        grad_q_x,grad_q_y = np.gradient(qfld.values)
        grad_u_x,grad_u_y = np.gradient(ufld.values)
        grad_v_x,grad_v_y = np.gradient(vfld.values)
        uvmet = wrf.getvar(ncfile, 'wspd_wdir', units="kt")        
        dx=uvmet.XLAT[1,0]-uvmet.XLAT[0,0]
        dx=dx.values*111100    
        MFC_advect=-1.* (ufld.values*grad_q_x/dx)+(vfld.values*grad_q_y/dx)
        MFC_conv=-1.*qfld.values*((grad_u_x/dx)+(grad_v_y/dx))
    
        VAR = qfld
        VAR.values = -1.*wrf.smooth2d(86400.*1000.*(MFC_advect + MFC_conv),30)
    
    elif ('10mwinds' in params['VAR']):
        VAR = wspd_wdir10[0,:,:]
        
    #--------------------------------------------------------------------------

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
    
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,220)   
        
    if ('convergence' in params['VAR']):
        colors1 = plt.cm.YlOrRd_r(np.linspace(0, 0.8, 120))    
        colors2 = plt.cm.PiYG(np.linspace(0.8, 0.7, 60))
        colors = np.vstack((colors1,colors2))#, colors2, colors1))
        cmap_conv= mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
        cmap_conv.set_over('white')
        cmap_conv.set_under('red')
        
    else:
        pcm     = ax.contourf(lon, lat,  VAR,  cmap=params['cmap'], levels=params['levels'], 
                          vmin=params['vmin'], vmax=params['vmax'], extend='both')
        cbar    = plt.colorbar(pcm, ax=ax, shrink=1, extend='both')
        cbar.cmap.set_under('white')
        
    ax.plot(prov[:,0],prov[:,1],color='k');     
    ax.set_xlim([-65.5,-62]); 
    ax.set_ylim([-34.5,-31])
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
    
    # agrego vertical cross section selection:
    #ax.plot(lon_points, lat_points, '-r', linewidth=1.2)

    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2)    
    ax.set_title(subtitle+' ('+params['plottitle']+')')
    ax.contour(lon, lat, REFL_10CM, levels=[40], colors=['darkred'], linewidths=1.2)
    
    return fig

#------------------------------------------------------------------------------
def plot_all_WRF_variables(time, domain, server):
    
    #EXPs    = ['WSM6_domain3_NoahMP', 'P3_3MOM_LF_domain3_NoahMP', 'initcond_fromwrf_domain3_WSM6_d01P3_54_test2']
    #domains = ['d02','d02','d01']
    #EXPtitle = ['WSM6', 'P3' ,'P3(WSM6init)']
    
    EXPs    = ['WSM6', 'P3mp54', ]
    domains = ['d02','d02']
    EXPtitle = ['WSM6', 'P3']
            
    #------------------
    colors1 = plt.cm.BrBG(np.linspace(0., 0.4, 120))
    colors2 = plt.cm.YlOrRd(np.linspace(0.0, 0.55, 70))
    colors3 = plt.cm.Greens(np.linspace(0.1, 0.7, 80))
    colors4 = plt.cm.Blues(np.linspace(0.2, 0.9, 80))
    colors = np.vstack((colors1,colors2,colors3,colors4))#, colors2, colors1))
    cmap_td= mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    cmap_td.set_over((0,0,0.69))
    cmap_td.set_under((0.2,0.1,0))
    
    import colormaps_new as cmapss
    colors1=np.array(cmapss._HomeyerRainbow_data[0:80])
    colors2=np.array(cmapss._HomeyerRainbow_data[100:230])
    colors3=np.array(cmapss._viridis_data[0:50])
    colors = np.vstack((colors3,colors1,colors2))
    cmap_temp = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
     
    # CAPE
    colors1 = plt.cm.YlOrRd(np.linspace(0.2, 0.65, 50))
    colors2 = plt.cm.Reds(np.linspace(0.6, 1.0, 64))
    colors3 = plt.cm.Purples(np.linspace(0.5, 0.9, 64))
    colors = np.vstack((colors1,colors2,colors3))#, colors2, colors1))
    cmap_cape= mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    cmap_cape.set_over((0.3,0,1))
    cmap_cape.set_under((1,1,1))
    
    #------------------
    folders   = config_folders.config_folders('cnrm_1312')

    save_dir_compare = folders['save_dir_compare']

    # #-------------------------------------------------------------------------
    # # Convergence
    # #-------------------------------------------------------------------------     
    colors1 = plt.cm.YlOrRd_r(np.linspace(0, 0.8, 120))    
    colors2 = plt.cm.PiYG(np.linspace(0.8, 0.7, 60))
    colors = np.vstack((colors1,colors2))#, colors2, colors1))
    cmap_conv= mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    cmap_conv.set_over('white')
    cmap_conv.set_under('red')
    params = { 'VAR': 'convergence', 
                'cmap': cmap_conv,
                'levels': np.arange(-600,20,100.), 
                'vmin': -600,
                'vmax': 0,
                'plottitle': '850 hPa moisture convergence '+time, 
                'filename': save_dir_compare+'/Comparison/Convergence_'+'Comparison_'+time+'.png' }    
    
    fig, ax = plt.subplots(nrows=1, ncols=3, constrained_layout=True, figsize=(8*3,8)) 
    counter=0
    for EXPi in EXPs:
        WRFfolder = folders[EXPi]
        fig = figure_stad_map(fig, ax[counter], domains[counter], time, params, EXPtitle[counter], server, WRFfolder)
        counter=counter+1
    fig.savefig(params['filename'], dpi=300,transparent=False,bbox_inches='tight')        
    plt.close(fig)    


    # #-------------------------------------------------------------------------
    # # Plot Td
    # #-------------------------------------------------------------------------
    params = { 'VAR': 'td_2m', 
               'cmap': cmap_td,
               'levels': np.arange(-12,28,2), 
               'vmin': -12,
               'vmax': 28,
               'plottitle': 'Td 2m '+time, 
               'filename': save_dir_compare+'/Comparison/td2m_'+'Comparison_'+time+'.png' }

    fig, ax = plt.subplots(nrows=1, ncols=3, constrained_layout=True, figsize=(8*3,8)) 
    counter=0
    for EXPi in EXPs:
        WRFfolder = folders[EXPi]
        fig = figure_stad_map(fig, ax[counter], domains[counter], time, params, EXPtitle[counter], server, WRFfolder)
        counter=counter+1
    fig.savefig(params['filename'], dpi=300,transparent=False,bbox_inches='tight')
    plt.close(fig)
    
    # #-------------------------------------------------------------------------
    # # Plot T
    # #-------------------------------------------------------------------------    
    params = { 'VAR': 't2', 
                'cmap': cmap_temp,
                'levels': np.arange(0,42,2), 
                'vmin': 0,
                'vmax': 40,
                'plottitle': 'T 2m '+time, 
                'filename': save_dir_compare+'/Comparison/t2m_'+'Comparison_'+time+'.png' }    
    
    fig, ax = plt.subplots(nrows=1, ncols=3, constrained_layout=True, figsize=(8*3,8)) 
    counter=0
    for EXPi in EXPs:
        WRFfolder = folders[EXPi]
        fig = figure_stad_map(fig, ax[counter], domains[counter], time, params, EXPtitle[counter], server, WRFfolder)
        counter=counter+1
    fig.savefig(params['filename'], dpi=300,transparent=False,bbox_inches='tight')
    plt.close(fig)
    
    # #-------------------------------------------------------------------------
    # # 0-6 km MLCAPE
    # #-------------------------------------------------------------------------        
    # params = { 'VAR': 'CAPE', 
    #             'cmap': cmap_cape,
    #             'levels': [100,500,1000,1500,2000,2500,3000,3500,4000], 
    #             'vmin': 100,
    #             'vmax': 4000,
    #             'plottitle': 'MLCAPE 0-6 km (J kg-1) '+time, 
    #             'filename': save_dir_compare+'/Comparison/mlCAPE_'+'Comparison_'+time+'.png' }    
    
    # fig, ax = plt.subplots(nrows=1, ncols=3, constrained_layout=True, figsize=(8*3,8)) 
    # counter=0
    # for EXPi in EXPs:
    #     WRFfolder = folders[EXPi]
    #     fig = figure_stad_map(fig, ax[counter], domains[counter], time, params, EXPtitle[counter], server, WRFfolder)
    #     counter=counter+1
    # fig.savefig(params['filename'], dpi=300,transparent=False,bbox_inches='tight')        
    # plt.close(fig)

    # #-------------------------------------------------------------------------
    # # Updraft_Helicity
    # #-------------------------------------------------------------------------  
    # params = { 'VAR': 'Updraft_Helicity', 
    #             'cmap':  get_cmap("gist_stern_r"),
    #             'levels':  np.arange(0,80,10), 
    #             'vmin': -10,
    #             'vmax': 80,
    #             'plottitle': 'Updraft Helicity (m2 s-2) '+time,  
    #             'filename': save_dir_compare+'/Comparison/Updraft_Helicity_'+'Comparison_'+time+'.png' }    
    
    # fig, ax = plt.subplots(nrows=1, ncols=3, constrained_layout=True, figsize=(8*3,8)) 
    # counter=0
    # for EXPi in EXPs:
    #     WRFfolder = folders[EXPi]
    #     fig = figure_stad_map(fig, ax[counter], domains[counter], time, params, EXPtitle[counter], server, WRFfolder)
    #     counter=counter+1
    # fig.savefig(params['filename'], dpi=300,transparent=False,bbox_inches='tight')        
    # plt.close(fig)
    
    #-------------------------------------------------------------------------
    # Storm relative Helicity
    #-------------------------------------------------------------------------  
    # params = { 'VAR': '0-1km_Storm_Relative_Helicity', 
    #             'cmap':  get_cmap("cubehelix_r"),
    #             'levels': np.arange(50,750,50), 
    #             'vmin': 50,
    #             'vmax': 750,
    #             'plottitle': '0-1 km AGL storm relative helicity (m2s-2) '+time, 
    #             'filename': save_dir_compare+'/Comparison/1km_Storm_Relative_Helicity_'+'Comparison_'+time+'.png' }    
    
    # fig, ax = plt.subplots(nrows=1, ncols=3, constrained_layout=True, figsize=(8*3,8)) 
    # counter=0
    # for EXPi in EXPs:
    #     WRFfolder = folders[EXPi]
    #     fig = figure_stad_map(fig, ax[counter], domains[counter], time, params, EXPtitle[counter], server, WRFfolder)
    #     counter=counter+1
    # fig.savefig(params['filename'], dpi=300,transparent=False,bbox_inches='tight')        
    # plt.close(fig)
    
    params = { 'VAR': '0-3km_Storm_Relative_Helicity', 
                 'cmap':  get_cmap("cubehelix_r"),
                 'levels': np.arange(50,950,50), 
                 'vmin': 50,
                 'vmax': 950,
                 'plottitle': '0-3 km AGL storm relative helicity (m2s-2) '+time, 
                 'filename': save_dir_compare+'/Comparison/3km_Storm_Relative_Helicity_'+'Comparison_'+time+'.png' }    
    
    fig, ax = plt.subplots(nrows=1, ncols=3, constrained_layout=True, figsize=(8*3,8)) 
    counter=0
    for EXPi in EXPs:
        WRFfolder = folders[EXPi]
        fig = figure_stad_map(fig, ax[counter], domains[counter], time, params, EXPtitle[counter], server, WRFfolder)
        counter=counter+1
    fig.savefig(params['filename'], dpi=300,transparent=False,bbox_inches='tight')        
    plt.close(fig)    
    
    # #-------------------------------------------------------------------------
    # # 1km Zh
    # #-------------------------------------------------------------------------  
    params = { 'VAR': '0-1km_Radar_Reflectivity', 
                'cmap':  get_cmap("gist_ncar"),
                'levels': np.arange(5,75,5.), 
                'vmin': 5,
                'vmax': 75,
                'plottitle': '1 km AGL radar reflectivity '+time, 
                'filename': save_dir_compare+'/Comparison/1kRadarReflectivity_'+'Comparison_'+time+'.png' }    
    
    fig, ax = plt.subplots(nrows=1, ncols=3, constrained_layout=True, figsize=(8*3,8)) 
    counter=0
    for EXPi in EXPs:
        WRFfolder = folders[EXPi]
        fig = figure_stad_map(fig, ax[counter], domains[counter], time, params, EXPtitle[counter], server, WRFfolder)
        counter=counter+1
    fig.savefig(params['filename'], dpi=300,transparent=False,bbox_inches='tight')        
    plt.close(fig)    
    
    # #-------------------------------------------------------------------------
    # # 850hPa Equiv. Potential temperature
    # #-------------------------------------------------------------------------      
    params = { 'VAR': 'potentialTemp850', 
                'cmap': matplotlib.cm.get_cmap("viridis_r"),
                'levels': np.arange(300,370,10.), 
                'vmin': 300,
                'vmax': 370,
                'plottitle': 'Equiv. Potential Temp. 850 hPa '+time, 
                'filename': save_dir_compare+'/Comparison/EquivPotentialTemp_'+'Comparison_'+time+'.png' }    
    
    fig, ax = plt.subplots(nrows=1, ncols=3, constrained_layout=True, figsize=(8*3,8)) 
    counter=0
    for EXPi in EXPs:
        WRFfolder = folders[EXPi]
        fig = figure_stad_map(fig, ax[counter], domains[counter], time, params, EXPtitle[counter], server, WRFfolder)
        counter=counter+1
    fig.savefig(params['filename'], dpi=300,transparent=False,bbox_inches='tight')        
    plt.close(fig)    

    # #-------------------------------------------------------------------------
    # # Winds
    # #-------------------------------------------------------------------------  
    params = { 'VAR': '10mwinds', 
                'cmap': matplotlib.cm.get_cmap("YlOrBr"),
                'levels': np.linspace(4.0, 24.0, 6), 
                'vmin': 4,
                'vmax': 24,
                'plottitle': '10m wind '+time, 
                'filename': save_dir_compare+'/Comparison/Winds_'+'Comparison_'+time+'.png' }    
    
    fig, ax = plt.subplots(nrows=1, ncols=3, constrained_layout=True, figsize=(8*3,8)) 
    counter=0
    for EXPi in EXPs:
        WRFfolder = folders[EXPi]
        fig = figure_stad_map(fig, ax[counter], domains[counter], time, params, EXPtitle[counter], server, WRFfolder)
        counter=counter+1
    fig.savefig(params['filename'], dpi=300,transparent=False,bbox_inches='tight')        
    plt.close(fig)  



    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
for h in range(17, 23):
    for m in range(0, 60, 30): 
          plot_all_WRF_variables(f"{h}:{m:02d}", 'd02', 'cnrm')

plot_all_WRF_variables('23:00', 'd02', 'cnrm')
plot_all_WRF_variables('23:30', 'd02', 'cnrm')


