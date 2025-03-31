import pyart 
import pyart 
import numpy as np 
import package_functions as vito_functions
import wrf
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib import ticker
import glob
import os 
from netCDF4 import Dataset
from PIL import Image
import config_folders
import netCDF4 as nc
from shapely.geometry import Polygon
import matplotlib

#------------------------------------------------------------------------------
def colormaps(variable):
    """
    Choose colormap for a radar variable

    variable : str
       Radar variable to define which colormap to use. (e.g. ref,
       dv, zdr..) More to come

    returns : matplotlib ColorList
       A Color List for the specified *variable*.

    """
    import matplotlib.colors as colors

    # Definicion de las disntitas paletas de colores:
    nws_ref_colors =([ [ 1.0/255.0, 159.0/255.0, 244.0/255.0],
                    [  3.0/255.0,   0.0/255.0, 244.0/255.0],
                    [  2.0/255.0, 253.0/255.0,   2.0/255.0],
                    [  1.0/255.0, 197.0/255.0,   1.0/255.0],
                    [  0.0/255.0, 142.0/255.0,   0.0/255.0],
                    [253.0/255.0, 248.0/255.0,   2.0/255.0],
                    [229.0/255.0, 188.0/255.0,   0.0/255.0],
                    [253.0/255.0, 149.0/255.0,   0.0/255.0],
                    [253.0/255.0,   0.0/255.0,   0.0/255.0],
                    #[212.0/255.0,   0.0/255.0,   0.0/255.0, 0.4],
                    [188.0/255.0,   0.0/255.0,   0.0/255.0],
                    [248.0/255.0,   0.0/255.0, 253.0/255.0],
                    [152.0/255.0,  84.0/255.0, 198.0/255.0]
                    ])
    # Definicion de las disntitas paletas de colores:
    nws_ref_colors_transparent =([ [ 1.0/255.0, 159.0/255.0, 244.0/255.0, 0.3],
                    [  3.0/255.0,   0.0/255.0, 244.0/255.0, 0.3],
                    [  2.0/255.0, 253.0/255.0,   2.0/255.0, 0.3],
                    [  1.0/255.0, 197.0/255.0,   1.0/255.0, 0.3],
                    [  0.0/255.0, 142.0/255.0,   0.0/255.0, 0.3],
                    [253.0/255.0, 248.0/255.0,   2.0/255.0, 0.3],
                    [229.0/255.0, 188.0/255.0,   0.0/255.0, 1],
                    [253.0/255.0, 149.0/255.0,   0.0/255.0, 1],
                    [253.0/255.0,   0.0/255.0,   0.0/255.0, 1],
                    #[212.0/255.0,   0.0/255.0,   0.0/255.0, 0.4],
                    [188.0/255.0,   0.0/255.0,   0.0/255.0, 1],
                    [248.0/255.0,   0.0/255.0, 253.0/255.0, 1],
                    [152.0/255.0,  84.0/255.0, 198.0/255.0, 1]
                    ])
    
    nws_zdr_colors = ([ [  1.0/255.0, 159.0/255.0, 244.0/255.0],
                    [  3.0/255.0,   0.0/255.0, 244.0/255.0],
                    [  2.0/255.0, 253.0/255.0,   2.0/255.0],
                    [  1.0/255.0, 197.0/255.0,   1.0/255.0],
                    [  0.0/255.0, 142.0/255.0,   0.0/255.0],
                    [253.0/255.0, 248.0/255.0,   2.0/255.0],
                    [229.0/255.0, 188.0/255.0,   2.0/255.0],
                    [253.0/255.0, 149.0/255.0,   0.0/255.0],
                    [253.0/255.0,   0.0/255.0,   0.0/255.0],
                    [188.0/255.0,   0.0/255.0,   0.0/255.0],
                    [152.0/255.0,  84.0/255.0, 198.0/255.0]
                    ])

    nws_dv_colors = ([  [0,  1,  1],
                    [0,  0.966666638851166,  1],
                    [0,  0.933333337306976,  1],
                    [0,  0.899999976158142,  1],
                    [0,  0.866666674613953,  1],
                    [0,  0.833333313465118,  1],
                    [0,  0.800000011920929,  1],
                    [0,  0.766666650772095,  1],
                    [0,  0.733333349227905,  1],
                    [0,  0.699999988079071,  1],
                    [0,  0.666666686534882,  1],
                    [0,  0.633333325386047,  1],
                    [0,  0.600000023841858,  1],
                    [0,  0.566666662693024,  1],
                    [0,  0.533333361148834,  1],
                    [0,  0.5,  1],
                    [0,  0.466666668653488,  1],
                    [0,  0.433333337306976,  1],
                    [0,  0.400000005960464,  1],
                    [0,  0.366666674613953,  1],
                    [0,  0.333333343267441,  1],
                    [0,  0.300000011920929,  1],
                    [0,  0.266666680574417,  1],
                    [0,  0.233333334326744,  1],
                    [0,  0.200000002980232,  1],
                    [0,  0.16666667163372,   1],
                    [0,  0.133333340287209,  1],
                    [0,  0.100000001490116,  1],
                    [0,  0.0666666701436043, 1],
                    [0,  0.0333333350718021, 1],
                    [0,  0,  1],
                    [0,  0,  0],
                    [0,  0,  0],
                    [0,  0,  0],
                    [0,  0,  0],
                    [1,  0,  0],
                    [1,  0.0322580635547638, 0],
                    [1,  0.0645161271095276, 0],
                    [1,  0.0967741906642914, 0],
                    [1,  0.129032254219055,  0],
                    [1,  0.161290317773819,  0],
                    [1,  0.193548381328583,  0],
                    [1,  0.225806444883347,  0],
                    [1,  0.25806450843811,   0],
                    [1,  0.290322571992874,  0],
                    [1,  0.322580635547638,  0],
                    [1,  0.354838699102402,  0],
                    [1,  0.387096762657166,  0],
                    [1,  0.419354826211929,  0],
                    [1,  0.451612889766693,  0],
                    [1,  0.483870953321457,  0],
                    [1,  0.516129016876221,  0],
                    [1,  0.548387110233307,  0],
                    [1,  0.580645143985748,  0],
                    [1,  0.612903237342834,  0],
                    [1,  0.645161271095276,  0],
                    [1,  0.677419364452362,  0],
                    [1,  0.709677398204803,  0],
                    [1,  0.74193549156189,   0],
                    [1,  0.774193525314331,  0],
                    [1,  0.806451618671417,  0],
                    [1,  0.838709652423859,  0],
                    [1,  0.870967745780945,  0],
                    [1,  0.903225779533386,  0],
                    [1,  0.935483872890472,  0],
                    [1,  0.967741906642914,  0],
                    [1,  1,  0]   ])


    cmap_nws_ref = colors.ListedColormap(nws_ref_colors)
    cmap_nws_zdr = colors.ListedColormap(nws_zdr_colors)
    cmap_nws_dv = colors.ListedColormap(nws_dv_colors)
    cmap_nws_ref_trans = colors.ListedColormap(nws_ref_colors_transparent)
    

    if variable == 'ref':
       return cmap_nws_ref
    if variable == 'ref2':
       return cmap_nws_ref_trans
        

    if variable == 'zdr':
       return cmap_nws_zdr

    if variable == 'dv':
       return cmap_nws_dv

#------------------------------------------------------------------------------
def pyplot_rings(lat_radar,lon_radar,radius):
    """
    Calculate lat-lon of the maximum range ring

    lat_radar : float32
       Radar latitude. Positive is north.

    lon_radar : float32
       Radar longitude. Positive is east.

    radius : float32
       Radar range in kilometers.

    returns : numpy.array
       A 2d array containing the 'radius' range latitudes (lat) and longitudes (lon)

    """
    import numpy as np

    R=12742./2.
    m=2.*np.pi*R/360.
    alfa=np.arange(-np.pi,np.pi,0.0001)

    nazim  = 360.0
    nbins  = 480.0
    binres = 0.5

    #azimuth = np.transpose(np.tile(np.arange(0,nazim,1), (int(nbins),1)))
    #rangos  = np.tile(np.arange(0,nbins,1)*binres, (int(nazim),1))
    #lats    = lat_radar + (rangos/m)*np.cos((azimuth)*np.pi/180.0)
    #lons    = lon_radar + (rangos/m)*np.sin((azimuth)*np.pi/180.0)/np.cos(lats*np.pi/180)

    lat_radius = lat_radar + (radius/m)*np.sin(alfa)
    lon_radius = lon_radar + ((radius/m)*np.cos(alfa)/np.cos(lat_radius*np.pi/180))

    return lat_radius, lon_radius

#------------------------------------------------------------------------------





#------------------------------------------------------------------------------
def plot_common_transect_MAP(EXP, EXP_P3, time, DOW6file, DOW7file, 
                             CSAPR2file, VAR, x_wsm6, y_wsm6, 
                             x_P3_3mom_LF, y_P3_3mom_LF):
    

    folders = config_folders.config_folders('yakaira')
    save_dir_compare = folders['save_dir_compare']
    
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    
    radar         = pyart.io.read(DOW6file) 
    radarLAT_DOW6 = radar.latitude['data'][0]
    radarLON_DOW6 = radar.longitude['data'][0]
    DOW6_range = round(np.nanmax(radar.range['data'])/1000)

    radar         = pyart.io.read(DOW7file) 
    radarLAT_DOW7 = radar.latitude['data'][0]
    radarLON_DOW7 = radar.longitude['data'][0]
    DOW7_range = round(np.nanmax(radar.range['data'])/1000)
    
    radar         = pyart.io.read(CSAPR2file)     
    radarLAT_CSPR2 = radar.latitude['data'][0]
    radarLON_CSPR2 = radar.longitude['data'][0]
    CSAPR2_range = round(np.nanmax(radar.range['data'])/1000)
    
    
    ncfile1 = Dataset( os.path.join(folders[EXP], 'wrfout_d02_2018-11-10_'+time+':00') ,'r')
    ncfile2 = Dataset( os.path.join(folders[EXP_P3], 'wrfout_d02_2018-11-10_'+time+':00') ,'r')
    
    
    [qr_WSM6, qs_WSM6, qi_WSM6, qc_WSM6, qg_WSM6, 
     qr_int_WSM6, qs_int_WSM6, qi_int_WSM6, qc_int_WSM6, qg_int_WSM6] = vito_functions.get_q_ints6(ncfile1)
    qi_sum_WSM6 = qs_int_WSM6 + qi_int_WSM6 + qg_int_WSM6

    #[qr_P3_3mom_LF, qi_P3_3mom_LF, qc_P3_3mom_LF, 
    # qr_int_P3_3mom_LF, qi_int_P3_3mom_LF, qc_int_P3_3mom_LF, qir, qib, qil] = vito_functions.get_q_ints3(ncfile2)
        
    
    lat_WSM6      =  wrf.getvar( ncfile1,"lat") 
    lon_WSM6      =  wrf.getvar( ncfile1,"lon")

    #lat_P3_3mom_LF      =  wrf.getvar( ncfile2,"lat") 
    #lon_P3_3mom_LF      =  wrf.getvar( ncfile2,"lon")    
    
    
    prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')
    
    fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
    
    
    
    fig, axes = plt.subplots(nrows=1, ncols=2, constrained_layout=True,figsize=[16,8])

    if 'QICE' in VAR:
        plot_data_WSM6 = qi_sum_WSM6
        plot_data_P3_3mom_LF = qi_int_P3_3mom_LF
        vminn = 0
        vmaxx = 12
        cmap  = get_cmap("viridis")
        plt.suptitle( 'q_i (tot) at ' + title + 'UTC')

    elif 'MAXRADAR_WRF' in VAR:
        data1 = wrf.getvar(ncfile1, "REFL_10CM")
        plot_data_WSM6 = np.nanmax(data1, axis=0)
        
        data2 = wrf.getvar(ncfile2, "REFL_10CM")  # (REFL_10CM) 
        plot_data_P3_3mom_LF = np.nanmax(data2, axis=0)
        vminn = 0
        vmaxx = 70
        cmap  = colormaps('ref')
        plt.suptitle( 'ZMAX (dbZ) at ' + time + 'UTC')
    
    elif 'RADAR_WRF_level4' in VAR:
        data1 = wrf.getvar(ncfile1, "REFL_10CM")
        plot_data_WSM6 = data1[4,:,:]
        
        data2 = wrf.getvar(ncfile2, "REFL_10CM")  # (REFL_10CM) 
        plot_data_P3_3mom_LF = data2[4,:,:]
        vminn = 0
        vmaxx = 70
        cmap  = colormaps('ref')
        plt.suptitle( 'Z(level 4) (dbZ) at ' + time + 'UTC')
        
    elif 'interp1km' in VAR:
        zh             = wrf.getvar(ncfile1, "REFL_10CM")
        z              = wrf.getvar(ncfile1, "z") 
        plot_data_WSM6 = wrf.interplevel(zh, z, 1000)

        zh2           = wrf.getvar(ncfile2, "REFL_10CM")  # (REFL_10CM) 
        z2            = wrf.getvar(ncfile2, "z") 
        data2         = wrf.interplevel(z2, z2, 1000)
        plot_data_P3_3mom_LF = data2.copy()
        vminn = 0
        vmaxx = 70
        cmap  = colormaps('ref')
        plt.suptitle( 'Z(1km) (dbZ) at ' + time + 'UTC')
        
        
    pcm0 = axes[0].pcolormesh(lon_WSM6, lat_WSM6, plot_data_WSM6, cmap=cmap, vmin=vminn,  vmax=vmaxx)
    cbar = plt.colorbar(pcm0, ax=axes[0], shrink=1)
    cbar.cmap.set_under('white')
    axes[0].grid()
    
    # agrego contorno de 500 y 1000m
    axes[0].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
    
    #pcm1 = axes[1].pcolormesh(lon_P3_3mom_LF, lat_P3_3mom_LF, plot_data_P3_3mom_LF, cmap=cmap, vmin=vminn,  vmax=vmaxx)

    #axes[1].plot(np.nan, np.nan, 'darkred',  linewidth=1.5, label='45dBZ cont.')
    #axes[1].legend(fontsize=9, loc='upper left')

    #cbar = plt.colorbar(pcm1, ax=axes[1], shrink=1)
    #cbar.cmap.set_under('white')
    
    for i in range(2):
        axes[i].grid()
        axes[i].set_xlim([-65.5,-62])
        axes[i].set_ylim([-33.5,-31.3])

    axes[0].set_title(EXP)                                                  
    axes[1].set_title(EXP_P3)                                                        


    # RMA1 
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)
    for iax in range(2): 
        axes[iax].plot(lon_radius, lat_radius, 'k', linewidth=0.8, label='RMA1 (120km)')
    # DOW6
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_DOW6,radarLON_DOW6,DOW6_range)
    for iax in range(2): 
        axes[iax].plot(lon_radius, lat_radius, 'darkred', linewidth=0.8, label='DOW6 ('+str(DOW6_range)+'km)')
    # DOW7
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_DOW7,radarLON_DOW7,DOW7_range)
    for iax in range(2): 
        axes[iax].plot(lon_radius, lat_radius, 'magenta', linewidth=0.8, label='DOW7 ('+str(DOW7_range)+'km)')
    # CSPR2
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_CSPR2,radarLON_CSPR2,CSAPR2_range)
    for iax in range(2): 
        axes[iax].plot(lon_radius, lat_radius, 'darkgreen', linewidth=0.8,  label='CSAPR-2 ('+str(CSAPR2_range)+'km)')
    for iax in range(2): 
        plt.plot(np.nan, np.nan, '-or' , linewidth=1.2, label= 'Transect')
    
    axes[0].legend(fontsize=10, loc='lower left')
     
    
    #
    xx = np.vstack([x_wsm6[[0]],x_wsm6[[1]]])
    yy = np.vstack([y_wsm6[[0]],y_wsm6[[1]]])
    axes[0].plot(xx, yy, '-or' , linewidth=1.2)
    #
    xx = np.vstack([x_P3_3mom_LF[[0]],x_P3_3mom_LF[[1]]])
    yy = np.vstack([y_P3_3mom_LF[[0]],y_P3_3mom_LF[[1]]])
    axes[1].plot(xx, yy, '-or' , linewidth=1.2)
    #
    fig.savefig(save_dir_compare+'/'+EXP+'_'+time+VAR+'_compare_contourmap.png', dpi=300,transparent=False)
    plt.show() 
    
    return

#------------------------------------------------------------------------------
def plot_common_transect_P3micro_ICE_derived(ncfile, x, y, title, savedir, fig, axes, ncol):
    
    # voy a querer qice, Nice, mui, Dm, DMax, rho
    
    # Need to interpolate first to common z_level grid 
    z_interp = 1e3*np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 
                 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20])
    
    Re       = 6.3781e6
    
    temp     = wrf.g_temp.get_tk(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL)
    z        = Re*geopo_p/(Re-geopo_p)

    # Create the start point and end point for the cross section for P3_3MOM_LF
    start_point = wrf.CoordPair(lat=y[0], lon=x[0])
    end_point   = wrf.CoordPair(lat=y[1], lon=x[1])

    [qr, qi, qc, qr_int, qi_int, qc_int, qir, qib, qil] = vito_functions.get_q_ints3(ncfile)
    QNICE   = wrf.getvar(ncfile, "QNICE")       # 'Ice Number concentration' ( kg-1)

    QICEKGKG = wrf.getvar(ncfile, "QICE") # kg/kg
    ZICE = wrf.getvar(ncfile, "QZI")
    RHOM = wrf.getvar(ncfile, "RHO_ICE")  # mass-weighted mean ice density cat 1
    DICE = wrf.getvar(ncfile, "D_ICE")    # mass-weighted mean ice size cat 1
    VICE = wrf.getvar(ncfile, "V_ICE")    # mass-weighted mean ice fallspedd cat 1

    VAR1 = (6*QICEKGKG)/((3.14*RHOM))**2
    mu = (QNICE/ZICE) / VAR1   
    
    REFL_10CM = wrf.getvar(ncfile, "REFL_10CM")

    # Initial guess for mu (you may need to adjust this based on your specific problem)
    #mu_guess = 5.0
    
    #mu = fsolve(mu_solveequation, mu_guess, args=((QNICE/ZICE) / VAR1 ,))


    # Compute the vertical cross-section interpolation.  Also, include the lat/lon points along the cross-section.
    qr_cross   = wrf.vertcross(qr, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True)
    qi_cross = wrf.vertcross(qi, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True)
    nice_cross   = wrf.vertcross(QNICE, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    mu_cross   = wrf.vertcross(mu, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    rho_cross   = wrf.vertcross(RHOM, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    v_cross   = wrf.vertcross(VICE, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    d_cross   = wrf.vertcross(DICE, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    Z10_cross   = wrf.vertcross(REFL_10CM, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 

    contour_x = np.arange(0, qi_cross.shape[-1], 1)
    contour_y = wrf.to_np(qi_cross.coords["vertical"])

    vminn = np.floor(np.log10(0.00002))-1; vmaxx = np.ceil(np.log10(0.005))+1; 
    lev_exp_qice = np.logspace(vminn, vmaxx, num=int(vmaxx-vminn+1))        
    lev_exp_ice   = np.logspace(-3, 10, num=10-(-3)+1)
    
    # MU SHOULD BE BETWEEN 0 AND 20? 
    
    contours_map_0 = axes[0,ncol].contourf(contour_x, contour_y, wrf.to_np(qr_cross), lev_exp_qice, locator=ticker.LogLocator(), cmap=get_cmap("viridis")) #, vmin=0.0002, vmax=0.01)
    contours_map_1 = axes[1,ncol].contourf(contour_x, contour_y, wrf.to_np(qi_cross), lev_exp_qice, locator=ticker.LogLocator(), cmap=get_cmap("viridis")) #, vmin=0.0002, vmax=0.01)
    contours_map_2 = axes[2,ncol].contourf(contour_x, contour_y, wrf.to_np(nice_cross), lev_exp_ice , locator=ticker.LogLocator(),  cmap=get_cmap("viridis")) #, vmin=0.00002, vmax=0.05)
    contours_map_3 = axes[3,ncol].contourf(contour_x, contour_y, wrf.to_np(mu_cross), cmap=get_cmap("viridis"))
    contours_map_4 = axes[4,ncol].contourf(contour_x, contour_y, wrf.to_np(rho_cross), cmap=get_cmap("viridis")) #, vmin=0.0002, vmax=0.005)
    contours_map_5 = axes[5,ncol].contourf(contour_x, contour_y, wrf.to_np(v_cross),   cmap=get_cmap("viridis"))
    contours_map_6 = axes[6,ncol].contourf(contour_x, contour_y, wrf.to_np(d_cross),   cmap=get_cmap("viridis"))
    contours_map_7 = axes[7,ncol].contourf(contour_x, contour_y, wrf.to_np(Z10_cross),   cmap=colormaps('ref'))
    
    axes[0,ncol].set_title('qr (kg/m3)')
    axes[1,ncol].set_title('qi (kg/m3)')
    axes[2,ncol].set_title('qnice (1/m3)')
    axes[3,ncol].set_title(r'$/mu$ice')
    axes[4,ncol].set_title(r'$\rho$ice (kg/m3)')
    axes[5,ncol].set_title('vice (m/s)')
    axes[6,ncol].set_title('dm (m)')
    axes[7,ncol].set_title('Z (dBZ)')
    
    cbar =plt.colorbar(contours_map_0, ax=axes[0,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_1, ax=axes[1,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_2, ax=axes[2,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_3, ax=axes[3,ncol]); cbar.cmap.set_under('white'); 
    cbar =plt.colorbar(contours_map_4, ax=axes[4,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_5, ax=axes[5,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_6, ax=axes[6,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_7, ax=axes[7,ncol]); cbar.cmap.set_under('white')

    #--------------------------------------------------------------------------     # now add 0C level    
    [output_array_T, output_array_Z] = vito_functions.find_00(temp, 273, z)
    HGT = wrf.getvar(ncfile, "HGT", timeidx=-1)
    HGT.data = output_array_Z
    isoline_hgt = wrf.interpline(HGT, wrfin=ncfile, start_point=start_point, end_point=end_point)
    xs = np.arange(0, qi_cross.shape[-1], 1)
    for ii in range(8):
        axes[ii, ncol].plot(xs, wrf.to_np(isoline_hgt), '-r', linewidth=1.2)
        axes[ii, ncol].set_ylabel("Height (m)", fontsize=12)

    #--------------------------------------------------------------------------
    # Set the X-ticks to be LATITUDE.
    coord_pairs = wrf.to_np(qi_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [round(pair.lat,2) for pair in wrf.to_np(coord_pairs)]
    for ii in range(8):
        axes[ii,ncol].set_xticks(x_ticks[::10])
        axes[ii,ncol].set_xticklabels(x_labels[::10],fontsize=12)                                 
        axes[ii,ncol].set_xlabel("Latitude", fontsize=12)    

    # Set the X-ticks to be LONGITUDE 
    coord_pairs = wrf.to_np(qi_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [round(pair.lon,2) for pair in wrf.to_np(coord_pairs)]
    for ii in range(8):
        axes[ii,ncol].set_xticks(x_ticks[::20])
        axes[ii,ncol].set_xticklabels(x_labels[::20],fontsize=12)   
        axes[ii,ncol].set_xlabel("Longitude", fontsize=12)  
        
    return 

#------------------------------------------------------------------------------
def return_qxs_WSM6(ncfile, domainlons, domainlats, z_interp):

    # Need to interpolate first to common z_level grid 
    z_interp =  [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 
                 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20]    
    
    temp     = wrf.g_temp.get_tk(ncfile)
    pressure = wrf.g_pressure.get_pressure(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL
    Re       = 6.3781e6
    z_level  = Re*geopo_p/(Re-geopo_p)
    lats      =  wrf.getvar( ncfile,"lat") 
    lons      =  wrf.getvar( ncfile,"lon")

    qr = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QRAIN"][0,:,:,:]  ),  "ght_msl", z_interp)          
    qs = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QSNOW"][0,:,:,:]  ),  "ght_msl", z_interp)          
    qi = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QICE"][0,:,:,:]   ),  "ght_msl", z_interp)          
    qc = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QCLOUD"][0,:,:,:] ),  "ght_msl", z_interp)         
    qg = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QGRAUP"][0,:,:,:] ),  "ght_msl", z_interp)   
    
    
    # Create a mask for the domain
    lat_mask = (lats >= domainlats[0]) & (lats <= domainlats[1])
    lon_mask = (lons >= domainlons[0]) & (lons <= domainlons[1])
    mask = lat_mask & lon_mask 
    
    # Apply the mask to the variable q
    qr_masked = np.where(mask, qr, np.nan)  # Replace values outside the domain with NaN
    qi_masked = np.where(mask, qi, np.nan)  # Replace values outside the domain with NaN
    qc_masked = np.where(mask, qc, np.nan)  # Replace values outside the domain with NaN
    qs_masked = np.where(mask, qs, np.nan)  # Replace values outside the domain with NaN
    qg_masked = np.where(mask, qg, np.nan)  # Replace values outside the domain with NaN
    z_level_masked = np.where(mask, z_level, np.nan)    
    
    return qr_masked, qs_masked, qi_masked, qc_masked, qg_masked, qg_masked+qi_masked+qs_masked, z_level_masked

#------------------------------------------------------------------------------
def return_qxs_p3(ncfile, domainlons, domainlats, z_interp):
    
    temp     = wrf.g_temp.get_tk(ncfile)
    pressure = wrf.g_pressure.get_pressure(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL
    Re       = 6.3781e6
    z_level  = Re*geopo_p/(Re-geopo_p)
    lats      =  wrf.getvar( ncfile,"lat") 
    lons      =  wrf.getvar( ncfile,"lon")

    qr = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QRAIN"][0,:,:,:]), "ght_msl", z_interp)   
    qi = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QICE"][0,:,:,:]),  "ght_msl", z_interp)    
    qc = wrf.vinterp(ncfile, np.squeeze(ncfile.variables["QCLOUD"][0,:,:,:]), "ght_msl", z_interp)   
    
    # Create a mask for the domain
    lat_mask = (lats >= domainlats[0]) & (lats <= domainlats[1])
    lon_mask = (lons >= domainlons[0]) & (lons <= domainlons[1])
    
    mask = lat_mask & lon_mask
    
    # Apply the mask to the variable q
    qr_masked = np.where(mask, qr, np.nan)  # Replace values outside the domain with NaN
    qi_masked = np.where(mask, qi, np.nan)  # Replace values outside the domain with NaN
    qc_masked = np.where(mask, qc, np.nan)  # Replace values outside the domain with NaN
    z_level_masked = np.where(mask, z_level, np.nan)
    
    return qr_masked, qi_masked, qc_masked, z_level_masked

#------------------------------------------------------------------------------
def plot_domain_qxs(wrf_file_wsm6, wrf_file_p3, title, savedir):
    
    # Domain    
    domainlons = [-65.5,-62]
    domainlats = [-33.5,-31.3] 
    
    # Need to interpolate first to common z_level grid 
    z_interp =  [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 
                 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20] 
    
    [qr, qs, qi, qc, qg, qitot, z_lev] = return_qxs_WSM6(wrf_file_wsm6, domainlons, domainlats, z_interp)
    [qr_p3, qi_p3, qc_p3, z_lev_p3] = return_qxs_p3(wrf_file_p3, domainlons, domainlats, z_interp)
    
    
    # Calculate the average inside the domain, ignoring NaN values
    #average_q = np.nanmean(q_masked)
   
    fig, axes = plt.subplots(nrows=1, ncols=3, constrained_layout=True,figsize=[12,6])
    axes[0].plot(np.nanmean(qr*1000, axis = (1, 2)), z_interp, color='darkred', linewidth=1.2, label='qr')    
    axes[0].plot(np.nanmean(qs*1000, axis = (1, 2)), z_interp, color='darkblue', linewidth=1.2, label='qs')    
    axes[0].plot(np.nanmean(qg*1000, axis = (1, 2)), z_interp, color='darkgreen', linewidth=1.2, label='qg')    
    axes[0].plot(np.nanmean(qi*1000, axis = (1, 2)), z_interp, color='indigo', linewidth=1.2, label='qi')    
    axes[0].plot(np.nanmean(qc*1000, axis = (1, 2)), z_interp, color='dodgerblue', linewidth=1.2, label='qc')    
    axes[0].plot(np.nanmean(qitot*1000, axis = (1, 2)), z_interp, color='k', linewidth=1.2, label='qg+qs+qi')  
    axes[0].set_xlabel('Mixing ratio (g/kg)')
    axes[0].set_ylabel('Height (km)')
    axes[0].set_title('WSM6')
    axes[0].legend()
    axes[0].grid(True)
    axes[0].set_xlim([0, 0.4])
    axes[0].set_ylim([0,20])

    axes[1].plot(np.nanmean(qr_p3*1000, axis = (1, 2)), z_interp, color='darkred', linewidth=1.2, label='qr')    
    axes[1].plot(np.nanmean(qi_p3*1000, axis = (1, 2)), z_interp, color='k', linewidth=1.2, label='qi (tot)')    
    axes[1].plot(np.nanmean(qc_p3*1000, axis = (1, 2)), z_interp, color='dodgerblue', linewidth=1.2, label='qc')    
    axes[1].set_xlabel('Mixing ratio (g/kg)')
    axes[1].set_ylabel('Height (km)')
    axes[1].set_title('P3 3MOM LF')
    axes[1].legend()
    axes[1].set_xlim([0, 0.35])
    axes[1].grid(True)
    axes[1].set_ylim([0,20])
                
    plt.suptitle('Domain average vertical profiles at' +title)    
    plt.show()    
    fig.savefig(savedir+'compareDomain_qxs'+title+'.png', dpi=300,transparent=False)

    return



    






#------------------------------------------------------------------------------
def plot_common_transect_WSM6micro_ICE_derived(EXP, x, y, time):
    
    # voy a querer qice, Nice, mui, Dm, DMax, rho
    
    fig, axes = plt.subplots(nrows=6, ncols=2, constrained_layout=True,figsize=[16,30])
    ncol = 0

    folders = config_folders.config_folders('yakaira')
    save_dir_compare = folders['save_dir_compare']
    
    ncfile = Dataset( os.path.join(folders[EXP], 'wrfout_d02_2018-11-10_'+time+':00') ,'r')
    
    # Need to interpolate first to common z_level grid 
    z_interp = 1e3*np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 
                 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20])
    
    Re       = 6.3781e6
    
    temp     = wrf.g_temp.get_tk(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL)
    z        = Re*geopo_p/(Re-geopo_p)

    # Create the start point and end point for the cross section for P3_3MOM_LF
    start_point = wrf.CoordPair(lat=y[0], lon=x[0])
    end_point   = wrf.CoordPair(lat=y[1], lon=x[1])

    [qr, qs, qi, qc, qg, qr_int, qs_int, qi_int, qc_int, qg_int] = vito_functions.get_q_ints6(ncfile)
    QICEKGKG = wrf.getvar(ncfile, "QICE") # kg/kg
    REFL_10CM = wrf.getvar(ncfile, "REFL_10CM")
    
    # Compute the vertical cross-section interpolation.  Also, include the lat/lon points along the cross-section.
    qr_cross   = wrf.vertcross(qr, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True)
    qi_cross   = wrf.vertcross(qi, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True)
    qs_cross   = wrf.vertcross(qs, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    qg_cross   = wrf.vertcross(qg, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    qitot_cross = wrf.vertcross( (qg+qs+qi), z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    Z10_cross   = wrf.vertcross( REFL_10CM, z, levels=z_interp, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, meta=True) 
    
    contour_x = np.arange(0, qi_cross.shape[-1], 1)
    contour_y = wrf.to_np(qi_cross.coords["vertical"])

    vminn = np.floor(np.log10(0.00002))-1; vmaxx = np.ceil(np.log10(0.005))+1; 
    lev_exp_qice = np.logspace(vminn, vmaxx, num=int(vmaxx-vminn+1))        
    
    # MU SHOULD BE BETWEEN 0 AND 20? 
    
    contours_map_0 = axes[0,ncol].contourf(contour_x, contour_y, wrf.to_np(qr_cross), lev_exp_qice, locator=ticker.LogLocator(), cmap=get_cmap("viridis")) #, vmin=0.0002, vmax=0.01)
    contours_map_1 = axes[1,ncol].contourf(contour_x, contour_y, wrf.to_np(qi_cross), lev_exp_qice, locator=ticker.LogLocator(), cmap=get_cmap("viridis")) #, vmin=0.0002, vmax=0.01)
    contours_map_2 = axes[2,ncol].contourf(contour_x, contour_y, wrf.to_np(qs_cross), lev_exp_qice , locator=ticker.LogLocator(),  cmap=get_cmap("viridis")) #, vmin=0.00002, vmax=0.05)
    contours_map_3 = axes[3,ncol].contourf(contour_x, contour_y, wrf.to_np(qg_cross), cmap=get_cmap("viridis"))
    contours_map_4 = axes[4,ncol].contourf(contour_x, contour_y, wrf.to_np(qitot_cross), cmap=get_cmap("viridis")) #, vmin=0.0002, vmax=0.005)
    contours_map_5 = axes[5,ncol].contourf(contour_x, contour_y, wrf.to_np(Z10_cross),   cmap=colormaps('ref'))
    
    axes[0,ncol].set_title('qr (kg/m3)')
    axes[1,ncol].set_title('qi (kg/m3)')
    axes[2,ncol].set_title('qs (kg/m3)')
    axes[3,ncol].set_title('qg (kg/m3)')
    axes[4,ncol].set_title('qs+qg+qi (kg/m3)')
    axes[5,ncol].set_title('Z (dBZ)')

    cbar =plt.colorbar(contours_map_0, ax=axes[0,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_1, ax=axes[1,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_2, ax=axes[2,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_3, ax=axes[3,ncol]); cbar.cmap.set_under('white'); 
    cbar =plt.colorbar(contours_map_4, ax=axes[4,ncol]); cbar.cmap.set_under('white')
    cbar =plt.colorbar(contours_map_5, ax=axes[5,ncol]); cbar.cmap.set_under('white')

    #--------------------------------------------------------------------------     # now add 0C level    
    [output_array_T, output_array_Z] = vito_functions.find_00(temp, 273, z)
    HGT = wrf.getvar(ncfile, "HGT", timeidx=-1)
    HGT.data = output_array_Z
    isoline_hgt = wrf.interpline(HGT, wrfin=ncfile, start_point=start_point, end_point=end_point)
    xs = np.arange(0, qi_cross.shape[-1], 1)
    for ii in range(6):
        axes[ii, ncol].plot(xs, wrf.to_np(isoline_hgt), '-r', linewidth=1.2)
        axes[ii, ncol].set_ylabel("Height (m)", fontsize=12)

    #--------------------------------------------------------------------------
    # Set the X-ticks to be LATITUDE.
    coord_pairs = wrf.to_np(qi_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [round(pair.lat,2) for pair in wrf.to_np(coord_pairs)]
    for ii in range(6):
        axes[ii,ncol].set_xticks(x_ticks[::10])
        axes[ii,ncol].set_xticklabels(x_labels[::10],fontsize=12)                                 
        axes[ii,ncol].set_xlabel("Latitude", fontsize=12)    

    # Set the X-ticks to be LONGITUDE 
    coord_pairs = wrf.to_np(qi_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [round(pair.lon,2) for pair in wrf.to_np(coord_pairs)]
    for ii in range(6):
        axes[ii,ncol].set_xticks(x_ticks[::20])
        axes[ii,ncol].set_xticklabels(x_labels[::20],fontsize=12)   
        axes[ii,ncol].set_xlabel("Longitude", fontsize=12)  
        
    fig.savefig(save_dir_compare+'/'+EXP+'_'+time+'_compare_contourmap.png', dpi=300,transparent=False)
        
        
    return 

#------------------------------------------------------------------------------
def plot_MAXRADAR_WRF_evolution(WRFfolder, iHH, title, save_dir_compare, other_WRFfolder, OTHER_NAME, elev):
    
    #start_time = (19*60) #(19*60)+20  #19:20
    #end_time = 20*60+40 #20:40    
    #start_time_1 = 1850
    #timereso = 10 #minutes
    
    start_time_1 = 1700
    start_time   = 17*60
    end_time = 20*60+30 #20:40    
    timereso = 30 #minutes
    
    
    # RMA1 
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    TH_name       = 'TH'
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)
    all_times_minutes = np.arange(start_time, end_time + timereso, timereso) 
    time1infile   = 82
    radar_folder =  '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RMA1/'
    wildfiles = 'cfrad.20181110*01.nc'  # try also with .02.nc maybe? 

    
    # RADAR PLOTS
    file_pattern = os.path.join(radar_folder, wildfiles)
    file_list    = sorted(glob.glob(file_pattern))
    timesr        = [filename[time1infile:time1infile+4] for filename in file_list]

    # Step 1: Round each time to the nearest 10 minutes
    rounded_minutes = [int(round(vito_functions.hhmm_to_minutes(t) / 10) * 10) for t in timesr]

    # Step 2: Create an array with NaN values for all intervals
    time_array = np.full(len(all_times_minutes), np.nan, dtype=object)
    filename_array = np.full(len(all_times_minutes), np.nan, dtype=object)
    
    # Step 3: Fill in times (and filenames) where data is available
    for i, rounded_time in enumerate(rounded_minutes):
        index = np.where(all_times_minutes == rounded_time)[0]
        if index.size > 0:
            time_array[index[0]] = vito_functions.minutes_to_hhmm(rounded_time)
            filename_array[index[0]] = file_list[i]

    # Convert all_times_minutes back to HHMM format for easy reading
    all_times_hhmm = [vito_functions.minutes_to_hhmm(m) for m in all_times_minutes]
    
    # for WSM6 iHH = 107
    # for P3_3MOM_LF iHH = 109

    all_files = sorted(glob.glob(os.path.join(WRFfolder, 'wrfout_d02_2018-11-10*')))

    # Filter filenames based on HH:MM > '19:30' and < 2040
    filtered_files = [
        file for file in all_files
        if start_time_1-10 < int(os.path.basename(file).split('_')[3].split(':')[0]) * 100 +
           int(os.path.basename(file).split('_')[3].split(':')[1]) < 2050 ]
        
    fig, axes = plt.subplots(nrows=3, ncols=4, constrained_layout=True,figsize=[16,16])
    counter = 0
    for ax, filename in zip(axes.flat, filtered_files):
        print(filename)
        times = filename[iHH:iHH+5]   
        ncfile       = Dataset(filename,'r')

        # Find the coincident radar name
        filename_radar = filename_array[counter]
        print(filename_radar)
        counter = counter +1 
        
        print(times)
        
        #---- select the actual observation
        radar       = pyart.io.read(filename_radar)         
        # Configure a gatefilter to filter out copolar correlation coefficient values > 0.9
        gatefilter = pyart.filters.GateFilter(radar)
        gatefilter.exclude_transition()
        gatefilter.exclude_equal('RHOHV', 0.9)
        
        #lons_,lats_,colmax = vito_functions.get_colmax(radar, TH_name, gatefilter)
        ZHelev18 = radar.get_field(elev, TH_name, copy=False)
        [lats_, lons_, _] = radar.get_gate_lat_lon_alt(elev, reset_gate_coords=False, filter_transitions=False)
        cont3 = ax.contour(lons_, lats_, ZHelev18, levels=[40], colors=['magenta'], linewidths=2)


        REFL_10CM = np.nanmax(wrf.getvar(ncfile, "REFL_10CM"), axis=0)
        lat       =  wrf.getvar( ncfile,"lat") 
        lon       =  wrf.getvar( ncfile,"lon")

        pcm = ax.pcolormesh(lon, lat,  REFL_10CM, cmap=colormaps('ref'), vmin=0,  vmax=70)
        cbar = plt.colorbar(pcm, ax=ax, shrink=1)
        cbar.cmap.set_under('white')
        ax.grid()
        ax.set_xlim([-65.5,-62])
        ax.set_ylim([-33.5,-31.3])
        ax.set_title(title+' '+times)      
        ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
        cont = ax.contour(lon, lat, REFL_10CM, levels=[50], colors=['darkblue'], linewidths=1.2)

    
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
        
        
    ax.plot(np.nan, np.nan, 'k', linewidth=0.8, label='RMA1 (120 km)')
    ax.plot(np.nan, np.nan, 'darkblue', linewidth=1.2, label=f'50 dBZ {OTHER_NAME}')
    ax.plot(np.nan, np.nan, 'magenta', linewidth=2, label=f'50 dBZ RMA1')
    ax.legend()
        
    plt.show()
    fig.savefig(save_dir_compare+'MAXZWRF_evolution_'+title+'long.png', dpi=300,transparent=False)
    
    return#------------------------------------------------------------------------------
    def plot_MAXRADAR_WRF_evolution(WRFfolder, iHH, title, save_dir_compare, other_WRFfolder, OTHER_NAME, elev):
        
        #start_time = (19*60) #(19*60)+20  #19:20
        #end_time = 20*60+40 #20:40    
        #start_time_1 = 1850
        #timereso = 10 #minutes
        
        start_time_1 = 1700
        start_time   = 17*60
        end_time = 20*60+30 #20:40    
        timereso = 30 #minutes
        
        
        # RMA1 
        radarLAT_RMA1 = -31.441389
        radarLON_RMA1 = -64.191944
        TH_name       = 'TH'
        [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)
        all_times_minutes = np.arange(start_time, end_time + timereso, timereso) 
        time1infile   = 82
        radar_folder =  '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RMA1/'
        wildfiles = 'cfrad.20181110*01.nc'  # try also with .02.nc maybe? 

        
        # RADAR PLOTS
        file_pattern = os.path.join(radar_folder, wildfiles)
        file_list    = sorted(glob.glob(file_pattern))
        timesr        = [filename[time1infile:time1infile+4] for filename in file_list]

        # Step 1: Round each time to the nearest 10 minutes
        rounded_minutes = [int(round(vito_functions.hhmm_to_minutes(t) / 10) * 10) for t in timesr]

        # Step 2: Create an array with NaN values for all intervals
        time_array = np.full(len(all_times_minutes), np.nan, dtype=object)
        filename_array = np.full(len(all_times_minutes), np.nan, dtype=object)
        
        # Step 3: Fill in times (and filenames) where data is available
        for i, rounded_time in enumerate(rounded_minutes):
            index = np.where(all_times_minutes == rounded_time)[0]
            if index.size > 0:
                time_array[index[0]] = vito_functions.minutes_to_hhmm(rounded_time)
                filename_array[index[0]] = file_list[i]

        # Convert all_times_minutes back to HHMM format for easy reading
        all_times_hhmm = [vito_functions.minutes_to_hhmm(m) for m in all_times_minutes]
        
        # for WSM6 iHH = 107
        # for P3_3MOM_LF iHH = 109

        all_files = sorted(glob.glob(os.path.join(WRFfolder, 'wrfout_d02_2018-11-10*')))

        # Filter filenames based on HH:MM > '19:30' and < 2040
        filtered_files = [
            file for file in all_files
            if start_time_1-10 < int(os.path.basename(file).split('_')[3].split(':')[0]) * 100 +
               int(os.path.basename(file).split('_')[3].split(':')[1]) < 2050 ]
            
        fig, axes = plt.subplots(nrows=3, ncols=4, constrained_layout=True,figsize=[16,16])
        counter = 0
        for ax, filename in zip(axes.flat, filtered_files):
            print(filename)
            times = filename[iHH:iHH+5]   
            ncfile       = Dataset(filename,'r')

            # Find the coincident radar name
            filename_radar = filename_array[counter]
            print(filename_radar)
            counter = counter +1 
            
            print(times)
            
            #---- select the actual observation
            radar       = pyart.io.read(filename_radar)         
            # Configure a gatefilter to filter out copolar correlation coefficient values > 0.9
            gatefilter = pyart.filters.GateFilter(radar)
            gatefilter.exclude_transition()
            gatefilter.exclude_equal('RHOHV', 0.9)
            
            #lons_,lats_,colmax = vito_functions.get_colmax(radar, TH_name, gatefilter)
            ZHelev18 = radar.get_field(elev, TH_name, copy=False)
            [lats_, lons_, _] = radar.get_gate_lat_lon_alt(elev, reset_gate_coords=False, filter_transitions=False)
            cont3 = ax.contour(lons_, lats_, ZHelev18, levels=[40], colors=['magenta'], linewidths=2)


            REFL_10CM = np.nanmax(wrf.getvar(ncfile, "REFL_10CM"), axis=0)
            lat       =  wrf.getvar( ncfile,"lat") 
            lon       =  wrf.getvar( ncfile,"lon")

            pcm = ax.pcolormesh(lon, lat,  REFL_10CM, cmap=colormaps('ref'), vmin=0,  vmax=70)
            cbar = plt.colorbar(pcm, ax=ax, shrink=1, label='Zhmax WRF WSM6')
            cbar.cmap.set_under('white')
            ax.grid()
            ax.set_xlim([-65.5,-62])
            ax.set_ylim([-33.5,-31.3])
            ax.set_title(title+' '+times)      
            ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
            cont = ax.contour(lon, lat, REFL_10CM, levels=[50], colors=['darkblue'], linewidths=1.2)

        
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
            
            
        ax.plot(np.nan, np.nan, 'k', linewidth=0.8, label='RMA1 (120 km)')
        ax.plot(np.nan, np.nan, 'darkblue', linewidth=1.2, label=f'50 dBZ {OTHER_NAME}')
        ax.plot(np.nan, np.nan, 'magenta', linewidth=2, label=f'50 dBZ RMA1')
        ax.legend()
            
        plt.show()
        fig.savefig(save_dir_compare+'MAXZWRF_evolution_'+title+'long.png', dpi=300,transparent=False)
        
        return
    
#------------------------------------------------------------------------------
def plot_MAXRADAR_WRF_evolution_radar(WRFfolder, iHH, title, save_dir_compare, other_WRFfolder, OTHER_NAME, elev):
    
    # same as above but plot radar reflectivity and WRF contour 
    #start_time = (19*60) #(19*60)+20  #19:20
    #end_time = 20*60+40 #20:40    
    #start_time_1 = 1850
    #timereso = 10 #minutes
    
    start_time_1 = 1700
    start_time   = 17*60
    end_time = 20*60+30 #20:40    
    timereso = 30 #minutes
    
    
    # RMA1 
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    TH_name       = 'TH'
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)
    all_times_minutes = np.arange(start_time, end_time + timereso, timereso) 
    time1infile   = 82
    radar_folder =  '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RMA1/'
    wildfiles = 'cfrad.20181110*01.nc'  # try also with .02.nc maybe? 

    
    # RADAR PLOTS
    file_pattern = os.path.join(radar_folder, wildfiles)
    file_list    = sorted(glob.glob(file_pattern))
    timesr        = [filename[time1infile:time1infile+4] for filename in file_list]

    # Step 1: Round each time to the nearest 10 minutes
    rounded_minutes = [int(round(vito_functions.hhmm_to_minutes(t) / 10) * 10) for t in timesr]

    # Step 2: Create an array with NaN values for all intervals
    time_array = np.full(len(all_times_minutes), np.nan, dtype=object)
    filename_array = np.full(len(all_times_minutes), np.nan, dtype=object)
    
    # Step 3: Fill in times (and filenames) where data is available
    for i, rounded_time in enumerate(rounded_minutes):
        index = np.where(all_times_minutes == rounded_time)[0]
        if index.size > 0:
            time_array[index[0]] = vito_functions.minutes_to_hhmm(rounded_time)
            filename_array[index[0]] = file_list[i]

    # Convert all_times_minutes back to HHMM format for easy reading
    all_times_hhmm = [vito_functions.minutes_to_hhmm(m) for m in all_times_minutes]
    
    # for WSM6 iHH = 107
    # for P3_3MOM_LF iHH = 109

    all_files = sorted(glob.glob(os.path.join(WRFfolder, 'wrfout_d02_2018-11-10*')))

    # Filter filenames based on HH:MM > '19:30' and < 2040
    filtered_files = [
        file for file in all_files
        if start_time_1-10 < int(os.path.basename(file).split('_')[3].split(':')[0]) * 100 +
           int(os.path.basename(file).split('_')[3].split(':')[1]) < 2050 ]
        
    fig, axes = plt.subplots(nrows=3, ncols=4, constrained_layout=True,figsize=[16,16])
    counter = 0
    for ax, filename in zip(axes.flat, filtered_files):
        print(filename)
        times = filename[iHH:iHH+5]   
        ncfile       = Dataset(filename,'r')

        # Find the coincident radar name
        filename_radar = filename_array[counter]
        print(filename_radar)
        counter = counter +1 
        
        print(times)
        
        #---- select the actual observation
        radar       = pyart.io.read(filename_radar)         
        # Configure a gatefilter to filter out copolar correlation coefficient values > 0.9
        gatefilter = pyart.filters.GateFilter(radar)
        gatefilter.exclude_transition()
        gatefilter.exclude_equal('RHOHV', 0.9)
        
        #lons_,lats_,colmax = vito_functions.get_colmax(radar, TH_name, gatefilter)
        ZHelev18 = radar.get_field(elev, TH_name, copy=False)
        [lats_, lons_, _] = radar.get_gate_lat_lon_alt(elev, reset_gate_coords=False, filter_transitions=False)
        pcm = ax.pcolormesh(lons_, lats_, ZHelev18, cmap=colormaps('ref'), vmin=0,  vmax=70)
        cbar = plt.colorbar(pcm, ax=ax, shrink=1, label='Zh RMA1 elev 3')
        cbar.cmap.set_under('white')
        ax.grid()        
        
        REFL_10CM = np.nanmax(wrf.getvar(ncfile, "REFL_10CM"), axis=0)
        lat       =  wrf.getvar( ncfile,"lat") 
        lon       =  wrf.getvar( ncfile,"lon")
        cont3 = ax.contour(lon, lat, REFL_10CM, levels=[40], colors=['magenta'], linewidths=2)

        ax.set_xlim([-65.5,-62])
        ax.set_ylim([-33.5,-31.3])
        ax.set_title(title+' '+times)      
        ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)

    
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
        
        
    #ax.plot(np.nan, np.nan, 'k', linewidth=0.8, label='RMA1 (120 km)')
    #ax.plot(np.nan, np.nan, 'darkblue', linewidth=1.2, label=f'50 dBZ {OTHER_NAME}')
    #ax.plot(np.nan, np.nan, 'magenta', linewidth=2, label=f'50 dBZ RMA1')
    #ax.legend()
        
    plt.show()
    fig.savefig(save_dir_compare+'radarevolution_MAXZWRF_evolutioncont_'+title+'long.png', dpi=300,transparent=False)
    
    return

#------------------------------------------------------------------------------





def make_wrfout_gif(WRFfolder, iHH, title, output_gif_path):

    tempsavedir =  '/home/vito.galligani/Work/Studies/HAILCASE_10112018/Plots/tmp/'
    
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)
    
    
    start_time_1 = 1700

    all_files = sorted(glob.glob(os.path.join(WRFfolder, 'wrfout_d02_2018-11-10*')))

    # Filter filenames based on HH:MM > '19:30' and < 2040
    filtered_files = [
        file for file in all_files
        if start_time_1 < int(os.path.basename(file).split('_')[3].split(':')[0]) * 100 +
           int(os.path.basename(file).split('_')[3].split(':')[1]) < 2200 ]
    
    
    images = []
    for filename in filtered_files:
        print(filename)
        times = filename[iHH:iHH+5]   
        ncfile       = Dataset(filename,'r')

        REFL_10CM = np.nanmax(wrf.getvar(ncfile, "REFL_10CM"), axis=0)
        lat       =  wrf.getvar( ncfile,"lat") 
        lon       =  wrf.getvar( ncfile,"lon")

        fig = plt.figure(figsize=(8,8)) 
        pcm = plt.pcolormesh(lon, lat,  REFL_10CM, cmap=colormaps('ref'), vmin=0,  vmax=70)
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.cmap.set_under('white')
        plt.grid()
        plt.xlim([-65.5,-62])
        plt.ylim([-33.5,-31.3])
        plt.title(title+' '+times)      
        plt.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
        

        # Save the figure as a temporary image file
        temp_filename = f"{times}.png"
        plt.savefig(tempsavedir+title+temp_filename) 
        images.append(Image.open(tempsavedir+title+temp_filename))
        plt.close(fig)  # Close the figure to free memory 
        
    # Step 3: Save as a GIF
    output_gif_path = output_gif_path+'WRFOUT_'+title+'.gif'
    images[0].save(output_gif_path, save_all=True, append_images=images[1:], duration=1000, loop=0)
    print(f"GIF saved at {output_gif_path}")

    return


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------  
def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    import matplotlib.pyplot as plt
    import numpy as np

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    
    return base.from_list(cmap_name, color_list, N)

#-------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_psuedo_HID_RMA1(test_transect, time, rfiles, ZDRoffset, xlim_range1, xlim_range2, freezing_lev1, freezing_lev2):
    

    folders      = config_folders.config_folders('yakaira')
    radar_folder = folders['rma1_dir']
    filename     = os.path.join(radar_folder, rfiles)
    radar        = pyart.io.read(filename) 

    #- Radar sweep
    nelev       = 0
    start_index = radar.sweep_start_ray_index['data'][3]
    end_index   = radar.sweep_end_ray_index['data'][3]
    lats0        = radar.gate_latitude['data'][start_index:end_index]
    lons0        = radar.gate_longitude['data'][start_index:end_index]
    azimuths     = radar.azimuth['data'][start_index:end_index]
    
    Ze_transect     = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); Ze_transect[:]=np.nan
    ZDR_transect    = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); ZDR_transect[:]=np.nan
    PHI_transect    = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); PHI_transect[:]=np.nan
    lon_transect    = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); lon_transect[:]=np.nan
    lat_transect    = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); lat_transect[:]=np.nan
    RHO_transect    = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); RHO_transect[:]=np.nan
    approx_altitude = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); approx_altitude[:]=np.nan
    color           = np.full((  len(radar.sweep_start_ray_index['data']), lats0.shape[1], 4), np.nan)
    gate_range      = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); gate_range[:]=np.nan

    azydims = lats0.shape[1]-1
    window = 9
    
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
        
        start_index = radar.sweep_start_ray_index['data'][nlev]
        end_index   = radar.sweep_end_ray_index['data'][nlev]       
        
        ZHZH    = radar.fields['TH']['data'][start_index:end_index]
        TV      = radar.fields['TV']['data'][start_index:end_index]
        ZDRZDR  = (ZHZH-TV)-ZDRoffset   
        RHORHO  = radar.fields['RHOHV']['data'][start_index:end_index]       
        #ZDRZDR[RHORHO<0.75] = np.nan
        #RHORHO[RHORHO<0.75] = np.nan
        
        lats        = radar.gate_latitude['data'][start_index:end_index]
        lons        = radar.gate_longitude['data'][start_index:end_index]
        
        # En verdad buscar azimuth no transecta ... 
        azimuths       = radar.azimuth['data'][start_index:end_index]
        target_azimuth = azimuths[test_transect]  #- target azimuth for nlev=0 test case is 301.5
        filas          = np.asarray(abs(azimuths-target_azimuth)<=0.1).nonzero()
        lon_transect[nlev,:] = lons[filas,:]
        lat_transect[nlev,:] = lats[filas,:]
        #
        gateZ    = radar.gate_z['data'][start_index:end_index]
        gateX    = radar.gate_x['data'][start_index:end_index]
        gateY    = radar.gate_y['data'][start_index:end_index]
        gates_range  = np.sqrt(gateX**2 + gateY**2 + gateZ**2)
        #
        Ze_transect[nlev,:]      = pyart.correct.phase_proc.smooth_and_trim( np.ravel(ZHZH[filas,:]), window)
        ZDR_transect[nlev,:]     = pyart.correct.phase_proc.smooth_and_trim( np.ravel(ZDRZDR[filas,:]), window)
        RHO_transect[nlev,:]     = pyart.correct.phase_proc.smooth_and_trim( np.ravel(RHORHO[filas,:]), window)
        # 
        [xgate, ygate, zgate]   = pyart.core.antenna_to_cartesian(gates_range[filas,:]/1e3, azimuths[filas],radar.get_elevation(nlev)[0]);
        approx_altitude[nlev,:] = zgate/1e3
        gate_range[nlev,:]      = gates_range[filas,:]/1e3;
                
    
    #---------------------------------------- PLOT 

    #---- REFLECTIVITY
    fig = plt.figure(figsize=[15,11])
    fig.add_subplot(111)
    mycolorbar = plt.pcolormesh(lon_transect, approx_altitude, Ze_transect, cmap=colormaps('ref'), vmin=0, vmax=70)
    
    #- De esta manera me guardo el color con el que rellenar los polygons (scatter plot para sacar el color de cada pixel)
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
         fig = plt.figure(figsize=[30,10])
         fig.add_subplot(221)
         sc = plt.scatter(lon_transect[nlev,:], approx_altitude[nlev,:],
                 s=1,c=Ze_transect[nlev,:],
                 cmap=colormaps('ref'), vmin=0, vmax=70)
         color[nlev,:,:] = sc.to_rgba(Ze_transect[nlev,:])
         plt.close()
    
    
    #- Try polygons
    fig2, axes = plt.subplots(nrows=3,ncols=1,constrained_layout=True,figsize=[8,6])  # 8,4 muy chiquito
    
    fig1 = plt.figure(figsize=(15,20))
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
         # Create the cone for each elevation IN TERMS OF RANGE. 
         # ===> ACA HABRIA QUE AGREGAR COMO CAMBIA LA ALTURA CON EL RANGE (?)
         ancho_haz_i0    = (np.pi/180*gate_range[nlev,0]/2)
         ancho_haz_i1099 = (np.pi/180*gate_range[nlev,azydims]/2)
         P1 = Polygon([( gate_range[nlev,0],    approx_altitude[nlev,0]-ancho_haz_i0      ),
                   ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]-ancho_haz_i1099),
                   ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]+ancho_haz_i1099),
                   ( gate_range[nlev,0],    approx_altitude[nlev,0]+ancho_haz_i0      )])
         ancho = 50/1E3
         # Location of gates? Every 100m?  
         LS = [Polygon([(gate_range[nlev,x]-ancho, 0),
                   (gate_range[nlev,x]+ancho, 0),
                   (gate_range[nlev,x]+ancho, 50),
                   (gate_range[nlev,x]-ancho, 50)]) for x in np.arange(approx_altitude.shape[1])]
         # Plot
         for i, l in enumerate(LS):
             # Get the polygon of the intersection between the cone and the space 
             #reserved for a specific point
             inter = l.intersection(P1)
             x,y = inter.exterior.xy    
             # Then plot it, filled by the color we want
             axes[0].fill(x, y, color = color[nlev,i,:], )
             x, y = P1.exterior.xy
    axes[0].set_ylim([0, 15])
    axes[0].set_ylabel('Altitude (km)')
    axes[0].grid()
    axes[0].set_xlim((xlim_range1, xlim_range2))
    norm = matplotlib.colors.Normalize(vmin=0.,vmax=70.)
    cax = matplotlib.cm.ScalarMappable(norm=norm, cmap=colormaps('ref'))
    cax.set_array(Ze_transect)
    cbar_z = fig2.colorbar(cax, ax=axes[0], shrink=1.1, ticks=np.arange(0,70.01,10), label='Zh (dBZ)')
    axes[0].axhline(y=freezing_lev1,color='k',linestyle='--', linewidth=1.2)
    axes[0].axhline(y=freezing_lev2,color='k',linestyle='--', linewidth=1.2)

    del mycolorbar, x, y, inter
    #---------------------------------------- ZDR
    N = (5+2)
    cmap_ZDR = discrete_cmap(int(N), 'jet') 
    #- Simple pcolormesh plot! 
    fig = plt.figure(figsize=[15,11])
    fig.add_subplot(221)
    mycolorbar = plt.pcolormesh(lon_transect, approx_altitude,
                ZDR_transect,
                cmap=cmap_ZDR, vmin=-2, vmax=5.)
    plt.close()
    
    #- De esta manera me guardo el color con el que rellenar los polygons
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
        # scatter plot para sacar el color de cada pixel 
        fig = plt.figure(figsize=[30,10])
        fig.add_subplot(221)
        sc = plt.scatter(lon_transect[nlev,:], approx_altitude[nlev,:],
                s=1,c=ZDR_transect[nlev,:],
                cmap=cmap_ZDR, vmin=-2, vmax=5.)
        color[nlev,:,:] = sc.to_rgba(ZDR_transect[nlev,:])
        plt.close()
    
    #- Try polygons
    #fig1.add_subplot(412)
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
        # Create the cone for each elevation IN TERMS OF RANGE. 
        ancho_haz_i0    = (np.pi/180*gate_range[nlev,0]/2)
        ancho_haz_i1099 = (np.pi/180*gate_range[nlev,azydims]/2)
        P1 = Polygon([( gate_range[nlev,0],   approx_altitude[nlev,0]-ancho_haz_i0      ),
                  ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]-ancho_haz_i1099),
                  ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]+ancho_haz_i1099),
                  ( gate_range[nlev,0],       approx_altitude[nlev,0]+ancho_haz_i0      )])
        ancho = 50/1E3
        # Location of gates? Every 100m?  
        LS = [Polygon([(gate_range[nlev,x]-ancho, 0),
                  (gate_range[nlev,x]+ancho, 0),
                  (gate_range[nlev,x]+ancho, 50),
                  (gate_range[nlev,x]-ancho, 50)]) for x in np.arange(approx_altitude.shape[1])]
        # Plot
        #ax1 = plt.gca()
        for i, l in enumerate(LS):
            inter = l.intersection(P1)
            x,y = inter.exterior.xy
            # Then plot it, filled by the color we want
            axes[1].fill(x, y, color = color[nlev,i,:], )
            x, y = P1.exterior.xy
    axes[1].set_ylim([0, 15])
    axes[1].set_ylabel('Altitude (km)')
    axes[1].grid()
    axes[1].set_xlim((xlim_range1, xlim_range2))
    norm = matplotlib.colors.Normalize(vmin=-2.,vmax=5.)
    cax = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_ZDR)
    cax.set_array(ZDR_transect)
    cbar_zdr = fig2.colorbar(cax, ax=axes[1], shrink=1.1, ticks=np.arange(-2.,5.01,1.), label='ZDR')     
    axes[1].axhline(y=freezing_lev1,color='k',linestyle='--', linewidth=1.2)
    axes[1].axhline(y=freezing_lev2,color='k',linestyle='--', linewidth=1.2)
    
    del mycolorbar, x, y, inter
    #---------------------------------------- RHOHV
    #- Simple pcolormesh plot! 
    fig = plt.figure(figsize=[15,11])
    fig.add_subplot(221)
    mycolorbar = plt.pcolormesh(lon_transect, approx_altitude,
                RHO_transect,
                cmap = pyart.graph.cm.RefDiff, vmin=0.7, vmax=1.)
    plt.close()
    
    #- De esta manera me guardo el color con el que rellenar los polygons
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
        # scatter plot para sacar el color de cada pixel 
        fig = plt.figure(figsize=[30,10])
        fig.add_subplot(221)
        sc = plt.scatter(lon_transect[nlev,:], approx_altitude[nlev,:],
                s=1,c=RHO_transect[nlev,:],
                cmap= pyart.graph.cm.RefDiff, vmin=0.7, vmax=1.)
        color[nlev,:,:] = sc.to_rgba(RHO_transect[nlev,:])   # pyart.graph.cm.RefDiff
        plt.close()
    
    #- Try polygons
    #fig1.add_subplot(412)
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
        # Create the cone for each elevation IN TERMS OF RANGE. 
        ancho_haz_i0    = (np.pi/180*gate_range[nlev,0]/2)
        ancho_haz_i1099 = (np.pi/180*gate_range[nlev,azydims]/2)
        P1 = Polygon([( gate_range[nlev,0],   approx_altitude[nlev,0]-ancho_haz_i0      ),
                  ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]-ancho_haz_i1099),
                  ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]+ancho_haz_i1099),
                  ( gate_range[nlev,0],       approx_altitude[nlev,0]+ancho_haz_i0      )])
        ancho = 50/1E3
        # Location of gates? Every 100m?  
        LS = [Polygon([(gate_range[nlev,x]-ancho, 0),
                  (gate_range[nlev,x]+ancho, 0),
                  (gate_range[nlev,x]+ancho, 50),
                  (gate_range[nlev,x]-ancho, 50)]) for x in np.arange(approx_altitude.shape[1])]
        # Plot
        #ax1 = plt.gca()
        for i, l in enumerate(LS):
            inter = l.intersection(P1)
            x,y = inter.exterior.xy
            # Then plot it, filled by the color we want
            axes[2].fill(x, y, color = color[nlev,i,:], )
            x, y = P1.exterior.xy
    axes[2].set_ylim([0, 15])
    axes[2].set_ylabel('Altitude (km)')
    axes[2].set_xlabel('Range (km)')
    
    axes[2].grid()
    axes[2].set_xlim((xlim_range1, xlim_range2))
    norm = matplotlib.colors.Normalize(vmin=0.7,vmax=1.)
    cax = matplotlib.cm.ScalarMappable(norm=norm, cmap=pyart.graph.cm.RefDiff)
    cax.set_array(RHO_transect)
    cbar_rho = fig2.colorbar(cax, ax=axes[2], shrink=1.1, ticks=np.arange(0.7,1.01,0.1), label=r'$\rho_{hv}$')     
    axes[2].axhline(y=freezing_lev1,color='k',linestyle='--', linewidth=1.2)
    axes[2].axhline(y=freezing_lev2,color='k',linestyle='--', linewidth=1.2)
    
    del mycolorbar, x, y, inter
    
    #- savefile
    fig2.suptitle('RMA1 Transect '+ time, fontweight='bold')
    #fig.savefig(fig_dir+'pseudo_RHI'+str(file)+'.png', dpi=300,transparent=False)
    save_dir_compare = folders['save_dir_compare']
    fig2.savefig(save_dir_compare+'/OBS/RMA1'+'/ZH_RMA1_obs_pseudo_RHI'+time+'.png', dpi=300,transparent=False, bbox_inches='tight')
    
    plt.show() 
       
        
    
    
    return


#------------------------------------------------------------------------------  
#------------------------------------------------------------------------------  
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#------------------------------------------------------------------------------  
#------------------------------------------------------------------------------
def plot_psuedo_HID_CSPR2(test_transect, time, rfiles, ZDRoffset, xlim_range1, xlim_range2, freezing_lev):
    

    folders      = config_folders.config_folders('yakaira')
    radar_folder = folders['csapr2_dir']
    filename     = os.path.join(radar_folder, rfiles)
    radar        = pyart.io.read(filename) 

    #- Radar sweep
    nelev       = 0
    start_index  = radar.sweep_start_ray_index['data'][nelev]
    end_index    = radar.sweep_end_ray_index['data'][nelev]
    lats0        = radar.gate_latitude['data'][start_index:end_index]
    lons0        = radar.gate_longitude['data'][start_index:end_index]
    azimuths     = radar.azimuth['data'][start_index:end_index] 
    Ze_transect     = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); Ze_transect[:]=np.nan
    ZDR_transect    = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); ZDR_transect[:]=np.nan
    PHI_transect    = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); PHI_transect[:]=np.nan
    lon_transect    = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); lon_transect[:]=np.nan
    lat_transect    = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); lat_transect[:]=np.nan
    RHO_transect    = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); RHO_transect[:]=np.nan
    approx_altitude = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); approx_altitude[:]=np.nan
    color           = np.full((  len(radar.sweep_start_ray_index['data']), lats0.shape[1], 4), np.nan)
    gate_range      = np.zeros( [len(radar.sweep_start_ray_index['data']), lats0.shape[1] ]); gate_range[:]=np.nan

    azydims = lats0.shape[1]-1
    window = 15 
    
    [lats_, lons_, _] = radar.get_gate_lat_lon_alt(1, reset_gate_coords=False, filter_transitions=False)
    
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
        start_index = radar.sweep_start_ray_index['data'][nlev]
        end_index   = radar.sweep_end_ray_index['data'][nlev]       
        ZHZH = radar.fields['attenuation_corrected_reflectivity_h']['data'][start_index:end_index]	
        ZDRZDR      = radar.fields['attenuation_corrected_differential_reflectivity']['data'][start_index:end_index]
        RHORHO      = radar.fields['copol_correlation_coeff']['data'][start_index:end_index]       
        #ZDRZDR[RHORHO<0.75]=np.nan
        #RHORHO[RHORHO<0.75]=np.nan
        lats        = radar.gate_latitude['data'][start_index:end_index]
        lons        = radar.gate_longitude['data'][start_index:end_index]
        # En verdad buscar azimuth no transecta ... 
        azimuths    = radar.azimuth['data'][start_index:end_index]
        filas = find_nearest(azimuths, test_transect)
        lon_transect[nlev,:]     = lons[filas,:]
        lat_transect[nlev,:]     = lats[filas,:]
        #
        gateZ    = radar.gate_z['data'][start_index:end_index]
        gateX    = radar.gate_x['data'][start_index:end_index]
        gateY    = radar.gate_y['data'][start_index:end_index]
        gates_range  = np.sqrt(gateX**2 + gateY**2 + gateZ**2)
        #
        Ze_transect[nlev,:]      = pyart.correct.phase_proc.smooth_and_trim( ZHZH[filas,:], window)
        ZDR_transect[nlev,:]     = pyart.correct.phase_proc.smooth_and_trim(ZDRZDR[filas,:], window)
        RHO_transect[nlev,:]     = pyart.correct.phase_proc.smooth_and_trim( RHORHO[filas,:], window)
        # 
        [xgate, ygate, zgate]   = pyart.core.antenna_to_cartesian(gates_range[filas,:]/1e3, azimuths[filas],radar.get_elevation(nlev)[0]);
        approx_altitude[nlev,:] = zgate/1e3
        gate_range[nlev,:]      = gates_range[filas,:]/1e3;
	                
    #---------------------------------------- REFLECTIVITY
    #- Simple pcolormesh plot! 
    fig = plt.figure(figsize=[15,11])
    fig.add_subplot(111)
    mycolorbar = plt.pcolormesh(lon_transect, approx_altitude, Ze_transect, cmap=colormaps('ref'), vmin=0, vmax=60)

    #- De esta manera me guardo el color con el que rellenar los polygons (scatter plot para sacar el color de cada pixel)
    print(len(radar.sweep_start_ray_index['data']))
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
         fig = plt.figure(figsize=[30,10])
         fig.add_subplot(221)
         sc = plt.scatter(lon_transect[nlev,:], approx_altitude[nlev,:],
                 s=1,c=Ze_transect[nlev,:],
                 cmap=colormaps('ref'), vmin=0, vmax=60)
         color[nlev,:,:] = sc.to_rgba(Ze_transect[nlev,:])
         plt.close()

    #- Try polygons
    fig2, axes = plt.subplots(nrows=3,ncols=1,constrained_layout=True,figsize=[8,6])  # 8,4 muy chiquito
    fig1 = plt.figure(figsize=(15,20))
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
         if nlev > 9: continue
         # Create the cone for each elevation IN TERMS OF RANGE. 
         # ===> ACA HABRIA QUE AGREGAR COMO CAMBIA LA ALTURA CON EL RANGE (?)
         ancho_haz_i0    = (np.pi/180*gate_range[nlev,0]/2)
         ancho_haz_i1099 = (np.pi/180*gate_range[nlev,azydims]/2)
         P1 = Polygon([( gate_range[nlev,0],    approx_altitude[nlev,0]-ancho_haz_i0      ),
                   ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]-ancho_haz_i1099),
                   ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]+ancho_haz_i1099),
                   ( gate_range[nlev,0],    approx_altitude[nlev,0]+ancho_haz_i0      )])
         ancho = 50/1E3
         # Location of gates? Every 100m?  
         LS = [Polygon([(gate_range[nlev,x]-ancho, 0),
                   (gate_range[nlev,x]+ancho, 0),
                   (gate_range[nlev,x]+ancho, 50),
                   (gate_range[nlev,x]-ancho, 50)]) for x in np.arange(approx_altitude.shape[1])]
         # Plot
         for i, l in enumerate(LS):
             # Get the polygon of the intersection between the cone and the space 
             #reserved for a specific point
             inter = l.intersection(P1)
             x,y = inter.exterior.xy    
             # Then plot it, filled by the color we want
             axes[0].fill(x, y, color = color[nlev,i,:], )
             x, y = P1.exterior.xy
    axes[0].set_ylim([0, 20])
    axes[0].set_ylabel('Altitude (km)')
    axes[0].grid()
    axes[0].set_xlim((xlim_range1, xlim_range2))
    norm = matplotlib.colors.Normalize(vmin=0.,vmax=60.)
    cax = matplotlib.cm.ScalarMappable(norm=norm, cmap=colormaps('ref'))
    cax.set_array(Ze_transect)
    cbar_z = fig2.colorbar(cax, ax=axes[0], shrink=1.1, ticks=np.arange(0,60.01,10), label='Zh (dBZ)')
    axes[0].axhline(y=freezing_lev,color='k',linestyle='--', linewidth=1.2)
    del mycolorbar, x, y, inter
    #---------------------------------------- ZDR
    N = (5+2)
    cmap_ZDR = discrete_cmap(int(N), 'jet') 
    #- Simple pcolormesh plot! 
    fig = plt.figure(figsize=[15,11])
    fig.add_subplot(221)
    mycolorbar = plt.pcolormesh(lon_transect, approx_altitude,
                ZDR_transect,
                cmap=cmap_ZDR, vmin=-2, vmax=5.)
    plt.close()

    #- De esta manera me guardo el color con el que rellenar los polygons
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
        # scatter plot para sacar el color de cada pixel 
        fig = plt.figure(figsize=[30,10])
        fig.add_subplot(221)
        sc = plt.scatter(lon_transect[nlev,:], approx_altitude[nlev,:],
                s=1,c=ZDR_transect[nlev,:],
                cmap=cmap_ZDR, vmin=-2, vmax=5.)
        color[nlev,:,:] = sc.to_rgba(ZDR_transect[nlev,:])
        plt.close()

    #- Try polygons
    #fig1.add_subplot(412)
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
        if nlev > 9: continue
        # Create the cone for each elevation IN TERMS OF RANGE. 
        ancho_haz_i0    = (np.pi/180*gate_range[nlev,0]/2)
        ancho_haz_i1099 = (np.pi/180*gate_range[nlev,azydims]/2)
        P1 = Polygon([( gate_range[nlev,0],   approx_altitude[nlev,0]-ancho_haz_i0      ),
                  ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]-ancho_haz_i1099),
                  ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]+ancho_haz_i1099),
                  ( gate_range[nlev,0],       approx_altitude[nlev,0]+ancho_haz_i0      )])
        ancho = 50/1E3
        # Location of gates? Every 100m?  
        LS = [Polygon([(gate_range[nlev,x]-ancho, 0),
                  (gate_range[nlev,x]+ancho, 0),
                  (gate_range[nlev,x]+ancho, 50),
                  (gate_range[nlev,x]-ancho, 50)]) for x in np.arange(approx_altitude.shape[1])]
        # Plot
        #ax1 = plt.gca()
        for i, l in enumerate(LS):
            inter = l.intersection(P1)
            x,y = inter.exterior.xy
            # Then plot it, filled by the color we want
            axes[1].fill(x, y, color = color[nlev,i,:], )
            x, y = P1.exterior.xy
    axes[1].set_ylim([0, 20])
    axes[1].set_ylabel('Altitude (km)')
    axes[1].grid()
    axes[1].set_xlim((xlim_range1, xlim_range2))
    norm = matplotlib.colors.Normalize(vmin=-2.,vmax=5.)
    cax = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_ZDR)
    cax.set_array(ZDR_transect)
    cbar_zdr = fig2.colorbar(cax, ax=axes[1], shrink=1.1, ticks=np.arange(-2.,5.01,1.), label='ZDR')     
    axes[1].axhline(y=freezing_lev,color='k',linestyle='--', linewidth=1.2)

    del mycolorbar, x, y, inter
    #---------------------------------------- RHOHV
    #- Simple pcolormesh plot! 
    fig = plt.figure(figsize=[15,11])
    fig.add_subplot(221)
    mycolorbar = plt.pcolormesh(lon_transect, approx_altitude,
                RHO_transect,
                cmap = pyart.graph.cm.RefDiff, vmin=0.7, vmax=1.)
    plt.close()

    #- De esta manera me guardo el color con el que rellenar los polygons
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
        # scatter plot para sacar el color de cada pixel 
        fig = plt.figure(figsize=[30,10])
        fig.add_subplot(221)
        sc = plt.scatter(lon_transect[nlev,:], approx_altitude[nlev,:],
                s=1,c=RHO_transect[nlev,:],
                cmap= pyart.graph.cm.RefDiff, vmin=0.7, vmax=1.)
        color[nlev,:,:] = sc.to_rgba(RHO_transect[nlev,:])   # pyart.graph.cm.RefDiff
        plt.close()

    #- Try polygons
    #fig1.add_subplot(412)
    for nlev in range(len(radar.sweep_start_ray_index['data'])):
        if nlev > 9: continue
        # Create the cone for each elevation IN TERMS OF RANGE. 
        ancho_haz_i0    = (np.pi/180*gate_range[nlev,0]/2)
        ancho_haz_i1099 = (np.pi/180*gate_range[nlev,azydims]/2)
        P1 = Polygon([( gate_range[nlev,0],   approx_altitude[nlev,0]-ancho_haz_i0      ),
                  ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]-ancho_haz_i1099),
                  ( gate_range[nlev,azydims], approx_altitude[nlev,azydims]+ancho_haz_i1099),
                  ( gate_range[nlev,0],       approx_altitude[nlev,0]+ancho_haz_i0      )])
        ancho = 50/1E3
        # Location of gates? Every 100m?  
        LS = [Polygon([(gate_range[nlev,x]-ancho, 0),
                  (gate_range[nlev,x]+ancho, 0),
                  (gate_range[nlev,x]+ancho, 50),
                  (gate_range[nlev,x]-ancho, 50)]) for x in np.arange(approx_altitude.shape[1])]
        # Plot
        #ax1 = plt.gca()
        for i, l in enumerate(LS):
            inter = l.intersection(P1)
            x,y = inter.exterior.xy
            # Then plot it, filled by the color we want
            axes[2].fill(x, y, color = color[nlev,i,:], )
            x, y = P1.exterior.xy
    axes[2].set_ylim([0, 20])
    axes[2].set_ylabel('Altitude (km)')
    axes[2].grid()
    axes[2].set_xlim((xlim_range1, xlim_range2))
    norm = matplotlib.colors.Normalize(vmin=0.7,vmax=1.)
    cax = matplotlib.cm.ScalarMappable(norm=norm, cmap=pyart.graph.cm.RefDiff)
    cax.set_array(RHO_transect)
    cbar_rho = fig2.colorbar(cax, ax=axes[2], shrink=1.1, ticks=np.arange(0.7,1.01,0.1), label=r'$\rho_{hv}$')     
    axes[2].axhline(y=freezing_lev,color='k',linestyle='--', linewidth=1.2)

    del mycolorbar, x, y, inter


    #- savefile
    fig2.suptitle('CSAPR2 Transect '+ time, fontweight='bold')
    #fig.savefig(fig_dir+'pseudo_RHI'+str(file)+'.png', dpi=300,transparent=False)
    save_dir_compare = folders['save_dir_compare']
    fig2.savefig(save_dir_compare+'/OBS/CSAPR2'+'/ZH_CSAPR2_obs_pseudo_RHI'+time+'.png', dpi=300,transparent=False, bbox_inches='tight')
    
    plt.show() 
       
    return




#------------------------------------------------------------------------------
def plot_precip_WRF(EXP, title, titleprevious, domain, servidor):

    
    folders   = config_folders.config_folders(servidor)
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
    

    prefix   = 'wrfout_'+domain+'_2018-11-10_'+title
    filename = os.path.join(WRFfolder, 'wrfout_'+domain+'_2018-11-10_'+title+':00')

    
    ncfile       = Dataset(filename,'r')        
    z            = wrf.getvar( ncfile,"z") 
    zh           = wrf.getvar(ncfile, "REFL_10CM")
    REFL_10CM    = wrf.interplevel(zh, z, 1000)
    lat          = wrf.getvar( ncfile,"lat") 
    lon          = wrf.getvar( ncfile,"lon")
    
    RAINNC       = wrf.getvar( ncfile,"RAINNC")
    
    # to get the previous hour and get the accumulated period of time (in this case 30mins)
    filename_preivous = os.path.join(WRFfolder, 'wrfout_'+domain+'_2018-11-10_'+titleprevious+':00')
    ncfile_preivous   = Dataset(filename_preivous,'r')        
    RAINNC_previous   = wrf.getvar( ncfile_preivous,"RAINNC")
    
    # 30MINS ACCUMULATED:
    precip = RAINNC - RAINNC_previous
   

    fig, ax = plt.subplots(figsize=(8,8)) 
    pcm = ax.pcolormesh(lon, lat,  precip, cmap=get_cmap("rainbow"), vmin=0,  vmax=50)   #MPL_gist_ncar
    cbar = plt.colorbar(pcm, ax=ax, shrink=1)
    ax.plot(prov[:,0],prov[:,1],color='k'); 
        
    if 'maite' in title: 
        ax.set_xlim([-66,-62]); 
        ax.set_ylim([-35,-31])
                
    else:
        ax.set_xlim([-65.5,-62]); 
        ax.set_ylim([-35,-31])
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
    # agrego contorno de 500 y 1000m
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['k','k'], linewidths=2)
    ax.contour(lon, lat, REFL_10CM, levels=[45], colors=['r'], linewidths=2)
        

    cbar.cmap.set_under('white')
    ax.grid()
    ax.set_title('Precip accumulated in 30 mins ('+title+')')
    
    ax.text(x=-65, y=-34.9, s='120 and 220 km radar rings')
    #plt.show()
    fig.savefig(save_dir_compare+'/'+EXP+'/WRF_PRECIPACCUM_30mins_'+domain+'_general_evolution_'+title+'.png', dpi=300,transparent=False,bbox_inches='tight')
    plt.close()
    
    return#------------------------------------------------------------------------------
