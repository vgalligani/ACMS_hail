import pyart 
from Plots4Analysis  import pyplot_rings
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import package_functions as funs
import Plots4Analysis as P4A
import config_folders
import netCDF4 as nc
import pandas as pd
import scipy.constants as spc
from csu_radartools import csu_fhc
import wradlib as wrl
import matplotlib.colors as colors


#------------------------------------------------------------------------------
def plot_Zh1km_radar_wTransect(rfiles, time, test_transect):

    prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')

    # RMA1 
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    TH_name       = 'TH'
    time1infile  = 82
     
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
    ZHelev18     = radar.get_field(3, TH_name, copy=False)
    [lats_, lons_, _] = radar.get_gate_lat_lon_alt(3, reset_gate_coords=False, filter_transitions=False)
    
    
    # --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----- ---- ---- ----        
    [temperature, pressure] = read_rsnd('/87344_20181110T170540Z.BUFR')  
    height                  = pressure2height(np.array(pressure), np.array(temperature) )
    freezing_lev            = calc_freezinglevel(temperature, pressure)    
    radar_T,radar_z         =  funs.interpolate_sounding_to_radar(np.array(temperature)-273, height, radar)
    radar = funs.add_field_to_radar_object(radar_T, radar, field_name='sounding_temperature')
    
    # --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----- ---- ---- ----        
    # CORRECION DE PHIDP#
    # Elimino cualquier mascara que venga en el dato
    PHIORIG = radar.fields['PHIDP']['data'].copy()
    mask = radar.fields['PHIDP']['data'].data.copy()
    mask[:] = False
    PHIORIG.mask = mask
    radar.add_field_like('PHIDP', 'PHIDP', PHIORIG, replace_existing=True)

    # Paso 1: Encontrar el sys_phase del radar
    sys_phase = funs.get_sys_phase_simple(radar)

    # Paso 2: correct phidp
    dphi, uphi, corr_phidp = funs.correct_phidp(radar, sys_phase, 280)

    # Paso 3: Calcular KDP. Probar el efecto de otros winlen
    calculated_KDP = wrl.dp.kdp_from_phidp(corr_phidp, winlen=7, dr=(radar.range['data'][1]-radar.range['data'][0])/1e3,
					   method='lanczos_conv', skipna=True)
    radar.add_field_like('RHOHV','corrKDP', calculated_KDP, replace_existing=True)
    radar.add_field_like('PHIDP','corrPHIDP', corr_phidp, replace_existing=True)

    # CORRECION ZH#
    a=0.08
    b=0.64884
    Atenua_espec, Z_correc_Gu = pyart.correct.attenuation.calculate_attenuation(radar, 0.0, refl_field='TH',
                                                                              ncp_field='RHOHV',rhv_field='RHOHV',
                                                                              phidp_field='corrPHIDP',
                                                                              a_coef=a,beta=b,rhv_min=0.5,debug=True)
    #rhv_min=0.85
    radar.add_field('dBZ_correc_Gu', Z_correc_Gu, replace_existing=True)
         
    do_plot = 0       
    if do_plot == 1:
        nelev = 0
        start_index = radar.sweep_start_ray_index['data'][nelev]
        end_index   = radar.sweep_end_ray_index['data'][nelev]
        lats0        = radar.gate_latitude['data'][start_index:end_index]
        lons0        = radar.gate_longitude['data'][start_index:end_index]
        field =  radar.fields['RHOHV']['data'][start_index:end_index]
        fig, ax = plt.subplots(figsize=[8,8])
        pcm = ax.pcolormesh(lons0, lats0, field, cmap=P4A.colormaps('ref'))
        cbar = plt.colorbar(pcm, ax=ax, shrink=1, label='Zh RMA1 elev 3')
        ax.grid()       
        ax.plot(prov[:,0],prov[:,1],color='k'); 
        ax.set_ylim( [-33.5,-30]  )
        ax.set_xlim( [-67,-61]  )
        
    # BEFORE HIDS FILTER TV FOR RHO<0.5 AND RHO<0.5 TO NANs
    rho_data = radar.fields['RHOHV']['data']
    zv_data  = radar.fields['TV']['data']
    corrKDP_data  = radar.fields['corrKDP']['data']
    th_data = radar.fields['dBZ_correc_Gu']['data']

    rho = rho_data.copy()

    zv_data[rho<0.7] = np.nan
    zv_data[rho<0.7] = np.nan
    corrKDP_data[rho<0.7] = np.nan
    th_data[rho<0.7] = np.nan
    rho_data[rho<0.7] = np.nan
    
    
    # HID CALC          
    fh          = csu_fhc.csu_fhc_summer(dz=th_data, zdr=radar.fields['dBZ_correc_Gu']['data']-zv_data,
                                         rho=rho_data, kdp=corrKDP_data, 
                                         use_temp=True, band='C', T=radar_T)
    
    radar = funs.add_field_to_radar_object(fh, radar)
    
    fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
                
    # que pasa aca??? con el cambio en calc_kdp_ray_fir.pyx ?! fix
    # https://github.com/CSU-Radarmet/CSU_RadarTools/issues/67
    # check el ejemplo? revisar si estan bien los fields del radar 1ue estoty ingresando
    
    #fig, ax = plt.subplots(figsize=[8,8])
    #--------------------------------------------
    nelev = 0
    savetitle = 'AllPolarimetrics_HIDS_time'+time+'_nlev'+str(nelev)
    fig, ax = plt.subplots(nrows=2, ncols=3, constrained_layout=True,figsize=[16,16])

    ylims_ranges = [-32.5,-31.8]
    xlims_ranges = [-64.5,-63.5]
    
    start_index = radar.sweep_start_ray_index['data'][nelev]
    end_index   = radar.sweep_end_ray_index['data'][nelev]
    lats0        = radar.gate_latitude['data'][start_index:end_index]
    lons0        = radar.gate_longitude['data'][start_index:end_index]
    
    # Zh
    ZHelev18 = th_data[start_index:end_index]
    pcm = ax[0,0].pcolormesh(lons0, lats0, ZHelev18, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
    cbar = plt.colorbar(pcm, ax=ax[0,0], shrink=1)
    ax[0,0].set_title('RMA1 Zh Gu corrected')
    ax[0,0].plot(prov[:,0],prov[:,1],color='k'); 

    # hid
    hid_colors = ['MediumBlue', 'DarkOrange', 'LightPink',
      'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    hidhid = fh.copy()
    hidhid = np.ma.array(hidhid)
    hidhid = np.ma.masked_where(rho < 0.7, hidhid)
    pcm = ax[0,1].pcolormesh(lons0, lats0, hidhid[start_index:end_index], cmap=cmaphid, vmin=1.8, vmax=10.4)
    cbar = plt.colorbar(pcm, ax=ax[0,1], shrink=1, label='CSU HID')
    cbar = funs.adjust_fhc_colorbar_for_pyart(cbar)
    ax[0,1].set_title('RMA1 HID')
    ax[0,1].plot(prov[:,0],prov[:,1],color='k'); 

    # RHO
    pcm = ax[1,0].pcolormesh(lons0, lats0, rho_data[start_index:end_index], cmap=P4A.colormaps('ref'))
    cbar = plt.colorbar(pcm, ax=ax[1,0], shrink=1)
    ax[1,0].set_title('RMA1 RHO')
    ax[1,0].plot(prov[:,0],prov[:,1],color='k'); 

    # KDP
    pcm = ax[1,1].pcolormesh(lons0, lats0, radar.fields['corrKDP']['data'][start_index:end_index],
                             cmap=P4A.colormaps('ref'), vmin=0, vmax=4)
    cbar = plt.colorbar(pcm, ax=ax[1,1], shrink=1)
    cbar.cmap.set_under('white')
    ax[1,1].set_title('RMA1 KDP')
    ax[1,1].plot(prov[:,0],prov[:,1],color='k'); 

    # ZDR
    zdr=(radar.fields['dBZ_correc_Gu']['data']-zv_data)
    pcm = ax[1,2].pcolormesh(lons0, lats0, zdr[start_index:end_index], cmap=P4A.colormaps('ref'), vmin=-5, vmax=10)
    cbar = plt.colorbar(pcm, ax=ax[1,2], shrink=1)
    ax[1,2].set_title('RMA1 ZDR')

    # TV
    tv = zv_data
    pcm = ax[0,2].pcolormesh(lons0, lats0, tv[start_index:end_index], cmap=P4A.colormaps('ref'))
    cbar = plt.colorbar(pcm, ax=ax[0,2], shrink=1)
    ax[0,2].set_title('RMA1 Zv')

    
    # En verdad buscar azimuth no transecta ... 
    azimuths       = radar.azimuth['data'][start_index:end_index]
    target_azimuth = azimuths[test_transect]  #- target azimuth for nlev=0 test case is 301.5
    filas          = np.asarray(abs(azimuths-target_azimuth)<=0.1).nonzero()
    lon_transect = lons_[filas,:]
    lat_transect = lats_[filas,:]
            
    samhi = pd.read_csv('/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/datos_base.csv')
    # Convert timestamp and extract only HH:MM
    samhi['formatted_time'] = pd.to_datetime(samhi['eventTime']).dt.strftime('%H:%M')

    # plot por severdad de granizo
    granizos_unk =samhi[samhi['hailMaxDiameterMM'].isna()]
    granizo_chico = samhi[(samhi['hailMaxDiameterMM']<20)] 
    granizo_chico_unk = pd.concat((granizos_unk, granizo_chico))
    granizo_sev = samhi[(samhi['hailMaxDiameterMM']>=20) & (samhi['hailMaxDiameterMM']<50)]
    granizo_grande = samhi[samhi['hailMaxDiameterMM']>=50]


    for ii in range(2):    
        ax[ii,2].plot(prov[:,0],prov[:,1],color='k'); 
        ax[ii,2].set_ylim( ylims_ranges )
        ax[ii,2].set_xlim( xlims_ranges )
        ax[ii,2].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
        ax[ii,2].plot(np.ravel(lon_transect), np.ravel(lat_transect), '-r', linewidth=1.5)
        ax[ii,2].scatter(granizo_chico_unk['longitude'].astype(float),granizo_chico_unk['latitude'].astype(float), s=200, marker = '^', color='grey', label='hail')
        ax[ii,2].scatter(granizo_sev['longitude'].astype(float),granizo_sev['latitude'].astype(float),             s=200, marker = '^', color='darkblue',  label='hail > 2 cm' )
        ax[ii,2].scatter(granizo_grande['longitude'].astype(float),granizo_grande['latitude'].astype(float),       s=200, marker = '^', color='darkred',  label='hail > 5 cm')
        ax[ii,2].scatter(-64.2177780,-31.9793750, s=200, marker = 'x', color='black', label='COW')
            
    for ii in range(2):    
        for jj in range(2):
            ax[ii,jj].set_ylim( ylims_ranges )
            ax[ii,jj].set_xlim( xlims_ranges )
            ax[ii,jj].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
            ax[ii,jj].plot(np.ravel(lon_transect), np.ravel(lat_transect), '-r', linewidth=1.5)
            ax[ii,jj].scatter(granizo_chico_unk['longitude'].astype(float),granizo_chico_unk['latitude'].astype(float), s=200, marker = '^', color='grey', label='hail')
            ax[ii,jj].scatter(granizo_sev['longitude'].astype(float),granizo_sev['latitude'].astype(float),             s=200, marker = '^', color='darkblue',  label='hail > 2 cm' )
            ax[ii,jj].scatter(granizo_grande['longitude'].astype(float),granizo_grande['latitude'].astype(float),       s=200, marker = '^', color='darkred',  label='hail > 5 cm')
            ax[ii,jj].scatter(-64.2177780,-31.9793750, s=200, marker = 'x', color='black', label='COW')
            
    leg = ax[1,1].legend(loc='lower left',  fontsize=12, framealpha=0.80)
    plt.show()
 
    save_dir_compare = folders['save_dir_compare']
    fig.savefig(save_dir_compare+'/OBS'+'/RMA1/'+'radarOBS_1km_wtransect_'+savetitle+'.png', dpi=300,transparent=False)
    

    return

#------------------------------------------------------------------------------
def plot_radar_cspr2_wcontour(filename, time, test_transect):
    
    prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')

    folders          = config_folders.config_folders('yakaira')
    radar_folder     = folders['csapr2_dir']
    save_dir_compare = folders['save_dir_compare']
       
    latrange = [-32.5,-31.8]  #[-33,-31.5] 
    lonrange = [-64.5,-63.5] #[-65.5,-63.5]

    #latrange =[-33,-31.5] 
    #lonrange = [-65.5,-63.5]
    
    radarLAT = -32.12641
    radarLON = -64.72837
    TH_name  = 'attenuation_corrected_reflectivity_h'
    
    [lat_radius, lon_radius] = pyplot_rings(radarLAT,radarLON,30)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT,radarLON,60) 
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT,radarLON,100) 
    
    #prefix = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/CSAPR2/corcsapr2cfrppiqcM1.b1.20181110'

    radar    = pyart.io.read(radar_folder+filename) 
    ZHelev18 = radar.get_field(0, TH_name, copy=False)
    [lats_, lons_, _] = radar.get_gate_lat_lon_alt(0, reset_gate_coords=False, filter_transitions=False)

    # --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----- ---- ---- ----        
    [temperature, pressure] = read_rsnd('/87344_20181110T170540Z.BUFR')  
    height                  = pressure2height(np.array(pressure), np.array(temperature) )
    freezing_lev            = calc_freezinglevel(temperature, pressure)    
    radar_T,radar_z         =  funs.interpolate_sounding_to_radar(np.array(temperature)-273, height, radar)
    # --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----- ---- ---- ----    
    
    fh    = csu_fhc.csu_fhc_summer(dz=radar.fields['attenuation_corrected_reflectivity_h']['data'], 
                                   zdr=radar.fields['attenuation_corrected_differential_reflectivity']['data'],
                                         rho=radar.fields['copol_correlation_coeff']['data'], 
                                         kdp=radar.fields['specific_differential_phase']['data'], 
                                         use_temp=True, band='C', T=radar_T)
    radar = funs.add_field_to_radar_object(fh, radar)
    
    fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
                
    # que pasa aca??? con el cambio en calc_kdp_ray_fir.pyx ?! fix
    # https://github.com/CSU-Radarmet/CSU_RadarTools/issues/67
    # check el ejemplo? revisar si estan bien los fields del radar 1ue estoty ingresando
    
    #fig, ax = plt.subplots(figsize=[8,8])
    #--------------------------------------------
    nelev = 0
    savetitle = 'AllPolarimetrics_HIDS_time'+time+'_nlev'+str(nelev)

    fig, ax = plt.subplots(nrows=2, ncols=3, constrained_layout=True,figsize=[16,16])

    ylims_ranges = [-32.5,-31.8]
    xlims_ranges = [-64.5,-63.5]
    
    start_index = radar.sweep_start_ray_index['data'][nelev]
    end_index   = radar.sweep_end_ray_index['data'][nelev]
    lats0        = radar.gate_latitude['data'][start_index:end_index]
    lons0        = radar.gate_longitude['data'][start_index:end_index]
    
    # Zh
    zhfield = 'attenuation_corrected_reflectivity_h'
    ZHelev18 = radar.fields[zhfield]['data'][start_index:end_index]
    pcm = ax[0,0].pcolormesh(lons0, lats0, ZHelev18, cmap=P4A.colormaps('ref'), vmin=0,  vmax=70)
    cbar = plt.colorbar(pcm, ax=ax[0,0], shrink=1)
    ax[0,0].set_title('CSPR2 Zh Gu corrected')
    ax[0,0].plot(prov[:,0],prov[:,1],color='k'); 

    # hid
    hid_colors = ['MediumBlue', 'DarkOrange', 'LightPink',
      'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    hidhid = fh.copy()
    hidhid = np.ma.array(hidhid)
    hidhid = np.ma.masked_where(radar.fields['copol_correlation_coeff']['data'] < 0.7, hidhid)
    pcm = ax[0,1].pcolormesh(lons0, lats0, hidhid[start_index:end_index], cmap=cmaphid, vmin=1.8, vmax=10.4)
    cbar = plt.colorbar(pcm, ax=ax[0,1], shrink=1, label='CSU HID')
    cbar = funs.adjust_fhc_colorbar_for_pyart(cbar)
    ax[0,1].set_title('CSPR2 HID')
    ax[0,1].plot(prov[:,0],prov[:,1],color='k'); 

    # RHO
    rho_data = radar.fields['copol_correlation_coeff']['data']
    pcm = ax[1,0].pcolormesh(lons0, lats0, rho_data[start_index:end_index], cmap=P4A.colormaps('ref'))
    cbar = plt.colorbar(pcm, ax=ax[1,0], shrink=1)
    ax[1,0].set_title('CSPR2 RHO')
    ax[1,0].plot(prov[:,0],prov[:,1],color='k'); 

    # KDP
    pcm = ax[1,1].pcolormesh(lons0, lats0, radar.fields['specific_differential_phase']['data'][start_index:end_index],
                             cmap=P4A.colormaps('ref'), vmin=0, vmax=4)
    cbar = plt.colorbar(pcm, ax=ax[1,1], shrink=1)
    cbar.cmap.set_under('white')
    ax[1,1].set_title('CSPR2 KDP')
    ax[1,1].plot(prov[:,0],prov[:,1],color='k'); 

    # ZDR
    zdr=radar.fields['attenuation_corrected_differential_reflectivity']['data']
    pcm = ax[1,2].pcolormesh(lons0, lats0, zdr[start_index:end_index], cmap=P4A.colormaps('ref'), vmin=-5, vmax=10)
    cbar = plt.colorbar(pcm, ax=ax[1,2], shrink=1)
    ax[1,2].set_title('CSPR2 ZDR')

    # TV
    tv = radar.fields['uncorrected_reflectivity_v']['data']
    pcm = ax[0,2].pcolormesh(lons0, lats0, tv[start_index:end_index], cmap=P4A.colormaps('ref'))
    cbar = plt.colorbar(pcm, ax=ax[0,2], shrink=1)
    ax[0,2].set_title('CSPR2 uncorrected Zv')
    
    # En verdad buscar azimuth no transecta ... 
    azimuths       = radar.azimuth['data'][start_index:end_index]
    filas          = P4A.find_nearest(azimuths, test_transect)
    lon_transect = lons_[filas,:]
    lat_transect = lats_[filas,:]
    
            
    samhi = pd.read_csv('/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/datos_base.csv')
    # Convert timestamp and extract only HH:MM
    samhi['formatted_time'] = pd.to_datetime(samhi['eventTime']).dt.strftime('%H:%M')

    # plot por severdad de granizo
    granizos_unk =samhi[samhi['hailMaxDiameterMM'].isna()]
    granizo_chico = samhi[(samhi['hailMaxDiameterMM']<20)] 
    granizo_chico_unk = pd.concat((granizos_unk, granizo_chico))
    granizo_sev = samhi[(samhi['hailMaxDiameterMM']>=20) & (samhi['hailMaxDiameterMM']<50)]
    granizo_grande = samhi[samhi['hailMaxDiameterMM']>=50]


    for ii in range(2):    
        ax[ii,2].plot(prov[:,0],prov[:,1],color='k'); 
        ax[ii,2].set_ylim( ylims_ranges )
        ax[ii,2].set_xlim( xlims_ranges )
        ax[ii,2].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
        ax[ii,2].plot(np.ravel(lon_transect), np.ravel(lat_transect), '-r', linewidth=1.5)
        ax[ii,2].scatter(granizo_chico_unk['longitude'].astype(float),granizo_chico_unk['latitude'].astype(float), s=200, marker = '^', color='grey', label='hail')
        ax[ii,2].scatter(granizo_sev['longitude'].astype(float),granizo_sev['latitude'].astype(float),             s=200, marker = '^', color='darkblue',  label='hail > 2 cm' )
        ax[ii,2].scatter(granizo_grande['longitude'].astype(float),granizo_grande['latitude'].astype(float),       s=200, marker = '^', color='darkred',  label='hail > 5 cm')
        ax[ii,2].scatter(-64.2177780,-31.9793750, s=200, marker = 'x', color='black', label='COW')
            
    for ii in range(2):    
        for jj in range(2):
            ax[ii,jj].set_ylim( ylims_ranges )
            ax[ii,jj].set_xlim( xlims_ranges )
            ax[ii,jj].contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
            ax[ii,jj].plot(np.ravel(lon_transect), np.ravel(lat_transect), '-r', linewidth=1.5)
            ax[ii,jj].scatter(granizo_chico_unk['longitude'].astype(float),granizo_chico_unk['latitude'].astype(float), s=200, marker = '^', color='grey', label='hail')
            ax[ii,jj].scatter(granizo_sev['longitude'].astype(float),granizo_sev['latitude'].astype(float),             s=200, marker = '^', color='darkblue',  label='hail > 2 cm' )
            ax[ii,jj].scatter(granizo_grande['longitude'].astype(float),granizo_grande['latitude'].astype(float),       s=200, marker = '^', color='darkred',  label='hail > 5 cm')
            ax[ii,jj].scatter(-64.2177780,-31.9793750, s=200, marker = 'x', color='black', label='COW')
            
    leg = ax[1,1].legend(loc='lower left',  fontsize=12, framealpha=0.80)
    plt.show()
 
    save_dir_compare = folders['save_dir_compare']
    fig.savefig(save_dir_compare+'/OBS'+'/CSAPR2/'+'radarOBS_1km_wtransect_'+savetitle+'.png', dpi=300,transparent=False)
    
    return



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def pressure2height(p, T):
    r"""Convert pressure to height based on the hydrostatic equilibrium.

    .. math::
       z = \int -\frac{\mathrm{d}p}{\rho g}

    Parameters:
        p (ndarray): Pressure [Pa].
        T (ndarray): Temperature [K].
            If ``None`` the standard atmosphere is assumed.

    See also:
        .. autosummary::
            :nosignatures:

            standard_atmosphere

    Returns:
        ndarray: Relative height above lowest pressure level [m].
    """

    layer_depth = np.diff(p)
    R   =  spc.gas_constant / 28.9645e-3  # J K^-1 kg^-1
    rho = p / (R * T)   
    rho_layer = 0.5 * (rho[:-1] + rho[1:])

    z = np.cumsum(-layer_depth / (rho_layer * spc.g ))

    return np.hstack([0, z])
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def calc_freezinglevel(temperature, pressure):

    zfield     = pressure2height(np.array(pressure), np.array(temperature) )

    tfield_ref = np.array(temperature) - 273 # convert to C
    
    freezing_lev_pres = np.array(pressure[P4A.find_nearest(tfield_ref, 0)]) # Pa
    freezing_lev      = zfield[P4A.find_nearest(tfield_ref, 0)] 
    
    temp_             = np.array(tfield_ref[P4A.find_nearest(tfield_ref, 0)])

    return freezing_lev


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def read_rsnd(file):
    
    # RE-CALCULATE FREEZING LEVEL AND DRIVE PDATE. 
    sndfolder = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/Radiosonde_Ar/' 
    sndfile   = sndfolder+file

    import eccodes 
    #import codes_bufr_new_from_file, codes_release, codes_get, codes_set, codes_unpack_bufr

    # Open the BUFR file
    file_path = sndfile
    
    pres = []
    temp = []
    
    with open(file_path, "rb") as f:
        while True:
            # Read BUFR message
            bufr = eccodes.codes_bufr_new_from_file(f)
            if bufr is None:
                break  # End of file
            
            # Decode the message
            eccodes.codes_set(bufr, "unpack", 1)
            
            try:
            
                # Get the units of pressure
                pressure_units = eccodes.codes_get(bufr, "pressure->units")
                print(f"Pressure units: {pressure_units}")
    
                # Extract radiosonde metadata
                station_id = eccodes.codes_get_array(bufr, "blockNumber")[0]  # Station ID
                latitude   = eccodes.codes_get_array(bufr, "latitude")[0]
                longitude  = eccodes.codes_get_array(bufr, "longitude")[0]                
                
                # Extract the full profile using sequencing
                num_subsets = eccodes.codes_get_array(bufr, "numberOfSubsets")[0]  # Number of levels
                pressure    = eccodes.codes_get_array(bufr, "pressure")            # Pressure levels
                temperature = eccodes.codes_get_array(bufr, "airTemperature")      # Temperature levels
    
                # Ensure we only display levels where both pressure & temperature exist
                min_length = min(len(pressure), len(temperature))  # Avoid mismatch issues
    
                for i in range(min_length):
                    pres.append(pressure[i])
                    temp.append(temperature[i])
                    
            except Exception as e:
                print(f"Error reading BUFR: {e}")
        
            # Release BUFR message
            eccodes.codes_release(bufr)

    return temp, pres


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main_rma1(rfile, time, test_transect):
    
    [temperature, pressure] = read_rsnd('/87344_20181110T170540Z.BUFR')    
    freezing_lev1 = calc_freezinglevel(temperature, pressure)
    del temperature, pressure 

    [temperature, pressure] = read_rsnd('/87344_20181110T200424Z.BUFR')    
    freezing_lev2 = calc_freezinglevel(temperature, pressure)
    del temperature, pressure 
    
    xlim_range1 = 0  
    xlim_range2 = 100
    ZDRoffset = 0
    plot_Zh1km_radar_wTransect(rfile, time, test_transect)
    #P4A.plot_psuedo_HID_RMA1(test_transect, time, rfile, ZDRoffset, 
    #                         xlim_range1, xlim_range2, freezing_lev1/1000 , freezing_lev2/1000)
        
    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main_cspr2(rfile, time, test_transect):
    
    [temperature, pressure] = read_rsnd('/87344_20181110T170540Z.BUFR')    
    freezing_lev1 = calc_freezinglevel(temperature, pressure)
    del temperature, pressure 

    [temperature, pressure] = read_rsnd('/87344_20181110T200424Z.BUFR')    
    freezing_lev2 = calc_freezinglevel(temperature, pressure)
    del temperature, pressure 
    
    xlim_range1 = 0  
    xlim_range2 = 110
    ZDRoffset = 0
    plot_radar_cspr2_wcontour(rfile, time, test_transect)
    #P4A.plot_psuedo_HID_CSPR2(test_transect, time, rfile, ZDRoffset, 
    #                         xlim_range1, xlim_range2, freezing_lev)
             
    return


# ============================================================================= RMA1
def run_rma1():
    
    # time = '20:00'
    rfiles = 'cfrad.20181110_200709.0000_to_20181110_201358.0000_RMA1_0301_01.nc'
    main_rma1(rfiles, '20:00', 180)
         
    # time = 20:15
    rfiles = 'cfrad.20181110_201524.0000_to_20181110_202213.0000_RMA1_0301_01.nc'
    main_rma1(rfiles, '20:15', 173)
    
    return



# ============================================================================= CSAPR2
def run_cspr2():

    filename = 'corcsapr2cfrppiqcM1.b1.20181110.200003.nc'
    time = '20:00'
    main_cspr2(filename, time, 73)
    
    filename = 'corcsapr2cfrppiqcM1.b1.20181110.201746.nc'
    time = '20:17'
    main_cspr2(filename, time, 85)

    return



# ================================= RUN 
run_cspr2()






