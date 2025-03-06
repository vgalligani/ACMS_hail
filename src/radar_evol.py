import pyart 
import config_folders
import os
import glob 
import numpy as np
import package_functions as funs
import matplotlib.pyplot as plt
import netCDF4 as nc
from Plots4Analysis  import pyplot_rings
import pandas as pd


def do_plot():
    
    
    # RMA1 
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    TH_name       = 'TH'
    time1infile  = 82
    
    radarLAT_RMA1 = -31.441389
    radarLON_RMA1 = -64.191944
    [lat_radius, lon_radius] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,120)   
    [lat_radius2, lon_radius2] = pyplot_rings(radarLAT_RMA1,radarLON_RMA1,220)  
    
    
    folders=config_folders.config_folders('yakaira')
    save_dir_compare = folders['save_dir_compare']

    radar_folder = folders['rma1_dir']
    wildfiles    = 'cfrad.20181110*01.nc'  # try also with .02.nc maybe?     
    file_pattern = os.path.join(radar_folder, wildfiles)
    file_list    = sorted(glob.glob(file_pattern))
    timesr       = [filename[time1infile:time1infile+4] for filename in file_list]    
    
    rounded_minutes = [int(round(funs.hhmm_to_minutes(t) / 10) * 10) for t in timesr]

    # Create an array with NaN values for all intervals
    timereso   = 30
    start_time = 17*60 # remplazar por 12
    end_time   = 21*60 #(21*60)+30
        
    all_times_minutes = np.arange(start_time, end_time + timereso, timereso) 
    time_array        = np.full(len(all_times_minutes), np.nan, dtype=object)
    filename_array    = np.full(len(all_times_minutes), np.nan, dtype=object)
    
    # Fill in times (and filenames) where data is available
    for i, rounded_time in enumerate(rounded_minutes):
        index = np.where(all_times_minutes == rounded_time)[0]
        if index.size > 0:
            time_array[index[0]] = funs.minutes_to_hhmm(rounded_time)
            filename_array[index[0]] = file_list[i]

    # Convert all_times_minutes back to HHMM format for easy reading
    all_times_hhmm = [funs.minutes_to_hhmm(m) for m in all_times_minutes]    
                
    ZH_array = []

    for filename in filename_array:
        
        times = filename[time1infile:time1infile+4]
        radar       = pyart.io.read(filename) 
        
        # Configure a gatefilter to filter out copolar correlation coefficient values > 0.9
        gatefilter = pyart.filters.GateFilter(radar)
        gatefilter.exclude_transition()
        gatefilter.exclude_equal('RHOHV', 0.9)
        
        #ZHelev18 = radar.get_field(3, TH_name, copy=False)
        lon_0,lat_0,out_compz = funs.get_colmax(radar, TH_name, gatefilter)


        ZH_array.append(np.array(out_compz))
        del out_compz
        #ZH_array.append(ZHelev18.data)


    # Thresholding (values > 45)
    threshold = 45
    t_idx, lat_idx, lon_idx = np.where(np.array(ZH_array) > threshold)  # Get indices where data > 45

    # Extract corresponding lat, lon, and time values
    [lats_, lons_, _] = radar.get_gate_lat_lon_alt(0, reset_gate_coords=False, filter_transitions=False)
    # Expand lat/lon to match (time, lat, lon) shape
    lons_3d = np.repeat(lons_[np.newaxis, :, :], len(filename_array), axis=0)  # Shape (8, 350, 400)
    lats_3d = np.repeat(lats_[np.newaxis, :, :], len(filename_array), axis=0)  # Shape (8, 350, 400)
    # Map
    contour_lons = lons_3d[t_idx, lat_idx, lon_idx]  # Correct lon values
    contour_lats = lats_3d[t_idx, lat_idx, lon_idx]  # Correct lat values
    contour_times = t_idx  # Assign correct time step

    #Other info for figure
    prov = np.genfromtxt("/home/vito.galligani/Work/Tools/Maps/provincias.txt", delimiter='')
    fn = '/home/vito.galligani/Work/Tools/etopo1_bedrock.nc'
    ds = nc.Dataset(fn)
    topo_lat = ds.variables['lat'][:]
    topo_lon = ds.variables['lon'][:]
    topo_dat = ds.variables['Band1'][:]/1e3
    lons_topo, lats_topo = np.meshgrid(topo_lon,topo_lat)
    
    
    fig, ax = plt.subplots(figsize=(8,8)) 
    sc = ax.scatter(contour_lons, contour_lats, c=contour_times, cmap="viridis", alpha=0.5, s=1)
    cbar = plt.colorbar(sc, orientation="vertical")
    cbar.set_label("Time (Hours)")  # Change label to your preferred unit
    cbar.set_ticks(np.arange(len(all_times_hhmm)))  # Set ticks at step indices (0 to 7)
    cbar.set_ticklabels(all_times_hhmm)  # Map ticks to actual time values

    ax.plot(prov[:,0],prov[:,1],color='k'); 
    ax.set_xlim([-65.5,-62]); 
    ax.set_ylim([-33.5,-31])
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
    
    ax.plot(lon_radius, lat_radius, 'k', linewidth=0.8)
    ax.plot(lon_radius2, lat_radius2, 'k', linewidth=0.8)
        
    # agrego contorno de 500 y 1000m
    ax.contour(lons_topo, lats_topo, topo_dat, levels=[0.5,1], colors=['gray','gray'], linewidths=2.5)
    
    # Agregar los reportes! 
    samhi = pd.read_csv('/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/datos_base.csv')

    # Convert timestamp and extract only HH:MM
    samhi['formatted_time'] = pd.to_datetime(samhi['eventTime']).dt.strftime('%H:%M')

    # Plot reports
    #ax.scatter(samhi['longitude'].astype(float),samhi['latitude'].astype(float), s=50,marker = '^', color='grey', zorder=10, label='Granizo: ')

    # plot por severdad de granizo
    granizos_unk =samhi[samhi['hailMaxDiameterMM'].isna()]
    granizo_chico = samhi[(samhi['hailMaxDiameterMM']<20)]
    granizo_chico_unk = pd.concat((granizos_unk, granizo_chico))
    granizo_sev = samhi[(samhi['hailMaxDiameterMM']>=20) & (samhi['hailMaxDiameterMM']<50)]
    granizo_grande = samhi[samhi['hailMaxDiameterMM']>=50]

    ax.scatter(granizo_chico_unk['longitude'].astype(float),granizo_chico_unk['latitude'].astype(float), s=200, marker = '^', color='grey', zorder=10, label='hail')
    ax.scatter(granizo_sev['longitude'].astype(float),granizo_sev['latitude'].astype(float),             s=200, marker = '^', color='darkblue', zorder=12, label='hail > 2 cm' )
    ax.scatter(granizo_grande['longitude'].astype(float),granizo_grande['latitude'].astype(float),       s=200, marker = '^', color='darkred', zorder=12, label='hail > 5 cm')
    
    # add time label
    samhi_lons = samhi['longitude'].astype(float)
    samhi_lats = samhi['latitude'].astype(float)
    samhi_labels = samhi['formatted_time']
    ii=0; ax.text(samhi_lons[ii]-0.05, samhi_lats[ii]+0.05, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    ii=1; ax.text(samhi_lons[ii]-0.05, samhi_lats[ii]+0.05, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    ii=2; ax.text(samhi_lons[ii]-0.05, samhi_lats[ii]+0.05, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')

    # AGRUPAR: 21:30 aprox. 
    #ii=3; ax.text(samhi_lons[ii]+0.1, samhi_lats[ii]+0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    #ii=4; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii], samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    ii=5; ax.text(samhi_lons[ii], samhi_lats[ii]+0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    
    # Agrupar tambien 21:00 (super sur)
    ii=6; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii]+0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    ii=7; ax.text(samhi_lons[ii]+0.1, samhi_lats[ii]+0.08, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    #ii=12; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii]-0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    #ii=13; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii]-0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    
    
    #Agrupar 21:00 (super celda)
    ii=8; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii]+0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    #ii=9; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii]-0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    #ii=10; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii]-0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    #ii=11; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii]-0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    
    # del otro lado de las sierras
    ii=14; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii]+0.08, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')

    # reporte a las 00:00 (MCS? o error?)
    #ii=15; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii]-0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    #ii=16; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii]-0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
    #ii=17; ax.text(samhi_lons[ii]-0.1, samhi_lats[ii]-0.1, samhi_labels[ii], fontsize=12, color='black', ha='left', va='center')
        
    leg = ax.legend(loc='upper right',  fontsize=12, framealpha=0.80)
    # Title and save
    ax.grid()
    ax.set_title('COLMAX RMA1 evolution + SAMHI hail reports')    
    plt.show()
    fig.savefig(save_dir_compare+'/COLMAXevol_SAMHI.png', dpi=300,transparent=False,bbox_inches='tight')
    
            
    return



do_plot()

