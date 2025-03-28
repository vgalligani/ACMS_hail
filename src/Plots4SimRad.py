import config_folders
from SimRadAR import phidp_corr
import numpy as np 
import wradlib as wrl    
import matplotlib.pyplot as plt
from SimRadAR import functions, readWRF
import os
from Plots4Analysis  import pyplot_rings
from scipy import integrate
import pyart 
import seaborn as sns
import wrf
from scipy.interpolate import RectBivariateSpline
from pandas import to_datetime
from netCDF4 import num2date
from pyart.core import Radar

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def to_dBZ(var):
    dbZ = 10.*np.log10(var)
    return dbZ

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_WRF_qxs(EXP, time, mp):
    
    folders=config_folders.config_folders('yakaira')
    save_dir_compare=folders['save_dir_compare']
    
    ncfile = os.path.join(folders[EXP], 'wrfout_d02_2018-11-10_'+time+':00')
    [z_level, lat, lon, u, v, w, qr, qs, qc, qg, qi, qtotal] = readWRF.readWRFvariables(ncfile, mp)

    
    radar_latitude = -31.441389
    radar_longitude = -64.191944
    
    #------------------------------------
    # Plot on single figure: 
    # row1 => qx from WRF
    # row2 => WRF Zh colmax (C-band)
    # row2 => WRF Zh colmax (S-band)
    # row3 => WRF Zh S-band - C-band for different WRF levels of interest
    # row3 => transformation to radar sweeps
    fig, axes = plt.subplots(nrows=1, ncols=3, constrained_layout=True,figsize=[8*3,8])
    #------------------------------------
    # Plot WRF grid qx (integrated on the vertical) 
    #------------------------------------
    qr_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qr)) , z_level, axis=0)
    qs_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qs)) , z_level, axis=0)
    qg_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qg)) , z_level, axis=0)
    #
    pcm0 = axes[0].pcolormesh(lon, lat, qr_int, vmin=0,  vmax=1.5)
    pcm1 = axes[1].pcolormesh(lon, lat, qs_int, vmin= 0,  vmax=1.5)
    pcm2 = axes[2].pcolormesh(lon, lat, qg_int, vmin=0,  vmax=1.5)
    cbar = plt.colorbar(pcm2, ax=axes[2], shrink=1, label='qx [kg/m$^2$]')
    axes[0].grid(); axes[1].grid(); axes[2].grid()                                                            
    axes[0].set_title('WRF-WSM6 column q_rain')                                                  
    axes[1].set_title('WRF-WSM6 column q_snow')                                                        
    axes[2].set_title('WRF-WSM6 column q_grau')         
    
    for i in [120, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar_latitude,radar_longitude,i)
        for ax in [0,1,2]: 
            axes[ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8) 
    for i in [0,1,2]: 
        axes[i].set_xlim([-65.5,-62]); 
        axes[i].set_ylim([-35,-31])
        
    plt.show()
    fig.savefig(save_dir_compare+'/RadarSim'+'/WRF_qxs_'+EXP+'2018-11-10_'+time+'.png', 
                dpi=300,transparent=False, bbox_inches='tight')
    
        
    return

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def running_median(X, Y, total_bins):
    
    bins = np.linspace(np.nanmin(X[X!=-np.inf]),np.nanmax(X[X!=-np.inf]), total_bins)
    delta = bins[1]-bins[0]
    idx  = np.digitize(X,bins)
    running_median = [np.nanmedian(Y[idx==k]) for k in range(total_bins)]

    return bins, delta, running_median

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_medianZhhvsPol(Zh_t_central, Zdr_t_central, KDP_t_central, bandID, radar, EXP, mp, time):
    
    folders=config_folders.config_folders('yakaira')
    save_dir_compare=folders['save_dir_compare']
    
    # Plot the median of Zhh. vs. ZDR, y Zhh vs. KDP for C-band, S-band and observations
    Zhh = np.ravel(to_dBZ(Zh_t_central[bandID,:,:,:]))
    ZDR = np.ravel(Zdr_t_central[bandID,:,:,:])
    KDP = np.ravel(KDP_t_central[bandID,:,:,:])             # [deg/km]
    obsZhh = np.ravel(radar.fields['TH']['data'].data)
    obsZDR = np.ravel(radar.fields['ZDR']['data'].data)
    obsKDP = np.ravel(radar.fields['corrKDP']['data'].data)
        
    fig, axes = plt.subplots(nrows=1, ncols=2, constrained_layout=True,figsize=[16,8])
    [bins, delta, running_median_ZDR] = running_median(Zhh, ZDR, 40)
    [obins, odelta, orunning_median_ZDR] = running_median(obsZhh[obsZDR!=-9999.0], obsZDR[obsZDR!=-9999.0], 40)
    axes[0].plot(bins-delta/2, running_median_ZDR,'r--', lw=1.5, alpha=.8)
    axes[0].plot(bins-delta/2, orunning_median_ZDR,'k--', lw=1.5, alpha=.8)
    [bins, delta, running_median_KDP] = running_median(Zhh, KDP, 40)
    [obins, odelta, orunning_median_KDP] = running_median(obsZhh[obsZhh!=-9999.0], obsKDP[obsZhh!=-9999.0], 40)
    axes[1].plot(bins-delta/2, running_median_KDP,'r--', lw=1.5, alpha=.8, label='Sim. C-band')
    axes[1].plot(bins-delta/2, orunning_median_KDP,'k--', lw=1.5, alpha=.8, label='Obs.')
    axes[0].set_xlabel('Zhh (dBZ)')
    axes[0].set_ylabel('ZDR (dBZ)')
    axes[1].set_xlabel('Zhh (dBZ)')
    axes[1].set_ylabel('kdp ($^o$/km)')
    axes[0].grid(True)
    axes[1].grid(True)
    
    # agrego a s-band pero tambien hacerlo en funcion prolijo:
    Zhh = np.ravel(to_dBZ(Zh_t_central[2,:,:,:]))
    ZDR = np.ravel(Zdr_t_central[2,:,:,:])
    KDP = np.ravel(KDP_t_central[2,:,:,:])             # [deg/km]
    obsZhh = np.ravel(radar.fields['TH']['data'].data)
    obsZDR = np.ravel(radar.fields['ZDR']['data'].data)
    obsKDP = np.ravel(radar.fields['corrKDP']['data'].data)    
    [bins, delta, running_median_ZDR] = running_median(Zhh, ZDR, 40)
    axes[0].plot(bins-delta/2, running_median_ZDR,'b--', lw=1.5, alpha=.8)
    [bins, delta, running_median_KDP] = running_median(Zhh, KDP, 40)
    axes[1].plot(bins-delta/2, running_median_KDP,'b--', lw=1.5, alpha=.8, label='Sim. S-band')
    
    axes[1].legend()
    plt.show()
    fig.savefig(save_dir_compare+'/RadarSim/'+'RadarSim_kdpZDR'+EXP+'2018-11-10_'+time+'.png', dpi=300,transparent=False, bbox_inches='tight')
         
    
    return





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def create_cappi(
    radar,
    fields=None,
    height=2000,
    gatefilter=None,
    vel_field="velocity",
    same_nyquist=True,
    nyquist_vector_idx=0,
):
    """
    Create a Constant Altitude Plan Position Indicator (CAPPI) from radar data.

    Parameters
    ----------
    radar : Radar
        Py-ART Radar object containing the radar data.
    fields : list of str, optional
        List of radar fields to be used for creating the CAPPI.
        If None, all available fields will be used. Default is None.
    height : float, optional
        The altitude at which to create the CAPPI. Default is 2000 meters.
    gatefilter : GateFilter, optional
        A GateFilter object to apply masking/filtering to the radar data.
        Default is None.
    vel_field : str, optional
        The name of the velocity field to be used for determining the Nyquist velocity.
        Default is 'velocity'.
    same_nyquist : bool, optional
        Whether to only stack sweeps with the same Nyquist velocity.
        Default is True.
    nyquist_vector_idx : int, optional
        Index for the Nyquist velocity vector if `same_nyquist` is True.
        Default is 0.

    Returns
    -------
    Radar
        A Py-ART Radar object containing the CAPPI at the specified height.

    Notes
    -----
    CAPPI (Constant Altitude Plan Position Indicator) is a radar visualization
    technique that provides a horizontal view of meteorological data at a fixed altitude.
    Reference: https://glossary.ametsoc.org/wiki/Cappi

    Author
    ------
    Hamid Ali Syed (@syedhamidali)
    """

    if fields is None:
        fields = list(radar.fields.keys())

    # Initialize the first sweep as the reference
    first_sweep = 0

    # Initialize containers for the stacked data and nyquist velocities
    data_stack = []
    nyquist_stack = []

    # Process each sweep individually
    for sweep in range(radar.nsweeps):
        sweep_slice = radar.get_slice(sweep)
        try:
            nyquist = radar.get_nyquist_vel(sweep=sweep)
            nyquist = np.round(nyquist)
        except LookupError:
            print(
                f"Nyquist velocity unavailable for sweep {sweep}. Estimating using maximum velocity."
            )
            nyquist = radar.fields[vel_field]["data"][sweep_slice].max()

        sweep_data = {}

        for field in fields:
            data = radar.get_field(sweep, field)

            # Apply gatefilter if provided
            if gatefilter is not None:
                data = np.ma.masked_array(
                    data, gatefilter.gate_excluded[sweep_slice, :]
                )
            time = radar.time["data"][sweep_slice]

            # Extract and sort azimuth angles
            azimuth = radar.azimuth["data"][sweep_slice]
            azimuth_sorted_idx = np.argsort(azimuth)
            azimuth = azimuth[azimuth_sorted_idx]
            data = data[azimuth_sorted_idx]

            # Store initial lat/lon for reordering
            if sweep == first_sweep:
                azimuth_final = azimuth
                time_final = time
            else:
                # Interpolate data for consistent azimuth ordering across sweeps
                interpolator = RectBivariateSpline(azimuth, radar.range["data"], data)
                data = interpolator(azimuth_final, radar.range["data"])

            sweep_data[field] = data[np.newaxis, :, :]

        data_stack.append(sweep_data)
        nyquist_stack.append(nyquist)

    nyquist_stack = np.array(nyquist_stack)

    # Filter for sweeps with similar Nyquist velocities
    if same_nyquist:
        nyquist_range = nyquist_stack[nyquist_vector_idx]
        nyquist_mask = np.abs(nyquist_stack - nyquist_range) <= 1
        data_stack = [
            sweep_data for i, sweep_data in enumerate(data_stack) if nyquist_mask[i]
        ]

    # Generate CAPPI for each field using data_stack
    fields_data = {}
    for field in fields:
        data_3d = np.concatenate(
            [sweep_data[field] for sweep_data in data_stack], axis=0
        )

        # Sort azimuth for all sweeps
        dim0 = data_3d.shape[1:]
        azimuths = np.linspace(0, 359, dim0[0])
        elevation_angles = radar.fixed_angle["data"][: data_3d.shape[0]]
        ranges = radar.range["data"]

        theta = (450 - azimuths) % 360
        THETA, PHI, R = np.meshgrid(theta, elevation_angles, ranges)
        Z = R * np.sin(PHI * np.pi / 180)

        # Extract the data slice corresponding to the requested height
        height_idx = np.argmin(np.abs(Z - height), axis=0)
        CAPPI = np.array(
            [
                data_3d[height_idx[j, i], j, i]
                for j in range(dim0[0])
                for i in range(dim0[1])
            ]
        ).reshape(dim0)

        # Retrieve units and handle case where units might be missing
        units = radar.fields[field].get("units", "").lower()

        # Determine valid_min and valid_max based on units
        if units == "dbz":
            valid_min, valid_max = -10, 80
        elif units in ["m/s", "meters per second"]:
            valid_min, valid_max = -100, 100
        elif units == "db":
            valid_min, valid_max = -7.9, 7.9
        else:
            # If units are not found or don't match known types, set default values or skip masking
            valid_min, valid_max = None, None

        # If valid_min or valid_max are still None, set them to conservative defaults or skip
        if valid_min is None:
            print(f"Warning: valid_min not set for {field}, using default of -1e10")
            valid_min = -1e10  # Conservative default
        if valid_max is None:
            print(f"Warning: valid_max not set for {field}, using default of 1e10")
            valid_max = 1e10  # Conservative default

        # Apply valid_min and valid_max masking
        if valid_min is not None:
            CAPPI = np.ma.masked_less(CAPPI, valid_min)
        if valid_max is not None:
            CAPPI = np.ma.masked_greater(CAPPI, valid_max)

        # Convert to masked array with the specified fill value
        CAPPI.set_fill_value(radar.fields[field].get("_FillValue", np.nan))
        CAPPI = np.ma.masked_invalid(CAPPI)
        CAPPI = np.ma.masked_outside(CAPPI, valid_min, valid_max)

        fields_data[field] = {
            "data": CAPPI,
            "units": radar.fields[field]["units"],
            "long_name": f"CAPPI {field} at {height} meters",
            "comment": f"CAPPI {field} calculated at a height of {height} meters",
            "_FillValue": radar.fields[field].get("_FillValue", np.nan),
        }

    # Set the elevation to zeros for CAPPI
    elevation_final = np.zeros(dim0[0], dtype="float32")

    # Since we are using the whole volume scan, report mean time
    try:
        dtime = to_datetime(
            num2date(radar.time["data"], radar.time["units"]).astype(str),
            format="ISO8601",
        )
    except ValueError:
        dtime = to_datetime(
            num2date(radar.time["data"], radar.time["units"]).astype(str)
        )
    dtime = dtime.mean()

    time = radar.time.copy()
    time["data"] = time_final
    time["mean"] = dtime

    # Create the Radar object with the new CAPPI data
    return Radar(
        time=radar.time.copy(),
        _range=radar.range.copy(),
        fields=fields_data,
        metadata=radar.metadata.copy(),
        scan_type=radar.scan_type,
        latitude=radar.latitude.copy(),
        longitude=radar.longitude.copy(),
        altitude=radar.altitude.copy(),
        sweep_number=radar.sweep_number.copy(),
        sweep_mode=radar.sweep_mode.copy(),
        fixed_angle=radar.fixed_angle.copy(),
        sweep_start_ray_index=radar.sweep_start_ray_index.copy(),
        sweep_end_ray_index=radar.sweep_end_ray_index.copy(),
        azimuth=radar.azimuth.copy(),
        elevation={"data": elevation_final},
        instrument_parameters=radar.instrument_parameters,
    )

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def ZH1km(zh, z, alt):
    
    return wrf.interplevel(zh, z, alt)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_WRF_grid_all_ZH1km(radar, EXP, time, mp, bandID, Zh_r_WRFGRID, Zh_s_WRFGRID, Zh_g_WRFGRID, Zh_tot_WRFGRID):
    
    folders=config_folders.config_folders('yakaira')
    save_dir_compare=folders['save_dir_compare']
    
    ncfile = os.path.join(folders[EXP], 'wrfout_d02_2018-11-10_'+time+':00')
    [z_level, lat, lon, u, v, w, qr, qs, qc, qg, qi, qtotal] = readWRF.readWRFvariables(ncfile, mp)
    
    #------------------------------------
    # Plot on single figure: 
    # row1 => qx from WRF
    # row2 => WRF Zh colmax (C-band)
    # row2 => WRF Zh colmax (S-band)
    # row3 => WRF Zh S-band - C-band for different WRF levels of interest
    # row3 => transformation to radar sweeps
    fig, axes = plt.subplots(nrows=4, ncols=4, constrained_layout=True,figsize=[26,28])
    #------------------------------------
     # Plot WRF grid qx (integrated on the vertical) 
    #------------------------------------
    qr_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qr)) , z_level, axis=0)
    qs_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qs)) , z_level, axis=0)
    qg_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qg)) , z_level, axis=0)
    #
    pcm0 = axes[0,0].pcolormesh(lon, lat, qr_int, vmin=0,  vmax=1.5)
    pcm1 = axes[0,1].pcolormesh(lon, lat, qs_int, vmin= 0,  vmax=1.5)
    pcm2 = axes[0,2].pcolormesh(lon, lat, qg_int, vmin=0,  vmax=1.5)
    cbar = plt.colorbar(pcm2, ax=axes[0,2], shrink=1, label='qx [kg/m$^2$]')
    axes[0,0].grid(); axes[0,1].grid(); axes[0,2].grid()                                                            
    axes[0,0].set_title('WRF-WSM6 column q_rain')                                                  
    axes[0,1].set_title('WRF-WSM6 column q_snow')                                                        
    axes[0,2].set_title('WRF-WSM6 column q_grau')         
    for i in [60, 120, 180, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [0,1,2]: 
            axes[0,ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    

    # Plot calculated WRF grid COLMAX(Zh) for each species at C-band 
    #------------------------------------
    pcm0 = axes[1,0].pcolormesh(lon, lat, to_dBZ(ZH1km(Zh_r_WRFGRID[bandID,:,:,:], z_level, 1000)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm1 = axes[1,1].pcolormesh(lon, lat, to_dBZ(ZH1km(Zh_s_WRFGRID[bandID,:,:,:], z_level, 1000)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm2 = axes[1,2].pcolormesh(lon, lat, to_dBZ(ZH1km(Zh_g_WRFGRID[bandID,:,:,:], z_level, 1000)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm3 = axes[1,3].pcolormesh(lon, lat, to_dBZ(ZH1km(Zh_tot_WRFGRID[bandID,:,:,:], z_level, 1000)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    cbar = plt.colorbar(pcm3, ax=axes[1,3], shrink=1, label='Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    axes[1,0].grid(); axes[1,1].grid(); axes[1,2].grid(); axes[1,3].grid()                                                                 
    axes[1,0].set_title('C-band 1km Zh_rain (at WRF grid)')                                                        
    axes[1,1].set_title('C-band 1km Zh_snow (at WRF grid)')                                                        
    axes[1,2].set_title('C-band 1km Zh_grau (at WRF grid)')    
    axes[1,3].set_title('C-band 1km Zh_tot (at WRF grid)')    
    for i in [60, 120, 180, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [0,1,2,3]: 
            axes[1,ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    

    # Plot calculated WRF grid COLMAX(Zh) for each species AT S-BAND
    #------------------------------------
    pcm0 = axes[2,0].pcolormesh(lon, lat, to_dBZ(ZH1km(Zh_r_WRFGRID[0,:,:,:], z_level, 1000)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm1 = axes[2,1].pcolormesh(lon, lat, to_dBZ(ZH1km(Zh_s_WRFGRID[0,:,:,:], z_level, 1000)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm2 = axes[2,2].pcolormesh(lon, lat, to_dBZ(ZH1km(Zh_g_WRFGRID[0,:,:,:], z_level, 1000)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm3 = axes[2,3].pcolormesh(lon, lat, to_dBZ(ZH1km(Zh_tot_WRFGRID[0,:,:,:], z_level, 1000)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    cbar = plt.colorbar(pcm3, ax=axes[2,3], shrink=1, label='Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    axes[2,0].grid(); axes[2,1].grid(); axes[2,2].grid(); axes[2,3].grid()                                                                  
    axes[2,0].set_title('S-band 1km Zh_rain (at WRF grid)')                                                        
    axes[2,1].set_title('S-band 1km Zh_snow (at WRF grid)')                                                        
    axes[2,2].set_title('S-band 1km Zh_grau (at WRF grid)')    
    axes[2,3].set_title('S-band 1km Zh_tot (at WRF grid)')    
    for i in [60, 120, 180, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [0,1,2,3]: 
            axes[2,ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    

        
    # Plot calculated WRF grid COLMAX(Zh) S-band - C-band diff for each species 
    #------------------------------------
    Zh1_qr = to_dBZ(ZH1km(Zh_r_WRFGRID[0,:,:,:], z_level, 1000)) - to_dBZ(ZH1km(Zh_r_WRFGRID[bandID,:,:,:], z_level, 1000))
    Zh1_qs = to_dBZ(ZH1km(Zh_s_WRFGRID[0,:,:,:], z_level, 6000)) - to_dBZ(ZH1km(Zh_s_WRFGRID[bandID,:,:,:], z_level, 6000))
    Zh1_qg = to_dBZ(ZH1km(Zh_g_WRFGRID[0,:,:,:], z_level, 6000)) - to_dBZ(ZH1km(Zh_g_WRFGRID[bandID,:,:,:], z_level, 6000))
    Zh1_tot = to_dBZ(ZH1km(Zh_tot_WRFGRID[0,:,:,:], z_level, 1000)) - to_dBZ(ZH1km(Zh_tot_WRFGRID[bandID,:,:,:], z_level, 1000))
    
    pcm0 = axes[3,0].pcolormesh(lon, lat, Zh1_qr, cmap=functions.colormaps('ref'), vmin=-10, vmax=10)
    pcm1 = axes[3,1].pcolormesh(lon, lat, Zh1_qs, cmap=functions.colormaps('ref'), vmin=-10, vmax=10)
    pcm2 = axes[3,2].pcolormesh(lon, lat, Zh1_qg, cmap=functions.colormaps('ref'), vmin=-10, vmax=10)
    pcm3 = axes[3,3].pcolormesh(lon, lat, Zh1_tot, cmap=functions.colormaps('ref'), vmin=-10, vmax=10)
    cbar = plt.colorbar(pcm3, ax=axes[3,3], shrink=1, label='Zh [dBZ]', ticks = np.arange(-10,10.01,2))
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    axes[3,0].grid(); axes[3,1].grid(); axes[3,2].grid(); axes[3,3].grid()                                                            
    axes[3,0].set_title('(S-band - C-band) Zh_rain (at WRF grid, zlevel=1km)')                                                        
    axes[3,1].set_title('(S-band - C-band) Zh_snow (at WRF grid, zlevel=6km)')                                                        
    axes[3,2].set_title('(S-band - C-band) Zh_grau (at WRF grid, zlevel=6klm)')    
    axes[3,3].set_title('(S-band - C-band) Zh_tot (at WRF grid, zlevel=1km)')    
    for i in [60, 120, 180, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [0,1,2,3]: 
            axes[3,ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)  
    
    #------------------------------------
    # Add for reference radar COLMAX observation
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()

    # Create CAPPI at 1,000 meters for the 'reflectivity' field
    radar_cappi = create_cappi(radar, fields=["TH"], vel_field="VRAD", height=1000, gatefilter=gatefilter)
    # compz       = pyart.retrieve.composite_reflectivity(radar, field="TH", gatefilter=gatefilter)
    start_index = radar.sweep_start_ray_index['data'][0]
    end_index   = radar.sweep_end_ray_index['data'][0]
    lats        = radar.gate_latitude['data'][start_index:end_index+1]
    lons        = radar.gate_longitude['data'][start_index:end_index+1]
    pcm0 = axes[0,3].pcolormesh(lons, lats, radar_cappi.fields['TH']['data'], cmap=functions.colormaps('ref'), vmin=0,  vmax=70)
    cbar = plt.colorbar(pcm0, ax=axes[0,3], shrink=1, label='Zh CAPPI 1km [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    for i in [60, 120, 180, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        axes[0,3].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    
    axes[0,3].set_title('RMA1 1km CAPPI')    
    
    plt.show()
    fig.savefig(save_dir_compare+'/RadarSim/'+'/RadarSim_WRFgrid_CAPPI'+EXP+'2018-11-10_'+time+'.png', 
                dpi=300,transparent=False, bbox_inches='tight')     
    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_WRF_grid_all(radar, EXP, time, mp, bandID, Zh_r_WRFGRID, Zh_s_WRFGRID, Zh_g_WRFGRID, Zh_tot_WRFGRID):
    
    folders=config_folders.config_folders('yakaira')
    save_dir_compare=folders['save_dir_compare']
    
    ncfile = os.path.join(folders[EXP], 'wrfout_d02_2018-11-10_'+time+':00')
    [z_level, lat, lon, u, v, w, qr, qs, qc, qg, qi, qtotal] = readWRF.readWRFvariables(ncfile, mp)
    
    #------------------------------------
    # Plot on single figure: 
    # row1 => qx from WRF
    # row2 => WRF Zh colmax (C-band)
    # row2 => WRF Zh colmax (S-band)
    # row3 => WRF Zh S-band - C-band for different WRF levels of interest
    # row3 => transformation to radar sweeps
    fig, axes = plt.subplots(nrows=4, ncols=4, constrained_layout=True,figsize=[26,28])
    #------------------------------------
     # Plot WRF grid qx (integrated on the vertical) 
    #------------------------------------
    qr_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qr)) , z_level, axis=0)
    qs_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qs)) , z_level, axis=0)
    qg_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qg)) , z_level, axis=0)
    #
    pcm0 = axes[0,0].pcolormesh(lon, lat, qr_int, vmin=0,  vmax=1.5)
    pcm1 = axes[0,1].pcolormesh(lon, lat, qs_int, vmin= 0,  vmax=1.5)
    pcm2 = axes[0,2].pcolormesh(lon, lat, qg_int, vmin=0,  vmax=1.5)
    cbar = plt.colorbar(pcm2, ax=axes[0,2], shrink=1, label='qx [kg/m$^2$]')
    axes[0,0].grid(); axes[0,1].grid(); axes[0,2].grid()                                                            
    axes[0,0].set_title('WRF-WSM6 column q_rain')                                                  
    axes[0,1].set_title('WRF-WSM6 column q_snow')                                                        
    axes[0,2].set_title('WRF-WSM6 column q_grau')         
    for i in [60, 120, 180, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [0,1,2]: 
            axes[0,ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    

    # Plot calculated WRF grid COLMAX(Zh) for each species at C-band 
    #------------------------------------
    pcm0 = axes[1,0].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_r_WRFGRID[bandID,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm1 = axes[1,1].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_s_WRFGRID[bandID,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm2 = axes[1,2].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_g_WRFGRID[bandID,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm3 = axes[1,3].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_tot_WRFGRID[bandID,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    cbar = plt.colorbar(pcm3, ax=axes[1,3], shrink=1, label='Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    axes[1,0].grid(); axes[1,1].grid(); axes[1,2].grid(); axes[1,3].grid()                                                                 
    axes[1,0].set_title('C-band COLMAX Zh_rain (at WRF grid)')                                                        
    axes[1,1].set_title('C-band COLMAX Zh_snow (at WRF grid)')                                                        
    axes[1,2].set_title('C-band COLMAX Zh_grau (at WRF grid)')    
    axes[1,3].set_title('C-band COLMAX Zh_tot (at WRF grid)')    
    for i in [60, 120, 180, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [0,1,2,3]: 
            axes[1,ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    

    # Plot calculated WRF grid COLMAX(Zh) for each species AT S-BAND
    #------------------------------------
    pcm0 = axes[2,0].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_r_WRFGRID[0,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm1 = axes[2,1].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_s_WRFGRID[0,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm2 = axes[2,2].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_g_WRFGRID[0,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm3 = axes[2,3].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_tot_WRFGRID[0,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    cbar = plt.colorbar(pcm3, ax=axes[2,3], shrink=1, label='Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    axes[2,0].grid(); axes[2,1].grid(); axes[2,2].grid(); axes[2,3].grid()                                                                  
    axes[2,0].set_title('S-band COLMAX Zh_rain (at WRF grid)')                                                        
    axes[2,1].set_title('S-band COLMAX Zh_snow (at WRF grid)')                                                        
    axes[2,2].set_title('S-band COLMAX Zh_grau (at WRF grid)')    
    axes[2,3].set_title('S-band COLMAX Zh_tot (at WRF grid)')    
    for i in [60, 120, 180, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [0,1,2,3]: 
            axes[2,ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    


    # Plot calculated WRF grid COLMAX(Zh) S-band - C-band diff for each species 
    #------------------------------------
    wrf_zlevelNr = 1
    Zh1_qr = to_dBZ(Zh_r_WRFGRID[0,wrf_zlevelNr,:,:]) - to_dBZ(Zh_r_WRFGRID[bandID,wrf_zlevelNr,:,:])
    Zh1_qs = to_dBZ(Zh_s_WRFGRID[0,18,:,:]) - to_dBZ( Zh_s_WRFGRID[bandID,18,:,:])
    Zh1_qg = to_dBZ(Zh_g_WRFGRID[0,18,:,:]) - to_dBZ(Zh_g_WRFGRID[bandID,18,:,:])
    Zh1_tot = to_dBZ(Zh_tot_WRFGRID[0,wrf_zlevelNr,:,:]) - to_dBZ(Zh_tot_WRFGRID[bandID,wrf_zlevelNr,:,:])
    pcm0 = axes[3,0].pcolormesh(lon, lat, to_dBZ(Zh1_qr), cmap=functions.colormaps('ref'), vmin=-10, vmax=10)
    pcm1 = axes[3,1].pcolormesh(lon, lat, to_dBZ(Zh1_qs), cmap=functions.colormaps('ref'), vmin=-10, vmax=10)
    pcm2 = axes[3,2].pcolormesh(lon, lat, to_dBZ(Zh1_qg), cmap=functions.colormaps('ref'), vmin=-10, vmax=10)
    pcm3 = axes[3,3].pcolormesh(lon, lat, to_dBZ(Zh1_tot), cmap=functions.colormaps('ref'), vmin=-10, vmax=10)
    cbar = plt.colorbar(pcm3, ax=axes[3,3], shrink=1, label='Zh [dBZ]', ticks = np.arange(-10,10.01,2))
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    axes[3,0].grid(); axes[3,1].grid(); axes[3,2].grid(); axes[3,3].grid()                                                            
    axes[3,0].set_title('(S-band - C-band) Zh_rain (at WRF grid, zlevel=1)')                                                        
    axes[3,1].set_title('(S-band - C-band) Zh_snow (at WRF grid, zlevel=18)')                                                        
    axes[3,2].set_title('(S-band - C-band) Zh_grau (at WRF grid, zlevel=18)')    
    axes[3,3].set_title('(S-band - C-band) Zh_tot (at WRF grid, zlevel=1)')    
    for i in [60, 120, 180, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [0,1,2,3]: 
            axes[3,ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)  
    
    #------------------------------------
    # Add for reference radar COLMAX observation
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    compz       = pyart.retrieve.composite_reflectivity(radar, field="TH", gatefilter=gatefilter)
    # compz.fields['composite_reflectivity']['data'] == (360, 521)
    start_index = radar.sweep_start_ray_index['data'][0]
    end_index   = radar.sweep_end_ray_index['data'][0]
    lats        = radar.gate_latitude['data'][start_index:end_index+1]
    lons        = radar.gate_longitude['data'][start_index:end_index+1]
    pcm0 = axes[0,3].pcolormesh(lons, lats, compz.fields['composite_reflectivity']['data'], cmap=functions.colormaps('ref'), vmin=0,  vmax=70)
    cbar = plt.colorbar(pcm0, ax=axes[0,3], shrink=1, label='ColMAX Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    for i in [60, 120, 180, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        axes[0,3].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    
    axes[0,3].set_title('RMA1 colmax')    
    
    plt.show()
    fig.savefig(save_dir_compare+'/RadarSim/'+'/RadarSim_WRFgrid_'+EXP+'2018-11-10_'+time+'.png', 
                dpi=300,transparent=False, bbox_inches='tight')     
    return


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_radar_grid(radar, EXP, time): 

    
    folders=config_folders.config_folders('yakaira')
    save_dir_compare=folders['save_dir_compare']
    
    # PHIDP CORRECTION 
    PHIORIG = radar.fields['PHIDP']['data'].copy() 
    mask = radar.fields['PHIDP']['data'].data.copy()    
    mask[:] = False
    PHIORIG.mask = mask
    radar.add_field_like('PHIDP', 'PHIDP', PHIORIG, replace_existing=True)     
    sys_phase = phidp_corr.get_sys_phase_simple(radar)
    dphi, uphi, corr_phidp = phidp_corr.correct_phidp(radar, sys_phase, 280)
    corr_phidp[corr_phidp==0] = np.nan
    calculated_KDP = wrl.dp.kdp_from_phidp(corr_phidp, winlen=7, dr=(radar.range['data'][1]-radar.range['data'][0])/1e3, 
    					   method='lanczos_conv', skipna=True)	
    radar.add_field_like('RHOHV','corrKDP', calculated_KDP, replace_existing=True)
    
    # Plot on single figure: 
    fig, axes = plt.subplots(nrows=1, ncols=4, constrained_layout=True,figsize=[8*4,8])
           
    # Plot REAL radar observatrions and last col is colmax
    #------------------------------------
    for nlev in [0,1,2]:
        start_index = radar.sweep_start_ray_index['data'][nlev]
        end_index   = radar.sweep_end_ray_index['data'][nlev]
        lats        = radar.gate_latitude['data'][start_index:end_index]
        lons        = radar.gate_longitude['data'][start_index:end_index]
        azimuths    = radar.azimuth['data'][start_index:end_index]
        TH          = radar.fields['TH']['data'][start_index:end_index]
        pcm1 = axes[nlev].pcolormesh(lons, lats, TH, cmap=functions.colormaps('ref'), vmin=0, vmax=70)
        axes[nlev].grid()
        axes[nlev].set_title('RMA1 Zh (elev='+str(radar.fixed_angle['data'][nlev])+'$^o$)')                                                        
    cbar = plt.colorbar(pcm1, ax=axes[nlev], shrink=1, label='Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    for i in [120, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [0,1,2]: 
            axes[ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    
    axes[0].text(-66.9,-29.5,'Radar rings at: \n120, 240 km')               
             
    # Add for reference radar COLMAX observation
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    compz       = pyart.retrieve.composite_reflectivity(radar, field="TH", gatefilter=gatefilter)
    # compz.fields['composite_reflectivity']['data'] == (360, 521)
    start_index = radar.sweep_start_ray_index['data'][0]
    end_index   = radar.sweep_end_ray_index['data'][0]
    lats        = radar.gate_latitude['data'][start_index:end_index+1]
    lons        = radar.gate_longitude['data'][start_index:end_index+1]
    pcm0 = axes[3].pcolormesh(lons, lats, compz.fields['composite_reflectivity']['data'], cmap=functions.colormaps('ref'), vmin=0,  vmax=70)
    cbar = plt.colorbar(pcm0, ax=axes[3], shrink=1, label='ColMAX Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    for i in [120, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        axes[3].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    
    axes[3].set_title('RMA1 colmax') 
    
    plt.show()
    fig.savefig(save_dir_compare+'/RadarSim/'+'/RadarOBS_'+EXP+'2018-11-10_'+time+'.png', 
                dpi=300,transparent=False, bbox_inches='tight')    
    
    return radar




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_radar_grid_all(lon, lat,radar, EXP, mp, Zh_r_central, Zh_s_central, Zh_g_central, Zh_t_central, bandID, time):
    
    theta_radar = (0.5, 0.9, 1.3, 1.9, 2.3, 3, 3.5, 5, 6.9, 9.1, 11.8, 15.1)                                       

    
    folders=config_folders.config_folders('yakaira')
    save_dir_compare=folders['save_dir_compare']
    
    # PHIDP CORRECTION 
    PHIORIG = radar.fields['PHIDP']['data'].copy() 
    mask = radar.fields['PHIDP']['data'].data.copy()    
    mask[:] = False
    PHIORIG.mask = mask
    radar.add_field_like('PHIDP', 'PHIDP', PHIORIG, replace_existing=True)     
    sys_phase = phidp_corr.get_sys_phase_simple(radar)
    dphi, uphi, corr_phidp = phidp_corr.correct_phidp(radar, sys_phase, 280)
    # A corr_phidp le tengo que poner el 0 == np.nan? 
    corr_phidp[corr_phidp==0] = np.nan
    calculated_KDP = wrl.dp.kdp_from_phidp(corr_phidp, winlen=7, dr=(radar.range['data'][1]-radar.range['data'][0])/1e3, 
    					   method='lanczos_conv', skipna=True)	
    radar.add_field_like('RHOHV','corrKDP', calculated_KDP, replace_existing=True)
    
    # Plot on single figure: 
    fig, axes = plt.subplots(nrows=5, ncols=4, constrained_layout=True,figsize=[26,28])
           
    # Plot REAL radar observatrions and last col is colmax
    #------------------------------------
    for nlev in [0,1,2]:
        start_index = radar.sweep_start_ray_index['data'][nlev]
        end_index   = radar.sweep_end_ray_index['data'][nlev]
        lats        = radar.gate_latitude['data'][start_index:end_index]
        lons        = radar.gate_longitude['data'][start_index:end_index]
        azimuths    = radar.azimuth['data'][start_index:end_index]
        TH          = radar.fields['TH']['data'][start_index:end_index]
        pcm1 = axes[0,nlev].pcolormesh(lons, lats, TH, cmap=functions.colormaps('ref'), vmin=0, vmax=70)
        axes[0,nlev].grid()
        axes[0,nlev].set_title('RMA1 Zh (elev='+str(radar.fixed_angle['data'][nlev])+'$^o$)')                                                        
    cbar = plt.colorbar(pcm1, ax=axes[0,nlev], shrink=1, label='Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    for i in [120, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [0,1,2]: 
            axes[0,ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    
    axes[0,0].text(-66.9,-29.5,'Radar rings at: \n120, 240 km')               
             
    # Add for reference radar COLMAX observation
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    compz       = pyart.retrieve.composite_reflectivity(radar, field="TH", gatefilter=gatefilter)
    # compz.fields['composite_reflectivity']['data'] == (360, 521)
    start_index = radar.sweep_start_ray_index['data'][0]
    end_index   = radar.sweep_end_ray_index['data'][0]
    lats        = radar.gate_latitude['data'][start_index:end_index+1]
    lons        = radar.gate_longitude['data'][start_index:end_index+1]
    pcm0 = axes[0,3].pcolormesh(lons, lats, compz.fields['composite_reflectivity']['data'], cmap=functions.colormaps('ref'), vmin=0,  vmax=70)
    cbar = plt.colorbar(pcm0, ax=axes[0,3], shrink=1, label='ColMAX Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    for i in [120, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        axes[0,3].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    
    axes[0,3].set_title('RMA1 colmax') 
    
    # Plot Zh transformed to radar grid for the first three elevations 
    #------------------------------------
    for i in [0,1,2]: 
        pcm0 = axes[1,i].pcolormesh(lon, lat, to_dBZ(Zh_r_central[bandID,i,:,:]), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
        pcm1 = axes[2,i].pcolormesh(lon, lat, to_dBZ(Zh_s_central[bandID,i,:,:]), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
        pcm2 = axes[3,i].pcolormesh(lon, lat, to_dBZ(Zh_g_central[bandID,i,:,:]), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
        axes[1,i].set_title('C-band Zh_rain (elev='+ str(theta_radar[i])+'$^o$)')                                                  
        axes[2,i].set_title('C-band Zh_snow (elev='+ str(theta_radar[i])+'$^o$)')    
        axes[3,i].set_title('C-band Zh_grau (elev='+ str(theta_radar[i])+'$^o$)')                                                        
        axes[1,i].grid(); axes[2,i].grid(); axes[3,i].grid() 
        for ii in [120, 240]: 
            [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],ii)
            axes[1,i].plot(lon_radius, lat_radius, 'k', linewidth=0.8)   
            axes[2,i].plot(lon_radius, lat_radius, 'k', linewidth=0.8)   
            axes[3,i].plot(lon_radius, lat_radius, 'k', linewidth=0.8)   
    
    cbar = plt.colorbar(pcm2, ax=axes[1,2], shrink=1, label='Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    
    cbar = plt.colorbar(pcm2, ax=axes[2,2], shrink=1, label='Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    
    cbar = plt.colorbar(pcm2, ax=axes[3,2], shrink=1, label='Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    
    
    #------------------------------------        
    # Plot Zh_tot transformed to radar grid for the first three elevations 
    #------------------------------------
    for i in [0,1,2]: 
        pcm0 = axes[4,i].pcolormesh(lon, lat, to_dBZ(Zh_t_central[bandID,i,:,:]), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
        axes[4,i].set_title('C-band Zh_tot (elev='+ str(theta_radar[i])+'$^o$)')                                                  
        axes[4,i].grid(); axes[4,i].grid(); axes[4,i].grid() 
        for ii in [120, 240]: 
            [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],ii)
            axes[4,i].plot(lon_radius, lat_radius, 'k', linewidth=0.8)   
    cbar = plt.colorbar(pcm0, ax=axes[4,2], shrink=1, label='Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    
    
    #------------------------------------        
    # Plot COLMAXZh transformed to radar grid for reference on the last column 
    #------------------------------------
    pcm0 = axes[1,3].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_r_central[bandID,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm1 = axes[2,3].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_s_central[bandID,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm2 = axes[3,3].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_g_central[bandID,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    pcm3 = axes[4,3].pcolormesh(lon, lat, to_dBZ(np.nanmax(Zh_t_central[bandID,:,:,:], axis=0)), cmap=functions.colormaps('ref'), vmin=0, vmax=70)
    for i in [1,2,3,4]:
        cbar = plt.colorbar(pcm2, ax=axes[i,3], shrink=1, label='Zh [dBZ]')
        cbar.cmap.set_under('white')
        cbar.cmap.set_over('white')
        axes[i,3].grid()                                                            
    axes[1,3].set_title('C-band COLMAX Zh_rain (at radar grid)')                                                        
    axes[2,3].set_title('C-band COLMAX Zh_snow (at radar grid)')                                                        
    axes[3,3].set_title('C-band COLMAX Zh_grau (at radar grid)')    
    axes[4,3].set_title('C-band COLMAX Zh_tot (at radar grid)')    
    for i in [120, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [1,2,3,4]: 
            axes[ax,3].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    
    
    plt.show()
    fig.savefig(save_dir_compare+'/RadarSim/'+'/RadarSim_RADARgrid_'+EXP+'2018-11-10_'+time+'.png', 
                dpi=300,transparent=False, bbox_inches='tight')     
    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_CvsS(lon, lat, radar, Zh_r_central, Zh_s_central, Zh_g_central, Zh_t_central, EXP, bandID, time ): 

    theta_radar = (0.5, 0.9, 1.3, 1.9, 2.3, 3, 3.5, 5, 6.9, 9.1, 11.8, 15.1)                                       

    folders=config_folders.config_folders('yakaira')
    save_dir_compare=folders['save_dir_compare']

    # Plot on single figure C-band vs. S-band on single grid
    fig, axes = plt.subplots(nrows=5, ncols=4, constrained_layout=True,figsize=[26,28])
           
    # Plot REAL radar observatrions and last col is colmax
    #------------------------------------
    for nlev in [0,1,2]:
        start_index = radar.sweep_start_ray_index['data'][nlev]
        end_index   = radar.sweep_end_ray_index['data'][nlev]
        lats        = radar.gate_latitude['data'][start_index:end_index]
        lons        = radar.gate_longitude['data'][start_index:end_index]
        azimuths    = radar.azimuth['data'][start_index:end_index]
        TH          = radar.fields['TH']['data'][start_index:end_index]
        pcm1 = axes[0,nlev].pcolormesh(lons, lats, TH, cmap=functions.colormaps('ref'), vmin=0, vmax=70)
        axes[0,nlev].grid()
        axes[0,nlev].set_title('RMA1 Zh (elev='+str(radar.fixed_angle['data'][nlev])+'$^o$)')                                                        
    cbar = plt.colorbar(pcm1, ax=axes[0,nlev], shrink=1, label='Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    for i in [120, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        for ax in [0,1,2]: 
            axes[0,ax].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    
    axes[0,0].text(-66.9,-29.5,'Radar rings at: \n120, 240 km')               
             
    # Add for reference radar COLMAX observation
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    compz       = pyart.retrieve.composite_reflectivity(radar, field="TH", gatefilter=gatefilter)
    # compz.fields['composite_reflectivity']['data'] == (360, 521)
    start_index = radar.sweep_start_ray_index['data'][0]
    end_index   = radar.sweep_end_ray_index['data'][0]
    lats        = radar.gate_latitude['data'][start_index:end_index+1]
    lons        = radar.gate_longitude['data'][start_index:end_index+1]
    pcm0 = axes[0,3].pcolormesh(lons, lats, compz.fields['composite_reflectivity']['data'], cmap=functions.colormaps('ref'), vmin=0,  vmax=70)
    cbar = plt.colorbar(pcm0, ax=axes[0,3], shrink=1, label='ColMAX Zh [dBZ]')
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    for i in [120, 240]: 
        [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],i)
        axes[0,3].plot(lon_radius, lat_radius, 'k', linewidth=0.8)    
    axes[0,3].set_title('RMA1 colmax')   
    
                    
    # Plot Zh transformed to radar grid for the first three elevations 
    #------------------------------------
    for i in [0,1,2]: 
        diff_r = to_dBZ(Zh_r_central[0,i,:,:]) -  to_dBZ(Zh_r_central[bandID,i,:,:])
        diff_s = to_dBZ(Zh_s_central[0,i,:,:]) -  to_dBZ(Zh_s_central[bandID,i,:,:])
        diff_g = to_dBZ(Zh_g_central[0,i,:,:]) -  to_dBZ(Zh_g_central[bandID,i,:,:]) 
        pcm0 = axes[1,i].pcolormesh(lon, lat, diff_r, cmap=functions.colormaps('ref'), vmin=-5, vmax=5)
        pcm1 = axes[2,i].pcolormesh(lon, lat, diff_s, cmap=functions.colormaps('ref'), vmin=-5, vmax=5)
        pcm2 = axes[3,i].pcolormesh(lon, lat, diff_g, cmap=functions.colormaps('ref'), vmin=-5, vmax=5)
        axes[1,i].set_title('(S-band - C-band) Zh_rain (elev='+ str(theta_radar[i])+'$^o$)')                                                  
        axes[2,i].set_title('(S-band - C-band) Zh_snow (elev='+ str(theta_radar[i])+'$^o$)')    
        axes[3,i].set_title('(S-band - C-band) Zh_grau (elev='+ str(theta_radar[i])+'$^o$)')                                                        
        axes[1,i].grid(); axes[2,i].grid(); axes[3,i].grid() 
        for ii in [120, 240]: 
            [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],ii)
            axes[1,i].plot(lon_radius, lat_radius, 'k', linewidth=0.8)   
            axes[2,i].plot(lon_radius, lat_radius, 'k', linewidth=0.8)   
            axes[3,i].plot(lon_radius, lat_radius, 'k', linewidth=0.8)   
    cbar = plt.colorbar(pcm2, ax=axes[1,2], shrink=1, label='Zh [dBZ]', ticks = np.arange(-5,5.01,1))
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    
    cbar = plt.colorbar(pcm2, ax=axes[2,2], shrink=1, label='Zh [dBZ]', ticks = np.arange(-5,5.01,1))
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    
    cbar = plt.colorbar(pcm2, ax=axes[3,2], shrink=1, label='Zh [dBZ]', ticks = np.arange(-5,5.01,1))
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    
    
    #------------------------------------        
    # Plot Zh_tot transformed to radar grid for the first three elevations 
    #------------------------------------
    for i in [0,1,2]: 
        diff_t = to_dBZ(Zh_t_central[0,i,:,:]) -  to_dBZ(Zh_t_central[bandID,i,:,:])
        pcm0 = axes[4,i].pcolormesh(lon, lat, diff_t, cmap=functions.colormaps('ref'), vmin=-5, vmax=5)
        axes[4,i].set_title('(S-band - C-band) Zh_tot (elev='+ str(theta_radar[i])+'$^o$)')                                                  
        axes[4,i].grid(); axes[4,i].grid(); axes[4,i].grid() 
        for ii in [120, 240]: 
            [lat_radius, lon_radius] = functions.pyplot_rings(radar.latitude['data'][0],radar.longitude['data'][0],ii)
            axes[4,i].plot(lon_radius, lat_radius, 'k', linewidth=0.8)   
    cbar = plt.colorbar(pcm0, ax=axes[4,2], shrink=1, label='Zh [dBZ]', ticks = np.arange(-5,5.01,1))
    cbar.cmap.set_under('white')
    cbar.cmap.set_over('white')
    
    axes[1,3].axis('off')
    axes[2,3].axis('off')
    axes[3,3].axis('off')
    axes[4,3].axis('off')

    plt.show()    
    fig.savefig(save_dir_compare+'/RadarSim/'+'/RadarSim_RADARgrid_CbandvsSband'+EXP+'2018-11-10_'+time+'.png', 
                dpi=300,transparent=False, bbox_inches='tight')        
    return


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def run_and_plot_stats(radar, Zh_r_central, Zh_s_central, Zh_g_central, Zh_t_central, EXP, time ): 

    
    theta_radar = (0.5, 0.9, 1.3, 1.9, 2.3, 3, 3.5, 5, 6.9, 9.1, 11.8, 15.1)                                       

    folders=config_folders.config_folders('yakaira')
    save_dir_compare=folders['save_dir_compare']


    # Filter obsdata w/ phidp and rhohv
    # Configure a gatefilter to filter out copolar correlation coefficient values > 0.9
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_below("RHOHV", 0.9)
    refl_array = np.ma.masked_where(gatefilter.gate_included == False, radar.fields['TH']['data'])        
    
    #for nlev in [0,1,2]:
        # fig = plt.figure(figsize=[15,11])
        # start_index = radar.sweep_start_ray_index['data'][nlev]
        # end_index   = radar.sweep_end_ray_index['data'][nlev]
        # lats        = radar.gate_latitude['data'][start_index:end_index]
        # lons        = radar.gate_longitude['data'][start_index:end_index]
        # plt.pcolormesh(lons, lats, refl_array[start_index:end_index], cmap=fun.colormaps('ref'), vmin=-80, vmax=80)
        # plt.colorbar()
        #fig = plt.figure(figsize=[15,11])
        #plt.pcolormesh(lon, lat, to_dBZ(Zh_t_central[1,nlev,:,:]), cmap=fun.colormaps('ref'), vmin=-80, vmax=80)
        #plt.colorbar()
    
    # for i in [0,1,2]:
    #     fig = plt.figure(figsize=[15,11])
    #     plt.pcolormesh(lon, lat, to_dBZ(Zh_s_central[1,i,:,:]), cmap=fun.colormaps('ref'),vmin=-20, vmax=0)
    
    # for i in [0,1,2]:
    #     fig = plt.figure(figsize=[15,11])
    #     plt.pcolormesh(lon, lat, to_dBZ(Zh_r_central[1,i,:,:]), cmap=fun.colormaps('ref'),vmin=-20, vmax=0)
    
    # for i in [0,1,2]:
    #     fig = plt.figure(figsize=[15,11])
    #     plt.pcolormesh(lon, lat, to_dBZ(Zh_g_central[1,i,:,:]), cmap=fun.colormaps('ref'),vmin=-20, vmax=0)
    
    
    fig, axes = plt.subplots(nrows=1, ncols=4, constrained_layout=True,figsize=[26,7])
    for i, ax in enumerate(axes):
        start_index = radar.sweep_start_ray_index['data'][i]
        end_index   = radar.sweep_end_ray_index['data'][i]
        TH          = refl_array[start_index:end_index]
        #sns.histplot(np.ravel(TH), color='red', ax=axes[i], label='obs')
        #sns.histplot(np.ravel(to_dBZ(Zh_r_central[1,i,:,:]) ), kde='True', color='blue', ax=axes[i], label='sim')  
        sns.kdeplot(np.ravel(TH), color='black', linewidth=1.2, ax=axes[i], label='obs')
        #sns.kdeplot(np.ravel(to_dBZ(Zh_t_central[0,i,:,:]) ),  linewidth=1.2, color='darkblue', linestyle='-', ax=axes[i], label='S-band sim')
        sns.kdeplot(np.ravel(to_dBZ(Zh_t_central[1,i,:,:]) ),  linewidth=1.2, color='darkred', linestyle='-', ax=axes[i], label='C-band sim')
    
        #sns.kdeplot(np.ravel(to_dBZ(Zh_r_central[1,i,:,:]) ), color='royalblue', ax=axes[i], label='sim(rain)')
        #sns.kdeplot(np.ravel(to_dBZ(Zh_s_central[1,i,:,:]) ), color='blue', ax=axes[i], label='sim(snow)')
        #sns.kdeplot(np.ravel(to_dBZ(Zh_g_central[1,i,:,:]) ), color='darkgreen', ax=axes[i], label='sim(graupel)')
    
        ax.set_title(f'C-band Zh_tot (elev='+ str(theta_radar[i])+'$^o$)')
        ax.set_xlabel('Zh_tot [dBZ]')
    ax.legend()
    plt.show() 
    #------------------------------------        
    fig.savefig(save_dir_compare+'/RadarSim/'+'/RadarSim_STATS1'+EXP+'2018-11-10_'+time+'.png', 
                dpi=300,transparent=False, bbox_inches='tight')      
    #------------------------------------        
    
    
    fig, axes = plt.subplots(nrows=3, ncols=4, constrained_layout=True,figsize=[26,28])
    for i in range(4):
        start_index = radar.sweep_start_ray_index['data'][i]
        end_index   = radar.sweep_end_ray_index['data'][i]
        TH          = refl_array[start_index:end_index]
        sns.kdeplot(np.ravel(TH), color='black', linewidth=1.2, ax=axes[0,i], label='obs')
        sns.kdeplot(np.ravel(TH), color='black', linewidth=1.2, ax=axes[1,i], label='obs')
        sns.kdeplot(np.ravel(TH), color='black', linewidth=1.2, ax=axes[2,i], label='obs')
    
        sns.kdeplot(np.ravel(to_dBZ(Zh_r_central[1,i,:,:]) ), color='darkred', ax=axes[0,i], label='sim(rain)')
        sns.kdeplot(np.ravel(to_dBZ(Zh_s_central[1,i,:,:]) ), color='darkblue', ax=axes[1,i], label='sim(snow)')
        sns.kdeplot(np.ravel(to_dBZ(Zh_g_central[1,i,:,:]) ), color='darkgreen', ax=axes[2,i], label='sim(graupel)')
    
        axes[0,i].set_title(f'C-band Zh_rain (elev='+ str(theta_radar[i])+'$^o$)')
        axes[1,i].set_title(f'C-band Zh_snow (elev='+ str(theta_radar[i])+'$^o$)')
        axes[2,i].set_title(f'C-band Zh_grau (elev='+ str(theta_radar[i])+'$^o$)')
    
    axes[2,0].set_xlabel('Zh [dBZ]')
    axes[0,3].legend()
    axes[1,3].legend()
    axes[2,3].legend()
    #------------------------------------  
    plt.show() 
      
    fig.savefig(save_dir_compare+'/RadarSim/'+'/RadarSim_STATS2'+EXP+'2018-11-10_'+time+'.png', 
                dpi=300,transparent=False, bbox_inches='tight')      
    #------------------------------------        
    
    return

