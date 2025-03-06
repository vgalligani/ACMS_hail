import numpy as np
from scipy.special import gamma
import wrf
from scipy import integrate
import copy
from scipy.interpolate import interp2d
from pyart.core.transforms import antenna_to_cartesian
import pyart
import math 
import wradlib as wrl

#------------------------------------------------------------------------------  
def adjust_fhc_colorbar_for_pyart(cb):
    
    # HID types:           Species #:  
    # -------------------------------
    # Drizzle                  1    
    # Rain                     2    
    # Ice Crystals             3
    # Aggregates               4
    # Wet Snow                 5
    # Vertical Ice             6
    # Low-Density Graupel      7
    # High-Density Graupel     8
    # Hail                     9
    # Big Drops                10
    
    
    cb.set_ticks(np.arange(2.4, 10, 0.9))
    cb.ax.set_yticklabels(['Rain', 'Ice Crystals', 'Aggregates',
                           'Wet Snow', 'Vertical Ice', 'LD Graupel',
                           'HD Graupel', 'Hail', 'Big Drops'])
    cb.ax.set_ylabel('')
    cb.ax.tick_params(length=0)
    return cb

#------------------------------------------------------------------------------
def get_z_from_radar(radar):
    """Input radar object, return z from radar (km, 2D)"""
    azimuth_1D = radar.azimuth['data']
    elevation_1D = radar.elevation['data']
    srange_1D = radar.range['data']
    sr_2d, az_2d = np.meshgrid(srange_1D, azimuth_1D)
    el_2d = np.meshgrid(srange_1D, elevation_1D)[1]
    xx, yy, zz = antenna_to_cartesian(sr_2d/1000.0, az_2d, el_2d) # Cartesian coordinates in meters from the radar.

    return zz + radar.altitude['data'][0]

#------------------------------------------------------------------------------
def interpolate_sounding_to_radar(snd_T, snd_z, radar):
    """Takes sounding data and interpolates it to every radar gate."""
    radar_z = get_z_from_radar(radar)
    radar_T = None
    shape   = np.shape(radar_z)
    rad_z1d = radar_z.ravel()
    rad_T1d = np.interp(rad_z1d, snd_z, snd_T)
    return np.reshape(rad_T1d, shape), radar_z

#------------------------------------------------------------------------------  
def add_field_to_radar_object(field, radar, field_name='FH', units='unitless', 
                              long_name='Hydrometeor ID', standard_name='Hydrometeor ID'):
    """
    Adds a newly created field to the Py-ART radar object. If reflectivity is a masked array,
    make the new field masked the same as reflectivity.
    """
    if 'TH' in radar.fields.keys():  
        dz_field='TH'
    elif 'DBZH' in radar.fields.keys():
        dz_field='DBZH'   
    elif 'attenuation_corrected_differential_reflectivity' in radar.fields.keys():
        dz_field='attenuation_corrected_differential_reflectivity'   
    fill_value = -32768
    masked_field = np.ma.asanyarray(field)
    masked_field.mask = masked_field == fill_value
    if hasattr(radar.fields[dz_field]['data'], 'mask'):
        setattr(masked_field, 'mask', 
                np.logical_or(masked_field.mask, radar.fields[dz_field]['data'].mask))
        fill_value = radar.fields[dz_field]['_FillValue']
    field_dict = {'data': masked_field,
                  'units': units,
                  'long_name': long_name,
                  'standard_name': standard_name,
                  '_FillValue': fill_value}
    radar.add_field(field_name, field_dict, replace_existing=True)

    return radar

#------------------------------------------------------------------------------
def airdensity(p,t): 
    """ Function to calculate the air density by trivial
    application of the ideal gas law
    """
    Re  = 287.04    
    rho = p / ( Re * t )       
    return rho


#------------------------------------------------------------------------------
def mixr2massconc(mixr,  pres, temp):
    """ Function to calculate the mass cocentration
    from mass mixing ratio assuming a mixture of an 
    atmospheric species and dry air. Output in kg/m3 
    """
    Re  = 287.04    
    rho = airdensity( pres, temp )          
    massconc = rho * mixr               
    return massconc

#------------------------------------------------------------------------------
def get_q_ints6(ncfile):
    
    Re       = 6.3781e6
    
    # ncfile1: WSM6
    temp     = wrf.g_temp.get_tk(ncfile)
    pressure = wrf.g_pressure.get_pressure(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL
    z_level  = Re*geopo_p/(Re-geopo_p)


    qr = mixr2massconc( np.squeeze(ncfile.variables["QRAIN"][0,:,:,:]  ), pressure, temp )        
    qs = mixr2massconc( np.squeeze(ncfile.variables["QSNOW"][0,:,:,:]  ), pressure, temp )        
    qi = mixr2massconc( np.squeeze(ncfile.variables["QICE"][0,:,:,:]   ), pressure, temp )        
    qc = mixr2massconc( np.squeeze(ncfile.variables["QCLOUD"][0,:,:,:] ), pressure, temp )       
    qg = mixr2massconc( np.squeeze(ncfile.variables["QGRAUP"][0,:,:,:] ), pressure, temp ) 

    qr_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qr)) , z_level, axis=0)
    qs_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qs)) , z_level, axis=0)
    qg_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qg)) , z_level, axis=0)
    qi_int = integrate.trapz(np.ma.array(qi, mask=np.isnan(qi)) , z_level, axis=0)
    qc_int = integrate.trapz(np.ma.array(qc, mask=np.isnan(qc)) , z_level, axis=0)
    
    qr_int[qr_int<0.0001] = np.nan
    qs_int[qs_int<0.0001] = np.nan
    qg_int[qg_int<0.0001] = np.nan
    qi_int[qi_int<0.0001] = np.nan
    qc_int[qc_int<0.0001] = np.nan
    
    
    return [qr, qs, qi, qc, qg, qr_int, qs_int, qi_int, qc_int, qg_int]

#------------------------------------------------------------------------------
def get_q_ints3(ncfile):
    
    Re       = 6.3781e6
    
    # ncfile1: WSM6
    temp     = wrf.g_temp.get_tk(ncfile)
    pressure = wrf.g_pressure.get_pressure(ncfile)
    geopo_p  = wrf.g_geoht.get_height(ncfile) # geopotential height as Mean Sea Level (MSL)
    z_level  = Re*geopo_p/(Re-geopo_p)


    qr = mixr2massconc( np.squeeze(ncfile.variables["QRAIN"][0,:,:,:]  ), pressure, temp )        
    qi = mixr2massconc( np.squeeze(ncfile.variables["QICE"][0,:,:,:]   ), pressure, temp )        
    qc = mixr2massconc( np.squeeze(ncfile.variables["QCLOUD"][0,:,:,:] ), pressure, temp )       

    qir = mixr2massconc( np.squeeze(ncfile.variables["QIR"][0,:,:,:] ), pressure, temp )       
    qib = mixr2massconc( np.squeeze(ncfile.variables["QIB"][0,:,:,:] ), pressure, temp )       
    qil = mixr2massconc( np.squeeze(ncfile.variables["QIL"][0,:,:,:] ), pressure, temp )       

    qr_int = integrate.trapz(np.ma.array(qr, mask=np.isnan(qr)) , z_level, axis=0)
    qi_int = integrate.trapz(np.ma.array(qi, mask=np.isnan(qi)) , z_level, axis=0)
    qc_int = integrate.trapz(np.ma.array(qc, mask=np.isnan(qc)) , z_level, axis=0)
    
    qr_int[qr_int<0.0001] = np.nan
    qi_int[qi_int<0.0001] = np.nan
    qc_int[qc_int<0.0001] = np.nan
    
    return [qr, qi, qc, qr_int, qi_int, qc_int, qir, qib, qil]    

#------------------------------------------------------------------------------
def find_0( arr, target ):
    
    # Find the indices and values closest to 273 in each column
    diff =  np.abs(arr.data - target)
    closest_indices = np.nanargmin(diff, axis=0)
    ### closest_values = arr.data[closest_indices, np.arange(arr.shape[1])]
    #closest_values = height_arr.data[closest_indices, np.arange(arr.shape[1])]

    return closest_indices

#------------------------------------------------------------------------------
def find_00( arr, target, arrZ):
    
    x_dim = arr.data.shape[1]
    y_dim = arr.data.shape[2]
    
    output_array_T = np.empty((x_dim, y_dim))
    output_array_Z = np.empty((x_dim, y_dim))
    
    # Iterate over each (x, y) position to find the closest value along the z-axis
    for i in range(x_dim):
        for j in range(y_dim):
            # Find the index in z-dimension closest to the target value
            closest_index = np.argmin(np.abs(arr.data[:, i, j] - target))
            # Store the value at this index
            output_array_T[i, j] = arr.data[closest_index, i, j]
            output_array_Z[i, j] = arrZ.data[closest_index, i, j]

    return output_array_T, output_array_Z

#------------------------------------------------------------------------------
# Helper function to convert HHMM to minutes since midnight
def hhmm_to_minutes(hhmm):
    hours = int(hhmm[:2])
    minutes = int(hhmm[2:])
    return hours * 60 + minutes

#------------------------------------------------------------------------------
# Helper function to convert minutes since midnight back to HHMM format
def minutes_to_hhmm(minutes):
    hours = minutes // 60
    mins = minutes % 60
    return f"{hours:02d}{mins:02d}"


#------------------------------------------------------------------------------
def get_colmax(radar, field, gatefilter):
    
    # Determine the lowest sweep (used for metadata and such)
    minimum_sweep = np.min(radar.sweep_number["data"])
    
    # loop over all measured sweeps
    for sweep in sorted(radar.sweep_number["data"]):
        # get start and stop index numbers
        sweep_slice = radar.get_slice(sweep)

        # grab radar data
        z = radar.get_field(sweep, field)
        z_dtype = z.dtype

        # Use gatefilter
        if gatefilter is not None:
           mask_sweep = gatefilter.gate_excluded[sweep_slice, :]
           z = np.ma.masked_array(z, mask_sweep)

        # extract lat lons
        lon = radar.gate_longitude["data"][sweep_slice, :]
        lat = radar.gate_latitude["data"][sweep_slice, :]   
        
        # get the range and time
        ranges = radar.range["data"]
        time = radar.time["data"]

        # get azimuth
        az = radar.azimuth["data"][sweep_slice]
        # get order of azimuths
        az_ids = np.argsort(az)

        # reorder azs so they are in order
        az = az[az_ids]
        z = z[az_ids]
        lon = lon[az_ids]
        lat = lat[az_ids]
        time = time[az_ids]       
        
        # if the first sweep, store re-ordered lons/lats
        if sweep == minimum_sweep:
            azimuth_final = az
            time_final = time
            lon_0 = copy.deepcopy(lon)
            lon_0[-1, :] = lon_0[0, :]
            lat_0 = copy.deepcopy(lat)
            lat_0[-1, :] = lat_0[0, :]
    
        else:
            # Configure the intperpolator
            z_interpolator = interp2d(ranges, az, z, kind="linear")

            # Apply the interpolation
            z = z_interpolator(ranges, azimuth_final)
            
        # if first sweep, create new dim, otherwise concat them up
        if sweep == minimum_sweep:
            z_stack = copy.deepcopy(z[np.newaxis, :, :])
        else:
            z_stack = np.concatenate([z_stack, z[np.newaxis, :, :]])

    
    # now that the stack is made, take max across vertical
    out_compz = np.nanmax(z_stack.data, axis=0) #.astype(z_dtype)
    

    
    
    return lon_0,lat_0,out_compz


#------------------------------------------------------------------------------
# Helper function to convert HHMM to minutes since midnight
def hhmm_to_minutes(hhmm):
    hours = int(hhmm[:2])
    minutes = int(hhmm[2:])
    return hours * 60 + minutes

#------------------------------------------------------------------------------
def maxHailSize(nit, lam, mu, pres, temp): 

    # !--------------------------------------------------------------------------
    # ! Computes the maximum hail size by estimating the maximum size that is
    # ! physically observable (and not just a numerical artifact of the complete
    # ! gamma size distribution).
    # !
    # ! Follows the method described in Milbrandt and Yau (2006a).
    # !
    # !--------------------------------------------------------------------------
    # ! Arguments:
    #  real, intent(in) :: rho        ! air density   [kg m-3]
    #  real, intent(in) :: nit        ! total num and total number mixing ratio
    #  real, intent(in) :: rhofaci    ! air density correction factor for ice fall speed
    #  real, intent(in) :: lam,mu     ! PSD slope and shape parameters
    
    # ! Local variables:
    #  real, parameter  :: dD       = 1.e-3    ! diameter bin width [m]
    #  real, parameter  :: Dmax_psd = 150.e-3  ! maximum diameter in PSD to compute integral  [m]
    #  real, parameter  :: Ncrit    = 5.e-4    ! threshold physically observable number concentration [# m-3]
    #  real, parameter  :: Rcrit    = 1.e-3    ! threshold physically observable number flux          [# m-2 s-1]
    #  real, parameter  :: ch       = 206.89   ! coefficient in V-D fall speed relation for hail (from MY2006a)
    #  real, parameter  :: dh       = 0.6384   ! exponent in V-D fall speed relation for hail (from MY2006a)
    #  double precision :: n0                  ! shape parameter in gamma distribution
    #  real             :: Di                  ! diameter  [m]
    #  real             :: N_tot               ! total number concentration  [# m-3]
    #  real             :: N_tail              ! number conc. from Di to infinity; i.e. trial for Nh*{D*} in MY2006a [# m-3]
    #  real             :: R_tail              ! number flux of large hail; i.e. trial for Rh*{D*} (corrected from MY2006a [# m-2 s-1]
    #  real             :: V_h                 ! fall speed of hail of size D     [m s-1]
    #  integer          :: nd                  ! maximum number of size bins for integral
    #  integer          :: i                   ! index for integration
    
    # !-----------------------------------------------------------------------
    # note that rhofaci = (rhrosui*inv_rho)**0.54
    rd     = 287.15
    rhosui = 60000/(rd*253.15)
    # calculate some time-varying atmospheric variables
    rho     = pres/(rd*temp)   # dry air density 
    inv_rho = 1/rho            # inverse
    rhofaci =  (rhosui*inv_rho)**0.54
    
    # !-----------------------------------------------------------------------    
    Dmax_psd = 150e-3  # maximum diameter in PSD to compute integral  [m]
    dD       = 1e-3    # diameter bin width [m]
    ch       = 206.89  # coefficient in V-D fall speed relation for hail (from MY2006a)
    Dh       = 0.6384  # exponent in V-D fall speed relation for hail (from MY2006a)
    Rcrit    = 1e-3    # threshold physically observable number flux          [# m-2 s-1]
                                                                                

    maxHailSize = 0
    nd  = int(Dmax_psd/dD)
    n0  = nit*lam**(mu+1.)/gamma(mu+1.)

    #-- method 1, based on Rh*crit:
    R_tail   = 0
    for i in range(nd, 0, -1):  
        Di  = i*dD
        V_h = rhofaci*(ch*Di**Dh)
        R_tail = R_tail + V_h*n0*Di**mu*np.exp(-lam*Di)*dD
        if (R_tail>Rcrit):
            maxHailSize = Di            
            return maxHailSize
        
 
    
#------------------------------------------------------------------------------
def get_sys_phase_simple(radar):
    """
    -------------------------------------------------------------------------------
    Estimate of the system phase.
    Uso los 60 gates en adelante de cada radial en RMA1 (probar con otros valores!)
    -------------------------------------------------------------------------------
    """
    Fromgates   = 60
    start_index = radar.sweep_start_ray_index['data'][0]
    end_index   = radar.sweep_end_ray_index['data'][0]

    phases_nlev = []

    for nlev in range(radar.nsweeps-1):
        start_index = radar.sweep_start_ray_index['data'][nlev]
        end_index   = radar.sweep_end_ray_index['data'][nlev]
        lats  = radar.gate_latitude['data'][start_index:end_index]
        lons  = radar.gate_longitude['data'][start_index:end_index]
        if 'TH' in radar.fields.keys():
         TH    = radar.fields['TH']['data'][start_index:end_index]
        elif 'DBZH' in radar.fields.keys():
          TH =radar.fields['DBZH']['data'][start_index:end_index]
        if 'TV' in radar.fields.keys():
          TV    = radar.fields['TV']['data'][start_index:end_index]
        elif 'DBZV' in radar.fields.keys():
          TV    = radar.fields['DBZV']['data'][start_index:end_index]
        RHOHV = radar.fields['RHOHV']['data'][start_index:end_index]
        PHIDP = np.array(radar.fields['PHIDP']['data'][start_index:end_index])
        PHIDP[np.where(PHIDP==radar.fields['PHIDP']['data'].fill_value)] = np.nan

        # Aca usamos algunos filtros:
        PHIDP = np.where( (RHOHV>0.8) & (TH>30), PHIDP, np.nan)

        # Por cada radial encontrar first non nan value:
        phases = []
        for radial in range(radar.sweep_end_ray_index['data'][0]):
            if firstNonNan(PHIDP[radial,Fromgates:]):
                phases.append(firstNonNan(PHIDP[radial,Fromgates:]))

        phases_nlev.append(np.median(phases))
    phases_out = np.nanmedian(phases_nlev)

    return phases_out
#%%
#------------------------------------------------------------------------------
def correct_phidp(radar, sys_phase, diferencia):
    """
    -------------------------------------------------------------------------------
    Correct phidp
    Ignoro los 1eros 60-1 de cada radial en RMA1 (probar con otros valores y ver efecto!)
    -------------------------------------------------------------------------------
    """
    Fromgates = 40
    rho_data  = radar.fields['RHOHV']['data']
    if 'TH' in radar.fields.keys():
      zh = radar.fields['TH']['data']
    elif 'DBZH' in radar.fields.keys():
      zh = radar.fields['DBZH']['data']
    phi  = radar.fields['PHIDP']['data']

    phiphi = phi.copy()
    rho = rho_data.copy()
    ni = phi.shape[0]
    nj = phi.shape[1]
    for i in range(ni):
        rho_h = rho[i,:]
        zh_h = zh[i,:]
        for j in range(nj):
            if (rho_h[j]<0.8) or (zh_h[j]<30):
                phiphi[i,j]  = np.nan
                rho[i,j]     = np.nan

    # No tomo en cuenta los primeros 60 gates de cada radial
    phiphi[:,0:Fromgates]  = np.nan
    rho[:,0:Fromgates]    = np.nan

    # Despeckle
    dphi = despeckle_phidp(phiphi, rho, zh)

    # Unfold
    uphi_i = unfold_phidp(dphi, rho, diferencia)

    # Accumu ------------------------------------------------------------------------------------->>>>>>>>>>
    uphi_accum = []
    for i in range(ni):
        phi_h = uphi_i[i,:]
        for j in range(1,nj-1,1):
            if phi_h[j] <= np.nanmax(np.fmax.accumulate(phi_h[0:j])):
              	uphi_i[i,j] = uphi_i[i,j-1]

    # Reemplazo nan por sys_phase para que cuando reste esos puntos queden en cero <<<<< ojo aca!
    uphi = uphi_i.copy()
    uphi = np.where(np.isnan(uphi), sys_phase, uphi)
    phi_cor = subtract_sys_phase(uphi, sys_phase)
    phi_cor[phi_cor < 0] = np.nan                                               # antes <= ?
    phi_cor[np.isnan(phi_cor)] = 0                                              # agregado para RMA1?

    # Smoothing final:
    for i in range(ni):
        phi_cor[i,:] = pyart.correct.phase_proc.smooth_and_trim(phi_cor[i,:], window_len=20,
                                            window='flat')

    return dphi, uphi_i, phi_cor
#%%

#------------------------------------------------------------------------------
def firstNonNan(listfloats):

  for item in listfloats:
    if math.isnan(item) == False:

        return item

#------------------------------------------------------------------------------------
def despeckle_phidp(phi, rho, zh):
    '''
    -------------------------------------------------------------------------------
    Elimina pixeles aislados de PhiDP
    -------------------------------------------------------------------------------
    '''
    # Unmask data and despeckle
    dphi = phi.copy()

    # Descartamos pixeles donde RHO es menor que un umbral (e.g., 0.7) o no está definido (e.g., NaN)
    dphi[np.isnan(rho)] = np.nan

    # Calculamos la textura de RHO (rhot) y descartamos todos los pixeles de PHIDP por encima
    # de un umbral de rhot (e.g., 0.25)
    rhot = wrl.dp.texture(rho)
    rhot_thr = 0.25
    dphi[rhot > rhot_thr] = np.nan

    # Eliminamos pixeles aislados rodeados de NaNs
    # https://docs.wradlib.org/en/stable/generated/wradlib.dp.linear_despeckle.html
    #dphi = wrl.dp.linear_despeckle(dphi, ndespeckle=5, copy=False)
    #https://docs.wradlib.org/en/latest/migration.html
    dphi = wrl.util.despeckle(dphi, n=5, copy=False)

    ni = phi.shape[0]
    nj = phi.shape[1]

    return dphi
#%%
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
def unfold_phidp(phi, rho, diferencia):
    '''
    -------------------------------------------------------------------------------
    Unfolding
    -------------------------------------------------------------------------------
    '''

    # Dimensión del PPI (elevaciones, azimuth, bins)
    nb = phi.shape[1]
    nr = phi.shape[0]

    phi_cor = np.zeros((nr, nb)) #Asigno cero a la nueva variable phidp corregida
    v1 = np.zeros(nb)  #Vector v1

    # diferencia = 200 # Valor que toma la diferencia entre uno y otro pixel dentro de un mismo azimuth
    for m in range(0, nr):
        v1 = phi[m,:]
        v2 = np.zeros(nb)
        for l in range(0, nb):
            a = v2[l-1] - v1[l]
            if np.isnan(a):
                v2[l] = v2[l-1]
            elif a > diferencia:  # np.abs(a) ?
                v2[l] = v1[l] + 360 # para -180to180 1ue le sumo 360, aca v1[l]-360
                if v2[l-1] - v2[l] > 100:  # Esto es por el doble folding cuando es mayor a 700
                    v2[l] = v1[l] + v2[l-1]
            else:
                v2[l] = v1[l]
        phi_cor[m,:] = v2

    return phi_cor
#%%
#------------------------------------------------------------------------------------
def subtract_sys_phase(phi, sys_phase):

    nb = phi.shape[1] #GATE
    nr = phi.shape[0] #AZYMUTH
    phi_final = np.copy(phi) * 0
    phi_err=np.ones((nr, nb)) * np.nan

    try:
        phi_final = phi-sys_phase
    except:
        phi_final = phi_err

    return phi_final


