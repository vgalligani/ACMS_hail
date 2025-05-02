#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 11:14:13 2025

@author: vito.galligani
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.colors as colors
from config import config_folders
import seaborn as sns 
# calculate the difference between footprints in the range for stratifom and 
# in the range for convective specifically for sieron,. ti think of conv. vs. stratiform

plt.matplotlib.rc('font', family='serif', size = 16)
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16  

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_maskedDomain(extent, d_cs, Tc):
    
    lon = d_cs['wrf_lon'].data
    lat = d_cs['wrf_lat'].data
    
    # Crop by extent
    mask = (
        (lon >= extent[0]) & (lon <= extent[1]) &
        (lat >= extent[2]) & (lat <= extent[3])
        )
    rows, cols = np.where(mask)
    
    # If any points are found, crop arrays to the bounding box
    if rows.size > 0 and cols.size > 0:
        row_min, row_max = rows.min(), rows.max()
        col_min, col_max = cols.min(), cols.max()
        
        # Crop all arrays accordingly
        Tc_crop = Tc[:, row_min:row_max+1, col_min:col_max+1]
        lat_crop = lat[row_min:row_max+1, col_min:col_max+1]
        lon_crop = lon[row_min:row_max+1, col_min:col_max+1]
        
    else:
        raise ValueError("No data points found within the specified extent.")
    
    
    return Tc_crop

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def cut_extent_liuliu(extent, d_cs, tb_rttov):
    
    # use extent to limit 
    var_cut    = []
    for isnow in range(11):
        rowi1 = []
        for igrau in range(11):
            var_cut1 = get_maskedDomain(extent, d_cs,  tb_rttov[isnow,igrau,:,:,:])
            rowi1.append(var_cut1)
        var_cut.append(rowi1)
    
    var_cut = np.array(var_cut)         
    
    return var_cut

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def cut_extent_liuliu_gaus(extent, d_cs, tb_rttov):
    
    var_cut    = []
    for isnow in range(11):
        rowi1 = []
        for igrau in range(11):
            var_cut1 = get_maskedDomain_gaussian(extent, d_cs,  tb_rttov[isnow,igrau,:,:,:])
            rowi1.append(var_cut1)
        var_cut.append(rowi1)
    
    var_cut = np.array(var_cut)         
    
    return var_cut

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_maskedDomain_gaussian(extent, d_cs, Tc):
    
    lon = d_cs['MHS_lon'].data
    lat = d_cs['MHS_lat'].data
    
    # Crop by extent
    mask = (
        (lon >= extent[0]) & (lon <= extent[1]) &
        (lat >= extent[2]) & (lat <= extent[3])
        )
    rows, cols = np.where(mask)
    
    # If any points are found, crop arrays to the bounding box
    if rows.size > 0 and cols.size > 0:
        row_min, row_max = rows.min(), rows.max()
        col_min, col_max = cols.min(), cols.max()
        
        # Crop all arrays accordingly
        Tc_crop = Tc[:, row_min:row_max+1, col_min:col_max+1]
        lat_crop = lat[row_min:row_max+1, col_min:col_max+1]
        lon_crop = lon[row_min:row_max+1, col_min:col_max+1]
        
    else:
        raise ValueError("No data points found within the specified extent.")
    
    return Tc_crop[:,2:,:]

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_liuliu_sieron_experiments_andcut(experiment, processedFolder):
    
    #extent = [ -64, -50, -40, -20]
    extent = [ -65, -50, -40, -20]

    instrument = 'MHS'
    outfile    = 'output_tb_'+instrument
    
    # rttov clearsky simulations
    d_cs         = xr.open_dataset(processedFolder+'/'+outfile+'rttov_processed_clearsky.nc')
    tb_cs        = d_cs['rttov_cs'].values 
    tb_cs_gaus   = d_cs['rttov_cs_Gaussianantennasigma_'].values
    tb_cs_mean   = d_cs['rttov_cs_footprintmean'].values
    tb_cs_interp = d_cs['rttov_cs_pointInterpNearest'].values
    
    # rttov allsky simulations
    tb_as        = [] 
    tb_as_gaus   = []
    tb_as_mean   = []
    tb_as_interp = []
    dtb_as        = [] 
    dtb_as_gaus   = []
    dtb_as_mean   = []
    dtb_as_interp = []


    for isnow in range(11):
        rowi  = []; rowig  = []; rowim  = []; rowiN  = []
        drowi = []; drowig = []; drowim = []; drowiN = []
        
        
        for  j in range(11):
            expname     = experiment+str(isnow)+'g'+str(j)+'.nc'
            d_liuliu    = xr.open_dataset(processedFolder+'/'+'output_tb_allsky_'+instrument+'sieron_sliu'+str(isnow)+'gliu'+str(j)+expname)    ##output_tb_allsky_MHSsieron_sliu9gliu9
            var         = d_liuliu['rttov_as'].values
            var_gaus    = d_liuliu['rttov_as_Gaussianantennasigma_'].values
            var_mean    = d_liuliu['rttov_as_footprintmean'].values
            var_interp  = d_liuliu['rttov_as_pointInterpNearest'].values
            
            rowi.append(var)
            rowig.append(var_gaus)
            rowim.append(var_mean)
            rowiN.append(var_interp)
            drowi.append(var-tb_cs)
            drowig.append(var_gaus-tb_cs_gaus)
            drowim.append(var_mean-tb_cs_mean)
            drowiN.append(var_interp-tb_cs_interp)
            
        tb_as.append(rowi)
        tb_as_gaus.append(rowig)
        tb_as_mean.append(rowim)
        tb_as_interp.append(rowiN)
        dtb_as.append(drowi)
        dtb_as_gaus.append(drowig)
        dtb_as_mean.append(drowim)
        dtb_as_interp.append(drowiN)
        
    tb_as        = np.array(tb_as)
    tb_as_gaus   = np.array(tb_as_gaus)
    tb_as_mean   = np.array(tb_as_mean)
    tb_as_interp = np.array(tb_as_interp)
    dtb_as        = np.array(dtb_as)
    dtb_as_gaus   = np.array(dtb_as_gaus)
    dtb_as_mean   = np.array(dtb_as_mean)
    dtb_as_interp = np.array(dtb_as_interp)
    
    # WRF data 
    # Id like to get all intg and ints for gaussian and interp. 
    WRFvars        = xr.open_dataset(processedFolder+'/'+'wrfdata_processed.nc')        
    #WRF_intTot_cut = get_maskedDomain2d(extent, d_cs, WRFvars['WRF_intTot'].data) 
    #cloudmask_WRF  = np.ma.masked_less_equal(WRF_intTot_cut, 0.1) 

    rttov_as =dict({'rttov_as':          cut_extent_liuliu(extent, d_cs, tb_as), 
                    'rttov_as_Gaussian': cut_extent_liuliu_gaus(extent, d_cs, tb_as_gaus), 
                    'rttov_as_mean':     cut_extent_liuliu_gaus(extent, d_cs, tb_as_mean),
                    'rttov_as_Interp':   cut_extent_liuliu_gaus(extent, d_cs, tb_as_interp),
                    'delta_rttov_as':    cut_extent_liuliu(extent, d_cs, dtb_as),
                    'delta_Gaussianan':  cut_extent_liuliu_gaus(extent, d_cs, dtb_as_gaus), 
                    'delta_mean':        cut_extent_liuliu_gaus(extent, d_cs, dtb_as_mean), 
                    'delta_Interp':      cut_extent_liuliu_gaus(extent, d_cs, dtb_as_interp) })     
    
    return rttov_as

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def return_data4stats_inputs(rttov): 
        
    # Return values for s9g3 y s9s9 and 89, 157 y 190GHz only 
    
    values     = rttov['rttov_as']    
    var        = values[0,0,0,:,:]
    nlon,nlat = var.shape
    
    values     = rttov['rttov_as_Gaussian']
    var        = values[0,0,0,:,:]
    nlon_obs,nlat_obs = var.shape

    delta_rttov  = np.zeros((11,11,5,nlon,nlat)); delta_rttov[:] = np.nan 
    delta_gaus   = np.zeros((11,11,5,nlon_obs,nlat_obs)); delta_gaus[:]  = np.nan 
    delta_mean   = np.zeros((11,11,5,nlon_obs,nlat_obs)); delta_mean[:]  = np.nan 
    delta_interp = np.zeros((11,11,5,nlon_obs,nlat_obs)); delta_interp[:] = np.nan 
    
    for isnow in range(11):
        for igrau in range(11):
            for ifreq in range(5):
                values     = rttov['rttov_as']
                var        = values[isnow,igrau,ifreq,:,:]
                #masked_var = np.where( (values[isnow,igrau,0,:,:]) <= 250, var, np.nan)
                delta_rttov[isnow,igrau,ifreq,:,:] = var  #masked_var

                values     = rttov['rttov_as_Gaussian']
                var        = values[isnow,igrau,ifreq,:,:]
                #masked_var = np.where( (values[isnow,igrau,0,:,:]) <= 250, var, np.nan)
                delta_gaus[isnow,igrau,ifreq,:,:] = var  #masked_var
                
                values     = rttov['rttov_as_mean']
                var        = values[isnow,igrau,ifreq,:,:]
                #masked_var = np.where( (values[isnow,igrau,0,:,:]) <= 250, var, np.nan)
                delta_mean[isnow,igrau,ifreq,:,:] = var  #masked_var
                
                values     = rttov['rttov_as_Interp']
                var        = values[isnow,igrau,ifreq,:,:]
                #masked_var = np.where( (values[isnow,igrau,0,:,:]) <= 250, var, np.nan)
                delta_interp[isnow,igrau,ifreq,:,:] = var  #masked_var
                
    ifreq = 0
    data_array_89 = {
        'wrf_grid': (delta_rttov[9,3,ifreq,:,:], delta_rttov[9,9,ifreq,:,:]  ),
        'gaussian': (delta_gaus[9,3,ifreq,:,:], delta_gaus[9,9,ifreq,:,:]  ),
        'delta_mean': (delta_mean[9,3,ifreq,:,:], delta_mean[9,9,ifreq,:,:]  ),
        'delta_Interp': (delta_interp[9,3,ifreq,:,:], delta_interp[9,9,ifreq,:,:] ),
        }

    ifreq = 1
    data_array_157 = {
        'wrf_grid': (delta_rttov[9,3,ifreq,:,:], delta_rttov[9,9,ifreq,:,:]  ),
        'gaussian': (delta_gaus[9,3,ifreq,:,:], delta_gaus[9,9,ifreq,:,:]  ),
        'delta_mean': (delta_mean[9,3,ifreq,:,:], delta_mean[9,9,ifreq,:,:]  ),
        'delta_Interp': ( delta_interp[9,3,ifreq,:,:], delta_interp[9,9,ifreq,:,:]  ),
        }

    ifreq = 4
    data_array_190 = {
        'wrf_grid': ( delta_rttov[9,3,ifreq,:,:], delta_rttov[9,9,ifreq,:,:]  ),
        'gaussian': (  delta_gaus[9,3,ifreq,:,:], delta_gaus[9,9,ifreq,:,:]  ),
        'delta_mean': ( delta_mean[9,3,ifreq,:,:], delta_mean[9,9,ifreq,:,:]  ),
        'delta_Interp': (  delta_interp[9,3,ifreq,:,:], delta_interp[9,9,ifreq,:,:]  ),
        }

    
    return data_array_89, data_array_157, data_array_190

#------------------------------------------------------------------------------
def get_diffin_15percentiles(percentile, data_array_0, data_array_1, data_array_4, experiment): 
        
    p89 = []
    p157 = []
    p190 = []

    
    #output deltaBT between footprint model and rttov_as for s9g9 y s9g3. meaning [3] y [4]
    for ii in ([0,1]): 
        p89.append(  np.round( np.nanpercentile( data_array_0['wrf_grid'][ii], percentile)  - np.nanpercentile( data_array_0[experiment][ii], percentile), 2) )
        p157.append(  np.round( np.nanpercentile( data_array_1['wrf_grid'][ii], percentile) - np.nanpercentile( data_array_0[experiment][ii], percentile), 2) )
        p190.append(  np.round( np.nanpercentile( data_array_4['wrf_grid'][ii], percentile) - np.nanpercentile( data_array_0[experiment][ii], percentile), 2) )

        #p25_89.append( np.nanpercentile( data_array_0['wrf_grid'][ii], 25) - np.nanpercentile( data_array_0[experiment][ii], 25) )
        #p25_157.append( np.nanpercentile( data_array_1['wrf_grid'][ii], 25) - np.nanpercentile( data_array_0[experiment][ii], 25) )
        #p25_190.append( np.nanpercentile( data_array_4['wrf_grid'][ii], 25) - np.nanpercentile( data_array_0[experiment][ii], 25) )

    print('-------- difference between wrf-grid and '+experiment)
    print(f'{p89[0]} and {p89[1]}')
    print(f'{p157[0]} and {p157[1]}')
    print(f'{p190[0]} and {p190[1]}')
    
    return  p89, p157, p190


#------------------------------------------------------------------------------
    
processedFolder = '/Users/vito.galligani/Work/Studies/HAIL_20181110/RTTOVout/Processed/'+'WRF-WSM6'       
rttov_sieron = get_liuliu_sieron_experiments_andcut('rttov_processed_allsky_sieron_rsg_s', processedFolder)
data_array_89, data_array_157, data_array_190 = return_data4stats_inputs(rttov_sieron)

get_diffin_15percentiles(15, data_array_89, data_array_157, data_array_190, 'gaussian')


isnow = 9; igrau=9
x_vars = rttov_sieron['rttov_as'][isnow,igrau,0,:,:]
deltavar_wrfgrid = rttov_sieron['delta_rttov_as']

# scatter plot de BT vs. deltaBt
ichan_title = ['89.0', '157.0', '190.311']
chan_indx   = [0,1,4]
x_bins = np.arange(0,310,5)
y_bins = np.arange(-250, -5, 5)
# mean dTB (rttov_as - rttov_cs) below 240K
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=[24,16])    #constrained_layout=True,
for index, i in enumerate(chan_indx):
    axes[0,index].set_title(ichan_title[index]+' GHz')        
    x1   = x_vars.flatten()
    var1 = deltavar_wrfgrid[isnow,igrau,i,:,:].flatten()
    mask1 = ~np.isnan(x1) & ~np.isnan(var1)
    h = axes[0,index].hist2d(x1[mask1], var1[mask1], bins=[x_bins, y_bins], cmap='viridis', norm=colors.LogNorm())
    fig.colorbar(h[3], ax=axes[0, index])   
    axes[0,index].grid(True)    
    axes[0, index].set_yticks(np.arange(-200,0,25))
    if index>0:
        axes[0, index].set_yticklabels([])

#del x_vars, x1, mask1
axes[1,0].scatter([],[],marker='o', label='Gaussian', color='darkblue'); 
axes[1,0].scatter([],[],marker='o', label='Interp.', color='darkred'); 
axes[1,0].legend(loc='upper left')

for index, i in enumerate(chan_indx):
    axes[1,index].set_title(ichan_title[index]+' GHz')        
    x1   = rttov_sieron['rttov_as_Gaussian'][isnow,igrau,0,:,:].flatten()
    var1 = rttov_sieron['delta_Gaussianan'][isnow,igrau,i,:,:].flatten()
    mask1 = ~np.isnan(x1) & ~np.isnan(var1)
    #h = axes[1,index].hist2d(x1[mask1], var1[mask1], bins=[x_bins, y_bins], cmap='viridis', norm=colors.LogNorm())
    #fig.colorbar(h[3], ax=axes[1, index])   
    axes[1,index].scatter(x1[mask1], var1[mask1], marker='o', color='darkblue'); 

    x1   = rttov_sieron['rttov_as_Interp'][isnow,igrau,0,:,:].flatten()
    var1 = rttov_sieron['delta_Interp'][isnow,igrau,i,:,:].flatten()
    mask1 = ~np.isnan(x1) & ~np.isnan(var1)  
    axes[1,index].scatter(x1[mask1], var1[mask1], marker='o', color='darkred'); 
    axes[1,index].grid(True)    

        
    #if index>0:
    #    axes[1, index].set_xticklabels([])
    #    axes[1, index].set_yticklabels([])    
    
for index, i in enumerate(chan_indx):
    axes[0, index].set_xlim([50, 300])
    axes[0, index].set_ylim([-200, 0])
    axes[1, index].set_xlim([50, 300])
    axes[1, index].set_ylim([-200, 0])
axes[0, 0].set_ylabel(r'$\Delta$BT [K] (BT(89clear)<250K) wrf grid')   
axes[1, 0].set_ylabel(r'$\Delta$BT [K] (BT(89clear)<250K) Footprint models')   

axes[0, 0].set_xlabel(r'BT [K] wrf grid')   
axes[1, 0].set_xlabel(r'BT [K] Footprint models')   
    
plt.show()
#fig.savefig(plotpath+'/poster/'+f'scattterplots_for_s{isnow}_g{igrau}.png', dpi=300,transparent=False, bbox_inches='tight')     


# still dont undestund the histogram vs the stacked ?


i = 0 
do_thisSSPs = [2,5,8,9,10]

base_colors = sns.color_palette('Paired')         
base_colors[10] = base_colors[10+1]  
     
fig, axes = plt.subplots(nrows=1, ncols=1, constrained_layout=True,figsize=[10,10])           
for igrau in do_thisSSPs: 
    varliu  = rttov_sieron['rttov_as'][3,igrau,i,:,:]
    varliu  = varliu.flatten() 
    counts, bin_edges = np.histogram(varliu, bins=np.arange(10,300,2),density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    axes.plot(bin_centers, counts, linewidth=2,linestyle='-', color=base_colors[igrau])

    varliu  = rttov_sieron['rttov_as'][9, igrau,i,:,:]
    varliu  = varliu.flatten() 
    counts, bin_edges = np.histogram(varliu, bins=np.arange(10,300,2),density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    axes.plot(bin_centers, counts, linewidth=2,linestyle='--', color=base_colors[igrau])
    
plt.show()
    

# - check delta bts. 
fig, axes = plt.subplots(nrows=1, ncols=1, constrained_layout=True,figsize=[10,10])           
for igrau in do_thisSSPs: 
    varliu  = rttov_sieron['delta_rttov_as'][3,igrau,i,:,:]
    masked_var = np.where((rttov_sieron['rttov_as'][3,igrau,0,:,:]) <= 250, varliu, 0)
    masked_var = masked_var.flatten() 
    counts, bin_edges = np.histogram(masked_var, bins=np.arange(-300,0,10),density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    axes.semilogy(bin_centers, counts, linewidth=2,linestyle='-', color=base_colors[igrau])

    varliu  = rttov_sieron['delta_rttov_as'][9, igrau,i,:,:]
    masked_var = np.where( (rttov_sieron['rttov_as'][9,igrau,0,:,:]) <= 250, varliu, 0)
    masked_var = masked_var.flatten() 
    counts, bin_edges = np.histogram(varliu, bins=np.arange(-300,0,10),density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    axes.semilogy(bin_centers, counts, linewidth=2,linestyle='--', color=base_colors[igrau])
    
plt.show()

for igrau in range(11):
    min_s3 = np.round( np.nanmedian(rttov_sieron['rttov_as'][3,igrau,0,:,:].flatten() ), 2)
    min_s9 = np.round( np.nanmedian(rttov_sieron['rttov_as'][9,igrau,0,:,:].flatten() ), 2)

    print(f'For 89GHz: min(s3): {str(min_s3)} vs. min(s9): {str(min_s9)}')
    
   

