#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 11:52:22 2025

@author: galliganiv
"""

def other(): 
    

    #
    #------------------------------------------------------------------------------------------
    # allskytesting exp1        
    # eqmassWSM6_rsg_s10g2: equal mass PSD consistency with WRF and 
    # default cloud overlap settings as above
    #d_asExp3   = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s3g2.nc')
    #d_asExp7   = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s7g2.nc')
    #d_asExp8   = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s8g2.nc')
    #d_asExp10  = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s10g2.nc')
    #_asExp11  = Dataset(processedFolder+'/'+outfile+'rttov_processed_allsky_rsg_s11g2.nc')
    
    #plt.pcolormesh(d_asExp3['wrf_lon'],d_asExp3['wrf_lat'], d_asExp3['rttov_as']) 
    
    
    #- Clear-sky pixels: For simulations I use cloudmask and for observations? 
    cloudmask_1 = np.ma.masked_greater_equal(WRFvars['MHSfootprintmean_intTot'], 1) # footprint mean
    cloudmask_2 = np.ma.masked_greater_equal(WRFvars['WRF_intTot'], 1)              # point model grid
    cloudmask_3 = np.ma.masked_greater_equal(WRFvars['MHSinterp_intTot'], 1)        # obs interp grid 
    cloudmask_4 = np.ma.masked_greater_equal(WRFvars['MHSGaussian_intTot'], 1)      # Gaussian interp grid 
    cloudobsmask =  np.ma.masked_greater_equal( d_cs['MHs_domain_obs'][3,:,:]-d_cs['MHs_domain_obs'][4,:,:], 0.2)     # very rought estimate: obs interp grid 
    # I also test a more restrictive condition to MHS observed: d_cs['MHs_domain_obs'][1,:,:]-280)<0
    cloudobsmask2 =  np.ma.masked_less_equal( d_cs['MHs_domain_obs'][1,:,:]-280, 0)
    
    cloudmasks = dict({'footprintmean':cloudmask_1, 'WRFgrid':cloudmask_2, 'MHSinterp':cloudmask_3, 
                     'footprintGaussian':cloudmask_4, 'roughOBS':cloudobsmask2})

    #- Histograms for clear sky (using 1) rough cloud masks and 2) all pixels
    filename = 'Stats_clearsky_allpixels_initHist_tests'
    T2P.simple_histograms(d_cs, 'rttov clear-sky', filename, plotpath)
    #
    filename = 'Stats_clearsky_cloudmasks_initHist_tests'
    T2P.simple_histograms_cloudmasks(d_cs, cloudmasks, 'rttov clear-sky (cloudmask)', filename, plotpath)

    #- Plot stats (mean, min, max, std for each simulation)
    EXPsTitles = ['obs', 'footprint mean', 'wrf grid', 'interp. nearest', 'Gaussian antenna']
    colorExp   = ['darkblue','darkred','red','magenta','darkgreen']
    cldmasks   = [cloudobsmask2, cloudmask_1, cloudmask_2, cloudmask_3, cloudmask_4]
    nchan      = 5 
    filename   = 'Stats_clearsky_cloudmasks_initstats'
    T2P.plot1ststats_cs(EXPsTitles, d_cs, cldmasks, nchan, colorExp, plotpath, filename)

    #----------------------------------------------------------------------
    #Look at min/max/std of clear-sky simulations with Cloudy Flags            
    EXPs        = ['MHs_domain_obs', 'rttov_cs_footprintmean', 'rttov_cs', 'rttov_cs_pointInterpNearest', 'rttov_cs_Gaussianantennasigma_']
    titles_axes = ['89.0', '157.0', '183.311$\pm$1', '183.311$\pm$3', '190.311']
    
    dsmean = np.zeros((nchan,nchan))
    dsstd  = np.zeros((nchan,nchan))
    dsmin  = np.zeros((nchan,nchan))
    dsmax  = np.zeros((nchan,nchan))
    ii = 0
    for EXP in EXPs:    
        dsmean[ii], dsstd[ii], dsmin[ii], dsmax[ii] = T2P.calcStats(d_cs[EXP], nchan, cldmasks[ii]) 
        ii=ii+1
    #----------------------------------------------------------------------
    name  = 'ClearSky_stats' 
    Stats = open(processedFolder+'/'+str(name)+'.txt','w')
    print('-----------------------------------', file = Stats)
    print('--- MHS obs w/ strong CloudFlag ---', file = Stats)
    print(' mean: ', np.round(dsmean[0],3), file = Stats)
    print(' std: ' , np.round(dsstd[0],3), file = Stats)
    print(' min: ' , np.round(dsmin[0],3), file = Stats)
    print(' max: ' , np.round(dsmin[0],3), file = Stats)
    print('--- footprintmean, wrfgrid, pointInterp, Gaussian ---', file = Stats)
    print(' mean',  np.round(dsmean[1],3), np.round(dsmean[2],3), np.round(dsmean[3],3), np.round(dsmean[4],3), file = Stats)
    print(' std',   np.round(dsstd[1],3), np.round(dsstd[2],3), np.round(dsstd[3],3), np.round(dsstd[4],3), file = Stats)
    print(' min: ', np.round(dsmin[1],3), np.round(dsmin[2],3), np.round(dsmin[3],3), np.round(dsmin[4],3), file = Stats)
    print(' max: ', np.round(dsmax[1],3), np.round(dsmax[2],3), np.round(dsmax[3],3), np.round(dsmax[4],3), file = Stats)
    print('-----------------------------------', file = Stats)
    Stats.close()
    #----------------------------------------------------------------------
    name  = 'ClearSky_diffstats' 
    Stats = open(processedFolder+'/'+str(name)+'.txt','w')
    print('-----------------------------------', file = Stats)
    print('--- MHS obs w/ strong CloudFlag diff with Gaussian ---', file = Stats)
    for i in range(nchan):
        obs = np.ma.array( d_cs['MHs_domain_obs'][i,:,:].data, mask=cloudobsmask2.mask) 
        sim = np.ma.array( d_cs['rttov_cs_Gaussianantennasigma_'][i,:,:].data, mask=cloudmask_4.mask) 
        print(r'CH ('+titles_axes[i]+'): '+str(np.nanmean(obs-sim)), file = Stats)
        #--------------------------------------------------- 
        # levels = np.arange(-10,10,1)
        # cmap = plt.get_cmap('viridis', len(levels)-1)
        # norm = mpl.colors.BoundaryNorm(boundaries=levels,ncolors=cmap.N)
        # plt.pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], 
        #                (obs-sim),
        #                cmap=cmap, norm=norm); 
        #--------------------------------------------------- 
    print('-----------------------------------', file = Stats)
    Stats.close()
    
    #----------------------------------------------------------------------
    # Std of the mean footprint size clear sky. 
    T2P.plot_MHSstd(d_cs['MHS_lon'], d_cs['MHS_lat'], d_cs['rttov_cs_footprintstd'], 'FootprintSTD_rttovclearsky', plotpath, 'cnrm')
    T2P.simple_footprint_Std_histograms(d_cs, 'clearsky Footprint std', 'FootprintSTD_rttovclearskyHIST', plotpath) 


    #----------------------------------------------------------------------
    #----------------------------------------------------------------------
    # All sky histogram checks look at Plots4Interc.py in harddrive
    T2P.plot_MHS(d_asTest['MHS_lon'], d_asTest['MHS_lat'], d_asTest['rttov_as_Gaussianantennasigma_'] ,'test', plotpath, 'cnrm')
    filename = 'Stats_allsky_initHist_tests'
    T2P.simple_histograms_as(d_asTest, 'rttov all-sky (test)', filename, plotpath)

    #----------------------------------------------------------------------
    # Std of the mean footprint size. 
    T2P.plot_MHSstd(d_asTest['MHS_lon'], d_asTest['MHS_lat'], d_asTest['rttov_as_footprintstd'], 'FootprintSTD_rttovallskytest', plotpath, 'cnrm')
    T2P.simple_footprint_allsky_Std_histograms(d_asTest, 'allsky Footprint std', 'FootprintSTD_rttovallskyHIST', plotpath) 

    #- Define CloudSignal dT
    dT_rttov_footprintmean = np.zeros((d_cs['rttov_cs_footprintmean'].shape)); dT_rttov_footprintmean[:]=np.nan
    dT_rttov_pointInterp   = np.zeros((d_cs['rttov_cs_footprintmean'].shape)); dT_rttov_pointInterp[:]=np.nan
    dT_rttov_WRFgrid       = np.zeros((d_cs['rttov_cs'].shape)); dT_rttov_WRFgrid[:]=np.nan
    dT_rttov_Gaussian      = np.zeros((d_cs['rttov_cs_footprintmean'].shape)); dT_rttov_Gaussian[:]=np.nan
    for i in range(nchan):
        dT_rttov_footprintmean[i,:,:] = d_cs['rttov_cs_footprintmean'][i,:,:] - d_asTest['rttov_as_footprintmean'][i,:,:]
        dT_rttov_pointInterp[i,:,:]   = d_cs['rttov_cs_pointInterpNearest'][i,:,:] - d_asTest['rttov_as_pointInterpNearest'][i,:,:]
        dT_rttov_Gaussian[i,:,:]      = d_cs['rttov_cs_Gaussianantennasigma_'][i,:,:] - d_asTest['rttov_as_Gaussianantennasigma_'][i,:,:]
        dT_rttov_WRFgrid[i,:,:]       = d_cs['rttov_cs'][i,:,:] - d_asTest['rttov_as'][i,:,:]

    #----------------------------------------------------------------------
    EXPs        = ['MHs_domain_obs', 'rttov_as_footprintmean', 'rttov_as', 'rttov_as_pointInterpNearest', 'rttov_as_Gaussianantennasigma_']
    titles_axes = ['89.0', '157.0', '183.311$\pm$1', '183.311$\pm$3', '190.311']
    
    dsmean = np.zeros((nchan,nchan))
    dsstd  = np.zeros((nchan,nchan))
    dsmin  = np.zeros((nchan,nchan))
    dsmax  = np.zeros((nchan,nchan))
    ii = 0
    for EXP in EXPs:    
        dsmean[ii], dsstd[ii], dsmin[ii], dsmax[ii] = T2P.calcStats_nomask(d_asTest[EXP], nchan) 
        ii=ii+1

    dTmean = np.zeros((nchan,nchan))
    dTstd  = np.zeros((nchan,nchan))
    dTmin  = np.zeros((nchan,nchan))
    dTmax  = np.zeros((nchan,nchan))
    dTmean[1], dTstd[1], dTmin[1], dTmax[1] = T2P.calcStats_nomask(dT_rttov_footprintmean, nchan) 
    dTmean[2], dTstd[2], dTmin[2], dTmax[2] = T2P.calcStats_nomask(dT_rttov_WRFgrid , nchan) 
    dTmean[3], dTstd[3], dTmin[3], dTmax[3] = T2P.calcStats_nomask(dT_rttov_pointInterp, nchan) 
    dTmean[4], dTstd[4], dTmin[4], dTmax[4] = T2P.calcStats_nomask(dT_rttov_Gaussian , nchan) 
        
    #----------------------------------------------------------------------
    # aca lo que importa es el max de los diferentes metodos son muy diferentes!
    #----------------------------------------------------------------------
    name  = 'AllSky_stats' 
    Stats = open(processedFolder+'/'+str(name)+'.txt','w')
    print('--------- BT -----------', file = Stats)
    print('--- footprintmean, wrfgrid, pointInterp, Gaussian ---', file = Stats)
    print(' mean',  np.round(dsmean[1],3), np.round(dsmean[2],3), np.round(dsmean[3],3), np.round(dsmean[4],3), file = Stats)
    print(' std',   np.round(dsstd[1],3), np.round(dsstd[2],3), np.round(dsstd[3],3), np.round(dsstd[4],3), file = Stats)
    print(' min: ', np.round(dsmin[1],3), np.round(dsmin[2],3), np.round(dsmin[3],3), np.round(dsmin[4],3), file = Stats)
    print(' max: ', np.round(dsmax[1],3), np.round(dsmax[2],3), np.round(dsmax[3],3), np.round(dsmax[4],3), file = Stats)
    print('-----------------------------------', file = Stats)
    print('-------- BTclear-BTallsky ---------', file = Stats)
    print('--- footprintmean, wrfgrid, pointInterp, Gaussian ---', file = Stats)
    print(' mean',  np.round(dTmean[1],3), np.round(dTmean[2],3), np.round(dTmean[3],3), np.round(dTmean[4],3), file = Stats)
    print(' std',   np.round(dTstd[1],3), np.round(dTstd[2],3), np.round(dTstd[3],3), np.round(dTstd[4],3), file = Stats)
    print(' min: ', np.round(dTmin[1],3), np.round(dTmin[2],3), np.round(dTmin[3],3), np.round(dTmin[4],3), file = Stats)
    print(' max: ', np.round(dTmax[1],3), np.round(dTmax[2],3), np.round(dTmax[3],3), np.round(dTmax[4],3), file = Stats)
    Stats.close()


    x = d_asTest['rttov_as_footprintmean']
    y = d_asTest['rttov_as_Gaussianantennasigma_'] 
    z = d_asTest['rttov_as_pointInterpNearest']
    #----------------------------------------------------------------------
    fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[30,7])  
    for i in range(5):
        axes[i].scatter(x[i,:,:].ravel(), y[i,:,:].ravel(), color='darkblue', alpha=0.8, label='Gaussian antenna')
        axes[i].scatter(x[i,:,:].ravel(), z[i,:,:].ravel(), color='darkred', alpha=0.3, label='Point Interp.')
        axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
        axes[i].grid(True)
    axes[0].set_xlabel('BT rttov footprint mean')
    axes[0].set_ylabel('BT rttov')
    axes[0].legend()
    fig.savefig(plotpath+'/RTTOV/'+'scatterplot_rttovac_test.png', dpi=300,transparent=False)     
    
    x   = d_asTest['rttov_as_footprintmean']
    y   = d_asTest['rttov_as_Gaussianantennasigma_'] 
    z   = d_asTest['rttov_as_pointInterpNearest']
    obs = d_asTest['MHs_domain_obs']
    #----------------------------------------------------------------------
    fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[30,7])  
    for i in range(5):
        axes[i].scatter(obs[i,:,:].ravel(), x[i,:,:].ravel(), color='darkblue', alpha=0.8, label='Footprint mean')
        axes[i].scatter(obs[i,:,:].ravel(), y[i,:,:].ravel(), color='darkgreen', alpha=0.3, label='Gaussian antenna')
        axes[i].scatter(obs[i,:,:].ravel(), z[i,:,:].ravel(), color='darkred', alpha=0.3, label='Point Interp.')
        axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
        axes[i].grid(True)
    axes[0].set_xlabel('BT obs')
    axes[0].set_ylabel('BT rttov')
    axes[0].legend()
    fig.savefig(plotpath+'/RTTOV/'+'scatterplotwOBS_rttovac_test.png', dpi=300,transparent=False)     


    #----------------------------------------------------------------------
    fig, axes = plt.subplots(nrows=4, ncols=5, constrained_layout=True,figsize=[32,14])  
    for i in range(5):
        axes[0,i].scatter(dT_rttov_footprintmean[i,:,:].ravel(), d_asTest['rttov_as_footprintmin'][i,:,:].ravel(), color='darkblue', alpha=0.8)
        #axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
        axes[0,i].grid(True)
    axes[0,0].set_ylabel('min BT in footprint')
    #
    for i in range(5):
        var = d_asTest['rttov_as_footprintmean'][i,:,:].ravel()  
        axes[1,i].scatter(dT_rttov_footprintmean[i,:,:].ravel(), ( var - d_asTest['rttov_as_footprintmin'][i,:,:].ravel()), color='darkred', alpha=0.8)
        #axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
        axes[1,i].grid(True)
    axes[1,0].set_ylabel('Diff in Footprint - min in footprint')        
    #
    for i in range(5):
        axes[2,i].scatter(dT_rttov_footprintmean[i,:,:].ravel(), d_asTest['rttov_as_footprintstd'][i,:,:].ravel(), color='red', alpha=0.8)
        #axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
        axes[2,i].grid(True)
    axes[2,0].set_ylabel('Footprint std')      
    #
    for i in range(5):
        var = d_asTest['rttov_as_Gaussianantennasigma_'][i,:,:].ravel()  
        axes[3,i].scatter(dT_rttov_footprintmean[i,:,:].ravel(), (var - d_asTest['rttov_as_footprintmean'][i,:,:].ravel()), color='blue', alpha=0.8)
        #axes[i].plot([50,300], [50,300], 'r-', linewidth=1.2)   
        axes[3,i].grid(True)
    axes[3,0].set_ylabel('diff Gaussian-Footprint')   
    axes[3,0].set_xlabel('rttov dT(Footsprintmean) (clear-allsky) ')
    #
    fig.savefig(plotpath+'/RTTOV/'+'scatterplot_FootprintDt_VS_footprintmintest.png', dpi=300,transparent=False)     
    #----------------------------------------------------------------------

    # mask = (dT_rttov_footprintmean[0,:,:] == 0) & ~np.isnan(dT_rttov_footprintmean[0,:,:])
    # indices_ = np.argwhere(mask)
    # dT_rttov_footprintmean[0,11,27]
    # CHECK THE ANTENNA PATTERN OF [11,27]
    

    # obs = d_asTest['MHs_domain_obs'][2,:,:].data
    # sim = d_asTest['rttov_as_Gaussianantennasigma_'][2,:,:].data
    # levels = np.arange(-10,10,1)
    # cmap = plt.get_cmap('viridis', len(levels)-1)
    # norm = mpl.colors.BoundaryNorm(boundaries=levels,ncolors=cmap.N)
    # plt.pcolormesh(d_cs['MHS_lon'], d_cs['MHS_lat'], 
    #                (obs-sim),
    #                cmap=cmap, norm=norm);  
    

    # meanrms = np.zeros((3,5)); meanrms[:]=np.nan
    # meansims = np.zeros((3,5)); meansims[:]=np.nan
    # meanobs = [] 
    # rmsobs = [] 
    # fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[30,7])    
    # for i in range(5):
    #     varsim1 = np.ma.array( ds['rttov_cs_footprintmean'].data[i,:,:], mask=cloudmask_1.mask) 
    #     varsim2 = np.ma.array( ds['rttov_cs'].data[i,:,:],               mask=cloudmask_2.mask) 
    #     varsim3 = np.ma.array( ds['rttov_cs_pointInterp'].data[i,:,:],   mask=cloudmask_3.mask)   # update to rttov_cs_pointInterpNearest
    #     varobs  = np.ma.array( ds['MHs_domain_obs'].data[i,:,:],         mask=cloudobsmask.mask) 
    #     sns.histplot( data=varobs.flatten(),  color='darkblue', fill=False, ax=axes[i], element='step', stat="density", bins=np.arange(200,320,1), label='MHS (BT89<270)')
    #     sns.histplot( data=varsim1.flatten(), color='darkred', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (footprint mean)')
    #     sns.histplot( data=varsim2.flatten(), color='red',     fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (point model grid)')
    #     sns.histplot( data=varsim3.flatten(), color='magenta', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (interp obs grid)')
    #     axes[i].set_title(titles_axes[i]+' GHz')
    #     axes[i].set_xlim([200,310])
    #     meansims[0,i] = np.round( np.nanmean( varsim1.flatten() ) ,2) 
    #     meansims[1,i] = np.round( np.nanmean( varsim2.flatten() ) ,2) 
    #     meansims[2,i] = np.round( np.nanmean( varsim3.flatten() ) ,2) 
    #     meanobs.append( np.round( np.nanmean( varobs.flatten() ), 2))

    #     meanrms[0,i] = np.round( np.nanstd( varsim1.flatten() ) ,2) 
    #     meanrms[1,i] = np.round( np.nanstd( varsim2.flatten() ) ,2) 
    #     meanrms[2,i] = np.round( np.nanstd( varsim3.flatten() ) ,2) 
    #     rmsobs.append( np.round( np.nanstd( varobs.flatten() ), 2))
        
    #     string = 'Mean (std) \n'
    #     string = string + 'rttov_cs footprint av.: ' + str( meansims[0,i] ) + 'K ('+ str( meanrms[0,i] ) + ') \n'
    #     string = string + 'rttov_cs (model grid): ' + str( meansims[1,i] ) + 'K ('+ str( meanrms[1,i] ) + ') \n'
    #     string = string + 'rttov_cs (nearest interp): ' + str( meansims[2,i] ) + 'K ('+ str( meanrms[2,i] ) + ') \n'
    #     string = string + 'obs: ' + str( np.nanmean( varobs.flatten() )) + 'K ('+ str( rmsobs[i] ) + ') \n'
        
    #     axes[i].text(.02, .98, string, ha='left', va='top', bbox=dict(facecolor='gray', 
    #                     edgecolor='black', boxstyle='round') ,  transform=axes[i].transAxes)

    # axes[0].set_xlabel('Brightness Temperature (K)')
    # plt.suptitle('Sims w/ CloudFlag (qtot_int < 1)')
    # fig.savefig(plotpath+'/RTTOV/'+'Stats_clearsky_qtotthresholds_init_tests_cloud1.png', dpi=300,transparent=False)    
    
    # #-----------------------------------------------------------------------------
    # #----- SIMPLE PLOTS for stats (cloudmask)        
    # cloudmask_1 = np.ma.masked_greater_equal(dwrf['MHSfootprintmean_intTot'], 0.1) # footprint mean
    # cloudmask_2 = np.ma.masked_greater_equal(dwrf['WRF_intTot'], 0.1) # point model grid
    # cloudmask_3 = np.ma.masked_greater_equal(dwrf['MHS_intTot'], 0.1) # obs interp grid 
    # cloudobsmask =  np.ma.masked_less_equal(ds['MHs_domain_obs'][0,:,:], 280) # obs interp grid 

    # meanrms = np.zeros((3,5)); meanrms[:]=np.nan
    # meansims = np.zeros((3,5)); meansims[:]=np.nan
    # meanobs = [] 
    # rmsobs = [] 
    # fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=True,figsize=[30,7])    
    # for i in range(5):
    #     varsim1 = np.ma.array( ds['rttov_cs_footprintmean'].data[i,:,:], mask=cloudmask_1.mask) 
    #     varsim2 = np.ma.array( ds['rttov_cs'].data[i,:,:],               mask=cloudmask_2.mask) 
    #     varsim3 = np.ma.array( ds['rttov_cs_pointInterp'].data[i,:,:],   mask=cloudmask_3.mask)   # update to rttov_cs_pointInterpNearest
    #     varobs  = np.ma.array( ds['MHs_domain_obs'].data[i,:,:],         mask=cloudobsmask.mask) 
    #     sns.histplot( data=varobs.flatten(),  color='darkblue', fill=False, ax=axes[i], element='step', stat="density", bins=np.arange(200,320,1), label='MHS (BT89<270)')
    #     sns.histplot( data=varsim1.flatten(), color='darkred', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (footprint mean)')
    #     sns.histplot( data=varsim2.flatten(), color='red',     fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (point model grid)')
    #     sns.histplot( data=varsim3.flatten(), color='magenta', fill=False, ax=axes[i], element='step',stat="density", bins=np.arange(200,320,1),label='rttov_cs (interp obs grid)')
    #     axes[i].set_title(titles_axes[i]+' GHz')
    #     axes[i].set_xlim([200,310])
    #     meansims[0,i] = np.round( np.nanmean( varsim1.flatten() ) ,2) 
    #     meansims[1,i] = np.round( np.nanmean( varsim2.flatten() ) ,2) 
    #     meansims[2,i] = np.round( np.nanmean( varsim3.flatten() ) ,2) 
    #     meanobs.append( np.round( np.nanmean( varobs.flatten() ), 2))

    #     meanrms[0,i] = np.round( np.nanstd( varsim1.flatten() ) ,2) 
    #     meanrms[1,i] = np.round( np.nanstd( varsim2.flatten() ) ,2) 
    #     meanrms[2,i] = np.round( np.nanstd( varsim3.flatten() ) ,2) 
    #     rmsobs.append( np.round( np.nanstd( varobs.flatten() ), 2))
        
    #     string = 'Mean (std) \n'
    #     string = string + 'rttov_cs footprint av.: ' + str( meansims[0,i] ) + 'K ('+ str( meanrms[0,i] ) + ') \n'
    #     string = string + 'rttov_cs (model grid): ' + str( meansims[1,i] ) + 'K ('+ str( meanrms[1,i] ) + ') \n'
    #     string = string + 'rttov_cs (nearest interp): ' + str( meansims[2,i] ) + 'K ('+ str( meanrms[2,i] ) + ') \n'
    #     string = string + 'obs: ' + str( np.nanmean( varobs.flatten() )) + 'K ('+ str( rmsobs[i] ) + ') \n'
        
    #     axes[i].text(.02, .98, string, ha='left', va='top', bbox=dict(facecolor='gray', 
    #                     edgecolor='black', boxstyle='round') ,  transform=axes[i].transAxes)

    # axes[0].set_xlabel('Brightness Temperature (K)')
    # plt.suptitle('Sims w/ CloudFlag (qtot_int < 0.1)')
    # fig.savefig(plotpath+'/RTTOV/'+'Stats_clearsky_qtotthresholds_init_tests_cloud01.png', dpi=300,transparent=False) 



    return