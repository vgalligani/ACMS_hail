# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns 
import matplotlib.ticker as ticker
import os
from datetime import datetime
from config import config_folders

plt.matplotlib.rc('font', family='serif', size = 16)
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16  

#------------------------------------------------
#------------------------------------------------
# define iwc [kg/m3]
def get_lwc(i_iwc):
    
    wc_slope = 0.01 # for version 13.0 which is the one I used 
    wc_offset = 301 # for version 13.0 which is the one I used 
 
    return 10**(wc_slope * (i_iwc - wc_offset))  #[g/3]
#    return 1e-3 * 10**(wc_slope * (i_iwc - wc_offset))  #[kg/3]

#------------------------------------------------
#------------------------------------------------
# define temperature
def get_temp(i_temp):
    # for snow or grau (n_type) =  177.0    
    temp_offset = 177
    return temp_offset + i_temp

#------------------------------------------------
#------------------------------------------------
def readf_txt(file):
    
    data = []
    with open(file) as f:
        for line in f:
            # Skip empty lines if any
            if line.strip():
                values = map(float, line.strip().split())
                data.extend(values)

    # size [80200] == 2*401*100 (n_dia=100)            
    array = np.array(data).reshape(2, 401, 100)       

    return array


#------------------------------------------------
#------------------------------------------------
def read_bulk_txt(file):
    
    data = []
    with open(file) as f:
        for line in f:
            # Skip empty lines if any
            if line.strip():
                values = map(float, line.strip().split())
                data.extend(values)

    # size [802] == 2*401            
    array = np.array(data).reshape(401)       

    return array


#------------------------------------------------
#------------------------------------------------
def readf_diam_txt(file):
    
    data = []
    with open(file) as f:
        for line in f:
            # Skip empty lines if any
            if line.strip():
                values = map(float, line.strip().split())
                data.extend(values)

    return data

#------------------------------------------------
#------------------------------------------------
def plot_nd(snow_diam, grau_diams, WSM6preNorm_snow, WSM6Norm_snow, WSM6preNorm_grau, WSM6Norm_grau, experiment):

    temperature = 263
    
    ssps_liu = np.arange(0,11,1)
    iwc     = [1,201,301]   # these indices correspond to 0.001, 0.1 
    iwc_val = [0.001, 0.1, 1]

    # Get the same colors as other plots:
    base_colors = sns.color_palette('Paired')         
        
    #-- plotn(d) norm and not-normed for all liu in one figure
    fig, axes = plt.subplots(nrows=3, ncols=2, constrained_layout=True,figsize=[16,8])        
    for issp in ssps_liu: 
        axes[0,1].plot([],[],color=base_colors[issp],label=f'Liu: {issp}') 
    axes[0,1].legend(ncol=2)
            
    ii = 0 
    for iiwc in iwc: 
        for issp in ssps_liu:             

            axes[ii,0].loglog(snow_diam[issp,:], WSM6preNorm_snow[issp,0,iiwc,:], linewidth=1.2, color=base_colors[issp]) 
            axes[ii,0].loglog(snow_diam[issp,:], WSM6Norm_snow[issp,0,iiwc,:], marker='o',  linestyle='none', markersize=3, color=base_colors[issp]) 

            axes[ii,1].loglog(grau_diam[issp,:], WSM6preNorm_grau[issp,0,iiwc,:], linewidth=1.2, color=base_colors[issp]) 
            axes[ii,1].loglog(grau_diam[issp,:], WSM6Norm_grau[issp,0,iiwc,:], linestyle='none', linewidth=1.2, markersize=3, marker='o', color=base_colors[issp]) 

            #axes[ii,0].set_ylim([1e-10,1e7])
            #axes[ii,1].set_ylim([1e-10,1e11])
            #axes[ii,0].set_ylim([1e-10,1e11])
            axes[ii,0].grid(True)
            axes[ii,1].grid(True)
            
            axes[ii,0].set_title(f'snow iwc: {iwc_val[ii]} g/m3')
            axes[ii,1].set_title(f'grau iwc: {iwc_val[ii]} g/m3')
            
        ii=ii+1
    
    axes[2,0].set_ylabel('nd(D) [1/m4]')
    axes[2,0].set_xlabel('D [m]')
    plt.suptitle(f'n(D) for {experiment} experiment (Temp: {temperature}K)', fontweight='bold', fontsize=14)
    plt.show()
    fig.savefig(plotpath+'/RTTOV/nD_all'+experiment+'.png', dpi=300,transparent=False)


    #-- plot the difference between n(d) norm and not-normed
    fig, axes = plt.subplots(nrows=3, ncols=2, constrained_layout=True,figsize=[16,8])        
    for issp in ssps_liu: 
        axes[0,1].plot([],[],color=base_colors[issp],label=f'Liu: {issp}') 
    axes[0,1].legend(ncol=2)
            
    ii = 0 
    for iiwc in iwc: 
        for issp in ssps_liu:             

            axes[ii,0].loglog(snow_diam[issp,:], ignore_diffzero(WSM6Norm_snow[issp,0,iiwc,:]-WSM6preNorm_snow[issp,0,iiwc,:]), marker='o', markersize=3, linewidth=1.2, color=base_colors[issp]) 
            axes[ii,1].loglog(grau_diam[issp,:], ignore_diffzero(WSM6Norm_grau[issp,0,iiwc,:]-WSM6preNorm_grau[issp,0,iiwc,:]), marker='o', markersize=3, linewidth=1.2, color=base_colors[issp]) 

            #axes[ii,1].set_ylim([1e-10,1e11])
            #axes[ii,0].set_ylim([1e-10,1e11])
            axes[ii,0].grid(True)
            axes[ii,1].grid(True)
            
            axes[ii,0].set_title(f'snow iwc: {iwc_val[ii]} g/m3')
            axes[ii,1].set_title(f'grau iwc: {iwc_val[ii]} g/m3')
            
            #axes[ii,0].yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
            #axes[ii,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
            #axes[ii,1].yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
            #axes[ii,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
            
            axes[ii,0].set_xlim([1e-4,1e-2])
            axes[ii,1].set_xlim([1e-4,1e-2])
            
        ii=ii+1
    
    axes[2,0].set_ylabel('$\Delta$ nd(D) due to normalization [1/m4]')
    axes[2,0].set_xlabel('D [m]')
    plt.suptitle(f'$\Delta$n(D) due to Norm. for {experiment} experiment (Temp: {temperature}K)', fontweight='bold', fontsize=14)
    plt.show()
    fig.savefig(plotpath+'/RTTOV/nD_all_deltaNorm'+experiment+'.png', dpi=300,transparent=False)
        
    #-- plotn(d) norm and not-normed one figure per liu ssp: 
    ii = 0 
    for issp in ssps_liu:    
        
        fig, axes = plt.subplots(nrows=3, ncols=1, constrained_layout=True,figsize=[12,10])
        
        axes[0].plot([],[], linestyle='-', linewidth=1.2, color='k', label='snow') 
        axes[0].plot([],[], linestyle='--', linewidth=1.2, color='k', label='grau') 
        axes[0].plot([],[], linestyle='none', marker='x', markersize=5, color='k', label='snow (renorm.)') 
        axes[0].plot([],[], linestyle='none', marker='o', markersize=5, color='k', label='grau (renorm.)') 
        axes[0].legend(ncol=2)
        
        ii=0
        for iiwc in iwc: 
         
            axes[ii].loglog(snow_diam[issp,:], WSM6preNorm_snow[issp,0,iiwc,:], linestyle='-', linewidth=1.2, color=base_colors[issp]) 
            axes[ii].loglog(snow_diam[issp,:], WSM6Norm_snow[issp,0,iiwc,:], marker='x',  linestyle='none', markersize=5, color=base_colors[issp]) 

            axes[ii].loglog(grau_diam[issp,:], WSM6preNorm_grau[issp,0,iiwc,:], linestyle='--', linewidth=1.2, color=base_colors[issp]) 
            axes[ii].loglog(grau_diam[issp,:], WSM6Norm_grau[issp,0,iiwc,:], linestyle='none', linewidth=1.2, markersize=5, marker='o', color=base_colors[issp]) 

            #axes[ii].set_ylim([1e-5,1e11])
            axes[ii].grid(True)
            axes[ii].grid(True)
            axes[ii].set_title(f'Grau and Snow (liu:{issp}), iwc: {iwc_val[ii]}')
    
            ii=ii+1
    
        axes[2].set_ylabel('nd(D) [1/m4]')
        axes[2].set_xlabel('D [m]')
        plt.suptitle(f'n(D) for {experiment} experiment (Temp: {temperature}K) and Liu: {issp}', fontweight='bold', fontsize=14)
        plt.show()
        #fig.savefig(plotpath+'/RTTOV/nD_perSSP'+experiment+f'liu{issp}.png', dpi=300,transparent=False)
        
    return
#------------------------------------------------
#------------------------------------------------
def ignore_diffzero(var):
    return np.where(var <= 0, np.nan, var)
    

#------------------------------------------------
#------------------------------------------------
def plot_nd_compareEqMass(WSM6path, eqMassWSM6path, sieronWSM6path):

    snow_diam, grau_diam, WSM6preNorm_snow, WSM6Norm_snow, WSM6preNorm_grau, WSM6Norm_grau = get_experiment(WSM6path)
    snow_diam1, grau_diam1, WSM6preNorm_snow1, WSM6Norm_snow1, WSM6preNorm_grau1, WSM6Norm_grau1 = get_experiment(eqMassWSM6path)
    snow_diam2, grau_diam2, WSM6preNorm_snow2, WSM6Norm_snow2, WSM6preNorm_grau2, WSM6Norm_grau2 = get_experiment(sieronWSM6path)
    
    temperature = 263    
    ssps_liu = np.arange(0,11,1)
    iwc     = [1,201,301]   # these indices correspond to 0.001, 0.1 
    iwc_val = [0.001, 0.1, 1]

    # Get the same colors as other plots:
    base_colors = sns.color_palette('Paired')         
        
    #-- plotn(d) norm and not-normed one figure per liu ssp: 
    ii = 0 
    for issp in ssps_liu:    
        
        fig, axes = plt.subplots(nrows=3, ncols=2, constrained_layout=True,figsize=[12,10])
        
        axes[2,0].plot([],[], linestyle='-', linewidth=1.2, color='k', label='WSM6') 
        axes[2,0].plot([],[], linestyle='--', linewidth=1.2, color='k', label='eqMassWSM6') 
        axes[2,0].plot([],[], linestyle='none', marker='o', markersize=5, color='k', label='WSM6 (renorm.)') 
        axes[2,0].plot([],[], linestyle='none', marker='x', markersize=5, color='k', label='eqMassWSM6 (renorm.)') 
        axes[2,0].legend(ncol=2, loc='lower left')
        
        ii=0
        for iiwc in iwc: 
         
            axes[ii,0].loglog(snow_diam[issp,:], WSM6preNorm_snow[issp,0,iiwc,:], linestyle='-', linewidth=1.2, color=base_colors[issp]) 
            axes[ii,0].loglog(snow_diam[issp,:], WSM6Norm_snow[issp,0,iiwc,:], marker='o',  linestyle='none', markersize=5, color=base_colors[issp]) 

            axes[ii,0].loglog(snow_diam1[issp,:], WSM6preNorm_snow1[issp,0,iiwc,:], linestyle='--', linewidth=1.2, color='k') 
            axes[ii,0].loglog(snow_diam1[issp,:], WSM6Norm_snow1[issp,0,iiwc,:], marker='x',  linestyle='none', markersize=5, color='k') 


            axes[ii,1].loglog(grau_diam[issp,:], WSM6preNorm_grau[issp,0,iiwc,:], linestyle='-', linewidth=1.2, color=base_colors[issp]) 
            axes[ii,1].loglog(grau_diam[issp,:], WSM6Norm_grau[issp,0,iiwc,:], linestyle='none', linewidth=1.2, markersize=5, marker='o', color=base_colors[issp]) 

            axes[ii,1].loglog(grau_diam1[issp,:], WSM6preNorm_grau1[issp,0,iiwc,:], linestyle='--', linewidth=1.2, color='k')  
            axes[ii,1].loglog(grau_diam1[issp,:], WSM6Norm_grau1[issp,0,iiwc,:], linestyle='none', linewidth=1.2, markersize=5, marker='x', color='k') 

            #axes[ii,0].set_ylim([1e-7,1e10])
            #axes[ii,1].set_ylim([1e-7,1e10])
            axes[ii,0].grid(True)
            axes[ii,1].grid(True)
            axes[ii,0].set_title(f'Snow (liu:{issp}), iwc: {iwc_val[ii]}')
            axes[ii,1].set_title(f'Grau (liu:{issp}), iwc: {iwc_val[ii]}')
    
            ii=ii+1
    
        axes[2,0].set_ylabel('nd(D) [1/m4]')
        axes[2,0].set_xlabel('D [m]')
        plt.suptitle(f'Compare experiments (Temp: {temperature}K) and Liu: {issp}', fontweight='bold', fontsize=14)
        plt.show()
        fig.savefig(plotpath+f'/RTTOV/nD_perSSP_comparison_liu{issp}.png', dpi=300,transparent=False)
        
    #-- plot delta n(d) for both experiments
    fig, axes = plt.subplots(nrows=3, ncols=2, constrained_layout=True,figsize=[16,8])        
    for issp in ssps_liu: 
        axes[2,0].plot([],[],color=base_colors[issp],label=f'Liu: {issp}') 
        
    axes[0,1].plot([],[], color='k', marker='o', markersize=5, label=f'WSM6') 
    axes[0,1].plot([],[], color='k', marker='x', markersize=3, label=f'eqMass WSM6') 

    axes[2,0].legend(ncol=3)
    axes[0,1].legend(ncol=2)
            
    ii = 0 
    for iiwc in iwc: 
        for issp in ssps_liu:

            axes[ii,0].loglog(snow_diam[issp,:],  ignore_diffzero(WSM6Norm_snow[issp,0,iiwc,:]-WSM6preNorm_snow[issp,0,iiwc,:]), marker='o', markersize=5, linewidth=1.2, color=base_colors[issp]) 
            axes[ii,0].loglog(snow_diam1[issp,:], ignore_diffzero(WSM6Norm_snow1[issp,0,iiwc,:]-WSM6preNorm_snow1[issp,0,iiwc,:]), marker='x', markersize=3, linewidth=1.2, color=base_colors[issp])
            
            axes[ii,1].loglog(grau_diam[issp,:],  ignore_diffzero(WSM6Norm_grau[issp,0,iiwc,:]-WSM6preNorm_grau[issp,0,iiwc,:]), marker='o', markersize=5, linewidth=1.2, color=base_colors[issp]) 
            axes[ii,1].loglog(grau_diam1[issp,:], ignore_diffzero(WSM6Norm_grau1[issp,0,iiwc,:]-WSM6preNorm_grau1[issp,0,iiwc,:]), marker='x', markersize=3, linewidth=1.2, color=base_colors[issp]) 

            #axes[ii,1].set_ylim([1e-4,1e2])
            #axes[ii,0].set_ylim([1e-4,1e2])
            axes[ii,0].grid(True)
            axes[ii,1].grid(True)
            
            axes[ii,0].set_title(f'snow iwc: {iwc_val[ii]} g/m3')
            axes[ii,1].set_title(f'grau iwc: {iwc_val[ii]} g/m3')
            
        ii=ii+1
    
    axes[2,0].set_ylabel(r'$\Delta$ nd(D) due to normalization [1/m4]')
    axes[2,0].set_xlabel('D [m]')
    plt.suptitle(f'$\Delta$n(D) due to Norm. for (Temp: {temperature}K)', fontweight='bold', fontsize=14)
    plt.show()
    fig.savefig(plotpath+'/RTTOV/nD_all_deltaNorm_comparison.png', dpi=300,transparent=False)    


    return
#------------------------------------------------
#------------------------------------------------
def plt_aDb(plotpath):
    
    
    diam = np.linspace(1e-6, 1e-2, num=100)
    # Get the same colors as other plots:
    base_colors = sns.color_palette('Paired')       
    
    a = [37.09, 116.12, 229.66, 122.66, 32.36, 0.32,  0.06,  0.07, 0.09, 0.002, 0.01]
    b = [3.00, 3.00, 3.00, 3.00, 3.00, 2.37, 2.12,  2.12,  2.13,  1.58, 1.90]
                         
    rhos = 300
    rhog = 500
    
    asnow = rhos*(3.14/6) 
    agrau = rhog*(3.14/6) 
    
    fig = plt.figure(figsize=(6, 6))
    
    for issp in range(11):
        plt.loglog(diam/1e-3,  a[issp]*diam ** b[issp] ,linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'Liu: {issp}') 
    
    plt.loglog(diam/1e-3, asnow*diam**3 ,linestyle='-', linewidth=1.2, color='k', label=r'$\rho_s$=0.1')
    plt.loglog(diam/1e-3, agrau*diam**3 ,linestyle='--', linewidth=1.2, color='k', label=r'$\rho_g$=0.5')

    plt.ylabel('mass [kg]')
    plt.xlabel('diameter [mm]')
    plt.grid(True)
    plt.legend(ncol=2)
    plt.show()
    fig.savefig(plotpath+'/RTTOV/mDb_general.png', dpi=300,transparent=False)
    
    return


#------------------------------------------------
def plt_aDb_density(plotpath):
    
    isspname = ['Long hex col.','Short hex col.','Block hex col.','Thick hex col','Thin hex col.','3b ros.','4b ros.','5b ros.','6b ros.','Sector','Dendrite']     

    diam = np.linspace(1e-6, 1e-2, num=100)
    # Get the same colors as other plots:
    base_colors = sns.color_palette('Paired')       
    
    a = [37.09, 116.12, 229.66, 122.66, 32.36, 0.32,  0.06,  0.07, 0.09, 0.002, 0.01]
    b = [3.00, 3.00, 3.00, 3.00, 3.00, 2.37, 2.12,  2.12,  2.13,  1.58, 1.90]
                         
    rhos = 300
    rhog = 500
    
    asnow = rhos*(3.14/6) 
    agrau = rhog*(3.14/6) 
    
    fig = plt.figure(figsize=(6, 6))
    
    do_thisSSPs = [2,5,9,10]
    
    for issp in do_thisSSPs:
        if issp == 10:
            base_colors[issp]=base_colors[issp+1]
        plt.loglog(diam/1e-3,  a[issp]*diam ** b[issp] ,linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'{isspname[issp]}')  
    
    plt.loglog(diam/1e-3, asnow*diam**3 ,linestyle='-', linewidth=1.2, color='k', label=r'$\rho_s$=0.1')
    plt.loglog(diam/1e-3, agrau*diam**3 ,linestyle='--', linewidth=1.2, color='k', label=r'$\rho_g$=0.5')

    plt.ylabel('mass [kg]')
    plt.xlabel('diameter [mm]')
    plt.grid(True)
    
    plt.legend(ncol=1, fontsize=10)
    plt.show()
    fig.savefig(plotpath+'/RTTOV/rho_4wg_general.png', dpi=300,transparent=False)
    
    return


#------------------------------------------------
#------------------------------------------------
def get_experiment(WSM6path): 

    WSM6preNorm_snow = []
    WSM6Norm_snow    = []
    WSM6preNorm_grau = []
    WSM6Norm_grau    = []
        
    grau_diam = []
    snow_diam = []
    
    for issp in range(11):    
        WSM6preNorm_snow.append(  readf_txt(WSM6path+f'nd_WSM6snow_preNormalization_liu{issp}.txt') )
        WSM6Norm_snow.append( readf_txt(WSM6path+f'nd_WSM6snow_wNormalization_liu{issp}.txt') )
        WSM6preNorm_grau.append(  readf_txt(WSM6path+f'nd_WSM6grau_preNormalization_liu{issp}.txt') )
        WSM6Norm_grau.append( readf_txt(WSM6path+f'nd_WSM6grau_wNormalization_liu{issp}.txt') )
    
        snow_diam.append( readf_diam_txt( WSM6path+f'diam_WSM6snow_liu{issp}.txt' ))
        grau_diam.append( readf_diam_txt( WSM6path+f'diam_WSM6grau_liu{issp}.txt' ))
    
    WSM6preNorm_snow = np.array(WSM6preNorm_snow) 
    WSM6Norm_snow    = np.array(WSM6Norm_snow) 
    WSM6preNorm_grau = np.array(WSM6preNorm_grau) 
    WSM6Norm_grau    = np.array(WSM6Norm_grau) 
    
    grau_diam = np.array(grau_diam) 
    snow_diam = np.array(snow_diam) 
    
    WSM6preNorm_snow[WSM6preNorm_snow == 0] = np.nan
    
    WSM6preNorm_snow[WSM6preNorm_snow == 0] = np.nan
    WSM6preNorm_grau[WSM6preNorm_grau == 0] = np.nan
    WSM6Norm_snow[WSM6Norm_snow == 0] = np.nan
    WSM6Norm_grau[WSM6Norm_grau == 0] = np.nan
    
    return snow_diam, grau_diam, WSM6preNorm_snow, WSM6Norm_snow, WSM6preNorm_grau, WSM6Norm_grau

#------------------------------------------------
#------------------------------------------------
def get_bulk(path, itype): 

    bulk_ext = []
    bulk_asm    = []
    bulk_ssa = []
    
    for issp in range(11):    
        bulk_ext.append(  read_bulk_txt(path+'bulk_'+itype+f'_ext_WSM6_liu{issp}.txt') )
        bulk_ssa.append(  read_bulk_txt(path+'bulk_'+itype+f'_ssa_WSM6_liu{issp}.txt') )
        bulk_asm.append(  read_bulk_txt(path+'bulk_'+itype+f'_asm_WSM6_liu{issp}.txt') )
        
    bulk_ext = np.array(bulk_ext) 
    bulk_ssa = np.array(bulk_ssa) 
    bulk_asm = np.array(bulk_asm) 
        
    return bulk_ext, bulk_ssa, bulk_asm

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_nd_poster(plotpath):
    
    isspname = ['Long hex col.','Short hex col.','Block hex col.','Thick hex col','Thin hex col.','3b ros.','4b ros.','5b ros.','6b ros.','Sector','Dendrite']     
    do_thisSSPs = [2,3,8,7,9,10]

    plt.matplotlib.rc('font', family='serif', size = 18)
    plt.rcParams['xtick.labelsize']=18
    plt.rcParams['ytick.labelsize']=18 

    n0  =  2e6*np.exp(0.12*(273.15-263))
    lam = (3.14*100*n0/0.1*1e3)**(1/4)
    diams_ = snow_diam[0,:]

    ng0  =  4e6
    lamg = (3.14*500*ng0/0.1*1e3)**(1/4)

    do_thisSSPs = [2,3,8,7,9,10]
    base_colors = sns.color_palette('Paired')   
        
    fig, axes = plt.subplots(nrows=2, ncols=1, constrained_layout=True,figsize=[6,10])

    for issp in do_thisSSPs:     
        if issp == 10:
            base_colors[issp]=base_colors[issp+1]
        axes[1].plot([],[],linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'{isspname[issp]}')     
    axes[1].plot([],[],linestyle='-', marker='o', markersize=3, linewidth=1.2, color='gray', label='eqMass_Mass')                  
    axes[1].plot([], [], '-k',  linewidth=1.2, label='WSM6')   
    axes[1].legend(ncol=2, loc='lower left', fontsize=14)

    # snow wsm6 
    for issp in do_thisSSPs: 

        # snow wsm6    
        axes[0].loglog(snow_diam[issp,:], WSM6Norm_snow[issp,0,201,:], linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'Liu: {isspname}')     
        # snow eq mass wsm6 
        axes[0].loglog(snow_diam1[0,:], WSM6Norm_snow1[0,0,201,:], linestyle='-', marker='o', markersize=3, linewidth=1.2, color='gray', label='eqMass_WSM6')                  
        axes[0].loglog(snow_diam[issp,:], WSM6Norm_snow2[issp,0,201,:], linestyle='--', linewidth=1.2, color=base_colors[issp], label=f'Liu: {isspname}')     

        axes[0].loglog(diams_, n0*np.exp(-lam*diams_), '-k',  linewidth=1.2, label='WSM6')   

        axes[0].set_title(f'WSM6 Snow (0.1 g/m3, 263K)')
        axes[0].set_ylabel('nd(D) [1/m4]')
        axes[0].set_ylim([1e-1,1e10])
        
        # grau wsm6    
        axes[1].loglog(grau_diam[issp,:], WSM6Norm_grau[issp,0,201,:], linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'Liu: {issp}')       
        # grau eq mass wsm6 
        axes[1].loglog(grau_diam1[0,:], WSM6Norm_grau1[0,0,201,:], linestyle='-', marker='o', markersize=5, linewidth=1.2, color='gray', label='eqMass_WSM6')         
        axes[1].loglog(diams_, ng0*np.exp(-lamg*diams_), '-k',  linewidth=1.2, label='WSM6')   
        axes[1].loglog(grau_diam[issp,:], WSM6Norm_grau2[issp,0,201,:], linestyle='--', linewidth=1.2, color=base_colors[issp], label=f'Liu: {issp}')       

        axes[1].set_title(f'WSM6 Graupel (0.1 g/m3, 263K)')
        axes[1].set_ylabel('nd(D) [1/m4]')
        axes[1].set_xlabel('D [m]')
        axes[1].set_ylim([1,1e7])
        
        axes[0].grid(True)
        axes[1].grid(True)
    plt.show()
    fig.savefig(plotpath+'/RTTOV/nd_summary_general.png', dpi=300,transparent=False)


    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_nd_poster_lessspecies(plotpath):
    
    isspname = ['Long hex col.','Short hex col.','Block hex col.','Thick hex col','Thin hex col.','3b ros.','4b ros.','5b ros.','6b ros.','Sector','Dendrite']     
    do_thisSSPs = [2,5,9,10]

    plt.matplotlib.rc('font', family='serif', size = 18)
    plt.rcParams['xtick.labelsize']=18
    plt.rcParams['ytick.labelsize']=18 

    n0  =  2e6*np.exp(0.12*(273.15-263))
    lam = (3.14*100*n0/0.1*1e3)**(1/4)
    diams_ = snow_diam[0,:]

    ng0  =  4e6
    lamg = (3.14*500*ng0/0.1*1e3)**(1/4)

    base_colors = sns.color_palette('Paired')   
        
    fig, axes = plt.subplots(nrows=2, ncols=1, constrained_layout=True,figsize=[12,10])

    for issp in do_thisSSPs:     
        if issp == 10:
            base_colors[issp]=base_colors[issp+1]
        axes[0].plot([],[],linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'{isspname[issp]}')    
    fig.legend(loc='lower left', ncol=4, bbox_to_anchor=(0.1,-0.08))
    
    axes[1].plot([], [], 'k',  linestyle=':', linewidth=2, label=r'n$_0$exp(-$\lambda$D)')   
    axes[1].plot([], [], '-k',  linewidth=1.2, label='iwc rescaling')   
    axes[1].plot([],[],linestyle='-', marker='o', markersize=3, linewidth=1.2, color='gray', label='eqMass')                  
    axes[1].plot([], [], '--k',  linewidth=1.2, label='M(D) Sieron et al. 2018') 
    axes[1].legend(ncol=1, loc='lower left', fontsize=14)
    axes[1].legend(loc='lower left')

    # snow wsm6 
    for issp in do_thisSSPs: 

        # snow wsm6    
        axes[0].loglog(snow_diam[issp,:], WSM6Norm_snow[issp,0,201,:], linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'Liu: {isspname}')     
        # snow eq mass wsm6 
        axes[0].loglog(snow_diam1[0,:], WSM6Norm_snow1[0,0,201,:], linestyle='-', marker='o', markersize=3, linewidth=1.2, color='gray', label='eqMass_WSM6')                  
        axes[0].loglog(snow_diam[issp,:], WSM6Norm_snow2[issp,0,201,:], linestyle='--', linewidth=1.2, color=base_colors[issp], label=f'Liu: {isspname}')     

        axes[0].loglog(diams_, n0*np.exp(-lam*diams_), 'k',  linestyle=':', linewidth=2, label='WSM6')   

        axes[0].set_title(f'WSM6 consistency for Snow (0.1 g/m3, 263K)')
        axes[0].set_ylabel('nd(D) [1/m4]')
        axes[0].set_ylim([1,1e10])
        
        # grau wsm6    
        axes[1].loglog(grau_diam[issp,:], WSM6Norm_grau[issp,0,201,:], linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'Liu: {issp}')       
        # grau eq mass wsm6 
        axes[1].loglog(grau_diam1[0,:], WSM6Norm_grau1[0,0,201,:], linestyle='-', marker='o', markersize=5, linewidth=1.2, color='gray', label='eqMass_WSM6')         
        axes[1].loglog(diams_, ng0*np.exp(-lamg*diams_),'k',  linestyle=':', linewidth=2, label='WSM6')   
        axes[1].loglog(grau_diam[issp,:], WSM6Norm_grau2[issp,0,201,:], linestyle='--', linewidth=1.2, color=base_colors[issp], label=f'Liu: {issp}')       

        axes[1].set_title(f'WSM6 consistency for Graupel (0.1 g/m3, 263K)')
        axes[1].set_ylabel('nd(D) [1/m4]')
        axes[1].set_xlabel('D [m]')
        axes[1].set_ylim([1,1e7])
        
        axes[0].grid(True)
        axes[1].grid(True)
        axes[0].set_xlim([1e-4,1e-2])
        axes[1].set_xlim([1e-4,1e-2])


    plt.show()
    fig.savefig(plotpath+'/RTTOV/nd_summary_general_lesspecies.png', dpi=300,transparent=False)

    return


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def plot_nd_poster_lessspecies_3subplots(plotpath):
    
    isspname = ['Long hex col.','Short hex col.','Block hex col.','Thick hex col','Thin hex col.','3b ros.','4b ros.','5b ros.','6b ros.','Sector','Dendrite']     
    do_thisSSPs = [2,5,9,10]

    plt.matplotlib.rc('font', family='serif', size = 18)
    plt.rcParams['xtick.labelsize']=18
    plt.rcParams['ytick.labelsize']=18 

    n0  =  2e6*np.exp(0.12*(273.15-263))
    lam = (3.14*100*n0/0.1*1e3)**(1/4)
    diams_ = snow_diam[0,:]

    ng0  =  4e6
    lamg = (3.14*500*ng0/0.1*1e3)**(1/4)

    base_colors = sns.color_palette('Paired')   
        
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=[15,10])
    fig.subplots_adjust(hspace=0.3)
    
    for issp in do_thisSSPs:     
        if issp == 10:
            base_colors[issp]=base_colors[issp+1]
        axes[0,0].plot([],[],linestyle='-', linewidth=2, color=base_colors[issp], label=f'{isspname[issp]}')    
    axes[0,0].plot([], [], color='k',  linestyle=':', linewidth=2, label=r'n$_0$exp(-$\lambda$D)')   
    fig.legend(loc='lower left', ncol=5, bbox_to_anchor=(0.1,-0.05))

    # axes[1].legend(ncol=1, loc='lower left', fontsize=14)
    # axes[1].legend(loc='lower left')

    # snow wsm6 
    for issp in do_thisSSPs: 

        # snow wsm6    
        axes[0,0].loglog(snow_diam[issp,:], WSM6Norm_snow[issp,0,201,:], linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'{isspname[issp]}')    
        # snow eq mass wsm6 
        axes[0,1].loglog(snow_diam1[0,:], WSM6Norm_snow1[0,0,201,:], linestyle='-', markersize=3, linewidth=1.2, color=base_colors[issp], label=f'{isspname[issp]}')                   
        axes[0,2].loglog(snow_diam[issp,:], WSM6Norm_snow2[issp,0,201,:], linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'{isspname[issp]}')    
        axes[0,2].loglog(snow_diam1[0,:], WSM6Norm_snow1[0,0,201,:], linestyle='-', marker='o', markersize=3, linewidth=1.2, alpha=0.2, color='gray', label=f'{isspname[issp]}')                   


        axes[0,0].loglog(diams_, n0*np.exp(-lam*diams_), 'k',  linestyle=':', linewidth=2, label='WSM6 n(D)')   
        axes[0,1].loglog(diams_, n0*np.exp(-lam*diams_), 'k',  linestyle=':', linewidth=2, label='WSM6 n(D)')   
        axes[0,2].loglog(diams_, n0*np.exp(-lam*diams_), 'k',  linestyle=':', linewidth=2, label='WSM6 n(D)')   
    
        axes[0,0].set_title(f'iwc rescaling')
        axes[0,1].set_title(f'eqMass')
        axes[0,2].set_title(f'M(D) Sieron et al. (2018)')

        axes[0,0].set_ylabel('nd(D) [1/m4]')
        
        axes[0,0].set_ylim([1,1e10])
        axes[0,1].set_ylim([1,1e10])
        axes[0,2].set_ylim([1,1e10])
        axes[0,0].grid(True)
        axes[0,1].grid(True)
        axes[0,2].grid(True)
        axes[0,0].set_xlim([1e-4,1e-2])
        axes[0,1].set_xlim([1e-4,1e-2])
        axes[0,2].set_xlim([1e-4,1e-2])
                
        fig.text(0.5, 0.94, r'WSM6 Snow ($\rho$=0.1g/m3, T=263K)', ha='center', va='center', fontsize=16,fontweight='bold')
        
        # # grau wsm6    
        # axes[1,0].loglog(grau_diam[issp,:], WSM6Norm_grau[issp,0,201,:], linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'Liu: {issp}')       
        # # grau eq mass wsm6 
        # axes[1,0].loglog(grau_diam1[0,:], WSM6Norm_grau1[0,0,201,:], linestyle='-', marker='o', markersize=5, linewidth=1.2, color='gray', label='eqMass_WSM6')         
        # axes[1,0].loglog(diams_, ng0*np.exp(-lamg*diams_),'k',  linestyle=':', linewidth=2, label='WSM6')   
        # axes[1,0].loglog(grau_diam[issp,:], WSM6Norm_grau2[issp,0,201,:], linestyle='--', linewidth=1.2, color=base_colors[issp], label=f'Liu: {issp}')       

        axes[1,0].loglog(grau_diam[issp,:], WSM6Norm_grau[issp,0,201,:], linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'{isspname[issp]}')     
        axes[1,1].loglog(grau_diam1[0,:], WSM6Norm_grau1[0,0,201,:], linestyle='-',  markersize=3, linewidth=1.2, color=base_colors[issp], label=f'{isspname[issp]}')                   
        axes[1,2].loglog(grau_diam[issp,:], WSM6Norm_grau2[0,0,201,:], linestyle='-', markersize=3, linewidth=1.2, color=base_colors[issp], label=f'{isspname[issp]}')                      
        #axes[1,2].loglog(grau_diam1[0,:], WSM6Norm_grau1[0,0,201,:], linestyle='-', marker='o', markersize=3, linewidth=1.2, alpha=0.2, color='gray', label=f'{isspname[issp]}')                   

        axes[1,0].loglog(diams_, ng0*np.exp(-lamg*diams_), 'k',  linestyle=':', linewidth=2, label='WSM6 n(D)')   
        axes[1,1].loglog(diams_, ng0*np.exp(-lamg*diams_), 'k',  linestyle=':', linewidth=2, label='WSM6 n(D)')   
        axes[1,2].loglog(diams_, ng0*np.exp(-lamg*diams_), 'k',  linestyle=':', linewidth=2, label='WSM6 n(D)')   
    
        axes[1,0].set_ylabel('nd(D) [1/m4]')
        axes[1,0].set_xlabel('D [m]')

        axes[1,0].set_ylim([1,1e10])
        axes[1,1].set_ylim([1,1e10])
        axes[1,2].set_ylim([1,1e10])
        axes[1,0].grid(True)
        axes[1,1].grid(True)
        axes[1,2].grid(True)
        axes[1,0].set_xlim([1e-4,1e-2])
        axes[1,1].set_xlim([1e-4,1e-2])
        axes[1,2].set_xlim([1e-4,1e-2])
    
        fig.text(0.5, 0.47, r'WSM6 grau ($\rho$=0.4g/m3, T=263K)', ha='center', va='center', fontsize=16,fontweight='bold')

    plt.show()
    fig.savefig(plotpath+'/RTTOV/nd4_wgmeeting.png', dpi=300,transparent=False, bbox_inches='tight')   

    return


#--------------------------------------------------------------------------------------------
#                       BULKS
#--------------------------------------------------------------------------------------------
def plot_bulks_together(plotpath):

    do_thisSSPs = [2,3,8,7,9,10]
    isspname = ['Long hex col.','Short hex col.','Block hex col.','Thick hex col','Thin hex col.','3b ros.','4b ros.','5b ros.','6b ros.','Sector','Dendrite']     

    plt.matplotlib.rc('font', family='serif', size = 18)
    plt.rcParams['xtick.labelsize']=18
    plt.rcParams['ytick.labelsize']=18 
    
    iwc = []
    for i_lwc in range(n_lwc):     
        iwc.append(get_lwc(i_lwc)) 
        
    temperature = 263
    
    # Get the same colors as other plots:
    ssps_liu = np.arange(0,11,1)
    base_colors = sns.color_palette('Paired')         
        
    #-- plotn(d) norm and not-normed one figure per liu ssp: 
    fig, axes = plt.subplots(nrows=3, ncols=2, constrained_layout=True,figsize=[12,10])

    for issp in do_thisSSPs:    
        
        if issp == 10:
            base_colors[issp]=base_colors[issp+1]
            
        axes[0,0].plot([],[], linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'{isspname[issp]}')    
    axes[0,0].plot([],[], linestyle='-', linewidth=1.2, color='gray', label='WSM6') 
    axes[0,0].plot([],[], linestyle='--', linewidth=1.2, color='gray', label='eqMassWSM6') 

    fig.legend(ncol=3, loc='lower center', bbox_to_anchor=(0.5,-0.13), fontsize=14)


    for issp in do_thisSSPs:    

        axes[0,0].loglog(iwc, bulk_snow_ext[issp,:],  linestyle='-', linewidth=1.2, color=base_colors[issp])          
        axes[1,0].semilogx(iwc, bulk_snow_ssa[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp])          
        axes[2,0].semilogx(iwc, bulk_snow_asm[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp])  
        

        axes[0,0].loglog(iwc, bulk_snow_ext1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp])          
        axes[1,0].semilogx(iwc, bulk_snow_ssa1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp])          
        axes[2,0].semilogx(iwc, bulk_snow_asm1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp])          
        axes[0,0].set_title(f'Snow')


        axes[0,1].loglog(iwc, bulk_grau_ext[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp]) 
        axes[1,1].semilogx(iwc, bulk_grau_ssa[issp,:],  linestyle='-', linewidth=1.2, color=base_colors[issp]) 
        axes[2,1].semilogx(iwc, bulk_grau_asm[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp]) 
        

        axes[0,1].loglog(iwc, bulk_grau_ext1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp]) 
        axes[1,1].semilogx(iwc, bulk_grau_ssa1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp]) 
        axes[2,1].semilogx(iwc, bulk_grau_asm1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp]) 
        axes[0,1].set_title(f'Graupel')
            
        axes[0,0].set_ylabel('bulk extinction [1/m]')
        axes[1,0].set_ylabel('bulk ssa [0-1] ')     # single scattering albedo
        axes[2,0].set_ylabel('bulk asm [0-1]')   # asymetry parameter

        axes[2,0].set_xlabel('iwc [g/m3]')
        for i in range(3):
            for j in range(2):
                axes[i,j].grid(True)
                axes[i,j].grid(True)
        
    plt.suptitle(f'Compare experiments (Temp: {temperature}K) and Liu: {issp} for 89GHz', fontweight='bold')
    #plt.tight_layout()
    plt.show()
    fig.savefig(plotpath+f'/RTTOV/bulk_scattering_iwc_someliu.png', dpi=300,transparent=False, bbox_inches='tight')

    return

def plot_bulks_together_lessspecies(plotpath):

    isspname = ['Long hex col.','Short hex col.','Block hex col.','Thick hex col','Thin hex col.','3b ros.','4b ros.','5b ros.','6b ros.','Sector','Dendrite']     
    do_thisSSPs = [2,5,9,10]
    
    iwc = []
    for i_lwc in range(n_lwc):     
        iwc.append(get_lwc(i_lwc)) 
    temperature = 263
   
    # Get the same colors as other plots:
    ssps_liu = np.arange(0,11,1)
    base_colors = sns.color_palette('Paired')    
    for issp in do_thisSSPs:    
        if issp == 10:
            base_colors[issp]=base_colors[issp+1]
                     
    #-- plotn(d) norm and not-normed one figure per liu ssp: 
    fig, axes = plt.subplots(nrows=3, ncols=2, constrained_layout=True,figsize=[14,10])
    ssp_handles = []
    for issp in do_thisSSPs:                
        line, = axes[0,0].plot([],[], linestyle='-', linewidth=1.2, color=base_colors[issp], label=f'{isspname[issp]}')    
        ssp_handles.append(line)
    #axes[0,0].plot([],[], linestyle='-', linewidth=1.2, color='gray', label='WSM6') 
    #axes[0,0].plot([],[], linestyle='--', linewidth=1.2, color='gray', label='eqMassWSM6') 
    method_handles = []
    method_handles.append(axes[0,0].plot([], [], color='gray',  linestyle='-', linewidth=2, label='iwc rescaling')[0])
    method_handles.append(axes[0,0].plot([], [], color='gray',  linestyle='-',  marker='o', markersize=3, linewidth=1.2, label='eqMass')[0])
    method_handles.append(axes[0,0].plot([], [], color='gray',  linestyle='--', linewidth=2, label='M(D) Sieron et al. 2018')[0])            
    #fig.legend(ncol=4, loc='lower center', bbox_to_anchor=(0.4,-0.08), fontsize=14)

    # -- First legend: SSPs
    legend1 = fig.legend(handles=ssp_handles, ncol=1, loc='lower center', bbox_to_anchor=(0.1, -0.16))
    
    # -- Second legend: Methods
    legend2 = fig.legend(handles=method_handles, ncol=1, loc='lower center', bbox_to_anchor=(0.5, -0.16))
    
    # -- Optional: add both legends to the figure manually
    fig.add_artist(legend1)
    fig.add_artist(legend2)

    

    for issp in do_thisSSPs:    
        axes[0,0].loglog(iwc, bulk_snow_ext[issp,:],  linestyle='-', linewidth=1.2, color=base_colors[issp])          
        axes[1,0].semilogx(iwc, bulk_snow_ssa[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp])          
        axes[2,0].semilogx(iwc, bulk_snow_asm[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp])  

        axes[0,0].loglog(iwc, bulk_snow_ext1[issp,:], linestyle='-', marker='o', markersize=3, linewidth=1.2, color=base_colors[issp])          
        axes[1,0].semilogx(iwc, bulk_snow_ssa1[issp,:], linestyle='-', marker='o',  markersize=3, linewidth=1.2, color=base_colors[issp])          
        axes[2,0].semilogx(iwc, bulk_snow_asm1[issp,:], linestyle='-',marker='o',  markersize=3, linewidth=1.2, color=base_colors[issp])          

        axes[0,0].loglog(iwc, bulk_snow_ext2[issp,:],  linestyle='--', linewidth=1.2, color=base_colors[issp])          
        axes[1,0].semilogx(iwc, bulk_snow_ssa2[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp])          
        axes[2,0].semilogx(iwc, bulk_snow_asm2[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp])  


        axes[0,1].loglog(iwc, bulk_grau_ext[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp]) 
        axes[1,1].semilogx(iwc, bulk_grau_ssa[issp,:],  linestyle='-', linewidth=1.2, color=base_colors[issp]) 
        axes[2,1].semilogx(iwc, bulk_grau_asm[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp]) 
        
        
        axes[0,1].loglog(iwc, bulk_grau_ext1[issp,:], linestyle='--', marker='o',  markersize=3, linewidth=1.2, color=base_colors[issp]) 
        axes[1,1].semilogx(iwc, bulk_grau_ssa1[issp,:], linestyle='--',marker='o',  markersize=3, linewidth=1.2, color=base_colors[issp]) 
        axes[2,1].semilogx(iwc, bulk_grau_asm1[issp,:], linestyle='--',marker='o',  markersize=3, linewidth=1.2, color=base_colors[issp]) 
            

        axes[0,1].loglog(iwc, bulk_grau_ext2[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp]) 
        axes[1,1].semilogx(iwc, bulk_grau_ssa2[issp,:],  linestyle='--', linewidth=1.2, color=base_colors[issp]) 
        axes[2,1].semilogx(iwc, bulk_grau_asm2[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp])         
        
        axes[0,0].set_ylabel('bulk extinction [1/m]')
        axes[1,0].set_ylabel('bulk ssa [0-1] ')     # single scattering albedo
        axes[2,0].set_ylabel('bulk asm [0-1]')   # asymetry parameter

        axes[2,0].set_xlabel('iwc [g/m3]')
        axes[0,0].set_title(f'Snow')
        axes[0,1].set_title(f'Graupel')

        for i in range(3):
            for j in range(2):
                axes[i,j].grid(True)
                axes[i,j].grid(True)
        
    plt.suptitle(f'Bulk optical properties under different consistency approaches ({temperature}K and 89GHz)', fontweight='bold')
    #plt.tight_layout()
    plt.show()
    fig.savefig(plotpath+f'/RTTOV/bulk_scattering_iwc_someliu_poster_legend.png', dpi=300,transparent=False, bbox_inches='tight')

    return

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
def plot_bulks():
    
    plt.matplotlib.rc('font', family='serif', size = 18)
    plt.rcParams['xtick.labelsize']=18
    plt.rcParams['ytick.labelsize']=18 
    
    iwc = []
    for i_lwc in range(n_lwc):     
        iwc.append(get_lwc(i_lwc)) 
        
    temperature = 263
    
    # Get the same colors as other plots:
    ssps_liu = np.arange(0,11,1)
    base_colors = sns.color_palette('Paired')         
        
    #-- plotn(d) norm and not-normed one figure per liu ssp: 
    for issp in ssps_liu:    
        
        if issp == 10:
            base_colors[issp]=base_colors[issp+1]
            
        fig, axes = plt.subplots(nrows=3, ncols=2, constrained_layout=True,figsize=[12,10])
        
        # axes[2,0].plot([],[], linestyle='-', linewidth=1.2, color='k', label='WSM6') 
        # axes[2,0].plot([],[], linestyle='--', linewidth=1.2, color='k', label='eqMassWSM6') 
        # axes[2,0].plot([],[], linestyle='none', marker='o', markersize=5, color='k', label='WSM6 (renorm.)') 
        # axes[2,0].plot([],[], linestyle='none', marker='x', markersize=5, color='k', label='eqMassWSM6 (renorm.)') 
        # axes[2,0].legend(ncol=2, loc='lower left')

        axes[0,0].loglog(iwc, bulk_snow_ext[issp,:],  linestyle='-', linewidth=1.2, color=base_colors[issp])          
        axes[1,0].semilogx(iwc, bulk_snow_ssa[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp])          
        axes[2,0].semilogx(iwc, bulk_snow_asm[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp])  
        

        axes[0,0].loglog(iwc, bulk_snow_ext1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp])          
        axes[1,0].semilogx(iwc, bulk_snow_ssa1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp])          
        axes[2,0].semilogx(iwc, bulk_snow_asm1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp])          

        axes[0,0].loglog(iwc, bulk_snow_ext2[issp,:], linestyle='-', marker='x', linewidth=1.2, color=base_colors[issp])          
        axes[1,0].semilogx(iwc, bulk_snow_ssa2[issp,:], linestyle='-',marker='x',  linewidth=1.2, color=base_colors[issp])          
        axes[2,0].semilogx(iwc, bulk_snow_asm2[issp,:], linestyle='-', marker='x', linewidth=1.2, color=base_colors[issp])   

        axes[0,0].set_title(f'Snow (liu:{issp})')


        axes[0,1].loglog(iwc, bulk_grau_ext[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp])
        axes[1,1].semilogx(iwc, bulk_grau_ssa[issp,:],  linestyle='-', linewidth=1.2, color=base_colors[issp]) 
        axes[2,1].semilogx(iwc, bulk_grau_asm[issp,:], linestyle='-', linewidth=1.2, color=base_colors[issp]) 
        

        axes[0,1].loglog(iwc, bulk_grau_ext1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp]) 
        axes[1,1].semilogx(iwc, bulk_grau_ssa1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp]) 
        axes[2,1].semilogx(iwc, bulk_grau_asm1[issp,:], linestyle='--', linewidth=1.2, color=base_colors[issp]) 
        axes[0,1].set_title(f'Graupel (liu:{issp})')

        axes[0,1].loglog(iwc, bulk_grau_ext2[issp,:], linestyle='-', marker='x', linewidth=1.2, color=base_colors[issp])          
        axes[1,1].semilogx(iwc, bulk_grau_ssa2[issp,:], linestyle='-',marker='x',  linewidth=1.2, color=base_colors[issp])          
        axes[2,1].semilogx(iwc, bulk_grau_asm2[issp,:], linestyle='-', marker='x', linewidth=1.2, color=base_colors[issp])   



        axes[0,0].set_ylabel('bulk extinction [1/m]')
        axes[1,0].set_ylabel('bulk ssa [0-1] ')     # single scattering albedo
        axes[2,0].set_ylabel('bulk asm [0-1]')   # asymetry parameter

        axes[2,0].set_xlabel('iwc [g/m3]')
        for i in range(3):
            for j in range(2):
                axes[i,j].grid(True)
                axes[i,j].grid(True)
        
        plt.show()
        plt.suptitle(f'Compare experiments (Temp: {temperature}K) and Liu: {issp} for 89GHz', fontweight='bold', fontsize=14)
        plt.show()
        fig.savefig(plotpath+f'/RTTOV/bulk_scattering_iwc_liu{issp}.png', dpi=300,transparent=False)

    return


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
n_temp = 1   # n_temp = 4 leads to 177, 178, 179, 180. quite high value 177! for ice actually 
n_lwc  = 401    # n_lwc  = 401 leads to values between 9.77e-7 and 0.0098 [kg/m3]             
temp   = get_temp(n_temp) 

#------------------------------------------------
# Get the array of lwc 
lwc = []
for i_lwc in range(n_lwc):     
    lwc.append(get_lwc(i_lwc)) 
    
#------------------------------------------------
# mas = aD**b plots: 
do_this = 0
if do_this == 1: 
    plt_aDb('/Users/vito.galligani/Work/Studies/HAIL_20181110/SCATANAL')

#------------------------------------------------
mpath           = '/Users/vito.galligani/Work/Studies/HAIL_20181110/RTTOVinout/mw_scatt_anal/'
WSM6path        = mpath+'rttov_hydro_scatt_WSM6/'
eqMassWSM6path  = mpath+'rttov_hydro_scatt_eqMassWSM6/'
sieronpath      = mpath+'rttov_hydro_scatt_sieron/'

#------------------------------------------------
# 1a) Plot n(D) for snow and grau. Look at the effect of re-normalization on the PSD
snow_diam, grau_diam, WSM6preNorm_snow, WSM6Norm_snow, WSM6preNorm_grau, WSM6Norm_grau = get_experiment(WSM6path)
#plot_nd(snow_diam, grau_diam, WSM6preNorm_snow, WSM6Norm_snow, WSM6preNorm_grau, WSM6Norm_grau, 'WSM6')

#------------------------------------------------
# 1b) Plot n(D) for snow and grau. Look at the effect of re-normalization on the PSD
snow_diam1, grau_diam1, WSM6preNorm_snow1, WSM6Norm_snow1, WSM6preNorm_grau1, WSM6Norm_grau1 = get_experiment(eqMassWSM6path)
#plot_nd(snow_diam1, grau_diam1, WSM6preNorm_snow1, WSM6Norm_snow1, WSM6preNorm_grau1, WSM6Norm_grau1, 'eqMassWSM6')

snow_diam2, grau_diam2, WSM6preNorm_snow2, WSM6Norm_snow2, WSM6preNorm_grau2, WSM6Norm_grau2 = get_experiment(sieronpath)
#plot_nd(snow_diam2, grau_diam2, WSM6preNorm_snow2, WSM6Norm_snow2, WSM6preNorm_grau2, WSM6Norm_grau2, 'sieron')

#------------------------------------------------
# Think of a plot for the poster for n(D)
plotpath, folders = config_folders('laptop')
    
#plotpath = '/Users/vito.galligani/Work/Studies/HAIL_20181110/SCATANAL'    
#plotpath = os.path.join(plotpath,datetime.now().strftime('%d%m%Y'))
#if not os.path.exists(plotpath):
#    os.makedirs(plotpath)        

do_this = 0
if do_this == 1: 
    #plot_nd_poster(plotpath)
    plot_nd_poster_lessspecies(plotpath)

plot_nd_poster_lessspecies_3subplots(plotpath)
    
    

#------------------------------------------------
# 3) BULKS
WSM6path       = mpath+'rttov_hydro_scatt_WSM6/'
eqMassWSM6path = mpath+'rttov_hydro_scatt_eqMassWSM6/'
sieronpath = mpath+'rttov_hydro_scatt_sieron/'

bulk_grau_ext, bulk_grau_ssa, bulk_grau_asm = get_bulk(WSM6path,'grau') 
bulk_snow_ext, bulk_snow_ssa, bulk_snow_asm = get_bulk(WSM6path,'snow') 

bulk_grau_ext1, bulk_grau_ssa1, bulk_grau_asm1 = get_bulk(eqMassWSM6path,'grau') 
bulk_snow_ext1, bulk_snow_ssa1, bulk_snow_asm1 = get_bulk(eqMassWSM6path,'snow') 

bulk_grau_ext2, bulk_grau_ssa2, bulk_grau_asm2 = get_bulk(sieronpath,'grau') 
bulk_snow_ext2, bulk_snow_ssa2, bulk_snow_asm2 = get_bulk(sieronpath,'snow') 
    
do_this = 0
if do_this == 1: 
    plot_bulks_together_lessspecies(plotpath)            # missing soft sphere?     
    
    

do_this = 0
if do_this == 1: 
    #------------------------------------------------
    # 2) compare plots
    plot_nd_compareEqMass(WSM6path, eqMassWSM6path, sieronpath)
    
    #------------------------------------------------
    # 3) BULKS
    WSM6path       = mpath+'rttov_hydro_scatt_WSM6/'
    eqMassWSM6path = mpath+'rttov_hydro_scatt_eqMassWSM6/'
    sieronpath = mpath+'rttov_hydro_scatt_sieron/'
    
    
    bulk_grau_ext, bulk_grau_ssa, bulk_grau_asm = get_bulk(WSM6path,'grau') 
    bulk_snow_ext, bulk_snow_ssa, bulk_snow_asm = get_bulk(WSM6path,'snow') 
    
    bulk_grau_ext1, bulk_grau_ssa1, bulk_grau_asm1 = get_bulk(eqMassWSM6path,'grau') 
    bulk_snow_ext1, bulk_snow_ssa1, bulk_snow_asm1 = get_bulk(eqMassWSM6path,'snow') 
    
    bulk_grau_ext2, bulk_grau_ssa2, bulk_grau_asm2 = get_bulk(sieronpath,'grau') 
    bulk_snow_ext2, bulk_snow_ssa2, bulk_snow_asm2 = get_bulk(sieronpath,'snow') 
    
         
    plot_bulks()                    # missing soft sphere? 
    plot_bulks_together()            # missing soft sphere? 

  
    
  
    
        








        
