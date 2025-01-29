#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------
@purpose : Analyse rttov 13.4 outputs  
@author  : Vito Galligani
@email   : vito.galligani@gmail.com
@condaenv: conda activate RTTOV14
-----------------------------------------------------------------
@main    : Functions to read and handle data files 
          
           Remember to open remote session:
               ssh -L 6000:yakaira:22 vito.galligani@portal.cima.fcen.uba.ar
               nohup python -m spyder_kernels.console -matplotlib='inline' -f=./remotemachine.json &
        
-----------------------------------------------------------------
"""
# Load libraries
import glob
import numpy as np
import os
from datetime import datetime
import pandas as pd
import Tools2Read as T2R


#-----------------------------------------------------------------------------
def LoadColAndResh(fname,nlev,nchan):
    """
    -------------------------------------------------------------
    Loads data stored column-wise by RTTOV-SCATT and reshapes it
    to (nlev,nchan)
    -------------------------------------------------------------
    OUT   var    Variable in nlev,nchan
    IN    fname  Variable in nlev x nchan
          nlev   Number of levels to read in. Based on ARTS
          nchan  Number of channels
    -------------------------------------------------------------
    """
    fname = sorted(glob.glob(str(fname)))
    print('paok',fname)
    var   = SingleCol(fname[0])
    var   = np.reshape(var,(nlev,nchan))
    return var

#-----------------------------------------------------------------------------
def read_scatt_aux(fnames,nlev,nchan):
    """
    -------------------------------------------------------------

    -------------------------------------------------------------

    -------------------------------------------------------------
    """
    out = np.empty((0,nlev,nchan))
    for i in fnames:
        d2 = LoadColAndResh(i,nlev,nchan)
        d2 = np.expand_dims(d2,0)
        out = np.append(out,d2,axis=0)
    out = np.array(out)
    return out


#-----------------------------------------------------------------------------
def SingleLayerTau_ltr(nlev,nchan,ds2):                         
    """
    -------------------------------------------------------------
    Loads total tau per surface zenith angle and returns the full
    dataset
    -------------------------------------------------------------
    OUT   ds_u   Layer tau at RTTOV levels
    IN    ds2    Datset for layer tau at RTTOV layers (files)
          nlev   The number of levels to read in.
          nchan  Number of channels
    -------------------------------------------------------------
    """
    out = np.empty((0,nchan,nlev))
    for i in ds2:
        d2 = SingleCol(i)
        d2 = np.reshape(d2,(nchan,nlev))
        d2 = np.expand_dims(d2,0)
        out = np.append(out,d2,axis=0)
    return out


#-----------------------------------------------------------------------------
def MultiCol(nlev,nchan,fdat):
    """
    -------------------------------------------------------------
    For all earth incident angles, loads and appends multi-column
    output from RTTOV. Example: transmission. It has nchan + 1
    columns. First column contains the number of model layers.
    Final row contains the # of channels.
    -------------------------------------------------------------
    OUT   out    Full dataset
    IN    fdat   List of filenames per sza
          nlev   The number of levels to read in
          nchan  Number of channels
    -------------------------------------------------------------
    """
    out = np.empty((0,nlev,nchan))
    for fi in fdat:
        with open(fi,"r") as f :
            lines   = f.readlines()
            data    = pd.DataFrame(l.strip().split() for l in lines)
            data    = np.array(data)
            data    = data.astype(float)
            data    = data[0:nlev,1:]    # From 2nd column
            data    = np.expand_dims(data,0)
            out     = np.append(out,data,axis=0)
            f.close()
    out = np.array(out)
    return out

#-----------------------------------------------------------------------------
def MultiRow(fdat):
    """
    -------------------------------------------------------------
    Loads multiple ascii files stored in a single row. Example:
    brightness temperature from the modified fortran scripts. Each
    file contains one row, one value per channel.
    -------------------------------------------------------------
    OUT   dset    Dataset
    IN    fdat    Filename
    -------------------------------------------------------------
    """
    dset = []
    for fi in fdat:
        dsi = np.genfromtxt(fi)
        dset.append(dsi)
    dset = np.array(dset)
    return dset

#-----------------------------------------------------------------------------
def SingleCol(fdat):
    """
    -------------------------------------------------------------
    Loads a single column ascii file. Note here that these files
    are stored as follows:

    Do ichan = 1, nchannels
      Do ilayer = 1, nlevels
      Enddo
    Enddo
    -------------------------------------------------------------
    OUT   dset    Dataset
    IN    fdat    Filename
    -------------------------------------------------------------
    """
    with open(fdat,"r") as f :
        lines   = f.readlines()
        dset    = pd.DataFrame(l.strip().split() for l in lines)
    dset        = np.array(dset)
    dset        = dset[:,0].astype(float)
    f.close()
    return dset


#-----------------------------------------------------------------------------
def read_outouts(ncfile, Aiprofile, nlev, nchan, sensor):
    """
    -------------------------------------------------------------
    Loads clear sky simulations
    -------------------------------------------------------------
    OUT   Tbs for both clearsky rttov13 and rttov14
    IN    ncfile    
          Aiprofile
    -------------------------------------------------------------
    """
    Tb = {}
    Trans = {}
    
    #----------------------------------------------------------------------
    # ---------------------- RTTOV V13.2 ----------------------------------
    rttovdir = '/home/vito.galligani/Work/RTTOV/rttov13.2/'
    rttov_flag = 'rttov132_Profile_'
    flag_name = rttov_flag +str( ncfile[-5:-3] )+'_'+str(Aiprofile)        

    # Path of rttov clear sky (cs) output
    prttov_cs  = rttovdir+'rttov_test/beta_test/'+sensor+'/TestClearSky_Nprof'+str(Aiprofile)+'/'
    prttov_as  = rttovdir+'rttov_test/beta_test/'+sensor+'/TestAllSky_'+str(Aiprofile)+'/'
    prttov_dummy_as  = rttovdir+'rttov_test/beta_test/'+sensor+'/TestAllSky_dummy_'+str(Aiprofile)+'/'
    # nlevels in rttov13 == 137 (full pressure levels)
    nlev_rttov = 137
    
    #----------------------
    frttov_tb  = sorted(glob.glob(str(prttov_cs)+'output_tb_*_atm_em1.0'+str(flag_name)+'.dat'))
    rttov_tb   = MultiRow(frttov_tb)
    Tb['rttov13'] = rttov_tb
    print(rttov_tb.shape)    

    #    
    frttov_tb      = sorted(glob.glob(str(prttov_dummy_as)+'output_tb_dummy_*_atm_em1.0'+str(flag_name)+'.dat'))
    Tb['rttov13_dummyscatt'] = MultiRow(frttov_tb)
    
    frttov_tb      = sorted(glob.glob(str(prttov_as)+'output_*atm_em1.0'+str(flag_name)+'.dat'))
    Tb['rttov13_scatt'] = MultiRow(frttov_tb)
    

    #----------------------
    # 2/ transmittance: level to space transmittance in main example_fw (CLEARSKY) 
    # Level to space transmittances: output_transm >> output_taulevel 
    frttov_tr  = sorted(glob.glob(str(prttov_cs)+'output_taulevel_*_atm_em1.0'+str(flag_name)+'.dat'))
    Trans['rttov13_cs_level2spaceTr'] = MultiCol(nlev_rttov, nchan, frttov_tr)

    # 3/ transmittance: Surface to space transmittance in main example_fw (CLEARSKY) 
    # Surface to space transmittances: surface2space_transm >> output_tautotal_$
    frttov_tr  = sorted(glob.glob(str(prttov_cs)+'output_tautotal_*_atm_em1.0'+str(flag_name)+'.dat'))
    Trans['rttov13_cs_Surface2spaceTr'] = MultiRow(frttov_tr)
    print(Trans['rttov13_cs_Surface2spaceTr'].shape)

    # 4/ layer OD (full levels)
    # Here interested only in nadir because in rttov-scatt scatt_aux%ext inputs the nadir only! 
    # in RTTOV-SCATT full levels are taken to nlev+1 so that od (half pressure levels == nlev)
    # This is in half pressure levels (nlevels+1)
    # nlevels == full pressure leves input into rttov clear. 
    # nlevels + 1 
    # od_rttov(1) == TOP RTTOV Level to space
    # od_rttov(2:nlevels) 
    # od_rttov(nlevels+1) == surface to bottom RTTOV (full pressure) level
    frttov_ltr = SingleCol(prttov_as+'Prop/'+'SingleLayerOpticalDepthsRTTOVLevels_cs_180.0.txt')
    Trans['rttov13_cs_layerOD_nlev+1_fullp']  = np.reshape(frttov_ltr, (nchan, nlev_rttov+1))   
    
    #---- check also this in rttov_transmitt
    # od_level in rttov_transmit [nlayers] == nlevels-1 == 136
    # where od_level == trans_aux_path%od_singlelayer(:,lay,j) = opdp_path(lay+1,j) - opdp_path(lay,j)
    # od_level(:,j) = MAX(-opdp_path(:,j), -max_optical_depth)
    ODnlayer_rttov13 = SingleCol(str(prttov_cs)+'Prop/Level2SpaceTauUserLevels_CS_minus1_  0.0.txt')
    ODnlayer_rttov13 = np.reshape(ODnlayer_rttov13, (nchan, nlev_rttov-1))
    check_out = np.empty(ODnlayer_rttov13.shape)
    for i in range( check_out.shape[1]-1 ):
        if i==0:
            check_out[:,i] = 0 
        else:
            check_out[:,i] = -1*(ODnlayer_rttov13[:,i] - ODnlayer_rttov13[:,i-1]) 
    Trans['rttov13_cs_tauUser'] = check_out   # rttov14_layerOD? 
    #----
    
    # 5/ layer OD (RTTOV-SCATT half-levels)
    # rttov_ltu (od in rttov_iniscatt.F90) is the od used in rttov-scatt!! nlevels 
    # this is between half pressure levels    
    # This one es el que tiene el efecto de HE y dz
    rttov_ltu  = SingleCol(prttov_as+'Prop/'+'SingleLayerOpticalDepthsOurLevels_cs_180.0.txt')
    Trans['rttov13_as_layerOD_nlev_halfp'] = np.reshape(rttov_ltu, (nchan, nlev_rttov))   # This is 're-allocation in 1/2 levels'! nlevels 
    
    # NOTE TO VITO: necesito entender el cambio de Trans['rttov13_cs_tauUser'] in rttov_transmitt a 
    # Trans['rttov13_cs_layerOD_nlev+1_fullp'] y Trans['rttov13_as_layerOD_nlev_halfp'] inside rttov_iniscatt! 
    # the differnce si the way the surface is handled! surface OD !! see that the surface OD is higher for rttov-scatt 
    
    # 6/ transmittance: Surface to space transmittance in main example_fw ((RTTOV-SCATT)) 
    # Total Optical depth for RTTOV-SCATT (after OD re-allocation to half levels and HE)
    frttov_tr = sorted(glob.glob(prttov_as+'Prop/'+'FullTransmission_BeEdd_*.txt'))[::-1]                      
    Trans['rttov13_as_Surface2spaceTr'] =  MultiRow(frttov_tr)         
    print(prttov_as+'atm_em1.0'+flag_name+'_Prop/'+'FullTransmission_BeEdd_*.txt')
    
    # 7/ nadir OD derivation issue in RTTOV-13 rttov-scatt
    # Ext from mw/rttov_iniscatt before hydro-mieprox-iniedd. Input nadir OD to rttov_scatt
    fextPro_angles = sorted(glob.glob(prttov_as+'Prop/'+'ExtProHydro_*'))[::-1]  
    Trans['rttov13_as_extPro_angles'] = read_scatt_aux(fextPro_angles,nlev_rttov,nchan)
    print(Trans['rttov13_as_extPro_angles'].shape)    
    print(prttov_as+'Prop/'+'ExtProHydro_*')
    print('-----------------------')
    



    #----------------------------------------------------------------------
    # ---------------------- RTTOV V14.0 ----------------------------------
    rttovdir = '/home/vito.galligani/Work/RTTOV/rttov14.0_beta/'
    rttov_flag = 'rttov14_Profile_' 
    flag_name = rttov_flag +str( ncfile[-5:-3] )+'_'+str(Aiprofile)        
    # nlevels in rttov14 == 138 (half pressure levels)
    nlev_rttov = 138
    
    # Path of rttov clear sky (cs) output FOR OLD RTTOV13 COEFFS
    prttov_cs  = rttovdir+'rttov_test/beta_test/'+sensor+'/TestClearSky_Nprof'+str(Aiprofile)+'/'
    prttov_as  = rttovdir+'rttov_test/beta_test/'+sensor+'/TestAllSky_Nprof'+str(Aiprofile)+'/'
    
    #----------------------
    # 1/ Tbs
    frttov_tb  = sorted(glob.glob(str(prttov_cs)+'output_tb_*_atm_em1.0'+str(flag_name)+'.dat'))
    rttov_tb   =  MultiRow(frttov_tb)    
    Tb['rttov14'] = rttov_tb
    print(rttov_tb.shape)    

    frttov_tb      = sorted(glob.glob(str(prttov_as)+'output_tb_*_atm_em1.0'+str(flag_name)+'.dat_optsHydroScatt_TRUE'))
    Tb['rttov14_dummyscattTRUE'] =  MultiRow(frttov_tb)
    print(Tb['rttov14_dummyscattTRUE'].shape)    
    
    frttov_tb      = sorted(glob.glob(str(prttov_as)+'output_tb_*_atm_em1.0'+str(flag_name)+'.dat_optsHydroScatt_FALSE'))
    Tb['rttov14_dummyscattFALSE'] =  MultiRow(frttov_tb)    
    print(Tb['rttov14_dummyscattFALSE'].shape)    

    frttov_tb      = sorted(glob.glob(str(prttov_as)+'output_tb_*_atm_em1.0'+str(flag_name)+'.dat'))
    Tb['rttov14_scatt'] =  MultiRow(frttov_tb)    
    print(Tb['rttov14_scatt'].shape)    
    
    
    #----------------------
    # 2/ transmittance: level to space transmittance in main example_fw (CLEARSKY) 
    # Level to space transmittances: output_transm >> output_taulevel 
    frttov_tr  = sorted(glob.glob(str(prttov_cs)+'output_taulevel_*_atm_em1.0'+str(flag_name)+'.dat'))
    Trans['rttov14_cs_level2spaceTr'] =  MultiCol(nlev_rttov, nchan, frttov_tr)
   
    # 3/ transmittance: Surface to space transmittance in main example_fw (CLEARSKY) 
    # Surface to space transmittances: surface2space_transm >> output_tautotal_$
    frttov_tr  = sorted(glob.glob(str(prttov_cs)+'output_tautotal_*_atm_em1.0'+str(flag_name)+'.dat'))
    Trans['rttov14_cs_Surface2spaceTr'] =  MultiRow(frttov_tr)

    frttov_tr  = sorted(glob.glob(str(prttov_as)+'output_tautotal_*_atm_em1.0'+str(flag_name)+'.dat_optsHydroScatt_FALSE'))
    Trans['rttov14_as_Surface2spaceTr_FALSE'] = MultiRow(frttov_tr)
    print(str(prttov_as)+'output_tautotal_atm_*_em1.0'+str(flag_name)+'.dat_optsHydroScatt_FALSE')

    frttov_tr  = sorted(glob.glob(str(prttov_as)+'output_tautotal_*_atm_em1.0'+str(flag_name)+'.dat_optsHydroScatt_TRUE'))
    Trans['rttov14_as_Surface2spaceTr_TRUE'] =  MultiRow(frttov_tr)
    
    #----------------------
    # od_level in rttov_transmit_thermal [nlayers] == nlevels-1 == 137
    # where od_level == trans_aux_path%od_singlelayer(:,lay,j) = opdp_path(lay+1,j) - opdp_path(lay,j)
    # od_level(:,j) = MAX(-opdp_path(:,j), -max_optical_depth)
    # nlayers == FULL LEVELS and nlevels == HALF PRES. LEVELS
    ODnlayer_rttov14 =  SingleCol(str(prttov_cs)+'Prop/Level2SpaceTauUserLevels_CS_minus1_0.0.txt')
    ODnlayer_rttov14 = np.reshape(ODnlayer_rttov14, (nchan, nlev_rttov-1))

    check_out = np.empty(ODnlayer_rttov14.shape)
    for i in range( check_out.shape[1]-1 ):
        if i==0:
            check_out[:,i] = 0 
        else:
            check_out[:,i] = -1*(ODnlayer_rttov14[:,i] - ODnlayer_rttov14[:,i-1]) 
    Trans['rttov14_cs_tauUser'] = check_out   # rttov14_layerOD? 
         
    
    # Test nadir angles and hydro == 0 extension scatt_aux % ext_all(1,:,ichan) 
    # a) how does it compare with the clear sky stuff? 
    # SCATT_AUX%EXT_ALL == scatt_aux % ext_all(0,:,ichan) + ext_clear(:)
    #---------------------------------------------------------------------------------------------
    path = str(prttov_as)+'Prop_optsHydroScatt_TRUE/'      
    
    hydroscatt_f = sorted(glob.glob(path+'scatt_aux_extHydro_*'))[::-1]      # ---> SCATT_AUX%EXT_ALL
    Trans['rttov14_aux_extAll'] =  read_scatt_aux(hydroscatt_f,nchan,nlev_rttov-1)
    print( Trans['rttov14_aux_extAll'].shape )
    print(path+'scatt_aux_extHydro_*') 
    
    scattfile = sorted(glob.glob(path+'scatt_aux_extLayerOD_*'))[::-1]      # ---> SCATT_AUX%EXT_ALL
    Trans['rttov14_aux_extLayerOD'] =  read_scatt_aux(scattfile,nchan,nlev_rttov-1)

    scattfile = sorted(glob.glob(path+'scatt_aux_extClear_*'))[::-1]      # ---> SCATT_AUX%EXT_ALL
    Trans['rttov14_aux_extClear'] =  read_scatt_aux(scattfile,nchan,nlev_rttov-1)


    #---------------------------------------------------------------------------------------------
    # Path of rttov clear sky (cs) output rttov_coefs1
    prttov_cs  = rttovdir+'rttov_test/beta_test/'+sensor+'/TestClearSky_Nprof'+str(Aiprofile)+'_option1/'
    #----------------------
    # 1/ Tbs
    frttov_tb  = sorted(glob.glob(str(prttov_cs)+'output_tb_*_atm_em1.0'+str(flag_name)+'.dat'))
    Tb['rttov14_option1'] =  MultiRow(frttov_tb)    

    #----------------------
    # 2/ transmittance: level to space transmittance in main example_fw (CLEARSKY) 
    # Level to space transmittances: output_transm >> output_taulevel 
    frttov_tr  = sorted(glob.glob(str(prttov_cs)+'output_taulevel_*_atm_em1.0'+str(flag_name)+'.dat'))
    Trans['rttov14_cs_option1_level2spaceTr'] =  MultiCol(nlev_rttov, nchan, frttov_tr)
   
    # 3/ transmittance: Surface to space transmittance in main example_fw (CLEARSKY) 
    # Surface to space transmittances: surface2space_transm >> output_tautotal_$
    frttov_tr  = sorted(glob.glob(str(prttov_cs)+'output_tautotal_*_atm_em1.0'+str(flag_name)+'.dat'))
    Trans['rttov14_cs_option1_Surface2spaceTr'] =  MultiRow(frttov_tr)

    #----------------------
    # od_level in rttov_transmit_thermal [nlayers] == nlevels-1 == 137
    # where od_level == trans_aux_path%od_singlelayer(:,lay,j) = opdp_path(lay+1,j) - opdp_path(lay,j)
    # od_level(:,j) = MAX(-opdp_path(:,j), -max_optical_depth)
    # nlayers == FULL LEVELS and nlevels == HALF PRES. LEVELS
    ODnlayer_rttov14 =  SingleCol(str(prttov_cs)+'Prop/Level2SpaceTauUserLevels_CS_minus1_0.0.txt')
    ODnlayer_rttov14 = np.reshape(ODnlayer_rttov14, (nchan, nlev_rttov-1))

    check_out = np.empty(ODnlayer_rttov14.shape)
    check_out[:] = np.nan
    for i in range( check_out.shape[1]-1 ):
        if i==0:
            check_out[:,i] = 0 
        else:
            check_out[:,i] = -1*(ODnlayer_rttov14[:,i] - ODnlayer_rttov14[:,i-1]) 
    Trans['rttov14_cs_tauUser_option1'] = check_out   # rttov14_layerOD? 

    #---------------------------------------------------------------------------------------------
    # Path of rttov clear sky (cs) output rttov_coefs1
    prttov_cs  = rttovdir+'rttov_test/beta_test/'+sensor+'/TestClearSky_Nprof'+str(Aiprofile)+'_option2b/'

    #----------------------
    # 1/ Tbs
    frttov_tb  = sorted(glob.glob(str(prttov_cs)+'output_tb_*_atm_em1.0'+str(flag_name)+'.dat'))
    Tb['rttov14_option2b'] =  MultiRow(frttov_tb)    

    #----------------------
    # 2/ transmittance: level to space transmittance in main example_fw (CLEARSKY) 
    # Level to space transmittances: output_transm >> output_taulevel 
    frttov_tr  = sorted(glob.glob(str(prttov_cs)+'output_taulevel_*_atm_em1.0'+str(flag_name)+'.dat'))
    Trans['rttov14_cs_option2b_level2spaceTr'] =  MultiCol(nlev_rttov, nchan, frttov_tr)
   
    # 3/ transmittance: Surface to space transmittance in main example_fw (CLEARSKY) 
    # Surface to space transmittances: surface2space_transm >> output_tautotal_$
    frttov_tr  = sorted(glob.glob(str(prttov_cs)+'output_tautotal_*_atm_em1.0'+str(flag_name)+'.dat'))
    Trans['rttov14_cs_option2b_Surface2spaceTr'] =  MultiRow(frttov_tr)

    #----------------------
    # od_level in rttov_transmit_thermal [nlayers] == nlevels-1 == 137
    # where od_level == trans_aux_path%od_singlelayer(:,lay,j) = opdp_path(lay+1,j) - opdp_path(lay,j)
    # od_level(:,j) = MAX(-opdp_path(:,j), -max_optical_depth)
    # nlayers == FULL LEVELS and nlevels == HALF PRES. LEVELS
    ODnlayer_rttov14 =  SingleCol(str(prttov_cs)+'Prop/Level2SpaceTauUserLevels_CS_minus1_0.0.txt')
    ODnlayer_rttov14 = np.reshape(ODnlayer_rttov14, (nchan, nlev_rttov-1))

    check_out = np.empty(ODnlayer_rttov14.shape)
    for i in range( check_out.shape[1]-1 ):
        if i==0:
            check_out[:,i] = 0 
        else:
            check_out[:,i] = -1*(ODnlayer_rttov14[:,i] - ODnlayer_rttov14[:,i-1]) 
    Trans['rttov14_cs_tauUser_option2b'] = check_out   # rttov14_layerOD? 
    
    return Tb, Trans

#-----------------------------------------------------------------------------
def read_outouts_SCATT(ncfile, Aiprofile, nlev, nchan, sensor):
    """
    -------------------------------------------------------------
    Loads clear sky simulations
    -------------------------------------------------------------
    OUT   Tbs for both clearsky rttov13 and rttov14
    IN    ncfile    
          Aiprofile
    -------------------------------------------------------------
    """
    Tb  = {}
    dTb = {}
    SCATT = {}
    HYDROSCATT = {}
    
    #----------------------------------------------------------------------
    # ---------------------- RTTOV V13.2 ----------------------------------
    rttovdir = '/home/vito.galligani/Work/RTTOV/rttov13.2/'
    rttov_flag = 'rttov132_Profile_'
    flag_name = rttov_flag +str( ncfile[-5:-3] )+'_'+str(Aiprofile)        

    # Path of rttov clear sky (cs) output
    prttov_cs  = rttovdir+'rttov_test/beta_test/'+sensor+'/TestClearSky_Nprof'+str(Aiprofile)+'/'
    prttov_as  = rttovdir+'rttov_test/beta_test/'+sensor+'/TestAllSky_'+str(Aiprofile)+'/'
        
    # nlevels in rttov13 == 137 (full pressure levels)
    nlev_rttov = 137
        
    #----------------------
    frttov_tb  = sorted(glob.glob(str(prttov_cs)+'output_tb_*_atm_em1.0'+str(flag_name)+'.dat'))
    rttov_tb   =  MultiRow(frttov_tb)
    Tb['rttov13'] = rttov_tb
    
    frttov_tb      = sorted(glob.glob(str(prttov_as)+'output_*atm_em1.0'+str(flag_name)+'.dat'))
    Tb['rttov13_scatt'] =  MultiRow(frttov_tb)
    print(Tb['rttov13_scatt'].shape)
    print(str(prttov_as)+'output_*atm_em1.0'+str(flag_name)+'.dat')
    
    dTb['rttov13'] = Tb['rttov13_scatt'] - Tb['rttov13']
    
    #----------------------  in mieproc
    # - Read RTTOV-SCATT SSA (SCATT-AUX%SSA)
    SCATT['rttov13_ssa']  = np.reshape( SingleCol(prttov_as+'Prop/'+'ssa_full_180.0.txt'),(nchan,nlev_rttov))
    SCATT['rttov13_ext']  = np.reshape( SingleCol(prttov_as+'Prop/'+'ext_full_180.0.txt'),(nchan,nlev_rttov))
    SCATT['rttov13_asm']  = np.reshape( SingleCol(prttov_as+'Prop/'+'asm_full_180.0.txt'),(nchan,nlev_rttov))

    # OJO ACA ES LAYER, NCHAN, HYDRO ...
    HYDROSCATT['rttov13_ssa']  = np.reshape( SingleCol(prttov_as+'Prop/'+'itypessa_angle_180.0.txt'),(nlev_rttov,nchan,5))
    HYDROSCATT['rttov13_ext']  = np.reshape( SingleCol(prttov_as+'Prop/'+'itypeext_angle_180.0.txt'),(nlev_rttov,nchan,5))
    HYDROSCATT['rttov13_asm']  = np.reshape( SingleCol(prttov_as+'Prop/'+'itypeasm_angle_180.0.txt'),(nlev_rttov,nchan,5))

    #----------------------------------------------------------------------
    # This comes from rttov_eddington_setup.F90 
    # scatt_aux % ext_all, asm_all, aam_all. ext hydro + gas 
    # SCATT_AUXS HYDRO IS AN INPUT >> 
    # a) rttov_optp_interpolate_setup
    # b) rttov_optp_interpolate
    # c) rttov_scatt_optp
    # ---------------------- RTTOV V14.0 ----------------------------------
    rttovdir = '/home/vito.galligani/Work/RTTOV/rttov14.0_beta/'
    rttov_flag = 'rttov14_Profile_' 
    flag_name = rttov_flag +str( ncfile[-5:-3] )+'_'+str(Aiprofile)        
    # nlevels in rttov14 == 138 (half pressure levels)
    nlev_rttov = 138
    
    # Path of rttov clear sky (cs) output FOR OLD RTTOV13 COEFFS
    prttov_cs  = rttovdir+'rttov_test/beta_test/'+sensor+'/TestClearSky_Nprof'+str(Aiprofile)+'/'
    prttov_as  = rttovdir+'rttov_test/beta_test/'+sensor+'/TestAllSky_Nprof'+str(Aiprofile)+'/'
    
    #----------------------
    # 1/ Tbs
    frttov_tb  = sorted(glob.glob(str(prttov_cs)+'output_tb_*_atm_em1.0'+str(flag_name)+'.dat'))
    rttov_tb   =  MultiRow(frttov_tb)    
    Tb['rttov14'] = rttov_tb

    frttov_tb      = sorted(glob.glob(str(prttov_as)+'output_tb_*_atm_em1.0'+str(flag_name)+'.dat_SCATT'))
    Tb['rttov14_scatt'] =  MultiRow(frttov_tb)    
    
    dTb['rttov14'] = Tb['rttov14_scatt'] - Tb['rttov14']

    #----------------------  in mieproc
    # - Read RTTOV-SCATT SSA (SCATT-AUX%SSA)
    SCATT['rttov14_ssa']  = np.reshape( SingleCol(prttov_as+'Prop_SCATT/'+'scatt_aux_ssaHydro_0.0.txt'),(nchan,nlev_rttov-1))
    SCATT['rttov14_ext']  = np.reshape( SingleCol(prttov_as+'Prop_SCATT/'+'scatt_aux_extHydro_0.0.txt'),(nchan,nlev_rttov-1))
    SCATT['rttov14_asm']  = np.reshape( SingleCol(prttov_as+'Prop_SCATT/'+'scatt_aux_asmHydro_0.0.txt'),(nchan,nlev_rttov-1))


    HYDROSCATT['rttov14_ssa']  = np.reshape( SingleCol(prttov_as+'Prop_SCATT/'+'optp_hydro_ssa_0.0.txt'),(nchan,nlev_rttov-1,5))
    HYDROSCATT['rttov14_ext']  = np.reshape( SingleCol(prttov_as+'Prop_SCATT/'+'optp_hydro_ext_0.0.txt'),(nchan,nlev_rttov-1,5))
    HYDROSCATT['rttov14_asm']  = np.reshape( SingleCol(prttov_as+'Prop_SCATT/'+'optp_hydro_asy_0.0.txt'),(nchan,nlev_rttov-1,5))



    return Tb, dTb, SCATT, HYDROSCATT



