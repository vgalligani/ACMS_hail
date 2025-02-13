import os 
import datetime

def config_folders(server): 

    folders ={} 
    
    if 'yakaira' in server: 
        # Leer wrfout y variables de interes: control-ndg
        #folder_WSM6  = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFout_WSM6_v4.5.2_all/'
        # re-run with correct namelist.input TWICE: uno con el correct ratio, y 2 con el correct domain: domain2
        
        # ---- WRFOUT files
        folders['WSM6_domain2']= '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFOUT/WSM6_domain2/'
        folders['WSM6_domain3']= '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFOUT/WSM6_domain3/'
        folders['WSM6_domain3_NoahMP']= '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFOUT/WSM6_domain3_NoahMP/'        
        folders['WSM6_domain3_NoahMP_10min']= '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFOUT/WSM6_domain3_NoahMP_10min'
        folders['WSM6_domain4_NoahMP']= '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFOUT/WSM6_domain4_NoahMP/'
        folders['P3_3MOM_LF_domain3_NoahMP'] = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFOUT/P3_3MOM_LF_domain3_NoahMP' 
        folders['P3_3MOM_LF_domain3_NoahMP_10min']='/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFOUT/P3_3MOM_LF_domain3_NoahMP_10min'
        folders['P3_3MOM_LF_domain3_NoahMP_highres']='/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/WRFOUT/P3_3MOM_LF_domain3_NoahMP_highres/'
        
        # ------ main savedir 
        # Within this folder, define the name of a sub-folder according to date
        folders['savedir'] = '/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/Plots/'
        folders['savedir']       = os.path.join(folders['savedir'],datetime.datetime.now().strftime('%d%m%Y'))
        if not os.path.exists(folders['savedir']):
            os.makedirs(folders['savedir'])

        folders['tempsavedir'] = '/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/Plots/tmp/'
        folders['output_gif_path'] = '/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/Plots/gifs/'
        
        folders['csapr2_dir'] = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/CSAPR2/'
        folders['rma1_dir'] = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RMA1/'
        folders['cswr_dir'] = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/CSWRdata/'
        
        folders['save_dir_compare'] = '/home/vito.galligani/datosmunin3/Work/Studies/HAILCASE_10112018/Plots/'
        folders['save_dir_compare'] = os.path.join(folders['save_dir_compare'],datetime.datetime.now().strftime('%d%m%Y'))


        csapr2_dir = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/CSAPR2/'
        rma1_dir = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/RMA1/'
        cswr_dir = '/home/vito.galligani/datosmunin3/Work/HAILCASE_10112018_datos/CSWRdata/'

        # Datos radar [CSAPR-2]
        folders['csapr2_str'] = [csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.170003.nc',
                      csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.174503.nc',
                      csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.194503.nc',
                      csapr2_dir+'corcsapr2cfrppiM1.a1.20181110.200003.nc',
                      csapr2_dir+'corcsapr2cfrppiqcM1.b1.20181110.201746.nc']
        
        # Datos de radar [RMA1]
        folders['rma1_str_02'] = [rma1_dir+'cfrad.20181110_170419.0000_to_20181110_170542.0000_RMA1_0301_02.nc',
                      rma1_dir+'cfrad.20181110_174533.0000_to_20181110_174657.0000_RMA1_0301_02.nc',
                      rma1_dir+'cfrad.20181110_194101.0000_to_20181110_194221.0000_RMA1_0301_02.nc',
                      rma1_dir+'cfrad.20181110_200546.0000_to_20181110_200709.0000_RMA1_0301_02.nc',
                      rma1_dir+'cfrad.20181110_202213.0000_to_20181110_202332.0000_RMA1_0301_02.nc']
        
        folders['rma1_str'] = [rma1_dir+'cfrad.20181110_170542.0000_to_20181110_171223.0000_RMA1_0301_01.nc',
                      rma1_dir+'cfrad.20181110_173842.0000_to_20181110_174523.0000_RMA1_0301_01.nc',
                      rma1_dir+'cfrad.20181110_194224.0000_to_20181110_194913.0000_RMA1_0301_01.nc',
                      rma1_dir+'cfrad.20181110_200709.0000_to_20181110_201358.0000_RMA1_0301_01.nc',
                      rma1_dir+'cfrad.20181110_202339.0000_to_20181110_203020.0000_RMA1_0301_01.nc']
        
        # Datos de radar [DOW6]
        folders['dow6_str'] = [cswr_dir+'cfrad.20181110_170011_DOW6high_v215_s01_el0.49_SUR.nc',
                      cswr_dir+'cfrad.20181110_174011_DOW6high_v227_s01_el0.49_SUR.nc',
                      cswr_dir+'cfrad.20181110_194512_DOW6high_v260_s02_el0.63_SUR.nc',    
                      cswr_dir+'cfrad.20181110_200014_DOW6high_v266_s02_el0.92_SUR.nc',
                      cswr_dir+'cfrad.20181110_202015_DOW6high_v274_s02_el0.90_SUR.nc'] # FALTA
        
        # Datos de radar [DOW7]
        folders['dow7_str'] = [cswr_dir+'cfrad.20181110_170011_DOW7high_v216_s01_el0.95_SUR.nc',
                    cswr_dir+'cfrad.20181110_172011_DOW7high_v222_s01_el0.97_SUR.nc',
                    cswr_dir+'cfrad.20181110_194054_DOW7high_v256_s01_el0.70_SUR.nc',  #LE LLUEVE ENCIMA! 
                    cswr_dir+'cfrad.20181110_200023_DOW7high_v264_s03_el1.14_SUR.nc',  #LE LLUEVE ENCIMA! 
                    cswr_dir+'cfrad.20181110_202012_DOW7high_v271_s02_el0.99_SUR.nc']  #JUSTO PASÓ
    
       # If the latter sub-folder does not exist, create it.
        if not os.path.exists(folders['save_dir_compare']):
            os.makedirs(folders['save_dir_compare'])        
        if not os.path.exists( os.path.join(folders['save_dir_compare'],'RadarSim')):
            os.makedirs(os.path.join(folders['save_dir_compare'],'RadarSim'))        
            

        if not os.path.exists( os.path.join(folders['save_dir_compare'],'OBS/RMA1')):
                os.makedirs(os.path.join(folders['save_dir_compare'],'OBS/RMA1'))    
        if not os.path.exists( os.path.join(folders['save_dir_compare'],'OBS/CSAPR2')):
                os.makedirs(os.path.join(folders['save_dir_compare'],'OBS/CSAPR2'))       
            
            
        EXPs = ['WSM6_domain2', 'WSM6_domain3', 'WSM6_domain3_NoahMP', 'WSM6_domain4_NoahMP', 'P3_3MOM_LF_domain3_NoahMP', 
                'P3_3MOM_LF_domain3_NoahMP_10min', 'WSM6_domain3_NoahMP_10min', 'P3_3MOM_LF_domain3_NoahMP_highres']
        for EXP in EXPs:    
            if not os.path.exists( os.path.join(folders['save_dir_compare'],EXP)):
                os.makedirs(os.path.join(folders['save_dir_compare'],EXP))                    
            if not os.path.exists( os.path.join(folders['save_dir_compare'],EXP+'/winds')):
                os.makedirs(os.path.join(folders['save_dir_compare'],EXP+'/winds'))      
            if not os.path.exists( os.path.join(folders['save_dir_compare'],EXP+'/convergence')):
                os.makedirs(os.path.join(folders['save_dir_compare'],EXP+'/convergence'))   
            if not os.path.exists( os.path.join(folders['save_dir_compare'],EXP+'/helicity')):
                os.makedirs(os.path.join(folders['save_dir_compare'],EXP+'/helicity'))                           

            
        # SimRad.AR data
        wdir = '/home/vito.galligani/datosmunin3/Work/SimRad.AR/data/'
        folders['LUT_WSM6'] = wdir+'UNIX_WSM6_LOOKUPTABLE_airmatrix_graupelAR1_ALLBANDS.pckl'
        
            
        
    else:
        # Leer wrfout y variables de interes: control-ndg
        print('CPU escritorio OFICINA')

    return folders 