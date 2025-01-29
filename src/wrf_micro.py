import numpy as np



#------------------------------------------------------------------------------
# WSM6 SNOW
def WSM6_get_psd_param_snow(q):

    rho_s = 100.;
    n0s = 2e6;
    lamdasmax = 1E5;
    Lambda = (np.pi*rho_s*n0s / (q))**(1./4.);
    if Lambda  > lamdasmax:
        Lambda = lamdasmax;
    N0 = n0s
    
    return [N0/1000., Lambda/1000.]

#------------------------------------------------------------------------------
# WSM6 GRAU
def WSM6_get_psd_param_grau(q):

    rho_g = 500.;
    n0g = 4e6;
    #rho = 0.8;
    lamdagmax = 6e4;
    Lambda = (np.pi*rho_g*n0g / (q))**(1./4.);
    if Lambda  > lamdagmax:
        Lambda = lamdagmax;
    N0 = n0g

    return [N0/1000., Lambda/1000.]

#------------------------------------------------------------------------------
# WSM6 ICE
def WSM6_get_psd_param_ice(q):

    rho_i = 900.;
    c     = 5.38E7;
    d     = 0.75;
    N     = c*((rho_i*q)**(d))

    return [N/1000.]
