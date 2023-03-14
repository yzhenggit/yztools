def find_logmstar_const_munshi21(do_print=False):
    import numpy as np 
    # Munshi+2021's SMHM is a broken power law
    # but they didn't specify the reference stellar mass
    # so here we are going to calculate it by making sure that at the high mass end
    # the mass is the same at the Behroozi+2013 relation 
    # and the constant is the same as used in Munshi+2021's github: 
    # https://github.com/emapple/smhm-toy-model/blob/main/toy-model.ipynb
    
    logm1 = 10 # characteristic halo mass in Munshi+2021 
    alpha_high = 1.9
    
    # use logmstar_b13's value at logmhalo = 1e11.514 to anchor the constant in Munshi since they don't give one 
    from yztools.mhalo2mstar import mhalo2mstar_behroozi13
    logmhalo_b13_M1 = 11.514 # characteristic M1 peak mass in Behroozi+2013
    logmstar_b13_ref = mhalo2mstar_behroozi13(logmhalo_b13_M1)[0] # get the corresponding stellar mass 

    # use that smhm at b13 at the ref mhalo to get a constant 
    logmstar_const_m21 = logmstar_b13_ref - alpha_high*np.log10(10**logmhalo_b13_M1 / 10**logm1)
    logmstar_const_m21 = np.around(logmstar_const_m21, decimals=2)
    if do_print == True: 
        print('Munshi+2021, logmstar_const = {:.2f}'.format(logmstar_const_m21))
        print("For Mstar = mstar_const * (Mhalo/M1)^alpha")
        
    return logmstar_const_m21
    
def m200c2mstar_munshi21(logmhalo, do_print=False, return_sigma=False):
    '''
    add Munshi+2021's SMHM based on parameters shown in their Figure 2 
    designed from lower mass dwarfs from LMC-mass to UFDs from mhalo=1e7 to mhalo = 1e11 Msun
    or mstar < 1e4 to 8.2e8 Msun

    Note that the Munshi+2021's SMHM relation is based on 200c definition.
    History: Yong Zheng, Dec 2022, RPI
    '''
    import numpy as np
    logmhalo_arr = np.asarray([logmhalo]).flatten()
    logmstar_arr = np.zeros(len(logmhalo_arr))
    logmstar_sigma_arr = np.zeros(len(logmhalo_arr))
    
    from yztools.mhalo2mstar import mhalo2mstar_behroozi13
    logmhalo_b13_M1 = 11.514 # characteristic M1 peak mass in Behroozi+2013
    
    # if there are values in logmhalo_arr higher than logmhalo_b13_M1, we use B13's relation 
    ind_logmhalo_high = logmhalo_arr >= logmhalo_b13_M1
    if len(logmhalo_arr[ind_logmhalo_high]) >0: 
        if do_print==True: 
            print("There are values higher than M1={}, gonna use Behroozi+2013's SMHM for those".format(logmhalo_b13_M1))
        logmstar_arr[ind_logmhalo_high] = mhalo2mstar_behroozi13(logmhalo_arr[ind_logmhalo_high])
        logmstar_sigma_arr[ind_logmhalo_high] = 0.3 # Munshi+2021 , page 4
    
    # for the rest of the values, use Munshi+2021's relation 
    logmstar_const_m21 = find_logmstar_const_munshi21()
    
    # fitted values from Figure 2 of Munshi+2021 
    sigma0 = 0.3
    gamma = -0.39

    alpha_high = 1.9
    alpha_low = 2.8
    
    logm1 = 10 # characteristic halo mass in Munshi+2021 

    # now generate the SMHM relation for Munshi+2021, only use it for low mass halo 
    ind_high = np.all([logmhalo_arr >= logm1, logmhalo_arr<logmhalo_b13_M1], axis=0)
    logmstar_arr[ind_high] = logmstar_const_m21 + alpha_high*np.log10(10**logmhalo_arr[ind_high] / 10**logm1)
    logmstar_sigma_arr[ind_high] = sigma0

    ind_low = logmhalo_arr<= logm1
    logmstar_arr[ind_low] = logmstar_const_m21 + alpha_low*np.log10(10**logmhalo_arr[ind_low] / 10**logm1)
    logmstar_sigma_arr[ind_low] = sigma0+gamma*(logmhalo_arr[ind_low] - logm1)
   
    if do_print == True:
        print(">> Note that Munshi+2021 derivation is for M*=4-8.2e8 Msunish, or Mhalo=1e7-1e11\n")
        print(">> At Mhalo>1e11.514, we use Behroozi+2013 here.")
    
    if return_sigma == True: 
        return logmstar_arr, logmstar_sigma_arr
    else: 
        return logmstar_arr

def m200c2mstar_moster13(log_mhalo, do_print=False):
    """
    Moster+2013's halo abundance matching.
    Based on Chabrier IMF
    Calculated over stellar mass of logM*=7.4--11.7 Msun.
    All virial masses are based on 200c 

    History:
    YZ, Dec 2022, RPI, update notes to m200c from mhalo
    YZ, 05/06/2022, UCB
    """

    import astropy.constants as const
    from astropy.cosmology import Planck15
    import astropy.units as u
    import numpy as np

    # from Table 1 in Moster+2013, 0 mean z=0
    log_M10 = 11.59
    M10 = 10**log_M10 # Msun,
    sig_logM10 = 0.236
    sig_M10 = sig_logM10*M10*np.log(10)

    N10 = 0.0351
    sig_N10 = 0.0058

    beta10 = 1.376
    sig_beta10 = 0.153

    gamma10 = 0.608
    sig_gamma10 = 0.059

    # from moster equation 2
    mhalo = 10**log_mhalo
    mstar = mhalo * (2.*N10*((mhalo/M10)**(-beta10)+(mhalo/M10)**(gamma10))**(-1.))

    if do_print == True:
        print(">> Note that Moster+2013 derivation is for logM*=7.4-11.7, Chabrier IMF")
    return np.log10(mstar)

def mhalo2mstar_moster10(log_mhalo, do_print=False):
    """
    Moster+2010's halo abundance matching.
    Based on Kroupa 2001 IMF
    Calculated over stellar mass of logM*=8.5--11.85 Msun.

    YZ. 05/10/2022, UCB.
    """
    import numpy as np

    log_mhalo = np.asarray([log_mhalo]).flatten()

    # from Table 6 in Moster+2010
    logM1 = 11.88
    sig_logM1 = 0.02
    M1 = 10**logM1 # Msun
    sig_M1 = sig_logM1*M1*np.log(10)

    mM0 = 0.0282
    sig_mM0 = 0.0005

    beta = 1.06
    sig_beta = 0.05

    gamma = 0.56
    sig_gamma = 0.0

    # from moster equation 2
    mhalo = 10**log_mhalo
    mstar = mhalo * (2.*mM0*((mhalo/M1)**(-beta)+(mhalo/M1)**(gamma))**(-1.))

    if do_print == True:
        print(">> Note that Moster+2010 derivation is for logM*=8.5-11.85, Kroupa IMF\n")
    return np.log10(mstar)

def b13_fx(x, alpha, delta, gamma):
    """
    Eq 3 in Behroozi+2013, where x=log10(Mhalo/M1), and M1 is the characteristic
    halo mass, as shown in Section 5.

    YZ. 05/10/2022. UCB
    """
    #
    # f(x)
    import numpy as np
    parta = -np.log10(10**(alpha*x)+1)
    partb = delta*(np.log10(1+np.exp(x)))**gamma / (1+np.exp(10**-x))

    return parta + partb

def mhalo2mstar_behroozi13(log_mhalo, do_print=False):
    """
    halo mass to stellar mass fuction as shown in Eq. 3 by Behroozi+2013
    Based on Chabrier IMF, which is similar to Kroupa IMF.
    Fit for logM*=7.25-11.85.

    YZ. 05/10/2022. UCB
    """
    import numpy as np

    log_mhalo = np.asarray([log_mhalo]).flatten()

    Mhalo = 10**log_mhalo

    # intrinsic parameter values as listed in Sec. 5 in Behroozi+2013
    log10_epsilon0 = -1.777
    epsilon0 = 10**log10_epsilon0

    log10_M10 = 11.514
    M10 = 10**log10_M10

    alpha0 = -1.412
    delta0 = 3.508
    gamma0 = 0.316

    # Eq 3
    x = np.log10(Mhalo/M10)
    log10_Mstar = np.log10(epsilon0*M10) + b13_fx(x, alpha0, delta0, gamma0) - b13_fx(0, alpha0, delta0, gamma0)

    if do_print == True:
        print(">> Note that Behroozi+2013 derivation is for logM*=7.25-11.85, Chabrier IMF\n")
    return log10_Mstar

def mhalo2mstar_behroozi13_gk14modified(log_mhalo, do_print=False, threshold_MhaloM10=0.1):
    """
    halo mass to stellar mass fuction as shown in Eq. 3 by Behroozi+2013
    Based on Chabrier IMF, which is similar to Kroupa IMF.
    Fit for logM*=7.25-11.85.

    The lower mass end (default Mhalo/M1<=0.1) is based on Garrison-Kimmel+2014, which found that
    a steeper slope (alpha=1.92 as compared to alpha=1.412) is better at reproducing the
    dwarf galaxy counts in the Local Group.

    YZ. 05/10/2022. UCB
    """

    import numpy as np
    log_mhalo = np.asarray([log_mhalo]).flatten()

    Mhalo = 10**log_mhalo

    # intrinsic parameter values as listed in Sec. 5 in Behroozi+2013
    log10_epsilon0 = -1.777
    epsilon0 = 10**log10_epsilon0

    log10_M10 = 11.514
    M10 = 10**log10_M10

    alpha0 = -1.412
    alpha0_lowmass = -1.92 # as suggested by Garrison-Kimmel-2014 for the lower mass end

    delta0 = 3.508
    gamma0 = 0.316

    # Eq 3
    x = np.log10(Mhalo/M10)
    ind_lowmass = x<threshold_MhaloM10
    ind_highmass = np.logical_not(ind_lowmass)

    # for low and high mass separately
    log10_Mstar = np.zeros(log_mhalo.size)
    log10_Mstar[ind_lowmass] = np.log10(epsilon0*M10) + b13_fx(x[ind_lowmass], alpha0_lowmass, delta0, gamma0) - \
                               b13_fx(0, alpha0_lowmass, delta0, gamma0)
    log10_Mstar[ind_highmass]= np.log10(epsilon0*M10) + b13_fx(x[ind_highmass], alpha0, delta0, gamma0) - \
                               b13_fx(0, alpha0, delta0, gamma0)

    if do_print == True:
        print(">> Note that Behroozi+2013 derivation is for logM*=7.25-11.85, Chabrier IMF\n")
        print(">> Lower mass end (Mhalo/M1<={}) use a steeper slope (alpha=1.92) based on Garrison-Kimmel+2014".format(threshold_MhaloM10))
    return log10_Mstar

def gk17_alpha_scatter(scatter):
    """
    Based on Garrison-Kimmel+2017, where alpha is a function of scatter in the SMHM relation

    YZ. 5/10/2022, UCB.
    """
    alpha = 0.24*scatter**2 + 0.16*scatter + 1.99
    return alpha

def mhalo2mstar_gk17(log_mhalo, scatter=0.2, do_print=False):
    """
    Based on Barrison-Kimmel+2017. The power lawer is anchored at M1=11.75 (see Behroozi+2013)

    YZ. 5/10/2022, UCB.
    """
    import numpy as np

    log_mhalo = np.asarray([log_mhalo]).flatten()
    Mhalo = 10**log_mhalo

    log10_M10 = 11.514
    M10 = 10**log10_M10

    # assuming constant scatter, field dwarfs, Eq 6 in Garrison-Kimmel+2017
    alpha_LF = gk17_alpha_scatter(scatter)

    x = np.log10(Mhalo/M10)
    # log10_Mstar = np.log10(epsilon0*M10) + b13_fx(x, alpha0, delta0, gamma0) - b13_fx(0, alpha0, delta0, gamma0)
    # 9.392 = np.log10(epsilon0*M10) - b13_fx(0, alpha0, delta0, gamma0)
    # when x=0, doesn't matter what input alpha is in b13_fx
    log10_Mstar = 9.392 + alpha_LF*x

    if do_print == True:
        print(">> Note that GK17 derivation is for logM*=6-8, Kroupa IMF in FIRE\n")
        print(">> this relation matches with Behroozi+2013 at the characteristic mass M1. ")

    return log10_Mstar

def mhalo2mstar_gk14(log_mhalo, do_print=False):
    """Garrison-Kimmel+2014a halo abundance matching for dwarfs logM*=6-8 """
    import numpy as np

    ## from Garrison-Kimmel+2014a, eq4
    # Mhalo = 10**np.arange(6, 17,0.001)
    alpha = 1.92
    mstar = 3e6*(10**log_mhalo/1e10)**alpha

    if do_print == True:
        print(">> Note that Garrison-Kimmel+2014a derivation is for logM*=6-8, FIRE use Kroupa IMF\n")

    return np.log10(mstar) 
