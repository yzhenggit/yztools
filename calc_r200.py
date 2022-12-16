import astropy.constants as const
from astropy.cosmology import Planck15
import astropy.units as u
import numpy as np
import sys

def calc_percentile(para_range, ndecimals=2):
    import numpy as np 
    percentile = np.percentile(para_range, [16, 50, 84])
    best_val = np.around(percentile[1], decimals=ndecimals)
    sig_lo = -np.around(np.diff(percentile)[0], decimals=ndecimals)
    sig_hi = np.around(np.diff(percentile)[1], decimals=ndecimals)
    # print(best_val, sig_lo, sig_hi)
    
    return best_val, sig_lo, sig_hi

def r200c(mhalo, delta=200): 
    # calculate r200 with respect to critical density
    # delta = 200 ## rho/rho_crit=200
    rho_c = Planck15.critical_density0
    r200c = ((3*mhalo/(4.*np.pi*delta*rho_c))**(1./3.)).to(u.kpc)
    return r200c.value

def r200m(mhalo, delta=200): 
    # calculate r200 with respect to matter density
    # delta = 200 ## rho/rho_crit=200
    rho_m = Planck15.critical_density0 * Planck15.Om0
    r200m = ((3*mhalo/(4.*np.pi*delta*rho_m))**(1./3.)).to(u.kpc)
    return r200m.value 

def calc_r200(logmhalo, sig_logmhalo=0., label='200m', do_print=False, Ntrials=4000):

    """
    Assuming an isothermal sphere. Based on Joo's code.
    Note we required log values as input!

    label: need to calrify whehther logmhalo is based on critical density (c) or mean matter density (m)

    Ntrials: default to 4000 mcmc trials to estimate uncertainties 
    """

    # print("calc_r200: Note! requried logmhalo instead mhalo input!")
    #if mstar != 0:
    #    if use_who == 'M10':
    #        from yztools.mstar2mhalo import mstar2mhalo
    #        mhalo = mstar2mhalo(mstar, do_print=do_print)
    #    else: # 'GK14'
    #        from yztools.mstar2mhalo_dwarfs import mstar2mhalo_dwarfs
    #        mhalo = mstar2mhalo_dwarfs(mstar, do_print=do_print)
    # print(mhalo)

    if label not in ['200m', '200c']: 
        print("label must be either 200m or 200c")
        sys.exit(0)

    if sig_logmhalo == 0.: 
        mhalo = 10**logmhalo * u.Msun 

        # calculate r200 with respect to critical density
        if label == '200c': 
            best_r200_kpc = np.around(r200c(mhalo), decimals=1)
        else# label == '200m': 
        # calculate r200 with respect to matter density
            best_r200_kpc = np.around(r200m(mhalo), decimals=1)
        sig_low_r200, sig_hi_r200 = 0., 0.

    else: # with error bars, use MCMC to estimate uncertainties
        # e.g., https://astrofrog.github.io/py4sci/_static/Practice%20Problem%20-%20Monte-Carlo%20Error%20Propagation%20-%20Sample%20Solution.html
        logmhalo_range = np.random.normal(logmhalo, sig_logmhalo, Ntrials)
        r200_range = np.zeros(Ntrials)+np.nan
        for i in range(Ntrials): 
            mhalo = 10**(logmhalo_range[i]) * u.Msun
            if label == '200c': 
                r200_range[i] = r200c(mhalo)
            else: 
                r200_range[i] = r200m(mhalo)
        
        best_r200_kpc, sig_low_r200, sig_hi_r200 = calc_percentile(r200_range, ndecimals=1)

    if do_print == True:
        if label == '200c'
            print(">> Use critical density, delta_c=rho/rho_c=200")
            print(">> r200c = {}+{}{} kpc\n".format(best_r200c_kpc, sig_hi_r200c, sig_low_r200c))
        else: 
            print(">> Use mean matter density, delta_c=rho/(rho_c*Omega_m)=200")
            print(">> r200m = {}+{}{} kpc\n".format(best_r200m_kpc, sig_hi_r200m, sig_low_r200m))

    return best_r200_kpc, sig_low_r200, sig_hi_r200

if __name__ == "__main__":
    import sys
    import numpy as np
    logmhalo = float(sys.argv[1])
    sig_logmhalo = float(sys.argv[2])
    label = sys.argv[3].lower() # either 200c or 200m 
 
    calc_r200(logmhalo, sig_logmhalo=sig_logmhalo, label=label, do_print=True)
