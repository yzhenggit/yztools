def mstar2mhalo(mstar):

    """Moster+2010's halo abundance matching. Based on Joo's code. """

    import astropy.constants as const
    from astropy.cosmology import Planck15
    import astropy.units as u
    import numpy as np

    M1 = 10**11.884 # Msun, from Table 4 or 6 in Moster+2010
    mM0 = 0.0282
    beta = 1.06
    gamma = 0.556
    Mhalo = 10**np.arange(7,13,0.001)

    # from moster equation 2
    Mstar = Mhalo * (2.*mM0*((Mhalo/M1)**(-beta)+(Mhalo/M1)**(gamma))**(-1.))

    # interpolate this to mstar
    from scipy import interpolate
    func = interpolate.interp1d(Mstar, Mhalo)
    mhalo = func(mstar)
    print(">> logMh = %.2f (%.1e) for logMstar = %.2f (%.1e)"%(np.log10(mhalo), mhalo,
                                                            np.log10(mstar), mstar))
    print(">> Note that Moster+2012 derivation is for logM*=8.5-11.85, see their p4\n")
    return mhalo

if __name__ == "__main__":
    import sys
    import numpy as np
    mstar = np.float(sys.argv[1])
    mstar2mhalo(mstar)
