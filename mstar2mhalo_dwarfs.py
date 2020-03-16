def mstar2mhalo_dwarfs(mstar):

    """Garrison-Kimmel+2014a halo abundance matching for dwarfs logM*=6-8 """

    import astropy.constants as const
    from astropy.cosmology import Planck15
    import astropy.units as u
    import numpy as np

    ## from Garrison-Kimmel+2014a, eq4
    Mhalo = 10**np.arange(6, 12,0.001)
    alpha = 1.92
    Mstar = 3e6*(Mhalo/1e10)**alpha

    # interpolate this to mstar
    from scipy import interpolate
    func = interpolate.interp1d(Mstar, Mhalo)
    mhalo = func(mstar)
    print(">> logMh = %.2f (%.1e) for logMstar = %.2f (%.1e)"%(np.log10(mhalo), mhalo,
                                                            np.log10(mstar), mstar))
    print(">> Note that Garrison-Kimmel+2014a derivation is for logM*=6-8, see their eq4\n")
    return mhalo

if __name__ == "__main__":
    import sys
    import numpy as np
    mstar = np.float(sys.argv[1])
    mstar2mhalo_dwarfs(mstar)
