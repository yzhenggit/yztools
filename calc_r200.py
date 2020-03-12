def calc_r200(mstar = 0., mhalo = 1e12):

    """
    Assuming an isothermal sphere. Based on Joo's code.
    """

    import astropy.constants as const
    from astropy.cosmology import Planck15
    import astropy.units as u
    import numpy as np
    import sys

    if mstar != 0:
        from yztools.mstar2mhalo import mstar2mhalo
        mhalo = mstar2mhalo(mstar)

    # print(mhalo)

    mhalo = mhalo*u.Msun
    delta = 200 ## rho/rho_crit=200

    # calculate r200 with respect to critical density
    rho_c = Planck15.critical_density0
    r200_c = ((3*mhalo/(4.*np.pi*delta*rho_c))**(1./3.)).to(u.kpc)
    print(">> If use critical density, delta_c=rho/rho_c=200")
    print(">> r200c = %.2f kpc\n"%(r200_c.value))

    # calculate r200 with respect to matter density
    rho_m = Planck15.critical_density0 * Planck15.Om0
    r200_m = ((3*mhalo/(4.*np.pi*delta*rho_m))**(1./3.)).to(u.kpc)
    print(">> If use critical MATTER density, delta_c=rho/(rho_c*Omega_m)=200")
    print(">> r200m = %.2f kpc\n"%(r200_m.value))

    return r200_c, r200_m, mhalo.value

if __name__ == "__main__":
    import sys
    import numpy as np
    mstar = np.float(sys.argv[1])
    calc_r200(mstar = mstar)
