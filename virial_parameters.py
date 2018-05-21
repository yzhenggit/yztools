def virial_radius(Mh, Delta_c=200):
    # Mh is the virial mass of the dark-matter halo 
    # Delta_c is the over-density with respect to the critical density of the Universe 
    
    from astropy.cosmology import Planck15
    import astropy.constants as const
    import numpy as np
    import astropy.units as u
    
    Mh = Mh * const.M_sun
    rho_c = 3*Planck15.H0**2 / (8*np.pi*const.G) # critial density 3H^2/8piG
    
    Rvir = (4*np.pi/3)**(-1/3) * Mh**(1/3) * Delta_c**(-1/3) * rho_c**(-1/3)
    
    return Rvir.to(u.kpc)

def virial_temperature(Mh, Delta_c=200):
    # Mh is the virial mass of the dark-matter halo 
    # Delta_c is the over-density with respect to the critical density of the Universe 
    
    import astropy.constants as const
    import numpy as np
    import astropy.units as u
    
    mu = 0.56 # for a primordial component of (X, Y, Z)=(0.75, 0.25, 0)
    Rvir = virial_radius(Mh, Delta_c=Delta_c)
    Tvir = mu*const.m_p/(5*const.k_B)*const.G*Mh*const.M_sun/Rvir
    
    return Tvir.cgs
