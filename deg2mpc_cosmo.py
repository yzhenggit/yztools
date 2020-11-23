#def deg2mpc_input_dL(theta_deg, dA_mpc, do_print=False):

#def deg2mpc_input_dA(theta_deg, dA_mpc, do_print=False):

def arcmin2mpc_input_z(theta_arcmin, redshift, do_print=False):
    """
    Calculation the linear size (mpc) of a subject wtih angular size of
    theta_arcminat redshift z
    """
    # impact and dist both in kpc
    from astropy.cosmology import Planck15 as cosmo
    import astropy.units as u

    dL = cosmo.luminosity_distance(redshift) # mpc
    dA = cosmo.angular_diameter_distance(redshift) # mpc
    linear_size_1arcsec_at_z = dA/206264.8
    impact_mpc = theta_arcmin*60*linear_size_1arcsec_at_z
    if do_print == True:
        print('%.4f arcmin (%.1f arcsec) at z=%.4f is'%(theta_arcin, theta_arcmin*60, redshift))
        print('%.6f mpc / %.2f kpc'%(impact_mpc.value, impact_mpc.value*1000))
    return impact_mpc, dL, dA

if __name__ == "__main__":
    import sys
    import numpy as np
    theta_arcmin = np.float(sys.argv[1])
    dist = np.float(sys.argv[2])
    note = sys.argv[3]
    if note.lower() in ['z' ,'redshift']:
        redshift = dist
        res = arcmin2mpc_input_z(theta_arcmin, redshift, do_print=True)
