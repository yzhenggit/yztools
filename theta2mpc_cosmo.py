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
        print('Input: ')
        print('    z = %.4f'%(redshift))
        print('    theta=%.2f arcmin, %.1f arcsec'%(theta_arcmin, theta_arcmin*60))
        print('Output: ')
        print('    dA = %.2f Mpc'%(dA.value))
        print('    dL = %.2f Mpc'%(dL.value))
        print('    %.2f arcmin = %.6f mpc / %.2f kpc'%(theta_arcmin, impact_mpc.value, impact_mpc.value*1000))
    return impact_mpc, dL, dA

if __name__ == "__main__":
    import sys
    import numpy as np
    theta = np.float(sys.argv[1])
    theta_unit = sys.argv[2]
    if theta_unit.lower() in ['deg', 'degree']:
        theta_arcmin = theta*60
    elif theta_unit.lower() in ['arcmin', 'arcminute', 'arc_minute', 'arc_min']:
        theta_arcmin = theta
    elif theta_unit.lower() in ['arcsec', 'arcsecond', 'arc_second', 'arc_sec']:
        theta_arcmin = theta/60.
    else:
        print('Do not recognize theta_unit %s.'%(theta_unit))
        sys.exit()

    dist_note = sys.argv[3]
    dist = np.float(sys.argv[4])
    if dist_note.lower() in ['z' ,'redshift']:
        redshift = dist
        res = arcmin2mpc_input_z(theta_arcmin, redshift, do_print=True)
