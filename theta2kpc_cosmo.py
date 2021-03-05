#def deg2mpc_input_dL(theta_deg, dA_mpc, do_print=False):

#def deg2mpc_input_dA(theta_deg, dA_mpc, do_print=False):
# example:
# for redshfit
# $ python theta2kpc_cosmo.py 7.12 arcmin 0.1 z 
# for distance
# $ python theta2kpc_cosmo.py 7.12 arcmin 48.2 kpc
def arcmin2kpc_input_z(theta_arcmin, redshift, do_print=False):
    """
    Calculation the linear size (mpc) of a subject wtih angular size of
    theta at redshift z

    Output: angular separation in Mpc,
            luminosity distance d_L in Mpc,
            angular diameter distance d_A in Mpc

    Example:
    $ python theta2kpc_cosmo.py 7.12 arcmin 48.2 kpc
    """
    # impact and dist both in kpc
    from astropy.cosmology import Planck15 as cosmo
    import astropy.units as u

    dL = cosmo.luminosity_distance(redshift) # mpc
    dA = cosmo.angular_diameter_distance(redshift) # mpc
    linear_size_1arcsec_at_z = dA/206264.8
    impact_mpc = theta_arcmin*60*linear_size_1arcsec_at_z
    impact_kpc = impact_mpc*1000
    if do_print == True:
        print('*'*60)
        print('Input: ')
        print('    z = %.4f'%(redshift))
        print('    theta=%.2f arcmin, %.1f arcsec'%(theta_arcmin, theta_arcmin*60))
        print('Output: ')
        print('    dA = %.2f Mpc'%(dA.value))
        print('    dL = %.2f Mpc'%(dL.value))
        print('    %.2f arcmin = %.6f mpc / %.2f kpc'%(theta_arcmin, impact_mpc.value, impact_mpc.value*1000))
        print('(Planck15: H0=67.8 km/s/Mpc, Omeba_b = 0.484, Omega_Lambda=0.692, Omega_m=0.308)')

    result = {'input_arcsec':theta_arcmin*60,
              'input_arcmin':theta_arcmin,
              'impact_kpc': impact_kpc.value,
              'input_z': redshift,
              'dL_mpc': dL,
              'dA_mpc': dA}
    return result

def arcsec2kpc_local(theta_arcsec, dist_kpc, do_print=True):
    """
    Input
    theta_arcsec: angular size of object in unit of arcsec
    dist_kpc: distance to the target in unit of kpc

    output:
    size_kpc

    Example:
    $ python theta2kpc_local.py 7.12 arcmin 48.2 kpc
    """

    size_kpc = dist_kpc/206265 * theta_arcsec
    size_pc = size_kpc*1000

    if do_print == True:
        print('*'*60)
        print('Input: ')
        print('    theta=%.2f arcmin, %.1f arcsec'%(theta_arcsec/60, theta_arcsec))
        print('    distance = %.1f kpc'%(dist_kpc))
        print('Output: ')
        print('    %.2f arcsec = %.5f kpc / %.2f pc'%(theta_arcsec, size_kpc, size_pc))

    result = {'theta_arcsec':theta_arcsec,
              'size_kpc': size_kpc,
              'input_d_kpc': dist_kpc}
    return result

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

    # then check the distance input
    dist = np.float(sys.argv[3])
    dist_unit = sys.argv[4]
    if dist_unit.lower() in ['kpc']:
        dist_kpc = dist
        theta_arcsec = theta_arcmin*60
        res = arcsec2kpc_local(theta_arcsec, dist_kpc, do_print=True)
    elif dist_unit.lower() in ['mpc']:
        dist_kpc = dist*1000
        theta_arcsec = theta_arcmin*60
        res = arcsec2kpc_local(theta_arcsec, dist_kpc, do_print=True)
    elif dist_unit.lower() in['z', 'redshift']:
        redshift = dist
        res = arcmin2kpc_input_z(theta_arcmin, redshift, do_print=True)
    else:
        print("Only know units of kpc/Kpc/MPC, mpc/Mpc/MPC, or z/redshift.")
        sys.exit()
