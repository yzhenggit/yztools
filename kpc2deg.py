def kpc2deg(impact_kpc, distance_kpc, yes_print=True):
    """
    Calculation the angular size of certain length (impact) at distance (dist)

    dist and impact both in kpc
    return value in degree
    """
    # impact and dist both in kpc
    theta_deg = impact_kpc/(distance_kpc/206265)/3600

    if yes_print == True:
        print('%.4f kpc at %.1f kpc is'%(impact_kpc, distance_kpc))
        print('%.4f deg'%(theta_deg))
        print('%.3f arcmin'%(theta_deg*60))
        print('%.1f arcsec'%(theta_deg*3600))
    return theta_deg


if __name__ == "__main__":
    import sys
    import numpy as np
    impact_kpc = np.float(sys.argv[1])
    distance_kpc = np.float(sys.argv[2])
    res = kpc2deg(impact_kpc, distance_kpc)
