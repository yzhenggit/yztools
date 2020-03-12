def deg2kpc(theta_deg, distance_kpc, yes_print=True):
    """
    Calculation the angular size of certain length (impact) at distance (dist)

    dist and impact both in kpc
    return value in degree
    """
    # impact and dist both in kpc
    impact_kpc = theta_deg * (distance_kpc/206265) * 3600
    if yes_print == True:
        print('%.4f deg at %.1f kpc is'%(theta_deg, distance_kpc))
        print('%.2f kpc'%(impact_kpc))
    return impact_kpc

if __name__ == "__main__":
    import sys
    import numpy as np
    theta_deg = np.float(sys.argv[1])
    distance_kpc = np.float(sys.argv[2])
    res = deg2kpc(theta_deg, distance_kpc)
