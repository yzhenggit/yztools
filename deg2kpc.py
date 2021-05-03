def deg2kpc(impact_deg, distance_kpc, do_print=False):
    """
    Calculation the physical size of certain length (impact) at distance (dist)

    dist in kpc and impact in deg
    return impact in kpc 
    Note: only good for very close stuff, and only good for rough calculation
    """
    # impact and dist both in kpc
    impact_kpc = impact_deg*(distance_kpc/206265)*3600

    if do_print == True:
        print('%.4f deg at %.1f kpc is'%(impact_deg, distance_kpc))
        print('%.4f kpc'%(impact_kpc))
    return impact_kpc


if __name__ == "__main__":
    import sys
    import numpy as np
    impact_deg = np.float(sys.argv[1])
    distance_kpc = np.float(sys.argv[2])
    res = deg2kpc(impact_deg, distance_kpc, do_print=True)
