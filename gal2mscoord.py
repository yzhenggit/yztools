
def gal2mscoord(l_deg, b_deg):
    """
    Transform from galactic coordinates l, b to magellanic coordinates

    History:
    03/12/2020, YZ, Berkeley

    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from gala.coordinates import MagellanicStreamNidever08

    gal_coord = SkyCoord(l=l_deg, b=b_deg, frame='galactic', unit=(u.deg, u.deg))

    # convert coordinates to mangellanic
    from gala.coordinates import MagellanicStreamNidever08
    ms_coords = gal_coord.transform_to(MagellanicStreamNidever08)
    ms_l = ms_coords.L.deg
    ms_b = ms_coords.B.deg

    print(">> l, b =%.2f, %.2f"%(l_deg, b_deg))
    print(">> ms_l, ms_b = %.2f, %.2f (Nidever+2008)"%(ms_l, ms_b))
    return ms_l, ms_b

if __name__ == '__main__':
    import sys
    import numpy as np
    l_deg = np.float(sys.argv[1])
    b_deg = np.float(sys.argv[2])
    ms_l, ms_b = gal2mscoord(l_deg, b_deg)
