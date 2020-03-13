def gal2radec(l_deg, b_deg):
    from astropy.coordinates import SkyCoord
    import numpy as np
    import astropy.units as u

    gal_coord = SkyCoord(l=l_deg, b=b_deg, frame='galactic', unit=(u.deg, u.deg))
    ra_deg = gal_coord.icrs.ra.deg
    dec_deg = gal_coord.icrs.dec.deg

    ra_hms = gal_coord.icrs.ra.to_string(u.hour)
    dec_dms = gal_coord.icrs.dec.to_string(u.deg)

    print(">> l, b =%.4f, %.4f"%(l_deg, b_deg))
    print(">> ra, dec = %.4f, %.4f "%(ra_deg, dec_deg))
    print(">> ra, dec = %s, %s"%(ra_hms, dec_dms))

    return l_deg, b_deg, ra_deg, dec_deg, ra_hms, dec_dms

if __name__ == '__main__':
    import sys
    import numpy as np
    l_deg = np.float(sys.argv[1])
    b_deg = np.float(sys.argv[2])

    res = gal2radec(l_deg, b_deg)
