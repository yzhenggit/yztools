def radec_deg2hms(ra_deg, dec_deg, do_print=False):
    """
    0316/2020,YZ.
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import numpy as np

    coord = SkyCoord(ra=ra_deg, dec=dec_deg, frame='icrs', unit=(u.deg, u.deg))
    ra_hms = coord.icrs.ra.to_string(u.hour)
    dec_dms = coord.icrs.dec.to_string(u.deg)

    if do_print == True:
        print(">>> Input: ra=%.4f deg, dec=%.4f deg"%(ra_deg, dec_deg))
        print(">>> Output: ra=%s, dec=%s"%(ra_hms, dec_dms))

    return ra_hms, dec_dms

if __name__ == '__main__':
    import sys
    import  numpy as np

    ra_deg = np.float(sys.argv[1])
    dec_deg = np.float(sys.argv[2])
    radec_deg2hms(ra_deg, dec_deg, do_print=True)
