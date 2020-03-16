def radec_hms2deg(ra_hms, dec_dms):
    """
    03/16/2020. YZ. 
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import numpy as np

    coord = SkyCoord(ra=ra_hms, dec=dec_dms, frame='icrs', unit=(u.hour, u.deg))
    ra_deg = coord.icrs.ra.deg
    dec_deg = coord.icrs.dec.deg

    print(">>> Input: ra=%s, dec=%s"%(ra_hms, dec_dms))
    print(">>> Output: ra=%.4f deg, dec=%.4f deg"%(ra_deg, dec_deg))

    return ra_deg, dec_deg

if __name__ == '__main__':
    import sys

    ra_hms = sys.argv[1]
    dec_dms = sys.argv[2]
    radec_hms2deg(ra_hms, dec_dms)
