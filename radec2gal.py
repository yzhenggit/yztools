def radec2gal(ra_deg, dec_deg):
    """
    input:
    $ python radec2gal.py 30 50
    or:
    from yztools.radec2gal import radec2gal
    radec2gal(ra_deg, dec_deg)

    will return:
    l_deg, b_deg, ra_deg, dec_deg, ra_hms, dec_dms, ms_l, ms_b
    """
    from astropy.coordinates import SkyCoord
    import numpy as np
    import astropy.units as u

    gal_coord = SkyCoord(ra=ra_deg, dec=dec_deg,
                         frame='icrs', unit=(u.deg, u.deg))
    l_deg = gal_coord.galactic.l.deg
    b_deg = gal_coord.galactic.b.deg

    ra_deg = gal_coord.icrs.ra.deg
    dec_deg = gal_coord.icrs.dec.deg

    ra_hms = gal_coord.icrs.ra.to_string(u.hour)
    dec_dms = gal_coord.icrs.dec.to_string(u.deg)

    print(">> ra, dec = %.4f  %.4f "%(ra_deg, dec_deg))
    print(">> ra, dec = %s  %s"%(ra_hms, dec_dms))
    print(">> l, b =%.4f  %.4f"%(l_deg, b_deg))

    ## bonus, also do ms transformation ##
    from yztools.gal2mscoord import gal2mscoord
    ms_l, ms_b = gal2mscoord(l_deg, b_deg)

    res = {'l_deg': l_deg, 'b_deg': b_deg,
           'ra_deg': ra_deg, 'dec_deg': dec_deg,
           'ra_hms': ra_hms, 'dec_dms': dec_dms,
           'ms_l': ms_l, 'ms_b': ms_b}

    # return l_deg, b_deg, ra_deg, dec_deg, ra_hms, dec_dms, ms_l, ms_b
    return res

if __name__ == '__main__':
    import sys
    import numpy as np
    ra_deg = np.float(sys.argv[1]) # deg
    dec_deg = np.float(sys.argv[2]) # deg

    radec2gal(ra_deg, dec_deg)
