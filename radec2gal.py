def radec2gal(ra, dec, unit_tag='deg'):

    from astropy.coordinates import SkyCoord
    import numpy as np
    import astropy.units as u

    if unit_tag.lower() == 'deg':
        units = (u.deg, u.deg)
    else: #  unit.lower() in ['hour', 'hourangle']:
        units = (u.hour, u.deg)

    gal_coord = SkyCoord(ra=ra, dec=dec, frame='icrs', unit=units)
    l_deg = gal_coord.galactic.l.deg
    b_deg = gal_coord.galactic.b.deg

    ra_deg = gal_coord.icrs.ra.deg
    dec_deg = gal_coord.icrs.dec.deg

    ra_hms = gal_coord.icrs.ra.to_string(u.hour)
    dec_dms = gal_coord.icrs.dec.to_string(u.deg)

    print(">> ra, dec = %.4f, %.4f "%(ra_deg, dec_deg))
    print(">> ra, dec = %s, %s"%(ra_hms, dec_dms))
    print(">> l, b =%.4f, %.4f"%(l_deg, b_deg))

    return l_deg, b_deg, ra_deg, dec_deg, ra_hms, dec_dms

if __name__ == '__main__':
    import sys
    import numpy as np
    ra = sys.argv[1]
    dec = sys.argv[2]
    unit_tag = sys.argv[3]
    if unit_tag.lower() == 'deg':
        ra = np.float(ra)
        dec = np.float(dec)
    res = radec2gal(ra, dec, unit_tag)
