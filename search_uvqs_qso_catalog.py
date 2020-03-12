def search_uvqs_qso_catalog(gal_name, gal_ra, gal_dec, gal_dist_kpc, within_radius_kpc=100):
    """
    Search new QSOs within certain radius of (gal_ra, gal_dec) from UVQS catalog.

    gal_ra: ra of host galaxy, in unit of degree
    gal_dec: dec of host galaxy, in unit of degree
    gal_dist_kpc: distance of host galaxy, in unit of kpc
    within_radius_kpc: search sightlines within this radius, in unit of kpc
    """

    import numpy as np
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery import mast
    from astropy.table import Table

    # coordinate setup for central galaxy
    gal_coord = SkyCoord(ra=gal_ra*u.deg, dec=gal_dec*u.deg, frame='icrs', distance=gal_dist_kpc*u.kpc)
    print("*"*90)
    print("Searching QSO (UVQS) within %.1f kpc of %s (RA=%.4f, DEC=%.4f, l=%.4f, b=%.4f)"%(within_radius_kpc,
                           gal_name, gal_ra, gal_dec, gal_coord.galactic.l.degree, gal_coord.galactic.b.degree))

    # read in the QSO catalog
    cat_dir = '/Users/Yong/Dropbox/Databucket'
    uvqs = Table.read(cat_dir+'/hlsp_uvqs_multi_multi_all_multi_v1_redshifts_monroe16.fits', format='fits')

    # search within certain radius
    found_id = []
    for j in range(len(uvqs)):
        qcoord = SkyCoord(ra=uvqs['RA'][j]*u.degree, dec=uvqs['DEC'][j]*u.degree,
                          distance=gal_dist_kpc*u.kpc, frame='icrs')
        impact = gal_coord.separation_3d(qcoord)
        if impact.kpc <= within_radius_kpc:
            found_id.append(j)

    if len(found_id) == 0:
        print("Found nothing.")
    else:
        print("%25s  %6s  %6s  %10s  %10s  %10s  %10s  %7s"%('QSO', 'FUVmag', 'z', 'ra', 'dec', 'l', 'b', 'impact'))
        c1, c2, c3, c4, c5, c6, c7, c8, c9 = [], [], [], [], [], [], [], [], []
        for i in found_id:
            qra, qdec = uvqs['RA'][i], uvqs['DEC'][i]
            qcoord = SkyCoord(ra=qra*u.degree, dec=qdec*u.degree, frame='icrs', distance=gal_dist_kpc*u.kpc)
            impact = gal_coord.separation_3d(qcoord)

            c1.append(uvqs['NAME'][i])
            c2.append(uvqs['FUV'][i])
            c3.append(uvqs['Z'][i])
            c4.append(qra)
            c5.append(qdec)
            c6.append(qcoord.galactic.l.degree)
            c7.append(qcoord.galactic.b.degree)
            c8.append(impact.kpc)
            c9.append(impact.kpc/within_radius_kpc)

        sortinds = np.argsort(c8)
        for k in sortinds:
            print('%25s  %6.2f  %6.3f  %10.4f  %10.4f  %10.4f  %10.4f  %7.1f(%.2f)'%(c1[k], c2[k], c3[k], c4[k], c5[k],
                                                                               c6[k], c7[k], c8[k], c9[k]))
    print("\n")
