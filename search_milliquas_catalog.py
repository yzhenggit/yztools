def search_milliquas_catalog(gal_name, gal_ra, gal_dec,
                             gal_dist_kpc, within_radius_kpc=100):
    """
    Search new QSOs within certain radius of (gal_ra, gal_dec)
    from the Million QSO catalog by Flesch 2017, version 5.2.

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
    from yzObs.kpc2deg import kpc2deg

    # coordinate setup for central galaxy
    gal_coord = SkyCoord(ra=gal_ra*u.deg, dec=gal_dec*u.deg,
                         frame='icrs', distance=gal_dist_kpc*u.kpc)

    print("*"*90)
    print("Searching the Million QSOs (Fleshch17; v5.2)")
    print("within %.1f kpc of %s (RA=%.4f, DEC=%.4f, l=%.4f, b=%.4f)"%(within_radius_kpc,
           gal_name, gal_ra, gal_dec, gal_coord.galactic.l.degree, gal_coord.galactic.b.degree))

    # read in the QSO catalog
    cat_dir = '/Users/Yong/Dropbox/Databucket'
    milliquas = Table.read(cat_dir+'/milliquas_v5.2_flesch2017_catalog_yzprocessed.fits', format='fits')

    # this catalog is toooooo large, let's do a rough filtering first
    sepdeg = kpc2deg(within_radius_kpc, gal_dist_kpc)
    rough_dist = np.sqrt((milliquas['RA']-gal_ra)**2 + (milliquas['DEC']-gal_dec)**2)
    close_ids = rough_dist <= 2*sepdeg
    close_milliquas = milliquas[close_ids]
    print("%d/%d QSOs are roughly within 2*%d kpc of the target."%(len(close_milliquas),
                                                                   len(milliquas),
                                                                   within_radius_kpc))

    # only look at targets with 100% possibility
    q100_ids = close_milliquas['is_qso?'] == 100
    final_milliquas = close_milliquas[q100_ids]
    print("Still too many, limit to those with is_qso? == 100. Found %d candidates."%(len(final_milliquas)))

    # search within certain radius
    found_id = []
    for j in range(len(final_milliquas)):
        qcoord = SkyCoord(ra=final_milliquas['RA'][j]*u.degree, dec=final_milliquas['DEC'][j]*u.degree,
                          distance=gal_dist_kpc*u.kpc, frame='icrs')
        impact = gal_coord.separation_3d(qcoord)
        if impact.kpc <= within_radius_kpc:
            found_id.append(j)

    if len(found_id) == 0:
        print("Found nothing.")
    else:
        print("%25s  %7s  %8s  %10s  %10s  %10s  %10s  %7s"%('QSO',  'is_qso?', 'z', 'ra', 'dec', 'l', 'b', 'impact'))

        # just to sort things using impact parameters
        c1, c2, c3, c4, c5, c6, c7, c8, c9 = [], [], [], [], [], [], [], [], []
        for i in found_id:
            qra, qdec = final_milliquas['RA'][i], final_milliquas['DEC'][i]
            qcoord = SkyCoord(ra=qra*u.degree, dec=qdec*u.degree,
                              frame='icrs', distance=gal_dist_kpc*u.kpc)
            impact = gal_coord.separation_3d(qcoord)
            c1.append(final_milliquas['NAME'][i])
            c2.append(final_milliquas['is_qso?'][i])
            c3.append(final_milliquas['z'][i])
            c4.append(qra)
            c5.append(qdec)
            c6.append(qcoord.galactic.l.degree)
            c7.append(qcoord.galactic.b.degree)
            c8.append(impact.kpc)
            c9.append(impact.kpc/within_radius_kpc)

        sortinds = np.argsort(c8)
        for k in sortinds:
            print('%25s  %7d  %8.3f  %10.4f  %10.4f  %10.4f  %10.4f  %7.1f(%.2f)'%(c1[k],
                                c2[k], c3[k], c4[k], c5[k], c6[k], c7[k], c8[k], c9[k]))
    print("\n")
