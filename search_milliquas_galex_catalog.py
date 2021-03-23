def search_milliquas_galex_catalog(gal_name, gal_coord1, gal_coord2,
                               gal_dist_kpc,
                               within_radius_kpc=100,
                               frame='icrs', fuvmag_limit=18.5):
    """
    Search new QSOs within certain radius of (gal_coord1, gal_coord2)
    from the cross-matched Million QSO catalog (Flesch 2017, v5.2) and Galex.

    when frame=icrs, gal_coord1 = ra, gal_coord2 = dec, in unit of deg
    when frame = galactic, gal_coord1 = gl, gal_coord_coord2 = gb, in unit of deg

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
    from yztools.kpc2deg import kpc2deg

    # coordinate setup for central galaxy
    if frame.lower() == 'icrs':
        gal_coord = SkyCoord(ra=gal_coord1*u.deg, dec=gal_coord2*u.deg, frame='icrs', distance=gal_dist_kpc*u.kpc)
    elif frame.lower() in ['galactic', 'gal', 'g']:
        gal_coord = SkyCoord(l=gal_coord1*u.deg, b=gal_coord2*u.deg, frame='galactic', distance=gal_dist_kpc*u.kpc)
    else:
        print("Do no recognize frame:%d. "%(frame))
        sys.exit(0)

    gal_ra = gal_coord.icrs.ra.degree
    gal_dec = gal_coord.icrs.dec.degree
    print("*"*90)
    print("Cross match the Million QSOs (Fleshch17; v5.2) + GALEX catalog")
    print("within %.1f kpc of %s (RA=%.4f, DEC=%.4f, l=%.4f, b=%.4f)"%(within_radius_kpc,
                           gal_name, gal_coord.icrs.ra.degree, gal_coord.icrs.dec.degree,
                           gal_coord.galactic.l.degree, gal_coord.galactic.b.degree))

    # read in the QSO catalog
    cat_dir = '/Users/Yong/Dropbox/Databucket'
    millgalex = Table.read(cat_dir+'/milliquas_v5.2_galex_jb032219.fits', format='fits')

    # this catalog is toooooo large, let's do a rough filtering first
    sepdeg = kpc2deg(within_radius_kpc, gal_dist_kpc)
    rough_dist = np.sqrt((millgalex['ra']-gal_ra)**2 + (millgalex['dec']-gal_dec)**2)
    close_ids = rough_dist <= 2*sepdeg
    close_millgalex = millgalex[close_ids]
    print("%d/%d QSOs are roughly within 2*%d kpc of the target."%(len(close_millgalex),
                                                                   len(millgalex),
                                                                   within_radius_kpc))

    # only look at targets with FUV mag within certian limits
    fuvbright_ids = close_millgalex['fuv'] <= fuvmag_limit
    final_millgalex = close_millgalex[fuvbright_ids]
    print("Still too many, limit to FUVmag<=%.2f. Found %d candidates."%(fuvmag_limit,
                                                                         len(final_millgalex)))

    # search within certain radius
    found_id = []
    for j in range(len(final_millgalex)):
        qcoord = SkyCoord(ra=final_millgalex['ra'][j]*u.degree, dec=final_millgalex['dec'][j]*u.degree,
                          distance=gal_dist_kpc*u.kpc, frame='icrs')
        impact = gal_coord.separation_3d(qcoord)
        if impact.kpc <= within_radius_kpc:
            found_id.append(j)

    if len(found_id) == 0:
        print("Found nothing.")
    else:
        print('%25s  %8s  %5s  %10s  %10s  %10s  %10s  %7s'%('QSO',
              'z', 'fuv', 'ra', 'dec', 'l', 'b', 'impact'))

        # just to sort things using impact parameters
        c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11 = [], [], [], [], [], [], [], [], [], [], []
        for i in found_id:
            qra, qdec = final_millgalex['ra'][i], final_millgalex['dec'][i]
            qcoord = SkyCoord(ra=qra*u.degree, dec=qdec*u.degree,
                              frame='icrs', distance=gal_dist_kpc*u.kpc)
            impact = gal_coord.separation_3d(qcoord)
            c1.append(final_millgalex['name'][i])
            c2.append(qra)
            c3.append(qdec)
            c4.append(qcoord.galactic.l.degree)
            c5.append(qcoord.galactic.b.degree)
            c6.append(final_millgalex['z'][i])
            c7.append(final_millgalex['Vmag'][i])
            c8.append(final_millgalex['fuv'][i])
            c9.append(final_millgalex['nuv'][i])
            c10.append(impact.kpc)
            c11.append(impact.kpc/within_radius_kpc)

        sortinds = np.argsort(c10)
        for k in sortinds:
            print('%25s  %8.3f  %5.2f  %10.4f  %10.4f  %10.4f  %10.4f  %7.1f(%.2f)'%(c1[k],
                  c6[k], c8[k], c2[k], c3[k], c4[k], c5[k], c10[k], c11[k]))
    print("\n")

if __name__ == '__main__':
    import sys
    import numpy as np

    gal_name = sys.argv[1]
    frame = sys.argv[2]  # 'icrs' or 'galactic'
    gal_coord1 = np.float(sys.argv[3]) # deg, ra, or gl
    gal_coord2 = np.float(sys.argv[4]) # deg, ra, or gb
    gal_dist_kpc = np.float(sys.argv[5])  # kpc
    within_radius_kpc = np.float(sys.argv[6]) # kpc

    search_milliquas_galex_catalog(gal_name, gal_coord1, gal_coord2,
                                   gal_dist_kpc,
                                   within_radius_kpc=within_radius_kpc,
                                   frame=frame)
