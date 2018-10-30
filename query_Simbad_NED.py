import numpy as np

def query_Simbad_NED(objname, z=np.nan):
    '''
    Make use to astroquery to extract basic coordiantes/redshift information 
    
    z: redshift. Normally NED query would give an answer. Provide z if NED answer is not satifying
    return: astropy SkyCoord object 
    '''
    
    import numpy as np
    from astroquery.simbad import Simbad
    from astroquery.ned import Ned 
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astropy.cosmology import Planck15 as cosmo

    sim = Simbad.query_object(objname)
    try:
        test1 = sim['RA']
    except: 
        print('Danger: %s cannot be recognized by Simbad, try other names.'%(objname))
        return False
    
    try:
        ned = Ned.query_object(objname)
    except:
        print('Danger: %s cannot be recognized by NED, try other names.'%(objname))
        return False

    galcoord = SkyCoord(ra=sim['RA'], dec=sim['DEC'], unit=(u.hour, u.deg), frame='icrs')
    if np.isfinite(z):
        gz = z
    else: 
        gz = ned['Redshift'][0]
    
    scale_kpc_arcmin = cosmo.kpc_comoving_per_arcmin(gz) 
    scale_pc_arcsec = scale_kpc_arcmin.to(u.pc/u.arcsec)
    
    print('Simbad: RA  = %.4f deg '%(galcoord.icrs.ra.deg), galcoord.icrs.ra)
    print('Simbad: DEC = %.4f deg '%(galcoord.icrs.dec.deg), galcoord.icrs.dec)
    print('Simbad: l = %.4f deg, b = %.4f deg'%(galcoord.galactic.l.deg, galcoord.galactic.b.deg))
    print('NED: z = %.4f; Input z = %.4f'%(ned['Redshift'], z))
    print('Scale at z = %.4f: %.1f kpc/arcmin, %.1f pc/arcsec'%(gz, abs(scale_kpc_arcmin.value), 
                                                                abs(scale_pc_arcsec.value)))
    print('Returning astropy SkyCoord object...')
    return galcoord
