import numpy as np
import astropy.io.fits as fits

def read_ionline(ql, qb, qname, ionline):
    
    specdir = '/Users/Yong/Dropbox/HSLA_Feb16/QSOSpec_YZ/QSOALS/'
    ionfile = '%s/%s/lines/%s_%s.fits.gz'%(specdir, qname, qname, ionline)
    
    spec = fits.open(ionfile)
    flux = spec[1].data['NORMFLUX']
    err = spec[1].data['NORMERR']
    
    from yztools.vhelio2vlsr import vhelio2vlsr_Westmeier
    vcorr = vhelio2vlsr_Westmeier(0, ql, qb, doradec=False)
    vlsr = spec[1].data['VELOCITY']+vcorr
    
    return vlsr, flux, err
    
def nv_aodm(flux, err, ionline):
    
    from yzSpec.find_line_data import find_line_data
    line_data = find_line_data(ionline)
    wave = line_data['wave']
    fval = line_data['fval']
    
    # This a modification from eqwrange.pro in YongIDL from Jess
    # decide the saturation level for the spectra
    sat_limit = 0.10  # COS tends to become saturated when flux <= 0.1
    ## sat_ind = np.any([flux<=sat_limit, flux<=error], axis=0) # not sure if f<=err would make weird profiles
    sat_ind = flux <= sat_limit
    flux[sat_ind] = err[sat_ind]

    # nv = 3.768e14*tau/(wave*fval)*delv  # GOD DAMN WRONG!! 
    tau = -np.log(flux) # the flux and error have been continuum normalized. 
    nv = 3.768e14*tau/(wave*fval)
    nverr = np.fabs(3.768e14/(wave*fval)*(err/flux))  # error proporgation 

    return nv, nverr


def cal_N_Vcent(ql, qb, qname, ionline, vmin=-100, vmax=100):

    vlsr, flux, err = read_ionline(ql, qb, qname, ionline)
    nv, nverr = nv_aodm(flux, err, ionline)
   
    # ok, this one is completely ok, no need to adjust for low or high velocity. 01/05/2018
    delv = np.fabs(np.mean(vlsr[1:] - vlsr[:-1]))

    ind = np.all([vlsr>=vmin, vlsr<=vmax], axis=0)
    finite = np.all([np.isfinite(nv), np.isfinite(nverr)], axis=0)
    vlsr = vlsr[np.all([ind, finite], axis=0)]
    nv = nv[np.all([ind, finite], axis=0)]
    nverr = nverr[np.all([ind, finite], axis=0)]

    colN = (nv*delv).sum()                   # checked
    colNerr = np.sqrt(np.sum((nverr*delv)**2))   # checked
    # lgN = np.log10(colN)
    # lgNerr = np.fabs(colNerr/colN/np.log(10))
    
    # centroid velocity                            # checked 
    vc = (vlsr*nv*delv).sum()/colN
    pA = (vlsr*nv*delv).sum()
    sigpA = np.sqrt(np.sum((vlsr*nverr*delv)**2))
    vcerr = np.fabs(vc)*np.sqrt((sigpA/pA)**2 + (colNerr/colN)**2)

    # second moment, Doppler width
    b_dp = np.sqrt(((vlsr-vc)**2*nv*delv).sum()/colN)
    b_dp_err = np.nan

    return colN, colNerr, vc, vcerr, b_dp, b_dp_err

#===============================================================================
# The following can be deleted later once finish testing the HLSP product.
def read_ionline_testhlsp(ql, qb, qname, ionline):
    import os 
    specdir = '/Users/Yong/Dropbox/HSLA_Feb16/QSOSpec_YZ/hlsp_qsoals/'
    tardir = '%s/%s/linedata_uv/'%(specdir, qname.lower())
    files = os.listdir(tardir)
    ifile = ''
    for ifile in files:
        if ionline.lower() in ifile.replace('-', ''): 
            break

    if len(ifile) == 0:
        print('Can\'t find the file. '+qname)
        import sys
        sys.exit(1)

    ionfile = '%s/%s'%(tardir, ifile)

    spec = fits.open(ionfile)
    flux = spec[1].data['NORMFLUX']
    err = spec[1].data['NORMERR']

    from yztools.vhelio2vlsr import vhelio2vlsr_Westmeier
    vcorr = vhelio2vlsr_Westmeier(0, ql, qb, doradec=False)
    vlsr = spec[1].data['VELOCITY']+vcorr

    return vlsr, flux, err


def cal_N_Vcent_testhlsp(ql, qb, qname, ionline, vmin=-100, vmax=100):

    vlsr, flux, err = read_ionline_testhlsp(ql, qb, qname, ionline)
    nv, nverr = nv_aodm(flux, err, ionline)
  
    # ok, this one is completely ok, no need to adjust for low or high velocity. 01/05/2018
    delv = np.fabs(np.mean(vlsr[1:] - vlsr[:-1]))

    ind = np.all([vlsr>=vmin, vlsr<=vmax], axis=0)
    finite = np.all([np.isfinite(nv), np.isfinite(nverr)], axis=0)
    vlsr = vlsr[np.all([ind, finite], axis=0)]
    nv = nv[np.all([ind, finite], axis=0)]
    nverr = nverr[np.all([ind, finite], axis=0)]

    colN = (nv*delv).sum()                   # checked
    colNerr = np.sqrt(np.sum((nverr*delv)**2))   # checked
    # lgN = np.log10(colN)
    # lgNerr = np.fabs(colNerr/colN/np.log(10))
   
    # centroid velocity                            # checked 
    vc = (vlsr*nv*delv).sum()/colN
    pA = (vlsr*nv*delv).sum()
    sigpA = np.sqrt(np.sum((vlsr*nverr*delv)**2))
    vcerr = np.fabs(vc)*np.sqrt((sigpA/pA)**2 + (colNerr/colN)**2)

    # second moment, Doppler width
    b_dp = np.sqrt(((vlsr-vc)**2*nv*delv).sum()/colN)
    b_dp_err = np.nan

    return colN, colNerr, vc, vcerr, b_dp, b_dp_err

