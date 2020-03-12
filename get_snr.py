import numpy as np

def get_snr(spec, wmin, wmax, has_continuum=False, segment='G130M', per_res_ele=False):
    '''
    To calculate the SNR within [wmin, wmax] range per pix

    wmin: the minimum wavelength
    wmax: the maximum wavelength 
    '''
    flux = spec['FLUX']
    wave = spec['WAVE']
    error = spec['ERROR']
        
    ind_wave = np.all([wave>wmin, wave<wmax], axis=0)
    s_flux = np.median(flux[ind_wave])
    n_err = np.median(error[ind_wave])
    snr_pix_flux = s_flux / n_err
    
    if has_continuum == True:
        cont = spec['CONTINUUM']
        s_cont = np.median(cont[ind_wave])
        snr_pix_cont = s_cont / n_err
    else: 
        snr_pix_cont = np.nan

    # SNR per resolution element
    if per_res_ele == True:
        npix_dict = {'G130M': 6, 'G160M':10}
        npix = npix_dict[segment]
        snr_res_flux = snr_pix_flux*np.sqrt(npix)
        snr_res_cont = snr_pix_cont*np.sqrt(npix)
        return snr_res_flux, snr_res_cont
    else: 
        return snr_pix_flux, snr_pix_cont

def get_siglevel(vel, flux, seg='G130M', per_res_ele=False):
    '''
    To calculate the SNR for the sliced lines of HSLA
    YZ noted on 07/23/18: this way of calculating the SNR is wrong, see email on Jan 29.  
    per_res_ele: set to True if want to compute SNR per resolution element
    '''
    
    isfinite = np.isfinite(flux)
    vel, flux = vel[isfinite], flux[isfinite]
    cen_flux = flux[np.all([vel>=-90, vel<=90], axis=0)]
        
    if cen_flux.size<=6: 
        snr = np.nan
    else:
        # the signal level
        indm = np.argmin(cen_flux)
        if indm<3: signal = 1-np.mean(cen_flux[0:6])
        else: signal = 1-np.mean(cen_flux[indm-3:indm+3])

        # the noise level, use -400-200 and +200+400
        l1, l2, r1, r2 = -400, -200, 200, 400
        # l1, l2, r1, r2 = -150, -90, 90, 150  # only for FeII 1143/1144
    
        indA = np.all([vel>=l1, vel<=l2], axis=0)
        if flux[indA].size<3: sigmaA = np.nan
        else: sigmaA = np.fabs(np.std(flux[indA]))
        
        indB = np.all([vel>=r1, vel<=r2], axis=0)
        if flux[indB].size<3: sigmaB = np.nan
        else: sigmaB = np.fabs(np.nanstd(flux[indB]))
        
        indAB = np.any([indA, indB], axis=0)
        if flux[indAB].size<3: sigmaAB = np.nan
        else: sigmaAB = np.fabs(np.nanstd(flux[indAB]))
        
        sigma = np.nanmin([sigmaA, sigmaB, sigmaAB])        
        if sigma == 0.: snr = np.nan
        else:
            # SNR per resolution element
            if per_res_ele == True:
                npix_dict = {'G130M': 6, 'G160M':10}
                npix = npix_dict[seg]
            else: npix = 1.
            snr = signal/sigma*np.sqrt(npix)
    return snr
