import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'

def calc_spec_snr(wave, flux, error, continuum=[],
                  dlambda=10, npix_per_res=6, figname='',
                  snr_sample_point=[1120, 1170, 1280,
                                     1350, 1490,
                                     1520, 1600, 1650, 1700]):
    """
    Calculate spectra signal to noise ratio with input of wave, flux, err of
    a spectra. Design to G130M and G160M HST/COS data. Kinda similar to the way
    HSLA calculate their SNR, to sample over the whole range while avoiding
    major strong absorption line region. And the SNR output is per resolution element.
    default to 6.

    Used to use npix = 6 /10 for G130M/G160M resolution element, but now cannot
    find the exact value for G160M anymore. From the COS cycle 27 instrument
    handbook, find that the ETC assume 6 pixels for G130M for sure. For G160M,
    gonna use npix=6 as well. Then compare the values from the value from HSLA.

    History:
    03/19/2020, now make it more generic for any COS spectra. double checked with
                HSLA SN result, values are consistent.
    04/10/2019, This is to calculate the SNR for coadd_x1d.pro(method=1) result
    and HSLA, bin1 per resolution element, npix=6 for both G130M, and G160M.
    check google note python calc_coadd_snr.py IC1613-A13
    """
    #if target == 'LBQS-0101+0009':
    #    snr_sample_point = np.array([1120, 1170, 1320, 1370, 1470, 1520, 1620])
    #if target == 'PG0044+030':
    #    ssnr_sample_point = np.array([1120, 1170, 1370, 1470, 1520, 1620])

    ###### now plot things! #####
    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(111)
    fs = 16

    ax.plot(wave, flux, color='k', lw=0.1, label='coadded spec')
    ax.plot(wave, error, color=plt.cm.Greys(0.3), lw=0.1, label='coadd error')
    if len(continuum) > 0:
        ax.plot(wave, continuum, color='r', lw=0.5, label='continuum')

    ymin=0
    ymax = np.nanmedian(flux)*5
    #ax.plot(snr_sample_point, all_snr, color=plt.cm.Reds(0.8), marker='o',
    #        linestyle='-', label='<SNR>=%.1f'%(np.nanmean(all_snr)))
    if 'bin3' not in figname:
        ### if bin3 filter is kinda corase, works for now, but might want to change later
        # YZ, 04/20/2020.
        all_snr = np.zeros(len(snr_sample_point))
        for i, left in enumerate(snr_sample_point):
            right = left+dlambda
            ind = np.all([wave>left, wave<right], axis=0)
            if len(flux[ind]) == 0 or np.nanmedian(error[ind])==0:
                all_snr[i] = np.nan
            else:
                all_snr[i] = np.nanmedian(flux[ind])/np.nanmedian(error[ind])*np.sqrt(npix_per_res)
            print("%d A:  %.1f"%(snr_sample_point[i], all_snr[i]))
        print("Mean: %.1f"%(np.nanmean(all_snr)))

        for i, left in enumerate(snr_sample_point):
            right = left+dlambda
            ax.fill_between([left, right], ymin, ymax, color=plt.cm.Blues(0.6), alpha=0.7, label=None)
            ax.text(left, ymax*0.7, 'SNR=%.1f'%(all_snr[i]), rotation=90)
            ax.set_title("Mean SNR: %.1f"%(np.nanmean(all_snr)))
    else:
        all_snr = []

    ax.set_xlim(np.nanmin(wave)-20, np.nanmax(wave)+20)
    ax.set_ylim(ymin, ymax)
    ax.legend(fontsize=fs-4)
    ax.minorticks_on()
    ax.grid('on', linestyle='--')
    ax.set_xlabel('Wavelength (A)', fontsize=fs)
    ax.set_ylabel('Flux', fontsize=fs)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fs-4)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fs-4)
    fig.tight_layout()

    if len(figname) != 0:
        fig.savefig(figname)
    else:
        import datetime
        now = datetime.datetime.now().strftime("%Y-%m-%d")
        fig.savefig('./spec_snr_%s.pdf'%(now))

    return all_snr

def get_snr(spec, wmin, wmax, has_continuum=False, segment='G130M', per_res_ele=False):
    '''
    03/19/2020, this is now abandoned. YZ.

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

if __name__ == '__main__':
    import numpy as np
    import sys
    import os
    # import astropy.io.fits as fits
    # e.g.;
    # python calc_spec_snr.py /Users/Yong/Dropbox/GitRepo/IC1613/targets/IC1613-A13/hsla/IC1613-A13_coadd_FUVM_final_lpALL.fits

    filename = sys.argv[1]
    from linetools.spectra import io as tio
    spec = tio.readspec(filename)
    wave = np.asarray(spec.wavelength)
    flux = np.asarray(spec.flux)
    err = np.asarray(spec.sig)

    # tarname = filename.split('/')[-1]
    # figname = './snr/snr_%s.pdf'%(tarname.split('_')[0])
    figname = filename+'_snr.pdf'
    print('>> Saving to ', figname)
    all_snr = calc_spec_snr(wave, flux, err, figname=figname)
