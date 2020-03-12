import datetime
import numpy as np
import astropy.io.fits as fits
from astropy.constants import c as speed_of_light_ms # speed of light, in m/s
from yzSpec.read_linelibrary import read_linelibrary

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_normspec(lt_spec, has_continuum=False):
    '''
    Read the linetool spec, if has continuum, then provide normalized profile, 
    otherwise return nan arrays
    '''

    if has_continuum == True:
        # to avoid 0 in the denominator
        continuum = lt_spec.co.value
        indnan = np.where(continuum == 0.)
        continuum[indnan] = 1.

        # normlized flux
        normflux = lt_spec.flux.value / continuum
        normflux[indnan] = np.nan

        # normlized error
        normsig = lt_spec.sig.value / continuum
        normsig[indnan] = np.nan  # still keep the nan since we want the velocity array to be continous

    else:
        arrsize = lt_spec.flux.value.size
        continuum = np.zeros(arrsize)+np.nan
        normflux = np.zeros(arrsize)+np.nan
        normsig = np.zeros(arrsize)+np.nan

    return continuum, normflux, normsig


def coadd_exposure_info(target, category):
    '''
    Extract the exposure information from the x1d header of each spec went into the coadd. 
    
    '''

    import os
    import numpy as np
    from astropy.table import Table
    
    datadir = os.path.expanduser('~')+'/Dropbox/HSLA_Feb16/datapile/'+category+'/'+target
    
    exptb = Table.read(datadir + '/all_exposures.txt', format='ascii')
    exptime1 = np.sum(exptb['Exptime'][exptb['Grating']=='G130M'])  # in s
    exptime2 = np.sum(exptb['Exptime'][exptb['Grating']=='G160M'])  # in s
    
    # get the exposure start/end time of each observation
    expstart = 100000
    expend = -99
    for rootname in exptb['Rootname']:
        spec = '%s/%s_x1d.fits.gz'%(datadir, rootname)
        header = fits.open(spec)[1].header
        ## find the ealiest start time 
        if expstart > header['EXPSTART']:  
            date_obs = header['DATE-OBS'] # UT date of start of observation (yyyy-mm-dd)  
            time_obs = header['TIME-OBS'] # UT time of start of observation (hh:mm:ss)  
            expstart = header['EXPSTART'] # exposure start time (Modified Julian Date)
        ## find the latest time 
        if expend < header['EXPEND']:
            expend = header['EXPEND'] # exposure end time (Modified Julian Date)
        
    return date_obs, time_obs, expstart, expend, exptime1, exptime2

def create_primary_header_UV(target_info, category):
    import astropy.io.fits as fits
    todate = str(datetime.datetime.now())
    hdulist = fits.open(target_info['DATAFILE'])

    ## method 1: get the exposure info by direclty reading the co-added spectra
    ##           but this one takes almost forever if we run it for every single target
    # expo_info = coadd_exposure_info(target_info['NAME'], category)

    ## method 2; pre-processed exposure for all targets, and save the info somewhere
    ##           now we just need to read that table in 
    import os
    expo_str = np.load(os.path.expanduser('~')+'/Dropbox/HSLA_Feb16/code/tables/exposure_structures.npy').item()
    try: 
        expo_info = expo_str[target_info['NAME']] 
    except KeyError:
        logger.info('Do not have exposure info for %s.'%(target_info['NAME']))
        return False

    # Add the following Keywords to the Primary Header
    # according to S. Flemming's instruction
    prihead = hdulist[0].header
    prihead['TELESCOP'] = 'HST'
    prihead['INSTRUME'] = 'COS'
    prihead['TARGNAME'] = (target_info['NAME'], 'Target name; HSLA (Peeples+2017)')
    prihead['RA_TARG']  = (round(target_info['RA'], 4), 'Right Ascension (deg); HSLA (Peeples+2017)')
    prihead['DEC_TARG'] = (round(target_info['DEC'], 4), 'Declination (deg); HSLA (Peeples+2017)') 
    prihead['DATE-OBS'] = (expo_info[0], 'UT date of start of the 1st obs (yyyy-mm-dd)')
    prihead['TIME-OBS'] =  (expo_info[1], 'UT time of start of the 1st obs (hh:mm:ss)')
    prihead['EXPSTART'] = (expo_info[2], 'Start time of the first exposure, MJD')
    prihead['EXPEND'] = (expo_info[3], 'End time of the last exposure, MJD')
    prihead['WAVEUNIT'] = ('vacuum', 'According to COS IHB for C25')
    prihead['EQUINOX'] = 2000.0

    if expo_info[4] != 0 and expo_info[5] != 0: 
        prihead['EXPTIME'] = (expo_info[4]+expo_info[5], 'Total effective exposure time in sec')
        prihead['FILTER'] = 'MULTI'  # if spectrum combines more than one grism, or the filter/grism if only one
        prihead['EXPTIME1'] = (expo_info[4], 'Total Effective exposure time in FILTER1 in sec')
        prihead['EXPTIME2'] = (expo_info[5], 'Total effective exposure time in FILTER2 in sec')
        prihead['FILTER1'] = 'G130M' # the first filter/grism if combines more than one']
        prihead['FILTER2'] = 'G160M' # the second filter/grism....etc.']
    else:
        if expo_info[4] != 0 and expo_info[5] == 0: # has G130M but not G160M 
            prihead['EXPTIME'] = (expo_info[4], 'Total effective exposure time in sec')
            prihead['FILTER'] = 'G130M' # the first filter/grism if combines more than one']
        elif expo_info[4] == 0 and expo_info[5] !=0:  # has G160M but not G130M 
            prihead['EXPTIME'] = (expo_info[5], 'Total effective exposure time in sec')
            prihead['FILTER'] = 'G160M' # the first filter/grism if combines more than one']

    prihead['HLSPNAME'] = ('COS-GAL', 'HLSP product name')
    prihead['HLSPLEAD'] = 'Yong Zheng'
    prihead['HISTORY'] = 'The QSO spectrum is coadded by HLSA(Peeples+2017).'
    prihead['HISTORY'] = 'The continuum normalization is performed on %s.'%(todate)
    prihead['HISTORY'] = 'with Linetools software (Prochaska+2016).'
    
    del prihead['CREATOR'] 
    del prihead['HIERARCH TIMESTAMP']
    del prihead['LEGACY']

    return prihead

def save_fullspec(lt_spec, target_info, category, has_continuum=False, filedir=''):
    '''
    Save the full continuum-normalized spectra to a fits file
    
    lt_spec: linetools.spectra.xspectrum1d.XSpectrum1D object
    target_info: a dict object that contains ('NAME', 'z', 'RA', 'DEC', 'l', 'b', 'S/N', 'DATAFILE')
    has_continuum: if the XSpectrum1D.fit_continuum is performed, this should be True.
    '''
   
    # save the data to fits, first create Primary header 
    prihead = create_primary_header_UV(target_info, category)
    if type(prihead) == bool:
        return False
    prihdu = fits.PrimaryHDU(header=prihead)
    
    # then, create the data array 
    continuum, normflux, normsig = get_normspec(lt_spec, has_continuum=has_continuum)
 
    col_wave = fits.Column(name='WAVE', format='D', array=lt_spec.wavelength.value)
    col_flux = fits.Column(name='FLUX', format='D', array=lt_spec.flux.value)
    col_sig = fits.Column(name='ERROR', format='D', array=lt_spec.sig.value)
    col_cont = fits.Column(name='CONTINUUM', format='D', array=continuum)
    col_normflux = fits.Column(name='NORMFLUX', format='D', array=normflux)
    col_normsig = fits.Column(name='NORMERR', format='D', array=normsig)
    col_arrs = [col_wave, col_flux, col_sig, col_cont, col_normflux, col_normsig]

    tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(col_arrs))
    tbhdu.header['TUNIT1'] = 'Angstroms' # (or whatever unit the "WAVE" column is).
    tbhdu.header['TUNIT2'] = 'ergs/s/cm^2/Ang' # for the flux array or whatever unit that column is, etc. for the other columns.
    tbhdu.header['TUNIT3'] = 'ergs/s/cm^2/Ang' # for the error array
    tbhdu.header['TUNIT4'] = 'ergs/s/cm^2/Ang' # for the continuum array
    # tbhdu.header['WAVEUNIT'] = '' # "air" or "vacuum" depending on what the wavelengths are measured in
    
    ## now put everything together, and write to fits
    thdulist = fits.HDUList([prihdu, tbhdu])

    ## name the file in HLSP format 
    if target_info['Grating'] == 'FUVM': grating = 'g130m-g160m'
    else: grating = target_info['Grating']
    filename = '%s/hlsp_cos-gal_hst_cos_%s_%s_v1_fullspec.fits.gz'%(filedir,
                                                                    target_info['NAME'].lower(),
                                                                    grating.lower())
    thdulist.writeto(filename, clobber=True)

    return filename

def save_linespec(lt_spec, target_info, category, has_continuum=False, filedir='', 
                  line='none', do_redshift=False, velwidth=1000):
    '''
    Save the each continuum-normalized individual line to a fits file
    
    lt_spec: linetools.spectra.xspectrum1d.XSpectrum1D object
    target_info: a dict object that contains ('NAME', 'z', 'RA', 'DEC', 'l', 'b', 'S/N', 'DATAFILE')
    has_continuum: if the XSpectrum1D.fit_continuum is performed, this should be True.
    '''
    # read the line library, the accepted line name, wavelength, (and oscillator strength)
    # line_info = read_linelibrary(lines=line, doprint=False)
    from yzSpec.find_line_data import find_line_data
    line_info = find_line_data(line)
    if len(line_info) == 0: 
        logger.info('Cannot find atomic data for %s'%(line))
        return False  # can't find this line in the library

    # save the data to fits, first create Primary header 
    prihead = create_primary_header_UV(target_info, category)
    if type(prihead) == bool:
        return False

    prihead['LINE'] = line
    prihead['LAMBDA'] = (line_info['wave'], 'Vacuum wavelength in Angstrom; Morton (2003)')
    prihead['FVAL'] = (line_info['fval'], 'Oscillator Strength; Morton (2003)')
    prihdu = fits.PrimaryHDU(header=prihead)

    # then, create the data array 
    wave = lt_spec.wavelength.value
    flux = lt_spec.flux.value
    sig = lt_spec.sig.value
    continuum, normflux, normsig = get_normspec(lt_spec, has_continuum=has_continuum) 

    rest_lambda = line_info['wave']
    if do_redshift == True: obs_lambda = rest_lambda*(1+target_info['z'])
    else: obs_lambda = rest_lambda # MW line

    # velocity vector 
    velocity = (wave - obs_lambda)/obs_lambda*(speed_of_light_ms.value/1e3)

    # slice +/- 1000 km/s from the line center, obs_lambda
    wave = lt_spec.wavelength.value
    lf_lambda = obs_lambda - velwidth/(speed_of_light_ms.value/1e3)*obs_lambda
    rt_lambda = obs_lambda + velwidth/(speed_of_light_ms.value/1e3)*obs_lambda
    slice_ind = np.all([wave>=lf_lambda, wave<=rt_lambda], axis=0)

    # in case there is very little data exist, we don't use it
    if wave[slice_ind].size < 3.: 
        logger.info('No enough data exist for %s'%(line))
        return False
    
    if velocity[slice_ind][0] > 350 or velocity[slice_ind][-1] < -350:
        logger.info('No enough data exist within [-400, 400] for %s'%(line))
        return False

    # for enough data, let's save the line spec
    col_wave = fits.Column(name='WAVE', format='D', array=lt_spec.wavelength.value[slice_ind])
    col_flux = fits.Column(name='FLUX', format='D', array=lt_spec.flux.value[slice_ind])
    col_sig = fits.Column(name='ERROR', format='D', array=lt_spec.sig.value[slice_ind])
    col_cont = fits.Column(name='CONTINUUM', format='D', array=continuum[slice_ind])
    col_normflux = fits.Column(name='NORMFLUX', format='D', array=normflux[slice_ind])
    col_normsig = fits.Column(name='NORMERR', format='D', array=normsig[slice_ind])
    col_velo = fits.Column(name='VELOCITY', format='D', array=velocity[slice_ind])
    col_arrs = [col_wave, col_flux, col_sig, col_cont, col_normflux, col_normsig, col_velo]

    tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(col_arrs))
    tbhdu.header['TUNIT1'] = 'Angstroms' # (or whatever unit the "WAVE" column is).
    tbhdu.header['TUNIT2'] = 'ergs/s/cm^2/Ang' # for the flux array or whatever unit that column is, etc. for the other columns.
    tbhdu.header['TUNIT3'] = 'ergs/s/cm^2/Ang' # for the error array
    tbhdu.header['TUNIT4'] = 'ergs/s/cm^2/Ang' # for the continuum array
    tbhdu.header['TUNIT7'] = 'km/s' # for the velocity array
    # tbhdu.header['WAVEUNIT'] = '' # "air" or "vacuum" depending on what the wavelengths are measured in
    tbhdu.header['VELFRAME'] = 'Heliocentric'  # still need to be double checked. as of 02/22/2018, yz

    ## now put everything together, and write to fits
    thdulist = fits.HDUList([prihdu, tbhdu])

    ## name the file in HLSP format 
    if target_info['Grating'] == 'FUVM': grating = 'g130m-g160m'
    else: grating = target_info['Grating']
    filename = '%s/hlsp_cos-gal_hst_cos_%s_%s_v1_%s-spec.fits.gz'%(filedir,
                                                                   target_info['NAME'].lower(),
                                                                   grating.lower(), 
                                                                   line_info['hlsp-name'].lower())
    thdulist.writeto(filename, clobber=True)
    return filename


def save_spec(lt_spec, target_info, has_continuum=False, line='none', filedir='', do_redshift=False, velwidth=1000):
    '''
    Save the spectra into fits. Could save the whole spectra, or a segment near a specific line. 
   
    lt_spec: linetools.spectra.xspectrum1d.XSpectrum1D object
    target_info: a dict object that contains ('NAME', 'z', 'RA', 'DEC', 'l', 'b', 'S/N', 'DATAFILE')
    has_continuum: if the XSpectrum1D.fit_continuum is performed, this should be True.
    line: default to 'none' if want to save the whole lt_spec; 
          otherwise set to, e.g., line='SII1250', to save a segment of the spec within 
          +/-1000 km/s of the line centroid. In this case, the fits header primary hdu
          would have the line info, i.e., line name, rest wavelength, and f value. 
    '''

    wave = lt_spec.wavelength.value
    flux = lt_spec.flux.value
    sig = lt_spec.sig.value

    if has_continuum == True:
        # to avoid 0 in the denominator
 
        continuum = lt_spec.co.value
        indnan = np.where(continuum == 0.)
        continuum[indnan] = 1.

        # normlized flux
        normflux = flux / continuum
        normflux[indnan] = np.nan

        # normlized error
        normsig = sig / continuum
        normsig[indnan] = np.nan  # still keep the nan since we want the velocity array to be continous

    else:
        continuum = np.zeros(wave.size)+np.nan
        normflux = np.zeros(wave.size)+np.nan
        normsig = np.zeros(wave.size)+np.nan

    if line != 'none': 
        # only save a specific line
        # read the line library, the accepted line name, wavelength, (and oscillator strength)
        line_info = read_linelibrary(lines=line, doprint=False)
        rest_lambda = line_info[1][0]
        if do_redshift == True:
            obs_lambda = rest_lambda*(1+target_info['z']) 
        else: 
            obs_lambda = rest_lambda*(1+0.0) # MW line
  
        # velocity vector 
        velocity = (wave-obs_lambda)/obs_lambda*(speed_of_light_ms.value/1e3)
    
        # slice +/- 1000 km/s from the line center, obs_lambda
        lf_lambda = obs_lambda-velwidth/(speed_of_light_ms.value/1e3)*obs_lambda
        rt_lambda = obs_lambda+velwidth/(speed_of_light_ms.value/1e3)*obs_lambda
        slice_ind = np.where(np.all([wave>=lf_lambda, wave<=rt_lambda], axis=0)==True)[0]
    else: 
        slice_ind = np.arange(wave.size) 

    # check if the data exist. 
    if wave[slice_ind].size <= 3: 
        return 'none'  
    else:  
        # save the data; first generate the header
        todate = str(datetime.datetime.now())
        hdulist = fits.open(target_info['DATAFILE'])
        head0 = hdulist[0].header
        head0['HISTORY'] = 'YZ added continuum using linetools. %s'%(todate)
   
        col_wave = fits.Column(name='WAVE', format='D', array=wave[slice_ind])
        col_flux = fits.Column(name='FLUX', format='D', array=flux[slice_ind])
        col_sig = fits.Column(name='ERROR', format='D', array=sig[slice_ind])
        col_cont = fits.Column(name='CONTINUUM', format='D', array=continuum[slice_ind])
        col_normflux = fits.Column(name='NORMFLUX', format='D', array=normflux[slice_ind])
        col_normsig = fits.Column(name='NORMERR', format='D', array=normsig[slice_ind])
        col_arrs = [col_wave, col_flux, col_sig, col_cont, col_normflux, col_normsig]
       
        # consider the case if just one line is needed 
        if line != 'none': 
            head0['LINE'] = line_info[0][0]
            head0['LAMBDA'] = line_info[1][0]
            head0['FVAL'] = line_info[2][0]
            head0['LINEREF'] = line_info[3]
            
            col_velo = fits.Column(name='VELOCITY', format='D', array=velocity[slice_ind])
            col_arrs.append(col_velo)

            filename = '%s/%s_%s.fits.gz'%(filedir, target_info['NAME'], line.replace(' ', ''))
        else: 
            filename = '%s/%s.fits.gz'%(filedir, target_info['NAME'])
    
        cols = fits.ColDefs(col_arrs)
        tbhdu = fits.BinTableHDU.from_columns(cols)
        prihdu = fits.PrimaryHDU(header=head0)
        thdulist = fits.HDUList([prihdu, tbhdu])
        thdulist.writeto(filename, clobber=True)
        
        return filename
