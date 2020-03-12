import numpy as np
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('TkAgg')  # backend for graphical interfaces
from yzSpec.yes_or_no import yes_or_no
from yzSpec.read_knots import read_knots

import os
homedir = os.path.expanduser('~')

from tqdm import *

import logging 
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def linetools_UVHI(target, lines='All', extract_HI=False):
    '''
    Call the linetools continuum fitting function to fit the spectra 
    for the input target. It will also find the corresponding HI 21cm spectra.
   
    '''
   
    from yzSpec.find_targetinfo import find_targetinfo
    category, target_info = find_targetinfo(target)

    if category == 'none':
        logger.info(target+' is not in QSOALS, GALAXY, STAR_EARLY.')
        return False
    elif len(target_info) == 0: 
        logger.info(target+' does not have either G130M or G160M observation.')
        return False
    else: 
        # OK, now fit continuum 
        ## (if the QSO has been processed before, the knots should be saved in knots/) 
        knotdir = '%s/Dropbox/HSLA_Feb16/code/knots'%(homedir)
        knotname = '%s/%s_coadd_%s_knots.jsn'%(knotdir, target_info['NAME'], target_info['Grating'])
        print(knotname)
        knotlist = read_knots(knotname)
 
        print(target_info['DATAFILE']) 
        logger.info('Start fitting %s ...'%(target)) 
        import linetools.spectra.xspectrum1d as lsx 
        lt_spec = lsx.XSpectrum1D.from_file(target_info['DATAFILE'])
        try:
            lt_spec.fit_continuum(knots=knotlist)
        except RuntimeError:
            logger.info('RuntimeError: Problem generating continuum spline knots. Use random knots.')
            tempknots = '%s/Dropbox/HSLA_Feb16/code/knots/tempknots.jsn'%(homedir)
            knotlist = read_knots(tempknots)
            lt_spec.fit_continuum(knots=knotlist)
   
        from shutil import copyfile 
        copyfile('_knots.jsn', knotname)
        os.remove('_knots.jsn')

        ## check if lt_spec has the same data size in each array, sanity check 
        ## sanitycheck = 
        '''
        ### this is the old way of doing the spectral saving.  
        ## save the whole continuum 
        savedir = homedir+'/Dropbox/HSLA_Feb16/QSOSpec_YZ/'+category
        filedir = savedir+'/'+target_info['NAME']
        if os.path.isdir(filedir) == False:
            os.makedirs(filedir)
            os.makedirs(filedir+'/lines')
            os.makedirs(filedir+'/figs')

        from yzSpec.save_spec import save_spec
        save_spec(lt_spec, target_info, has_continuum=True, filedir=filedir)

        from yzSpec.find_line_data import import_lines
        lines = import_lines()
        for iline in lines:
            savename = save_spec(lt_spec, target_info, has_continuum=True, 
                                 line=iline, filedir=filedir+'/lines')

        from yzSpec.plot_spec import stack_spec
        stack_spec(target_info, filedir+'/lines')
        ''' 
        
        
        ## add on 02/20/2018, to add HLSP-formatted product 
        hlsp_savedir = homedir+'/Dropbox/HSLA_Feb16/QSOSpec_YZ/hlsp_qsoals'
        hlsp_filedir = hlsp_savedir+'/'+target_info['NAME'].lower()
        if os.path.isdir(hlsp_filedir) == False:
            os.makedirs(hlsp_filedir)
            os.makedirs(hlsp_filedir+'/linedata_uv')
            os.makedirs(hlsp_filedir+'/linedata_21cm')
            os.makedirs(hlsp_filedir+'/figs')
        
        logger.info('Saving, ploting, stacking spec for %s ...'%(target)) 
        import time
        t1 = time.time()
        
        from yzSpec.save_spec import save_fullspec
        name = save_fullspec(lt_spec, target_info, category, 
                             has_continuum=True, filedir=hlsp_filedir)
        if type(name) == bool: 
            logger.info('%s does not have full spectra, odd. '%(target_info['NAME']))
            return False
        
        t2 = time.time()
        logger.info('Saving the whole spec: %.1f seconds.'%(t2-t1))
        
        ### ADD by YZ, 02/21/2021, to add HLSP-formmatted product
        from yzSpec.save_spec import save_linespec
        from yzSpec.plot_spec import plot_OneLine
        from yzSpec.find_line_data import import_lines
        if lines[0].lower() == 'all': 
            lines = import_lines()
            ## basically rewrite the linedata_uv directory 
            for fn in os.listdir(hlsp_filedir+'/linedata_uv'):    
                os.remove(os.path.join(hlsp_filedir, 'linedata_uv', fn))
            for fn in os.listdir(hlsp_filedir+'/figs'):
                if fn.lower().endswith('spec.pdf'):  # to replot the UV line as well
                    os.remove(os.path.join(hlsp_filedir, 'figs', fn))
                
        logger.info('Slicing/plotting each line: ... ')
        for iline in tqdm(lines):
            savename = save_linespec(lt_spec, target_info, category, 
                                     has_continuum=True, line=iline, 
                                     filedir=hlsp_filedir+'/linedata_uv')
            if type(savename) != bool:
                plot_OneLine(target_info, savename, filedir=hlsp_filedir+'/figs')
        t3 = time.time()
 
        ## check if HI data exist; if not, extract it from HI4PI
        if extract_HI == True:
            from yzGALFAHI.extract_HI21cm import save_HIspec_fits
            savedir = hlsp_filedir+'/linedata_21cm'
            lab_good = save_HIspec_fits(target_info, beam=1.0, savedir=savedir, observation='LAB',
                                        datadir='/Users/Yong/Dropbox/databucket/LAB/labh_glue.fits')
            hi4_good = save_HIspec_fits(target_info, beam=1.0, savedir=savedir, observation='HI4PI',
                                        datadir='/Volumes/YongData2TB/HI4PI')
            ghi_good = save_HIspec_fits(target_info, beam=1.0, savedir=savedir, observation='GALFA-HI',
                                        datadir='/Volumes/YongData2TB/GALFAHI_DR2/RC5/Wide')
            hi4_good2 = save_HIspec_fits(target_info, beam=16.2/60, savedir=savedir, observation='HI4PI',
                                         datadir='/Volumes/YongData2TB/HI4PI')
            ghi_good2 = save_HIspec_fits(target_info, beam=4/60., savedir=savedir, observation='GALFA-HI',
                                         datadir='/Volumes/YongData2TB/GALFAHI_DR2/RC5/Wide')
            if lab_good == False or hi4_good == False or hi4_good2 == False:
                logger.info('BAD HI::::: ', target_info['NAME'])

        t4 = time.time()
        logger.info('Finding HI 21cm: %.1f seconds.'%(t4-t3))
        
        ## stacking the spectra
        from yzSpec.plot_spec import stack_allline
        stack_allline(target_info, hlsp_filedir, pltrange=[-400, 400], 
                      savedir=hlsp_filedir+'/figs', plt_HI=True)

        t5 = time.time()
        logger.info('Stacking: %.1f seconds.'%(t5-t4))
        
        return target_info['NAME']
