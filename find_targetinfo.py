import os
import numpy as np
from astropy.table import Table
homedir = os.path.expanduser('~')

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def locate_category(target):
    '''
    Find the category in HSLA for the input target. For now, it can
    only find target in QSOALS, GALAXY, or STAR_EARLY. 
 
    Return: category name (if founded)
    '''

    category = 'none'
    thislist = 'none'
    for ltag in ['QSOALS', 'GALAXY', 'STAR_EARLY', 'STAR_POST_AGB']:
        lname = '%s/Dropbox/HSLA_Feb16/targetlists/%s_sample.txt'%(homedir, ltag)
        objs = Table.read(lname, format='ascii')
        if target in objs['ID']: 
            category = ltag
            thislist = lname
            break
    return category, thislist

def find_targetinfo(target):
    '''
    Provided the target name, this func would find which category (QSOALS, 
    GALAXY, or STAR_EARLY) it belongs to, and the gratings/coordinates info
    for this target. 
    '''
 
    category, objlist = locate_category(target)
    if category == 'none': 
        target_info = dict() 
    else: 
        objs = Table.read(objlist, format='ascii')
        ind = np.where(target == objs['ID'])[0]
        g130m = objs['G130M'][ind][0]
        g160m = objs['G160M'][ind][0]

        if g130m == 1: 
            if g160m == 1: gtag = 'FUVM'
            else: gtag = 'G130M'
        else: 
            if g160m == 1: gtag = 'G160M'
            else: gtag = 'none'

        # 101518, Yong, add this following line for only those 7 targets that HSLA store G130M/G160M wrong
        # 2XMM-J141348.3+440014/ H1821+643/ MRK304/ PG1126-041/ PG1435-067/ QSO-B1124+271/ SDSSJ085259.22+031320.6
        # gtag = 'G160M'
       
        datadir = homedir+'/Dropbox/HSLA_Feb16/datapile/'+category
        filename = '%s/%s/%s_coadd_%s_final_all.fits.gz'%(datadir, target, target, gtag) 
        if os.path.isfile(filename) == False: filename='none'
        
        if objs['z'][ind][0] == '.': redshift = np.nan
        else: redshift = float(objs['z'][ind][0])
        target_info = {'NAME': target, 
                       'z': redshift, 
                       'RA': float(objs['RA'][ind][0]), 
                       'DEC': float(objs['DEC'][ind][0]), 
                       'l': float(objs['l'][ind][0]),  
                       'b': float(objs['b'][ind][0]), 
                       'S/N': float(objs['S/N'][ind][0]),
                       'Grating': gtag, 
                       'DATAFILE': filename
                      }
    return category, target_info

