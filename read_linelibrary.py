import astropy.io.fits as fits
import sys, re
import numpy as np

from yzSpec.find_line_data import import_lines

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def read_linelibrary(lines='All', doprint=True):
    '''
    For the input lines, find the wavelength and fval. 
    '''

    # find the path, read in the line library 
    # xidl line library is really outdated I think, don't use this one. 
    import sys
    for ipath in sys.path:
        if 'GitRepo' in ipath:
            linepath = ipath
            break
    all_lin = fits.getdata(linepath+'/yzSpec/files/xidl_all_lin.fits')
    # Lines that covered in G130 M, should really expand this later on. 
    if lines == 'All':
        dolines = import_lines()
    else:
        # Do spec slicing for some particular lines
        # lines could be just one: lines = ['FeII 1142']
        #              or several: lines = ['SiII 1250', 'SiII 1253', 'SiII 1259']

        # reformat the lines, in case it is not in forms, like, 'FeII 1142', but 'Fe II 1142'
        dolines = []
        defaultlines = import_lines()
        if type(lines) is str: lines = [lines]  # make sure the for loop will work
        for il in range(len(lines)):
            eles = re.split('(\d+)', lines[il].replace(' ', ''))
            newline = '%s %s'%(eles[0], eles[1])
            if newline in defaultlines: dolines.append(newline)
            else: logger.info(lines[il]+' not in the right format. skip.')

    # read line wavelength and oscillator strength
    line_lambda = []
    line_fval = []
    liblines = []
    for i in range(len(dolines)):
        for j in range(len(all_lin)):
            if dolines[i] == all_lin[j][0]:
                line_lambda.append(all_lin[j][1])
                liblines.append(dolines[i])

                # # FeII 1144. Check emails exchange with X.  
                if liblines[-1] == 'FeII 1144': line_fval.append(0.083)
                else: line_fval.append(all_lin[j][2])
                
                break

    if doprint == True: logger.info("Found these in our library: ",liblines)
    line_ref = 'Morton (2003)'  
    return liblines, line_lambda, line_fval, line_ref

