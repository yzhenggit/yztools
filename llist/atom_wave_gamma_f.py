from astropy.table import Table
import numpy as np

def atom_wave_gamma_f(ionline, do_print=False):
    # ionline in the form of 'HI1215', or 'SiIII1206', nospace

    # ion_tb = Table.read('/Users/Yong/Dropbox/atom_wave_gamma_f.dat', format='ascii')
    line_dir = '/Users/Yong/Dropbox/GitRepo/yztools'
    ion_tb = Table.read(line_dir+'/atom_wave_gamma_f.dat', format='ascii')

    ind = np.nan
    for i in range(len(ion_tb)):
        name = '%s%d'%(ion_tb['Ion'][i], ion_tb['Wavelength[A]'][i])
        if name == ionline:
            ind = i
            break
    if np.isfinite(ind):
        line_gamma = ion_tb['gamma'][ind]
        line_fvalue = ion_tb['f_value'][ind]
        line_wrest = ion_tb['Wavelength[A]'][ind]
        if do_print == True: 
            print(ion_tb[ind])
    else:
        print('Did not find this line: ', ionline)

    lineres = {'gamma': line_gamma,
                'f': line_fvalue,
                'wrest': line_wrest}
    return lineres

if __name__ == '__main__':
    import sys
    ionline = sys.argv[1]
    res = atom_wave_gamma_f(ionline, do_print=True)
