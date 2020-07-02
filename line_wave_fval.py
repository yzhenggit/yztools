import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def line_wave_fval(line):
    from astropy.table import Table
    import numpy as np
    import re

    if line == 'FeII1144': 
        line = 'FeII1145'
    # line = 'SiIV1393'
    line_ele = re.split('(\d+)',line.replace(' ', ''))[0]
    line_wave = int(re.split('(\d+)',line.replace(' ', ''))[1])

    all_line_info = Table.read('/Users/Yong/Dropbox/atom_wave_gamma_f.dat', format='ascii')
    sub_info = all_line_info[all_line_info['Ion'] == line_ele]

    ind = np.argmin(abs(sub_info['Wavelength[A]'] - line_wave))

    line_acc_wave = sub_info['Wavelength[A]'][ind]
    line_fval = sub_info['f_value'][ind]

    line_info = {'wave': line_acc_wave,
                 'fval': line_fval,
                 'Ref': 'atom_wave_gamma_f.dat'}

    print('Input    ele  acc_wave fval')
    print(line, sub_info['Ion'][0], line_acc_wave, line_fval)

    return line_info

def import_lines():
    '''
    Load the lines of interest for COS-GAL
    '''

    defaultlines = ['PII 1152', 'FeII 1142', 'FeII 1143', 'FeII 1144', 'FeII 1608',
                    'CII 1334', 'CIV 1548', 'CIV 1550',
                    'SiII 1190', 'SiII 1193', 'SiII 1260', 'SiII 1304', 'SiII 1526',
                    'SiIII 1206', 'SiIV 1393', 'SiIV 1402',
                    'SII 1250', 'SII 1253', 'SII 1259', 'NV 1238', 'NV 1242']
    return defaultlines

def find_line_data_old(ionline):
    # this code has been abondoned. see find_line_data for a more extended lib.
    ## this needs to be more sophisticated
    ## Morton 2003
    ionline = ionline.replace(' ', '')
    line_struct = {'SiII1190': {'wave': 1190.4158, 'fval': 2.92E-01,
                                'Ref': 'Morton2003', 'hlsp-name': 'Si-II-1190'},
                   'SiII1193': {'wave': 1193.2897, 'fval': 5.82E-01,
                                'Ref': 'Morton2003', 'hlsp-name': 'Si-II-1193'},
                   'SiII1260': {'wave': 1260.4221, 'fval': 1.18E+00,
                                'Ref': 'Morton2003', 'hlsp-name': 'Si-II-1260'},
                   'SiII1304': {'wave': 1304.3702, 'fval': 8.63E-02,
                                'Ref': 'Morton2003', 'hlsp-name': 'Si-II-1304'},
                   # Linetools give f=0.127 for SiII1526, with ref of Morton2003, seems wrong
                   'SiII1526': {'wave': 1526.7070, 'fval': 1.33E-01,
                                'Ref': 'Morton2003', 'hlsp-name': 'Si-II-1526'},
                   'SiIII1206':{'wave': 1206.500,  'fval': 1.63E+00,
                                'Ref': 'Morton2003', 'hlsp-name': 'Si-III-1206'},
                   # Linetools use Verner1994 for SiIV doublets
                   'SiIV1393': {'wave': 1393.7602, 'fval': 5.13E-01,
                                'Ref': 'Morton2003', 'hlsp-name': 'Si-IV-1393'},
                   'SiIV1402': {'wave': 1402.7729, 'fval': 2.54E-01,
                                'Ref': 'Morton2003', 'hlsp-name': 'Si-IV-1402'},
                   'CII1334':  {'wave': 1334.5323, 'fval': 1.28E-01,
                                'Ref': 'Morton2003', 'hlsp-name': 'C-II-1334'},
                   # Linetools use Verner1994 for CIV doublets
                   'CIV1548':  {'wave': 1548.204,  'fval': 1.899E-01,
                                'Ref': 'Morton2003', 'hlsp-name': 'C-IV-1548'},
                   'CIV1550':  {'wave': 1550.781,  'fval': 9.475E-02,
                                'Ref': 'Morton2003', 'hlsp-name': 'C-IV-1550'},
                   'SII1250':  {'wave': 1250.578,  'fval': 5.43E-03,
                                'Ref': 'Morton2003', 'hlsp-name': 'S-II-1250'},
                   'SII1253':  {'wave': 1253.805,  'fval': 1.09E-02,
                                'Ref': 'Morton2003', 'hlsp-name': 'S-II-1253'},
                   'SII1259':  {'wave': 1259.518,  'fval': 1.66E-02,
                                'Ref': 'Morton2003', 'hlsp-name': 'S-II-1259'},
                   'FeII1142': {'wave': 1142.3656, 'fval': 4.01E-03,
                                'Ref': 'Morton2003', 'hlsp-name': 'Fe-II-1142'},
                   # there is another 1142 line in FeII, but seems like this one is the usual one
                   'FeII1143': {'wave': 1143.2260, 'fval': 1.92E-02,
                                'Ref': 'Morton2003', 'hlsp-name': 'Fe-II-1143'},
                   'FeII1144': {'wave': 1144.9379, 'fval': 8.30E-02,
                                'Ref': 'Morton2003', 'hlsp-name': 'Fe-II-1144'},
                   'FeII1608': {'wave': 1608.4511, 'fval': 5.77E-02,
                                'Ref': 'Morton2003', 'hlsp-name': 'Fe-II-1608'},
                   'OI1302':   {'wave': 1302.1685, 'fval': 4.80E-02,
                                'Ref': 'Morton2003', 'hlsp-name': 'O-I-1302'},
                   'OVI1031':  {'wave': 1031.9261, 'fval': 1.325E-01,
                                'Ref': 'Morton2003', 'hlsp-name': 'O-VI-1031'},
                   'OVI1037':  {'wave': 1037.6167, 'fval': 6.580E-02,
                                'Ref': 'Morton2003', 'hlsp-name': 'O-VI-1037'},
                   'PII1152':  {'wave': 1152.8180, 'fval': 2.45E-01,
                                'Ref': 'Morton2003', 'hlsp-name': 'P-II-1152'},
                   'NV1238':   {'wave': 1238.821,  'fval': 1.560E-01,
                                'Ref': 'Morton2003', 'hlsp-name': 'N-V-1238'},
                   'NV1242':   {'wave': 1242.804,  'fval': 7.770E-02,
                                'Ref': 'Morton2003', 'hlsp-name': 'N-V-1242'},
                   'MgII1239': {'wave': 1239.9253, 'fval': 6.330E-04,
  				                'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'Mg-II-1239'},
                   'MgII1240': {'wave': 1240.3947, 'fval': 3.550E-04,
				                'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'Mg-II-1240'},
                   'AlII1670': {'wave': 1670.7874, 'fval': 1.880E+00,
                                'Ref':'atom_wave_gamma_f.dat', 'hlsp-name': 'Al-II-1670'},
		           'SiII1808': {'wave': 1808.0130, 'fval': 2.186E-03,
                                'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'Si-II-1808'},
                'PII1532': {'wave': 1532.5330, 'fval': 7.610E-03,
                            'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'P-II-1532'},
                'FeII1901': {'wave': 1901.7730, 'fval': 1.009E-04,
                            'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'Fe-II-1901'},
                'FeII1611': {'wave': 1611.2005, 'fval': 1.360E-03,
                            'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'Fe-II-1611'},
                'FeII1608': {'wave': 1608.4511, 'fval': 5.800E-02,
                            'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'Fe-II-1608'},
                'FeII1588': {'wave': 1588.6876, 'fval': 1.48e-4,
                            'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'Fe-II-1588'},
                'FeII1144': {'wave': 1144.9379, 'fval': 1.060E-01,
                            'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'Fe-II-1144'},
                'FeII1143': {'wave': 1143.2260, 'fval': 1.770E-02,
                            'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'Fe-II-1143'},
                'FeII1142': {'wave': 1142.3656, 'fval': 4.200E-03,
                            'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'Fe-II-1142'},
                'HI1215': {'wave': 1215.6701, 'fval': 0.416400,
                           'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'H-I-1215'},
                'NV1238': {'wave': 1238.8210, 'fval': 1.570E-01,
                           'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'N-V-1238'},
                'NV1242': {'wave': 1242.8040,  'fval': 7.823E-02,
                           'Ref': 'atom_wave_gamma_f.dat', 'hlsp-name': 'N-V-1242'}

                  }

    try:
        line_info = line_struct[ionline]
    except KeyError:
        logger.info('KeyError: %s does not exist, or is not included in the library.'%(ionline))
        line_info = {}
    return line_info

if __name__ == "__main__":
    import sys
    line = sys.argv[1]
    line_info = line_wave_fval(line)
