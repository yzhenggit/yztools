def load_atomic_data_morton03(ion):
    
    ion = ion.replace(' ', '')
    ### relevant atomic data
    # based on Morton 2003, table 2
    # https://iopscience.iop.org/article/10.1086/377639/fulltext/
    
    ion_hlsp_names = {'CIV': 'c-iv',
                      'SiIV': 'si-iv',
                      'NV': 'n-v', 
                      'OVI': 'o-vi', 
                      'SiII': 'si-ii', 
                      'SiIII': 'si-iii', 
                      'CII': 'c-ii', 
                      #'CII*': 'c-ii-s'
                      'SII': 's-ii', 
                      'FeII': 'fe-ii', 
                      'PII': 'p-ii'}

    ### atomic data, morton 2003
    lines_names = {'SiIV': ['SiIV1393', 'SiIV1402'],
                   'CIV': ['CIV1548', 'CIV1550'],
                   'NV': ['NV1238', 'NV1242'], 
                   'OVI': ['OVI1031', 'OVI1037'], 
                   'SiII': ['SiII1190', 'SiII1193', 'SiII1260', 
                            'SiII1526', 'SiII1808'], 
                   'SiIII': ['SiIII1206'], 
                   'CII': ['CII1334'], 
                   # 'CII*':['CII*1335'] , 
                   'SII': ['SII1250', 'SII1253', 'SII1259'], 
                   }

    lines_lambda = {'SiIV': [1393.7602, 1402.7729],
                    'CIV': [1548.204, 1550.781],
                    'NV': [1238.821, 1242.804], 
                    'OVI': [1031.9261, 1037.6167], 
                    'SiII': [1190.4158, 1193.2897, 1260.4221, 1526.7070, 1808.0129], 
                    'SiIII': [1206.500], 
                    'CII': [1334.5323], 
                    #'CII*': 
                    'SII': [1250.578, 1253.805, 1259.518]
                    }

    lines_f = {'SiIV': [5.13E-01, 2.54E-01],
               'CIV': [1.899E-01, 9.475E-02],
               'NV': [1.560E-01, 7.770E-02], 
               'OVI': [1.325E-01, 6.580E-02], 
               'SiII': [2.92E-01, 5.82E-01, 1.18E+00, 1.33E-01, 2.08E-03], 
               'SiIII': [1.63E+00], 
               'CII': [1.28E-01], 
               'SII': [5.43E-03, 1.09E-02, 1.66E-02]
              }

    line_info = {'ion_name': [ion_hlsp_names[ion], ion_hlsp_names[ion],],
                 'line_name': lines_names[ion],
                 'lambda': lines_lambda[ion],
                 'f': lines_f[ion]}
    return line_info

