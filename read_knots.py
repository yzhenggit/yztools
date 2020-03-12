import os

def read_knots(knotname):
    '''
    This is to read the **_knots.jsn file for the continuum fitting func.
    '''
    if os.path.isfile(knotname) is True:
        knotfile = open(knotname).readline().split('],')
        knotlist = []
        for iknot in knotfile:
            a = float(iknot.split(',')[0].replace('[', ''))
            b = float(iknot.split(',')[1].replace(']', ''))
            knotlist.append([a, b])
    else: knotlist = None  # the default

    return knotlist

