from yztools.mhalo2mstar import *

def mstar2mhalo(log_mstar, ref='B13', do_print=False):
    """
    ref: 'B13' = Behroozi13
         'B13-GK14'=Behoorzi13 with lower end modified with Garrison-Kimmel+2014
         'M10' = Moster+2010
         'M13' = Moster+2013
         'GK14' = Garrison-Kimmel+2014
         'GK17' = Garrison-Kimmel+2017

    YZ. 05/09/2022
    """

    import numpy as np
    log_mhalo_arr = np.arange(8,15,0.001)

    # interpolate this to mstar
    if ref.upper() == 'B13':
        log_mstar_arr = mhalo2mstar_behroozi13(log_mhalo_arr)
    elif ref.upper() in ['B13-GK14', 'B13_GK14']:
        log_mstar_arr = mhalo2mstar_behroozi13_gk14modified(log_mhalo_arr)
    elif ref.upper() == 'M10':
        log_mstar_arr = mhalo2mstar_moster10(log_mhalo_arr)
    elif ref.upper() == 'M13':
        log_mstar_arr = mhalo2mstar_moster13(log_mhalo_arr)
    elif ref.upper() == 'GK14':
        log_mstar_arr = mhalo2mstar_gk14(log_mhalo_arr)
    elif ref.upper() == 'GK17':
        log_mstar_arr = mhalo2mstar_gk17(log_mhalo_arr)
    else:
        print('Cannot find the corresponding ref equation, please check')
        import sys
        sys.exit()

    # interpolate the result
    from scipy import interpolate
    func = interpolate.interp1d(log_mstar_arr, log_mhalo_arr)
    log_mhalo = func(log_mstar)

    if do_print==True:
        print('Input: logMstar={:.2f}'.format(log_mstar))
        print('Output: logMhalo={:.2f}'.format(log_mhalo))

    return log_mhalo

if __name__ == "__main__":
    import sys
    import numpy as np
    log_mstar = float(sys.argv[1])
    ref = sys.argv[2]
    mstar2mhalo(log_mstar, ref=ref, do_print=True)
