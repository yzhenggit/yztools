from yztools.mhalo2mstar import *

def mstar2mhalo_moster10(log_mstar):

    import numpy as np
    log_mhalo_arr = np.arange(8,15,0.001)

    # interpolate this to mstar
    log_mstar_arr = mhalo2mstar_moster10(log_mhalo_arr)

    # interpolate the result
    from scipy import interpolate
    func = interpolate.interp1d(log_mstar_arr, log_mhalo_arr)
    log_mhalo = func(log_mstar)

    return log_mhalo

def mstar2mhalo_moster13(log_mstar):

    import numpy as np
    log_mhalo_arr = np.arange(8,15,0.001)

    # interpolate this to mstar
    log_mstar_arr = mhalo2mstar_moster13(log_mhalo_arr)

    # interpolate the result
    from scipy import interpolate
    func = interpolate.interp1d(log_mstar_arr, log_mhalo_arr)
    log_mhalo = func(log_mstar)

    return log_mhalo

def mstar2mhalo_gk14(log_mstar):

    import numpy as np
    log_mhalo_arr = np.arange(8,15,0.001)

    # interpolate this to mstar
    log_mstar_arr = mhalo2mstar_gk14(log_mhalo_arr)

    # interpolate the result
    from scipy import interpolate
    func = interpolate.interp1d(log_mstar_arr, log_mhalo_arr)
    log_mhalo = func(log_mstar)

    return log_mhalo

def mstar2mhalo_dwarfcgm(log_mstar, do_print=False):
    """
    specifically derived for the Dwarf CGM project.
    At low mass end with logM*=6-8, we use Garrison-Kimmel+2017
    At high mass end with logM*>8, we use Behroozi+2013

    YZ. 05/09/2022
    """

    import numpy as np
    if log_mstar < 8:
        print('Using Garrison-Kimmel+2017...good for logMstar=[6,8]')
        mhalo2mstar_func = mhalo2mstar_gk17
        log_mhalo_arr = np.arange(9.8, 11,0.001)
    else:
        print('Using Behroozi+2013...good for logMstar=[7.25, 11.5]')
        mhalo2mstar_func = mhalo2mstar_behroozi13
        log_mhalo_arr = np.arange(10, 15, 0.001)

    log_mstar_arr = mhalo2mstar_func(log_mhalo_arr)

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
    mstar2mhalo_dwarfcgm(log_mstar, do_print=True)
