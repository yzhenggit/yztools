from yztools.mhalo2mstar import *

def calc_percentile(para_range, ndecimals=2):
    import numpy as np 
    percentile = np.percentile(para_range, [16, 50, 84])
    best_val = np.around(percentile[1], decimals=ndecimals)
    sig_lo = -np.around(np.diff(percentile)[0], decimals=ndecimals)
    sig_hi = np.around(np.diff(percentile)[1], decimals=ndecimals)
    # print(best_val, sig_lo, sig_hi)
    
    return best_val, sig_lo, sig_hi

def mstar2mhalo(logmstar, sig_logmstar=0., ref='B13-GK14', do_print=False, Ntrials=4000):
    """
    ref: 'B13' = Behroozi13
         'B13-GK14'=Behoorzi13 with lower end modified with Garrison-Kimmel+2014
         'M10' = Moster+2010
         'M13' = Moster+2013
         'GK14' = Garrison-Kimmel+2014
         'GK17' = Garrison-Kimmel+2017
         'M21' = Munshi+2021; for logM200c>11.514, use Behroozi+2013

    Ntrials: default to 4000 trials of MCMC 

    Example: 
    $ python mstar2mhalo.py logmstar sig_logmstar B13-GK14

    YZ. 05/09/2022
    Updated on 06/12, adding error bars 
    Updated on 07/28, added M21 method
    """

    import numpy as np
    from scipy import interpolate
    logmhalo_arr = np.arange(8,16,0.001)

    # interpolate this to mstar
    if ref.upper() == 'B13':
        logmstar_arr = mhalo2mstar_behroozi13(logmhalo_arr)
    elif ref.upper() in ['B13-GK14', 'B13_GK14']:
        logmstar_arr = mhalo2mstar_behroozi13_gk14modified(logmhalo_arr)
    elif ref.upper() in ['M21']: 
        #if do_print == True: 
        logmstar_arr = m200c2mstar_munshi21(logmhalo_arr)
    elif ref.upper() == 'M10':
        logmstar_arr = mhalo2mstar_moster10(logmhalo_arr)
    elif ref.upper() == 'M13':
        logmstar_arr = mhalo2mstar_moster13(logmhalo_arr)
    elif ref.upper() == 'GK14':
        logmstar_arr = mhalo2mstar_gk14(logmhalo_arr)
    elif ref.upper() == 'GK17':
        logmstar_arr = mhalo2mstar_gk17(logmhalo_arr)
    else:
        print('Cannot find the corresponding ref equation, please check')
        import sys
        sys.exit()

    # interpolate the result
    func = interpolate.interp1d(logmstar_arr, logmhalo_arr)

    if sig_logmstar == 0.: 
        best_logmhalo = np.around(float(func(logmstar)), decimals=2)
        sig_lo, sig_hi = 0., 0.
    else: # if errors in logmstar is given, we run mcmc to estimate the range
        # e.g., https://astrofrog.github.io/py4sci/_static/Practice%20Problem%20-%20Monte-Carlo%20Error%20Propagation%20-%20Sample%20Solution.html
        logmstar_range = np.random.normal(logmstar, sig_logmstar, Ntrials) 
        logmhalo_range = np.zeros(Ntrials)+np.nan
        for i in range(Ntrials): 
            logmhalo_range[i] = float(func(logmstar_range[i]))       
        
        # calculate the 16, 50, and 84 percentiles 
        best_logmhalo, sig_lo, sig_hi = calc_percentile(logmhalo_range)
        
    if do_print==True:
        print('Input: logMstar={}+/-{}'.format(logmstar, sig_logmstar))
        print('Output: logMhalo={}+{}{}'.format(best_logmhalo, sig_hi, sig_lo))

    return best_logmhalo, sig_lo, sig_hi

if __name__ == "__main__":
    import sys
    import numpy as np
    logmstar = float(sys.argv[1])
    sig_logmstar = float(sys.argv[2])
    ref = sys.argv[3]
    mstar2mhalo(logmstar, sig_logmstar=sig_logmstar, ref=ref, do_print=True)
