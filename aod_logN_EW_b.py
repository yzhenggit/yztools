# Have already double checked this code with COS-HALOs eqwrange.pro, yield
# consistent result of logNaod and EW. - Y. Zheng, UCB, 2019/07/25

def aod_logN_EW_b(wave_arr, flux_arr, err_arr, linename, vmin, vmax,
                  target_z=0.0, vel_arr=[]):
    """
    Calculate column density, EW, and b for an input line arrays

    Input:
    wave_arr: wavelength array, in unit of A
    flux_arr: normalized flux array
    err_arr: normalized error array
    linename: line to study, can be, e.g., 'SiII 1190' or 'SiII1190'
    vmin, vmax: the velocity range to make the integration
    target_z: the redshift of the target, default to 0 for MW stuff, if for
              other target, need to be target_z = vhelio_target/speed_of_light if
              wave_arr is in helio frame.
    vel_arr: if [], then we will calculate the vel_arr from the input wave_arr,
             and consider the input target_z; if has values (same structure
             as wave_arr), then we will use that directly.

    Output:
    A dict with keys of v_cent, sigma_v, b, N, Nerr, logN, logNerr, EW_mA,
    EWerr_MA

    How to use:
    from yztools.aod_logN_EW_b import aod_logN_EW_b
    > aod_logN_EW_b(wave, flux, err, linename, vmin, vmax, target_z=vhelio/speed_of_light)

    History:
    03/17/2020, put together into a module, YZ. Most of the stuff are from the COS
    Halos IDL package, see Lines/eqwrange.pro

    """
    import numpy as np

    # get line info
    from yztools.line_wave_fval import line_wave_fval
    lineinfo = line_wave_fval(linename.replace(' ', ''))
    linewave = lineinfo['wave']
    linefval = lineinfo['fval']
    linewave_at_z = linewave*(1+target_z)
    print('>> %s  lambda=%.4f  fval=%.4f'%(linename, linewave, linefval))
    print('>> Line shift to lambda=%.4f at z=%.4f'%(linewave_at_z, target_z))

    # if velocity array is imported, we use that directly, otherwise
    # we convert wavelength to velocity for the particular line
    if len(vel_arr) == 0:
        print("there is no imported velocity array, so I am going to calculate that from wave_arr")
        from astropy import constants as const
        import astropy.units as u
        vel_arr = ((wave_arr-linewave_at_z)/linewave_at_z * const.c).to(u.km/u.s).value # km/s

    # we only care about stuff within some velocity range as indicated by (vmin, vmax)
    vel_arr[np.isnan(vel_arr)] = -10000 # in case there are nan values, we don't include that
    indv = np.all([vel_arr>=vmin, vel_arr<=vmax], axis=0)
    vel_arr = vel_arr[indv]
    wave_arr = wave_arr[indv]
    flux_arr = flux_arr[indv]
    err_arr = err_arr[indv]

    # get the velocity and wave_arr intevals
    dv = vel_arr[1:]-vel_arr[:-1]
    dv = np.concatenate([dv, [dv[-1]]])
    dlambda = wave_arr[1:] - wave_arr[:-1]
    dlambda = np.concatenate([dlambda, [dlambda[-1]]])

    # get the column density using apparent optical depth method
    # line 145 and 163 in COS-Halos IDL package, Lines/eqwrange.pro
    flux_arr[flux_arr<=0.] = 0.0001
    tau_v = -np.log(flux_arr) # iflux_arr = I_obs / I_continuum; line 134
    N_v = 3.768e14*tau_v/linewave/linefval # atom cm-2/(kms-1); line 145
    N = np.fabs(np.sum(N_v*dv)) # line 163

    # get column density err_arror
    # line 160 in COS-Halos IDL package, Lines/eqwrange.pro
    tauerr_v = err_arr/flux_arr # line 158
    Nerr_v = 3.768e14*tauerr_v/linewave/linefval  # line 160
    Nerr = np.sqrt(np.sum((Nerr_v*dv)**2))  # line 164

    # get logN and sig_logN
    # line 172 in COS-Halos IDL package, Lines/eqwrange.pro
    logN = np.log10(N)
    logNerr = np.log10(N+Nerr)-logN # line 172

    # get centroid velocity, the velocity where cumulative tau reaches 50%
    # line 136 in COS-Halos IDL package, Lines/eqwrange.pro
    cum_tau_v = np.zeros(tau_v.size)
    for i in range(tau_v.size):
        cum_tau_v[i] = tau_v[:(i+1)].sum()/tau_v.sum()
    ind50 = np.argmin(np.fabs(cum_tau_v-0.5))
    v_cent = vel_arr[ind50]

    # get sigma v and line width from the profile
    # based on line 138-143; and page 2 in Heckman et al. 2002, ApJ, 577:691
    top_part = np.sum((vel_arr-v_cent)**2 * tau_v * dv)
    bot_part = np.sum(tau_v * dv)
    sigma_v = np.sqrt(top_part/bot_part)
    doppler_b = np.sqrt(2) * sigma_v

    # get Equivalent width over the same velocity range
    # See Eq. 1 and 5 in ISM review by Savage+1996
    ew_mA = np.sum(dlambda*(1-flux_arr))*1000 # mA, line 111 in YongIDL/Lines/eqwrange.pro
    ewerr_mA = np.sqrt(np.sum((dlambda*err_arr)**2))*1000 # mA, line 120 in YongIDL/Lines/eqwrange.pro

    res = {'v_cent': v_cent,
           'sigma_v': sigma_v,
           'b': doppler_b,
           'N': N,
           'Nerr': Nerr,
           'logN': logN,
           'logNerr': logNerr,
           'EW_mA': ew_mA,
           'EWerr_mA': ewerr_mA
           }

    print('*'*60)
    print(">> v range in vsys: %.1f %.1f (vmin, vmax) km/s"%(vmin, vmax))
    print(">> v_cent: %.2f km/s"%(v_cent))
    print(">> sigma_v: %.2f km/s"%(sigma_v))
    print(">> b   : %.2f km/s = sqrt(sigv)"%(doppler_b))
    print(">> EW  : %.1f,  %.1f,  %.1f sigma (EW, eEW, EW/eEW) mA"%(ew_mA, ewerr_mA, ew_mA/ewerr_mA))
    print('>> N   : %.2e cm-2'%(N))
    print(">> Nerr: %.2e cm-2"%(Nerr))
    print('>> logN: %.2f'%(logN))
    print('>> elogN: %.2f = np.log10(N+Nerr)-logN'%(logNerr))
    print(">> logN_2sig: %.2f = np.log10(N+2*Nerr)"%(np.log10(N+2*Nerr)))
    print('*'*60)

    return res
