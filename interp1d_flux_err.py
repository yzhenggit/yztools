import numpy as np

def interp1d_flux_err(interp_wave_grid, wave, flux, error, fill_values=np.nan): 
    """
    Interpolate spectra linearly, and do error proporgation properly 
    
    input: 
    interp_wave_grid: the wavelength grid you want to interpolate the spectra too 
    wave: original wave array 
    flux: original flux array 
    error: original error array
    fill_values: if interpolation is attempted beyond the wave range, fill in this value, 
              could be np.nan or 0 
              
    reutrn: 
    interp_flux_grid
    interp_err_grid
    
    History: 
    First written for M33/Keck DEIMOS project, in discussion with AS. 
    09/08/2021, YZ. UCB. 
    """

    interp_flux_grid = np.zeros(interp_wave_grid.size) 
    interp_err_grid = np.zeros(interp_wave_grid.size)
    
    in_range = np.all([interp_wave_grid<=wave.max(), interp_wave_grid>=wave.min()], axis=0)
    
    # fill the out_range with fill_values
    out_range = np.logical_not(in_range)
    if len(interp_flux_grid[out_range]) != 0: 
        interp_flux_grid[out_range] = fill_values
        interp_err_grid[out_range] = fill_values

    # Return the indices of the bins to which each value in input array belongs
    wave_grid = interp_wave_grid[in_range]
    
    bin_idx = np.digitize(wave_grid, wave, right=False)  # bins[i-1] <= x < bins[i]
    coeff_grid = (wave[bin_idx]-wave_grid)/(wave[bin_idx]-wave[bin_idx-1])

    # variance 
    var_grid = ((1-coeff_grid)*error[bin_idx])**2 + (coeff_grid*error[bin_idx-1])**2
    error_grid = np.sqrt(var_grid)
    interp_err_grid[in_range] = error_grid

    # flux (to cross match with the interp1d to double check)
    flux_grid = (1-coeff_grid)*flux[bin_idx] + (coeff_grid*flux[bin_idx-1])
    interp_flux_grid[in_range] = flux_grid

    # cross checked, all good. 
    from scipy.interpolate import interp1d
    interp1d_func = interp1d(wave, flux, kind='linear', bounds_error=False, fill_value=np.nan)
    interp1d_newflux = interp1d_func(interp_wave_grid)
    
    return interp1d_newflux, interp_err_grid
