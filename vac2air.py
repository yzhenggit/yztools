import astropy.units as u
import numpy as np

def vac2air_Edlen66(lambda_vac):
    """
    convert vacuum wavelength to air/obs wavelength using formula from Edlen1966.
    Tested, same as using
    from pypeit.core.wave import airtovac

    input:
    lambda_vac: input vaccum wavelength arrays

    output:
    the corresponding lambda in obs/air for the input vac wave
    """

    #let's make a func from vac to air based on equation 1 from Edlen66 and Hanuschik03
    lambda_vac = np.asarray(lambda_vac)

    wmin = np.nanmin(lambda_vac)-5
    wmax = np.nanmax(lambda_vac)+5
    dw = 0.0001
    wave_vac_arr = np.mgrid[wmin:wmax:dw] # in unit of Ang
    sig = 1/(wave_vac_arr*u.AA.to(u.um)) # sigma is 1/lambda where lambda is in um
    n = 1 + (1e-8)*(8342.13+ 2406030/(130-sig**2) + 15997/(38.9-sig**2))
    wave_air_arr = wave_vac_arr/n

    from scipy.interpolate import interp1d
    func_vac2air = interp1d(wave_vac_arr, wave_air_arr)

    # now get the air wavelength we want
    lambda_air = func_vac2air(np.asarray(lambda_vac))
    return np.around(lambda_air, decimals=3)

if __name__ == "__main__":
    import sys
    import numpy as np
    lambda_vac = np.float(sys.argv[1])
    lambda_air = vac2air_Edlen66(lambda_vac)
    print("(Edlen1966)")
    print("Input:  vac wave = {:.3f}".format(lambda_vac))
    print("Output: air wave = {:.3f}".format(lambda_air))
