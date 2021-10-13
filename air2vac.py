import astropy.units as u
import numpy as np

def air2vac_Edlen66(lambda_air):
    """
    convert air/obs wavelength to vacuum wavelength using formula from Edlen1966.
    Tested, same as using
    from pypeit.core.wave import airtovac

    input:
    lambda_air: input air wavelength arrays

    output:
    the corresponding lambda in vaccum for the input air wave
    """

    #let's make a func from vac to air based on equation 1 from Edlen66 and Hanuschik03
    lambda_air = np.asarray(lambda_air)

    wmin = np.nanmin(lambda_air)-5
    wmax = np.nanmax(lambda_air)+5
    dw = 0.0001
    wave_vac_arr = np.mgrid[wmin:wmax:dw] # in unit of Ang
    sig = 1/(wave_vac_arr*u.AA.to(u.um)) # sigma is 1/lambda where lambda is in um
    n = 1 + (1e-8)*(8342.13+ 2406030/(130-sig**2) + 15997/(38.9-sig**2))
    wave_air_arr = wave_vac_arr/n

    from scipy.interpolate import interp1d
    func_air2vac = interp1d(wave_air_arr, wave_vac_arr)

    # now get the vac wavelength we want
    lambda_vac = func_air2vac(np.asarray(lambda_air))
    return np.around(lambda_vac, decimals=3)

if __name__ == "__main__":
    import sys
    import numpy as np
    lambda_air = np.float(sys.argv[1])
    lambda_vac = air2vac_Edlen66(lambda_air)
    print("(Edlen1966)")
    print("Input:  air wave = {:.3f}".format(lambda_air))
    print("Output: vac wave = {:.3f}".format(lambda_vac))
