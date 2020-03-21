
def reformat_spec_wave_flux_error(input_filename, output_filename):
    import astropy.io.fits as fits
    import numpy as np
    from astropy.table import Table

    from linetools.spectra import io as tio
    spec = tio.readspec(input_filename)
    wave = np.asarray(spec.wavelength)
    flux = np.asarray(spec.flux)
    error = np.asarray(spec.sig)

    data = Table([wave, flux, error], names=['WAVE', 'FLUX', 'ERROR'])
    data.write(output_filename, format='ascii', overwrite=True)

    return wave, flux, error

if __name__ == '__main__':
    import sys
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]

    res = reformat_spec_wave_flux_error(input_filename, output_filename)
