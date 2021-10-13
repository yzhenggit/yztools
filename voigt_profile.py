# History:
# 09/03/21, update voigt profile using functions from linetools, YZ.
# first version, based on Draine 2011 Chapter 6 & 9
# 03/03/2021, Yong Zheng

import numpy as np
from astropy import constants as const
import astropy.units as u
from astropy.table import Table

def lt_Voigt_profile(line, wave, logN, b, z):
    """
    line: e.g., 'SiIV 1393'
    wave_AA: wavelength array values in unit of AA
    """

    import astropy.units as u
    from linetools.spectralline import AbsLine
    line_comp = AbsLine(trans=line)

    line_comp.attrib['N'] = 10**logN*u.cm**(-2)
    line_comp.attrib['b'] = b*u.km/u.s
    line_comp.setz(z)

    line_comp_voigt = line_comp.generate_voigt(wave=wave)

    return line_comp_voigt

def lt_Voigt_LSF(wave, flux, instrument='COS', grating='G130M'):
    """
    line: e.g., 'SiIV 1393'
    wave: wavelength array values in unit of AA
    instrument: only for 'COS', 'STIS'
    grating: only for 'G130M', 'G160M', 'E140M'
    """

    from linetools.spectra.lsf import LSF

    voigt_flux = flux # line_comp_voigt.flux
    voigt_flux_lsf = np.zeros_like(voigt_flux)

    # first determine the LSF
    for ll,w0 in enumerate(wave):
        lolim_idx = np.max([0,ll-20])
        uplim_idx = np.min([len(wave)-1,ll+20])

        # decide which instrument configuration to use
        if (instrument == 'COS') & (grating == 'G130M'):
            instr_config = LSF({'name':'COS','grating':'G130M','life_position':'1'})
        elif (instrument == 'COS') & (grating == 'G160M'):
            instr_config = LSF({'name':'COS','grating':'G160M','life_position':'1'})
        else:
            print("I don't recognize this instrument/grating: %s/%s. stop."%(instrument, grating))
            return
        # NEED to add STIS here

        lsf_tab = instr_config.interpolate_to_wv_array(wave[lolim_idx:uplim_idx+1])
        # LSF, why do this?
        #if w0 < 1800.*u.AA:
        #    lsf_tab = instr_config.interpolate_to_wv_array(wave[lolim_idx:uplim_idx+1])
        #else:
        #    new_kernel_wv = np.min([line_spec_wave[lolim_idx:uplim_idx+1].value,
        #                        np.zeros_like(line_spec['Wave'][lolim_idx:uplim_idx+1].value)+1799.99],axis=0)
        #    lsf_tab = instr_config.interpolate_to_wv_array(new_kernel_wv*u.AA)

        # convolve
        voigt_flux_lsf[lolim_idx:uplim_idx+1] += voigt_flux[ll]*lsf_tab['kernel']

    return voigt_flux_lsf

def phi_nu_voigt(nu, sigv, r_ul, nu_ul):
    # nu, frequency array
    # sigv, thermal broadening
    # r_ul, see above
    # nu_ul, frequency at line center

    dv = 0.005*u.km/u.s # km/s
    v = np.mgrid[-5*sigv.value:5*sigv.value:dv.value]*u.km/u.s
    v = v.reshape(1, v.size)

    gauss_part = dv/sigv*np.exp(-v**2/2/sigv**2)/np.sqrt(2*np.pi)
    lorentz_part = 4*r_ul / (16*np.pi**2 * (nu-(1-v/const.c)*nu_ul)**2 + r_ul**2)

    phi_nu = np.sum(gauss_part * lorentz_part, axis=1).cgs
    return phi_nu

def tau_nu(f_lu, N_l, phi_nu):
    # f_lu, oscillator strength of the line
    # N_l column density of ion at the low levels, (which is where most ions are at in the ISM/CGM environment)
    # also this ignores the stimulated emission

    tau_nu = (np.pi*const.e.esu**2/const.m_e/const.c * f_lu * N_l * phi_nu).cgs
    return tau_nu

def I_nu_voigt(nu, sigv, r_ul, nu_ul, f_lu, N_l):
    # nu, frequency array
    # sigv, thermal broadening
    # r_ul, see above
    # nu_ul, frequency at line center
    # f_lu, oscillator strength of the line
    # N_l column density of ion at the low levels, (which is where most ions are at in the ISM/CGM environment)


    phi_nu = phi_nu_voigt(nu, sigv, r_ul, nu_ul)
    line_tau = tau_nu(f_lu, N_l, phi_nu)
    I_nu = np.exp(-line_tau)

    return I_nu


def voigt_profile(ionline, N, b, wave_halfwidth=5, dwave=0.005):
    # wave_halfwidth: how broad the line region you want to generate Voigt profile for, Angstrom
    # dwave: Angstrom, the wavelength width you want to sample the function
    # N: column density, unit of cm-2 # only value input though
    # b: doppler width, unit of km/s # only value input though

    from yztools.atom_wave_gamma_f import atom_wave_gamma_f
    lineres = atom_wave_gamma_f(ionline)
    line_gamma = lineres['gamma']
    line_wrest = lineres['wrest']
    line_fvalue = lineres['f']

    # add unit
    r_ul = line_gamma/u.s
    lambda_ul = line_wrest*u.Angstrom
    nu_ul = (const.c/lambda_ul).cgs
    f_lu = line_fvalue

    # wavelength array
    wave_wd = wave_halfwidth # Ang
    dw = dwave # Ang
    wave = np.mgrid[line_wrest-wave_wd:line_wrest+wave_wd:dw]*u.Angstrom

    nu = (const.c/wave).cgs
    nu = nu.reshape(nu.size, 1)

    from yztools.voigt_profile import I_nu_voigt
    sigv = b/np.sqrt(2)*u.km/u.s
    N_l = N/u.cm**2
    I_nu = I_nu_voigt(nu, sigv, r_ul, nu_ul, f_lu, N_l)

    return nu.flatten(), wave, I_nu
