from astropy.io import fits
from astropy import units as u
from astropy import constants as c
import numpy as np
import os


def flux_planck_law(wavelength, temperature=4500*u.K, radius=2*u.R_sun, distance=4.32e+18*u.m):
    planck = (2.0 * c.h * c.c**2 / wavelength**5) / (np.exp(c.h * c.c / (wavelength * c.k_B * temperature)) - 1.0)
    flux = np.pi * planck * radius**2 / distance**2
    return flux.to(u.Jy, equivalencies=u.spectral_density(wavelength))


def read_data(sed_fits_file, stokes='I'):
    fits_header = fits.getheader(sed_fits_file)
    fits_data = fits.getdata(sed_fits_file)

    nr_wave = np.shape(fits_data)[2]
    sed_wavelengths = np.zeros(nr_wave) * u.m
    for i_wave in range(nr_wave):
        sed_wavelengths[i_wave] = float(fits_header[f'HIERARCH WAVELENGTH{i_wave+1}']) * u.m

    _stokes = ['I', 'Q', 'U', 'V']
    sed_data = {}
    for i_s, i_stokes in enumerate(_stokes):
        sed_data[i_stokes] = fits_data[i_s,:,:] * u.Jy

    return sed_wavelengths, sed_data[stokes]


def compare():
    sed_wavelengths, sed_data = read_data(
        os.path.join('projects', 'test', 'stellar_sed', 'dust_mc', 'data', 'polaris_detector_nr0001_sed.fits.gz'))
    reference = flux_planck_law(sed_wavelengths)

    for sed, ref in zip(sed_data[0], reference):
        if abs(sed / ref - 1.0) > 1e-4:
            return False

    return True


if __name__ == '__main__':
    res = compare()
    if res:
        print('Test passed')
    else:
        print('Test failed')
