from astropy.io import fits
from astropy import units as u
import numpy as np


"""
----------------------------------------------
Stellar scattering at a sphere
----------------------------------------------

Test whether scattering in a homogenous sphere
looks the same from different viewing angles.
We test four wavelengths, and the sphere is
optically thick at the two shortest wavelengths
and optically thin at the two longer ones.

In order to pass this test, all detectors
(with different orientations) have to yield
the same stellar + scattered flux for each
wavelength. This is tested by calculating
the relative difference of every detector
to the results of the first detector
(Should get better with more photons)
"""


def read_data(sed_fits_file):
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

    return sed_data


def compare():
    sed_data = read_data('projects/test/stellar_scattering_sphere/dust_mc/data/polaris_detector_nr0001_sed.fits.gz')
    for det in range(2, 6):
        tmp_data = read_data(f'projects/test/stellar_scattering_sphere/dust_mc/data/polaris_detector_nr{det:04d}_sed.fits.gz')

        max_rel_diff = np.max(np.abs( sed_data['I'] / tmp_data['I'] - 1.0 ))
        if max_rel_diff > 1e-2:
            raise Exception(f'Test failed: Stokes I does not match (max. relative difference = {max_rel_diff})')
        
        sed_polarization = np.sqrt(sed_data['Q']**2 + sed_data['U']**2 + sed_data['V']**2) / sed_data['I']
        tmp_polarization = np.sqrt(tmp_data['Q']**2 + tmp_data['U']**2 + tmp_data['V']**2) / tmp_data['I']
        max_abs_diff = np.max(np.abs( sed_polarization - tmp_polarization ))
        if max_abs_diff > 1e-2:
            raise Exception(f'Test failed: Polarization does not match (max. absolute difference = {max_abs_diff})')
    
    max_polarization = np.max(sed_polarization)
    if max_polarization > 1e-2:
        raise Exception(f'Test failed: Polarization is too large (max. value = {max_polarization})')

    return True


if __name__ == '__main__':
    res = compare()
    if res:
        print('Test passed')
