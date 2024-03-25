from astropy.io import fits
from astropy import units as u
from astropy import constants as c
import numpy as np
import glob
import os


def flux_planck_law(wavelength, temperature=6000*u.K, radius=1*u.R_sun, distance=3.0856775814671917e+17*u.m):
    planck = (2.0 * c.h * c.c**2 / wavelength**5) / (np.exp(c.h * c.c / (wavelength * c.k_B * temperature)) - 1.0)
    flux = np.pi * planck * radius**2 / distance**2
    return flux.to(u.Jy, equivalencies=u.spectral_density(wavelength))


def lambertian(phase, albedo=1.0):
    return albedo * 2.0 / (3.0 * np.pi) * (np.sin(phase) + (np.pi - phase) * np.cos(phase))


def lommelseeliger(phase, albedo=1.0):
    res = 0.125 * albedo
    return res * (1.0 + np.sin(0.5 * phase) * np.tan(0.5 * phase) * np.log(np.tan(0.25 * phase), out=np.zeros_like(phase), where=phase > 0.0))


def read_data(sed_fits_file):
    fits_header = fits.getheader(sed_fits_file)
    fits_data = fits.getdata(sed_fits_file)

    sed_data = {}

    phase_angle = float(fits_header['RANGLE2'])
    wavelength = float(fits_header['HIERARCH WAVELENGTH1']) * u.m

    _stokes = ['I', 'Q', 'U', 'V']
    for i_s, i_stokes in enumerate(_stokes):
        sed_data[i_stokes] = fits_data[i_s,:,:] * u.Jy
        planet_radius = 7.0e7 * u.m
        planet_distance = 0.1 * u.au
        sed_data[i_stokes] /= flux_planck_law(wavelength) * planet_radius**2 / planet_distance**2
        sed_data[i_stokes] = sed_data[i_stokes].to(u.dimensionless_unscaled)

    return phase_angle, sed_data


def compare():
    benchmarks = [
        "lambertian_exten", "lommelseeliger_exten",
        "rayleigh_tau0.3_surface1.0_exten", "rayleigh_tau5.0_surface0.0_exten", "rayleigh_tau10.0_surface0.3_exten",
        "lambertian_plane", "lommelseeliger_plane",
        "rayleigh_tau0.3_surface1.0_plane", "rayleigh_tau5.0_surface0.0_plane", "rayleigh_tau10.0_surface0.3_plane"]

    for benchmark in benchmarks:
        phase_angles = []
        sed_data = {'I': [], 'P': []}

        files = sorted(glob.glob(
            os.path.join('projects', 'test', 'planet', benchmark, 'data', 'polaris_detector_nr????_sed.fits.gz')))

        for file in files:
            _phase_angle, _sed_data = read_data(file)
            phase_angles.append(_phase_angle)
            sed_data['I'].append(_sed_data['I'][0][0])
            sed_data['P'].append(np.sqrt(_sed_data['Q'][0][0]**2 + _sed_data['U'][0][0]**2) / _sed_data['I'][0][0])

        if benchmark.replace('_exten', '').replace('_plane', '') == 'lambertian':
            I_ref = lambertian(np.deg2rad(phase_angles))
            P_ref = np.zeros_like(I_ref)
        elif benchmark.replace('_exten', '').replace('_plane', '') == 'lommelseeliger':
            I_ref = lommelseeliger(np.deg2rad(phase_angles))
            P_ref = np.zeros_like(I_ref)
        else:
            ref_file = os.path.join('projects', 'test', 'planet', 'reference', benchmark.replace('_exten', '').replace('_plane', '') + '.txt')
            ref_data = np.loadtxt(ref_file, comments="#", delimiter=",")
            I_ref = ref_data[:,1] * u.dimensionless_unscaled
            Q_ref = ref_data[:,2] * u.dimensionless_unscaled
            P_ref = [np.abs(Q / I) if I > 0.0 else 0.0 for I, Q in zip(I_ref, Q_ref)]

        for phase, sed, ref in zip(phase_angles, sed_data['I'], I_ref):
            if abs(sed - ref) > 1e-2:
                raise Exception(f'{benchmark} Test FAILED (intensity, {sed} != {ref} at {phase} deg)')
        
        for phase, pol, ref in zip(phase_angles, sed_data['P'], P_ref):
            if abs(pol - ref) > 1e-2:
                raise Exception(f'{benchmark} Test FAILED (polarization, {pol} != {ref} at {phase} deg)')

    return True


if __name__ == '__main__':
    res = compare()
    if res:
        print('Tests OK')
