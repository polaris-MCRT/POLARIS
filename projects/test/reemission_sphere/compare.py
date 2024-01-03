from astropy.io import fits
from astropy import units as u
from astropy import constants as c
import numpy as np


"""
----------------------------------------------
Re-emission of a homogeneous sphere
----------------------------------------------

We test whether the re-emission radiation is
a modified Blackbody described by Planck law
with a wavelength-dependent emission coefficient.

We test whether both detector geometries of the
raytracer (cartesian and polar) yield the same results.

In order to pass this test, the simulated radiation
from POLARIS should match the analytical solution.
"""


def planck_law(wavelength, temperature):
    return (2.0 * c.h * c.c**2 / wavelength**5) / (np.exp(c.h * c.c / (wavelength * c.k_B * temperature)) - 1.0)


def flux_planck_law(wavelength, temperature=4500*u.K, radius=2*u.R_sun, distance=4.32e+18*u.m):
    planck = planck_law(wavelength, temperature)
    flux = np.pi * planck * radius**2 / distance**2
    return flux.to(u.Jy, equivalencies=u.spectral_density(wavelength))


def calc_flux_ana():
    # Dust grain radius, default 1µm; compare with dust_mixture_001.dat
    s = 1 *u.micron

    # Bulk density, default 3800kg/m^3; compare with dust_mixture_001.dat
    rho_bulk = 3800 * u.kg / (u.m)**3

    # Wavelength for SED
    wavelength = np.loadtxt('projects/test/reemission_sphere/dust/data/dust_mixture_001.dat', skiprows=27, usecols=0) * u.m

    # Extinction cross-section
    C_ext = np.loadtxt('projects/test/reemission_sphere/dust/data/dust_mixture_001.dat', skiprows=27, usecols=1) * u.m**2

    # Absorption (= emission) cross-section
    C_abs = np.loadtxt('projects/test/reemission_sphere/dust/data/dust_mixture_001.dat', skiprows=27, usecols=3) * u.m**2

    # Dust temperature, default 300K
    T_dust = 300 *u.K

    # Inner radius of slab, default 0.5au
    R_inner = 0.5 * u.au

    # Outer radius of slab, default 10au
    R_outer = 10 * u.au

    # Dust mass of sphere, default 1e-7 M_sun
    M_sphere = 1e-7 * c.M_sun

    # Distance to observer, default 4.32e18m (~140pc)
    distance = 4.32e18 * u.m

    volume_sphere = 4./3. * np.pi * (R_outer**3 - R_inner**3)
    mass_density = M_sphere / volume_sphere
    m_grain = 4./3. * np.pi * s**3 * rho_bulk
    number_density = mass_density / m_grain

    # Dust emisison per m**3
    # Factor pi to account for each radiation direction
    # Factor 4 is the ratio of surface to cross-section of a sphere
    emission = np.pi * 4 * C_abs * planck_law(wavelength, T_dust) * number_density

    # Number of steps through the sphere (rho direction)
    n_step_rho = int(1e5)

    # Path length through a half-sphere with inner hole for
    # different rho (cylindrical distance from center)
    # h = sin(x) * R_outer, sin(x) is evenly distributed
    path_length_half_sphere = R_outer * np.linspace(1,0,n_step_rho)
    # Inner part starts at boundary of inner hole, outer part at 0 (center line)
    rho = R_outer * np.sqrt( 1 - np.linspace(1,0,n_step_rho)**2 )
    path_length_half_sphere -= np.sqrt( np.maximum(0,R_inner**2 - rho**2) )

    # Number of steps through the sphere (z direction)
    n_step_z = int(1e5)

    # calculate the optical depth through sphere (with inner hole)
    # for different offsets from the center
    # Factor 2: path length is only for a half-sphere
    # tau_total_path: 1st dim = wavelength, 2nd dim = rho
    tau_total_path = (C_ext[:,None] * number_density * 2*path_length_half_sphere[None,:]).decompose()

    # The radiation of each step is reduced by np.exp( -tau_total/n_step )
    # Thus, we need to sum over these reductions from i=1 to n of np.exp( -tau/n*i )
    # Wolfram Alpha says that this is (1 - exp(-tau)) / (exp(tau/n) - 1)
    # exp_tau_sum: 1st dim = wavelength, 2nd dim = rho
    exp_tau_sum = np.divide(
        1 - np.exp(-tau_total_path),
        np.exp(tau_total_path/n_step_z) - 1,
        where=tau_total_path/n_step_z > 1e-15)

    # Total flux of a column, i.e. per projected detector area, is emission per volume
    # times step width reduced by the optical depth
    # flux_tmp: 1st dim = wavelength, 2nd dim = rho
    flux_tmp = emission[:,None] * 2*path_length_half_sphere[None,:] / n_step_z * exp_tau_sum

    # integrate over phi
    flux_tmp *= 2*np.pi
    # integrate over rho*d(rho)
    flux_sphere_ana = np.trapz(flux_tmp * rho, x=rho, axis=1)
    # account for the distance of the object
    flux_sphere_ana /= 4*np.pi*distance**2

    # convert to Jy
    return flux_sphere_ana.to(u.Jy, equivalencies=u.spectral_density(wavelength))


def read_data(sed_fits_file, map_fits_file, stokes='I'):
    sed_fits_header = fits.getheader(sed_fits_file)
    sed_fits_data = fits.getdata(sed_fits_file)

    nr_wave = np.shape(sed_fits_data)[2]
    sed_wavelengths = np.zeros(nr_wave) * u.m
    for i_wave in range(nr_wave):
        sed_wavelengths[i_wave] = float(sed_fits_header[f'HIERARCH WAVELENGTH{i_wave+1}']) * u.m

    _stokes = ['I', 'Q', 'U', 'V', 'TAU']
    sed_data = {}
    for i_s, i_stokes in enumerate(_stokes):
        if i_stokes == 'TAU':
            sed_data[i_stokes] = sed_fits_data[i_s,0,:] * u.dimensionless_unscaled
        else:
            sed_data[i_stokes] = sed_fits_data[i_s,0,:] * u.Jy

    map_fits_header = fits.getheader(map_fits_file)
    map_fits_data = fits.getdata(map_fits_file)

    nr_wave = np.shape(sed_fits_data)[1]
    map_wavelengths = np.zeros(nr_wave) * u.m
    for i_wave in range(nr_wave):
        map_wavelengths[i_wave] = float(map_fits_header[f'HIERARCH WAVELENGTH{i_wave+1}']) * u.m

    _stokes = ['I', 'Q', 'U', 'V']
    map_data = {}
    for i_s, i_stokes in enumerate(_stokes):
        map_data[i_stokes] = map_fits_data[i_s,:,:,:] * u.Jy / u.pix

    return sed_data[stokes], map_data[stokes]


def compare():
    sed_data_pol, map_data_pol = read_data(
        'projects/test/reemission_sphere/dust/data/polaris_detector_nr0001_sed.fits.gz',
        'projects/test/reemission_sphere/dust/data/polaris_detector_nr0001.fits.gz')
    sed_data_car, map_data_car = read_data(
        'projects/test/reemission_sphere/dust/data/polaris_detector_nr0002_sed.fits.gz',
        'projects/test/reemission_sphere/dust/data/polaris_detector_nr0002.fits.gz')
    reference = calc_flux_ana()

    max_rel_diff = np.max(np.abs( sed_data_pol / reference - 1.0 ))
    if max_rel_diff > 1e-3:
        raise Exception(f'Test failed: Polar detector and reference do not match (max. relative difference = {max_rel_diff})')

    max_rel_diff = np.max(np.abs( sed_data_car / sed_data_pol - 1.0 ))
    if max_rel_diff > 1e-3:
        raise Exception(f'Test failed: Cartesian and polar detector do not match (max. relative difference = {max_rel_diff})')
    
    max_rel_diff = np.max(np.abs( np.sum(map_data_pol, axis=(1,2)) * u.pix / reference - 1.0 ))
    if max_rel_diff > 1e-3:
        raise Exception(f'Test failed: Sum of polar map detector and reference do not match (max. relative difference = {max_rel_diff})')

    max_rel_diff = np.max(np.abs( np.sum(map_data_car, axis=(1,2)) / np.sum(map_data_pol, axis=(1,2)) - 1.0 ))
    if max_rel_diff > 1e-3:
        raise Exception(f'Test failed: Sum of cartesian and polar map detector do not match (max. relative difference = {max_rel_diff})')

    return True


if __name__ == '__main__':
    res = compare()
    if res:
        print('Test passed')
