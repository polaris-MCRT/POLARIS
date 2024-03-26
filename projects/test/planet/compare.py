from astropy.io import fits
from astropy import units as u
from astropy import constants as c
from scipy.special import erfc
import numpy as np
import glob
import os


def flux_planck_law(wavelength, temperature=6000*u.K, radius=1*u.R_sun, distance=3.0856775814671917e+17*u.m):
    planck = (2.0 * c.h * c.c**2 / wavelength**5) / (np.exp(c.h * c.c / (wavelength * c.k_B * temperature)) - 1.0)
    flux = np.pi * planck * radius**2 / distance**2
    return flux.to(u.Jy, equivalencies=u.spectral_density(wavelength))


def lambertian(phases, albedo=1.0):
    return albedo * 2.0 / (3.0 * np.pi) * (np.sin(phases) + (np.pi - phases) * np.cos(phases))


def lommelseeliger(phases, albedo=1.0):
    res = 0.125 * albedo
    return res * (1.0 + np.sin(0.5 * phases) * np.tan(0.5 * phases) * np.log(np.tan(0.25 * phases), out=np.zeros_like(phases), where=phases > 0.0))


def ocean(phases, foam_albedo=0.2, wind_speed=3.0, refractive_index=1.33):
    def getFresnelR(cosi, ri):
        sino_sq = (1.0 - cosi**2) / ri**2
        coso = np.sqrt(1.0 - sino_sq)
        rpll = (ri * cosi - coso) / (ri * cosi + coso)
        rper = (cosi - ri * coso) / (cosi + ri * coso)
        return 0.5 * (np.abs(rpll)**2 + np.abs(rper)**2)

    def getFresnelP(cosi, ri):
        sino_sq = (1.0 - cosi**2) / ri**2
        coso = np.sqrt(1.0 - sino_sq)
        rpll = (ri * cosi - coso) / (ri * cosi + coso)
        rper = (cosi - ri * coso) / (cosi + ri * coso)
        return 0.5 * (np.abs(rpll)**2 - np.abs(rper)**2)

    def getProb(cosn, sigma_sq):
        return np.exp(-(1.0 - cosn**2) / (sigma_sq * cosn**2)) / (np.pi * sigma_sq * cosn**3)

    def getShadow(cosi, coso, sigma_sq):
        sini = np.sqrt(1.0 - cosi**2)
        etai = cosi / (np.sqrt(sigma_sq) * sini)
        sino = np.sqrt(1.0 - coso**2)
        etao = coso / (np.sqrt(sigma_sq) * sino)
        lambdai = 0.5 * (np.exp(-etai**2) / (np.sqrt(np.pi) * etai) - erfc(etai))
        lambdao = 0.5 * (np.exp(-etao**2) / (np.sqrt(np.pi) * etao) - erfc(etao))
        return 1.0 / (1.0 + lambdai + lambdao)
    
    def getCoxMunkInt(sigma_sq, size):
        res_int = np.zeros(size)
        res_mu0 = np.zeros(size)
        step = 90.0 / (size - 1)
        points, weights = np.polynomial.legendre.leggauss(40)

        for i in range(size):
            cos_0 = np.cos(np.deg2rad(i * step))
            sin_0 = np.sqrt(1.0 - cos_0**2)
            res_mu0[i] = cos_0
            a1 = 0.0
            b1 = 0.5 * np.pi
            a2 = 0.0
            b2 = 2.0 * np.pi

            for (x1, w1) in zip(points, weights):
                xt = 0.5 * x1 * (b1 - a1) + 0.5 * (a1 + b1)
                prob = getProb(np.cos(xt), sigma_sq)
                for (x2, w2) in zip(points, weights):
                    xp = 0.5 * x2 * (b2 - a2) + 0.5 * (a2 + b2)
                    cos_i_n = sin_0 * np.tan(xt) * np.cos(xp) + cos_0
                    if cos_i_n > 0.0:
                        res_int[i] += w1 * w2 * cos_i_n * prob * np.sin(xt)
            res_int[i] *= 0.25 * (b1 - a1) * (b2 - a2)

        return res_int, res_mu0

    sigma_sq = 0.003 + 0.00512 * wind_speed
    frac = 2.95e-6 * wind_speed**3.52
    cm_norm, cm_mu = getCoxMunkInt(sigma_sq, 901)

    flux = np.zeros_like(phases)
    polarization = np.zeros_like(phases)
    points, weights = np.polynomial.legendre.leggauss(40)

    for a, alpha in enumerate(phases):
        a1 = 0.0
        b1 = np.pi
        a2 = alpha - 0.5 * np.pi
        b2 = 0.5 * np.pi

        for (x1, w1) in zip(points, weights):
            xt = 0.5 * x1 * (b1 - a1) + 0.5 * (a1 + b1)
            for (x2, w2) in zip(points, weights):
                xp = 0.5 * x2 * (b2 - a2) + 0.5 * (a2 + b2)
                mu_0 = np.sin(xt) * np.cos(xp)
                mu_r = np.sin(xt) * np.cos(alpha - xp)
                mu_i = np.cos(0.5 * alpha)
                mu_n = 0.5 * (mu_0 + mu_r) / mu_i

                prob = getProb(mu_n, sigma_sq)
                refl = getFresnelR(mu_i, refractive_index)
                pol = getFresnelP(mu_i, refractive_index)

                brdf_f = foam_albedo / np.pi
                brdf_o = 0.25 * prob * refl / (mu_0 * mu_r * mu_n)
                bpdf_o = 0.25 * prob * pol / (mu_0 * mu_r * mu_n)

                const = mu_i / (np.interp(mu_0, cm_mu[::-1], cm_norm[::-1]) * mu_n)
                brdf_o *= const
                bpdf_o *= const

                const = getShadow(mu_0, mu_r, sigma_sq)
                brdf_o *= const
                bpdf_o *= const

                brdf = frac * brdf_f + (1.0 - frac) * brdf_o
                bpdf = (1.0 - frac) * bpdf_o
                flux[a] += w1 * w2 * brdf * np.sin(xt) * mu_0 * mu_r
                polarization[a] += w1 * w2 * bpdf * np.sin(xt) * mu_0 * mu_r

        flux[a] *= 0.25 * (b1 - a1) * (b2 - a2)
        polarization[a] *= 0.25 * (b1 - a1) * (b2 - a2)

    return flux, np.abs(polarization / flux)


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
        "lambertian_exten", "lommelseeliger_exten", "ocean_exten",
        "rayleigh_tau0.3_surface1.0_exten", "rayleigh_tau5.0_surface0.0_exten", "rayleigh_tau10.0_surface0.3_exten",
        "lambertian_plane", "lommelseeliger_plane", "ocean_plane",
        "rayleigh_tau0.3_surface1.0_plane", "rayleigh_tau5.0_surface0.0_plane", "rayleigh_tau10.0_surface0.3_plane"
        ]

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
        elif benchmark.replace('_exten', '').replace('_plane', '') == 'ocean':
            # since polaris includes multiple scattering on the surface, only small wind speeds should match
            I_ref, P_ref = ocean(np.deg2rad(phase_angles))
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

        print(f'{benchmark} Test ... OK')


if __name__ == '__main__':
    compare()
