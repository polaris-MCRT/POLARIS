from re import I
from astropy.io import fits
from astropy import units as u
from astropy import constants as c
import numpy as np
import matplotlib.pyplot as plt
import glob
import os, sys


def flux_planck_law(wavelength, temperature=6000*u.K, radius=1*u.R_sun, distance=10*u.pc):
    planck = (2.0 * c.h * c.c**2 / wavelength**5) / (np.exp(c.h * c.c / (wavelength * c.k_B * temperature)) - 1.0)
    flux = np.pi * planck * radius**2 / distance**2
    return flux.to(u.Jy, equivalencies=u.spectral_density(wavelength))

def read_data(fits_file, model):
    fits_header = fits.getheader(fits_file)
    fits_data = fits.getdata(fits_file)

    sed_data = {}
    _stokes = ['I', 'Q', 'U', 'V']

    if model in ['cloudy', 'rayleigh', 'ocean']:
        phase_angle = float(fits_header['RANGLE2'])
        wavelength = float(fits_header['HIERARCH WAVELENGTH1']) * u.m
        for i_s, i_stokes in enumerate(_stokes):
            sed_data[i_stokes] = fits_data[i_s,0,0] * u.Jy
            planet_radius = 7.0e7 * u.m
            planet_distance = 0.1 * u.au
            sed_data[i_stokes] /= flux_planck_law(wavelength) * planet_radius**2 / planet_distance**2
            sed_data[i_stokes] = sed_data[i_stokes].to(u.dimensionless_unscaled)

        return phase_angle, sed_data

    elif model == 'ocean_res':
        wavelength = float(fits_header['HIERARCH WAVELENGTH1']) * u.m
        for i_s, i_stokes in enumerate(_stokes):
            sed_data[i_stokes] = fits_data[i_s,0,:,:] * u.Jy
            sed_data[i_stokes] = sed_data[i_stokes].to(u.W / u.m**2 / u.um, equivalencies=u.spectral_density(wavelength))
        
        return sed_data

    elif model == 'methane':
        nr_wave = int(fits_header['NAXIS1'])
        wavelengths = np.zeros(nr_wave) * u.m
        for iw in range(nr_wave):
            wavelengths[iw] = float(fits_header[f'HIERARCH WAVELENGTH{iw+1}']) * u.m
        for i_s, i_stokes in enumerate(_stokes):
            sed_data[i_stokes] = fits_data[i_s,0,:] * u.Jy
            planet_radius = 7.0e7 * u.m
            planet_distance = 0.1 * u.au
            sed_data[i_stokes] /= flux_planck_law(wavelengths) * planet_radius**2 / planet_distance**2
            sed_data[i_stokes] = sed_data[i_stokes].to(u.dimensionless_unscaled)

        return wavelengths.to(u.um), sed_data

    elif model == 'ringed':
        wavelength = float(fits_header['HIERARCH WAVELENGTH1']) * u.m
        for i_s, i_stokes in enumerate(_stokes):
            sed_data[i_stokes] = fits_data[i_s,0,:,:] * u.Jy
            sed_data[i_stokes] = sed_data[i_stokes].to(u.W / u.m**2 / u.um, equivalencies=u.spectral_density(wavelength))
        
        return sed_data
    
    elif model == 'venus':
        phase_angle = float(fits_header['RANGLE2'])
        wavelength = float(fits_header['HIERARCH WAVELENGTH1']) * u.m
        for i_s, i_stokes in enumerate(_stokes):
            sed_data[i_stokes] = fits_data[i_s,0,0] * u.Jy
            planet_radius = 6.05e6 * u.m
            planet_distance = 0.7282 * u.au
            sed_data[i_stokes] /= flux_planck_law(wavelength, temperature=5772*u.K, distance=1*u.pc) * planet_radius**2 / planet_distance**2
            sed_data[i_stokes] = sed_data[i_stokes].to(u.dimensionless_unscaled)

        return phase_angle, sed_data
    

def load_obs_data():
    with open(os.path.join('projects', 'venus', 'coffeen_gehrels_1969.txt')) as file:
        phase_angles = {'340nm': [], '365nm': [], '445nm': [], '520nm': [], '550nm': [], '655nm': [], '685nm': [], '740nm': [], '875nm': [], '990nm': []}
        polarization = {'340nm': [], '365nm': [], '445nm': [], '520nm': [], '550nm': [], '655nm': [], '685nm': [], '740nm': [], '875nm': [], '990nm': []}
        wavel_dict = {'N': '340nm', 'U': '365nm', 'B': '445nm', 'GU': '520nm', 'GI': '550nm', 'R4': '655nm', 'R': '685nm', 'R1': '740nm', 'R2': '875nm', 'I': '990nm'}
        for line in file.readlines():
            if line.startswith('#'):
                continue
            alpha, w_filt, pol, error, angle = line.split()
            wavel = wavel_dict[w_filt]
            phase_angles[wavel].append(float(alpha))
            sign = 1
            try:
                angle = float(angle)
                if angle > 180: angle -= 180    
                if 45 <= angle < 135: sign = -1
            except:
                if angle == '-': sign = -1
            polarization[wavel].append(float(pol) * sign)
        
    return phase_angles, polarization

def plot(model):
    if model in ['cloudy', 'rayleigh', 'ocean']:
        files = sorted(glob.glob(
            os.path.join('projects', model, 'data', 'polaris_detector_nr????_sed.fits.gz')))

        sed_data = {'I': [], 'P': []}
        phase_angles = []
        for file in files:
            _phase_angle, _sed_data = read_data(file, model)
            phase_angles.append(_phase_angle)
            sed_data['I'].append(_sed_data['I'])
            sed_data['P'].append(np.sqrt(_sed_data['Q']**2 + _sed_data['U']**2) / _sed_data['I'])

        fig, ax = plt.subplots(layout='constrained')
        ax.plot(phase_angles, sed_data['I'], label='Normalized flux')
        ax.plot(phase_angles, sed_data['P'], label='Degree of polarization')
        ax.set_xlabel('Phase angle [deg]')
        ax.legend()
        fig.savefig(os.path.join('projects', model, f'{model}.png'), bbox_inches='tight')

        if model == 'ocean':
            nr = len(files) // 2
            file = os.path.join('projects', model, 'data', f'polaris_detector_nr{nr:04d}.fits.gz')
            sed_data = read_data(file, f'{model}_res')
            sed_data['P'] = np.sqrt(sed_data['Q']**2 + sed_data['U']**2) / sed_data['I']

            fig1, ax1 = plt.subplots(layout='constrained')
            im1 = ax1.imshow(sed_data['I'].value, origin='lower')
            fig1.colorbar(im1, label='Surface brightness [W m$^{-2}$ µm$^{-1}$ px$^{-1}$]')
            fig1.savefig(os.path.join('projects', model, f'{model}_I.png'), bbox_inches='tight')

            fig2, ax2 = plt.subplots(layout='constrained')
            im2 = ax2.imshow(sed_data['P'].value, origin='lower')
            fig2.colorbar(im2, label='Degree of Polarization')
            fig2.savefig(os.path.join('projects', model, f'{model}_P.png'), bbox_inches='tight')

    elif model == 'methane':
        file = os.path.join('projects', model, 'data', 'polaris_detector_nr0001_sed.fits.gz')

        wavelengths, sed_data = read_data(file, model)
        sed_data['P'] = np.sqrt(sed_data['Q']**2 + sed_data['U']**2) / sed_data['I']

        fig, ax = plt.subplots(layout='constrained')
        ax.plot(wavelengths, sed_data['I'], label='Normalized flux')
        ax.plot(wavelengths, sed_data['P'], label='Degree of polarization')
        ax.set_xlabel('Wavelength [µm]')
        ax.legend()
        fig.savefig(os.path.join('projects', model, f'{model}.png'), bbox_inches='tight')

    elif model == 'ringed':
        file = os.path.join('projects', model, 'data', 'polaris_detector_nr0001.fits.gz')
        sed_data = read_data(file, model)
        sed_data['P'] = np.sqrt(sed_data['Q']**2 + sed_data['U']**2) / sed_data['I']

        fig1, ax1 = plt.subplots(layout='constrained')
        im1 = ax1.imshow(sed_data['I'].value, origin='lower')
        fig1.colorbar(im1, label='Surface brightness [W m$^{-2}$ µm$^{-1}$ px$^{-1}$]')
        fig1.savefig(os.path.join('projects', model, f'{model}_I.png'), bbox_inches='tight')

        fig2, ax2 = plt.subplots(layout='constrained')
        im2 = ax2.imshow(sed_data['P'].value, origin='lower')
        fig2.colorbar(im2, label='Degree of Polarization')
        fig2.savefig(os.path.join('projects', model, f'{model}_P.png'), bbox_inches='tight')
    
    elif model == 'venus':
        phase_obs, pol_obs = load_obs_data()
        for wavel in ['340nm', '365nm', '445nm', '520nm', '550nm', '655nm', '685nm', '740nm', '875nm', '990nm']:
            files = sorted(glob.glob(
                os.path.join('projects', model, wavel, 'data', 'polaris_detector_nr????_sed.fits.gz')))

            sed_data = {'P': []}
            phase_angles = []
            for file in files:
                _phase_angle, _sed_data = read_data(file, model)
                phase_angles.append(_phase_angle)
                pol = np.sqrt(_sed_data['Q']**2 + _sed_data['U']**2) / _sed_data['I']
                sign = -np.sign(_sed_data['Q'])
                sed_data['P'].append(100 * pol * sign)

            fig, ax = plt.subplots(layout='constrained')
            ax.plot(phase_angles, sed_data['P'], label='POLARIS')
            ax.scatter(phase_obs[wavel], pol_obs[wavel], label='Coffeen & Gehrels (1969)')
            ax.set_title(wavel)
            ax.legend()
            ax.set_ylabel('Degree of Polarization [%]')
            ax.set_xlabel('Phase angle [deg]')
            fig.savefig(os.path.join('projects', model, f'{model}_{wavel}.png'), bbox_inches='tight')

    else:
        print('Model name unknown.')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Wrong amount of arguments. Needs exactly 1 argument: python plot.py model_name (cloudy, methane, ocean, rayleigh, ringed, venus)')
    else:
        plot(sys.argv[1])
