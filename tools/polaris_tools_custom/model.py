#!/usr/bin/env python
# -*- coding: utf-8 -*-

from polaris_tools_modules.base import Model
from polaris_tools_modules.atmosphere import AtmosphereRoutines
import numpy as np

"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_model_dict(dictionary):
    model_dict = {
        'jupiter_like': JupiterLike,
        'earth_like': EarthLike,
        'custom': CustomModel,
        'custom_disk': CustomDisk,
    }
    dictionary.update(model_dict)


class JupiterLike(Model):
    def __init__(self):
        Model.__init__(self)
        self.ar = AtmosphereRoutines()

        self.planetary_mass = self.math.const['M_jup']
        self.planetary_radius = self.math.const['R_jup']
        self.gravity = self.math.const['G'] * self.planetary_mass / self.planetary_radius**2

        # ==================================================
        ### Jupiter atmospheric profile
        ### Seiff et al. 1998, J. Geophysical Research, 103, 22857

        # pressure profile [mbar]
        self.pressure_profile = np.array([
            2.200e+04, 2.100e+04, 2.000e+04, 1.900e+04, 1.800e+04, 1.700e+04, 1.600e+04, 1.500e+04, 1.400e+04, 1.300e+04,
            1.200e+04, 1.100e+04, 1.000e+04, 9.000e+03, 8.000e+03, 7.000e+03, 6.000e+03, 5.000e+03, 4.800e+03, 4.600e+03,
            4.400e+03, 4.200e+03, 4.000e+03, 3.800e+03, 3.600e+03, 3.400e+03, 3.200e+03, 3.000e+03, 2.800e+03, 2.600e+03,
            2.400e+03, 2.200e+03, 2.000e+03, 1.800e+03, 1.600e+03, 1.400e+03, 1.200e+03, 1.000e+03, 9.000e+02, 8.000e+02,
            7.000e+02, 6.000e+02, 5.000e+02, 4.680e+02, 3.515e+02, 1.358e+02, 4.374e+01, 1.640e+01, 7.177e+00, 3.079e+00,
            1.342e+00, 6.192e-01, 2.824e-01, 1.257e-01, 5.593e-02, 2.475e-02, 1.094e-02, 4.800e-03, 2.177e-03, 1.152e-03,
            6.036e-04, 3.367e-04, 2.107e-04, 1.463e-04, 1.072e-04, 8.149e-05, 6.518e-05, 5.265e-05, 4.256e-05, 3.484e-05,
            2.903e-05, 2.445e-05, 2.067e-05, 1.749e-05, 1.481e-05, 1.256e-05, 1.068e-05, 9.117e-06, 7.827e-06, 6.765e-06,
            5.887e-06, 5.152e-06, 4.527e-06, 3.990e-06, 3.522e-06, 3.112e-06, 2.749e-06, 2.427e-06, 2.141e-06, 1.888e-06,
            1.664e-06, 1.468e-06, 1.297e-06, 1.148e-06, 1.018e-06, 9.633e-07])
        # mbar -> Pa
        self.pressure_profile *= 1e2

        # temperature_profile [K]
        self.temperature_profile = np.array([
            427.71, 421.83, 415.75, 409.46, 402.93, 396.13, 389.04, 381.64, 373.87, 365.71,
            357.08, 347.92, 338.16, 327.67, 316.31, 303.88, 290.10, 274.54, 271.16, 267.68,
            264.09, 260.38, 256.53, 252.55, 248.40, 244.09, 239.59, 234.87, 229.92, 224.71,
            219.19, 213.32, 207.06, 200.31, 193.00, 184.99, 176.11, 166.10, 160.56, 154.57,
            148.04, 140.84, 132.79, 130.00, 122.9,  113.2,  122.6,  143.8,  158.1,  149.8,
            160.5,  168.6,  157.4,  158.2,  157.2,  155.7,  151.0,  152.8,  177.6,  194.2,
            198.6,  231.7,  289.2,  370.2,  392.9,  483.1,  535.7,  535.8,  545.9,  594.3,
            642.4,  664.0,  670.3,  669.8,  674.6,  680.2,  692.7,  712.8,  741.0,  776.1,
            812.0,  840.0,  862.2,  877.2,  884.4,  885.2,  881.4,  875.2,  869.1,  865.4,
            866.0,  872.2,  883.9,  897.0,  903.4,  899.9])

        # boundaries of the grid [m]
        # first value is the inner model radius
        # last value is the outer model radius
        self.spherical_parameter['radius_list'] = np.array([
            -1.324e+05, -1.294e+05, -1.263e+05, -1.230e+05, -1.197e+05, -1.162e+05, -1.125e+05, -1.087e+05, -1.047e+05, -1.005e+05,
            -9.608e+04, -9.137e+04, -8.635e+04, -8.095e+04, -7.512e+04, -6.874e+04, -6.169e+04, -5.376e+04, -5.204e+04, -5.028e+04,
            -4.845e+04, -4.657e+04, -4.463e+04, -4.261e+04, -4.052e+04, -3.835e+04, -3.609e+04, -3.373e+04, -3.125e+04, -2.865e+04,
            -2.591e+04, -2.300e+04, -1.991e+04, -1.660e+04, -1.302e+04, -9.120e+03, -4.820e+03,  0.000e+00,  2.660e+03,  5.530e+03,
             8.650e+03,  1.209e+04,  1.595e+04,  1.729e+04,  2.330e+04,  4.000e+04,  6.000e+04,  8.000e+04,  1.000e+05,  1.200e+05,
             1.400e+05,  1.600e+05,  1.800e+05,  2.000e+05,  2.200e+05,  2.400e+05,  2.600e+05,  2.800e+05,  3.000e+05,  3.200e+05,
             3.400e+05,  3.600e+05,  3.800e+05,  4.000e+05,  4.200e+05,  4.400e+05,  4.600e+05,  4.800e+05,  5.000e+05,  5.200e+05,
             5.400e+05,  5.600e+05,  5.800e+05,  6.000e+05,  6.200e+05,  6.400e+05,  6.600e+05,  6.800e+05,  7.000e+05,  7.200e+05,
             7.400e+05,  7.600e+05,  7.800e+05,  8.000e+05,  8.200e+05,  8.400e+05,  8.600e+05,  8.800e+05,  9.000e+05,  9.200e+05,
             9.400e+05,  9.600e+05,  9.800e+05,  1.000e+06,  1.020e+06,  1.029e+06])
        self.spherical_parameter['radius_list'] += self.planetary_radius

        # ==================================================

        # Use spherical coordinate system
        self.parameter['grid_type'] = 'spherical'
        # sf_r = 0: user defined radial boundaries
        self.spherical_parameter['sf_r'] = 0

        # set inner and outer grid radius (bottom and top of atmosphere)
        self.parameter['inner_radius'] = self.spherical_parameter['radius_list'][0]
        self.parameter['outer_radius'] = self.spherical_parameter['radius_list'][-1]

        # number of radial grid cells
        self.spherical_parameter['n_r'] = len(self.spherical_parameter['radius_list']) - 1

        # sf_th = -1: equally distributed theta cells
        self.spherical_parameter['sf_th'] = -1
        # number of theta cells
        self.spherical_parameter['n_th'] = 1

        # sf_ph = -1: equally distributed phi cells
        self.spherical_parameter['sf_ph'] = -1
        # number of phi cells
        self.spherical_parameter['n_ph'] = 1

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters and update
        disk parameter that depend on other parameter."""
        # Use extra_parameter to adjust the model without changing the model.py file

    def dust_density_distribution(self):
        density = np.zeros(2, dtype=float)
        # density[0]: molecular hydrogen (input/cross_sections/molecular_hydrogen.dat)
        # density[1]: ammonia ice haze/cloud (input/refractive_indices/methane_ice_mo94.nk)

        pos = np.linalg.norm(self.position)
        idx = np.searchsorted(self.spherical_parameter['radius_list'], pos)

        # calculate the number density assuming hydrostatic equilibrium and an ideal gas
        density[0] = self.ar.getNumberDensityFromPressureAltitude(
            self.pressure_profile[idx-1] - self.pressure_profile[idx],
            self.spherical_parameter['radius_list'][idx] - self.spherical_parameter['radius_list'][idx-1],
            self.ar.getMolarMass('1H2'),
            self.gravity)

        # cloud particles from 1 bar to 0.1358 bar
        # np.interp assumes increasing x data points -> reverse arrays
        idx_bot = np.searchsorted(self.pressure_profile[::-1], 1.0 * 1e5)
        idx_top = np.searchsorted(self.pressure_profile[::-1], 0.1358 * 1e5)
        if self.spherical_parameter['radius_list'][idx_bot] < pos < self.spherical_parameter['radius_list'][idx_top]:
            tau = 5.0             # optical depth
            c_ext = 1.45262e-12   # precalculated cross section [m^2] using r_eff = 0.5 µm, veff = 0.05 (550 nm)
            density[1] = self.ar.getNumberDensityFromOpticalDepth(
                tau,
                self.spherical_parameter['radius_list'][idx_top] - self.spherical_parameter['radius_list'][idx_bot],
                c_ext)

        return density

    def dust_temperature(self):
        # return temperature at self.position
        pos = np.linalg.norm(self.position)
        # interpolate temperature profile based on current position
        return np.interp(pos, self.spherical_parameter['radius_list'], self.temperature_profile)


class EarthLike(Model):
    def __init__(self):
        Model.__init__(self)
        self.ar = AtmosphereRoutines()

        self.planetary_mass = self.math.const['M_earth']
        self.planetary_radius = self.math.const['R_earth']
        self.gravity = self.math.const['G'] * self.planetary_mass / self.planetary_radius**2

        # ==================================================
        ### Earth profile
        ### McClatchey (1972), Optical Properties of the Atmosphere
        ### Midlatitude Summer

        # pressure profile [mbar]
        self.pressure_profile = np.array([
            1.013e+03, 9.020e+02, 8.020e+02, 7.100e+02, 6.280e+02, 5.540e+02, 4.870e+02, 4.260e+02, 3.720e+02, 3.240e+02,
            2.810e+02, 2.430e+02, 2.090e+02, 1.790e+02, 1.530e+02, 1.300e+02, 1.110e+02, 9.500e+01, 8.120e+01, 6.950e+01,
            5.950e+01, 5.100e+01, 4.370e+01, 3.760e+01, 3.220e+01, 2.770e+01, 1.320e+01, 6.520e+00, 3.330e+00, 1.760e+00,
            9.510e-01, 6.710e-02, 3.000e-04])
        # mbar -> Pa
        self.pressure_profile *= 1e2

        # temperature profile [K]
        self.temperature_profile = np.array([
            294, 290, 285, 279, 273, 267, 261, 255, 248, 242,
            235, 229, 222, 216, 216, 216, 216, 216, 216, 217,
            218, 219, 220, 222, 223, 224, 234, 245, 258, 270,
            276, 218, 210])

        # boundaries [m]
        self.boundaries = np.array([
            0.00e+00, 1.00e+03, 2.00e+03, 3.00e+03, 4.00e+03, 5.00e+03, 6.00e+03, 7.00e+03, 8.00e+03, 9.00e+03,
            1.00e+04, 1.10e+04, 1.20e+04, 1.30e+04, 1.40e+04, 1.50e+04, 1.60e+04, 1.70e+04, 1.80e+04, 1.90e+04,
            2.00e+04, 2.10e+04, 2.20e+04, 2.30e+04, 2.40e+04, 2.50e+04, 3.00e+04, 3.50e+04, 4.00e+04, 4.50e+04,
            5.00e+04, 7.00e+04, 1.00e+05]) + self.planetary_radius
        # ==================================================

        # Use spherical coordinate system
        self.parameter['grid_type'] = 'spherical'
        # sf_r = 0: user defined radial boundaries
        self.spherical_parameter['sf_r'] = 0

        # set inner and outer grid radius (bottom and top of atmosphere)
        self.parameter['inner_radius'] = self.spherical_parameter['radius_list'][0]
        self.parameter['outer_radius'] = self.spherical_parameter['radius_list'][-1]

        # number of radial grid cells
        self.spherical_parameter['n_r'] = len(self.spherical_parameter['radius_list']) - 1

        # sf_th = -1: equally distributed theta cells
        self.spherical_parameter['sf_th'] = -1
        # number of theta cells
        self.spherical_parameter['n_th'] = 1

        # sf_ph = -1: equally distributed phi cells
        self.spherical_parameter['sf_ph'] = -1
        # number of phi cells
        self.spherical_parameter['n_ph'] = 1

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters and update
        disk parameter that depend on other parameter."""
        # Use extra_parameter to adjust the model without changing the model.py file

    def dust_density_distribution(self):
        density = np.zeros(2, dtype=float)
        # density[0]: air (input/dust/air.dat)
        # density[1]: clouds (input/dust/water_s81.nk)

        pos = np.linalg.norm(self.position)
        idx = np.searchsorted(self.spherical_parameter['radius_list'], pos)

        # calculate the number density assuming hydrostatic equilibrium and an ideal gas
        density[0] = self.ar.getNumberDensityFromPressureAltitude(
            self.pressure_profile[idx-1] - self.pressure_profile[idx],
            self.spherical_parameter['radius_list'][idx] - self.spherical_parameter['radius_list'][idx-1],
            self.ar.getMolarMass('air'),
            self.gravity)

        # water cloud particles from 802 mbar to 628 mbar (2 km and 4 km)
        idx_bot = np.searchsorted(self.spherical_parameter['radius_list'], self.planetary_radius + 2e3)
        idx_top = np.searchsorted(self.spherical_parameter['radius_list'], self.planetary_radius + 4e3)
        if self.spherical_parameter['radius_list'][idx_bot] < pos < self.spherical_parameter['radius_list'][idx_top]:
            tau = 2.0               # optical depth
            c_ext = 2.0579466e-11   # using r_eff = 2 µm, veff = 0.1, (550 nm)
            density[1] = self.ar.getNumberDensityFromOpticalDepth(
                tau,
                self.spherical_parameter['radius_list'][idx_top] - self.spherical_parameter['radius_list'][idx_bot],
                c_ext)

        return density
    
    def dust_temperature(self):
        # return temperature at self.position
        pos = np.linalg.norm(self.position)
        # interpolate temperature profile based on current position
        return np.interp(pos, self.spherical_parameter['radius_list'], self.temperature_profile)


class CustomModel(Model):
    """Change this to the model you want to use.
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        # Set parameters of the custom model (see parent Model class for all available options)
        self.parameter['grid_type'] = 'spherical'
        self.parameter['inner_radius'] = 0.1 * self.math.const['au']
        self.parameter['outer_radius'] = 100.0 * self.math.const['au']
        self.parameter['gas_mass'] = 1e-2 * self.math.const['M_sun']

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters and update
        model parameter that depend on other parameter.
        """
        # Use extra_parameter to adjust the model without changing the model.py file

    def gas_density_distribution(self):
        """Define here your routine to calculate the density at a given position
        in the model space.

        Notes:
            Use 'self.position' to calculate the quantity depending on position.
            self.position[0] = x coordinate
            self.position[1] = y coordinate
            self.position[2] = z coordinate

            Define also the following routines if necessary:
                dust_density_distribution(self), gas_temperature(self),
                dust_temperature(self), velocity_field(
                    self), magnetic_field(self),
                dust_id(self), dust_min_size(self), dust_max_size(
                    self), dust_size_param(self)

            xyz_density_distribution can return a density or 2D list of densities.
                - The first dimension is used to define multiple density distributions
                    for different dust compositions (see CustomDust in dust.py for explanation)
                - With the second dimension, multiple regions of the density distribution
                    of the same dust composition can be normalized individually to different total masses.
                - The self.parameter['gas_mass'] needs to have the same dimension and size as the return of this

        Returns:
            float: Gas density at a given position.
        """
        gas_density = 1.0
        # Or a function that depends on the position!
        # See other models to find prewritten distributions (shakura & sunyaev, ...)
        return gas_density


class CustomDisk(Model):
    """Power-law or viscous accretion disk model.
    Andrews et al. (2009), ApJ 700, 1502
    Kwon et al. (2015), ApJ 808, 102
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        #: Set parameters of the disk model
        self.parameter['gas_mass'] = 1e-3 * self.math.const['M_sun']
        self.parameter['grid_type'] = 'cylindrical'
        self.parameter['inner_radius'] = 0.1 * self.math.const['au']
        self.parameter['outer_radius'] = 100. * self.math.const['au']
        
        # In the case of a spherical grid
        self.spherical_parameter['n_r'] = 100
        self.spherical_parameter['n_th'] = 181
        self.spherical_parameter['n_ph'] = 1
        self.spherical_parameter['sf_r'] = 1.03
        # sf_th = 1 is sinus; sf_th > 1 is exp with step width sf_th; rest is linear
        self.spherical_parameter['sf_th'] = 1.0
        # In the case of a cylindrical grid
        self.cylindrical_parameter['n_r'] = 100
        self.cylindrical_parameter['n_z'] = 181
        self.cylindrical_parameter['n_ph'] = 1
        self.cylindrical_parameter['sf_r'] = 1.03
        # sf_z = -1 is using scale height; sf_z = 1 is sinus;
        # sf_z > 1 is exp with step width sf_z and rest is linear
        self.cylindrical_parameter['sf_z'] = -1
        
        # Default disk parameter
        self.parameter['ref_radius'] = 100. * self.math.const['au']
        self.parameter['ref_scale_height'] = 10. * self.math.const['au']
        self.parameter['beta'] = 1.1
        self.parameter['alpha'] = 3.0 * (self.parameter['beta'] - 0.5)
        self.parameter['tapered_gamma'] = None

        self.custom_parameter_list = [
            'ref_radius',
            'ref_scale_height',
            'beta',
            'alpha',
            'tapered_gamma',
        ]

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters.
        """
        # Use extra parameter to vary the disk structure
        if extra_parameter is None:
            return

        for i, param in enumerate(extra_parameter[::2]):
            if param not in self.custom_parameter_list:
                print(f'  - Warning: invalid parameter: {param}')
                continue
            try:
                self.parameter[param] = float(eval(extra_parameter[2*i+1]))
            except:
                self.parameter[param] = extra_parameter[2*i+1]
            print(f'  - Info: {param} updated to {self.parameter[param]}')

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Returns:
            float: Gas density at a given position.
        """
        gas_density = self.math.default_disk_density(self.position,
                                                     inner_radius=self.parameter['inner_radius'],
                                                     outer_radius=self.parameter['outer_radius'],
                                                     ref_radius=self.parameter['ref_radius'],
                                                     ref_scale_height=self.parameter['ref_scale_height'],
                                                     alpha=self.parameter['alpha'], beta=self.parameter['beta'],
                                                     tapered_gamma=self.parameter['tapered_gamma'])
        return gas_density

    def get_scale_height(self, radius):
        """Calculates the scale height at a certain position.

        Args:
            radius (float) : Cylindrical radius of current position

        Returns:
            float: Scale height.
        """
        scale_height = self.math.default_disk_scale_height(radius,
                                                           ref_radius=self.parameter['ref_radius'],
                                                           ref_scale_height=self.parameter['ref_scale_height'],
                                                           beta=self.parameter['beta'])
        return scale_height
