#!/usr/bin/env python
# -*- coding: utf-8 -*-

from polaris_tools_modules.math import Math
from polaris_tools_modules.base import Model
from polaris_tools_custom.model import *
from polaris_tools_modules.atmosphere import AtmosphereRoutines
import numpy as np


class ModelChooser:
    """The ModelChooser class provides the chosen model.
    """

    def __init__(self, parse_args):
        """Initialisation of all usable options.

        Notes:
            To create your own model, add its name to the dictionary
            and write a class with its options as a derived class of class Model.

        Args:
            parse_args (ArgumentParser) : Provides all parameters chosen
            by user when executing PolarisTools.
        """
        self.parse_args = parse_args

        # Get math module
        self.math = Math()

        # dict: Dictionary with all usable models
        self.model_dict = {
            'default': Model,
            'cloudy': Cloudy,
            'rayleigh': Rayleigh,
            'ringed': Ringed,
            'venus': Venus,
            'disk': Disk,
            'sphere': Sphere,
        }
        update_model_dict(self.model_dict)

    def get_module(self):
        """Chooses model class from user input

            Note:
                Parameters set by PolarisTools overwrite preset values in the
                separate model classes.

            Returns:
                Instance of chosen model.
        """
        if self.parse_args.model_name in self.model_dict.keys():
            model = self.model_dict[self.parse_args.model_name]()
        elif self.parse_args.model_name is not None:
            raise ValueError('Model name not known! You can add a new model in model.py.')
        else:
            model = self.model_dict['default']()

        # Set user input variables
        if 'grid_type' in vars(self.parse_args).keys():
            model.update_parameter(self.parse_args.extra_parameter)
            if self.parse_args.grid_type is not None:
                model.parameter['grid_type'] = self.parse_args.grid_type
            if self.parse_args.gas_mass is not None:
                model.parameter['gas_mass'] = self.math.parse(
                    self.parse_args.gas_mass, 'mass')
            if self.parse_args.inner_radius is not None:
                model.parameter['inner_radius'] = self.math.parse(
                    self.parse_args.inner_radius, 'length')
            if self.parse_args.outer_radius is not None:
                model.parameter['outer_radius'] = self.math.parse(
                    self.parse_args.outer_radius, 'length')
            if self.parse_args.z_max is not None:
                model.cylindrical_parameter['z_max'] = self.math.parse(
                    self.parse_args.z_max, 'length')
            if self.parse_args.n_r is not None:
                model.spherical_parameter['n_r'] = self.parse_args.n_r
                model.cylindrical_parameter['n_r'] = self.parse_args.n_r
            if self.parse_args.n_ph is not None:
                model.spherical_parameter['n_ph'] = self.parse_args.n_ph
                model.cylindrical_parameter['n_ph'] = self.parse_args.n_ph
            if self.parse_args.n_th is not None:
                model.spherical_parameter['n_th'] = self.parse_args.n_th
            if self.parse_args.n_z is not None:
                model.cylindrical_parameter['n_z'] = self.parse_args.n_z
            if self.parse_args.sf_r is not None:
                model.spherical_parameter['sf_r'] = self.parse_args.sf_r
                model.cylindrical_parameter['sf_r'] = self.parse_args.sf_r
            if self.parse_args.sf_ph is not None:
                model.spherical_parameter['sf_ph'] = self.parse_args.sf_ph
                model.cylindrical_parameter['sf_ph'] = self.parse_args.sf_ph
            if self.parse_args.sf_th is not None:
                model.spherical_parameter['sf_th'] = self.parse_args.sf_th
            if self.parse_args.sf_z is not None:
                model.cylindrical_parameter['sf_z'] = self.parse_args.sf_z
        elif 'distance' in vars(self.parse_args).keys():
            if self.parse_args.distance is not None:
                model.parameter['distance'] = self.math.parse(
                    self.parse_args.distance, 'length')
    
        # Set the grid extent if global extent is set
        if model.parameter['grid_type'] == 'octree' and model.parameter['outer_radius'] is not None:
            model.octree_parameter['sidelength'] = 2. * \
                model.parameter['outer_radius']
        elif model.parameter['grid_type'] == 'spherical':
            if model.parameter['inner_radius'] is not None:
                model.spherical_parameter['inner_radius'] = model.parameter['inner_radius']
            if model.parameter['outer_radius'] is not None:
                model.spherical_parameter['outer_radius'] = model.parameter['outer_radius']
        elif model.parameter['grid_type'] == 'cylindrical':
            if model.parameter['inner_radius'] is not None:
                model.cylindrical_parameter['inner_radius'] = model.parameter['inner_radius']
            if model.parameter['outer_radius'] is not None:
                model.cylindrical_parameter['outer_radius'] = model.parameter['outer_radius']
                if model.cylindrical_parameter['z_max'] is None:
                    model.cylindrical_parameter['z_max'] = model.parameter['outer_radius']

        return model


class Rayleigh(Model):
    def __init__(self):
        Model.__init__(self)
        self.ar = AtmosphereRoutines()

        # boundaries of the grid [m]
        # first value is the inner model radius
        # last value is the outer model radius
        self.spherical_parameter['radius_list'] = np.array([7.0e7, 7.01e7])

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

        self.parameter['optical_depth'] = 1

        self.parameter_list = [
            'optical_depth', # vertical optical depth of the atmosphere
        ]

    def update_parameter(self, extra_parameter):
        if extra_parameter is None:
            return

        for i, param in enumerate(extra_parameter[::2]):
            if param not in self.parameter_list:
                print(f'  - Warning: invalid parameter: {param}')
                continue
            try:
                self.parameter[param] = float(eval(extra_parameter[2*i+1]))
            except:
                self.parameter[param] = extra_parameter[2*i+1]
            print(f"  - Info: {param} updated to {self.parameter[param]}")

    def dust_density_distribution(self):
        # density: molecular hydrogen (input/cross_sections/molecular_hydrogen.dat)
        # calculate the number density based on a given optical depth, radial boundaries, and cross section
        density = self.ar.getNumberDensityFromOpticalDepth(
            self.parameter['optical_depth'],
            self.spherical_parameter['radius_list'][1] - self.spherical_parameter['radius_list'][0],
            self.ar.getRayleighCrossSection(5.5e-7, '1H2'))

        return density

    def dust_temperature(self):
        # return constant temperature
        return 300.0


class Cloudy(Model):
    def __init__(self):
        Model.__init__(self)
        self.ar = AtmosphereRoutines()

        self.planetary_radius = 7.0e7 # [m]
        self.gravity = 25 # [m/s^2]

        # pressure profile of the atmosphere (bottom to top)
        self.pressure_profile = np.geomspace(1e+1, 1e-5, 7) * 1e5  # bar -> Pa

        # constant temperature_profile [K]
        self.temperature_profile = 300 * np.ones_like(self.pressure_profile)

        # boundaries of the grid [m]
        # first value is the inner model radius
        # last value is the outer model radius
        self.spherical_parameter['radius_list'] = self.ar.getAltitudeFromPressure(
            self.pressure_profile,
            self.temperature_profile,
            self.ar.getMolarMass('1H2'),
            self.gravity)
        self.spherical_parameter['radius_list'] += self.planetary_radius

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

        self.parameter['optical_depth'] = 1

        self.parameter_list = [
            'optical_depth', # vertical optical depth of the cloud layer
        ]

    def update_parameter(self, extra_parameter):
        if extra_parameter is None:
            return

        for i, param in enumerate(extra_parameter[::2]):
            if param not in self.parameter_list:
                print(f'  - Warning: invalid parameter: {param}')
                continue
            try:
                self.parameter[param] = float(eval(extra_parameter[2*i+1]))
            except:
                self.parameter[param] = extra_parameter[2*i+1]
            print(f"  - Info: {param} updated to {self.parameter[param]}")

    def dust_density_distribution(self):
        density = np.zeros(2, dtype=float)
        # density[0]: molecular hydrogen (input/cross_sections/molecular_hydrogen.dat)
        # density[1]: clouds (input/refractive_indices/water_s81.nk)

        pos = np.linalg.norm(self.position)
        idx = np.searchsorted(self.spherical_parameter['radius_list'], pos)

        # calculate the number density assuming hydrostatic equilibrium and an ideal gas
        density[0] = self.ar.getNumberDensityFromPressureAltitude(
            self.pressure_profile[idx-1] - self.pressure_profile[idx],
            self.spherical_parameter['radius_list'][idx] - self.spherical_parameter['radius_list'][idx-1],
            self.ar.getMolarMass('1H2'),
            self.gravity)

        # add water clouds between 1 bar and 0.1 bar
        c_ext = 2.05786e-11  # precalculated cross section [m^2] using r_eff = 2.0e-06 m, veff = 0.1 (550 nm)
        if self.spherical_parameter['radius_list'][1] < pos < self.spherical_parameter['radius_list'][2]:
            density[1] = self.ar.getNumberDensityFromOpticalDepth(
                self.parameter['optical_depth'],
                self.spherical_parameter['radius_list'][2] - self.spherical_parameter['radius_list'][1],
                c_ext)

        return density
    
    def dust_temperature(self):
        # return temperature at self.position
        pos = np.linalg.norm(self.position)
        # interpolate temperature profile based on current position
        return np.interp(pos, self.spherical_parameter['radius_list'], self.temperature_profile)


class Ringed(Model):
    def __init__(self):
        Model.__init__(self)
        self.ar = AtmosphereRoutines()

        self.planetary_radius = 7.0e7 # [m]
        self.ring_inner_radius = 1.2 * self.planetary_radius
        self.ring_outer_radius = 2.3 * self.planetary_radius
        self.ring_opening_angle = 0.2 / 3600 # [deg]

        # boundaries of the grid [m]
        # first value is the inner model radius
        # last value is the outer model radius
        self.spherical_parameter['radius_list'] = np.array([self.planetary_radius, 7.01e7, self.ring_inner_radius, self.ring_outer_radius])

        # Use spherical coordinate system
        self.parameter['grid_type'] = 'spherical'
        # sf_r = 0: user defined radial boundaries
        self.spherical_parameter['sf_r'] = 0

        # set inner and outer grid radius (bottom and top of atmosphere)
        self.parameter['inner_radius'] = self.spherical_parameter['radius_list'][0]
        self.parameter['outer_radius'] = self.spherical_parameter['radius_list'][-1]

        # number of radial grid cells
        self.spherical_parameter['n_r'] = len(self.spherical_parameter['radius_list']) - 1

        # sf_th = 0: user defined theta cells
        self.spherical_parameter['sf_th'] = 0
        # custom theta cells [rad]
        self.spherical_parameter['theta_list'] = np.deg2rad( np.array([0, 90 - 0.5*self.ring_opening_angle, 90 + 0.5*self.ring_opening_angle, 180]) )
        # number of theta cells
        self.spherical_parameter['n_th'] = len(self.spherical_parameter['theta_list']) - 1

        # sf_ph = -1: equally distributed phi cells
        self.spherical_parameter['sf_ph'] = -1
        # number of phi cells
        self.spherical_parameter['n_ph'] = 1

        self.parameter['optical_depth_gas'] = 1
        self.parameter['optical_depth_ring'] = 1

        self.parameter_list = [
            'optical_depth_gas', # vertical optical depth of the atmosphere
            'optical_depth_ring', # vertical optical depth of the ring
        ]

    def update_parameter(self, extra_parameter):
        if extra_parameter is None:
            return

        for i, param in enumerate(extra_parameter[::2]):
            if param not in self.parameter_list:
                print(f'  - Warning: invalid parameter: {param}')
                continue
            try:
                self.parameter[param] = float(eval(extra_parameter[2*i+1]))
            except:
                self.parameter[param] = extra_parameter[2*i+1]
            print(f"  - Info: {param} updated to {self.parameter[param]}")

    def dust_density_distribution(self):
        density = np.zeros(2, dtype=float)
        # density[0]: molecular hydrogen (input/cross_sections/molecular_hydrogen.dat)
        # density[1]: ring particles (input/refractive_indices/silicate_d03.nk)
        
        pos = np.linalg.norm(self.position)
        theta = 0.0
    
        if pos > 0.0:
            theta = np.rad2deg( np.arccos(self.position[2] / pos) )

        if self.spherical_parameter['radius_list'][0] < pos < self.spherical_parameter['radius_list'][1]:
            # calculate the number density based on a given optical depth, radial boundaries, and cross section
            density[0] = self.ar.getNumberDensityFromOpticalDepth(
                self.parameter['optical_depth_gas'],
                self.spherical_parameter['radius_list'][1] - self.spherical_parameter['radius_list'][0],
                self.ar.getRayleighCrossSection(5.5e-7, '1H2'))

        if self.ring_inner_radius < pos < self.ring_outer_radius and np.abs(theta - 90.0) < 0.5 * self.ring_opening_angle:
            # vertical height of the ring at outer radius
            length = 2.0 * self.ring_outer_radius * np.tan(0.5 * np.deg2rad(self.ring_opening_angle))
            c_ext = 7.13263e-13  # precalculated cross section using q = -3, s_min = 1e-7 m, s_max = 1e-5 m (550 nm)
            density[1] = self.parameter['optical_depth_ring'] / (c_ext * length)

        return density
    
    def dust_temperature(self):
        # return constant temperature
        return 300.0
    

class Venus(Model):
    def __init__(self):
        Model.__init__(self)
        self.ar = AtmosphereRoutines()

        self.planetary_radius = 6.05e6 # [m]
        self.gravity = 8.9 # [m/s^2]

        # pressure profile [bar]
        # Venus Global Reference Atmospheric Model sample output
        self.pressure_profile = np.array([
            9.30e+01, 6.65e+01, 4.64e+01, 3.20e+01, 2.30e+01, 1.52e+01, 9.82e+00, 6.01e+00, 3.55e+00, 1.98e+00,
            1.06e+00, 5.34e-01, 2.30e-01, 8.42e-02, 3.16e-02, 1.18e-02, 4.28e-03, 1.38e-03, 3.89e-04, 1.03e-04,
            2.85e-05, 7.35e-06, 1.66e-06, 3.27e-07, 6.69e-08, 1.46e-08, 2.93e-09, 7.57e-10, 2.08e-10, 7.90e-11])
        # bar -> Pa
        self.pressure_profile *= 1e5

        # temperature_profile [K]
        # Venus Global Reference Atmospheric Model sample output
        self.temperature_profile = np.array([
            737, 703, 664, 626, 586, 545, 504, 461, 418, 382,
            347, 301, 248, 226, 231, 227, 211, 190, 175, 170,
            169, 160, 146, 140, 140, 139, 141, 142, 144, 146])

        # boundaries of the grid [m]
        # first value is the inner model radius
        # last value is the outer model radius
        self.spherical_parameter['radius_list'] = np.arange(0.0, 150e3, 5e3)
        self.spherical_parameter['radius_list'] += self.planetary_radius

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

        self.parameter['optical_depth'] = 30

        self.parameter_list = []

    def update_parameter(self, extra_parameter):
        if extra_parameter is None:
            return

        for i, param in enumerate(extra_parameter[::2]):
            if param not in self.parameter_list:
                print(f'  - Warning: invalid parameter: {param}')
                continue
            try:
                self.parameter[param] = float(eval(extra_parameter[2*i+1]))
            except:
                self.parameter[param] = extra_parameter[2*i+1]
            print(f"  - Info: {param} updated to {self.parameter[param]}")

    def dust_density_distribution(self):
        density = np.zeros(2, dtype=float)
        # density[0]: carbon dioxide (input/cross_sections/carbon_dioxide.dat)
        # density[1]: venus-like clouds (input/refractive_indices/venus_clouds_hh74.nk)

        pos = np.linalg.norm(self.position)
        idx = np.searchsorted(self.spherical_parameter['radius_list'], pos)

        # calculate the number density assuming hydrostatic equilibrium and an ideal gas
        density[0] = self.ar.getNumberDensityFromPressureAltitude(
            self.pressure_profile[idx-1] - self.pressure_profile[idx],
            self.spherical_parameter['radius_list'][idx] - self.spherical_parameter['radius_list'][idx-1],
            self.ar.getMolarMass('12C-16O2'),
            self.gravity)

        # add clouds between 50 km and 70 km
        idx_bot = np.searchsorted(self.spherical_parameter['radius_list'], self.planetary_radius + 50e3)
        idx_top = np.searchsorted(self.spherical_parameter['radius_list'], self.planetary_radius + 70e3)
        c_ext = 6.631955598e-12  # precalculated cross section [m^2] using r_eff = 1.05e-06 m, veff = 0.07 (550 nm)
        if self.spherical_parameter['radius_list'][idx_bot] < pos < self.spherical_parameter['radius_list'][idx_top]:
            density[1] = self.ar.getNumberDensityFromOpticalDepth(
                self.parameter['optical_depth'],
                self.spherical_parameter['radius_list'][idx_top] - self.spherical_parameter['radius_list'][idx_bot],
                c_ext)

        return density
    
    def dust_temperature(self):
        # return temperature at self.position
        pos = np.linalg.norm(self.position)
        # interpolate temperature profile based on current position
        return np.interp(pos, self.spherical_parameter['radius_list'], self.temperature_profile)


class Disk(Model):
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


class Sphere(Model):
    """A sphere model with constant density
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        #: Set parameters of the sphere model
        self.parameter['grid_type'] = 'spherical'
        self.parameter['inner_radius'] = 0.1 * self.math.const['au']
        self.parameter['outer_radius'] = 100. * self.math.const['au']
        
        self.spherical_parameter['n_r'] = 100
        self.spherical_parameter['n_th'] = 91
        self.spherical_parameter['n_ph'] = 91
        self.spherical_parameter['sf_r'] = 1.03
        self.parameter['gas_mass'] = 1e-4 * self.math.const['M_sun']
        
        self.parameter['mag_field_geometry'] = 'toroidal'
        self.parameter['mag_field_strength'] = 1e-10

        self.custom_parameter_list = [
            'mag_field_geometry',
            'mag_field_strength'
        ]

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Returns:
            float: Gas density at a given position.
        """
        gas_density = self.math.const_sphere_density(self.position,
                                                     outer_radius=self.parameter['outer_radius'],
                                                     inner_radius=self.parameter['inner_radius'])
        return gas_density

    def magnetic_field(self):
        """Calculates the magnetic field strength at a given position.

        Returns:
            List[float, float, float]: Magnetic field strength at a given position.
        """
        if self.parameter['mag_field_geometry'] == 'toroidal':
            magnetic_field = self.math.toroidal_mag_field(
                position=self.position, mag_field_strength=self.parameter['mag_field_strength'])
        elif self.parameter['mag_field_geometry'] == 'vertical':
            magnetic_field = self.math.simple_mag_field(
                mag_field_strength=self.parameter['mag_field_strength'], axis='z')
        elif self.parameter['mag_field_geometry'] == 'radial':
            magnetic_field = self.math.radial_mag_field(
                mag_field_strength=self.parameter['mag_field_strength'], position=self.position)
        else:
            magnetic_field = [0, 0, 0]
            print(f'  - Warning: invalid magnetic field geometry: {self.parameter["mag_field_geometry"]}')
        return magnetic_field

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters and update 
        model parameter that depend on other parameter.
        """
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
