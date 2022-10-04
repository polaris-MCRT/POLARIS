#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from polaris_tools_modules.math import Math
from polaris_tools_modules.base import Model

"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_model_dict(dictionary):
    model_dict = {
        'custom': CustomModel,
        'pascucci': Pascucci,
        'subdisk': SubDisk,
    }
    dictionary.update(model_dict)


class CustomModel(Model):
    """Change this to the model you want to use.
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        # Set parameters of the custom model (see parent Model class for all available options)
        self.parameter['distance'] = 140.0 * self.math.const['pc']
        self.parameter['grid_type'] = 'spherical'
        self.parameter['inner_radius'] = 0.1 * self.math.const['au']
        self.parameter['outer_radius'] = 100.0 * self.math.const['au']
        self.parameter['gas_mass'] = 1e-2 * self.math.const['M_sun']
        # Define which other choise are default for this model
        self.parameter['background_source'] = 'bg_plane'
        self.parameter['radiation_source'] = 't_tauri'
        self.parameter['dust_composition'] = 'mrn'
        self.parameter['gas_species'] = 'oh'
        self.parameter['detector'] = 'cartesian'

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters and update 
        disk parameter that depend on other parameter.
        """
        # Use extra_parameter to adjust the model without changing the model.py file

    def gas_density_distribution(self):
        """Define here your routine to calculate the density at a given position
        in the model space.

        Notes:
            Use 'self.position' to calculate the quantity depending on position.

            Define also the following routines if necessary:
                dust_density_distribution(self), gas_temperature(self),
                dust_temperature(self), velocity_field(self), magnetic_field(self),
                dust_id(self), dust_min_size(self), dust_max_size(self), dust_size_param(self)

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

class Pascucci(Model):
    """Pascucci benchmark disk 2004
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        # Set parameters of the pascucci model (see parent Model class for all available options)
        self.parameter['distance'] = 140.0 * self.math.const['pc']
        self.parameter['grid_type'] = 'spherical'
        self.parameter['inner_radius'] = 1 * self.math.const['au']
        self.parameter['outer_radius'] = 1000.0 * self.math.const['au']
        self.parameter['gas_mass'] = 1.1e-2 * self.math.const['M_sun']
        # Define which other choise are default for this model
        self.parameter['radiation_source'] = 't_tauri'
        self.parameter['dust_composition'] = 'mrn'
        self.parameter['gas_species'] = 'oh'
        self.parameter['detector'] = 'cartesian'
        
        # In the case of a spherical grid
        self.spherical_parameter['n_r'] = 1000
        self.spherical_parameter['n_th'] = 121
        self.spherical_parameter['n_ph'] = 1
        self.spherical_parameter['sf_r'] = 1.03

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters and update 
        disk parameter that depend on other parameter.
        """
        # Use extra_parameter to adjust the model without changing the model.py file

    def gas_density_distribution(self):
        """Define here your routine to calculate the density at a given position
        in the model space.

        Notes:
            Use 'self.position' to calculate the quantity depending on position.

        Returns:
            float: Gas density at a given position.
        """
        
        #: float: Parameters from Pascucci 2004
        r_d = self.parameter['outer_radius']/2
        z_d = self.parameter['outer_radius']/8
        
        #: float: Radial distance from center and z
        radius = np.sqrt(self.position[0] ** 2 + self.position[1]
                         ** 2 + self.position[2] ** 2)
        z = self.position[2]
        
        #: float:More parameters from Pascucci 2004
        if radius > 0:
            f1 = (radius/r_d) ** -1.0
            h = z_d * ((radius/r_d) ** 1.125)
            f2 = np.exp(-np.pi/4 * ((z/h) ** 2))
            gas_density = f1 * f2
        else:
            gas_density = 0
        
        return gas_density

class SubDisk(Model):
    """
    Model for disks with sublimated inner radii. No normalization needed.
    Dust mass is decreasing with increasing sublimation radius.
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        # Set parameters of the custom model (see parent Model class for all available options)
        self.parameter['distance'] = 140.0 * self.math.const['pc']
        self.parameter['grid_type'] = 'cylindrical'
        self.parameter['inner_radius'] = 0.06 * self.math.const['au']
        self.parameter['outer_radius'] = 20.0 * self.math.const['au']
        self.parameter['stable_radius'] = 0.1 * self.math.const['au']
        self.parameter['gas_mass'] = 1e-3 * self.math.const['M_sun']
        # Define which other choise are default for this model
        self.parameter['radiation_source'] = 't_tauri'
        self.parameter['dust_composition'] = 'mrn'
        self.parameter['gas_species'] = 'oh'
        self.parameter['detector'] = 'cartesian'
        # In the case of a spherical grid
        self.spherical_parameter['n_r'] = 200
        self.spherical_parameter['n_th'] = 181
        self.spherical_parameter['n_ph'] = 1
        self.spherical_parameter['sf_r'] = 1.03
        # sf_th = 1 is sinus; sf_th > 1 is exp with step width sf_th; rest is linear
        self.spherical_parameter['sf_th'] = 1.0
        # In the case of a cylindrical grid
        self.cylindrical_parameter['n_r'] = 200
        self.cylindrical_parameter['n_z'] = 181
        self.cylindrical_parameter['n_ph'] = 1
        self.cylindrical_parameter['sf_r'] = 1.03
        # sf_z = -1 is using scale height; sf_z = 1 is sinus;
        # sf_z > 1 is exp with step width sf_z and rest is linear
        self.cylindrical_parameter['sf_z'] = -1
        # Default disk parameter
        self.parameter['ref_radius'] = 100. * self.math.const['au']
        self.parameter['ref_scale_height'] = 10. * self.math.const['au']
        self.parameter['alpha'] = 1.8 # Andrews 2010 with alpha-beta and formula below
        # self.parameter['beta'] = self.parameter['alpha'] / 3 + 0.5
        self.parameter['beta'] = 1.1 # Woitke 2019
        self.parameter['f_dg'] = 0.01
        self.parameter['zoom'] = 5

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters and update 
        disk parameter that depend on other parameter.
        """
        # Set smallest stable radius (e.g. Akeson 2005 12.5Rsun)
        if extra_parameter is not None:
            self.parameter['stable_radius'] = self.math.parse(extra_parameter[0], 'length')

    def gas_density_distribution(self):
        """Define here your routine to calculate the density at a given position
        in the model space.
        Args:
            position (List[float, float, float]): Position in model space.
            inner_radius (float): Inner radius of the disk.
            outer_radius (float): Outer radius of the disk.
            sub_radius (float): Sublimation radius of the disk.
            ref_scale_height (float): Reference scale height.
            ref_radius (float): Reference radius.
            gas_mass (float): Total gas mass.
            alpha (float): Exponent for radial density decrease.
            beta (float): Exponent for disk flaring.
            f_dg (float): Dust to gas ratio.
        Notes:
            Use 'self.position' to calculate the quantity depending on position.
        Returns:
            float: Gas density at a given position.
        """
        #: float: Total dust mass
        m_dust = self.parameter['gas_mass'] * self.parameter['f_dg']
        #: float: Cylindrical radius
        radius_cy = np.sqrt(self.position[0] ** 2 + self.position[1] ** 2)
        if self.parameter['outer_radius']*self.parameter['zoom'] >= radius_cy >= self.parameter['inner_radius']:
            #: float: Vertical height
            vert_height = abs(self.position[2])
            #: float: Vertical scale height
            scale_height = self.scale_height(radius_cy)
            #: float gas surface density scaling factor
            sig_0 = (m_dust * ( 2 - self.parameter['alpha'])) / (self.parameter['f_dg'] * 2 * np.pi * self.parameter['ref_radius'] ** self.parameter['alpha'] * \
                    ((self.parameter['outer_radius']*self.parameter['zoom']) ** (2 - self.parameter['alpha']) - self.parameter['stable_radius'] ** (2 - self.parameter['alpha'])))
            #: float gas surface density
            sigma = sig_0 * (radius_cy / self.parameter['ref_radius']) ** (-self.parameter['alpha'])
            #: float: Shakura and Sunyaev like density distribution
            density =  sigma / (np.sqrt(2 * np.pi) * scale_height ) * \
                np.exp(-0.5 * (vert_height / scale_height) ** 2)
        else:
            density = 0.

        return density
    
    def scale_height(self, radius):
        """Calculates the scale height at a certain position.

        Args:
            radius (float) : Cylindrical radius of current position

        Returns:
            float: Scale height.
        """
        scale_height = self.parameter['ref_scale_height'] * (radius / self.parameter['ref_radius']) ** self.parameter['beta']
        
        return scale_height
