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
