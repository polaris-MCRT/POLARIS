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
        'cube': Cube,
        'custom': CustomModel,
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
        # Define which other choice are default for this model
        self.parameter['background_source'] = 'bg_plane'
        self.parameter['radiation_source'] = 't_tauri'
        self.parameter['dust_composition'] = 'mrn'
        self.parameter['gas_species'] = 'oh'
        self.parameter['detector'] = 'cartesian'

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


class Cube(Model):
    """A cube model with constant density
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        #: Set parameters of the sphere model
        self.parameter['distance'] = 140.0 * self.math.const['pc']
        # 2.8e-14 * self.math.const['M_sun']
        self.parameter['gas_mass'] = np.array(
            [[0.67 * 1e-6 * self.math.const['M_sun']], [0.33 * 1e-6 * self.math.const['M_sun']]])
        self.parameter['outer_radius'] = 100.0 * \
            self.math.const['au']  # 0.5 * self.math.const['au']
        #self.parameter['radiation_source'] = 'isrf'
        self.parameter['dust_composition'] = 'silicate_oblate'
        self.parameter['detector'] = 'cartesian'

    def dust_temperature(self):
        """Calculates the dust temperature at a given position.

        Returns:
            float: Dust temperature at a given position.
        """
        dust_temperature = 20.
        return dust_temperature

    def gas_temperature(self):
        """Calculates the dust temperature at a given position.

        Returns:
            float: Dust temperature at a given position.
        """
        gas_temperature = 10.
        return gas_temperature

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Returns:
            float: Gas density at a given position.
        """
        # gas_density = self.math.bonor_ebert_density(self.position, outer_radius=self.parameter['outer_radius'],
        #                                        truncation_radius=1 * self.math.const['au'])
        # gas_density = self.math.random_density_distribution(
        #    self.position, d_exp=2)
        # gas_density = 1.0
        return np.array([[1.], [1.]])

    def magnetic_field(self):
        """Calculates the magnetic field strength at a given position.

        Returns:
            List[float, float, float]: Magnetic field strength at the given
            position.
        """
        return self.math.simple_mag_field(mag_field_strength=1e-10, axis='z')

