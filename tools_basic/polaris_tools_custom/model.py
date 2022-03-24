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
        # 'cube': Cube,
        'custom': CustomModel,
        'custom_disk': CustomDisk,
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
        # self.parameter['distance'] = 140.0 * self.math.const['pc']
        self.parameter['grid_type'] = 'spherical'
        self.parameter['inner_radius'] = 0.1 * self.math.const['au']
        self.parameter['outer_radius'] = 100.0 * self.math.const['au']
        self.parameter['gas_mass'] = 1e-2 * self.math.const['M_sun']
        # Define which other choice are default for this model
        # self.parameter['background_source'] = 'bg_plane'
        # self.parameter['radiation_source'] = 't_tauri'
        # self.parameter['dust_composition'] = 'mrn'
        # self.parameter['gas_species'] = 'oh'
        # self.parameter['detector'] = 'cartesian'

    # def update_parameter(self, extra_parameter):
    #     """Use this function to set model parameter with the extra parameters and update
    #     model parameter that depend on other parameter.
    #     """
    #     # Use extra_parameter to adjust the model without changing the model.py file

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


class CustomDisk(Model):
    """The disk model with the Shakura and Sunyaev disk density profile.
    """

    def __init__(self):
        """Initialisation of the model parameters.

        Notes:
            Shakura and Sunyaev (1973)
            Link: http://adsabs.harvard.edu/abs/1973A&A....24..337S
        """
        Model.__init__(self)

        #: Set parameters of the disk model
        # self.parameter['distance'] = 140.0 * self.math.const['pc']
        self.parameter['gas_mass'] = 1e-3 * self.math.const['M_sun']
        self.parameter['grid_type'] = 'cylindrical'
        self.parameter['inner_radius'] = 0.1 * self.math.const['au']
        self.parameter['outer_radius'] = 100. * self.math.const['au']
        # Define the used sources, dust composition and gas species
        # self.parameter['radiation_source'] = 't_tauri'
        # self.parameter['dust_composition'] = 'mrn'
        # self.parameter['gas_species'] = 'co'
        # self.parameter['detector'] = 'cartesian'
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
        self.parameter['alpha'] = 2.625
        self.parameter['beta'] = 1.125

    # def update_parameter(self, extra_parameter):
    #     """Use this function to set model parameter with the extra parameters.
    #     """
    #     # Use extra parameter to vary the disk structure
    #     if extra_parameter is not None:
    #         if len(extra_parameter) == 4:
    #             self.parameter['ref_radius'] = self.math.parse(
    #                 extra_parameter[0], 'length')
    #             self.parameter['ref_scale_height'] = self.math.parse(
    #                 extra_parameter[1], 'length')
    #             self.parameter['alpha'] = float(extra_parameter[2])
    #             self.parameter['beta'] = float(extra_parameter[3])

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
                                                     alpha=self.parameter['alpha'], beta=self.parameter['beta'])
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

# class Cube(Model):
#     """A cube model with constant density
#     """

#     def __init__(self):
#         """Initialisation of the model parameters.
#         """
#         Model.__init__(self)

#         #: Set parameters of the sphere model
#         self.parameter['distance'] = 140.0 * self.math.const['pc']
#         # 2.8e-14 * self.math.const['M_sun']
#         self.parameter['gas_mass'] = np.array(
#             [[0.67 * 1e-6 * self.math.const['M_sun']], [0.33 * 1e-6 * self.math.const['M_sun']]])
#         self.parameter['outer_radius'] = 100.0 * \
#             self.math.const['au']  # 0.5 * self.math.const['au']
#         # self.parameter['radiation_source'] = 'isrf'
#         # self.parameter['dust_composition'] = 'silicate_oblate'
#         # self.parameter['detector'] = 'cartesian'

#     def dust_temperature(self):
#         """Calculates the dust temperature at a given position.

#         Returns:
#             float: Dust temperature at a given position.
#         """
#         dust_temperature = 20.
#         return dust_temperature

#     def gas_temperature(self):
#         """Calculates the dust temperature at a given position.

#         Returns:
#             float: Dust temperature at a given position.
#         """
#         gas_temperature = 10.
#         return gas_temperature

#     def gas_density_distribution(self):
#         """Calculates the gas density at a given position.

#         Returns:
#             float: Gas density at a given position.
#         """
#         # gas_density = self.math.bonor_ebert_density(self.position, outer_radius=self.parameter['outer_radius'],
#         #                                        truncation_radius=1 * self.math.const['au'])
#         # gas_density = self.math.random_density_distribution(
#         #    self.position, d_exp=2)
#         # gas_density = 1.0
#         return np.array([[1.], [1.]])

#     def magnetic_field(self):
#         """Calculates the magnetic field strength at a given position.

#         Returns:
#             List[float, float, float]: Magnetic field strength at the given
#             position.
#         """
#         return self.math.simple_mag_field(mag_field_strength=1e-10, axis='z')
