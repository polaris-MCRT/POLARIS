#!/usr/bin/env python
# -*- coding: utf-8 -*-

from polaris_tools_modules.math import Math
from polaris_tools_modules.base import Model
from polaris_tools_custom.model import *


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
            'constantCylinder': ConstantCylinder,
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
            raise ValueError(
                'Model name not known! You can add a new model in model.py.')
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
        # Convert radius in various units
        if model.parameter['grid_type'] == 'octree':
            model.tmp_parameter['radius_x_m'] = model.octree_parameter['sidelength'] / 2.
            model.tmp_parameter['radius_y_m'] = model.octree_parameter['sidelength'] / 2.
        elif model.parameter['grid_type'] == 'spherical':
            model.tmp_parameter['radius_x_m'] = model.spherical_parameter['outer_radius']
            model.tmp_parameter['radius_y_m'] = model.spherical_parameter['outer_radius']
        elif model.parameter['grid_type'] == 'cylindrical':
            model.tmp_parameter['radius_x_m'] = model.cylindrical_parameter['outer_radius']
            model.tmp_parameter['radius_y_m'] = model.cylindrical_parameter['outer_radius']
        else:
            raise ValueError('Grid type not known!')
        model.tmp_parameter['radius_x_arcsec'] = self.math.length_conv(
            model.tmp_parameter['radius_x_m'], 'arcsec', model.parameter['distance'])
        model.tmp_parameter['radius_y_arcsec'] = self.math.length_conv(
            model.tmp_parameter['radius_y_m'], 'arcsec', model.parameter['distance'])
        model.tmp_parameter['radius_x_pc'] = self.math.length_conv(
            model.tmp_parameter['radius_x_m'], 'pc')
        model.tmp_parameter['radius_y_pc'] = self.math.length_conv(
            model.tmp_parameter['radius_y_m'], 'pc')
        model.tmp_parameter['radius_x_au'] = self.math.length_conv(
            model.tmp_parameter['radius_x_m'], 'au')
        model.tmp_parameter['radius_y_au'] = self.math.length_conv(
            model.tmp_parameter['radius_y_m'], 'au')
        model.tmp_parameter['radius_x_ae'] = model.tmp_parameter['radius_x_au']
        model.tmp_parameter['radius_y_ae'] = model.tmp_parameter['radius_y_au']
        return model


class ConstantCylinder(Model):
    """The model of the polarimetrix experiment with a cylindrical dust cloud.
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        #: Set parameters of the disk model (SI units)
        self.parameter['distance'] = 0.5
        self.parameter['grid_type'] = 'cylindrical'
        self.parameter['inner_radius'] = 1e-5
        self.parameter['outer_radius'] = 0.03
        self.parameter['density'] = 1e13  # [m-3]
        self.parameter['num_dens'] = 1
        # Define the grid geometry
        self.cylindrical_parameter['n_r'] = 3
        self.cylindrical_parameter['n_z'] = 3
        self.cylindrical_parameter['n_ph'] = 1
        self.cylindrical_parameter['sf_r'] = -1
        self.cylindrical_parameter['sf_z'] = 0.5
        self.cylindrical_parameter['z_max'] = 0.015

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters.
        """
        # Use extra parameter to vary the disk structure
        if extra_parameter is not None and len(extra_parameter) == 1:
            self.parameter['density'] = float(extra_parameter[0])
        else:
            print('Only one extra parameter (constant mass density) expected.')

    def dust_density_distribution(self):
        """Calculates the gas mass density at a given position.
        The dust mass density is set with the same distribution but scaled
        with self.parameter['dust_gas_ratio'].

        Returns:
            float: Gas density at a given position.
        """
        gas_density = self.parameter['density']
        return gas_density