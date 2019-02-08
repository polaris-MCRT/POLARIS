#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from custom.source import update_sources_dict
from modules.base import StellarSource
from modules.math import Math


class SourceChooser:
    """The StellarSourceChooser class provides the chosen stellar source
    such as a single star in the centre of the model space.
    """

    def __init__(self, model, file_io, parse_args):
        """Initialisation of all usable options.

        Notes:
            To create your own stellar source, add its name to the dictionary
            and write a class with its options as a derived class of class StellarSource.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
            file_io : Handles file input/output and all
                necessary paths.
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
        """
        self.model = model
        self.file_io = file_io
        self.parse_args = parse_args

        # Get math module
        self.math = Math()

        #: dict: Dictionary with all usable stellar sources
        self.sources_dict = {
            'isrf': ISRF,
            'mathis_isrf': MathisISRF,
            't_tauri': TTauri,
            'herbig_ae': HerbigAe,
            'binary': BinaryStar,
            'sun': Sun,
        }
        update_sources_dict(self.sources_dict)

    def get_module(self):
        """Chooses radiation source from user input

        Note:
            Parameters set by PolarisTools overwrite preset values in the
            separate radiation source classes.

        Returns:
            Instance of chosen stellar source.
        """
        if self.parse_args.radiation_source in self.sources_dict.keys():
            source_name = self.parse_args.radiation_source
        elif self.parse_args.radiation_source is None:
            if self.model is not None and self.model.parameter['radiation_source'] is not None:
                source_name = self.model.parameter['radiation_source']
            else:
                return None
        else:
            raise ValueError(
                'stellar source not known! You can add a new source in source.py')
        radiation_source = self.sources_dict[source_name](
            self.file_io, self.parse_args)
        # Overwrite default values with user input
        if self.parse_args.nr_photons is not None:
            radiation_source.parameter['nr_photons'] = int(
                self.parse_args.nr_photons)
        if self.parse_args.rad_source_position is not None:
            radiation_source.parameter['position'] = [self.math.parse(self.parse_args.rad_source_position[i], 'length')
                                                      for i in range(3)]
        if self.parse_args.rad_source_temperature is not None:
            radiation_source.parameter['temperature'] = self.parse_args.rad_source_temperature
        if self.parse_args.rad_source_luminosity is not None:
            radiation_source.parameter['luminosity'] = self.math.parse(
                self.parse_args.rad_source_luminosity, 'luminosity')
            radiation_source.parameter['radius'] = 0.
        elif self.parse_args.rad_source_radius is not None:
            radiation_source.parameter['radius'] = self.math.parse(
                self.parse_args.rad_source_radius, 'length')
        if self.parse_args.rad_source_mass is not None:
            radiation_source.parameter['mass'] = self.math.parse(
                self.parse_args.rad_source_mass, 'mass')
        radiation_source.update_parameter(
            self.parse_args.rad_source_extra_parameter)
        return radiation_source

    def get_module_from_name(self, radiation_source_name):
        """Chooses radiation source from name string

        Args:
            radiation_source_name (str): Name of the stellar source.

        Returns:
            Instance of chosen radiation source .
        """
        if radiation_source_name in self.sources_dict.keys():
            radiation_source = self.sources_dict[radiation_source_name](
                self.file_io, self.parse_args)
        else:
            raise ValueError('No source with the name ' +
                             str(radiation_source_name) + '!')
        return radiation_source


class ISRF:
    """The ISRF class is the interstellar radiation field, usable for temperature
    and rat simulations. The data is read from file in the input directory.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        self.file_io = file_io
        self.parse_args = parse_args

        #: dict: Parameters of the interstellar radiation field class
        self.parameter = {
            'filename': 'interstellar_radiation_field.dat',
            'nr_photons': 1e6,
            'kepler_usable': False,
        }

    def get_command_line(self):
        """Provides ISRF command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the interstellar radiation field.
        """
        return '\t<source_isrf nr_photons = "' \
               + str(int(self.parameter['nr_photons'])) + '">\t''\"' \
               + self.file_io.path['input'] \
               + self.parameter['filename'] + '\"\n'

    def get_command(self):
        """Provides ISRF command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the interstellar radiation field.
        """
        return self.get_command_line()

    def update_parameter(self, extra_parameter):
        """Use this function to set isrf source parameter with the extra parameters and update isrf source parameter that depend on other parameter.
        """
        # Use extra_parameter to adjust the radiation source without changing the source.py file


class MathisISRF(ISRF):
    """The ISRF class is the interstellar radiation field, usable for temperature
    and rat simulations. The data is read from file in the input directory.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the radiation source parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        ISRF.__init__(self, file_io, parse_args)

        self.file_io = file_io
        self.parse_args = parse_args

        #: dict: Parameters of the interstellar radiation field class
        mathis_isrf_parameter = {
            'g_zero': 1,
        }

        # Updates the parameter dictionary
        self.parameter.update(mathis_isrf_parameter)

    def get_command_line(self):
        """Provides ISRF command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the interstellar radiation field.
        """
        return '\t<source_isrf nr_photons = "' \
               + str(int(self.parameter['nr_photons'])) + '">\t' \
               + str(self.parameter['g_zero']) + '\n'

    def get_command(self):
        """Provides ISRF command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the interstellar radiation field.
        """
        return self.get_command_line()

    def update_parameter(self, extra_parameter):
        """Use this function to set isrf source parameter with the extra parameters and update isrf source parameter that depend on other parameter.
        """
        # Use extra_parameter to adjust the radiation source without changing the source.py file


class Sun(StellarSource):
    """The Sun class is a sun-like star in the center.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the radiation source parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        StellarSource.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        star_parameter = {
            # in SI units (m)
            'position': [0, 0, 0],
            'temperature': 5778,
            # In R_sun
            'radius': 1. * self.math.const['R_sun'],
            'mass': 1. * self.math.const['M_sun'],
        }

        # Updates the parameter dictionary
        self.parameter.update(star_parameter)


class TTauri(StellarSource):
    """Change this to the star you want to use.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the radiation source parameters.

        Args:
            file_io : Handles file input/output and all necessary paths.
        """
        StellarSource.__init__(self, file_io, parse_args)

        # Position of the star [m, m, m]
        self.parameter['position'] = [0, 0, 0]
        # Effective temperature of the star [K]
        self.parameter['temperature'] = 4000
        # Radius of the star [R_sun]
        self.parameter['radius'] = 0.9 * self.math.const['R_sun']
        # Mass of the star [M_sun] (for Keplerian rotation)
        self.parameter['mass'] = 0.7 * self.math.const['M_sun']
        # Number of photons if no number is chosen via --photons
        self.parameter['nr_photons'] = 1e6


class HerbigAe(StellarSource):
    """Change this to the star you want to use.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the radiation source parameters.

        Args:
            file_io : Handles file input/output and all necessary paths.
        """
        StellarSource.__init__(self, file_io, parse_args)

        # Position of the star [m, m, m]
        self.parameter['position'] = [0, 0, 0]
        # Effective temperature of the star [K]
        self.parameter['temperature'] = 8500
        # Radius of the star [R_sun]
        self.parameter['radius'] = 2.0 * self.math.const['R_sun']
        # Mass of the star [M_sun]
        self.parameter['mass'] = 5.0 * self.math.const['M_sun']
        # Number of photons if no number is chosen via --photons
        self.parameter['nr_photons'] = 1e6


class BinaryStar(StellarSource):
    """The BinaryStar class is two pre-main sequence stars in the center.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the radiation source parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        StellarSource.__init__(self, file_io, parse_args)

        self.parameter['nr_photons'] = 1e5

        from math import fmod
        a = 5.0  # AU
        omega = a ** (-1.5)

        # Number of the MHD-simulations
        num = 340.0  # 260.0

        g_time = 561.985178 * num
        angle_1 = fmod(omega * g_time, np.pi * 2.0)
        angle_2 = fmod(omega * g_time + np.pi, np.pi * 2.0)

        #: dict: Parameters which are different to the default values
        star_parameter = {
            'temperature': 4500.0,
            'radius': 0.4123886893502353 * self.math.const['R_sun'],
            'mass': 0.5 * self.math.const['M_sun'],
            'position_star_1': [a * self.math.const['au'] * np.cos(angle_1),
                                a * self.math.const['au'] * np.sin(angle_1),
                                0.0],
            'position_star_2': [a * self.math.const['au'] * np.cos(angle_2),
                                a * self.math.const['au'] * np.sin(angle_2),
                                0.0],
        }

        # Updates the parameter dictionary
        self.parameter.update(star_parameter)

    def get_command(self):
        new_command_line = str()
        self.parameter['position'] = self.parameter['position_star_1']
        new_command_line += self.get_command_line()
        self.parameter['position'] = self.parameter['position_star_2']
        new_command_line += self.get_command_line()
        return new_command_line
