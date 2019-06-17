#!/usr/bin/env python
# -*- coding: utf-8 -*-

from polaris_tools_modules.base import Dust, NoDust
from polaris_tools_custom.dust import *
import numpy as np


class DustChooser:
    """The DustChooser class provides the composition of the chosen dust
    species.
    """

    def __init__(self, file_io, parse_args, model=None):
        """Initialisation of all usable options.

        Notes:
            To create your own dust composition, add its name to the dictionary
            and write a class with its options as a derived class of class Dust.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
            model: Handles the model space including various
                quantities such as the density distribution.
        """
        self.file_io = file_io
        self.parse_args = parse_args
        self.model = model

        #: dict: Dictionary with all usable dust compositions
        self.dust_dict = {
            'none': NoDust,
            'graphite_oblate': GraphiteOblate,
            'graphite_perpend': GraphitePerpend,
            'graphite_parallel': GraphiteParallel,
            'silicate_oblate': SilicateOblate,
            'silicate': Silicate,
            'mrn_oblate': MrnOblate,
            'mrn': Mrn,
            'CM20': CM20,
            'aPyM5': aPyM5,
            'aOlM5': aOlM5,
            'themis': Themis,
            'pah_neutral': PAH0,
            'pah_ion': PAH1,
        }
        update_dust_dict(self.dust_dict)

    def get_module(self):
        """Chooses dust class from user input

        Note:
            Parameters set by PolarisTools overwrite preset values in the
            separate dust classes.

        Returns:
            Instance of chosen dust composition.
        """
        if self.parse_args.dust_composition in self.dust_dict.keys():
            dust_composition = self.parse_args.dust_composition
        elif self.parse_args.dust_composition is None:
            if self.parse_args.simulation_type in ['line', 'zeeman']:
                dust_composition = 'none'
            elif self.model is not None and self.model.parameter['dust_composition'] is not None:
                dust_composition = self.model.parameter['dust_composition']
            else:
                dust_composition = 'none'
        else:
            raise ValueError(
                'dust component not known! You can add a new component in dust.py')

        dust = self.dust_dict[dust_composition](self.file_io, self.parse_args)
        # Overwrite preset parameters from user input
        self.update_user_parameters(dust)
        return dust

    def get_module_from_name(self, dust_name, user_input=True):
        """Chooses dust class from user input

        Args:
            user_input (bool): Update parameters with user input?

        Note:
            Parameters set by PolarisTools overwrite preset values in the
            separate dust classes.

        Returns:
            Instance of chosen dust composition.
        """
        if dust_name not in self.dust_dict.keys():
            raise ValueError(
                'dust component not known! You can add a new component in dust.py')
        dust = self.dust_dict[dust_name](self.file_io, self.parse_args)
        # Overwrite preset parameters from user input
        if user_input:
            self.update_user_parameters(dust)
        return dust

    def update_user_parameters(self, dust):
        """Overwrite preset parameters from user input.

        Args:
            dust : Instance of chosen dust composition.
        """
        if self.parse_args.dust_size is not None:
            dust.parameter['amin'] = self.math.parse(
                self.parse_args.dust_size[0], 'length')
            dust.parameter['amax'] = self.math.parse(
                self.parse_args.dust_size[1], 'length')
        if self.parse_args.dust_size_distribution is not None:
            if len(self.parse_args.dust_size_distribution) < 2:
                raise ValueError('Defined dust size distribution not usable!')
            else:
                dust.parameter['size_keyword'] = self.parse_args.dust_size_distribution[0]
                if 'plaw' in dust.parameter['size_keyword']:
                    min_len = 1
                    if '-ed' in dust.parameter['size_keyword']:
                        min_len += 3
                    if '-cv' in dust.parameter['size_keyword']:
                        min_len += 3
                elif 'logn' in dust.parameter['size_keyword']:
                    min_len = 2
                else:
                    min_len = 0
                if min_len != len(self.parse_args.dust_size_distribution) - 1:
                    raise ValueError(
                        'Defined dust size distribution not usable!')
                else:
                    dust.parameter['size_parameter'] = []
                    for i in range(min_len):
                        dust.parameter['size_parameter'].append(self.math.parse(
                            self.parse_args.dust_size_distribution[i + 1], 'length'))


class GraphiteOblate(Dust):
    """Dust class for graphite grains.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        # Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'graphite_oblate.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 2250
        # Minimum dust grain size
        self.parameter['amin'] = 5e-9
        # Maximum dust grain size
        self.parameter['amax'] = 2e-6


class GraphitePerpend(Dust):
    """Dust class for perpendicular oriented graphite grains.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'graphite_perpend.nk'
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 2250
        # Minimum dust grain size
        self.parameter['amin'] = 5e-9
        # Maximum dust grain size
        self.parameter['amax'] = 250e-9


class GraphiteParallel(Dust):
    """Dust class for  parallel oriented graphite grains.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'graphite_parallel.nk'
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 2250
        # Minimum dust grain size
        self.parameter['amin'] = 5e-9
        # Maximum dust grain size
        self.parameter['amax'] = 250e-9


class SilicateOblate(Dust):
    """Dust class for silicate grains."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        # Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'silicate_oblate.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 3800
        # Minimum dust grain size
        self.parameter['amin'] = 5e-9
        # Maximum dust grain size
        self.parameter['amax'] = 2e-6


class Silicate(Dust):
    """Dust class for silicate grains."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        # For creation of a dust catalog
        self.parameter['dust_cat_file'] = 'silicate.nk'
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 3500.
        # Minimum dust grain size
        self.parameter['amin'] = 5e-9
        # Maximum dust grain size
        self.parameter['amax'] = 250e-9


class MrnOblate(Dust):
    """Dust class for MRN dust grains.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Note:
            Mathis, J. S. 1972, ApJ, 176, 651
            Link: http://adsabs.harvard.edu/doi/10.1086/151667

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        # Set the name of one component to allow print of sizes and wavelengths
        self.parameter['dust_cat_file'] = 'silicate_oblate.dat'

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the MRN dust composition.
         """
        new_command_line = str()
        dust = self.dust_chooser.get_module_from_name('silicate_oblate')
        dust.parameter['fraction'] = 0.625
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('graphite_oblate')
        dust.parameter['fraction'] = 0.375
        new_command_line += dust.get_command_line()
        return new_command_line


class Mrn(Dust):
    """Dust class for MRN dust grains.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Note:
            Mathis, J. S. 1972, ApJ, 176, 651
            Link: http://adsabs.harvard.edu/doi/10.1086/151667

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['scattering'] = 'MIE'
        # Set the name of one component to allow print of sizes and wavelengths
        self.parameter['dust_cat_file'] = 'silicate.nk'
        # Minimum dust grain size
        self.parameter['amin'] = 5e-9
        # Maximum dust grain size
        self.parameter['amax'] = 250e-9

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the MRN dust composition.
         """
        new_command_line = str()
        dust = self.dust_chooser.get_module_from_name('silicate')
        dust.parameter['fraction'] = 0.625
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('graphite_perpend')
        dust.parameter['fraction'] = 0.25
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('graphite_parallel')
        dust.parameter['fraction'] = 0.125
        new_command_line += dust.get_command_line()
        return new_command_line


class CM20(Dust):
    """Dust class for carbon grains from Themis model."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'CM20.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 1600
        self.parameter['subl_temp'] = 2100
        self.parameter['size_keyword'] = 'plaw-ed'
        self.parameter['size_parameter'] = [-5, 10e-9, 50e-9,  1.0]
        # Minimum dust grain size
        self.parameter['amin'] = 0.4e-9
        # Maximum dust grain size
        self.parameter['amax'] = 4.9e-06


class aPyM5(Dust):
    """Dust class for carbon grains from Themis model."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'aPyM5.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 2190
        self.parameter['subl_temp'] = 1200
        self.parameter['size_keyword'] = 'logn'
        self.parameter['size_parameter'] = [8e-9, 1.0]
        # Minimum dust grain size (5 nm)
        self.parameter['amin'] = 1e-9
        # Maximum dust grain size
        self.parameter['amax'] = 4.9e-06


class aOlM5(Dust):
    """Dust class for carbon grains from Themis model."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'aOlM5.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 2190
        self.parameter['subl_temp'] = 1200
        self.parameter['size_keyword'] = 'logn'
        self.parameter['size_parameter'] = [8e-9, 1.0]
        # Minimum dust grain size
        self.parameter['amin'] = 1e-9
        # Maximum dust grain size
        self.parameter['amax'] = 4.9e-06


class Themis(Dust):
    """Dust class for themis dust grains.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Note:
            Link: https://www.ias.u-psud.fr/themis/

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 0
        # Set the name of one component to allow print of sizes and wavelengths
        self.parameter['dust_cat_file'] = 'CM20.dat'

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the MRN dust composition.
        """
        new_command_line = str()
        dust = self.dust_chooser.get_module_from_name('CM20')
        dust.parameter['fraction'] = 0.229
        # If more than the fraction is changed, overwrite parameter with user input
        self.dust_chooser.update_user_parameters(dust)
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('CM20')
        dust.parameter['material_density'] = 1570
        dust.parameter['size_keyword'] = 'logn'
        dust.parameter['size_parameter'] = [7e-9, 1.0]
        dust.parameter['amin'] = 0.5e-9
        dust.parameter['fraction'] = 0.085
        # If more than the fraction is changed, overwrite parameter with user input
        self.dust_chooser.update_user_parameters(dust)
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('aPyM5')
        dust.parameter['fraction'] = 0.343
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('aOlM5')
        dust.parameter['fraction'] = 0.343
        new_command_line += dust.get_command_line()
        return new_command_line


class PAH0(Dust):
    """Dust class for neutral PAH grains (see Draine & Lee 2001)."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'pah_neutral.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 2240.
        # Add _pah to the keyword to consider the pah mass correctly
        self.parameter['size_keyword'] = 'plaw_pah'
        # Minimum dust grain size
        self.parameter['amin'] = 3e-10
        # Maximum dust grain size
        self.parameter['amax'] = 5e-9
        # For creation of a dust catalog
        self.parameter['input_file'] = 'PAH0_DL07'
        self.parameter['input_type'] = 'dustem'


class PAH1(Dust):
    """Dust class for ionized PAH grains (see Draine & Lee 2001)."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'pah_ion.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 2240.
        # Add _pah to the keyword to consider the pah mass correctly
        self.parameter['size_keyword'] = 'plaw_pah'
        # Minimum dust grain size
        self.parameter['amin'] = 3e-10
        # Maximum dust grain size
        self.parameter['amax'] = 5e-9
        # For creation of a dust catalog
        self.parameter['input_file'] = 'PAH1_DL07'
        self.parameter['input_type'] = 'dustem'
