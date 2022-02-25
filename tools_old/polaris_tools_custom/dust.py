#!/usr/bin/env python
# -*- coding: utf-8 -*-

from polaris_tools_modules.base import Dust


"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_dust_dict(dictionary):
    dust_dict = {
        'trust_silicate': TrustSilicate,
        'trust_graphite': TrustGraphite,
        'trust_pah': TrustPah,
        'custom': CustomDust,
    }
    dictionary.update(dust_dict)


class CustomDust(Dust):
    """Change this to the dust component you want to use.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        # Define the dust catalog file in the POLARIS standard file format
        # (relative to the polaris/input/ directory)
        self.parameter['dust_cat_file'] = 'custom.dat'
        # Relative fraction of this dust composition to mix multiple dust compositions
        self.parameter['fraction'] = 1.0
        # Material density of the custom composition [kg/m^3]
        self.parameter['material_density'] = 2500
        # Minimum dust grain size
        self.parameter['amin'] = 5e-9
        # Maximum dust grain size
        self.parameter['amax'] = 250e-9
        # Possible dust size distributions 'plaw', 'plaw-ed', 'logn'
        self.parameter['size_keyword'] = 'plaw'
        # List of size parameter for dust size distribution (plaw -> [exponent])
        self.parameter['size_parameter'] = [-3.5]

    def get_command(self):
        """Provides dust component command line for POLARIS .cmd file.

        Note:
            This demonstrates how to mix multiple dust components together and use different
                dust grain compositions throughout the model.

        Returns:
            str: Command line to consider the custom dust component.
        """
        # This shows how to mix multiple dust components and use them as one
        new_command_line = str()
        dust = self.dust_chooser.get_module_from_name('silicate_oblate')
        dust.parameter['fraction'] = 0.625
        dust.parameter['choice_id'] = 1
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('graphite_oblate')
        dust.parameter['fraction'] = 0.375
        dust.parameter['choice_id'] = 1
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('silicate_oblate')
        dust.parameter['fraction'] = 1.
        dust.parameter['choice_id'] = 2
        new_command_line += dust.get_command_line()
        return new_command_line


class TrustSilicate(Dust):
    """Dust class for silicate grains from Trust benchmark."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'trust_silicate.dat'
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 3500
        self.parameter['size_keyword'] = 'zda'
        # c0, b0, m1, b1, a1, m2, b2, a2, m3, b3, a3, m4, b4, a4,
        self.parameter['size_parameter'] = [-8.47091, -3.68708, 22.5489, 2.37316e-5, 7.64943e-3,
                                            0, 0, 0, 12.1717, 2961.28, 0.480229, 0, 0, 0]
        # Minimum dust grain size
        self.parameter['amin'] = 0.00035e-6
        # Maximum dust grain size
        self.parameter['amax'] = 0.37e-6
        # For creation of a dust catalog
        self.parameter['input_file'] = 'suvSil_121_1201'
        self.parameter['input_type'] = 'trust'
        self.parameter['subl_temp'] = 1200


class TrustGraphite(Dust):
    """Dust class for graphite grains from Trust benchmark."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'trust_graphite.dat'
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 2240
        self.parameter['size_keyword'] = 'zda'
        # c0, b0, m1, b1, a1, m2, b2, a2, m3, b3, a3, m4, b4, a4,
        self.parameter['size_parameter'] = [-9.86, -5.02082, 4.63229, 5.81215e-3, 0.415861, 0, 0, 0,
                                            3.69897, 1125.02, 0.160344, 3.69967, 1126.02, 0.160501]
        # Minimum dust grain size
        self.parameter['amin'] = 0.00035e-6
        # Maximum dust grain size
        self.parameter['amax'] = 0.33e-6
        # For creation of a dust catalog
        self.parameter['input_file'] = 'Gra_121_1201'
        self.parameter['input_type'] = 'trust'
        self.parameter['subl_temp'] = 2100


class TrustPah(Dust):
    """Dust class for PAH grains from Trust benchmark."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'trust_pah.dat'
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 2240
        self.parameter['size_keyword'] = 'zda'
        # c0, b0, m1, b1, a1, m2, b2, a2, m3, b3, a3, m4, b4, a4,
        self.parameter['size_parameter'] = [-8.02895, -3.45764, -8.20551, 1.18396e3, 1.0, 0, 0, 0,
                                            12.0146, 1.0e24, -5.29496e-3, 0, 0, 0]
        # Minimum dust grain size
        self.parameter['amin'] = 0.00035e-6
        # Maximum dust grain size
        self.parameter['amax'] = 0.005e-6
        # For creation of a dust catalog
        self.parameter['input_file'] = 'PAH_28_1201_neu'
        self.parameter['input_type'] = 'trust'
        self.parameter['subl_temp'] = 2100


