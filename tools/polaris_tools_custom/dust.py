#!/usr/bin/env python
# -*- coding: utf-8 -*-

from polaris_tools_modules.base import Dust
import numpy as np


"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_dust_dict(dictionary):
    dust_dict = {
        'vincent_mrn_oblate': VincentMrnOblate,
        'silicate_oblate_vg': SilicateOblateVG,
        'multi_sil': MultiSilicate,
        'multi_mrn': MultiMRN,
        'multi_themis': MultiThemis,
        'themis_sil': ThemisLargeSilicate,
        'olivine': Olivine,
        'silicate_pah': SilicatePAH,
        'olivine_pah': OlivinePAH,
        'olivine_themis_CM20': OlivineThemisCM20,
        'olivine_pah_carbon': OlivinePAHCarbon,
        'trust_silicate': TrustSilicate,
        'trust_graphite': TrustGraphite,
        'trust_pah': TrustPah,
        'thomas_themis': ThomasThemis,
        'thomas_CM20': ThomasCM20,
        'thomas_aPyM5': ThomasPyM5,
        'thomas_aOlM5': ThomasOlM5,
        'pebbles_typ1': PebblesTyp1,
        'pebbles_typ2': PebblesTyp2,
        'pebbles_typ3': PebblesTyp3,
        'pebbles_typ4': PebblesTyp4,
        'pebbles_typ5': PebblesTyp5,
        'pebbles_typ6': PebblesTyp6,
        'pebbles_typ7': PebblesTyp7,
        'pebbles_typ8': PebblesTyp8,
        'pebbles_typ9': PebblesTyp9,
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


class VincentMrnOblate(Dust):
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
        dust.parameter['fraction'] = 0.67
        dust.parameter['amin'] = 8.0e-9
        dust.parameter['amax'] = 4.0e-7
        dust.parameter['size_keyword'] = 'plaw'
        dust.parameter['size_parameter'] = [-3.5]
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('graphite_oblate')
        dust.parameter['fraction'] = 0.33
        dust.parameter['amin'] = 1.0e-8
        dust.parameter['amax'] = 1.7e-7
        dust.parameter['size_keyword'] = 'plaw'
        dust.parameter['size_parameter'] = [-3.9]
        new_command_line += dust.get_command_line()
        return new_command_line


class SilicateOblateVG(Dust):
    """Dust class for silicate grains."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        # Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'test_silicate_prolate_vg.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 3800
        # Minimum dust grain size
        self.parameter['amin'] = 3e-9
        # Maximum dust grain size
        self.parameter['amax'] = 2e-6


class MultiSilicate(Dust):
    """Dust class for silicate grains with different sizes in different
    regions of the model.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'silicate.nk'
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 3500

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the MRN dust composition.
        """
        amax_list = [1e-6, 2.5e-6]
        new_command_line = str()
        for i in range(len(amax_list)):
            self.parameter['amin'] = 5e-9
            self.parameter['amax'] = amax_list[i]
            self.parameter['choice_id'] = i
            new_command_line += self.get_command_line()
        return new_command_line


class MultiMRN(Dust):
    """Dust class for silicate and PAH grains.
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
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 0

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the MRN dust composition.
        """
        # amax_list = [0.25e-6, 2.5e-6]  # [0.25e-6, 0.5e-6]
        #amax_list = [2.5e-07, 0.3e-6, 0.4e-6, 0.5e-6, 2.0e-6, 3.5e-6, 300e-6]
        amax_list = [0.25e-6, 3e-6, 100e-6]
        new_command_line = str()
        for i in range(len(amax_list)):
            dust = self.dust_chooser.get_module_from_name('silicate')
            dust.parameter['fraction'] = 0.625
            dust.parameter['amin'] = 5e-9
            dust.parameter['amax'] = amax_list[i]
            dust.parameter['choice_id'] = i
            new_command_line += dust.get_command_line()
            dust = self.dust_chooser.get_module_from_name('graphite_perpend')
            dust.parameter['fraction'] = 0.25
            dust.parameter['amin'] = 5e-9
            dust.parameter['amax'] = amax_list[i]
            dust.parameter['choice_id'] = i
            new_command_line += dust.get_command_line()
            dust = self.dust_chooser.get_module_from_name('graphite_parallel')
            dust.parameter['fraction'] = 0.125
            dust.parameter['amin'] = 5e-9
            dust.parameter['amax'] = amax_list[i]
            dust.parameter['choice_id'] = i
            new_command_line += dust.get_command_line()
        return new_command_line


class MultiThemis(Dust):
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
        # List of dust choice IDs
        self.choice_id = [0, 1]

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the MRN dust composition.
        """
        new_command_line = str()
        # First dust choice is the full themis model
        dust = self.dust_chooser.get_module_from_name('CM20_1')
        dust.parameter['fraction'] = 0.229
        dust.parameter['choice_id'] = self.choice_id[0]
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('CM20_2')
        dust.parameter['fraction'] = 0.085
        dust.parameter['choice_id'] = self.choice_id[0]
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('aPyM5')
        dust.parameter['fraction'] = 0.343
        dust.parameter['choice_id'] = self.choice_id[0]
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('aOlM5')
        dust.parameter['fraction'] = 0.343
        dust.parameter['choice_id'] = self.choice_id[0]
        new_command_line += dust.get_command_line()
        # Adding only PAHs as a second dust choice
        dust = self.dust_chooser.get_module_from_name('CM20_1')
        dust.parameter['fraction'] = 1.0
        dust.parameter['choice_id'] = self.choice_id[1]
        new_command_line += dust.get_command_line()
        return new_command_line


class ThemisLargeSilicate(Dust):
    """Dust class for themis dust grains plus larger silicates.
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
        # First dust choice is the full themis model
        dust = self.dust_chooser.get_module_from_name('CM20')
        dust.parameter['choice_id'] = 0
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('CM20')
        dust.parameter['material_density'] = 1570
        dust.parameter['size_keyword'] = 'logn'
        dust.parameter['size_parameter'] = [7e-9, 1.0]
        dust.parameter['amin'] = 0.5e-9
        dust.parameter['choice_id'] = 1
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('aPyM5')
        dust.parameter['choice_id'] = 2
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('aOlM5')
        dust.parameter['choice_id'] = 3
        new_command_line += dust.get_command_line()
        # Adding larger silicates to the THEMIS model
        dust = self.dust_chooser.get_module_from_name('silicate')
        dust.parameter['amin'] = 10e-6
        dust.parameter['amax'] = 1e-3
        dust.parameter['choice_id'] = 4
        new_command_line += dust.get_command_line()
        return new_command_line


class Olivine(Dust):
    """Dust class for olivine silicate grains."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        # For creation of a dust catalog
        self.parameter['dust_cat_file'] = 'olivine.nk'
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 1500.
        # Minimum dust grain size
        self.parameter['amin'] = 0.03e-6
        # Maximum dust grain size
        self.parameter['amax'] = 1000.0e-6


class SilicatePAH(Dust):
    """Dust class for silicate and PAH grains.
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
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 0

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the MRN dust composition.
        """
        new_command_line = str()
        # First dust choice is the full themis model
        dust = self.dust_chooser.get_module_from_name('silicate')
        dust.parameter['size_keyword'] = 'plaw'
        dust.parameter['size_parameter'] = [-3.5]
        dust.parameter['amin'] = 0.03e-6
        dust.parameter['amax'] = 1000.e-6
        dust.parameter['choice_id'] = 0
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('pah_ion')
        # Add _pah to the keyword to consider the pah mass correctly
        dust.parameter['size_keyword'] = 'plaw_pah'
        dust.parameter['size_parameter'] = [-3.5]
        dust.parameter['amin'] = 4e-10
        dust.parameter['amax'] = 4e-10
        dust.parameter['choice_id'] = 1
        new_command_line += dust.get_command_line()
        return new_command_line


class OlivinePAH(Dust):
    """Dust class for olivine and PAH grains.
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
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 0

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the MRN dust composition.
        """
        new_command_line = str()
        # First dust choice is the full themis model
        dust = self.dust_chooser.get_module_from_name('olivine')
        # For HD169142
        dust.parameter['size_parameter'] = [-3.7]
        dust.parameter['amin'] = 0.04637e-6
        dust.parameter['amax'] = 1e-3  # 6506.4e-6
        # ---
        dust.parameter['choice_id'] = 0
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('pah_ion')
        # Add _pah to the keyword to consider the pah mass correctly
        dust.parameter['size_keyword'] = 'plaw_pah'
        dust.parameter['size_parameter'] = [-3.5]
        dust.parameter['amin'] = 4e-10
        dust.parameter['amax'] = 4e-10
        dust.parameter['choice_id'] = 1
        new_command_line += dust.get_command_line()
        return new_command_line


class OlivineThemisCM20(Dust):
    """Dust class for olivine and PAH grains.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Note:
            Link: https: // www.ias.u-psud.fr/themis/

        Args:
            file_io: Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 0

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the MRN dust composition.
        """
        new_command_line = str()
        # First dust choice is the full themis model
        dust = self.dust_chooser.get_module_from_name('olivine')
        # For HD169142
        dust.parameter['size_parameter'] = [-3.7]
        dust.parameter['amin'] = 0.04637e-6
        dust.parameter['amax'] = 1e-3  # 6506.4e-6
        # ---
        dust.parameter['choice_id'] = 0
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('thomas_CM20')
        # Add _pah to the keyword to consider the pah mass correctly
        dust.parameter['choice_id'] = 1
        new_command_line += dust.get_command_line()
        return new_command_line


class OlivinePAHCarbon(Dust):
    """Dust class for olivine and olivine + PAH grains.
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
        self.parameter['scattering'] = 'MIE'
        self.parameter['material_density'] = 0

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the MRN dust composition.
        """
        new_command_line = str()
        # First dust choice is the full themis model
        dust = self.dust_chooser.get_module_from_name('olivine')
        dust.parameter['size_keyword'] = 'plaw'
        dust.parameter['size_parameter'] = [-3.5]
        dust.parameter['amin'] = 0.03e-6
        dust.parameter['amax'] = 1000.0e-6
        dust.parameter['choice_id'] = 0
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('pah_ion')
        # Add _pah to the keyword to consider the pah mass correctly
        dust.parameter['size_keyword'] = 'plaw_pah'
        dust.parameter['size_parameter'] = [-3.5]
        dust.parameter['amin'] = 4e-10
        dust.parameter['amax'] = 4e-10
        dust.parameter['choice_id'] = 1
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('trust_graphite')
        dust.parameter['size_keyword'] = 'plaw'
        dust.parameter['size_parameter'] = [-3.5]
        dust.parameter['amin'] = 0.001e-6
        dust.parameter['amax'] = 0.01e-6
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


class ThomasCM20(Dust):
    """Dust class for carbon grains from Themis model."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'thomas_CM20.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 1600
        self.parameter['subl_temp'] = 2100
        self.parameter['size_keyword'] = 'plaw-ed'
        self.parameter['size_parameter'] = [-5, 10e-9, 50e-9,  1.0]
        # Minimum dust grain size (5 nm)
        self.parameter['amin'] = 0.4e-09
        # Maximum dust grain size
        self.parameter['amax'] = 4.9e-06
        # For creation of a dust catalog
        self.parameter['input_file'] = 'CM20_0.10'
        self.parameter['input_type'] = 'dustem'
        self.parameter['dustem_wl_file'] = 'LAMBDA_THOMAS'


class ThomasPyM5(Dust):
    """Dust class for silicate grains from Themis model."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'thomas_aPyM5.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 2190
        self.parameter['subl_temp'] = 1200
        self.parameter['size_keyword'] = 'logn'
        self.parameter['size_parameter'] = [8e-9, 1.0]
        # Minimum dust grain size (5 nm)
        self.parameter['amin'] = 1e-09
        # Maximum dust grain size
        self.parameter['amax'] = 4.9e-06
        # For creation of a dust catalog
        self.parameter['input_file'] = 'thomas_aPyM5'
        self.parameter['input_type'] = 'dustem'
        self.parameter['dustem_wl_file'] = 'LAMBDA_THOMAS'


class ThomasOlM5(Dust):
    """Dust class for silicate grains from Themis model."""

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'thomas_aOlM5.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 2190
        self.parameter['subl_temp'] = 1200
        self.parameter['size_keyword'] = 'logn'
        self.parameter['size_parameter'] = [8e-9, 1.0]
        # Minimum dust grain size (5 nm)
        self.parameter['amin'] = 1e-09
        # Maximum dust grain size
        self.parameter['amax'] = 4.9e-06
        # For creation of a dust catalog
        self.parameter['input_file'] = 'thomas_aOlM5'
        self.parameter['input_type'] = 'dustem'
        self.parameter['dustem_wl_file'] = 'LAMBDA_THOMAS'


class ThomasThemis(Dust):
    """Dust class for themis dust grains (thomas version 1).
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
        # Set the name of one component to allow print of sizes and wavelengths
        self.parameter['dust_cat_file'] = 'thomas_CM20.dat'

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the MRN dust composition.
        """
        new_command_line = str()
        dust = self.dust_chooser.get_module_from_name('thomas_CM20')
        dust.parameter['material_density'] = 1596
        dust.parameter['choice_id'] = 0
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('thomas_CM20')
        dust.parameter['material_density'] = 1596
        dust.parameter['size_keyword'] = 'logn'
        dust.parameter['size_parameter'] = [7e-9, 1.0]
        dust.parameter['amin'] = 0.5e-9
        dust.parameter['choice_id'] = 1
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('thomas_aPyM5')
        dust.parameter['choice_id'] = 2
        new_command_line += dust.get_command_line()
        dust = self.dust_chooser.get_module_from_name('thomas_aOlM5')
        dust.parameter['choice_id'] = 3
        new_command_line += dust.get_command_line()
        return new_command_line


class PebblesTyp1(Dust):
    """Dust class for large silicates from THEMIS model.
    See: From grain to pebble (Ysard et al.)
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        self.parameter['dust_cat_file'] = 'pebbles_typ1.dat'
        self.parameter['scattering'] = 'HG'
        self.parameter['material_density'] = 2950
        self.parameter['subl_temp'] = 1200
        self.parameter['size_keyword'] = 'logn'
        self.parameter['size_parameter'] = [0.1, 7e-1]
        # Minimum dust grain size
        self.parameter['amin'] = 1e-8
        # Maximum dust grain size
        self.parameter['amax'] = 9.99e-2
        # For creation of a dust catalog
        self.parameter['input_file'] = 'typ1'
        self.parameter['input_type'] = 'dustem'
        self.parameter['dustem_wl_file'] = 'LAMBDA_PEBBLES'


class PebblesTyp2(PebblesTyp1):
    def __init__(self, file_io, parse_args):
        PebblesTyp1.__init__(self, file_io, parse_args)
        self.parameter['dust_cat_file'] = 'pebbles_typ2.dat'
        self.parameter['input_file'] = 'typ2'
        self.parameter['material_density'] = 1600


class PebblesTyp3(PebblesTyp1):
    def __init__(self, file_io, parse_args):
        PebblesTyp1.__init__(self, file_io, parse_args)
        self.parameter['dust_cat_file'] = 'pebbles_typ3.dat'
        self.parameter['input_file'] = 'typ3'
        self.parameter['material_density'] = 1300


class PebblesTyp4(PebblesTyp1):
    def __init__(self, file_io, parse_args):
        PebblesTyp1.__init__(self, file_io, parse_args)
        self.parameter['dust_cat_file'] = 'pebbles_typ4.dat'
        self.parameter['input_file'] = 'typ4'
        self.parameter['material_density'] = 2500


class PebblesTyp5(PebblesTyp1):
    def __init__(self, file_io, parse_args):
        PebblesTyp1.__init__(self, file_io, parse_args)
        self.parameter['dust_cat_file'] = 'pebbles_typ5.dat'
        self.parameter['input_file'] = 'typ5'
        self.parameter['material_density'] = 2400


class PebblesTyp6(PebblesTyp1):
    def __init__(self, file_io, parse_args):
        PebblesTyp1.__init__(self, file_io, parse_args)
        self.parameter['dust_cat_file'] = 'pebbles_typ6.dat'
        self.parameter['input_file'] = 'typ6'
        self.parameter['material_density'] = 1200


class PebblesTyp7(PebblesTyp1):
    def __init__(self, file_io, parse_args):
        PebblesTyp1.__init__(self, file_io, parse_args)
        self.parameter['dust_cat_file'] = 'pebbles_typ7.dat'
        self.parameter['input_file'] = 'typ7'
        self.parameter['material_density'] = 1680


class PebblesTyp8(PebblesTyp1):
    def __init__(self, file_io, parse_args):
        PebblesTyp1.__init__(self, file_io, parse_args)
        self.parameter['dust_cat_file'] = 'pebbles_typ8.dat'
        self.parameter['input_file'] = 'typ8'
        self.parameter['material_density'] = 1880


class PebblesTyp9(PebblesTyp1):
    def __init__(self, file_io, parse_args):
        PebblesTyp1.__init__(self, file_io, parse_args)
        self.parameter['dust_cat_file'] = 'pebbles_typ9.dat'
        self.parameter['input_file'] = 'typ9'
        self.parameter['material_density'] = 1340
