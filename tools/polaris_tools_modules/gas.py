#!/usr/bin/env python
# -*- coding: utf-8 -*-

from polaris_tools_modules.math import Math
from polaris_tools_modules.base import Gas, NoGas
from polaris_tools_custom.gas import *


class GasChooser:
    """The GasChooser class provides the chosen gas species.
    """

    def __init__(self, file_io, parse_args, model=None):
        """Initialisation of all usable options.

        Notes:
            To create your own gas species, add its name to the dictionary
            and write a class with its options as a derived class of class Gas.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
            file_io : Handles file input/output and all
                necessary paths.
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
        """
        self.file_io = file_io
        self.parse_args = parse_args
        self.model = model

        # Get math module
        self.math = Math()

        #: dict[str: int]: Index of different level population  methods
        self.lvl_pop = {
            'LTE': 1,
            'FEP': 2,
            'LVG': 3,
        }

        #: dict: Dictionary with all usable gas species
        self.gas_dict = {
            'none': NoGas,
            'CO': GasCO,
            '13CO': Gas13CO,
            'C18O': GasC18O,
            'HCO': GasHCO,
            'SO': GasSO,
            'OH': GasOH,
            'CN': GasCN,
            'H1': GasH1,
        }
        update_gas_dict(self.gas_dict)

    def get_module(self):
        """Chooses gas class from user input

        Note:
            Parameters set by PolarisTools overwrite preset values in the
            separate gas classes.

        Returns:
            Instance of chosen gas species.
        """
        if any(simu_type in self.parse_args.simulation_type for simu_type in ['temp', 'rat', 'dust_mc', 'dust']):
            gas_species_name = 'none'
        elif self.parse_args.gas_species is not None and \
                self.parse_args.gas_species.upper() in self.gas_dict.keys():
            gas_species_name = self.parse_args.gas_species.upper()
        elif self.model is not None and self.model.parameter['gas_species'] is not None:
            gas_species_name = self.model.parameter['gas_species'].upper()
        else:
            raise ValueError(
                'Gas species not known! You can add a new species in gas.py!')
        gas = self.gas_dict[gas_species_name](self.file_io, self.parse_args)
        # Overwrite default values with user input
        if self.parse_args.level_pop_type is not None:
            gas.parameter['level_pop_type'] = self.lvl_pop[self.parse_args.level_pop_type]
        if self.parse_args.abundance is not None:
            gas.parameter['abundance'] = self.parse_args.abundance

        return gas

    def get_module_from_name(self, species_name):
        """Chooses gas class from name string

        Args:
            species_name (str): Name of the gas species.

        Returns:
            Instance of chosen gas species.
        """
        if species_name in self.gas_dict.keys():
            gas = self.gas_dict[species_name](self.file_io, self.parse_args)
        else:
            raise ValueError('No gas species with the name ' +
                             str(species_name) + '!')
        return gas


class GasCustom(Gas):
    """Change this to the gas species you want to use.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Gas.__init__(self, file_io, parse_args)

        # Set parameters of the custom gas species
        # (see parent Gas class for available options)
        self.parameter['filename'] = 'gas_database_file.dat'
        self.parameter['abundance'] = 1e-4
        self.parameter['zeeman_usable'] = False

    def shift_2_mag(self, velocity_shift, f_0, i_trans):
        """Calculates the magnetic field strength related to a Zeeman shift of
            a given spectral line transition.

        Args:
            velocity_shift (float) : Zeeman shift in [m/s]
            f_0 (float) : Rest frequency of spectral line [Hz]
            i_trans (int) : transition index of spectral line

        Returns:
            float: Magnetic field strength in [T]
        """
        return 1.0

    def mag_2_shift(self, magnetic_field, f_0, i_trans):
        """Calculates the Zeeman frequency shift of a spectral line
        transition related to a given magnetic field strength.

        Args:
            magnetic_field (float) : Magnetic field strength [T]
            f_0 (float) : Rest frequency of spectral line [Hz]
            i_trans (int) : transition index of spectral line

        Returns:
            float: Zeeman shift in [m/s]
        """
        return 1.0


class GasCO(Gas):
    """Gas class for carbon monoxide.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Gas.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        gas_species_co_parameter = {
            'filename': 'co.dat',
            'abundance': 1e-4,
            'nr_velocity_channels': 35,
        }

        # Updates the parameter dictionary
        self.parameter.update(gas_species_co_parameter)


class Gas13CO(Gas):
    """Gas class for carbon monoxide (13CO isotope).
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Gas.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        gas_species_13co_parameter = {
            'filename': '13co.dat',
        }

        # Updates the parameter dictionary
        self.parameter.update(gas_species_13co_parameter)


class GasC18O(Gas):
    """Gas class for carbon monoxide (C18O isotope).
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Gas.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        gas_species_c18o_parameter = {
            'filename': 'c18o.dat',
        }

        # Updates the parameter dictionary
        self.parameter.update(gas_species_c18o_parameter)


class GasHCO(Gas):
    """Gas class for HCO+.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Gas.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        gas_species_hco_parameter = {
            'filename': 'hco+@xpol.dat',
        }

        # Updates the parameter dictionary
        self.parameter.update(gas_species_hco_parameter)


class GasSO(Gas):
    """Gas class for sulfur monoxide.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Gas.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        gas_species_so_parameter = {
            'filename': 'so.dat',
            'abundance': 1e-10,  # Pacheco-Vázquez et al. 2016
            'nr_velocity_channels': 35,
            'zeeman_usable': True,
        }

        #: dict: Zeeman shifts for each supported transition 2*nu/B [Hz/muG]
        self.zeeman_shift = {
            3: 2. * 0.51848,  # 1.0
            8: 1.7,
            10: 0.8,
            28: 1.0,
            34: 0.5,
        }

        # Updates the parameter dictionary
        self.parameter.update(gas_species_so_parameter)

    def shift_2_mag(self, velocity_shift, f_0, i_trans):
        """Calculates the magnetic field strength related to a Zeeman shift of
            a given spectral line transition.

        Args:
            velocity_shift (float) : Zeeman shift in [m/s]
            f_0 (float) : Rest frequency of spectral line [Hz]
            i_trans (int) : transition index of spectral line

        Returns:
            float: Magnetic field strength in [T]
        """
        try:
            frequency_shift = self.math.velocity_2_frequency(
                velocity_shift, f_0)
            magnetic_field = 1e-10 * frequency_shift / \
                (0.5 * self.zeeman_shift[i_trans])
        except ValueError:
            print('The chosen species has a known Zeeman splitting in some of its transitions. '
                  'However, you have not chosen one of them!')
        return magnetic_field

    def mag_2_shift(self, magnetic_field, f_0, i_trans):
        """Calculates the Zeeman frequency shift of a spectral line
        transition related to a given magnetic field strength.

        Args:
            magnetic_field (float) : Magnetic field strength [T]
            f_0 (float) : Rest frequency of spectral line [Hz]
            i_trans (int) : transition index of spectral line

        Returns:
            float: Zeeman shift in [m/s]
        """
        try:
            frequency_shift = (magnetic_field / 1e-10) * \
                (0.5 * self.zeeman_shift[i_trans])
            velocity_shift = self.math.frequency_2_velocity(
                frequency_shift, f_0)
        except ValueError:
            print('The chosen species has a known Zeeman splitting in some of its transitions. '
                  'However, you have not chosen one of them!')
        return velocity_shift


class GasOH(Gas):
    """Gas class for hydroxyl radical.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Gas.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        gas_species_oh_parameter = {
            'filename': 'oh@hfs.dat',
            'abundance': 1e-7,
            'nr_velocity_channels': 35,
            'zeeman_usable': True,
            'lande_factor': [0, 0, 1.16828926744, 0.700973560462],
        }
        # Lande factors for Zeeman splitting of OH
        # (zero level: not existent!, 1st level: not available for Zeeman!)

        # Updates the parameter dictionary
        self.parameter.update(gas_species_oh_parameter)

    def shift_2_mag(self, velocity_shift, f_0, i_trans):
        """Calculates the magnetic field strength related to a Zeeman shift of
            a given spectral line transition.

        Args:
            velocity_shift (float) : Zeeman shift in [m/s]
            f_0 (float) : Rest frequency of spectral line [Hz]
            i_trans (int) : transition index of spectral line

        Returns:
            float: Magnetic field strength in [T]
        """
        try:
            frequency_shift = self.math.velocity_2_frequency(
                velocity_shift, f_0)
            constant = self.math.const['h'] / self.math.const['mu_B']
            magnetic_field = frequency_shift * constant / \
                float(self.parameter['lande_factor'][i_trans])
        except ValueError:
            print('The chosen species has a known Zeeman splitting in some of its transitions. '
                  'However, you have not chosen one of them!')
        return magnetic_field

    def mag_2_shift(self, magnetic_field, f_0, i_trans):
        """Calculates the Zeeman frequency shift of a spectral line
        transition related to a given magnetic field strength.

        Args:
            magnetic_field (float) : Magnetic field strength [T]
            f_0 (float) : Rest frequency of spectral line [Hz]
            i_trans (int) : transition index of spectral line

        Returns:
            float: Zeeman shift in [m/s]
        """
        try:
            constant = self.math.const['mu_B'] / self.math.const['h']
            frequency_shift = magnetic_field * constant * \
                float(self.parameter['lande_factor'][i_trans])
            velocity_shift = self.math.frequency_2_velocity(
                frequency_shift, f_0)
        except ValueError:
            print('The chosen species has a known Zeeman splitting in some of its transitions. '
                  'However, you have not chosen one of them!')
        return velocity_shift


class GasCN(Gas):
    """Gas class for cyanide radical.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Gas.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        gas_species_so_parameter = {
            'filename': 'cn.dat',
            'abundance': 4e-9,  # Falgarone et al. 2008
            'nr_velocity_channels': 35,
            'zeeman_usable': True,
        }

        #: dict: Zeeman shifts for each supported transition 2*nu/B [Hz/muG]
        self.zeeman_shift = {
            11: 2.17719367,
            12: -0.31102767,
            13: 0.62205534,
            14: 2.17719367,
            15: 0.5598498,
            16: 0.62205534,
            17: 1.61734387
        }

        # Updates the parameter dictionary
        self.parameter.update(gas_species_so_parameter)

    def shift_2_mag(self, velocity_shift, f_0, i_trans):
        """Calculates the magnetic field strength related to a Zeeman shift of
            a given spectral line transition.

        Args:
            velocity_shift (float) : Zeeman shift in [m/s]
            f_0 (float) : Rest frequency of spectral line [Hz]
            i_trans (int) : transition index of spectral line

        Returns:
            float: Magnetic field strength in [T]
        """
        try:
            frequency_shift = self.math.velocity_2_frequency(
                velocity_shift, f_0)
            magnetic_field = frequency_shift / \
                (0.5 * self.zeeman_shift[i_trans]) * 1e-10
        except ValueError:
            print('The chosen species has a known Zeeman splitting in some of its transitions. '
                  'However, you have not chosen one of them!')
        return magnetic_field

    def mag_2_shift(self, magnetic_field, f_0, i_trans):
        """Calculates the Zeeman frequency shift of a spectral line
        transition related to a given magnetic field strength.

        Args:
            magnetic_field (float) : Magnetic field strength [T]
            f_0 (float) : Rest frequency of spectral line [Hz]
            i_trans (int) : transition index of spectral line

        Returns:
            float: Zeeman shift in [m/s]
        """
        try:
            frequency_shift = (magnetic_field / 1e-10) * \
                (0.5 * self.zeeman_shift[i_trans])
            velocity_shift = self.math.frequency_2_velocity(
                frequency_shift, f_0)
        except ValueError:
            print('The chosen species has a known Zeeman splitting in some of its transitions. '
                  'However, you have not chosen one of them!')
        return velocity_shift


class GasH1(Gas):
    """Gas class for atomic hydrogen.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Gas.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        gas_species_h1_parameter = {
            'filename': 'h1.dat',
            'abundance': 0.8,  # Source missing
            'nr_velocity_channels': 35,
            'zeeman_usable': True,
        }

        #: dict: Zeeman shifts for each supported transition 2*nu/B [Hz/muG]
        self.zeeman_shift = {
            1: 2.799,
        }

        # Updates the parameter dictionary
        self.parameter.update(gas_species_h1_parameter)

    def shift_2_mag(self, velocity_shift, f_0, i_trans):
        """Calculates the magnetic field strength related to a Zeeman shift of
            a given spectral line transition.

        Args:
            velocity_shift (float) : Zeeman shift in [m/s]
            f_0 (float) : Rest frequency of spectral line [Hz]
            i_trans (int) : transition index of spectral line

        Returns:
            float: Magnetic field strength in [T]
        """
        try:
            frequency_shift = self.math.velocity_2_frequency(
                velocity_shift, f_0)
            magnetic_field = frequency_shift / \
                (0.5 * self.zeeman_shift[i_trans]) * 1e-10
        except ValueError:
            print('The chosen species has a known Zeeman splitting in some of its transitions. '
                  'However, you have not chosen one of them!')
        return magnetic_field

    def mag_2_shift(self, magnetic_field, f_0, i_trans):
        """Calculates the Zeeman frequency shift of a spectral line
        transition related to a given magnetic field strength.

        Args:
            magnetic_field (float) : Magnetic field strength [T]
            f_0 (float) : Rest frequency of spectral line [Hz]
            i_trans (int) : transition index of spectral line

        Returns:
            float: Zeeman shift in [m/s]
        """
        try:
            frequency_shift = (magnetic_field / 1e-10) * \
                (0.5 * self.zeeman_shift[i_trans])
            velocity_shift = self.math.frequency_2_velocity(
                frequency_shift, f_0)
        except ValueError:
            print('The chosen species has a known Zeeman splitting in some of its transitions. '
                  'However, you have not chosen one of them!')
        return velocity_shift
