#!/usr/bin/env python
# -*- coding: utf-8 -*-

from modules.math import Math
from modules.base import Gas

"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_gas_dict(dictionary):
    gas_dict = {
        'custom': GasCustom,
    }
    dictionary.update(gas_dict)


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
