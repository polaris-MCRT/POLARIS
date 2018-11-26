#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from modules.base import ExternalInput

"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_ext_input_dict(dictionary):
    ext_input_dict = {
        'custom': CustomInput,
    }
    dictionary.update(ext_input_dict)


class CustomInput(ExternalInput):
    """Change this to the external input you want to use.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the external input parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        ExternalInput.__init__(self, file_io, parse_args)

        # Index of the MHD simulation
        self.model_index = 320

    def init_data(self):
        """Reads the data of e.g. an MHD simulation
        """
        # Save the data in this list
        self.data = [1.0]

    def gas_density_distribution(self):
        """Define here your routine to get the density at a given position
        in the model space from your read input data.

        Notes:
            Define also the following routines if necessary:
                - dust_density_distribution(self)
                - gas_temperature(self)
                - dust_temperature(self)
                - velocity_field(self)
                - magnetic_field(self)
                - dust_id(self)

        Returns:
            float: Gas density at a given position.
        """
        return self.data[0]
