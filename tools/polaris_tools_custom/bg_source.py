#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from polaris_tools_modules.base import BGSource


"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_bg_source_dict(dictionary):
    bg_source_dict = {
        'custom': CustomBGsource,
    }
    dictionary.update(bg_source_dict)


class CustomBGsource(BGSource):
    """Change this to the background source you want to use.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the background source parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """

        BGSource.__init__(self, file_io, parse_args)

        #: dict: Parameters which are different to the default values
        bg_source_parameter = {
            'A': 1.0e-3,
            'T': 1000.,
            'Q': 1.,
            'U': 0.,
            'V': 0.,
        }

        # Updates the parameter dictionary
        self.parameter.update(bg_source_parameter)
