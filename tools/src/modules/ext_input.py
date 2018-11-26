#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from modules.base import ExternalInput
from custom.ext_input import *


class ExtChooser:
    """The ExtChooser class provides external input data.
    """

    def __init__(self, model, file_io, parse_args):
        """Initialisation of all usable options.

        Notes:
            To create your own external input, add its name to the dictionary
            and write a class with its options as a derived class of class ExternalInput.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
            file_io : Handles file input/output and all
                necessary paths.
            parse_args (ArgumentParser) : Provides all parameters chosen
                by user when executing PolarisTools.
        """
        self.model = model
        self.file_io = file_io
        self.parse_args = parse_args

        # dict: Dictionary with all usable external inputs
        self.ext_input_dict = {
        }
        update_ext_input_dict(self.ext_input_dict)

    def get_module(self):
        """Chooses model class from user input

            Note:
                Parameters set by PolarisTools overwrite preset values in the
                separate model classes.

            Returns:
                Instance of chosen model.
        """
        if self.parse_args.external_input_name is not None and \
                self.parse_args.external_input_name in self.ext_input_dict.keys():
            ext_name = self.parse_args.external_input_name
        elif self.model.parameter['external_input_name'] is not None:
            ext_name = self.model.parameter['external_input_name']
        else:
            return None
        ext_input = self.ext_input_dict[ext_name](
            self.file_io, self.parse_args)
        ext_input.init_data()
        return ext_input
