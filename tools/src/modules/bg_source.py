#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from modules.base import BGSource
from custom.bg_source import *


class BGSourceChooser:
    """The BGSourceChooser class provides the chosen background source
    such as background stars.
    """

    def __init__(self, model, file_io, parse_args):
        """Initialisation of all usable options.

        Notes:
            To create your own background source, add its name to the dictionary
            and write a class with its options as a derived class of class BGSource.

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

        #: dict: Dictionary with all usable background source types
        self.bg_source_dict = {
            'bg_plane': BackgroundPlane,
        }
        update_bg_source_dict(self.bg_source_dict)

    def get_module(self):
        """Chooses background source from user input

        Note:
            Parameters set by PolarisTools overwrite preset values in the
            separate stellar source classes.

        Returns:
            Instance of chosen background source.
        """

        if self.parse_args.bg_source in self.bg_source_dict.keys():
            bg_source_name = self.parse_args.bg_source
        elif self.model.parameter['background_source'] is not None:
            bg_source_name = self.model.parameter['background_source']
        else:
            raise ValueError('Background source not known! '
                             'You can add a new background source in bg_source.py')
        bg_source = self.bg_source_dict[bg_source_name](
            self.file_io, self.parse_args)
        # Overwrite default values with user input
        if self.parse_args.nr_photons is not None:
            bg_source.parameter['nr_photons'] = int(self.parse_args.nr_photons)
        return bg_source


class BackgroundPlane(BGSource):
    """The BackgroundPlane class is a radiating plane behind the model space.
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
        }

        # Updates the parameter dictionary
        self.parameter.update(bg_source_parameter)
