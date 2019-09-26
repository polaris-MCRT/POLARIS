#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from polaris_tools_modules.base import Detector


"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_detector_dict(dictionary):
    detector_dict = {
        'gg_tau': GGTauDetector,
        'hd97048': HD97048Detector,
        'hd169142': HD169142Detector,
        'hd100546': HD100546Detector,
        'custom': CustomDetector,
    }
    dictionary.update(detector_dict)


class CustomDetector(Detector):
    """Change this to the detector you want to use.
    """

    def __init__(self, model, parse_args):
        """Initialisation of the detector configuration.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
        """
        Detector.__init__(self, model, parse_args)
        # First wavelength of the observing wavelengths
        self.parameter['wavelength_min'] = 1e-6
        # Last wavelength of the observing wavelengths
        self.parameter['wavelength_max'] = 1e-6
        # Number of logarithmically distributed wavelengths
        # between wavelength_min and wavelength_max
        self.parameter['nr_of_wavelength'] = 1
        # Rotation angle around the first rotation axis
        self.parameter['rot_angle_1'] = 0.
        # Rotation angle around the second rotation axis
        self.parameter['rot_angle_2'] = 0.
        # Number of pixel per axis of the detector
        self.parameter['nr_pixel_x'] = 256
        # Number of pixel per axis of the detector
        self.parameter['nr_pixel_y'] = 256
        # Index of the related background source (for dust/line simulations)
        self.parameter['source_id'] = 1
        # Index of the related gas_species (for line simulations)
        self.parameter['gas_species_id'] = 1
        # Index of the related transition (for line simulations)
        self.parameter['transition_id'] = 1
        # Factor to zoom onto the observing object
        self.parameter['sidelength_zoom_x'] = 1
        self.parameter['sidelength_zoom_y'] = 1
        # Offset of the detector in x- and y-direction
        self.parameter['map_shift_x'] = 0.
        self.parameter['map_shift_y'] = 0.
        # Acceptance angle for dust_mc simulations without peel-off technique
        self.parameter['acceptance_angle'] = 1.0
        # Detector shape (defining the background grid for raytracing simulations)
        self.parameter['shape'] = 'cartesian'

    def get_dust_scattering_command(self):
        """Provides detector configuration command line for Monte-Carlo
        simulations for POLARIS .cmd file.

        Returns:
            str: Command line to consider the detector configuration.
        """
        '''To add multiple detectors, use the following:
        new_command_line = str()
        self.parameter['rot_angle_1'] = 0.0
        new_command_line += self.get_dust_scattering_command_line()
        self.parameter['rot_angle_1'] = 90.0
        new_command_line += self.get_dust_scattering_command_line()
        return new_command_line
        '''
        return self.get_dust_scattering_command_line()

    def get_dust_emission_command(self):
        """Provides detector configuration command line for raytrace
        simulations for POLARIS .cmd file.

        Returns:
            str: Command line to consider the detector configuration.
        """
        '''To add multiple detectors, use the following:
        new_command_line = str()
        self.parameter['rot_angle_1'] = 0.0
        new_command_line += self.get_dust_emission_command_line()
        self.parameter['rot_angle_1'] = 90.0
        new_command_line += self.get_dust_emission_command_line()
        return new_command_line
        '''
        return self.get_dust_emission_command_line()

    def get_line_command(self):
        """Provides detector configuration command line for spectral line
        simulations for POLARIS .cmd file.

        Returns:
            str: Command line to consider the detector configuration.
        """
        '''To add multiple detectors, use the following:
        new_command_line = str()
        self.parameter['rot_angle_1'] = 0.0
        new_command_line += self.get_line_command_line()
        self.parameter['rot_angle_1'] = 90.0
        new_command_line += self.get_line_command_line()
        return new_command_line
        '''
        return self.get_line_command_line()


class GGTauDetector(Detector):
    """Change this to the detector you want to use.
    """

    def __init__(self, model, parse_args):
        """Initialisation of the detector configuration.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
        """
        Detector.__init__(self, model, parse_args)
        # ------ Rotation angles -----
        # Cite: inclination (Duchêne et al. 2004)
        self.parameter['rot_angle_1'] = -37.0
        # Cite: position angle of Aa and Ab (Yang et al. 2017)
        self.parameter['rot_angle_2'] = 25.
        # Cite: rotation axis of inclination (Yang et al. 2017)
        inc_PA = (360. - 270. - 7.) / 180. * np.pi  # (360. - 270. - 7.)
        if self.parameter['shape'] == 'slice' or \
                self.parse_args.raytracing_shape == 'slice':
            inc_PA += (39.) / 180. * np.pi  # +39
        self.parameter['rot_axis_1'] = [np.sin(inc_PA), np.cos(inc_PA), 0]
        self.parameter['rot_axis_2'] = [0, 0, 1]
        self.parameter['nr_pixel_x'] = 512
        self.parameter['nr_pixel_y'] = 512
        self.parameter['wavelength_list'] = np.array(
            [1.6, 2.2, 3.5, 4.0, 5.0, 7.7, 10.0, 11.3, 12.0, 15.5, 23.0]) * 1e-6


class HD97048Detector(Detector):
    """Change this to the detector you want to use.
    """

    def __init__(self, model, parse_args):
        """Initialisation of the detector configuration.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
        """
        Detector.__init__(self, model, parse_args)
        # Rotation angle around the first rotation axis
        self.parameter['rot_angle_1'] = 0.
        # Rotation angle around the second rotation axis
        self.parameter['rot_angle_2'] = 43.
        # Number of pixel per axis of the detector
        self.parameter['nr_pixel_x'] = 201
        self.parameter['nr_pixel_y'] = 201
        # Wavelengths
        self.parameter['wavelength_list'] = np.array([1.25, 8.6, 17.8]) * 1e-6


class HD169142Detector(Detector):
    """Change this to the detector you want to use.
    """

    def __init__(self, model, parse_args):
        """Initialisation of the detector configuration.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
        """
        Detector.__init__(self, model, parse_args)
        # Rotation angle around the first rotation axis
        self.parameter['rot_angle_1'] = 0.
        # Rotation angle around the second rotation axis
        self.parameter['rot_angle_2'] = 13.
        # Number of pixel per axis of the detector
        self.parameter['nr_pixel_x'] = 201
        self.parameter['nr_pixel_y'] = 201
        # Wavelengths
        self.parameter['wavelength_list'] = np.array([1.25, 8.6, 17.8]) * 1e-6


class HD100546Detector(Detector):
    """Detector for the radiation field disk for Thomas.
    """

    def __init__(self, model, parse_args):
        """Initialisation of the detector configuration.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
        """
        Detector.__init__(self, model, parse_args)
        # Set the background RT grid
        if self.parse_args.raytracing_shape is None:
            self.parameter['shape'] = 'polar'
        # Set the wavelengths
        self.parameter['wavelength_list'] = np.array(
            [3, 3.28, 3.38, 3.42, 3.7]) * 1e-6
        # Rotation angle around the first rotation axis
        self.parameter['rot_angle_1'] = 0.
        # Rotation angle around the second rotation axis
        self.parameter['rot_angle_2'] = 42.
