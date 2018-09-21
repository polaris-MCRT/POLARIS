#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from modules.base import Detector


"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""
def update_detector_dict(dictionary):
    detector_dict = {
        'gg_tau': GGTauDetector,
        'hd97048': HD97048Detector,
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
        # Number of logaritmically distributed wavelengths 
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
        self.wavelengths = None
        if parse_args.wavelength is None: 
            #self.wavelengths = np.array([3.5, 4.5, 5.6, 7.7, 10.5, 15.5]) * 1e-6 
            self.wavelengths = np.array([10., 450., 1300.]) * 1e-6
        self.parameter['rot_angle_1'] = -37.0

    def get_dust_emission_command(self):
        """Provides detector configuration command line for raytrace
        simulations for POLARIS .cmd file.

        Returns:
            str: Command line to consider the detector configuration.
        """
        if self.wavelengths is not None:
            new_command_line = str()
            for wl in self.wavelengths:
                self.parameter['wavelength_min'] = wl
                self.parameter['wavelength_max'] = wl
                new_command_line += self.get_dust_emission_command_line()
            return new_command_line 
        return self.get_dust_emission_command_line()

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
        self.wavelengths = None
        if parse_args.wavelength is None:  
            self.wavelengths = [1.25e-6, 8.6e-6, 17.8e-6]            
        # Rotation angle around the first rotation axis
        self.parameter['rot_angle_1'] = 0.
        # Rotation angle around the second rotation axis
        self.parameter['rot_angle_2'] = 43.
        # Number of pixel per axis of the detector
        self.parameter['nr_pixel_x'] = 201
        self.parameter['nr_pixel_y'] = 201

    def get_dust_emission_command(self):
        """Provides detector configuration command line for raytrace
        simulations for POLARIS .cmd file.

        Returns:
            str: Command line to consider the detector configuration.
        """
        if self.wavelengths is not None:
            new_command_line = str()
            for wl in self.wavelengths:
                self.parameter['wavelength_min'] = wl
                self.parameter['wavelength_max'] = wl
                new_command_line += self.get_dust_emission_command_line()
            return new_command_line 
        return self.get_dust_emission_command_line()