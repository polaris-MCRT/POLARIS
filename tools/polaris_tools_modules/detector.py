#!/usr/bin/env python
# -*- coding: utf-8 -*-

from polaris_tools_custom.detector import *
from polaris_tools_modules.math import Math


class DetectorChooser:
    """The DetectorChooser class provides the chosen detector configuration.
    """

    def __init__(self, model, parse_args):
        """Initialisation of all usable options.

        Notes:
            To create your own detector configuration, add its name to the dictionary
            and write a class with its options as a derived class of class Detector.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
        """

        self.model = model
        self.parse_args = parse_args

        # Get math module
        self.math = Math()

        #: dict: Dictionary with all usable detector configurations
        self.detector_dict = {
            'cartesian': Detector,
            'polar': PolarDetector,
            'allsky': AllSkyMap,
        }
        update_detector_dict(self.detector_dict)

    def get_module(self):
        """Chooses detector class from user input

        Note:
            Parameters set by PolarisTools overwrite preset values in the
            separate detector classes.

        Returns:
            Instance of chosen detector configuration.
        """

        if self.parse_args.detector in self.detector_dict.keys():
            detector_name = self.parse_args.detector
        elif self.model.parameter['detector'] is not None:
            detector_name = self.model.parameter['detector']
        else:
            raise ValueError(
                'Detector name not known! You can add a new detector in detector.py')
        detector = self.detector_dict[detector_name](
            self.model, self.parse_args)
        # Overwrite default values with user input
        if self.parse_args.wavelength is not None:
            if len(self.parse_args.wavelength) == 1:
                detector.parameter['wavelength_min'] = float(
                    self.math.parse(self.parse_args.wavelength[0], 'length'))
                detector.parameter['wavelength_max'] = float(
                    self.math.parse(self.parse_args.wavelength[0], 'length'))
                detector.parameter['nr_of_wavelength'] = 1
            elif len(self.parse_args.wavelength) == 2:
                detector.parameter['wavelength_min'] = float(
                    self.math.parse(self.parse_args.wavelength[0], 'length'))
                detector.parameter['wavelength_max'] = float(
                    self.math.parse(self.parse_args.wavelength[1], 'length'))
                detector.parameter['nr_of_wavelength'] = 2
            elif len(self.parse_args.wavelength) == 3:
                detector.parameter['wavelength_min'] = float(
                    self.math.parse(self.parse_args.wavelength[0], 'length'))
                detector.parameter['wavelength_max'] = float(
                    self.math.parse(self.parse_args.wavelength[1], 'length'))
                detector.parameter['nr_of_wavelength'] = int(
                    self.parse_args.wavelength[2])
            else:
                raise ValueError('Too many wavelengths IDs!')
        if self.parse_args.rot_angle_1 is not None:
            detector.parameter['rot_angle_1'] = self.math.parse(
                self.parse_args.rot_angle_1, 'angle')
        if self.parse_args.rot_angle_2 is not None:
            detector.parameter['rot_angle_2'] = self.math.parse(
                self.parse_args.rot_angle_2, 'angle')

        if self.parse_args.obs_pos is not None:
            detector.parameter['obs_position_x'] = self.math.parse(
                self.parse_args.obs_pos[0], 'length')
            detector.parameter['obs_position_y'] = self.math.parse(
                self.parse_args.obs_pos[1], 'length')
            detector.parameter['obs_position_z'] = self.math.parse(
                self.parse_args.obs_pos[2], 'length')
        if self.parse_args.obs_vel is not None:
            detector.parameter['obs_velocity_x'] = self.math.parse(
                self.parse_args.obs_vel[0], 'velocity')
            detector.parameter['obs_velocity_y'] = self.math.parse(
                self.parse_args.obs_vel[1], 'velocity')
            detector.parameter['obs_velocity_z'] = self.math.parse(
                self.parse_args.obs_vel[2], 'velocity')
        if self.parse_args.nr_pixel is not None:
            detector.parameter['nr_pixel_x'] = self.parse_args.nr_pixel
            detector.parameter['nr_pixel_y'] = self.parse_args.nr_pixel
        else:
            if self.parse_args.nr_pixel_x is not None:
                detector.parameter['nr_pixel_x'] = self.parse_args.nr_pixel_x
            if self.parse_args.nr_pixel_y is not None:
                detector.parameter['nr_pixel_y'] = self.parse_args.nr_pixel_y
        if self.parse_args.nr_sides is not None:
            detector.parameter['nr_sides'] = self.parse_args.nr_sides
        if self.parse_args.transition_id is not None:
            detector.parameter['transition_id'] = self.parse_args.transition_id
        if self.parse_args.max_velocity is not None:
            detector.parameter['max_velocity'] = self.math.parse(
                self.parse_args.max_velocity, 'velocity')
        if self.parse_args.nr_velocity_channels is not None:
            detector.parameter['nr_velocity_channels'] = self.parse_args.nr_velocity_channels
        if self.parse_args.sidelength_zoom is not None:
            detector.parameter['sidelength_zoom_x'] = self.parse_args.sidelength_zoom
            detector.parameter['sidelength_zoom_y'] = self.parse_args.sidelength_zoom
        if self.parse_args.sidelength_zoom_x is not None:
            detector.parameter['sidelength_zoom_x'] = self.parse_args.sidelength_zoom_x
        if self.parse_args.sidelength_zoom_y is not None:
            detector.parameter['sidelength_zoom_y'] = self.parse_args.sidelength_zoom_y
        if self.parse_args.map_shift_x is not None:
            detector.parameter['map_shift_x'] = self.math.parse(
                self.parse_args.map_shift_x, 'length')
        if self.parse_args.map_shift_y is not None:
            detector.parameter['map_shift_y'] = self.math.parse(
                self.parse_args.map_shift_y, 'length')
        if self.parse_args.acceptance_angle is not None:
            detector.parameter['acceptance_angle'] = self.parse_args.acceptance_angle
        if self.parse_args.raytracing_shape is not None:
            detector.parameter['shape'] = self.parse_args.raytracing_shape
        # elif self.model.parameter['grid_type'] in ['spherical', 'cylindrical'] and \
        #        'mc' not in self.parse_args.simulation_type:
        #    print('--- HINT: Raytracing simulation of spherical or cylindrical grid detected:\n'
        #        '    --- Use polar raytracing grid by default!\n'
        #        '    --- If this is causing problems, use \"--rt_grid cartesian\" instead!')
        #    detector.parameter['shape'] = 'polar'
        return detector


class PolarDetector(Detector):
    """This is the detector class for simulations with a polar background grid.
    """

    def __init__(self, model, parse_args):
        """Initialisation of the detector configuration.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
        """
        Detector.__init__(self, model, parse_args)

        # Set detector shape to polar
        self.parameter['shape'] = 'polar'


class AllSkyMap(Detector):
    """This is the detector class for simulations with the observer inside of the model.
    """

    def __init__(self, model, parse_args):
        """Initialisation of the detector configuration.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
        """
        Detector.__init__(self, model, parse_args)

        # Set detector shape to polar
        self.parameter['shape'] = 'healpix'
        #self.parameter['obs_position_z'] = self.math.const['au'] * 10.
