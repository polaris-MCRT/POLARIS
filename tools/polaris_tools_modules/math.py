# -*- coding: utf-8 -*-

import numpy as np

# Import curve_fit from scipy package
# from scipy.optimize import curve_fit


class Math:
    """Constants and math functions.
    """

    def __init__(self):
        """Initialisation of all constants and conversion factors.
        """
        #: dict: Conversion factor from (Kruegel et. al 2008)
        # self.conversion_factor = {
        #     'U': 1810,
        #     'B': 4260,
        #     'V': 3640,
        #     'R': 3100,
        #     'I': 2500,
        #     'J': 1635,
        #     'H': 1090,
        #     'K': 665,
        #     'L': 277,
        #     'M': 164,
        #     'N': 37,
        #     'Q': 10,
        # }

        #: dict: Constants taken from astropy (reference: CODATA 2014, IAU 2012 Resolution B2)
        self.const = {
            'M_sun':        1.9884754153381438e+30,  # Solar mass [kg]
            'M_jup':        1.8981871658715508e+27,  # Jupiter mass [kg]
            'R_sun':        695700000.0,             # Nominal solar radius [m]
            # Nominal solar luminosity [W]
            'L_sun':        3.828e+26,
            'au':           149597870700.0,          # Astronomical Unit [m]
            'pc':           3.0856775814671916e+16,  # Parsec [m]
            'amu':          1.66053904e-27,          # Atomic mass unit [kg]
            # Gravitational constant [m^3 / (kg * s^2)]
            'G':            6.67408e-11,
            'h':            6.62607004e-34,          # Planck constant [Js]
            # Reduced Planck constant [Js]
            'hbar':         1.0545718e-34,
            # Speed of light in vacuum [m / s]
            'c':            299792458,
            'e':            1.6021766208e-19,        # Electron charge [C]
            'mu_B':         9.274009994e-24,         # Bohr magneton [J / T]
            # Boltzmann constant [J / K]
            'k_B':          1.38064852e-23,
            'Rydberg':      10973731.568508,         # Rydberg constant [1 / m]
            # Stefan-Boltzmann constant [W / (m^2 * K^4)]
            'sigma_sb':     5.670367e-08,
            'm_e':          9.10938356e-31,          # Electron mass [kg]
            'm_p':          1.672621898e-27,         # Proton mass [kg]
            # Vacuum permittivity [F / m]
            'epsilon_0':    8.854187817620389e-12,
            'avg_gas_mass': 2.,                      # Average atomic mass unit per gas particle
        }

        #: dict: Covalent radii of various atoms
        # self.covalent_radii = {
        #     'H': 31e-12,
        #     'C': 76e-12,
        #     'N': 71e-12,
        #     'O': 66e-12,
        #     'S': 105e-12,
        # }

    # @staticmethod
    # def latex_float(number):
    #     """Convert number to a string that can be used by latex.

    #     Args:
    #         number (float): Number which should be converted.

    #     Returns:
    #         str: Number converted to a string that can be used by latex (label, ...).
    #     """
    #     float_str = "{:1.2f}".format(number)
    #     if "e" in float_str:
    #         base, exponent = float_str.split("e")
    #         return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    #     else:
    #         return float_str

    def length_conv(self, length, unit, distance=None):
        """Converted the length to various units if given in meters.

        Args:
            length (float or List): A given length [m].
            unit (str): Unit to convert length into.
            distance (float): Distance to object used for conversion into arc seconds.

        Returns:
            float or List: Length in another unit ['arcsec', 'au', 'pc'].
        """
        if unit == 'arcsec':
            if distance is not None:
                new_length = np.divide(
                    length, self.const['au'] * (distance / self.const['pc']))
            else:
                raise ValueError('The distance is not set. Without the distance, '
                                 'no conversion to arcseconds is available!')
        elif unit == 'au':
            new_length = np.divide(length, self.const['au'])
        elif unit == 'pc':
            new_length = np.divide(length, self.const['pc'])
        elif unit == 'microns':
            new_length = np.multiply(length, 1e6)
        elif unit == 'mm':
            new_length = np.multiply(length, 1e3)
        elif unit == 'm':
            new_length = length
        else:
            raise ValueError('There is no unit type available for conversion!')
        return new_length

    def conv_length_factor(self, unit, distance=None):
        """Calculates the conversion factor to convert a length from unit to m.

        Args:
            unit (str): Unit of the original length.
            distance (float): Distance to object used for conversion into arc seconds

        Returns:
            float: Conversion factor
        """
        if unit == 'arcsec':
            if distance is not None:
                return (self.const['au'] * (distance / self.const['pc']))
            else:
                raise ValueError('The distance is not set. Without the distance, '
                                 'no conversion to arcseconds is available!')
        elif unit == 'au':
            return self.const['au']
        elif unit == 'pc':
            return self.const['pc']
        elif unit == 'm':
            return 1.
        else:
            raise ValueError('There is no unit type available for conversion!')

    def parse(self, value, unit_type, unit=None):
        """Convert a value of a certain type (length, angle, ...) to SI or unit.

        Args:
            value (float): Value to parse/convert.
            unit_type (str): Type of the input value
                ('length', 'mass', 'angle', 'velocity', 'luminosity')
            unit (str): Unit to convert to (SI if None).
        """
        # Make lower case to ensure comparability
        if unit is not None:
            unit = unit.lower()
        # Check input unit for different quantity types
        if unit_type == 'length':
            parsed_value = self.parse_length(value)
            if unit == 'au':
                return parsed_value / self.const['au']
            elif unit == 'pc':
                return parsed_value / self.const['pc']
            elif unit == 'km':
                return parsed_value / 1e3
            elif unit == 'cm':
                return parsed_value / 1e-2
            elif unit == 'mm':
                return parsed_value / 1e-3
            elif unit == 'microns':
                return parsed_value / 1e-6
            elif unit == 'nm':
                return parsed_value / 1e-9
            elif unit == 'm' or unit is None:
                return parsed_value
            else:
                raise ValueError(
                    'There is no unit type available for conversion!')
        elif unit_type == 'velocity':
            parsed_value = self.parse_velocity(value)
            if unit == 'km/s':
                return parsed_value * 1e-3
            elif unit == 'km/h':
                return parsed_value * 3600 / 1e3
            elif unit == 'm/s' or unit is None:
                return parsed_value
            else:
                raise ValueError(
                    'There is no unit type available for conversion!')
        elif unit_type == 'mass':
            parsed_value = self.parse_mass(value)
            if unit == 'm_jup' or unit == 'mjup':
                return parsed_value / self.const['M_jup']
            elif unit == 'm_sun' or unit == 'msun':
                return parsed_value / self.const['M_sun']
            elif unit == 'kg' or unit is None:
                return parsed_value
            else:
                raise ValueError(
                    'There is no unit type available for conversion!')
        elif unit_type == 'luminosity':
            parsed_value = self.parse_luminosity(value)
            if unit == 'l_sun' or unit == 'lsun':
                return parsed_value / self.const['L_sun']
            elif unit == 'w' or unit is None:
                return parsed_value
            else:
                raise ValueError(
                    'There is no unit type available for conversion!')
        elif unit_type == 'angle':
            parsed_value = self.parse_angle(value)
            if unit == 'arcsec' or unit == 'arc_sec':
                return parsed_value * 3600
            elif unit == 'rad':
                return parsed_value / 180. * np.pi
            elif unit == 'degree' or unit is None:
                return parsed_value
            else:
                raise ValueError(
                    'There is no unit type available for conversion!')
        else:
            raise ValueError(
                'There is no unit type available for conversion!')

    def parse_length(self, length):
        """Convert input length string into length in meters.

        Args:
            length (str): String that includes the length and possibly
                the unit (e.g. 13au).
        """
        conv = 1
        if 'km' in length.lower():
            conv = 1e3
            length = length.lower().replace('km', '')
        elif 'cm' in length.lower():
            conv = 1e-2
            length = length.lower().replace('cm', '')
        elif 'mm' in length.lower():
            conv = 1e-3
            length = length.lower().replace('mm', '')
        elif 'microns' in length.lower():
            conv = 1e-6
            length = length.lower().replace('microns', '')
        elif 'nm' in length.lower():
            conv = 1e-9
            length = length.lower().replace('nm', '')
        elif 'm' in length.lower():
            length = length.lower().replace('m', '')
        elif 'au' in length.lower():
            conv = self.const['au']
            length = length.lower().replace('au', '')
        elif 'pc' in length.lower():
            conv = self.const['pc']
            length = length.lower().replace('pc', '')
        elif 'r_sun' in length.lower() or 'rsun' in length.lower():
            conv = self.const['R_sun']
            length = length.lower().replace('r_sun', '').replace('rsun', '')
        try:
            return float(length) * conv
        except ValueError:
            print('Length string ' + length + ' cannot be interpreted!')

    def parse_velocity(self, velocity):
        """Convert input velocity string into velocity in meter/seconds.

        Args:
            velocity (str): String that includes the velocity and possibly
                the unit (e.g. 13m/s).
        """
        conv = 1
        if 'km/s' in velocity.lower():
            conv = 1e3
            velocity = velocity.lower().replace('km/s', '')
        elif 'km/h' in velocity.lower():
            conv = 1e3 / 3600.
            velocity = velocity.lower().replace('km/h', '')
        elif 'm/s' in velocity.lower():
            velocity = velocity.lower().replace('m/s', '')
        try:
            return float(velocity) * conv
        except ValueError:
            print('Velocity string ' + velocity + ' cannot be interpreted!')

    def parse_mass(self, mass):
        """Convert input mass string into mass in kg.

        Args:
            mass (str): String that includes the mass and possibly
                the unit (e.g. 13m_sun).
        """
        conv = 1
        if 'kg' in mass.lower():
            mass = mass.lower().replace('kg', '')
        elif 'mg' in mass.lower():
            conv = 1e-6
            mass = mass.lower().replace('mg', '')
        elif 'g' in mass.lower():
            conv = 1e-3
            mass = mass.lower().replace('g', '')
        elif 'm_jup' in mass.lower() or 'mjup' in mass.lower():
            conv = self.const['M_jup']
            mass = mass.lower().replace('m_jup', '').replace('mjup', '')
        elif 'm_sun' in mass.lower() or 'msun' in mass.lower():
            conv = self.const['M_sun']
            mass = mass.lower().replace('m_sun', '').replace('msun', '')
        try:
            return float(mass) * conv
        except ValueError:
            print('Mass string ' + mass + ' cannot be interpreted!')

    def parse_angle(self, angle):
        """Convert input angle string into angle in degree.

        Args:
            angle (str): String that includes the angle and possibly
                the unit (e.g. 13arcsec).
        """
        conv = 1
        if 'degree' in angle.lower():
            angle = angle.lower().replace('degree', '')
        elif 'rad' in angle.lower():
            conv = 180. / np.pi
            angle = angle.lower().replace('rad', '')
        elif 'arcsec' in angle.lower():
            conv = 1 / 3600.
            angle = angle.lower().replace('arcsec', '')
        try:
            return float(angle) * conv
        except ValueError:
            print('Angle string ' + angle + ' cannot be interpreted!')

    def parse_luminosity(self, luminosity):
        """Convert input luminosity string into luminosity in Watt.

        Args:
            luminosity (str): String that includes the luminosity and possibly
                the unit (e.g. 13L_sun).
        """
        conv = 1
        if 'mW' in luminosity:
            conv = 1e-3
            luminosity = luminosity.lower().replace('mw', '')
        elif 'kW' in luminosity:
            conv = 1e3
            luminosity = luminosity.lower().replace('kw', '')
        elif 'MW' in luminosity:
            conv = 1e6
            luminosity = luminosity.lower().replace('mw', '')
        elif 'w' in luminosity.lower():
            luminosity = luminosity.lower().replace('w', '')
        elif 'l_sun' in luminosity.lower() or 'lsun' in luminosity.lower():
            conv = self.const['L_sun']
            luminosity = luminosity.lower().replace('l_sun', '').replace('lsun', '')
        try:
            return float(luminosity) * conv
        except ValueError:
            print('Luminosity string ' + luminosity + ' cannot be interpreted!')

    # @staticmethod
    # def angle_conv(angle, unit):
    #     """Converted the angle in various units into radians.

    #     Args:
    #         angle (float or List): A given angle.
    #         unit (str): Unit of angle.

    #     Returns:
    #         float or List: Angle in radians.
    #     """
    #     if unit == 'arcsec':
    #         new_angle = np.divide(np.multiply(angle, np.pi), 3600. * 180.)
    #     elif unit == 'arcmin':
    #         new_angle = np.divide(np.multiply(angle, np.pi), 60. * 180.)
    #     elif unit == 'degree':
    #         new_angle = np.divide(np.multiply(angle, np.pi), 180.)
    #     else:
    #         raise ValueError('Unit of angle not understood!')
    #     return new_angle

    # @staticmethod
    # def get_velocity(vch, nr_channel, max_velocity):
    #     return vch / (nr_channel - 1) * (2 * max_velocity) - max_velocity

    # @staticmethod
    # def get_vector_color(tbldata, cmap, cmap_scaling=None, vec_color=None):
    #     """
    #     Args:
    #         tbldata: Numpy array of the imshow plot.
    #         cmap(str): Colormap of of the imshow plot
    #         cmap_scaling (str): Scaling of the colormap.
    #         vec_color (str): Overwrite automatic choice by user.

    #     Returns:
    #         str: Name of the color that has to be used by the vector plots
    #     """
    #     # If the user input sets the vector color, use it
    #     if vec_color is not None:
    #         return vec_color
    #     # White has a good contrast for viridis
    #     '''
    #     if cmap in ['viridis', 'magma']:
    #         return 'white'
    #     '''
    #     # Return black as the vector color if no or negative data is used
    #     if cmap_scaling == 'log' and np.min(tbldata) < 0.:
    #         return 'black'
    #     # Get color that fits well with the colorplot
    #     import matplotlib.cm as cm
    #     import copy
    #     colormap = copy.copy(cm.get_cmap(cmap))
    #     if cmap_scaling == 'log':
    #         mean_value = 10 ** ((np.log10(np.max(tbldata)) +
    #                              np.log10(np.min(tbldata[np.where(tbldata > 0.)]))) / 2.)
    #     else:
    #         mean_value = (np.max(tbldata) + np.min(tbldata)) / 2.
    #     if np.mean(tbldata) >= mean_value:
    #         return colormap(0), colormap(256)
    #     else:
    #         return colormap(256), colormap(0)

    # def flux_2_mag(self, flux, filter_system='V'):
    #     """Calculates the magnitude of a given flux in a given filter system.

    #     Args:
    #         flux (float): Flux that will be converted [Jy].
    #         filter_system (str): Name of the filter system used to convert
    #             the flux.

    #     Returns:
    #         float: Magnitude in the chosen filter system.
    #     """
    #     if filter_system in self.conversion_factor.keys():
    #         if float(flux) != 0:
    #             magnitude = 2.5 * \
    #                 np.log10(
    #                     self.conversion_factor[filter_system] / float(flux))
    #         else:
    #             raise ValueError(
    #                 'Flux is zero and cannot be converted to a magnitude!')
    #     else:
    #         raise ValueError('The chosen filter system does not exists!')
    #     return magnitude

    # def mag_2_flux(self, magnitude, filter_system='V'):
    #     """Calculates the flux of a given magnitude in a given filter system.

    #     Args:
    #         magnitude (float): Magnitude that will be converted [mag].
    #         filter_system (str): Name of the filter system used to convert
    #             the magnitude.

    #     Returns:
    #         float: flux in the chosen filter system.
    #     """
    #     if filter_system in self.conversion_factor.keys():
    #         flux = 10 ** (magnitude / -2.5) * \
    #             self.conversion_factor[str(filter_system)]
    #     else:
    #         raise ValueError('The chosen filter system does not exists!')
    #     return flux

    # @staticmethod
    # def mbol_2_lum(mbol):
    #     """Calculates the luminosity of an object from its bolometric magnitude.

    #     Args:
    #         mbol (float): Bolometric magnitude [mag].

    #     Returns:
    #         float: Luminosity [L_sun].
    #     """
    #     luminosity = 10 ** ((4.7554 - mbol) / 2.5) * 3.826e26
    #     return luminosity

    # @staticmethod
    # def luminosity_to_radius(luminosity, temperature):
    #     """Calculates the radius for a certain luminosity of a black body.

    #     Args:
    #         luminosity (float): Luminosity of a black body [W].
    #         temperature (float): Temperature of a black body [K].

    #     Returns:
    #         float: Radius [m].
    #     """
    #     radius = np.sqrt(luminosity / (4. * np.pi *
    #                                    5.670367e-08 * temperature ** 4.))
    #     return radius

    # @staticmethod
    # def best_vel_map_distribution(nr_channels, vel_map_plots):
    #     """Get a reasonable distribution of channels for a vel_map plot.

    #     Args:
    #         nr_channels (int): Number of channels of vel_map.
    #         vel_map_plots (int): No more than this should be the amount of channels.
    #     """
    #     channel_skip = int((nr_channels - 1) / (vel_map_plots - 1))
    #     leftover = (nr_channels - 1) % (vel_map_plots - 1)
    #     if leftover % 2 == 0:
    #         start_channel = int(leftover / 2.0)
    #         vch_list = list(range(start_channel, nr_channels, channel_skip))
    #     else:
    #         if nr_channels % 2 == 1:
    #             start_channel = int((leftover - 1) / 2.0)
    #             vch_list = list(range(start_channel, int(
    #                 (nr_channels - 1) / 2.0), channel_skip))
    #             vch_list = list(np.append(vch_list, np.subtract(
    #                 nr_channels - 1, vch_list[::-1])))
    #         else:
    #             raise ValueError(
    #                 'Odd number of vel_map_plots does not work with even number of channels!')

    #     if vel_map_plots > 5:
    #         nr_y_images = int(np.sqrt(len(vch_list)))
    #         nr_x_images = int(np.ceil(len(vch_list) / nr_y_images))
    #     else:
    #         nr_y_images = 1
    #         nr_x_images = int(len(vch_list))
    #     return vch_list, nr_x_images, nr_y_images

    # def planet_lum(self, M, M_acc, R):
    #     """Calculates the luminosity of a circumplanetary disk
    #     (according to Zhu 2015).

    #     Args:
    #         M (float): Mass of the planet [kg].
    #         M_acc(float): Accretion rate of the circumplanetary disk [M_sun/yr].
    #         R (float): Radius of the planet [m].

    #     Returns:
    #         float: Luminosity of the circumplanetary disk [L_sun].
    #     """
    #     M_acc_si = M_acc * self.const['M_sun'] * (1/31536000.)
    #     res = self.const['G'] * M * M_acc_si / (2. * R)
    #     return res

    # def planet_temp(self, M, M_acc, R):
    #     """Calculates the effective temperature of a circumplanetary disk.

    #     Args:
    #         M (float): Mass of the planet [kg].
    #         M_acc(float): Accretion rate of the circumplanetary disk [M_sun/yr].
    #         R (float): Radius of the planet [m].

    #     Returns:
    #         float: Effective temperature of the circumplanetary disk [K].
    #     """
    #     M_acc_si = M_acc * self.const['M_sun'] * (1/31536000.)
    #     res = 3 * self.const['G'] * M * M_acc_si / \
    #         (8 * np.pi * self.const['sigma_sb'] ** 1 * R ** 3)
    #     return 0.488 * res ** (1. / 4.)

    # @staticmethod
    # def apply_inclination(pos, inclination, inc_PA, inc_offset=0, inv=False):
    #     """Calculate ellipse caused by inclination.
    #     """
    #     inc_offset *= np.linalg.norm(pos)
    #     dir_inc = [np.cos(inc_PA), np.sin(inc_PA)]
    #     dir_perp = [-np.sin(inc_PA), np.cos(inc_PA)]
    #     if not inv:
    #         return np.multiply(np.dot(pos, dir_perp), dir_perp) + \
    #             np.multiply(np.dot(pos, dir_inc) * np.cos(inclination), dir_inc) + \
    #             np.multiply(inc_offset, dir_inc)
    #     else:
    #         shifted_pos = pos - np.multiply(inc_offset, dir_inc)
    #         return np.multiply(np.dot(shifted_pos, dir_perp), dir_perp) + \
    #             np.multiply(np.dot(shifted_pos, dir_inc) /
    #                         np.cos(inclination), dir_inc)

    # @staticmethod
    # def tbldata_2_radial_profile(tbldata, bins=50):
    #     """Calculates a radial brightness profile from 2D array.

    #     Args:
    #         tbldata : A 2D array with a size of n x n pixel.
    #         bins (int): Number of bins for radial profile.

    #     Returns:
    #         Radial profile as 1D array.
    #     """
    #     # Check is 2D array is symmetric
    #     if len(tbldata[0, :]) != len(tbldata[:, 0]):
    #         raise ValueError('2D array has not a symmetric shape!')
    #     # Get length and half length of each axis of the 2D array
    #     n = len(tbldata[0, :])
    #     n_half = len(tbldata[0, :]) / 2.
    #     # Create radial profile for output
    #     radial_profile = np.zeros(bins)
    #     # Create array to save the number of pixel per bin
    #     number_of_pixel_per_bin = np.zeros(bins)
    #     # For each Pixel get information about radial position and array value
    #     for i_x in range(n):
    #         for i_y in range(n):
    #             #: float: Radial position in [m]
    #             radius = np.sqrt((i_x - n_half) ** 2 +
    #                              (i_y - n_half) ** 2) / n_half
    #             #: float: Radial position in bins [0, bins - 1]
    #             radius_bin = int(radius * (bins - 1))
    #             if radius_bin < bins:
    #                 radial_profile[radius_bin] += tbldata[i_x, i_y]
    #                 number_of_pixel_per_bin[radius_bin] += 1
    #     # Loop over each bin to normalize it on the value per pixel
    #     for i_bin in range(bins):
    #         radial_profile[i_bin] /= number_of_pixel_per_bin[i_bin]
    #     return radial_profile

    # @staticmethod
    # def tbldata_2_polarization_intensity(tbldata, bins=50, i_index=0, p_index=7):
    #     """Calculates radial profile of the degree of polarization P
    #     in dependence on the normalized Intensity.

    #     Args:
    #         tbldata : A 3D array with a size of q x n x n pixel.
    #             (q is number of quantities and needs to include I at i_index and P at p_index)
    #         bins (int): Number of bins for radial profile.
    #         i_index (int): Index of tbldata where the 2D intensity data is.
    #         p_index (int): Index of tbldata where the 2D degree of polarization data is.

    #     Returns:
    #         Radial profile as 1D array.
    #     """
    #     # Check is 2D array is symmetric
    #     if len(tbldata[0, 0, :]) != len(tbldata[0, :, 0]):
    #         raise ValueError('Array has not a symmetric shape!')
    #     # Get Length of each axis of the 2D array
    #     n = len(tbldata[0, :])
    #     # Create radial profile for output
    #     radial_profile = np.zeros(bins)
    #     # Create array to save the number of pixel per bin
    #     number_of_pixel_per_bin = np.zeros(bins)
    #     # Calculate the maximum intensity for normalization
    #     max_intensity = np.max(tbldata[0, :, :]).item()
    #     # For each Pixel get information about radial position and array value
    #     for i_x in range(n):
    #         for i_y in range(n):
    #             if tbldata[p_index, i_x, i_y] > 0:
    #                 #: float: Radial position in [m]
    #                 intensity = tbldata[i_index, i_x, i_y] / max_intensity
    #                 #: float: Radial position in bins [0, bins - 1]
    #                 radius_bin = int(intensity * (bins - 1))
    #                 radial_profile[radius_bin] += tbldata[p_index, i_x, i_y]
    #                 number_of_pixel_per_bin[radius_bin] += 1
    #     # Loop over each bin to normalize it on the value per pixel
    #     for i_bin in range(bins):
    #         if radial_profile[i_bin] == 0 or number_of_pixel_per_bin[i_bin] == 0:
    #             radial_profile[i_bin] = radial_profile[i_bin - 1] \
    #                 + (radial_profile[i_bin - 1] - radial_profile[i_bin - 2])
    #         else:
    #             radial_profile[i_bin] /= number_of_pixel_per_bin[i_bin]
    #     return radial_profile

    # @staticmethod
    # def velocity_2_frequency(velocity, f_0):
    #     """Calculates the frequency corresponding to a given velocity.

    #     Args:
    #         velocity (float): Velocity for conversion [m/s].
    #         f_0 (float) : Rest frequency [Hz].

    #     Returns:
    #         float: Frequency [Hz].
    #     """
    #     frequency = velocity * (f_0 / 299792458.0)
    #     return frequency

    # @staticmethod
    # def frequency_2_velocity(frequency, f_0):
    #     """Calculates the velocity corresponding to a given frequency.

    #     Args:
    #         frequency (float): Frequency for conversion [Hz].
    #         f_0 (float) : Rest frequency [Hz].

    #     Returns:
    #         float: Velocity [m/s].
    #     """
    #     velocity = frequency * (299792458.0 / f_0)
    #     return velocity

    # @staticmethod
    # def scuba_2_obs(f, sigma, wavelength):
    #     """Calculates the observation time of the SCUBA-2 instrument
    #     with the Daisy mapping mode.

    #     Notes:
    #         Link: www.eaobservatory.org/jcmt/instrumentation/continuum/scuba-2/time-and-sensitivity/

    #     Args:
    #         f (int): Sampling factor.
    #             f = ( pixel size requested / default pixel size )^2
    #         sigma (float): 1-sigma depth at the given wavelength [Jy].
    #         wavelength: The wavelength at which the sensitivity is calculated [m].

    #     Returns:
    #         float: Observation time to reach the sensitivity.
    #     """
    #     if wavelength == 4.50E-4:
    #         return 1. / f * ((689. / 0.197 - 118.) * 1. / (sigma * 1e3)) ** 2
    #     elif wavelength == 8.50E-4:
    #         return 1. / f * ((189. / 0.720 - 48.) * 1. / (sigma * 1e3)) ** 2
    #     else:
    #         raise ValueError(
    #             'The chosen wavelength is not provided by the SCUBA-2 instrument!')

    @staticmethod
    def spherical_to_cartesian(spherical_coord):
        """Calculates cartesian coordinates from spherical ones.

        Args:
            spherical_coord (List[float, float, float]): Spherical coordinates.
                radius, theta, phi

        Returns:
            List[float, float, float]: Cartesian coordinates
        """
        cartesian_coord = np.zeros(3)
        cartesian_coord[0] = spherical_coord[0] * \
            np.sin(spherical_coord[1]) * np.cos(spherical_coord[2])
        cartesian_coord[1] = spherical_coord[0] * \
            np.sin(spherical_coord[1]) * np.sin(spherical_coord[2])
        cartesian_coord[2] = spherical_coord[0] * \
            np.cos(spherical_coord[1])
        return cartesian_coord

    @staticmethod
    def cartesian_to_spherical(cartesian_coord):
        """Calculates spherical coordinates from cartesian ones.

        Args:
            cartesian_coord (List[float, float, float]): Cartesian coordinates.
                x, y, z

        Returns:
            List[float, float, float]: Spherical coordinates
        """
        spherical_coord = np.zeros(3)
        if np.linalg.norm(cartesian_coord[:]) != 0:
            spherical_coord[0] = np.linalg.norm(cartesian_coord[:])
            spherical_coord[1] = np.arccos(cartesian_coord[2] / spherical_coord[0])
            spherical_coord[2] = np.arctan2(
                                    cartesian_coord[1], cartesian_coord[0])
        return spherical_coord

    @staticmethod
    def cylindrical_to_cartesian(cylindrical_coord):
        """Calculates cartesian coordinates from cylindrical ones.

        Args:
            cylindrical_coord (List[float, float, float]): Cylindrical coordinates.
                radius, phi, z

        Returns:
            List[float, float, float]: Cartesian coordinates
        """
        cartesian_coord = np.zeros(3)
        cartesian_coord[0] = cylindrical_coord[0] * \
            np.cos(cylindrical_coord[1])
        cartesian_coord[1] = cylindrical_coord[0] * \
            np.sin(cylindrical_coord[1])
        cartesian_coord[2] = cylindrical_coord[2]
        return cartesian_coord

    @staticmethod
    def cartesian_to_cylindrical(cartesian_coord):
        """Calculates cylindrical coordinates from cartesian ones.

        Args:
            cartesian_coord (List[float, float, float]): Cartesian coordinates.
                x, y, z

        Returns:
            List[float, float, float]: Cylindrical coordinates
        """
        cylindrical_coord = np.zeros(3)
        cylindrical_coord[0] = np.linalg.norm(cartesian_coord[0:2])
        cylindrical_coord[1] = np.arctan2(
            cartesian_coord[1], cartesian_coord[0]) + np.pi
        cylindrical_coord[2] = cartesian_coord[2]
        return cylindrical_coord

    # @staticmethod
    # def spherical_to_cartesian_direction(spherical_coord, spherical_direction):
    #     """Calculates cartesian coordinates from spherical ones.

    #     Args:
    #         spherical_coord (List[float, float, float]): Spherical coordinates.
    #             radius, theta, phi
    #         spherical_direction (List[float, float, float]): Direction in Spherical coordinates.

    #     Returns:
    #         List[float, float, float]: Cartesian coordinates
    #     """
    #     e_r = np.zeros(3)
    #     e_r[0] = np.sin(spherical_coord[1]) * np.cos(spherical_coord[2])
    #     e_r[1] = np.sin(spherical_coord[1]) * np.sin(spherical_coord[2])
    #     e_r[2] = np.cos(spherical_coord[1])

    #     e_t = np.zeros(3)
    #     e_t[0] = np.cos(spherical_coord[1]) * np.cos(spherical_coord[2])
    #     e_t[1] = np.cos(spherical_coord[1]) * np.sin(spherical_coord[2])
    #     e_t[2] = -np.sin(spherical_coord[1])

    #     e_p = np.zeros(3)
    #     e_p[0] = -np.sin(spherical_coord[2])
    #     e_p[1] = np.cos(spherical_coord[2])
    #     e_p[2] = 0

    #     cartesian_direction = np.zeros(3)
    #     cartesian_direction[:] += spherical_direction[0] * e_r[:]
    #     cartesian_direction[:] += spherical_direction[1] * e_t[:]
    #     cartesian_direction[:] += spherical_direction[2] * e_p[:]
    #     return cartesian_direction

    # @staticmethod
    # def rotate_coord_system(position, rotation_axis, rotation_angle, inv=False):
    #     """Converts the postion coordinates into a coordinate system that is
    #     rotated around an rotation_axis.

    #     Args:
    #         position (List[float, float, float]): Position in cartesian coordinates.
    #         rotation_axis (List[float, float, float]): Rotation axis in cartesian coordinates.
    #         rotation_angle (float): Angle to rotate around [rad].
    #         inv (bool): Invert rotation?

    #     Note:
    #         Source: https://en.wikipedia.org/wiki/Rotation_matrix

    #     Returns:
    #         List[float, float, float]: Rotated cartesian coordinates
    #     """
    #     # Ignore if angle is zero
    #     if rotation_angle == 0:
    #         return position

    #     rotation_axis /= np.linalg.norm(rotation_axis)
    #     (u_x, u_y, u_z) = rotation_axis
    #     rot_cos = np.cos(rotation_angle)
    #     rot_sin = np.sin(rotation_angle)
    #     rotation_matrix = np.array([[rot_cos + u_x ** 2 * (1 - rot_cos),
    #                                  u_x * u_y * (1 - rot_cos) - u_z * rot_sin,
    #                                  u_x * u_z * (1 - rot_cos) + u_y * rot_sin],
    #                                 [u_y * u_x * (1 - rot_cos) + u_z * rot_sin,
    #                                  rot_cos + u_y ** 2 * (1 - rot_cos),
    #                                  u_y * u_z * (1 - rot_cos) - u_x * rot_sin],
    #                                 [u_z * u_x * (1 - rot_cos) - u_y * rot_sin,
    #                                  u_z * u_y * (1 - rot_cos) + u_x * rot_sin,
    #                                  rot_cos + u_z ** 2 * (1 - rot_cos)]])
    #     if inv:
    #         rotation_matrix = rotation_matrix.T
    #     rotated_position = np.dot(position, rotation_matrix)
    #     return rotated_position

    @staticmethod
    def sin_list(start, stop, total_number):
        """Calculates sinus distribution between two values.

        Args:
            start (float): starting value
            stop (float): last value
            total_number (int): total amount of distributed values.

        Returns:
             number_list (list): Distributed numbers
        """
        number_list = np.zeros(total_number + 1)
        number_list[0] = start
        number_list[total_number] = stop

        inter = (stop - start)
        dang = np.pi / total_number
        mid = (total_number - 0.5) / 2.0

        for i_x in range(1, total_number):
            if i_x <= mid:
                number_list[i_x] = start + inter * (0.5 * np.sin(i_x * dang))
            else:
                number_list[i_x] = start + inter * \
                    (1 - 0.5 * np.sin(i_x * dang))
        return number_list

    @staticmethod
    def lin_list(start, stop, total_number):
        """Calculates linear distribution between two values.

        Args:
            start (float): starting value
            stop (float): last value
            total_number (int): total amount of distributed values.

        Returns:
             number_list (list): Distributed numbers
        """
        number_list = np.zeros(total_number + 1)
        number_list[0] = start
        number_list[total_number] = stop

        dx = (stop - start) / total_number

        for i_x in range(1, total_number):
            number_list[i_x] = start + i_x * dx
        return number_list

    @staticmethod
    def exp_list(start, stop, total_number, base):
        """Calculates exponential distribution between two values.

        Args:
            start (float): starting value
            stop (float): last value
            total_number (int): total amount of distributed values.
            base (float): distribution factor

        Returns:
             number_list (list): Distributed numbers
        """
        number_list = np.zeros(total_number + 1)
        number_list[0] = start
        number_list[total_number] = stop

        if base > 1:
            dx = (stop - start) * (base - 1.0) / (pow(base, total_number) - 1)

            for i_x in range(0, total_number):
                number_list[i_x] = start + dx * \
                    (pow(base, i_x) - 1) / (base - 1.0)
        else:
            raise ValueError('only positive exp bases are allowed!')
        return number_list

    def exp_list_sym(self, start, stop, total_number, base):
        """Calculates exponential distribution between two values.

        Args:
            start (float): starting value
            stop (float): last value
            total_number (int): total amount of distributed values.
            base (float): distribution factor

        Returns:
             number_list (list): Distributed numbers
        """
        number_list = np.zeros(total_number + 1)
        number_list[0] = start
        number_list[total_number] = stop
        tmp_mid = start + 0.5 * (stop - start)
        if total_number % 2 == 0:
            midN = int(total_number / 2)
            tmp_list = self.exp_list(tmp_mid, stop, midN, base)
            for i_x in range(midN, total_number + 1):
                number_list[i_x] = tmp_list[i_x - midN]
            diff = 0
            for i_x in range(midN):
                diff += tmp_list[i_x + 1] - tmp_list[i_x]
                number_list[midN - 1 - i_x] = tmp_mid - diff
        else:
            midN = int((total_number - 1) / 2)
            tmpN = midN + 1
            tmp_list = self.exp_list(tmp_mid, stop, tmpN, base)
            for i_x in range(tmpN, total_number + 1):
                number_list[i_x] = tmp_list[i_x - midN]
            diff = 0
            for i_x in range(midN):
                diff += tmp_list[i_x + 1] - tmp_list[i_x]
                number_list[midN - i_x] = tmp_mid - diff
        return number_list

    # @staticmethod
    # def angle_from_stokes(stokes_q, stokes_u):
    #     """Calculates the polarization angle from Q and U Stokes component.

    #     Args:
    #         stokes_q (float): Q-Stokes component [Jy].
    #         stokes_u (float): U-Stokes component [Jy].

    #     Returns:
    #         float: Polarization angle.
    #     """
    #     #: float: Polarization angle from Stokes Q component
    #     q_angle = 0.
    #     if stokes_q >= 0:
    #         q_angle = np.pi / 2.
    #     #: float: Polarization angle from Stokes U component
    #     u_angle = 0.
    #     if stokes_u >= 0:
    #         u_angle = np.pi / 4.
    #     elif stokes_u < 0:
    #         if stokes_q >= 0:
    #             u_angle = np.pi * 3. / 4.
    #         elif stokes_q < 0:
    #             u_angle = -np.pi / 4.
    #     #: float: x vector components from both angles
    #     x = abs(stokes_q) * np.sin(q_angle)
    #     x += abs(stokes_u) * np.sin(u_angle)
    #     #: float: y vector components from both angles
    #     y = abs(stokes_q) * np.cos(q_angle)
    #     y += abs(stokes_u) * np.cos(u_angle)
    #     # Define a global direction of the polarization vector since polarization vectors
    #     # are ambiguous in both directions.
    #     if x < 0:
    #         x *= -1.0
    #         y *= -1.0
    #     #: float: Polarization angle calculated from Q and U components
    #     pol_angle = np.arctan2(y, x)
    #     return pol_angle

    def kepler_rotation(self, position, stellar_mass):
        """Calculates Kepler rotation velocity.

        Args:
            position (List[float, float, float]): Position in model space.
            stellar_mass (float): Mass of central stellar object [M_sun].

        Returns:
            List[float, float, float]: Velocity at the given position.
        """
        #: float: Cylindrical radius
        radius_cy = np.sqrt(position[0] ** 2 + position[1] ** 2)
        #: float: Kepler constant ( v=sqrt(GM/a) )
        kepler_const = (self.const['G'] * stellar_mass *
                        self.const['M_sun'] / radius_cy) ** 0.5
        velocity = [-1.0 * position[1] / radius_cy * kepler_const,
                    position[0] / radius_cy * kepler_const, 0.]
        return velocity

    # @staticmethod
    # def velocity_shift_from_profile(profile_i, profile_v, vel_channel_width):
    #     """Calculate magnetic field strength from profile.
    #     Args:
    #         profile_i (list): Intensity profile.
    #         profile_v (list): Circular polarization profile.
    #         vel_channel_width (float): Width between channels in [m/s].

    #     Notes:
    #         See Crutcher et al. (1993)

    #     Returns:
    #         Float: Velocity shift in [m/s] calculated from I and V profiles.
    #     """
    #     # Derivative of the intensity I with respect to the velocity
    #     intensity_derivative = np.gradient(
    #         profile_i, vel_channel_width, axis=0)

    #     velocity_shift = np.zeros(intensity_derivative.shape[1:])
    #     for i in np.ndindex(intensity_derivative.shape[1:]):
    #         # Define fit function to get proportionality factor
    #         def fit_function(vch, prop_factor): return prop_factor * \
    #             intensity_derivative[(vch,) + i]

    #         # A total zero line profile corresponds to a magnetic field of zero
    #         if (intensity_derivative[(slice(None),) + i] == 0.).all():
    #             velocity_shift[i] = 0.
    #         else:
    #             # Calculate best-fit parameter
    #             popt, pcov = curve_fit(f=fit_function, xdata=range(profile_i.shape[0]), ydata=profile_v[(slice(None),) + i],
    #                                    sigma=1e-3 * np.ones(profile_i.shape[0]))
    #             velocity_shift[i] = popt[0]
    #     # Get magnetic field strength from Zeeman information of used species
    #     return velocity_shift

    # @staticmethod
    # def relative_line_strength_zero(j, m, transition):
    #     """Calculates the relative line strength of a transition between two Zeeman sublevels
    #     (if Delta j = 0).

    #     Args:
    #         j (float): Quantum number of the effective angular momentum.
    #         m (float): Quantum number of the projection of j on the magnetic field direction.
    #         transition (int): Type of transition between the Zeeman sublevels.

    #     Returns:
    #         Float: The relative line strength of a transition between two Zeeman sublevels
    #     """
    #     if transition == 0:
    #         return 3 * m ** 2 / (j * (j + 1) * (2 * j + 1))
    #     elif transition == +1:
    #         return 3 * (j - m) * (j + 1 + m) / (4 * j * (j + 1) * (2 * j + 1))
    #     elif transition == -1:
    #         return 3 * (j + m) * (j + 1 - m) / (4 * j * (j + 1) * (2 * j + 1))

    # @staticmethod
    # def relative_line_strength_plus(j, m, transition):
    #     """Calculates the relative line strength of a transition between two Zeeman sublevels
    #     (if Delta j = +1).

    #     Args:
    #         j (float): Quantum number of the effective angular momentum.
    #         m (float): Quantum number of the projection of j on the magnetic field direction.
    #         transition (int): Type of transition between the Zeeman sublevels.

    #     Returns:
    #         Float: The relative line strength of a transition between two Zeeman sublevels
    #     """
    #     if transition == 0:
    #         return 3 * ((j + 1) ** 2 - m ** 2) / (2 * (j + 1) * (2 * j + 1) * (2 * j + 3))
    #     elif transition == +1:
    #         return 3 * ((j + 1 + m) * (j + 2 + m)) / (4 * (j + 1) * (2 * j + 1) * (2 * j + 3))
    #     elif transition == -1:
    #         return 3 * ((j + 1 - m) * (j + 2 - m)) / (4 * (j + 1) * (2 * j + 1) * (2 * j + 3))

    # @staticmethod
    # def relative_line_strength_minus(j, m, transition):
    #     """Calculates the relative line strength of a transition between two Zeeman sublevels
    #     (if Delta j = -1).

    #     Args:
    #         j (float): Quantum number of the effective angular momentum.
    #         m (float): Quantum number of the projection of j on the magnetic field direction.
    #         transition (int): Type of transition between the Zeeman sublevels.

    #     Returns:
    #         Float: The relative line strength of a transition between two Zeeman sublevels
    #     """
    #     if transition == 0:
    #         return 3 * (j ** 2 - m ** 2) / (2 * j * (2 * j - 1) * (2 * j + 1))
    #     elif transition == +1:
    #         return 3 * ((j - m) * (j - 1 - m)) / (4 * j * (2 * j - 1) * (2 * j + 1))
    #     elif transition == -1:
    #         return 3 * ((j + m) * (j - 1 + m)) / (4 * j * (2 * j - 1) * (2 * j + 1))

    # @staticmethod
    # def lande_g_cn(n, j, f):
    #     """Calculates the lande g-factor for CN.

    #     Args:
    #         n (float): Quantum number of the total orbital angular momentum.
    #         j (float): Quantum number of the effective angular momentum.
    #         f (float): Quantum number of the effective atomic angular momentum.

    #     Returns:
    #         Float: Lande g-factor of the chosen hyperfine sublevel.
    #     """
    #     #: float: Spin angular momentum quantum number of CN
    #     s = 0.5
    #     #: float: Atomic angular momentum quantum number of CN
    #     i = 1
    #     #: float: Orbital Lande factor of CN
    #     g_j = (j * (j + 1) + s * (s + 1) - n * (n + 1)) / (j * (j + 1))
    #     #: float: Atomic Lande factor of CN
    #     g_f = g_j * (f * (f + 1) + j * (j + 1) -
    #                  i * (i + 1)) / (2 * f * (f + 1))
    #     return g_f

    # @staticmethod
    # def lande_g_co(n, j):
    #     """Calculates the lande g-factor for CO.

    #     Args:
    #         n (float): Quantum number of the total orbital angular momentum.
    #         j (float): Quantum number of the effective angular momentum.

    #     Returns:
    #         Float: Lande g-factor of the chosen hyperfine sublevel.
    #     """
    #     return 0

    # @staticmethod
    # def lande_g_so(n, j):
    #     """Calculates the lande g-factor for SO.

    #     Args:
    #         n (float): Quantum number of the total orbital angular momentum.
    #         j (float): Quantum number of the effective angular momentum.

    #     Returns:
    #         Float: Lande g-factor of the chosen zeeman sublevel.
    #     """
    #     if j == 0 or n == 0:
    #         return 0.0
    #     delta = j - n
    #     if delta == 0:
    #         #: float: Orbital Lande factor of CN
    #         g_jn = 2 / (j * (j + 1))
    #     else:
    #         #: float: Clark & Johnson 1974 / Tiemann 1974
    #         const = 7.3528
    #         # Calculate the parts of the Atomic Lande factor of CN
    #         part1 = 1 / (2 * j * (j + 1))
    #         part2 = (2 * j + 1) ** 2 - const
    #         part3 = (1 - const) ** 2 + 4 * j * (j + 1)
    #         #: float: Orbital Lande factor of CN
    #         g_jn = 2 * part1 * (1 + delta * part2 / part3 ** 0.5)
    #     return g_jn

    # @staticmethod
    # def lande_g_ccs(n, j):
    #     """Calculates the lande g-factor for CCS.

    #     Args:
    #         n (float): Quantum number of the total orbital angular momentum.
    #         j (float): Quantum number of the effective angular momentum.

    #     Returns:
    #         Float: Lande g-factor of the chosen hyperfine sublevel.
    #     """
    #     #: list: Eigenvector coefficients for the F_1 and F_3 components for CCS
    #     # (Shinnaga & Yamamoto 2000)
    #     p_list = [0, 0.6558, 0.7534, 0.8097, 0.8490, 0.8782, 0.9003, 0.9175, 0.9308, 0.9414, 0.9499, 0.9568, 0.9623,
    #               0.9670, 0.9708, 0.9741, 0.9768, 0.9792, 0.9812, 0.9829, 0.9844, 0.9858, 0.9869, 0.9880, 0.9889,
    #               0.9897, 0.9904, 0.9911, 0.9917, 0.9922, 0.9927]
    #     q_list = [1, 0.7549, 0.6576, 0.5868, 0.5283, 0.4783, 0.4352, 0.3978, 0.3654, 0.3372, 0.3126, 0.2909, 0.2718,
    #               0.2548, 0.2397, 0.2262, 0.2141, 0.2031, 0.1931, 0.1841, 0.1758, 0.1682, 0.1612, 0.1547, 0.1488,
    #               0.1432, 0.1381, 0.1333, 0.1288, 0.1246, 0.1207]
    #     #: float: Spin Lande factor of CCS
    #     g_s = 2.00232
    #     #: float: Orbital Lande factor of CCS
    #     g_j = g_s
    #     if (j - 1) < 0:
    #         g_j = 0.
    #     else:
    #         g_j *= p_list[int(j)] ** 2 / ((j - 1) + 1) - \
    #             q_list[int(j)] ** 2 / (j + 1)
    #     return g_j

    # @staticmethod
    # def lande_g_h1(l, f):
    #     """Calculates the lande g-factor for atomic hydrogen.

    #     Args:
    #         l (float): Quantum number of the total orbital angular momentum.
    #         f (float): Quantum number of the effective atomic angular momentum.

    #     Returns:
    #         Float: Lande g-factor of the chosen hyperfine sublevel.
    #     """
    #     #: float: Spin angular momentum quantum number of H1
    #     s = 0.5
    #     #: float: Effective angular momentum quantum number of H1
    #     j = l + s
    #     #: float: Atomic angular momentum quantum number of H1
    #     i = 0.5
    #     #: float: Orbital Lande factor of H1
    #     g_j = 3. / 2. + (s * (s + 1) - l * (l + 1)) / (2 * j * (j + 1))
    #     if f > 0:
    #         #: float: Atomic Lande factor of H1
    #         g_f = g_j * (f * (f + 1) - i * (i + 1) +
    #                      j * (j + 1)) / (2 * f * (f + 1))
    #     else:
    #         g_f = 0.0
    #     return g_f

    # @staticmethod
    # def get_weight_from_qn(database_code, qn_1, qn_2, qn_3, qn_4, qn_5, qn_6):
    #     """Calculates the degeneracy/weight from quantum numbers of the given species.

    #     Args:
    #         database_code (int): Code of the species for which the degeneracy/weight will be calculated.
    #         qn_1 (int): Quantum number 1 (N).
    #         qn_2 (int): Quantum number 2 (K).
    #         qn_3 (int): Quantum number 3 (v).
    #         qn_4 (int): Quantum number 4 (F).
    #         qn_5 (int): Quantum number 5.
    #         qn_6 (int): Quantum number 6.

    #     Returns:
    #         float: Degeneracy/weight calculated from quantum numbers of the given species.
    #     """
    #     if str(database_code).zfill(6) in ['004001', '004501']:
    #         if qn_2 % 2 == 0:
    #             weight = ((qn_1 * 2) + 1) * 1
    #         else:
    #             weight = ((qn_1 * 2) + 1) * 3
    #     elif str(database_code).zfill(6) in ['004101', '004201']:
    #         weight = ((qn_1 * 2) + 1) * 1
    #     else:
    #         if qn_6 == -1:
    #             if qn_5 == -1:
    #                 if qn_4 == -1:
    #                     if qn_3 == -1:
    #                         if qn_2 == -1:
    #                             if qn_1 == -1:
    #                                 raise ValueError(
    #                                     'No valid quantum numbers!')
    #                             else:
    #                                 weight = (qn_1 * 2) + 1
    #                         else:
    #                             weight = (qn_2 * 2) + 1
    #                     else:
    #                         weight = (qn_3 * 2) + 1
    #                 else:
    #                     weight = (qn_4 * 2) + 1
    #             else:
    #                 weight = (qn_5 * 2) + 1
    #         else:
    #             weight = (qn_6 * 2) + 1
    #     return weight

    @staticmethod
    def default_disk_density(position, inner_radius, outer_radius, ref_scale_height=10. * 149597870700.0,
                             ref_radius=100. * 149597870700.0, alpha=0.9, beta=0.8, tapered_gamma=None,
                             column_dens_exp=None, real_zero=True):
        """Shakura and Sunyaev disk density profile.

        Args:
            position (List[float, float, float]): Position in model space.
            inner_radius (float): Inner radius of the disk.
            outer_radius (float): Outer radius of the disk.
            ref_scale_height (float): Reference scale height.
            ref_radius (float): Reference radius.
            alpha (float): Exponent for radial density decrease.
            beta (float): Exponent for disk flaring.
            column_dens_exp (float): If set, calculate alpha from the surface density exponent.
                (Defined positively)
            real_zero (bool): No minimum value for the density.

        Returns:
            Float: Density at the given position.
        """
        #: float: Cylindrical radius
        radius_cy = np.sqrt(position[0] ** 2 + position[1] ** 2)
        if outer_radius >= radius_cy >= inner_radius:
            if column_dens_exp is not None:
                alpha = beta + column_dens_exp
            #: float: Vertical height
            vert_height = abs(position[2])
            #: float: Vertical scale height
            scale_height = ref_scale_height * (radius_cy / ref_radius) ** beta
            #: float: Shakura and Sunyaev density distribution
            density = (radius_cy / ref_radius) ** (-alpha) * \
                np.exp(-0.5 * (vert_height / scale_height) ** 2)
        else:
            density = 0.

        if tapered_gamma is not None:
            density *= np.exp(-(radius_cy / ref_radius)
                              ** (2 - tapered_gamma))
        if not real_zero:
            density = max(density, 1e-200)

        return density

    @staticmethod
    def default_disk_scale_height(radius, beta=0.8,
                                  ref_scale_height=10. * 149597870700.0, ref_radius=100. * 149597870700.0):
        """Shakura and Sunyaev disk density profile.

        Args:
            radius (float): Distance from the center in the midplane.
            ref_scale_height (float): Reference scale height.
            ref_radius (float): Reference radius.
            beta (float): Exponent for disk flaring.

        Returns:
            Float: Scale Height at the given position.
        """
        scale_height = ref_scale_height * (radius / ref_radius) ** beta
        return scale_height

    @staticmethod
    def const_sphere_density(position, outer_radius, inner_radius=None):
        """Density profile with a sphere of constant density.

        Args:
            position (List[float, float, float]): Position in model space.
            outer_radius (float): Outer radius of the sphere.
            inner_radius (float): Inner radius of the sphere.

        Returns:
            Float: Density at the given position.
        """
        #: float: Radial distance from center
        radius = np.sqrt(position[0] ** 2 + position[1]
                         ** 2 + position[2] ** 2)
        #: float: Constant density inside the Sphere
        density = 0.
        if inner_radius is None:
            if radius <= outer_radius:
                density = 1.
        else:
            if outer_radius >= radius >= inner_radius:
                density = 1.
        return density

    # @staticmethod
    # def bonor_ebert_density(position, outer_radius, truncation_radius, exponent=-2.):
    #     """Density profile with a Bonnor-Ebert sphere.

    #     Notes:
    #         Link: http://arxiv.org/abs/1401.5064

    #     Args:
    #         position (List[float, float, float]): Position in model space.
    #         outer_radius (float): Radius of the sphere.
    #         truncation_radius (float): Radius of inner regions of constant density.
    #         exponent (float): Exponential decrease of he density.

    #     Returns:
    #         Float: Density at the given position.
    #     """
    #     #: float: Radial distance
    #     radius = np.sqrt(position[0] ** 2 + position[1]
    #                      ** 2 + position[2] ** 2)
    #     #: float: Bonnor-Ebert sphere density distribution
    #     density = 0.0
    #     if radius <= truncation_radius:
    #         density = truncation_radius ** exponent
    #     elif radius <= outer_radius:
    #         density = radius ** exponent
    #     return density

    # @staticmethod
    # def random_density_distribution(position, d_exp=3):
    #     """Density profile with random clumps.

    #     Args:
    #         position (List[float, float, float]): Position in model space.
    #         outer_radius (float): Radius of the sphere.

    #     Returns:
    #         Float: Density at the given position.
    #     """
    #     density = 10 ** (np.random.random() * d_exp)
    #     return density

    @staticmethod
    def simple_mag_field(mag_field_strength, axis='z',
                         random_variations=False, rnd_b_min=0.):
        """Magnetic field pointing in one direction.

        Args:
            mag_field_strength (float): Amplitude of the magnetic field strength.
            axis (str): Axis name of the magnetic field direction.
                [x, y, z]
            random_variations (bool): Instead of a constant magnetic field strength,
                the field strength is randomly chosen between rnd_b_min and mag_field_strength.
            rnd_b_min (float): Minimum magnetic field strength for random_variations.

        Returns:
            List[float, float, float]: Magnetic field strength at any position.
        """
        #: List: Magnetic field strength
        mag = np.array([0., 0., 0.])
        if random_variations:
            mag_field_strength = rnd_b_min + \
                np.random.random() * (mag_field_strength - rnd_b_min)
        if axis == 'x':
            mag[0] += mag_field_strength
        elif axis == 'y':
            mag[1] += mag_field_strength
        elif axis == 'z':
            mag[2] += mag_field_strength
        else:
            raise ValueError(
                'Chosen axis direction for the magnetic field strength is not valid!')
        return mag

    @staticmethod
    def radial_mag_field(mag_field_strength, position):
        """Magnetic field pointing in one direction.

        Args:
            mag_field_strength (float): Amplitude of the magnetic field strength.
            position ([float, float, float]): Position in the grid

        Returns:
            List[float, float, float]: Magnetic field strength at any position.
        """
        #: List: Magnetic field strength
        if np.linalg.norm(position) == 0:
            return np.zeros(3)
        return np.ones(3) * mag_field_strength * position[:] / np.linalg.norm(position)

    @staticmethod
    def disturbed_mag_field(mag_field_strength, main_axis='z', rel_strength=0.1):
        """Magnetic field pointing in one direction, but with a small
        random disturbance in the perpendicular direction. The strength of the disturbance
        is randomized.

        Args:
            mag_field_strength (float): Amplitude of the magnetic field strength.
            main_axis (str): Axis name of the main (stronger) magnetic field direction.
            rel_strength ( float): Relative strength between secondary and main magnetic
                field component.

        Returns:
            List[float, float, float]: Magnetic field strength at any position.
        """
        #: list: Magnetic field strength
        mag = [0., 0., 0.]
        #: float: Disturbing magnetic field strength
        sec_mag_strength_1 = mag_field_strength * \
            rel_strength * (2. * np.random.random() - 1.)
        sec_mag_strength_2 = mag_field_strength * \
            rel_strength * (2. * np.random.random() - 1.)
        if main_axis == 'x':
            mag[0] += mag_field_strength
            mag[1] += sec_mag_strength_1
            mag[2] += sec_mag_strength_2
        elif main_axis == 'y':
            mag[0] += sec_mag_strength_1
            mag[1] += mag_field_strength
            mag[2] += sec_mag_strength_2
        elif main_axis == 'z':
            mag[0] += sec_mag_strength_1
            mag[1] += sec_mag_strength_2
            mag[2] += mag_field_strength
        else:
            raise ValueError(
                'Chosen main axis direction for the magnetic field strength is not valid!')
        return mag

    @staticmethod
    def disturbed_mag_field_2(mag_field_strength, main_axis='z', max_angle=20):
        """Magnetic field with small offset in their angle around a fixed direction.

        Args:
            mag_field_strength (float): Amplitude of the magnetic field strength.
            main_axis (str): Axis name of the main (stronger) magnetic field direction.
            max_angle ( float): Maximum angle of the disturbed magnetic field vector [deg].

        Returns:
            List[float, float, float]: Magnetic field strength at any position.
        """
        #: list: Magnetic field strength
        mag = [0., 0., 0.]
        #: float: Disturbed field angles
        disturbed_theta_angle = np.arccos(np.random.random()) * max_angle / 90.
        disturbed_phi_angle = np.random.random() * 2. * np.pi

        mag_strength_perp_1 = mag_field_strength * \
            np.cos(disturbed_phi_angle) * np.sin(disturbed_theta_angle)
        mag_strength_perp_2 = mag_field_strength * \
            np.sin(disturbed_phi_angle) * np.sin(disturbed_theta_angle)
        mag_strength_main = mag_field_strength * np.cos(disturbed_theta_angle)
        if main_axis == 'x':
            mag[0] += mag_strength_main
            mag[1] += mag_strength_perp_1
            mag[2] += mag_strength_perp_2
        elif main_axis == 'y':
            mag[0] += mag_strength_perp_1
            mag[1] += mag_strength_main
            mag[2] += mag_strength_perp_2
        elif main_axis == 'z':
            mag[0] += mag_strength_perp_1
            mag[1] += mag_strength_perp_2
            mag[2] += mag_strength_main
        else:
            raise ValueError(
                'Chosen main axis direction for the magnetic field strength is not valid!')
        return mag

    @staticmethod
    def two_simple_mag_field(position, mag_field_strength):
        """Magnetic field pointing in the z-direction in the first half of the
        model space and pointing in the y-direction in the other half.

        Args:
            position (List[float, float, float]): position in model space.
            mag_field_strength (float): Amplitude of the magnetic field strength.

        Returns:
            List[float, float, float]: Magnetic field strength at the given
            position.
        """
        if position[2] < 0:
            mag = [0, 0, mag_field_strength]
        else:
            mag = [0, mag_field_strength, 0]
        return mag

    def toroidal_mag_field(self, position, mag_field_strength):
        """Magnetic field pointing in toroidal (phi) direction.

        Notes:
            Link: https://en.wikipedia.org/wiki/Toroidal_and_poloidal

        Args:
            position (List[float, float, float]): position in model space.
            mag_field_strength (float): Amplitude of the magnetic field strength.

        Returns:
            List[float, float, float]: Magnetic field strength at the given
            position.
        """
        #: List(float, float, float): Spherical coordinates
        spherical_coord = self.cartesian_to_spherical(position)
        # Only a field component in the xy-plane
        mag = [mag_field_strength, mag_field_strength, 0]
        # Multiplication with the phi direction unit vectors
        mag[0] *= -np.sin(spherical_coord[2])
        mag[1] *= np.cos(spherical_coord[2])
        return mag

    def poloidal_mag_field(self, position, mag_field_strength, torus_r_distance):
        """Magnetic field pointing a the poloidal direction.

        Notes:
            Link: https://en.wikipedia.org/wiki/Toroidal_and_poloidal

        Args:
            position (List[float, float, float]): position in model space.
            mag_field_strength (float): Amplitude of the magnetic field strength.
            torus_r_distance (float): Radial distance of the centre of the poloidal torus

        Returns:
            List[float, float, float]: Magnetic field strength at the given
            position.
        """
        #: List(float, float, float): Spherical coordinates
        spherical_coord = self.cartesian_to_spherical(position)
        #: float: Cylindrical radius
        radius_cy = spherical_coord[0] * np.cos(spherical_coord[1])
        #: float: Theta angle related to the radial distance
        theta = np.arctan2(position[2], (radius_cy - torus_r_distance))
        mag = np.zeros(3)
        # Magnetic field is rotating around the radial ring
        mag[0] = -np.sin(theta) * -np.cos(spherical_coord[2])
        mag[1] = -np.sin(theta) * -np.sin(spherical_coord[2])
        mag[2] = np.cos(theta)
        mag *= mag_field_strength
        return mag

    def hourglass_mag_field(self, position, mag_field_strength, radius):
        """Hourglass magnetic field.

        Args:
            position (List[float, float, float]): position in model space.
            mag_field_strength (float): Amplitude of the magnetic field strength.
            radius (float): Radial extent of the model space.

        Returns:
            List[float, float, float]: Magnetic field strength at the given
            position.
        """
        #: List(float, float, float): Spherical coordinates
        spherical_coord = self.cartesian_to_spherical(position)
        #: float: Weighting factor
        gamma = 5
        #: float: Radial component of the magnetic field
        mag_r = mag_field_strength * \
            (gamma * radius ** 2 / (radius + spherical_coord[0]) ** 2)
        #: float: Z component of the magnetic field
        mag_z = mag_field_strength
        mag = np.zeros(3)
        # Conversion to cartesian coordinates
        mag[0] += mag_r * np.cos(spherical_coord[1]) * \
            np.cos(spherical_coord[2])
        mag[1] += mag_r * np.cos(spherical_coord[1]) * \
            np.sin(spherical_coord[2])
        mag[2] += mag_r * np.sin(spherical_coord[1])
        # The field should point into the same direction above and below the xy-plane
        if spherical_coord[1] < 0:
            mag *= -1
        mag[2] += mag_z
        return mag

    # @staticmethod
    # def logarithmic_bin_distribution(position, length, bins, min_exponent=0., max_exponent=3.):
    #     """Distribute values logarithmically over a number of bins.

    #     Args:
    #         position (float): 1 dimensional position.
    #         length (float): Length over which the bins are distributed.
    #         bins (int): Number of bins.
    #         min_exponent (float): Lowest exponent.
    #         max_exponent (float): Highest exponent.

    #     Returns:
    #         Float: Value related to the bin at the given position.
    #     """
    #     #: int: Bin id corresponding to the x-position
    #     pos_bin = int(position / length * bins)
    #     #: float: Exponent to scale the temperature logarithmically
    #     exponent = min_exponent - \
    #         (min_exponent - max_exponent) * (pos_bin / (bins - 1))
    #     value = 10 ** exponent
    #     return value

    # def crutcher_mag_field_strength(self, density, correction_factor=1):
    #     """Magnetic field strength from Crutcher et al. (2010)

    #     Args:
    #         density (float): Hydrogen mass density.
    #         correction_factor (float): Normalization factor from grid creation.

    #     Returns:
    #         List[float, float, float]: Magnetic field strength at any position.
    #     """
    #     #: Float: Hydrogen number density [m^-3]
    #     n = correction_factor * density / \
    #         (self.const['avg_gas_mass'] * self.const['amu'])
    #     #: Float: Reference hydrogen number density [m^-3]
    #     n_0 = 300e6
    #     #: Float: Reference magnetic field strength [T]
    #     b_0 = 10e-10
    #     #: Float: Exponential factor
    #     alpha = 0.65

    #     #: Float: Magnetic field strength from equation of Crutcher et al. (2010)
    #     mag_field_strength = b_0 * (n / n_0) ** alpha
    #     return mag_field_strength

    # def crutcher_mag_field(self, density, correction_factor=1, axis='z'):
    #     """Magnetic field pointing in one direction.

    #     Args:
    #         density (float): Hydrogen mass density.
    #         correction_factor (float): Normalization factor from grid creation.
    #         axis (str): Axis name of the magnetic field direction.

    #     Returns:
    #         List[float, float, float]: Magnetic field strength at any position.
    #     """
    #     #: Float: Magnetic field strength from equation of Crutcher et al. (2010)
    #     mag_field_strength = self.crutcher_mag_field_strength(
    #         density, correction_factor)

    #     #: List: Magnetic field strength
    #     mag = [0., 0., 0.]

    #     if axis == 'x':
    #         mag[0] += mag_field_strength
    #     elif axis == 'y':
    #         mag[1] += mag_field_strength
    #     elif axis == 'z':
    #         mag[2] += mag_field_strength
    #     else:
    #         raise ValueError(
    #             'Chosen axis direction for the magnetic field strength is not valid!')
    #     return mag

    # @staticmethod
    # def dustem_size_distributions(dist_type, nr_of_dust_species, a_min, a_max, dist_parameter):
    #     """Calculates various grain size distributions.
    #         From DustEM Code "https://www.ias.u-psud.fr/DUSTEM/"

    #     Args:
    #         dist_type (str): Type of the size distribution (from GRAIN.DAT).
    #             (LOGN -> log normal or PLAW -> power law)
    #         nr_of_dust_species (int): Number of dust grain sizes.
    #         a_min (int): Minimum dust grain size. [m]
    #         a_max (int): Maximum dust grain size. [m]
    #         dist_parameter (list): List of size distribution parameters (from GRAIN.DAT).
    #     """
    #     # Convert to cgs
    #     a_min *= 1e2
    #     a_max *= 1e2

    #     if not ('PLAW' in dist_type.upper() or 'LOGN' in dist_type.upper()):
    #         raise ValueError(
    #             'Dust grains size distribution dist_type not known!')

    #     if nr_of_dust_species != 1:
    #         da = (np.log(a_max) - np.log(a_min)) / \
    #             float(nr_of_dust_species - 1)
    #     else:
    #         da = np.log(a_max) - np.log(a_min)
    #     ava = np.zeros(nr_of_dust_species)
    #     size_ava = ava.copy()

    #     if 'PLAW' in dist_type.upper():
    #         if ('-ED' in dist_type.upper() and '-CV' in dist_type.upper() and len(dist_parameter) != 7) or \
    #                 (('-ED' in dist_type.upper() or '-CV' in dist_type.upper()) and len(dist_parameter) != 4):
    #             raise ValueError(
    #                 'Wrong number of parameters for size distribution!')
    #         for i_dust in range(nr_of_dust_species):
    #             aux = np.log(a_min) + float(i_dust) * da
    #             # argu = (4.0 + dist_parameter[0]) * aux
    #             argu = (1.0 + dist_parameter[0]) * aux
    #             if argu > -350.0:
    #                 ava[i_dust] = np.exp(argu)
    #             size_ava[i_dust] = np.exp(aux)
    #         if '-ED' in dist_type.upper():
    #             for i_dust in range(nr_of_dust_species):
    #                 if size_ava[i_dust] >= dist_parameter[1]:
    #                     ava[i_dust] *= np.exp(
    #                         -((size_ava[i_dust] - (dist_parameter[1])) /
    #                           (dist_parameter[2])) ** dist_parameter[3])
    #         if '-CV' in dist_type.upper():
    #             if '-ED' in dist_type.upper():
    #                 au = dist_parameter[4]
    #                 zeta = abs(dist_parameter[5])
    #                 zxp = np.sign(1.0, dist_parameter[5])
    #                 gama = dist_parameter[6]
    #             else:
    #                 au = dist_parameter[1]
    #                 zeta = abs(dist_parameter[2])
    #                 zxp = np.sign(1.0, dist_parameter[2])
    #                 gama = dist_parameter[3]
    #             for i_dust in range(nr_of_dust_species):
    #                 ava[i_dust] *= (1.0 + zeta *
    #                                 (size_ava[i_dust] / au) ** gama) ** zxp
    #     elif 'LOGN' in dist_type.upper():
    #         if dist_parameter[0] == 0. or dist_parameter[1] == 0.:
    #             raise ValueError(
    #                 'Centroid or sigma of log-normal cannot be 0!')
    #         if nr_of_dust_species != 1:
    #             da = (np.log(a_max) - np.log(a_min)) / \
    #                 float(nr_of_dust_species - 1)
    #         else:
    #             da = np.log(a_max) - np.log(a_min)
    #         for i_dust in range(nr_of_dust_species):
    #             aux = np.log(a_min) + float(i_dust) * da
    #             # argu = 3.0 * aux - 0.5 * ((aux - np.log(dist_parameter[0])) / dist_parameter[1]) ** 2
    #             argu = - 0.5 * \
    #                 ((aux - np.log(dist_parameter[0])
    #                   ) / dist_parameter[1]) ** 2
    #             if argu > -350.0:
    #                 ava[i_dust] = np.exp(argu)
    #             size_ava[i_dust] = np.exp(aux)
    #     # Convert to SI
    #     size_ava *= 1e-2
    #     ava *= 1e2
    #     return size_ava, ava

    # @staticmethod
    # def trust_size_distributions(aeff, a1, a2, a3, a4, b0, b1, b2, b3, b4, c0, m1, m2, m3, m4):
    #     log_g = c0 + b0 * np.log10(aeff)
    #     if b1 != 0 or a1 != 0 or m1 != 0:
    #         log_g -= b1 * abs(np.log10(aeff / a1)) ** m1
    #     if b2 != 0 or a2 != 0 or m2 != 0:
    #         log_g -= b2 * abs(np.log10(aeff / a2)) ** m2
    #     if b3 != 0 or a3 != 0 or m3 != 0:
    #         log_g -= b3 * abs(aeff - a3) ** m3
    #     if b4 != 0 or a4 != 0 or m4 != 0:
    #         log_g -= b4 * abs(aeff - a4) ** m4
    #     return 10**log_g
