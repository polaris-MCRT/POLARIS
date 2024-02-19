# -*- coding: utf-8 -*-

import numpy as np


class Math:
    """Constants and math functions.
    """

    def __init__(self):
        """Initialisation of all constants and conversion factors.
        """

        #: dict: Constants taken from astropy (reference: CODATA 2014, IAU 2012 Resolution B2)
        self.const = {
            'M_sun':        1.9884754153381438e+30,  # Solar mass [kg]
            'M_jup':        1.8981871658715508e+27,  # Jupiter mass [kg]
            'R_sun':        695700000.0,             # Nominal solar radius [m]
            'R_jup':        71492000.0,              # Nominal Jupiter equatorial radius [m]
            'L_sun':        3.828e+26,               # Nominal solar luminosity [W]
            'au':           149597870700.0,          # Astronomical Unit [m]
            'pc':           3.0856775814671916e+16,  # Parsec [m]
            'u':            1.66053904e-27,          # Atomic mass unit [kg]
            'G':            6.67408e-11,             # Gravitational constant [m^3 / (kg * s^2)]
            'h':            6.62607004e-34,          # Planck constant [J * s]
            'hbar':         1.0545718e-34,           # Reduced Planck constant [J * s]
            'c':            299792458.0,             # Speed of light in vacuum [m / s]
            'e':            1.6021766208e-19,        # Electron charge [C]
            'b_wien':       0.00289777196,           # Wien wavelength displacement law constant [m * K]
            'muB':          9.274009994e-24,         # Bohr magneton [J / T]
            'k_B':          1.38064852e-23,          # Boltzmann constant [J / K]
            'N_A':          6.02214129e+23,          # Avogadro’s number [1 / mol]
            'R':            8.31446262,              # Gas constant [J / (K * mol)]
            'Ryd':          10973731.568508,         # Rydberg constant [1 / m]
            'sigma_sb':     5.670367e-08,            # Stefan-Boltzmann constant [W / (m^2 * K^4)]
            'm_e':          9.10938356e-31,          # Electron mass [kg]
            'm_p':          1.672621898e-27,         # Proton mass [kg]
            'eps0':         8.854187817620389e-12,   # Vacuum permittivity [F / m]
            'avg_gas_mass': 2.,                      # Average atomic mass unit per gas particle
        }

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
                raise ValueError('The distance is not set. Without the distance, no conversion to arcseconds is available!')
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
                raise ValueError('The distance is not set. Without the distance, no conversion to arcseconds is available!')
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
                raise ValueError('There is no unit type available for conversion!')
        elif unit_type == 'velocity':
            parsed_value = self.parse_velocity(value)
            if unit == 'km/s':
                return parsed_value * 1e-3
            elif unit == 'km/h':
                return parsed_value * 3600 / 1e3
            elif unit == 'm/s' or unit is None:
                return parsed_value
            else:
                raise ValueError('There is no unit type available for conversion!')
        elif unit_type == 'mass':
            parsed_value = self.parse_mass(value)
            if unit == 'm_jup' or unit == 'mjup':
                return parsed_value / self.const['M_jup']
            elif unit == 'm_sun' or unit == 'msun':
                return parsed_value / self.const['M_sun']
            elif unit == 'kg' or unit is None:
                return parsed_value
            else:
                raise ValueError('There is no unit type available for conversion!')
        elif unit_type == 'luminosity':
            parsed_value = self.parse_luminosity(value)
            if unit == 'l_sun' or unit == 'lsun':
                return parsed_value / self.const['L_sun']
            elif unit == 'w' or unit is None:
                return parsed_value
            else:
                raise ValueError('There is no unit type available for conversion!')
        elif unit_type == 'angle':
            parsed_value = self.parse_angle(value)
            if unit == 'arcsec' or unit == 'arc_sec':
                return parsed_value * 3600
            elif unit == 'rad':
                return parsed_value / 180. * np.pi
            elif unit == 'degree' or unit is None:
                return parsed_value
            else:
                raise ValueError('There is no unit type available for conversion!')
        else:
            raise ValueError('There is no unit type available for conversion!')

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
        elif 'r_jup' in length.lower() or 'rjup' in length.lower():
            conv = self.const['R_jup']
            length = length.lower().replace('r_jup', '').replace('rjup', '')
        elif 'r_sun' in length.lower() or 'rsun' in length.lower():
            conv = self.const['R_sun']
            length = length.lower().replace('r_sun', '').replace('rsun', '')
        try:
            return float(length) * conv
        except ValueError:
            print(f'Length string {length} cannot be interpreted!')

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
            print(f'Velocity string {velocity} cannot be interpreted!')

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
            print(f'Mass string {mass} cannot be interpreted!')

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
            print(f'Angle string {angle} cannot be interpreted!')

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
            print(f'Luminosity string {luminosity} cannot be interpreted!')

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
        cylindrical_coord[1] = np.arctan2(cartesian_coord[1], cartesian_coord[0])
        if cylindrical_coord[1] < 0.0:
            cylindrical_coord[1] += 2.0 * np.pi
        cylindrical_coord[2] = cartesian_coord[2]
        return cylindrical_coord

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

    @staticmethod
    def default_disk_density(position, inner_radius, outer_radius, ref_scale_height=10. * 149597870700.0,
                             ref_radius=100. * 149597870700.0, alpha=1.8, beta=1.1, tapered_gamma=None,
                             column_dens_exp=None, real_zero=True):
        """Power-law or viscous accretion disk model (exponential taper).
        Andrews et al. (2009), ApJ 700, 1502
        Kwon et al. (2015), ApJ 808, 102

        Args:
            position (List[float, float, float]): Position in model space.
            inner_radius (float): Inner radius of the disk.
            outer_radius (float): Outer radius of the disk.
            ref_scale_height (float): Reference scale height.
            ref_radius (float): Reference radius.
            alpha (float): Exponent for radial density decrease.
            beta (float): Exponent for disk flaring.
            tapered_gamma (float): If set, add an exponential taper at large radii characterized by tapered_gamma.
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
            density *= np.exp(-(radius_cy / ref_radius) ** (2 - tapered_gamma))
        if not real_zero:
            density = max(density, 1e-200)

        return density

    @staticmethod
    def default_disk_scale_height(radius, beta=1.1,
                                  ref_scale_height=10. * 149597870700.0, ref_radius=100. * 149597870700.0):
        """Width of the Gaussian vertical density profile.

        Args:
            radius (float): Distance from the center in the midplane.
            beta (float): Exponent for disk flaring.
            ref_scale_height (float): Reference scale height.
            ref_radius (float): Reference radius.

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
            raise ValueError('Chosen axis direction for the magnetic field strength is not valid!')
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
            raise ValueError('Chosen main axis direction for the magnetic field strength is not valid!')
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
            raise ValueError('Chosen main axis direction for the magnetic field strength is not valid!')
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
