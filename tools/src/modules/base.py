"""base.py includes the base versions of each Polaris component to be loaded by other scripts.
"""

import numpy as np
from modules.math import Math


# ----- Background source class -----
class BGSource:
    """The BGSource class is the base version for each background sources.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the background source parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """

        self.file_io = file_io
        self.parse_args = parse_args

        #: dict: Includes parameters of a specific background source
        self.parameter = {
            'nr_photons': 1e5,
            'A': 0.,
            'T': 0.,
            'Q': 0.,
            'U': 0.,
            'V': 0.,
            'theta': 0.,
            'phi': 0.,
        }

    def get_command_line(self):
        """Provides background source command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the background source.
        """

        if np.sqrt(self.parameter['Q'] ** 2 + self.parameter['U'] ** 2) > 1. or self.parameter['V'] > 1.:
            raise ValueError('Error: Polarization degree is larger than 100%!')

        return '\t<source_background nr_photons = "' \
               + str(self.parameter['nr_photons']) + '">\t' \
               + str(self.parameter['A']) + '\t' \
               + str(self.parameter['T']) + '\t' \
               + str(self.parameter['Q']) + '\t' \
               + str(self.parameter['U']) + '\t' \
               + str(self.parameter['V']) + '\t' \
               + str(self.parameter['theta']) + '\t' \
               + str(self.parameter['phi']) + '\n'

    def get_command(self):
        """Provides background source command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the background source.
        """
        return self.get_command_line()


# ----- Detector class -----
class Detector:
    """This is the default class for detector configurations.
    """

    def __init__(self, model, parse_args):
        """Initialisation of all usable options.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
        """

        self.model = model
        self.parse_args = parse_args

        # Get math module
        from modules.math import Math
        self.math = Math()

        #: dict: Includes parameters of a specific detector configuration
        self.parameter = {
            'shape': 'cartesian',
            'wavelength_min': 1e-6,
            'wavelength_max': 1e-6,
            'nr_of_wavelength': 1,
            'wavelength_list': None,
            'distance': model.parameter['distance'],
            'nr_pixel_x': 201,
            'nr_pixel_y': 201,
            'max_subpixel_lvl': 0,
            'source_id': 1,
            'gas_species_id': 1,
            'transition_id': 1,
            'max_velocity': 3000,
            'nr_velocity_channels': 35,
            'sidelength_zoom_x': 1,
            'sidelength_zoom_y': 1,
            'map_shift_x': 0.,
            'map_shift_y': 0.,
            'acceptance_angle': 1.0,
            'rot_angle_1': 0.,
            'rot_angle_2': 0.,
            'rot_axis_1': [1, 0, 0],
            'rot_axis_2': [0, 1, 0],
            # Parameter for healpix all-sky maps
            'nr_sides': 32,
            'obs_position_x': 0.,
            'obs_position_y': 0.,
            'obs_position_z': 0.,
            'obs_velocity_x': None,
            'obs_velocity_y': None,
            'obs_velocity_z': None,
            'all_sky_l_min': None,
            'all_sky_l_max': None,
            'all_sky_b_min': None,
            'all_sky_b_max': None,
        }

    def get_zoom_cmd(self):
        """Provides the command strings for zooming the raytracing map.
        """
        cmd_string = ''
        # Add zoom factor if chosen
        for direction in ['_x', '_y']:
            zoom_string = 'sidelength_zoom' + direction
            if self.parameter[zoom_string] != 1.:
                sidelength = 2. * \
                    self.model.tmp_parameter['radius' +
                                             direction + '_m'] / self.parameter[zoom_string]
                cmd_string += '\t' + str(sidelength)
            elif self.parameter['sidelength_zoom_x'] != 1. or self.parameter['sidelength_zoom_y'] != 1.:
                cmd_string += '\t' + \
                    str(2. *
                        self.model.tmp_parameter['radius' + direction + '_m'])
            elif self.parameter['map_shift_x'] != 0 or self.parameter['map_shift_y'] != 0:
                cmd_string += '\t' + \
                    str(2. *
                        self.model.tmp_parameter['radius' + direction + '_m'])
        return cmd_string

    def get_shift_cmd(self):
        """Provides the command strings for shifting the raytracing map.
        """
        cmd_string = ''
        # Add detector map shift if chosen
        for direction in ['_x', '_y']:
            shift_string = 'map_shift' + direction
            if self.parameter['map_shift_x'] != 0 or self.parameter['map_shift_y'] != 0:
                cmd_string += '\t' + str(self.parameter[shift_string])
        return cmd_string

    def get_nr_pixel_cmd(self):
        cmd_string = str(int(self.parameter['nr_pixel_x'])) + \
            '*' + str(int(self.parameter['nr_pixel_y']))
        return cmd_string

    def get_dust_scattering_command_line(self):
        """Provides detector configuration command line for Monte-Carlo
        simulations for POLARIS .cmd file.

        Returns:
            str: Command line to consider the detector configuration.
        """
        param_name_list = ['wavelength_min', 'wavelength_max', 'nr_of_wavelength',
                           'rot_angle_1', 'rot_angle_2', 'distance']
        cmd_string = '\t<detector_dust_mc nr_pixel = "' + self.get_nr_pixel_cmd() + \
            '">'
        for param_name in param_name_list:
            cmd_string += '\t' + str(self.parameter[param_name])
        # Zooming but not shifting is possible with the MC detector
        cmd_string += self.get_zoom_cmd() + '\n'
        return cmd_string

    def get_dust_scattering_command(self):
        """Provides detector configuration command line for Monte-Carlo
        simulations for POLARIS .cmd file.

        Returns:
            str: Command line to consider the detector configuration.
        """
        if self.parse_args.wavelength is None and self.parameter['wavelength_list'] is not None:
            new_command_line = str()
            for wl in self.parameter['wavelength_list']:
                self.parameter['wavelength_min'] = wl
                self.parameter['wavelength_max'] = wl
                new_command_line += self.get_dust_scattering_command_line()
            return new_command_line
        return self.get_dust_scattering_command_line()

    def get_dust_emission_command_line(self):
        """Provides detector configuration command line for raytrace
        simulations for POLARIS .cmd file.

        Returns:
            str: Command line to consider the detector configuration.
        """
        param_name_list = ['wavelength_min',
                           'wavelength_max', 'nr_of_wavelength', 'source_id']
        if self.parameter['shape'] == 'healpix' and self.parameter['obs_position_x'] is not None and \
                self.parameter['obs_position_y'] is not None and self.parameter['obs_position_z'] is not None:
            cmd_string = '\t<detector_dust_healpix nr_sides = "' + \
                str(int(self.parameter['nr_sides'])) + '">'
            param_name_list.extend(
                ['obs_position_x', 'obs_position_y', 'obs_position_z'])
            if self.parameter['all_sky_l_min'] is not None and self.parameter['all_sky_l_max'] is not None and \
                    self.parameter['all_sky_b_min'] is not None and self.parameter['all_sky_b_max'] is not None:
                param_name_list.extend(
                    ['all_sky_l_min', 'all_sky_l_max', 'all_sky_b_min', 'all_sky_b_max'])
        elif self.parameter['shape'] == 'polar':
            cmd_string = '\t<detector_dust_polar nr_pixel = "' + self.get_nr_pixel_cmd() + \
                '">'
            param_name_list.extend(['rot_angle_1', 'rot_angle_2', 'distance'])
        elif self.parameter['shape'] == 'slice':
            cmd_string = '\t<detector_dust_slice nr_pixel = "' + self.get_nr_pixel_cmd() + \
                '">'
            param_name_list.extend(['rot_angle_1', 'rot_angle_2', 'distance'])
        else:
            cmd_string = '\t<detector_dust nr_pixel = "' + self.get_nr_pixel_cmd() + '">'
            param_name_list.extend(['rot_angle_1', 'rot_angle_2', 'distance'])
        # Create cmd string for detector
        for param_name in param_name_list:
            cmd_string += '\t' + str(self.parameter[param_name])
        cmd_string += self.get_zoom_cmd() + self.get_shift_cmd() + '\n'
        return cmd_string

    def get_dust_emission_command(self):
        """Provides detector configuration command line for raytrace
        simulations for POLARIS .cmd file.

        Returns:
            str: Command line to consider the detector configuration.
        """
        if self.parse_args.wavelength is None and self.parameter['wavelength_list'] is not None:
            new_command_line = str()
            for wl in self.parameter['wavelength_list']:
                self.parameter['wavelength_min'] = wl
                self.parameter['wavelength_max'] = wl
                new_command_line += self.get_dust_emission_command_line()
            return new_command_line
        return self.get_dust_emission_command_line()

    def get_line_command_line(self):
        """Provides detector configuration command line for spectral line
        simulations for POLARIS .cmd file.

        Returns:
            str: Command line to consider the detector configuration.
        """
        param_name_list = ['gas_species_id',
                           'transition_id', 'source_id', 'max_velocity']
        if self.parameter['shape'] == 'healpix' and self.parameter['obs_position_x'] is not None and \
                self.parameter['obs_position_y'] is not None and self.parameter['obs_position_z'] is not None:
            cmd_string = '\t<detector_line_healpix nr_sides = "' + str(int(self.parameter['nr_sides'])) + \
                '" vel_channels = "' + \
                str(int(self.parameter['nr_velocity_channels'])) + '">'
            param_name_list.extend(
                ['obs_position_x', 'obs_position_y', 'obs_position_z'])
            if self.parameter['all_sky_l_min'] is not None and self.parameter['all_sky_l_max'] is not None and \
                    self.parameter['all_sky_b_min'] is not None and self.parameter['all_sky_b_max'] is not None:
                param_name_list.extend(
                    ['all_sky_l_min', 'all_sky_l_max', 'all_sky_b_min', 'all_sky_b_max'])
                if self.parameter['obs_velocity_x'] is not None and self.parameter['obs_velocity_y'] is not None and \
                        self.parameter['obs_velocity_z'] is not None:
                    param_name_list.extend(
                        ['obs_velocity_x', 'obs_velocity_y', 'obs_velocity_z'])
        elif self.parameter['shape'] == 'polar':
            cmd_string = '\t<detector_line_polar nr_pixel = "' + self.get_nr_pixel_cmd() + \
                '" vel_channels = "' + \
                str(int(self.parameter['nr_velocity_channels'])) + '">'
            param_name_list.extend(['rot_angle_1', 'rot_angle_2', 'distance'])
        elif self.parameter['shape'] == 'slice':
            cmd_string = '\t<detector_line_slice nr_pixel = "' + self.get_nr_pixel_cmd() + \
                '" vel_channels = "' + \
                str(int(self.parameter['nr_velocity_channels'])) + '">'
            param_name_list.extend(['rot_angle_1', 'rot_angle_2', 'distance'])
        else:
            cmd_string = '\t<detector_line nr_pixel = "' + self.get_nr_pixel_cmd() + \
                '" vel_channels = "' + \
                str(int(self.parameter['nr_velocity_channels'])) + '">'
            param_name_list.extend(['rot_angle_1', 'rot_angle_2', 'distance'])
        # Create cmd string for detector
        for param_name in param_name_list:
            cmd_string += '\t' + str(self.parameter[param_name])
        cmd_string += self.get_zoom_cmd() + self.get_shift_cmd() + '\n'
        return cmd_string

    def get_line_command(self):
        """Provides detector configuration command line for spectral line
        simulations for POLARIS .cmd file.

        Returns:
            str: Command line to consider the detector configuration.
        """
        return self.get_line_command_line()


# ----- Dust class -----
class Dust:
    """The Dust class is the base version for each dust composition.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        self.file_io = file_io
        self.parse_args = parse_args

        # Init dust_chooser to mix dust components
        from modules.dust import DustChooser
        self.dust_chooser = DustChooser(file_io, parse_args)

        #: dict: Includes parameters related to a specific dust composition
        self.parameter = {
            'dust_cat_file': None,
            'fraction': 1.0,
            'amin': None,
            'amax': None,
            'size_keyword': 'plaw',
            'size_parameter': [-3.5],
            'choice_id': None,
            'scattering': 'ISO',
            # Parameter to create dust catalogs
            'input_file': None,
            'input_type': None,  # 'dustem' or 'trust'
            'dustem_wl_file': 'LAMBDA',  # Only if input_type == dustem
            'wavelength_list': None,
            'size_list': None,
            'theta_angle_list': np.linspace(0, 180, 181),
            'inc_angle_list': np.zeros(1),
            'phi_angle_list': np.zeros(1),
            'material_density': 0,
            'calorimetry_type': 'heat_capacity',
            'subl_temp': 9999,
            'align': 0,
            'aspect_ratio': 1,
            'rat_delta': 0,
        }

    def create_catalog(self):
        """Create a dust catalog file for POLARIS.
        """
        # Get dust creator module
        from modules.create_dust import DustCreator
        self.dust_creator = DustCreator(
            self.file_io, self.parse_args, self.parameter)

        # Check if filename is defined
        if self.parameter['input_file'] is None:
            raise ValueError('No input file for dust catalog defined!')

        if self.parameter['input_type'] == 'dustem':
            self.dust_creator.create_polaris_dust_from_dustem()
        elif self.parameter['input_type'] == 'trust':
            self.dust_creator.create_polaris_dust_from_trust()

    def print_info(self, microns=False):
        # Define some colors to show error and ok more clearly
        color_blue = '\033[0;34m'
        color_red = '\033[0;31m'
        color_std = '\033[0m'
        color_bold = '\033[1m'
        # If micron is True convert to wl and size to microns
        if microns:
            unit = 'microns'
            mult = 1e6
        else:
            unit = 'm'
            mult = 1
        # Get dust data from catalog file
        size_list, wavelength_list = self.file_io.read_dust_file(
            self.parameter)
        print(color_bold + 'index\t' +
              'dust grain size [' + unit + ']' + color_std)
        for i, size in enumerate(size_list):
            print(color_bold + str(i + 1) + '\t' +
                  color_blue + str(size * mult) + color_std)
        print(color_bold + 'index\t' + 'wavelength [' + unit + ']' + color_std)
        for i, wl in enumerate(wavelength_list):
            print(color_bold + str(i + 1) + '\t' +
                  color_red + str(wl * mult) + color_std)

    def get_command_line(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the dust composition.
        """
        if self.parameter['choice_id'] is not None:
            dust_string = '\t<dust_component id = \"' + \
                str(self.parameter['choice_id']) + '\">\t"'
        else:
            dust_string = '\t<dust_component>\t"'
        if self.parameter['amin'] is not None and self.parameter['amax'] is not None:
            dust_string += self.file_io.path['dust'] + self.parameter['dust_cat_file'] +  \
                '" "' + self.parameter['size_keyword'] + '" ' + str(self.parameter['fraction']) + ' ' + \
                str(self.parameter['material_density']) + ' ' + str(self.parameter['amin']) + ' ' + \
                str(self.parameter['amax'])
            for i in range(len(self.parameter['size_parameter'])):
                dust_string += ' ' + str(self.parameter['size_parameter'][i])
        else:
            raise ValueError(
                'Minimum AND maximum dust grain size need to be defined!')
        return dust_string + '\n'

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the dust composition.
        """
        return self.get_command_line()


# ----- No dust class -----
class NoDust(Dust):
    """Dust class for graphite grains.
        """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Dust.__init__(self, file_io, parse_args)

    def get_command(self):
        """Provides dust composition command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the dust composition.
        """
        return ""


# ----- Gas class -----
class Gas:
    """The Gas class is the base version for each gas species.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
        """
        self.file_io = file_io
        self.parse_args = parse_args

        # Get math module
        self.math = Math()

        #: dict: Includes parameters related to a specific gas species
        self.parameter = {
            'filename': None,
            'level_pop_type': 1,
            'abundance': 1e-5,
            'max_velocity': 3000,
            'zeeman_usable': False,
        }

    def get_command_line(self):
        """Provides gas species command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the gas species.
        """
        gas_species_string = '\t<gas_species>\t"' + str(self.file_io.path['gas']) + str(self.parameter['filename']) + \
            '"\t' + str(self.parameter['level_pop_type']) + \
            '\t' + str(self.parameter['abundance'])
        if self.parse_args.simulation_type == 'zeeman':
            if self.parameter['zeeman_usable'] is True:
                gas_species_string += '\t"' + self.file_io.path['gas'] + \
                    self.parameter['filename'].replace(
                        '.dat', '_zeeman.dat') + '"'
            else:
                raise ValueError(
                    'Chosen gas species does not support Zeeman splitting!')
        return gas_species_string + '\n'

    def get_command(self):
        """Provides gas species command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the gas species.
        """
        return self.get_command_line()

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
        if self.parameter['zeeman_usable']:
            raise ValueError(
                'Please define the \"shift_2_mag\" routine for the chosen species!')
        else:
            raise ValueError(
                'The chosen species has no Zeeman splitting information available!')

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
        if self.parameter['zeeman_usable']:
            raise ValueError(
                'Please define the \"mag_2_shift\" routine for the chosen species!')
        else:
            raise ValueError(
                'The chosen species has no Zeeman splitting information available!')


# ----- No gas class -----
class NoGas(Gas):
    """Gas class if no gas species is used
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the gas parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        Gas.__init__(self, file_io, parse_args)

    def get_command(self):
        """Provides gas species command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the gas species.
        """
        return ""


# ----- Model class -----
class Model:
    """The Model class is the base version for each model.
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        # Get math module
        self.math = Math()

        #: dict: Includes parameters of a specific model
        self.parameter = {
            'distance': 1. * self.math.const['pc'],
            'gas_mass': 1. * self.math.const['M_sun'],
            'dust_mass': None,
            # Global extent variables to set all grid types at once
            'inner_radius': 1. * self.math.const['au'],
            'outer_radius': 100. * self.math.const['au'],
            'mass_fraction': 0.01,
            'grid_type': 'octree',
            # Global parameters of gas phase
            'turbulent_velocity': 100.,
            # Define which other objects will be connected to this model
            'background_source': 'bg_plane',
            'detector': None,
            'radiation_source': None,
            'dust_composition': None,
            'gas_species': None,
            'external_input_name': None,
            'vel_is_speed_of_sound': False,
            'enforced_scattering': True,
            'peel_off': True,
            # Plot parameter
            'midplane_points': 256,
            'midplane_zoom': 1,
        }

        #: dict: Includes parameters for the octree grid
        self.octree_parameter = {
            'sidelength': None,
            'max_tree_level': 5,
        }

        #: dict: Includes parameters for the spherical grid
        self.spherical_parameter = {
            'inner_radius': None,
            'outer_radius': None,
            'n_r': 100,
            'n_th': 101,
            'n_ph': 1,
            'sf_r': 1.03,
            'sf_ph': 1.0,
            'sf_th': -1.0,
            # These list are used as cell borders if sf_r or sf_th is zero
            'radius_list': [],
            'phi_list': [],
            'theta_list': [],
            # Split the first radial cell into multiple
            'split_first_cell': 1
        }

        #: dict: Includes parameters for the cylindrical grid
        self.cylindrical_parameter = {
            'inner_radius': None,
            'outer_radius': None,
            'z_max': None,
            'n_r': 100,
            'n_ph': 101,
            'n_z': 5,
            'sf_r': 1.03,
            'sf_ph': 1.0,
            'sf_z': -1.0,
            # This list is used as cell borders if sf_r or sf_z is zero
            'radius_list': [],
            'phi_list': [],
            'z_list': [],
            # Split the first radial cell into multiple
            'split_first_cell': 1
        }

        #: dict: Includes conversion factors of different quantities
        self.conv_parameter = {
            'conv_dens': 1.,
            'conv_len': 1.,
            'conv_mag': 1.,
            'conv_vel': 1.,
        }

        #: dict: Parameters that are used internally to modify the model class
        self.tmp_parameter = {
            'relative_gas_densities': None,
            'relative_dust_densities': None,
            'ignored_gas_density': 0,
            'ignored_dust_density': 0,
        }

        self.position = None
        self.volume = None

    def init_position(self, node, cell_IDs=None):
        """Initialise the grid position to calculate the necessary cell data.

        Args:
            node (node class instance)): instance of node.
            cell_IDs (List): cell_IDs of the cells (alternative to the grid position).
                Spherical -> [i_r, i_t, i_p]
                Cylindrical -> [i_r, i_p, i_z]
        """
        self.position = node.parameter['position']
        self.volume = node.parameter['volume']

    def get_gas_temperature(self):
        """The gas temperature can be modified by the code here if neccessary.

        Returns:
            float: Gas temperature at a given position.
        """
        return self.gas_temperature()

    def get_dust_temperature(self):
        """The dust temperature can be modified by the code here if neccessary.

        Returns:
            float: Dust temperature at a given position.
        """
        return self.dust_temperature()

    def get_gas_density_distribution(self):
        """Modifies the mass density defined by get_gas_density_distribution() to match the preset gas masses.

        Returns:
            float: Gas density at a given position.
        """
        if self.tmp_parameter['relative_gas_densities'] is not None:
            return np.multiply(self.gas_density_distribution(), self.tmp_parameter['relative_gas_densities'])
        return self.gas_density_distribution()

    def get_dust_density_distribution(self):
        """Modifies the mass density defined by get_dust_density_distribution() to match the preset dust masses.

        Returns:
            float: Dust density at a given position.
        """
        if self.tmp_parameter['relative_dust_densities'] is not None:
            return np.multiply(self.dust_density_distribution(), self.tmp_parameter['relative_dust_densities'])
        return self.dust_density_distribution()

    def get_velocity_field(self):
        """The velocity field can be modified by the code here if necessary.

        Returns:
            List[float, float, float]: Velocity at a given position.
        """
        return self.velocity_field()

    def get_magnetic_field(self):
        """The magnetic field can be modified by the code here if neccessary.

        Returns:
            List[float, float, float]: Magnetic field strength at a given position.
        """
        return self.magnetic_field()

    def get_dust_id(self):
        """The dust choice id can be modified by the code here if neccessary.

        Returns:
            int: dust ID.
        """
        return self.dust_id()

    def get_dust_min_size(self):
        """The minimum dust grain size can be modified by the code here if neccessary.

        Returns:
            float: minimum grain size
        """
        return self.dust_min_size()

    def get_dust_max_size(self):
        """The maximum dust grain size can be modified by the code here if neccessary.

        Returns:
            float: minimum grain size
        """
        return self.dust_max_size()

    def get_dz(self, radius):
        """Calculates the width of each vertical cell border depending on the radial position.

        Args:
            radius (float) : Cylindrical radius of current position

        Returns:
            float: Width between two cell borders.
        """
        return 10. * self.scale_height(radius) / self.cylindrical_parameter['n_z']

    def adjust_extent(self, sidelength_x, sidelength_y):
        """Adjust the extent of the model.

        Args:
            sidelength_x (float): New x-axis extent.
            sidelength_y (float): New y-axis extent.
        """
        factor_x = (0.5 * sidelength_x) / self.tmp_parameter['radius_x_m']
        factor_y = (0.5 * sidelength_y) / self.tmp_parameter['radius_y_m']
        self.tmp_parameter['radius_x_m'] *= factor_x
        self.tmp_parameter['radius_y_m'] *= factor_y
        self.tmp_parameter['radius_x_arcsec'] *= factor_x
        self.tmp_parameter['radius_y_arcsec'] *= factor_y
        self.tmp_parameter['radius_x_au'] *= factor_x
        self.tmp_parameter['radius_y_au'] *= factor_y
        self.tmp_parameter['radius_x_pc'] *= factor_x
        self.tmp_parameter['radius_y_pc'] *= factor_y

    def adjust(self, value, name):
        """Adjust a parameter of the model.

        Args:
            value (float): New value.
            name (str): Name of model parameter.
        """
        self.parameter[name] = value

    def ignore_cell(self, node=None):
        """Ignore a cell for grid refinement, if necessary for a given model.

        Args:
            node: Current node to check if it should be ignored for grid refinement

        Returns:
            bool: True if ignore, False if considered for refinement.
        """
        return False

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters and update 
        model parameter that depend on other parameter.
        """
        # Use extra_parameter to adjust the model without changing the model.py file

    '''--------------------------------------------------------------------------------------------
    The following functions can be changed by the user and modified in each derived class of Model!
    --------------------------------------------------------------------------------------------'''

    def gas_temperature(self):
        """Calculates the gas temperature at a given position.

        Returns:
            float: Gas temperature at a given position.
        """
        return None

    def dust_temperature(self):
        """Calculates the dust temperature at a given position.

        Returns:
            float: Dust temperature at a given position.
        """
        return None

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Notes:
            Possible options:
                gas_density = 1.
                gas_density = self.math.shakura_disk(self.position)
                gas_density = self.math.sphere_density(self.position,inner_radius, outer_radius)
                gas_density = self.math.bonor_ebert_density(self.position, outer_radius, truncation_radius)

        Returns:
            float: Gas density at a given position.
        """
        return None

    def dust_density_distribution(self):
        """Calculates the dust density at a given position.

        Notes:
            Possible options:
                dust_density = 1.
                dust_density = self.math.shakura_disk(self.position)
                dust_density = self.math.sphere_density(self.position,inner_radius, outer_radius)
                dust_density = self.math.bonor_ebert_density(self.position, outer_radius, truncation_radius)

        Returns:
            float: Dust density at a given position.
        """
        return None

    def velocity_field(self):
        """Calculates the velocity at a given position.

        Notes:
            Possible options:
                velocity = [0., 0., 0.]
                velocity = self.math.kepler_rotation(self.position, stellar_mass=0.7)

        Returns:
            List[float, float, float]: Velocity at a given position.
        """
        return None

    def magnetic_field(self):
        """Calculates the magnetic field strength at a given position.

        Notes:
            Possible options:
                magnetic_field = self.math.simple_mag_field(mag_field_strength, axis)
                magnetic_field = self.math.two_simple_mag_field(self.position, mag_field_strength)
                magnetic_field = self.math.toroidal_mag_field(self.position, mag_field_strength)
                magnetic_field = self.math.poloidal_mag_field(self.position, mag_field_strength,
                torus_r_distance)
                magnetic_field = self.math.hourglass_mag_field(self.position, mag_field_strength, radius)

        Returns:
            List[float, float, float]: Magnetic field strength at a given position.
        """
        return None

    def dust_id(self):
        """Calculates the dust ID depending on the position in the grid.
        The dust ID is related to the dust composition. With this, one can
        change the dust properties inside the disk.

        Returns:
            int: dust ID.
        """
        return None

    def dust_min_size(self):
        """Calculates the minimum dust grain size depending on the position in the grid.
        This overwrites the global minimum grain size, but has no effect if it is smaller than it.

        Returns:
            float: minimum grain size
        """
        return None

    def dust_max_size(self):
        """Calculates the maximum dust grain size depending on the position in the grid.
        This overwrites the global maximum grain size, but has no effect if it is larger than it.

        Returns:
            float: maximum grain size
        """
        return None

    def scale_height(self, radius):
        """Calculates the scale height at a certain position.

        Args:
            radius (float) : Cylindrical radius of current position

        Returns:
            float: Scale height.
        """
        return None

# ----- Server class -----


class Server:
    """The Server class is the base version for each server/cluster.
    """

    def __init__(self, parse_args):
        """Initialisation of the server/cluster parameters.
        """
        self.parse_args = parse_args

        #: dict: Includes parameters related to a specific server/cluster
        self.parameter = {
            'simulation_directory': None,  # set in at file_io initialization
            'node_name': None,  # to identify the running system
            'address': '@server.uni-city.de:~/',
            'user_id': None,
            'server_polaris_dir': 'polaris/',
            'queue_system': 'PBS',
            'short_walltime': '2:00:00',
            'short_batch_class': 'express',
            'medium_walltime': '48:00:00',
            'medium_batch_class': 'medium',
            'long_walltime': '100:00:00',
            'long_batch_class': 'long',
            'walltime': None,
            'batch_class': None,
            'ram_usage': '16gb',
            'node_number': '1',
            'nr_threads': 16,
        }

    def check_for_remote(self):
        # Check if configure
        if self.parameter['user_id'] is None or self.parameter['server_polaris_dir'] is None:
            raise ValueError('polaris-remote cannot be used if userid and serverdir '
                             'are not set at the ./configure execution!')

    def check_for_walltime(self):
        if self.parse_args.queue_long:
            # Extra time for long simulations
            self.parameter['walltime'] = self.parameter['long_walltime']
            self.parameter['batch_class'] = self.parameter['long_batch_class']
        elif self.parse_args.queue_short:
            # Less time for short simulations
            self.parameter['walltime'] = self.parameter['short_walltime']
            self.parameter['batch_class'] = self.parameter['short_batch_class']
        else:
            self.parameter['walltime'] = self.parameter['medium_walltime']
            self.parameter['batch_class'] = self.parameter['long_batch_class']

    def get_command_line(self):
        """Provides server/cluster command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider a server/cluster.
        """
        new_command_line = str()
        return new_command_line

    def get_command(self):
        """Provides server/cluster command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider a server/cluster.
        """
        return self.get_command_line()


# ----- radiation source class -----
class StellarSource:
    """The StellarSource class is the base version for each stellar sources.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the radiation source parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        self.file_io = file_io
        self.parse_args = parse_args

        # Get math module
        from modules.math import Math
        self.math = Math()

        #: dict: Includes parameters of a specific stellar source
        self.parameter = {
            'position': [0., 0., 0.],
            'temperature': 0.,
            'radius': 0.,
            'luminosity': 0.,
            'mass': 0.,
            'nr_photons': 0,
            'kepler_usable': True,
        }

    def check_source_properties(self):
        """Checking for correct stellar parameters and calculating the radius,
        if only a luminosity is set.
        """
        if self.parameter['temperature'] == 0.:
            raise ValueError(
                'Effective temperature of radiation source is zero!')
        elif self.parameter['radius'] == 0.:
            if self.parameter['luminosity'] == 0.:
                raise ValueError(
                    'No radius or luminosity is set for the stellar source!')
            else:
                radius = self.math.luminosity_to_radius(self.parameter['luminosity'],
                                                        self.parameter['temperature'])
        else:
            radius = self.parameter['radius']
        return radius

    def update_parameter(self, extra_parameter):
        """Use this function to set radiation source parameter with the extra parameters and update
        radiation source parameter that depend on other parameter.
        """
        # Use extra_parameter to adjust the radiation source without changing the source.py file

    def get_command_line(self):
        """Provides radiation source command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the stellar source.
        """
        radius = self.check_source_properties()
        return '\t<source_star nr_photons = "' + str(int(self.parameter['nr_photons'])) + '">\t' \
               + str(self.parameter['position'][0]) + '\t' + str(self.parameter['position'][1]) + '\t' \
               + str(self.parameter['position'][2]) + '\t' + str(radius / self.math.const['R_sun']) + '\t' \
               + str(self.parameter['temperature']) + '\n'

    def get_command(self):
        """Provides radiation source command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the stellar source.
        """
        return self.get_command_line()


# ----- External input class -----
class ExternalInput:
    """The ExternalInput class is the base version for each external grid input.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the external grid parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        self.file_io = file_io
        self.parse_args = parse_args

        # Get math module
        from modules.math import Math
        self.math = Math()

        #: List: Data container
        self.data = None

        #: dict: Parameters that are used internally to modify the model class
        self.tmp_parameter = {
            'relative_gas_densities': None,
            'relative_dust_densities': None,
            'ignored_gas_density': 0,
            'ignored_dust_density': 0,
        }

    def external_data_loaded(self):
        if self.data is not None:
            return True
        else:
            return False

    def init_data(self):
        self.data = None

    def init_position(self, position, cell_IDs=None):
        """Initialise the grid position to calculate the nessecary cell data.

        Args:
            position (List[float, float, float]): position in model space.
            cell_IDs (List): cell_IDs of the cells (alternative to the grid position).
                Spherical -> [i_r, i_t, i_p]
                Cylindrical -> [i_r, i_p, i_z]
        """
        if cell_IDs is not None:
            self.i_r = cell_IDs[0]
            self.i_t = cell_IDs[1]
            self.i_p = cell_IDs[2]
        else:
            raise ValueError(
                'Error: External data is set, but cell_IDs of cell is not!')
        if self.i_t < 217 or self.i_t >= (217 + 128):
            self.pos_r = None
            self.pos_t = None
            self.pos_p = None
        else:
            self.i_t -= 217
            self.pos_r = self.data[0][0][self.i_r]
            self.pos_t = self.data[0][1][self.i_t]
            self.pos_p = self.data[0][2][self.i_p]

    def get_gas_temperature(self):
        """The gas temperature can be modified by the code here if neccessary.

        Returns:
            float: Gas temperature at a given position.
        """
        return self.gas_temperature()

    def get_dust_temperature(self):
        """The dust temperature can be modified by the code here if neccessary.

        Returns:
            float: Dust temperature at a given position.
        """
        return self.dust_temperature()

    def get_gas_density_distribution(self):
        """Modifies the density defined by get_gas_density_distribution() to match the preset gas masses.

        Returns:
            float: Gas density at a given position.
        """
        if self.tmp_parameter['relative_gas_densities'] is not None:
            return np.multiply(self.gas_density_distribution(), self.tmp_parameter['relative_gas_densities'])
        return self.gas_density_distribution()

    def get_dust_density_distribution(self):
        """Modifies the density defined by get_dust_density_distribution() to match the preset dust masses.

        Returns:
            float: Dust density at a given position.
        """
        if self.tmp_parameter['relative_dust_densities'] is not None:
            return np.multiply(self.dust_density_distribution(), self.tmp_parameter['relative_dust_densities'])
        return self.dust_density_distribution()

    def get_velocity_field(self):
        """The velocity field can be modified by the code here if neccessary.

        Returns:
            List[float, float, float]: Velocity at a given position.
        """
        return self.velocity_field()

    def get_magnetic_field(self):
        """The magnetic field can be modified by the code here if neccessary.

        Returns:
            List[float, float, float]: Magnetic field strength at a given position.
        """
        return self.magnetic_field()

    def get_dust_id(self):
        """The dust choice id can be modified by the code here if neccessary.

        Returns:
            int: dust ID.
        """
        return self.dust_id()

    def adjust(self, value, name):
        """Adjust a parameter of the model.

        Args:
            value (float): New value.
            name (str): Name of model parameter.
        """
        self.tmp_parameter[name] = value

    def ignore_cell(self, node=None):
        """Ignore a cell for grid refinement, if necessary for a given model.

        Args:
            node: Current node to check if it should be ignored for grid refinement

        Returns:
            bool: True if ignore, False if considered for refinement.
        """
        return False

    '''--------------------------------------------------------------------------------------------
    The following functions can be changed by the user and modified in each derived class of Model!
    --------------------------------------------------------------------------------------------'''

    def gas_temperature(self):
        """Calculates the gas temperature at a given position.

        Returns:
            float: Gas temperature at a given position.
        """
        gas_temperature = 0.
        return gas_temperature

    def dust_temperature(self):
        """Calculates the dust temperature at a given position.

        Returns:
            float: Dust temperature at a given position.
        """
        dust_temperature = 0.
        return dust_temperature

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Notes:
            Possible options:
                gas_density = 1.
                gas_density = self.math.shakura_disk(self.position)
                gas_density = self.math.sphere_density(self.position,inner_radius, outer_radius)
                gas_density = self.math.bonor_ebert_density(self.position, outer_radius, truncation_radius)

        Returns:
            float: Gas density at a given position.
        """
        gas_density = 0.
        return gas_density

    def dust_density_distribution(self):
        """Calculates the dust density at a given position.

        Notes:
            Possible options:
                dust_density = 1.
                dust_density = self.math.shakura_disk(self.position)
                dust_density = self.math.sphere_density(self.position,inner_radius, outer_radius)
                dust_density = self.math.bonor_ebert_density(self.position, outer_radius, truncation_radius)

        Returns:
            float: Dust density at a given position.
        """
        dust_density = 0.
        return dust_density

    def velocity_field(self):
        """Calculates the velocity at a given position.

        Notes:
            Possible options:
                velocity = [0., 0., 0.]
                velocity = self.math.kepler_rotation(self.position, stellar_mass=0.7)

        Returns:
            List[float, float, float]: Velocity at a given position.
        """
        velocity = [0., 0., 0.]
        return velocity

    def magnetic_field(self):
        """Calculates the magnetic field strength at a given position.

        Notes:
            Possible options:
                magnetic_field = self.math.simple_mag_field(mag_field_strength, axis)
                magnetic_field = self.math.two_simple_mag_field(self.position, mag_field_strength)
                magnetic_field = self.math.toroidal_mag_field(self.position, mag_field_strength)
                magnetic_field = self.math.poloidal_mag_field(self.position, mag_field_strength,
                torus_r_distance)
                magnetic_field = self.math.hourglass_mag_field(self.position, mag_field_strength, radius)

        Returns:
            List[float, float, float]: Magnetic field strength at a given position.
        """
        magnetic_field = [0., 0., 0.]
        return magnetic_field

    def dust_id(self):
        """Calculates the dust ID depending on the position in the grid.
        The dust ID is related to the dust composition. With this, one can
        change the dust properties inside the disk.

        Returns:
            int: dust ID.
        """
        dust_id = 0
        return dust_id
