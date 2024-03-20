"""base.py includes the base versions of each Polaris component to be loaded by other scripts.
"""

import numpy as np
from polaris_tools_modules.math import Math

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
            'gas_mass': None,
            'dust_mass': None,
            # Global extent variables to set all grid types at once
            'inner_radius': 1. * self.math.const['au'],
            'outer_radius': 100. * self.math.const['au'],
            'mass_fraction': 0.01,
            'grid_type': 'octree',
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
            'sf_z': -2.0,
            # This list is used as cell borders if sf_r or sf_z is zero
            'radius_list': [],
            'phi_list': [],
            'z_list': [],
            # Split the first radial cell into multiple
            'split_first_cell': 1
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

    def get_dust_size_param(self):
        """The size distribution parameter can be modified by the code here if neccessary.

        Returns:
            float: minimum grain size
        """
        return self.dust_size_param()

    def get_dz(self, radius):
        """Calculates the width of each vertical cell border depending on the radial position.

        Args:
            radius (float) : Cylindrical radius of current position

        Returns:
            float: Width between two cell borders.
        """
        return 10. * self.get_scale_height(radius) / self.cylindrical_parameter['n_z']

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
                gas_density = self.math.const_sphere_density(self.position,inner_radius, outer_radius)
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
                dust_density = self.math.const_sphere_density(self.position,inner_radius, outer_radius)
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

    def dust_size_param(self):
        """Calculates the size distribution parameter depending on the position in the grid.
        This overwrites the global size distribution parameter.

        Returns:
            float: minimum grain size
        """
        return None

    def get_scale_height(self, radius):
        """Calculates the scale height at a certain position.

        Args:
            radius (float) : Cylindrical radius of current position

        Returns:
            float: Scale height.
        """
        return None
