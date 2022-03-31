#!/usr/bin/env python
# -*- coding: utf-8 -*-

import struct
from sys import stdout

import numpy as np


class Grid:
    """This is the base class to create various grids based on the models defined in model.py.
    """

    def __init__(self, model, path, parse_args):
        """Initialisation of grid parameters.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
            path : Handles file input/output and all necessary paths.
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
        """
        self.model = model
        self.path = path
        self.parse_args = parse_args

        # Get math module
        from polaris_tools_modules.math import Math
        self.math = Math()

        # Define which class is responsible for the grid data
        # if ext_input is not None:
        #     self.data = ext_input
        #     self.data.init_data()
        # else:
        self.data = model

        # Set the position to an arbitrary value to check get_density functions
        self.data.position = [0, 0, 0]
        self.data.volume = 0.

        #: float: Total gas mass of the grid nodes
        self.total_gas_mass = None

        #: float: Total dust mass of the grid nodes
        self.total_dust_mass = None

        # Number of gas density distributions
        if self.model.parameter['gas_mass'] is not None:
            if isinstance(self.model.parameter['gas_mass'], (float, int)):
                self.nr_gas_densities = 1
            else:
                self.nr_gas_densities = len(self.model.parameter['gas_mass'])
        else:
            self.nr_gas_densities = 0
        # Number of dust density distributions
        if self.model.parameter['dust_mass'] is not None:
            if isinstance(self.model.parameter['dust_mass'], (float, int)):
                self.nr_dust_densities = 1
            else:
                self.nr_dust_densities = len(self.model.parameter['dust_mass'])
        else:
            self.nr_dust_densities = 0

        # Check if get_density functions provide fitting arrays
        self.check_density_arrays()

        #: int: Number of quantities per grid cell
        self.data_length = self.nr_gas_densities + self.nr_dust_densities + \
            (self.data.get_dust_temperature() is not None) + \
            (self.data.get_gas_temperature() is not None) + \
            (self.data.get_magnetic_field() is not None) * 3 + \
            (self.data.get_velocity_field() is not None) * 3 + \
            (self.data.dust_min_size() is not None) + \
            (self.data.dust_max_size() is not None) + \
            (self.data.dust_size_param() is not None)
        if self.data.get_dust_id() is not None:
            if self.nr_gas_densities > 1 or self.nr_dust_densities > 1:
                raise ValueError(
                    'Multiple density distributions and dust_id function defined')
            else:
                self.data_length += 1

    def update_mass_measurement(self, node, remove=False):
        """Updates the total mass for normalization and allows the definition
        of custom regions to track their mass too.

        Args:
            node: A grid cell.
            remove (bool): Remove the density instead of adding?
        """
        if self.nr_gas_densities > 0:
            gas_mass = np.multiply(
                node.parameter['gas_density'], node.parameter['volume'])
            if remove:
                gas_mass *= -1
            if self.total_gas_mass is None:
                self.total_gas_mass = gas_mass
            else:
                self.total_gas_mass += gas_mass

        if self.nr_dust_densities > 0:
            dust_mass = np.multiply(
                node.parameter['dust_density'], node.parameter['volume'])
            if remove:
                dust_mass *= -1
            if self.total_dust_mass is None:
                self.total_dust_mass = dust_mass
            else:
                self.total_dust_mass += dust_mass

    def check_density_arrays(self):
        """Check if get_density functions provide fitting arrays.
        """
        # Check for errors
        try:
            if self.nr_gas_densities == 1:
                if self.data.get_gas_density_distribution() is not None and \
                        not isinstance(self.data.get_gas_density_distribution(), (float, int)):
                    if len(self.data.get_gas_density_distribution()[0]) != \
                            len(self.model.parameter['gas_mass'][0]):
                        raise ValueError("gas_density_distribution does not provied the same array than "
                                         "defined in self.parameter['gas_mass']")
            elif self.nr_gas_densities > 1:
                for i_gas_dens in range(self.nr_gas_densities):
                    if len(self.data.get_gas_density_distribution()[i_gas_dens]) != \
                            len(self.model.parameter['gas_mass'][i_gas_dens]):
                        raise ValueError("gas_density_distribution does not provied the same array than "
                                         "defined in self.parameter['gas_mass']")
        except:
            raise ValueError(
                "the gas_density function and the defined gas_mass do not fit!")
        try:
            if self.nr_dust_densities == 1:
                if self.data.get_dust_density_distribution() is not None and \
                        not isinstance(self.data.get_dust_density_distribution(), (float, int)):
                    raise ValueError(
                        "get_dust_density_distribution provides not float!")
            elif self.nr_dust_densities > 1:
                for i_dust_dens in range(self.nr_dust_densities):
                    if len(self.data.get_dust_density_distribution()[i_dust_dens]) != \
                            len(self.model.parameter['dust_mass'][i_dust_dens]):
                        raise ValueError("dust_density_distribution does not provied the same array than "
                                         "defined in self.parameter['dust_mass']")
        except:
            raise ValueError(
                "the dust_density function and the defined dust_mass do not fit!")

    def write_header(self, grid_file, grid_type='', root=None):
        """Writes general header to binary file.

        Args:
            grid_file: Input grid file (tmp_grid).
            grid_type (str): Name of the grid type.
                (octree, spherical, cylindrical)
            root: Root node of the octree grid.

        Notes:
            GRIDgas_dens   0
            GRIDdust_dens  1
            GRIDdust_temp  2
            GRIDgas_temp   3
            GRIDmx         4
            GRIDmy         5
            GRIDmz         6
            GRIDvx         7
            GRIDvy         8
            GRIDvz         9
            GRIDpx         10
            GRIDpy         11
            GRIDpz         12
            GRIDa_alg      13
            GRIDa_min      14
            GRIDa_max      15
            GRIDq          16
            GRIDratio      17
            GRIv_turb      18
            GRIDPDA        19
            GRIDopiate     20
            GRIDdust_id    21
            GRIDn_th       22  //number density of thermal electrons
            GRIDT_e        23  //Temperature of thermal electrons
            GRIDn_cr       24  //number density of CR electrons
            GRIDg_min      25  //gamma min for power law distribution
            GRIDg_max      26  //gamma max for power law distribution
            GRIDp          27  //power law exponent
            GRIDgas_mdens  28
            GRIDdust_mdens 29
        """
        if grid_type == 'octree':
            grid_id = 20
        elif grid_type == 'spherical':
            grid_id = 30
        elif grid_type == 'cylindrical':
            grid_id = 40
        else:
            raise ValueError('Grid type: ' + str(grid_type) + 'is not known!')

        grid_file.write(struct.pack('H', grid_id))
        grid_file.write(struct.pack('H', self.data_length))
        # Gas density index: 0 or gas mass density index: 28
        for i_gas_dens in range(self.nr_gas_densities):
            grid_file.write(struct.pack('H', 28))
        # Dust density index: 1 or dust mass density index: 29
        for i_dust_dens in range(self.nr_dust_densities):
            grid_file.write(struct.pack('H', 29))
        if self.data.get_dust_temperature() is not None:
            # Dust temperature index: 2
            grid_file.write(struct.pack('H', 2))
        if self.data.get_gas_temperature() is not None:
            # Gas temperature index: 3
            grid_file.write(struct.pack('H', 3))
        if self.data.get_magnetic_field() is not None:
            # Magnetic field x-component index: 4
            grid_file.write(struct.pack('H', 4))
            # Magnetic field y-component index: 5
            grid_file.write(struct.pack('H', 5))
            # Magnetic field z-component index: 6
            grid_file.write(struct.pack('H', 6))
        if self.data.get_velocity_field() is not None:
            # Velocity field x-component index: 7
            grid_file.write(struct.pack('H', 7))
            # Velocity field y-component index: 8
            grid_file.write(struct.pack('H', 8))
            # Velocity field z-component index: 9
            grid_file.write(struct.pack('H', 9))
        # Add an index for the dust choice in a certain cell
        if self.data.get_dust_id() is not None:
            # Dust choice index: 21
            grid_file.write(struct.pack('H', 21))
        if self.data.dust_min_size() is not None:
            # Minimum dust grain size: 14
            grid_file.write(struct.pack('H', 14))
        if self.data.dust_max_size() is not None:
            # Maximum dust grain size: 15
            grid_file.write(struct.pack('H', 15))
        if self.data.dust_size_param() is not None:
            # Size distribution parameter: 16
            grid_file.write(struct.pack('H', 16))
        if grid_type == 'octree':
            if root is not None:
                grid_file.write(struct.pack(
                    'd', self.model.octree_parameter['sidelength']))
                grid_file.write(struct.pack('H', root.parameter['is_leaf']))
                grid_file.write(struct.pack('H', root.parameter['level']))
            else:
                raise ValueError(
                    'root node has to be defined for writing the grid header!')

    def write_node_data(self, grid_file, node, data_type='f', cell_IDs=None, rewrite=False):
        """Write data of each node.

        Args:
            grid_file: Input grid file (tmp_grid).
            node: Instance of octree node.
            data_type (str): Type of the node data ('f': float or 'd': double).
            cell_IDs (List): indices of the cells (used for external purpose).
                Spherical -> [i_r, i_t, i_p]
                Cylindrical -> [i_r, i_p, i_z]
            rewrite (bool): Remove the the last cell data to overwrite it.
        """
        # Set data type length in bytes depending on the data type
        if data_type == 'f':
            data_type_length = 4
        elif data_type == 'd':
            data_type_length = 8
        else:
            raise ValueError(
                'Do not understand the data type ' + data_type + ' in grid.py!')

        if rewrite:
            grid_file.seek(-(self.data_length * data_type_length), 1)

        # Transmit the cell position to the model to obtain the cell quantities
        self.data.init_position(node, cell_IDs)

        # Write gas density to grid for each defined distribution
        if self.nr_gas_densities > 0:
            node.parameter['gas_density'] = self.data.get_gas_density_distribution()
            for i_gas_dens in range(self.nr_gas_densities):
                if isinstance(node.parameter['gas_density'], (float, int)):
                    grid_file.write(struct.pack(
                        data_type, node.parameter['gas_density']))
                else:
                    grid_file.write(struct.pack(data_type, np.sum(
                        node.parameter['gas_density'][i_gas_dens])))
        # Write dust density to grid for each defined distribution
        if self.nr_dust_densities > 0:
            node.parameter['dust_density'] = self.data.get_dust_density_distribution()
            for i_dust_dens in range(self.nr_dust_densities):
                if isinstance(node.parameter['dust_density'], (float, int)):
                    grid_file.write(struct.pack(
                        data_type, node.parameter['dust_density']))
                else:
                    grid_file.write(struct.pack(data_type, np.sum(
                        node.parameter['dust_density'][i_dust_dens])))

        # Set additional quantities
        dust_temperature = self.data.get_dust_temperature()
        if dust_temperature is not None:
            grid_file.write(struct.pack(
                data_type, dust_temperature))
        gas_temperature = self.data.get_gas_temperature()
        if gas_temperature is not None:
            grid_file.write(struct.pack(
                data_type, gas_temperature))

        mag_field = self.data.get_magnetic_field()
        if mag_field is not None:
            grid_file.write(struct.pack(data_type, mag_field[0]))
            grid_file.write(struct.pack(data_type, mag_field[1]))
            grid_file.write(struct.pack(data_type, mag_field[2]))
        velocity = self.data.get_velocity_field()
        if velocity is not None:
            grid_file.write(struct.pack(data_type, velocity[0]))
            grid_file.write(struct.pack(data_type, velocity[1]))
            grid_file.write(struct.pack(data_type, velocity[2]))

        dust_id = self.data.get_dust_id()
        if dust_id is not None:
            grid_file.write(struct.pack(data_type, dust_id))
        a_min = self.data.get_dust_min_size()
        if a_min is not None:
            grid_file.write(struct.pack(data_type, a_min))
        a_max = self.data.get_dust_max_size()
        if a_max is not None:
            grid_file.write(struct.pack(data_type, a_max))
        size_param = self.data.get_dust_size_param()
        if size_param is not None:
            grid_file.write(struct.pack(data_type, size_param))

    def read_write_node_data(self, tmp_file, grid_file, data_type='f'):
        """Read node data to binary grid_file.

        Args:
            tmp_file: Input grid file (tmp_grid).
            grid_file: Output grid file (final_grid).
            data_type (str): Type of the node data ('f': float or 'd': double).
        """
        # Set data type length in bytes depending on the data type
        if data_type == 'f':
            data_type_length = 4
        elif data_type == 'd':
            data_type_length = 8
        else:
            raise ValueError(
                'Do not understand the data type ' + data_type + ' in grid.py!')

        # Calculate normalized density
        for i_gas_dens in range(self.nr_gas_densities):
            gas_density = tmp_file.read(data_type_length)
            normed_gas_density = 0.
            if isinstance(np.sum(self.total_gas_mass), (float, int)) and np.sum(self.total_gas_mass) > 0.:
                normed_gas_density = struct.unpack(data_type, gas_density)[0] * \
                    (np.sum(
                        self.model.parameter['gas_mass']) / np.sum(self.total_gas_mass))
            grid_file.write(struct.pack(data_type, normed_gas_density))
        for i_dust_dens in range(self.nr_dust_densities):
            dust_density = tmp_file.read(data_type_length)
            normed_dust_density = 0.
            if isinstance(np.sum(self.total_dust_mass), (float, int)) and np.sum(self.total_dust_mass) > 0.:
                normed_dust_density = struct.unpack(data_type, dust_density)[0] * \
                    (np.sum(
                        self.model.parameter['dust_mass']) / np.sum(self.total_dust_mass))
            grid_file.write(struct.pack(data_type, normed_dust_density))
        for i in range(self.data_length - self.nr_gas_densities - self.nr_dust_densities):
            grid_file.write(tmp_file.read(data_type_length))

    # def update_grid(self, grid_file, tmp_file, reverse):
    #     """Update grid to be in agreement with POLARIS newest version.

    #     Args:
    #         grid_file: Input grid file (previous grid).
    #         tmp_file: Output grid file (updated grid).
    #         reverse (bool): Downgrading instead?
    #     """
    #     grid_id = grid_file.read(2)
    #     tmp_file.write(grid_id)

    #     data_length = grid_file.read(2)
    #     tmp_file.write(data_length)

    #     for i in range(struct.unpack('H', data_length)[0]):
    #         tmp_file.write(grid_file.read(2))

    #     if struct.unpack('H', grid_id)[0] == 20:
    #         print('Nothing to do!')
    #     elif struct.unpack('H', grid_id)[0] == 30:
    #         tmp_file.write(grid_file.read(8))
    #         tmp_file.write(grid_file.read(8))
    #         n_r = grid_file.read(2)
    #         tmp_file.write(n_r)
    #         n_ph = grid_file.read(2)
    #         tmp_file.write(n_ph)
    #         n_th = grid_file.read(2)
    #         tmp_file.write(n_th)
    #         sf_r = grid_file.read(8)
    #         tmp_file.write(sf_r)
    #         if not reverse:
    #             # Add log_Phi value
    #             tmp_file.write(struct.pack('d', 1.0))
    #         else:
    #             # Ignore log_Phi value
    #             sf_ph = grid_file.read(8)
    #         byte = grid_file.read(1)
    #         while byte != b'':
    #             tmp_file.write(byte)
    #             byte = grid_file.read(1)
    #     elif struct.unpack('H', grid_id)[0] == 40:
    #         tmp_file.write(grid_file.read(8))
    #         tmp_file.write(grid_file.read(8))
    #         tmp_file.write(grid_file.read(8))
    #         n_r = grid_file.read(2)
    #         tmp_file.write(n_r)
    #         n_ph = grid_file.read(2)
    #         tmp_file.write(n_ph)
    #         n_z = grid_file.read(2)
    #         tmp_file.write(n_z)
    #         sf_r = grid_file.read(8)
    #         tmp_file.write(sf_r)
    #         if not reverse:
    #             # Add log_Phi value
    #             tmp_file.write(struct.pack('d', 1.0))
    #         else:
    #             # Ignore log_Phi value
    #             sf_ph = grid_file.read(8)
    #         byte = grid_file.read(1)
    #         while byte != b'':
    #             tmp_file.write(byte)
    #             byte = grid_file.read(1)
    #     else:
    #         raise ValueError(
    #             'Grid index: ' + str(struct.unpack('H', grid_id)[0]) + 'is not known!')

    # def set_quantiy_in_grid(self, grid_file, tmp_file, quantity_id, quantity_value):
    #     """Set the value of the desired quantity to value in each cell of the grid.

    #     Args:
    #         grid_file: Input grid file (previous grid).
    #         tmp_file: Output grid file (updated grid).
    #         quantity_id (int): Index of the quantity in the POLARIS grid.
    #         quantity_value (float): Value to which the quantity should be set.
    #     """
    #     grid_id = grid_file.read(2)
    #     tmp_file.write(grid_id)

    #     data_length = struct.unpack('H', grid_file.read(2))[0]
    #     pos = None
    #     found = False
    #     id_list = []
    #     for i in range(data_length):
    #         tmp_quantity_id = grid_file.read(2)
    #         if struct.unpack('H', tmp_quantity_id)[0] == quantity_id:
    #             if pos is not None:
    #                 raise ValueError(
    #                     'Desired quantity is set in grid multiple times!')
    #             pos = i
    #             found = True
    #             print('HINT: Quantity found in grid!')
    #         id_list.append(tmp_quantity_id)
    #     if not found:
    #         pos = data_length
    #         data_length += 1
    #         id_list.append(struct.pack('H', quantity_id))

    #     tmp_file.write(struct.pack('H', data_length))
    #     for i in range(data_length):
    #         tmp_file.write(id_list[i])
    #     if struct.unpack('H', grid_id)[0] == 20:
    #         tmp_file.write(grid_file.read(8))
    #         tmp_file.write(grid_file.read(2))
    #         tmp_file.write(grid_file.read(2))
    #         data_type = 'f'
    #     elif struct.unpack('H', grid_id)[0] == 30:
    #         tmp_file.write(grid_file.read(8))
    #         tmp_file.write(grid_file.read(8))
    #         # Get number of radial cells
    #         n_r = grid_file.read(2)
    #         tmp_file.write(n_r)
    #         # Get number of phi cells
    #         n_ph = grid_file.read(2)
    #         tmp_file.write(n_ph)
    #         # Get number of theta cells
    #         n_th = grid_file.read(2)
    #         tmp_file.write(n_th)
    #         # Get step width factors (zero for custom)
    #         sf_r = grid_file.read(8)
    #         tmp_file.write(sf_r)
    #         sf_ph = grid_file.read(8)
    #         tmp_file.write(sf_ph)
    #         sf_th = grid_file.read(8)
    #         tmp_file.write(sf_th)
    #         # Write the custom cell distribution if chosen
    #         if struct.unpack('d', sf_r)[0] == 0:
    #             for i_r in range(struct.unpack('H', n_r)[0] - 1):
    #                 tmp_file.write(grid_file.read(8))
    #         if struct.unpack('d', sf_ph)[0] == 0:
    #             for i_ph in range(struct.unpack('H', n_ph)[0]):
    #                 tmp_file.write(grid_file.read(8))
    #         if struct.unpack('d', sf_th)[0] == 0:
    #             for i_th in range(struct.unpack('H', n_th)[0]):
    #                 tmp_file.write(grid_file.read(8))
    #         data_type = 'd'
    #     elif struct.unpack('H', grid_id)[0] == 40:
    #         tmp_file.write(grid_file.read(8))
    #         tmp_file.write(grid_file.read(8))
    #         tmp_file.write(grid_file.read(8))
    #         # Get number of radial cells
    #         n_r = grid_file.read(2)
    #         tmp_file.write(n_r)
    #         # Get number of phi cells
    #         n_ph = grid_file.read(2)
    #         tmp_file.write(n_ph)
    #         # Get number of theta cells
    #         n_z = grid_file.read(2)
    #         tmp_file.write(n_z)
    #         # Get step width factors (zero for custom)
    #         sf_r = grid_file.read(8)
    #         tmp_file.write(sf_r)
    #         sf_ph = grid_file.read(8)
    #         tmp_file.write(sf_ph)
    #         sf_z = grid_file.read(8)
    #         tmp_file.write(sf_z)
    #         if struct.unpack('d', sf_r)[0] == 0:
    #             for i_r in range(struct.unpack('H', n_r)[0] - 1):
    #                 tmp_file.write(grid_file.read(8))
    #         if struct.unpack('d', sf_ph)[0] == 0:
    #             for i_ph in range(struct.unpack('H', n_ph)[0] - 1):
    #                 tmp_file.write(grid_file.read(8))
    #         elif struct.unpack('d', sf_ph)[0] == -1:
    #             for i_r in range(struct.unpack('H', n_r)[0]):
    #                 tmp_file.write(grid_file.read(2))
    #         if struct.unpack('d', sf_z)[0] == 0:
    #             for i_z in range(struct.unpack('H', n_z)[0] - 1):
    #                 tmp_file.write(grid_file.read(8))
    #         elif struct.unpack('d', sf_z)[0] == -1:
    #             for i_r in range(struct.unpack('H', n_r)[0]):
    #                 tmp_file.write(grid_file.read(8))
    #         data_type = 'd'
    #     else:
    #         raise ValueError(
    #             'Grid index: ' + str(struct.unpack('H', grid_id)[0]) + 'is not known!')

    #     if data_type == 'f':
    #         data_type_length = 4
    #     elif data_type == 'd':
    #         data_type_length = 8

    #     is_leaf = grid_file.read(2)
    #     while is_leaf != b'':
    #         tmp_file.write(is_leaf)
    #         level = grid_file.read(2)
    #         tmp_file.write(level)
    #         if struct.unpack('H', is_leaf)[0]:
    #             for i in range(data_length):
    #                 if i == pos:
    #                     tmp_file.write(struct.pack(data_type, quantity_value))
    #                     if found:
    #                         grid_file.read(data_type_length)
    #                 else:
    #                     tmp_file.write(grid_file.read(data_type_length))
    #         is_leaf = grid_file.read(2)


class OcTree(Grid):
    """This class creates OcTree grids based on the models defined in model.py.
    """

    def __init__(self, model, path, parse_args):
        """Initialization of grid parameters.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
            path : Handles file input/output and all necessary paths.
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
        """
        Grid.__init__(self, model, path, parse_args)

    def init_root(self):
        """Initialise the root node.

        Return:
            Root node of the octree grid.
        """
        # Create root node
        root = Node('octree')
        # Set root node position
        root.parameter['position'] = [0., 0., 0.]
        # Set width between node borders
        root.parameter['sidelength'] = self.model.octree_parameter['sidelength']
        # Set node volume
        root.parameter['volume'] = self.get_volume(node=root)
        return root

    @staticmethod
    def write_node_header(grid_file, node, rewrite=False):
        """Write header of each node.

        Args:
            grid_file: Input grid file (tmp_grid).
            node: Instance of octree node.
            rewrite (bool): Remove the last 4 bytes to rewrite the header.
        """
        if rewrite:
            grid_file.seek(-4, 1)
        grid_file.write(struct.pack('H', node.parameter['is_leaf']))
        grid_file.write(struct.pack('H', node.parameter['level']))

    def create_grid(self, grid_file, node, max_tree_level=None, refinement_limit=0.1):
        """Create an octree grid and calculate the total mass of the grid nodes.

        Args:
            grid_file: Input grid file (tmp_grid).
            node: Instance of grid node.
            max_tree_level (int): Maximum number of grid levels.
            refinement_limit (float): maximum of the density difference between the 8
                children nodes and the parent node.
        """
        # Set max tree level from user input.
        max_tree_level = self.model.octree_parameter['max_tree_level']
        # Show percentage based on the second and third grid level
        tmp_node = node
        while tmp_node.parameter['level'] > 2:
            tmp_node = tmp_node.parent
        percentage_count = 0
        if tmp_node.parameter['level'] == 2:
            percentage_count = tmp_node.parameter['index'] + \
                tmp_node.parent.parameter['index'] * 8
        elif tmp_node.parameter['level'] == 1:
            percentage_count = tmp_node.parameter['index'] * 8
        stdout.write('--- Generate cartesian grid: ' +
                     str(int(100 * percentage_count / (8 * 8 - 1))) + ' %      \r')
        stdout.flush()
        if node.parameter['level'] < max_tree_level and not self.model.ignore_cell(node):
            # Add 8 children to node
            self.add_level(node=node)
            for i_leaf in range(8):
                # Refer to current children node
                node = node.children[i_leaf]
                # Write the header of the node
                self.write_node_header(grid_file=grid_file, node=node)
                # Recursive execute this function to add another 8 children to node
                self.create_grid(grid_file=grid_file, node=node, max_tree_level=max_tree_level,
                                 refinement_limit=refinement_limit)
                # Go to parent node
                node = node.parent
                if i_leaf == 7:
                    # Calculate a difference between various quantities to do grid refinement
                    difference = self.grid_refinement(node=node)
                    if difference < refinement_limit and node.parameter['level'] > 3:
                        # If the difference is small enough and the level larger than 3,
                        # delete the children and use the parent node only.
                        for j_leaf in range(8):
                            # Subtract the children mass from total mass
                            self.remove_node_from_grid(
                                grid_file=grid_file, node=node.children[j_leaf])
                        # Parent node is now a leaf
                        node.parameter['is_leaf'] = True
                        # Rewrite the header of the parent node
                        self.write_node_header(
                            grid_file=grid_file, node=node, rewrite=True)
                        # Write data of the parent node
                        self.write_node_data(grid_file=grid_file, node=node)
                        # Add the mass of the parent node to the total mass
                        self.update_mass_measurement(node=node)
        else:
            # Initialize and write a leaf node including data
            node.parameter['is_leaf'] = True
            # Rewrite the node header to include is_leaf = True
            self.write_node_header(grid_file=grid_file,
                                   node=node, rewrite=True)
            # Write data of the node
            self.write_node_data(grid_file=grid_file, node=node)
            # Add the mass of the node to the total mass
            self.update_mass_measurement(node=node)

    @staticmethod
    def add_level(node):
        """Add 8 children nodes to an octree node.

        Args:
            node: Instance of octree node.
        """
        #: float: Distance to a children node on the x-axis
        d_x = node.parameter['sidelength'] / 4.
        #: float: Distance to a children node on the y-axis
        d_y = node.parameter['sidelength'] / 4.
        #: float: Distance to a children node on the z-axis
        d_z = node.parameter['sidelength'] / 4.
        #: Node is not a leaf anymore
        node.parameter['is_leaf'] = False
        for i_leaf in range(8):
            # Each children is a new node
            node.children[i_leaf] = Node('octree')
            # Calculate the position of the children nodes
            position = [0., 0., 0.]
            if i_leaf == 0:
                position = np.add(
                    node.parameter['position'], [-d_x, -d_y, -d_z])
            elif i_leaf == 1:
                position = np.add(
                    node.parameter['position'], [d_x, -d_y, -d_z])
            elif i_leaf == 2:
                position = np.add(
                    node.parameter['position'], [-d_x, d_y, -d_z])
            elif i_leaf == 3:
                position = np.add(node.parameter['position'], [d_x, d_y, -d_z])
            elif i_leaf == 4:
                position = np.add(
                    node.parameter['position'], [-d_x, -d_y, d_z])
            elif i_leaf == 5:
                position = np.add(node.parameter['position'], [d_x, -d_y, d_z])
            elif i_leaf == 6:
                position = np.add(node.parameter['position'], [-d_x, d_y, d_z])
            elif i_leaf == 7:
                position = np.add(node.parameter['position'], [d_x, d_y, d_z])
            # Set the position of the children nodes
            node.children[i_leaf].parameter['position'] = position
            # Set the extent of the children nodes
            node.children[i_leaf].parameter['sidelength'] = np.divide(
                node.parameter['sidelength'], 2.)
            # Set the volume of the children nodes
            node.children[i_leaf].parameter['volume'] = node.parameter['volume'] / 8.
            # Set the level of the children nodes
            node.children[i_leaf].parameter['level'] = node.parameter['level'] + 1
            # Set the index of the children nodes (from 0 to 7)
            node.children[i_leaf].parameter['index'] = i_leaf
            # Set the parent node of the children nodes
            node.children[i_leaf].parent = node

    def remove_node_from_grid(self, grid_file, node):
        """Remove node from binary grid file.

        Args:
            grid_file: Output grid file.
            node: Instance of octree node.

        Returns:
            Volume of the node
        """
        if node.parameter['is_leaf'] == True:
            # Subtract the children mass from total mass
            self.update_mass_measurement(node=node, remove=True)
            # Go some bytes back to remove the children from binary file
            grid_file.seek(-(self.data_length * 4 + 2 + 2), 1)
            # Remove written nodes from memory
            del node
        else:
            # Go some bytes back to remove the children from binary file
            grid_file.seek(-(2 + 2), 1)
            for j_leaf in range(8):
                # Subtract the children mass from total mass
                self.remove_node_from_grid(grid_file, node.children[j_leaf])

    @staticmethod
    def get_volume(node):
        """Calculate the volume of an octree node.

        Args:
            node: Instance of octree node.

        Returns:
            Volume of the node
        """
        volume = (node.parameter['sidelength'] *
                  node.parameter['sidelength'] * node.parameter['sidelength'])
        return volume

    def grid_refinement(self, node):
        """Calculates the maximum difference of the density in the
        center of the parent and in each of the 8 children.

        Args:
            node: Instance of octree node.

        Returns:
            Float: Maximum difference of the density in the
            center of the parent and in each of the 8 children.
        """
        # Set density of parent node for refinement calculation
        self.data.init_position(node)
        node.parameter['gas_density'] = self.model.get_gas_density_distribution()
        if isinstance(node.parameter['gas_density'], (float, int)):
            if np.sum(node.parameter['gas_density']) > 0:
                return 9999
        else:
            for i_gas_dens, gas_dens in enumerate(node.parameter['gas_density']):
                if np.sum(gas_dens) > 0:
                    return 9999
        return 0
        # node_mass = node.parameter['gas_density'] * node.parameter['volume']
        # node_leaf_mass = 0
        # for i_leaf in range(8):
        #    node_leaf_mass += node.children[i_leaf].parameter['gas_density'] * node.children[i_leaf].parameter['volume']
        # if (node_mass + node_leaf_mass) > 0.:
        #    difference = abs(node_mass - node_leaf_mass) / (node_mass + node_leaf_mass)
        # else:
        #    difference = 0.
        # self.grid_refinement_extra_mag()
        # self.grid_refinement_extra_t_gas()
        # difference = 99999.0
        # return difference

    # def grid_refinement_extra_mag(self, node):
    #     """Calculates grid refinement from the magnetic field strength.

    #     Args:
    #         node: Instance of octree node.

    #     Returns:
    #         Float: Maximum difference of the magnetic field strength in the
    #         center of the parent and in each of the 8 children.
    #     """
    #     difference = 0.
    #     self.data.init_position(node)
    #     parent_mag = self.model.get_magnetic_field()
    #     for i_leaf in range(8):
    #         self.data.init_position(node.children[i_leaf])
    #         children_mag = self.model.get_magnetic_field()
    #         for i in range(3):
    #             if abs(parent_mag[i] + children_mag[i]) > 0:
    #                 diff = abs(parent_mag[i] - children_mag[i]) / \
    #                     abs(parent_mag[i] + children_mag[i])
    #                 difference = max(difference, diff)
    #     return difference

    # def grid_refinement_extra_t_gas(self, node):
    #     """Calculates grid refinement from the temperature.

    #     Args:
    #         node: Instance of octree node.

    #     Returns:
    #         Float: Maximum difference of the gas temperature in the
    #         center of the parent and in each of the 8 children.
    #     """
    #     difference = 0.
    #     self.data.init_position(node)
    #     parent_t_gas = self.model.get_gas_temperature()
    #     for i_leaf in range(8):
    #         self.data.init_position(node.children[i_leaf])
    #         children_t_gas = self.model.get_gas_temperature()
    #         if abs(parent_t_gas + children_t_gas) > 0:
    #             diff = 10 * abs(parent_t_gas - children_t_gas) / \
    #                 abs(parent_t_gas + children_t_gas)
    #             difference = max(difference, diff)
    #     return difference

    def normalize_density(self, tmp_file, grid_file):
        """Read the temporary octree grid and normalize the
        density to the total model mass.

        Args:
            tmp_file: Input grid file (tmp_grid).
            grid_file: Output grid file (final_grid).
        """
        self.read_write_header(tmp_file=tmp_file, grid_file=grid_file)
        while True:
            is_leaf = self.read_write_node_header(
                tmp_file=tmp_file, grid_file=grid_file)
            if len(is_leaf) == 0:
                break
            is_leaf = struct.unpack('H', is_leaf)[0]
            if is_leaf == 1:
                self.read_write_node_data(
                    tmp_file=tmp_file, grid_file=grid_file)

    def read_write_header(self, tmp_file, grid_file):
        """Read and write octree header from binary file.

        Args:
            tmp_file: Input grid file (tmp_grid).
            grid_file: Output grid file (final_grid).
        """
        grid_file.write(tmp_file.read(2))
        grid_file.write(tmp_file.read(2))
        for i in range(self.data_length):
            grid_file.write(tmp_file.read(2))
        grid_file.write(tmp_file.read(8))
        grid_file.write(tmp_file.read(2))
        grid_file.write(tmp_file.read(2))

    @staticmethod
    def read_write_node_header(tmp_file, grid_file):
        """Read and write header of each node.

        Args:
            tmp_file: Input grid file (tmp_grid).
            grid_file: Output grid file (final_grid).

        Returns:
            Int: Binary value of the is_leaf parameter.
        """
        is_leaf = tmp_file.read(2)
        if len(is_leaf) == 2:
            grid_file.write(is_leaf)
            level = tmp_file.read(2)
            grid_file.write(level)
        elif len(is_leaf) != 0:
            raise ValueError('Problem with is_leaf in grid normalization!')
        return is_leaf


class Spherical(Grid):
    """This class creates spherical grids based on the models defined in model.py.
    """

    def __init__(self, model, path, parse_args):
        """Initialisation of grid parameters.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
            path : Handles file input/output and all necessary paths.
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
        """
        Grid.__init__(self, model, path, parse_args)

    def init_root(self):
        """Initialise the root node.

        Return:
            Root node of the spherical grid.
        """
        # Create root node
        root = Node('spherical')
        # Set root node position
        root.parameter['position'] = [0., 0., 0.]
        # Set width between node borders.
        root.parameter['extent'] = [self.model.spherical_parameter['inner_radius'],
                                    self.model.spherical_parameter['outer_radius'],
                                    0, np.pi, 0., 2. * np.pi]
        # Set node volume
        root.parameter['volume'] = self.get_volume(node=root)
        return root

    def create_grid(self, grid_file, root):
        """Create a spherical grid and calculate the total mass of the nodes.

        Args:
            grid_file: Input grid file (tmp_grid).
            root: Instance of spherical grid root node.
        """
        #: Parameter from the chosen model used for the grid creation
        sp_param = self.model.spherical_parameter

        #: Array of radius values
        if sp_param['sf_r'] == 0:
            if sp_param['radius_list'][0] != sp_param['inner_radius'] or \
                    sp_param['radius_list'][-1] != sp_param['outer_radius']:
                raise ValueError(
                    'radius_list does not agree with the inner and outer grid borders!')
            radius_list = sp_param['radius_list']
            sp_param['n_r'] = len(radius_list) - 1
        elif sp_param['sf_r'] == 1:
            radius_list = self.math.sin_list(
                sp_param['inner_radius'], sp_param['outer_radius'], sp_param['n_r'])
        elif sp_param['sf_r'] > 1:
            radius_list = self.math.exp_list(sp_param['inner_radius'],
                                             sp_param['outer_radius'], sp_param['n_r'], sp_param['sf_r'])
        else:
            radius_list = self.math.lin_list(
                sp_param['inner_radius'], sp_param['outer_radius'], sp_param['n_r'])

        if sp_param['split_first_cell'] > 1 and sp_param['sf_r'] != 0:
            sp_param['radius_list'] = np.hstack((np.linspace(radius_list[0], radius_list[1],
                                                             sp_param['split_first_cell'] + 1), radius_list[2:])).ravel()
            radius_list = sp_param['radius_list']
            sp_param['sf_r'] = 0
            sp_param['n_r'] = len(radius_list) - 1
            sp_param['split_first_cell'] = 1

        #: Array of phi values
        if sp_param['sf_ph'] == 0:
            phi_list = sp_param['phi_list']
            sp_param['n_ph'] = len(phi_list) - 1
        else:
            phi_list = self.math.lin_list(0., 2. * np.pi, sp_param['n_ph'])

        #: Array of theta values
        if sp_param['sf_th'] == 0:
            theta_list = sp_param['theta_list']
            sp_param['n_th'] = len(theta_list) - 1
        elif sp_param['sf_th'] == 1:
            theta_list = self.math.sin_list(0, np.pi, sp_param['n_th'])
        elif sp_param['sf_th'] > 1:
            theta_list = self.math.exp_list_sym(0, np.pi, sp_param['n_th'], sp_param['sf_th'])
        else:
            theta_list = self.math.lin_list(0, np.pi, sp_param['n_th'])

        grid_file.write(struct.pack('d', sp_param['inner_radius']))
        grid_file.write(struct.pack('d', sp_param['outer_radius']))
        grid_file.write(struct.pack('H', sp_param['n_r']))
        grid_file.write(struct.pack('H', sp_param['n_ph']))
        grid_file.write(struct.pack('H', sp_param['n_th']))
        # f_radius
        grid_file.write(struct.pack('d', sp_param['sf_r']))
        # f_phi
        grid_file.write(struct.pack('d', sp_param['sf_ph']))
        # f_theta
        grid_file.write(struct.pack('d', sp_param['sf_th']))
        # Write radius list if custom
        if sp_param['sf_r'] == 0:
            for tmp_rad in radius_list[1:-1]:
                grid_file.write(struct.pack('d', tmp_rad))
        # Write phi list if custom
        if sp_param['sf_ph'] == 0:
            for tmp_phi in phi_list[1:-1]:
                grid_file.write(struct.pack('d', tmp_phi))
        # Write theta list if custom
        if sp_param['sf_th'] == 0:
            for tmp_theta in theta_list[1:-1]:
                grid_file.write(struct.pack('d', tmp_theta))

        # Calculate the total number of cells
        nr_cells = 0
        for i_r in range(sp_param['n_r']):
            nr_cells += sp_param['n_ph'] * sp_param['n_th']

        i_node = 0
        for i_r in range(sp_param['n_r']):
            for i_p in range(sp_param['n_ph']):
                for i_t in range(sp_param['n_th']):
                    stdout.write('--- Generate spherical grid: ' +
                                 str(round(100.0 * i_node / nr_cells, 3)) + ' %      \r')
                    stdout.flush()
                    # Calculate the cell midpoint in spherical coordinates
                    spherical_coord = np.zeros(3)
                    spherical_coord[0] = (
                        radius_list[i_r] + radius_list[i_r + 1]) / 2.
                    spherical_coord[1] = (
                        theta_list[i_t] + theta_list[i_t + 1]) / 2.
                    spherical_coord[2] = (
                        phi_list[i_p] + phi_list[i_p + 1]) / 2.
                    # Convert the spherical coordinate into cartesian node position
                    position = self.math.spherical_to_cartesian(
                        spherical_coord)
                    node = Node('spherical')
                    node.parameter['position'] = position
                    node.parameter['extent'] = [radius_list[i_r], radius_list[i_r + 1],
                                                theta_list[i_t], theta_list[i_t + 1],
                                                phi_list[i_p], phi_list[i_p + 1]]
                    node.parameter['volume'] = self.get_volume(node=node)
                    self.write_node_data(grid_file=grid_file, node=node, data_type='d',
                                         cell_IDs=[i_r, i_t, i_p])
                    self.update_mass_measurement(node=node)
                    del node
                    i_node += 1
        node = Node('spherical')
        node.parameter['position'] = [0., 0., 0.]
        node.parameter['extent'] = [0., radius_list[0],
                                    0, np.pi,
                                    0, 2. * np.pi]
        node.parameter['volume'] = self.get_volume(node=node)
        self.write_node_data(grid_file=grid_file, node=node, data_type='d',
                             cell_IDs=[-1, -1, -1])
        self.update_mass_measurement(node=node)

    @staticmethod
    def get_volume(node):
        """Calculate the volume of a spherical node.

        Args:
            node: Instance of spherical node.

        Returns:
            Volume of the node
        """
        volume = (node.parameter['extent'][1] ** 3 - node.parameter['extent'][0] ** 3) * \
                 (np.cos(node.parameter['extent'][2]) - np.cos(node.parameter['extent'][3])) * \
                 (node.parameter['extent'][5] -
                  node.parameter['extent'][4]) / 3.
        return volume

    def normalize_density(self, tmp_file, grid_file):
        """Read the temporary spherical grid and normalize the
        density to the total model mass.

        Args:
            tmp_file: Input grid file (tmp_grid).
            grid_file: Output grid file (final_grid).
        """
        #: Parameter from the chosen model used for the grid creation
        sp_param = self.model.spherical_parameter

        self.read_write_header(tmp_file=tmp_file, grid_file=grid_file)
        for i_r in range(sp_param['n_r']):
            for i_t in range(sp_param['n_th']):
                for i_p in range(sp_param['n_ph']):
                    self.read_write_node_data(
                        tmp_file=tmp_file, grid_file=grid_file, data_type='d')
        self.read_write_node_data(
            tmp_file=tmp_file, grid_file=grid_file, data_type='d')

    def read_write_header(self, tmp_file, grid_file):
        """Read and write spherical header from binary file.

        Args:
            tmp_file: Input grid file (tmp_grid).
            grid_file: Output grid file (final_grid).
        """
        grid_file.write(tmp_file.read(2))
        grid_file.write(tmp_file.read(2))
        for i in range(self.data_length):
            grid_file.write(tmp_file.read(2))
        grid_file.write(tmp_file.read(8))
        grid_file.write(tmp_file.read(8))
        # Get number of radial cells
        n_r = tmp_file.read(2)
        grid_file.write(n_r)
        # Get number of phi cells
        n_ph = tmp_file.read(2)
        grid_file.write(n_ph)
        # Get number of theta cells
        n_th = tmp_file.read(2)
        grid_file.write(n_th)
        # Get step width factors (zero for custom)
        sf_r = tmp_file.read(8)
        grid_file.write(sf_r)
        sf_ph = tmp_file.read(8)
        grid_file.write(sf_ph)
        sf_th = tmp_file.read(8)
        grid_file.write(sf_th)
        # Write the custom cell distribution if chosen
        if struct.unpack('d', sf_r)[0] == 0:
            for i_r in range(struct.unpack('H', n_r)[0] - 1):
                grid_file.write(tmp_file.read(8))
        if struct.unpack('d', sf_ph)[0] == 0:
            for i_ph in range(struct.unpack('H', n_ph)[0] - 1):
                grid_file.write(tmp_file.read(8))
        if struct.unpack('d', sf_th)[0] == 0:
            for i_th in range(struct.unpack('H', n_th)[0] - 1):
                grid_file.write(tmp_file.read(8))

    def write_other_grid(self, tmp_file, code_name):
        """Generate a grid for other RT codes as well.

        Args:
            tmp_file: Input grid file (tmp_grid).
            code_name (str): Name of the other RT code.
        """
        #: Parameter from the chosen model used for the grid creation
        sp_param = self.model.spherical_parameter
        from astropy.io import fits

        if code_name == 'mcfost':
            if sp_param['n_ph'] == 1:
                tbldata = np.zeros(
                    (int(sp_param['n_th'] / 2), sp_param['n_r']))
            else:
                tbldata = np.zeros((sp_param['n_ph'], int(
                    sp_param['n_th'] / 2), sp_param['n_r']))
            tmp_file.read(2 + 2 + 2 * self.data_length +
                          8 + 8 + 2 + 2 + 2 + 8 + 8)
            for i_r in range(sp_param['n_r']):
                for i_p in range(sp_param['n_ph']):
                    for i_t in range(sp_param['n_th']):
                        density = tmp_file.read(8)
                        if np.sum(self.total_gas_mass) > 0. and i_t >= int(sp_param['n_th'] / 2.0):
                            i_t_tmp = i_t - int(sp_param['n_th'] / 2.0)
                            data = struct.unpack('d', density)[0] * \
                                (self.model.parameter['gas_mass'] /
                                 np.sum(self.total_gas_mass))
                            if sp_param['n_ph'] == 1:
                                tbldata[i_t_tmp, i_r] = data
                            else:
                                tbldata[i_p, i_t_tmp, i_r] = data
                        for i in range(self.data_length - 1):
                            tmp_file.read(8)
            hdu = fits.PrimaryHDU(tbldata)
            hdulist = fits.HDUList([hdu])
            hdulist.writeto(
                self.path['model'] + code_name + '_grid.fits', overwrite=True)


class Cylindrical(Grid):
    """This class creates cylindrical grids based on the models defined in model.py.
    """

    def __init__(self, model, path, parse_args):
        """Initialisation of grid parameters.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
            path : Handles file input/output and all necessary paths.
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
        """
        Grid.__init__(self, model, path, parse_args)

    def init_root(self):
        """Initialise the root node.

        Return:
            Root node of the cylindrical grid.
        """
        # Create root node
        root = Node('cylindrical')
        # Set root node position
        root.parameter['position'] = [0., 0., 0.]
        # Set width between node borders.
        root.parameter['extent'] = [self.model.cylindrical_parameter['inner_radius'],
                                    self.model.cylindrical_parameter['outer_radius'],
                                    0., 2. * np.pi,
                                    -self.model.cylindrical_parameter['z_max'],
                                    self.model.cylindrical_parameter['z_max']]
        # Set node volume
        root.parameter['volume'] = self.get_volume(node=root)
        return root

    def create_grid(self, grid_file, root):
        """Create an cylindrical grid and calculate the total mass of the nodes.

        Args:
            grid_file: Input grid file (tmp_grid).
            root: Instance of cylindrical grid root node.
        """
        #: Parameter from the chosen model used for the grid creation
        cy_param = self.model.cylindrical_parameter

        #: Array of radius values
        if cy_param['sf_r'] == 0:
            if cy_param['radius_list'][0] != cy_param['inner_radius'] or \
                    cy_param['radius_list'][-1] != cy_param['outer_radius']:
                raise ValueError(
                    'radius_list does not agree with the inner and outer grid borders!')
            radius_list = cy_param['radius_list']
            cy_param['n_r'] = len(radius_list) - 1
        elif cy_param['sf_r'] == 1:
            radius_list = self.math.sin_list(
                cy_param['inner_radius'], cy_param['outer_radius'], cy_param['n_r'])
        elif cy_param['sf_r'] > 1:
            radius_list = self.math.exp_list(cy_param['inner_radius'],
                                             cy_param['outer_radius'], cy_param['n_r'], cy_param['sf_r'])
        else:
            radius_list = self.math.lin_list(
                cy_param['inner_radius'], cy_param['outer_radius'], cy_param['n_r'])

        # Add refinement to innermost cell
        if cy_param['split_first_cell'] > 1 and cy_param['sf_r'] != 0:
            cy_param['radius_list'] = np.hstack((np.linspace(radius_list[0], radius_list[1],
                                                             cy_param['split_first_cell'] + 1), radius_list[2:])).ravel()
            radius_list = cy_param['radius_list']
            cy_param['sf_r'] = 0
            cy_param['n_r'] = len(radius_list) - 1
            cy_param['split_first_cell'] = 1

        # Prepare n_ph for grid creation
        if isinstance(cy_param['n_ph'], int):
            cy_param['n_ph'] = [cy_param['n_ph']] * cy_param['n_r']

        #: Array of phi values
        if cy_param['sf_ph'] == 0:
            if len(cy_param['phi_list']) > 0:
                if cy_param['phi_list'][0] != 0 or \
                        cy_param['phi_list'][-1] != 2. * np.pi:
                    raise ValueError('phi_list does not fullfil a full cicle!')
                phi_list = np.array([cy_param['phi_list']
                                     for i_r in range(cy_param['n_r'])])
                cy_param['n_ph'] = [len(phi_list[0])-1] * cy_param['n_r']
            else:
                raise ValueError(
                    'Cell distriution in phi-direction not understood!')
        elif cy_param['sf_ph'] == -1:
            phi_list = [self.math.lin_list(0., 2. * np.pi, n_ph)
                        for n_ph in cy_param['n_ph']]
        else:
            phi_list = np.array([self.math.lin_list(0., 2. * np.pi, cy_param['n_ph'][0])
                                 for i_r in range(cy_param['n_r'])])

        #: Array of z values
        if cy_param['sf_z'] == 0:
            if len(cy_param['z_list']) > 0:
                if cy_param['z_list'][0] != -cy_param['z_max'] or \
                        cy_param['z_list'][-1] != cy_param['z_max']:
                    raise ValueError(
                        'z_list does not agree with the inner and outer grid borders!')
                z_list = np.array([cy_param['z_list']
                                   for i_r in range(cy_param['n_r'])])
                cy_param['n_z'] = len(z_list[0]) - 1
            else:
                raise ValueError(
                    'Cell distribution in z-direction not understood!')
        elif cy_param['sf_z'] == -1:
            z_max_tmp = [self.model.get_dz(
                radius_list[i_r]) * cy_param['n_z'] / 2. for i_r in range(cy_param['n_r'])]
            z_list = np.array(
                [self.math.lin_list(-zmax, zmax, cy_param['n_z']) for zmax in z_max_tmp])
        elif cy_param['sf_z'] == 1.0:
            z_list = np.array([self.math.sin_list(-cy_param['z_max'], cy_param['z_max'], cy_param['n_z'])
                               for i_r in range(cy_param['n_r'])])
        elif cy_param['sf_z'] > 1.0:
            z_list = np.array([self.math.exp_list_sym(-cy_param['z_max'], cy_param['z_max'],
                                                      cy_param['n_z'], cy_param['sf_z']) for i_r in range(cy_param['n_r'])])
        else:
            z_list = np.array([self.math.lin_list(-cy_param['z_max'], cy_param['z_max'], cy_param['n_z'])
                               for i_r in range(cy_param['n_r'])])

        grid_file.write(struct.pack('d', cy_param['inner_radius']))
        grid_file.write(struct.pack('d', cy_param['outer_radius']))
        grid_file.write(struct.pack('d', cy_param['z_max']))
        grid_file.write(struct.pack('H', cy_param['n_r']))
        grid_file.write(struct.pack('H', cy_param['n_ph'][0]))
        grid_file.write(struct.pack('H', cy_param['n_z']))
        # f_rho
        grid_file.write(struct.pack('d', cy_param['sf_r']))
        # f_phi
        grid_file.write(struct.pack('d', cy_param['sf_ph']))
        # f_z
        grid_file.write(struct.pack('d', cy_param['sf_z']))
        # Write radius list if custom
        if cy_param['sf_r'] == 0:
            for tmp_rho in radius_list[1:-1]:
                grid_file.write(struct.pack('d', tmp_rho))
        # Write phi list if custom
        if cy_param['sf_ph'] == 0:
            for tmp_phi in phi_list[0][1:-1]:
                grid_file.write(struct.pack('d', tmp_phi))
        elif cy_param['sf_ph'] == -1:
            for i_r in range(cy_param['n_r']):
                grid_file.write(struct.pack('H', cy_param['n_ph'][i_r]))
        # Write dz or z, if chosen
        if cy_param['sf_z'] == 0:
            for tmp_z in z_list[0][1:-1]:
                grid_file.write(struct.pack('d', tmp_z))
        elif cy_param['sf_z'] == -1:
            for rho_tmp in radius_list[:-1]:
                grid_file.write(struct.pack('d', self.model.get_dz(rho_tmp)))

        # Calculate the total number of cells
        nr_cells = 0
        for i_r in range(cy_param['n_r']):
            nr_cells += cy_param['n_ph'][i_r] * cy_param['n_z']

        i_node = 0
        for i_r in range(cy_param['n_r']):
            for i_p in range(cy_param['n_ph'][i_r]):
                for i_z in range(cy_param['n_z']):
                    stdout.write('--- Generate cylindrical grid: ' +
                                 str(round(100.0 * i_node / nr_cells, 3)) + ' %      \r')
                    stdout.flush()
                    # Calculate the cell midpoint in cylindrical coordinates
                    cylindrical_coord = np.array([
                        (radius_list[i_r] + radius_list[i_r + 1]) / 2.,
                        (phi_list[i_r][i_p] + phi_list[i_r][i_p + 1]) / 2.,
                        (z_list[i_r][i_z] + z_list[i_r][i_z + 1]) / 2.
                    ])
                    # Convert the cylindrical coordinate into cartesian node position
                    node = Node('cylindrical')
                    node.parameter['position'] = self.math.cylindrical_to_cartesian(
                        cylindrical_coord)
                    node.parameter['extent'] = [radius_list[i_r], radius_list[i_r + 1],
                                                phi_list[i_r][i_p], phi_list[i_r][i_p + 1],
                                                z_list[i_r][i_z], z_list[i_r][i_z + 1]]
                    node.parameter['volume'] = self.get_volume(node=node)
                    self.write_node_data(grid_file=grid_file, node=node, data_type='d',
                                         cell_IDs=[i_r, i_p, i_z])
                    self.update_mass_measurement(node=node)
                    del node
                    i_node += 1

        for i_z in range(cy_param['n_z']):
            node = Node('cylindrical')
            node.parameter['position'] = [
                0., 0., (z_list[0][i_z] + z_list[0][i_z + 1]) / 2.]
            node.parameter['extent'] = [0., radius_list[0],
                                        phi_list[0][0], phi_list[0][cy_param['n_ph'][0]], z_list[0][i_z], z_list[0][i_z + 1]]
            node.parameter['volume'] = self.get_volume(node=node)
            self.write_node_data(grid_file=grid_file, node=node,
                                 data_type='d', cell_IDs=[-1, -1, i_z])
            self.update_mass_measurement(node=node)
            del node

    @staticmethod
    def get_volume(node):
        """Calculate the volume of a cylindrical node.

        Args:
            node: Instance of cylindrical node.

        Returns:
            Volume of the node
        """
        volume = (node.parameter['extent'][1] ** 2 - node.parameter['extent'][0] ** 2) * \
                 (node.parameter['extent'][3] - node.parameter['extent'][2]) * \
                 (node.parameter['extent'][5] -
                  node.parameter['extent'][4]) / 2.
        return volume

    def normalize_density(self, tmp_file, grid_file):
        """Read the temporary cylindrical grid and normalize the
        density to the total model mass.

        Args:
            tmp_file: Input grid file (tmp_grid).
            grid_file: Output grid file (final_grid).
        """
        #: Parameter from the chosen model used for the grid creation
        cy_param = self.model.cylindrical_parameter

        self.read_write_header(tmp_file=tmp_file, grid_file=grid_file)
        for i_r in range(cy_param['n_r']):
            for i_p in range(cy_param['n_ph'][i_r]):
                for i_z in range(cy_param['n_z']):
                    self.read_write_node_data(
                        tmp_file=tmp_file, grid_file=grid_file, data_type='d')
        for i_z in range(cy_param['n_z']):
            self.read_write_node_data(
                tmp_file=tmp_file, grid_file=grid_file, data_type='d')

    def read_write_header(self, tmp_file, grid_file):
        """Read and write spherical header from binary file.

        Args:
            tmp_file: Input grid file (tmp_grid).
            grid_file: Output grid file (final_grid).
        """
        grid_file.write(tmp_file.read(2))
        grid_file.write(tmp_file.read(2))
        for i in range(self.data_length):
            grid_file.write(tmp_file.read(2))
        grid_file.write(tmp_file.read(8))
        grid_file.write(tmp_file.read(8))
        grid_file.write(tmp_file.read(8))
        # Get number of radial cells
        n_r = tmp_file.read(2)
        grid_file.write(n_r)
        # Get number of phi cells
        n_ph = tmp_file.read(2)
        grid_file.write(n_ph)
        # Get number of theta cells
        n_z = tmp_file.read(2)
        grid_file.write(n_z)
        # Get step width factors (zero for custom)
        sf_r = tmp_file.read(8)
        grid_file.write(sf_r)
        sf_ph = tmp_file.read(8)
        grid_file.write(sf_ph)
        sf_z = tmp_file.read(8)
        grid_file.write(sf_z)
        if struct.unpack('d', sf_r)[0] == 0:
            for i_r in range(struct.unpack('H', n_r)[0] - 1):
                grid_file.write(tmp_file.read(8))
        if struct.unpack('d', sf_ph)[0] == 0:
            for i_ph in range(struct.unpack('H', n_ph)[0] - 1):
                grid_file.write(tmp_file.read(8))
        elif struct.unpack('d', sf_ph)[0] == -1:
            for i_r in range(struct.unpack('H', n_r)[0]):
                grid_file.write(tmp_file.read(2))
        if struct.unpack('d', sf_z)[0] == 0:
            for i_z in range(struct.unpack('H', n_z)[0] - 1):
                grid_file.write(tmp_file.read(8))
        elif struct.unpack('d', sf_z)[0] == -1:
            for i_r in range(struct.unpack('H', n_r)[0]):
                grid_file.write(tmp_file.read(8))


class Node:
    """The Node class includes the information of one node in the grid.
    """

    def __init__(self, grid_type):
        """Initialisation of node parameters.

        Args:
            grid_type (str): Name of the grid type.
        """
        #: dict: Includes parameters of a node
        self.parameter = {
            'position': [0., 0., 0.],
            'extent': [0., 0., 0.],
            'volume': 0.,
            'gas_density': 0.,
        }

        #: Parent node
        self.parent = None

        if grid_type == 'octree':
            #: Children nodes
            self.children = [None, None, None, None, None, None, None, None]
            # Index for the level depth of the node
            self.parameter['level'] = 0
            # Index of the node (from 0 to 7)
            self.parameter['index'] = 0
            # Is the node a leaf with data?
            self.parameter['is_leaf'] = False
            # Amount of bytes used by a node
            self.parameter['byte_count'] = 0
        elif grid_type in ['spherical', 'cylindrical']:
            #: Children nodes
            self.children = []
