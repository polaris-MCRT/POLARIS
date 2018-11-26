#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from modules.base import ExternalInput

"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_ext_input_dict(dictionary):
    ext_input_dict = {
        'mflock': MarioFlock,
        'custom': CustomInput,
    }
    dictionary.update(ext_input_dict)


class CustomInput(ExternalInput):
    """Change this to the external input you want to use.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the external input parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        ExternalInput.__init__(self, file_io, parse_args)

        # Index of the MHD simulation
        self.model_index = 320

    def init_data(self):
        """Reads the data of e.g. an MHD simulation
        """
        # Save the data in this list
        self.data = [1.0]

    def gas_density_distribution(self):
        """Define here your routine to get the density at a given position
        in the model space from your read input data.

        Notes:
            Define also the following routines if necessary:
                - dust_density_distribution(self)
                - gas_temperature(self)
                - dust_temperature(self)
                - velocity_field(self)
                - magnetic_field(self)
                - dust_id(self)

        Returns:
            float: Gas density at a given position.
        """
        return self.data[0]


class MarioFlock(ExternalInput):
    """Get model input data from Mario Flock MHD simulations.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the external input parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        ExternalInput.__init__(self, file_io, parse_args)

        # Index of the MHD simulation
        self.model_index = 320

    def init_data(self):
        """Reads the data of an MHD simulation from Mario Flock (Mario.Flock@jpl.nasa.gov)

        Notes:
            unit_dens = 2.987e-11 g/cm^3
            unit_Velocity=2.112e+06 cm/s
            unit_length = 1.496e+13 cm
            unit_mag = 4.092e+01 Gauss

        Returns:
            Output is [0]->Grid (x1,x2,x3)
            [1]->3D Density array
            [2]->3D Velocity array
        """

        # Convert grid density to kg/m^3
        unit_dens = 2.987e-11 * 1e-3 * 1e6
        # Convert grid velocity to m/s
        # unit_vel = 2.112e+06 * 1e-3
        # Convert grid length from AU to m
        unit_length = self.math.const['au']
        # Convert grid magnetic field to Tesla
        unit_mag = 4.092e+01 * 1e-4

        n1 = 256
        n2 = 128
        n3 = 512

        path = self.file_io.path['model'] + 'mhd_flock/'
        filename = path + 'grid.dat'
        file = open(filename, mode='rb')
        bin1 = array.array('f')
        bin2 = array.array('f')
        bin3 = array.array('f')

        bin1.fromfile(file, n1)
        radius = np.reshape(bin1, n1)
        bin2.fromfile(file, n2)
        theta = np.reshape(bin2, n2)
        bin3.fromfile(file, n3)
        phi = np.reshape(bin3, n3)
        file.close()

        filename = path + 'rho.%(number)04d.dbl' % {"number": model_index}
        file = open(filename, mode='rb')
        bin_array = array.array('d')
        bin_array.fromfile(file, n1 * n2 * n3)
        rho = np.reshape(bin_array, (n3, n2, n1))

        filename = path + 'vx1.%(number)04d.dbl' % {"number": model_index}
        file = open(filename, mode='rb')
        bin_array = array.array('d')
        bin_array.fromfile(file, n1 * n2 * n3)
        vr = np.reshape(bin_array, (n3, n2, n1))

        filename = path + 'vx2.%(number)04d.dbl' % {"number": model_index}
        file = open(filename, mode='rb')
        bin_array = array.array('d')
        bin_array.fromfile(file, n1 * n2 * n3)
        vt = np.reshape(bin_array, (n3, n2, n1))

        filename = path + 'vx3.%(number)04d.dbl' % {"number": model_index}
        file = open(filename, mode='rb')
        bin_array = array.array('d')
        bin_array.fromfile(file, n1 * n2 * n3)
        vp = np.reshape(bin_array, (n3, n2, n1))

        filename = path + 'prs.%(number)04d.dbl' % {"number": model_index}
        file = open(filename, mode='rb')
        bin_array = array.array('d')
        bin_array.fromfile(file, n1 * n2 * n3)
        prs = np.reshape(bin_array, (n3, n2, n1))

        cs = np.sqrt(prs / rho)

        filename = path + 'bx1.%(number)04d.dbl' % {"number": model_index}
        file = open(filename, mode='rb')
        bin_array = array.array('d')
        bin_array.fromfile(file, n1 * n2 * n3)
        br = np.reshape(bin_array, (n3, n2, n1))

        filename = path + 'bx2.%(number)04d.dbl' % {"number": model_index}
        file = open(filename, mode='rb')
        bin_array = array.array('d')
        bin_array.fromfile(file, n1 * n2 * n3)
        bt = np.reshape(bin_array, (n3, n2, n1))

        filename = path + 'bx3.%(number)04d.dbl' % {"number": model_index}
        file = open(filename, mode='rb')
        bin_array = array.array('d')
        bin_array.fromfile(file, n1 * n2 * n3)
        bp = np.reshape(bin_array, (n3, n2, n1))

        # Give information to user
        print('--- Read external grid data with index =',
              model_index, 'finished!')

        self.data = ((radius * unit_length, theta, phi), rho.T * unit_dens,
                     vr / cs, vt / cs, vp / cs, br * unit_mag, bt * unit_mag, bp * unit_mag)

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Args:
            position (List[float, float, float]): Position in model space.

        Returns:
            float: Gas density at the given position.
        """
        if None in [self.pos_r, self.pos_t, self.pos_p]:
            gas_density = self.data[1][self.i_r, self.i_t, self.i_p]
        else:
            gas_density = 0.
        return gas_density

    def velocity_field(self):
        """Calculates the velocity at a given position.

        Args:
            position (List[float, float, float]): Position in model space.

        Notes:
            Possible options:
                velocity = [0., 0., 0.]
                velocity = self.math.kepler_rotation(position, stellar_mass=0.7)

        Returns:
            List[float, float, float]: Velocity at the given position.
        """
        if None in [self.pos_r, self.pos_t, self.pos_p]:
            velocity = [0., 0., 0.]
        else:
            velocity_sp = [self.data[2][self.i_p, self.i_t, self.i_r],
                           self.data[3][self.i_p, self.i_t, self.i_r],
                           self.data[4][self.i_p, self.i_t, self.i_r]]
            velocity = self.math.spherical_to_carthesian_direction(
                [self.pos_r, self.pos_t, self.pos_p], velocity_sp)
        return velocity

    def magnetic_field(self):
        """Calculates the magnetic field strength at a given position.

        Args:
            position (List[float, float, float]): position in model space.
            cell_IDs (List): cell_IDs of the cells (alternative to the grid position).
                Spherical -> [i_r, i_t, i_p]
                Cylindrical -> [i_r, i_p, i_z]

        Returns:
            List[float, float, float]: Magnetic field strength at the given
            position.
        """
        if None in [self.pos_r, self.pos_t, self.pos_p]:
            mag_field = [0., 0., 0.]
        else:
            mag_field_sp = [self.data[5][self.i_p, self.i_t, self.i_r],
                            self.data[6][self.i_p, self.i_t, self.i_r],
                            self.data[7][self.i_p, self.i_t, self.i_r]]
            mag_field = self.math.spherical_to_carthesian_direction(
                [self.pos_r, self.pos_t, self.pos_p], mag_field_sp)
        return mag_field
