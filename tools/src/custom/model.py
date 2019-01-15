#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from modules.math import Math
from modules.base import Model

"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_model_dict(dictionary):
    model_dict = {
        'cube': Cube,
        'filament': Filament,
        'galaxy': Galaxy,
        'mhd_binary': MhdFlock,
        'gg_tau_disk': GGTauDisk,
        'hd97048': HD97048,
        'hd169142': HD169142,
        'themis_disk': ThemisDisk,
        'multi_disk': MultiDisk,
        'test': TestModel,
        'custom': CustomModel,
    }
    dictionary.update(model_dict)


class CustomModel(Model):
    """Change this to the model you want to use.
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        # Set parameters of the custom model (see parent Model class for all available options)
        self.parameter['distance'] = 140.0 * self.math.const['pc']
        self.parameter['grid_type'] = 'spherical'
        self.parameter['inner_radius'] = 0.1 * self.math.const['au']
        self.parameter['outer_radius'] = 100.0 * self.math.const['au']
        self.parameter['gas_mass'] = 1e-2 * self.math.const['M_sun']
        # Define which other choice are default for this model
        self.parameter['background_source'] = 'bg_plane'
        self.parameter['stellar_source'] = 't_tauri'
        self.parameter['dust_composition'] = 'mrn'
        self.parameter['gas_species'] = 'oh'
        self.parameter['detector'] = 'cartesian'

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters and update
        disk parameter that depend on other parameter.
        """
        # Use extra_parameter to adjust the model without changing the model.py file

    def gas_density_distribution(self):
        """Define here your routine to calculate the density at a given position
        in the model space.

        Notes:
            Use 'self.position' to calculate the quantity depending on position.

            Define also the following routines if necessary:
                dust_density_distribution(self), gas_temperature(self),
                dust_temperature(self), velocity_field(
                    self), magnetic_field(self),
                dust_id(self), dust_min_size(self), dust_max_size(self)

            xyz_density_distribution can return a density or 2D list of densities.
                - The first dimension is used to define multiple density distributions
                    for different dust compositions (see CustomDust in dust.py for explanation)
                - With the second dimension, multiple regions of the density distribution
                    of the same dust composition can be normalized individually to different total masses.
                - The self.parameter['gas_mass'] needs to have the same dimension and size as the return of this

        Returns:
            float: Gas density at a given position.
        """
        gas_density = 1.0
        # Or a function that depends on the position!
        # See other models to find prewritten distributions (shakura & sunyaev, ...)
        return gas_density


class Cube(Model):
    """A cube model with constant density
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        #: Set parameters of the sphere model
        self.parameter['distance'] = 140.0 * self.math.const['pc']
        # 2.8e-14 * self.math.const['M_sun']
        self.parameter['gas_mass'] = 1e-6 * self.math.const['M_sun']
        self.parameter['outer_radius'] = 100.0 * \
            self.math.const['au']  # 0.5 * self.math.const['au']
        self.parameter['stellar_source'] = 'isrf'
        self.parameter['dust_composition'] = 'silicate_oblate'
        self.parameter['detector'] = 'cartesian'

    def dust_temperature(self):
        """Calculates the dust temperature at a given position.

        Returns:
            float: Dust temperature at a given position.
        """
        dust_temperature = 20.
        return dust_temperature

    def gas_temperature(self):
        """Calculates the dust temperature at a given position.

        Returns:
            float: Dust temperature at a given position.
        """
        gas_temperature = 10.
        return gas_temperature

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Returns:
            float: Gas density at a given position.
        """
        # gas_density = self.math.bonor_ebert_density(self.position, outer_radius=self.parameter['outer_radius'],
        #                                        truncation_radius=1 * self.math.const['au'])
        gas_density = self.math.random_density_distribution(
            self.position, d_exp=2)
        # gas_density = 1.0
        return gas_density

    def magnetic_field(self):
        """Calculates the magnetic field strength at a given position.

        Returns:
            List[float, float, float]: Magnetic field strength at the given
            position.
        """
        return self.math.simple_mag_field(mag_field_strength=1e-10, axis='z')


class Galaxy(Model):
    """A galaxy model.
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        #: Set parameters of the Bok globule model
        self.parameter['distance'] = 100.0 * self.math.const['pc']
        self.parameter['dust_composition'] = 'mrn'
        self.parameter['detector'] = 'cartesian'
        # Density factor to convert MHD simulation from g/cm^3 to kg/m^3
        self.conv_parameter['conv_dens'] = 6195019.204559535
        # Lengths factor to convert MHD simulation from cm to m
        self.conv_parameter['conv_len'] = 1e-2
        # Magnetic field factor to convert MHD simulation from Gauss to Tesla
        self.conv_parameter['conv_mag'] = 1e-4
        # Velocity factor to convert MHD simulation from cm/s to m/s
        self.conv_parameter['conv_vel'] = 1e-2


class TestModel(Model):
    """The test grid model.
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        #: Set parameters of the disk model
        self.parameter['distance'] = 140.0 * self.math.const['pc']
        self.parameter['inner_radius'] = 0.1 * self.math.const['au']
        self.parameter['outer_radius'] = 100.0 * self.math.const['au']
        self.parameter['grid_type'] = 'spherical'
        self.spherical_parameter['n_r'] = 100
        self.spherical_parameter['n_th'] = 91
        self.spherical_parameter['n_ph'] = 1
        self.spherical_parameter['sf_r'] = 1.03
        # self.parameter['gas_mass'] = [[1.22138e-07 * self.math.const['M_sun']],
        #                              [8.77862e-07 * self.math.const['M_sun']]]
        self.parameter['gas_mass'] = [[1e-6 * self.math.const['M_sun']],
                                      [1e-5 * self.math.const['M_sun']]]
        self.parameter['stellar_source'] = 't_tauri'
        self.parameter['dust_composition'] = 'silicate'
        self.parameter['detector'] = 'cartesian'
        self.parameter['variable_dust'] = True
        self.parameter['variable_size_limits'] = True

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Returns:
            float: Gas density at a given position.
        """
        gas_density1 = self.math.sphere_density(self.position,
                                                inner_radius=self.parameter['inner_radius'],
                                                outer_radius=self.parameter['outer_radius'])
        # gas_density2 = self.math.sphere_density(self.position,
        #                                         inner_radius=100. *
        #                                         self.parameter['inner_radius'],
        #                                         outer_radius=self.parameter['outer_radius'])
        # return [[gas_density1], [gas_density2]]
        if np.linalg.norm(self.position) < 0.5 * self.spherical_parameter['outer_radius']:
            return [[gas_density1], [0]]
        else:
            return [[0], [gas_density1]]

    def dust_id(self):
        """Calculates the dust ID depending on the position in the grid.
        The dust ID is related to the dust composition. With this, one can
        change the dust properties inside the disk.

        Returns:
            int: dust ID.
        """
        if np.linalg.norm(self.position) < 0.5 * self.spherical_parameter['outer_radius']:
            dust_id = 0
        else:
            dust_id = 1
        return dust_id

    def dust_min_size(self):
        """Calculates the minimum dust grain size depending on the position in the grid.
        This overwrites the global minimum grain size, but has no effect if it is smaller than it.

        Returns:
            float: minimum grain size
        """
        dust_min_size = 5e-9
        return dust_min_size

    def dust_max_size(self):
        """Calculates the maximum dust grain size depending on the position in the grid.
        This overwrites the global maximum grain size, but has no effect if it is larger than it.

        Returns:
            float: maximum grain size
        """
        if np.linalg.norm(self.position) < 0.5 * self.spherical_parameter['outer_radius']:
            dust_max_size = 0.25e-6
        else:
            dust_max_size = 0.001
        return dust_max_size


class Filament(Model):
    """A sphere model with constant density
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        #: Set parameters of the sphere model
        self.parameter['distance'] = 140.0 * self.math.const['pc']
        self.parameter['outer_radius'] = 10 * self.math.const['pc']
        self.parameter['gas_mass'] = 100.0 * self.math.const['M_sun']
        self.parameter['detector'] = 'cartesian'
        # self.parameter['stellar_source'] = 't_tauri'
        self.parameter['dust_composition'] = 'silicate_oblate'
        self.parameter['detector'] = 'cartesian'

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Returns:
            float: Gas density at a given position.
        """
        density = np.exp(-abs(self.position[0]) / self.math.const['pc'])
        return density


class MhdFlock(Model):
    """The disk model with the Shakura and Sunyaev disk density profile.
    """

    def __init__(self):
        """Initialisation of the model parameters.

        Notes:
            Shakura and Sunyaev (1973)
            Link: http://adsabs.harvard.edu/abs/1973A&A....24..337S
        """
        Model.__init__(self)

        #: Set parameters of the disk model
        self.parameter['distance'] = 100.0 * self.math.const['pc']
        self.parameter['inner_radius'] = 20.0 * self.math.const['au']
        self.parameter['outer_radius'] = 100.0 * self.math.const['au']
        self.parameter['grid_type'] = 'spherical'
        self.spherical_parameter['n_r'] = 256
        self.spherical_parameter['n_th'] = 562
        self.spherical_parameter['n_ph'] = 512
        self.spherical_parameter['sf_r'] = 1.0063066707156978
        self.parameter['gas_mass'] = 1e-2 * self.math.const['M_sun']
        self.parameter['external_input_name'] = 350
        self.parameter['vel_is_speed_of_sound'] = True
        self.parameter['stellar_source'] = 'binary'
        self.parameter['dust_composition'] = 'mrn_oblate'
        self.parameter['detector'] = 'cartesian'


class GGTauDisk(Model):
    """The disk model for GG Tau.

    Notes:
        Based on Duchêne et al. 2004
    """

    def __init__(self):
        """Initialisation of the model parameters.

        Notes:
            - Shakura and Sunyaev 1973 (Disk density distribution)
                "http://adsabs.harvard.edu/abs/1973A&A....24..337S"
            - Duchêne et. al 2004 (Total dust mass and size constaints)
                "https://e-reports-ext.llnl.gov/pdf/304661.pdf"
        """
        Model.__init__(self)

        # ----------------------------------------------
        # ------ Set parameters of the disk model ------
        # ----------------------------------------------
        # Cite: distance (Bertout et al. 1999)
        self.parameter['distance'] = 140. * self.math.const['pc']
        self.parameter['grid_type'] = 'cylindrical'
        # ------ With circumstellar disks -----
        self.parameter['gas_mass'] = np.array(
            [[1.e-3, 1.e-5, 1.e-5, 1.3e-1]]) * self.math.const['M_sun']
        # ---- Without circumstellar disks ----
        # self.parameter['gas_mass'] = 1.3e-1 * self.math.const['M_sun']
        # -------------------------------------
        self.parameter['stellar_source'] = 'gg_tau_stars'
        # Cite: larger grains in cb disk (McCabe et al. 2002)
        self.parameter['dust_composition'] = 'multi_sil'
        self.parameter['detector'] = 'gg_tau'
        self.parameter['variable_dust'] = True
        # ----------------------------------------------
        # --- Parameter for the density distribution ---
        # ----------------------------------------------
        # Position angle of the stars (Ab12 is Ab1 and Ab2)
        self.angle_Aa = 3. / 2. * np.pi
        self.angle_Ab12 = self.angle_Aa + np.pi
        # Half-major axis of the stars (Ab12 is Ab1 and Ab2)
        # Cite: separation Aa and Ab (White et al. 1999)
        rot_angle_2 = 25. + 7.
        self.a_Aab = 36. / 2. * self.math.const['au'] * np.sqrt(
            (np.cos(rot_angle_2 / 180 * np.pi) / np.cos(37 / 180 * np.pi))**2 +
            np.sin(rot_angle_2 / 180 * np.pi)**2)
        self.a_Ab12 = 5. / 2. * self.math.const['au']
        # Inclination of the GG Tau Aa and Ab12 orbits
        self.orbit_inclination = 0 / 180. * np.pi
        # Inclination of the circumstellar disks around the stars
        self.inclination_Aa = 0 / 180. * np.pi
        self.inclination_Ab1 = 0 / 180. * np.pi
        self.inclination_Ab2 = 0. / 180. * np.pi
        self.inclination_rotation_axis = [
            np.cos(25 / 180. * np.pi),
            np.sin(25 / 180. * np.pi),
            0
        ]
        # Vertical shift of the disks
        self.vertical_shift_Aa = 0
        self.vertical_shift_Ab12 = 0

        # ------ Conversion of position angles -----
        # pos = [
        #     np.sin(-25 / 180. * np.pi),
        #     np.cos(-24.9 / 180. * np.pi),
        #     0
        # ]
        # rot_pos = self.math.rotate_coord_system(pos,
        #     rotation_axis=[
        #         np.sin(277. / 180. * np.pi),
        #         np.cos(277. / 180. * np.pi),
        #         0
        #     ],
        #     rotation_angle=37. / 180. * np.pi)
        # print(np.arctan2(rot_pos[0], rot_pos[1]) / np.pi  * 180.)

        # Extend of the circumstellar disks around the stars
        self.inner_radius = 0.15 * self.math.const['au']
        self.outer_radius_Aa = 7. * self.math.const['au']
        self.outer_radius_Ab12 = 2. * self.math.const['au']
        # Parameter of the circumbinary disk
        # Cite: scale height (McCabe et al. 2002)
        # Range: 16 AU, 21 AU, 26 AU, 31 AU
        self.ref_scale_height = 21. * self.math.const['au']
        self.ref_radius = 180. * self.math.const['au']
        self.beta = 1.05
        self.surf_dens_exp = -1.7
        self.cut_off = 2. * self.math.const['au']
        # ----------------------------------------------
        # -------- Distribution of cell borders --------
        # ----------------------------------------------
        self.cylindrical_parameter['n_z'] = 251
        self.cylindrical_parameter['sf_r'] = 0  # Custom radial cell borders
        # Custom number of phi-cells per ring
        self.cylindrical_parameter['sf_ph'] = -1
        # Custom width of z-cell borders per ring
        self.cylindrical_parameter['sf_z'] = -1
        # Radial cells
        r_list_cs_disks = np.linspace(self.a_Aab - 8. * self.math.const['au'],
                                      self.a_Aab + 8. * self.math.const['au'], 150)
        r_list_cb_disk = self.math.exp_list(180. * self.math.const['au'],
                                            260. * self.math.const['au'], 50, 1.03)
        # ------ With circumstellar disks -----
        self.cylindrical_parameter['radius_list'] = np.hstack(
            (r_list_cs_disks, 140 * self.math.const['au'], r_list_cb_disk)).ravel()
        # ---- Without circumstellar disks ----
        # self.cylindrical_parameter['radius_list'] = np.multiply(r_list_cb_disk, self.math.const['au'])
        # -------------------------------------
        # Cite: extent of circumbinary disk 180 AU - 260 AU (Dutrey et al. 2014)
        self.parameter['outer_radius'] = self.cylindrical_parameter['radius_list'][-1]
        self.parameter['inner_radius'] = self.cylindrical_parameter['radius_list'][0]
        # Phi cells
        n_ph_list_1 = [600] * 150
        n_ph_list_2 = [180] * 51
        # ------ With circumstellar disks -----
        self.cylindrical_parameter['n_ph'] = np.hstack(
            (n_ph_list_1, n_ph_list_2)).ravel()
        # ---- Without circumstellar disks ----
        # self.cylindrical_parameter['n_ph'] = n_ph_list_2
        # -------------------------------------

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters.
        """
        # Use extra parameter to vary the disk structure
        self.disk_Aa = False
        self.disk_Ab1 = False
        self.disk_Ab2 = False

        if extra_parameter is not None:
            if len(extra_parameter) == 3:
                if bool(int(extra_parameter[0])):
                    self.disk_Aa = True
                if bool(int(extra_parameter[1])):
                    self.disk_Ab1 = True
                if bool(int(extra_parameter[2])):
                    self.disk_Ab2 = True
            else:
                raise ValueError('Wrong number of extra parameters!')

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.

        Returns:
            float: Gas density at a given position.
        """
        # Calculate cylindrical radius
        radius_cy = np.sqrt(self.position[0] ** 2 + self.position[1] ** 2)

        # --- GG Tau Aa
        if self.disk_Aa and self.a_Aab - self.outer_radius_Aa <= radius_cy <= self.a_Aab + self.outer_radius_Aa:
            # Add inclination
            pos_Aa = self.math.rotate_coord_system([
                self.position[0],
                self.position[1] - self.a_Aab * np.sin(self.angle_Aa),
                self.position[2] + self.a_Aab *
                np.sin(self.orbit_inclination) + self.vertical_shift_Aa
            ], rotation_axis=self.inclination_rotation_axis, rotation_angle=self.inclination_Aa)
            # Calculate the density
            disk_density_Aa = self.math.default_disk_density(pos_Aa, outer_radius=self.outer_radius_Aa,
                                                             inner_radius=self.inner_radius)
        else:
            # Set to zero outside of the disk
            disk_density_Aa = 0.

        # --- GG Tau Ab1 and Ab2
        if self.a_Aab - self.a_Ab12 - self.outer_radius_Ab12 <= radius_cy <= self.a_Aab + self.a_Ab12 + self.outer_radius_Ab12:
            if self.disk_Ab1:
                # --- GG Tau Ab1
                # Add inclination
                pos_Ab1 = self.math.rotate_coord_system([
                    self.position[0] + self.a_Ab12,
                    self.position[1] - self.a_Aab * np.sin(self.angle_Ab12),
                    self.position[2] + self.a_Aab *
                    np.sin(self.orbit_inclination) + self.vertical_shift_Ab12
                ], rotation_axis=self.inclination_rotation_axis, rotation_angle=self.inclination_Ab1)
                # Calculate the density
                disk_density_Ab1 = self.math.default_disk_density(pos_Ab1, outer_radius=self.outer_radius_Ab12,
                                                                  inner_radius=self.inner_radius)
            else:
                disk_density_Ab1 = 0.

            if self.disk_Ab2:
                # --- GG Tau Ab2
                # Add inclination
                pos_Ab2 = self.math.rotate_coord_system([
                    self.position[0] - self.a_Ab12,
                    self.position[1] - self.a_Aab * np.sin(self.angle_Ab12),
                    self.position[2] + self.a_Aab *
                    np.sin(self.orbit_inclination) + self.vertical_shift_Ab12
                ], rotation_axis=self.inclination_rotation_axis, rotation_angle=self.inclination_Ab2)
                # Calculate the density
                disk_density_Ab2 = self.math.default_disk_density(pos_Ab2, outer_radius=self.outer_radius_Ab12,
                                                                  inner_radius=self.inner_radius)
            else:
                disk_density_Ab2 = 0.
        else:
            # Set to zero outside of the disks
            disk_density_Ab1 = 0.
            disk_density_Ab2 = 0.

        # --- GG Tau A CB disk
        if 180. * self.math.const['au'] <= radius_cy <= 260. * self.math.const['au']:
            # Calculate the density
            disk_density = self.math.default_disk_density(self.position, outer_radius=260. * self.math.const['au'],
                                                          inner_radius=180. * self.math.const['au'], ref_scale_height=self.ref_scale_height,
                                                          ref_radius=self.ref_radius, column_dens_exp=self.surf_dens_exp, beta=self.beta)
            # Add exponential decay ath the inner border
            if radius_cy < 190. * self.math.const['au']:
                disk_density *= np.exp(-0.5 * (
                    (190. * self.math.const['au'] - radius_cy) / self.cut_off) ** 2)
        else:
            # Set to zero outside of the disk
            disk_density = 0.

        # Return the densities of each region
        # ------ With circumstellar disks -----
        return [[disk_density_Aa, disk_density_Ab1, disk_density_Ab2, disk_density]]
        # ---- Without circumstellar disks ----
        # return disk_density
        # -------------------------------------

    def dust_id(self):
        """Calculates the dust ID depending on the position in the grid.
        The dust ID is related to the dust composition. With this, one can
        change the dust properties inside the disk.

        Returns:
            int: dust ID.
        """
        # Calculate cylindrical radius
        radius_cy = np.sqrt(self.position[0] ** 2 + self.position[1] ** 2)

        if 180. * self.math.const['au'] <= radius_cy <= 260. * self.math.const['au']:
            dust_id = 1
        else:
            dust_id = 0
        return dust_id

    def scale_height(self, radius):
        """Calculates the scale height at a certain position.

        Args:
            radius (float) : Cylindrical radius of current position

        Returns:
            float: Scale height.
        """
        if(radius <= self.a_Aab + self.outer_radius_Aa):
            scale_height = 0.3 * self.math.const['au']
        else:
            scale_height = self.math.default_disk_scale_height(
                radius, ref_radius=self.ref_radius, ref_scale_height=self.ref_scale_height, beta=self.beta)
        return scale_height


class HD97048(Model):
    """The disk model for HD97048.
    """

    def __init__(self):
        """Initialisation of the model parameters.

        Notes:
            Shakura and Sunyaev (1973)
            Link: http://adsabs.harvard.edu/abs/1973A&A....24..337S
        """
        Model.__init__(self)

        #: Set parameters of the disk model
        self.parameter['distance'] = 185.0 * self.math.const['pc']
        self.parameter['inner_radius'] = 0.3 * self.math.const['au']
        self.parameter['outer_radius'] = 400. * self.math.const['au']
        self.parameter['grid_type'] = 'cylindrical'
        # In the case of a cylindrical grid
        self.cylindrical_parameter['n_r'] = 250
        self.cylindrical_parameter['n_z'] = 61
        self.cylindrical_parameter['n_ph'] = 1
        self.cylindrical_parameter['sf_r'] = 0
        self.cylindrical_parameter['radius_list'] = np.logspace(
            np.log10(self.parameter['inner_radius']),
            np.log10(self.parameter['outer_radius']),
            self.cylindrical_parameter['n_r']
        )
        self.cylindrical_parameter['radius_list'][0] = self.parameter['inner_radius']
        self.cylindrical_parameter['radius_list'][-1] = self.parameter['outer_radius']
        self.cylindrical_parameter['sf_z'] = -1
        self.cylindrical_parameter['split_first_cell'] = 20
        # Define the used sources, dust composition and gas species
        self.parameter['detector'] = 'hd97048'
        self.parameter['stellar_source'] = 'hd97048'
        self.parameter['dust_composition'] = 'olivine_pah'
        # Use multiple dust compositionas depending on the region in the grid
        self.parameter['variable_dust'] = True
        self.use_cont = False

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters.
        """
        # Use the continuum or ring version of the model and set PAH to silicate ration
        if len(extra_parameter) == 1:
            self.use_cont = bool(int(extra_parameter[0]))

        # Radial cell border list
        # list_inner_disk = np.logspace(np.log10(0.3 * self.math.const['au']),
        #                              np.log10(2.6 * self.math.const['au']), 50)

        # Set the gas density
        if self.use_cont:
            self.parameter['gas_mass'] = np.array([
                [1e-4, 4e-3, 5e-3, 1e-1, 0],
                [0, 0, 0, 0, 0.2 * 1e-3],
                [0, 0, 0, 0, 0.8 * 1e-3]
            ]) * self.math.const['M_sun']
        else:
            mf_pah = 1e-3
            self.parameter['gas_mass'] = np.array([
                [1e-4, (1 - mf_pah) * 4e-3, (1 - mf_pah)
                 * 5e-3, (1 - mf_pah) * 1e-1],
                [0, mf_pah * 4e-3, mf_pah * 5e-3, mf_pah * 1e-1]
            ]) * self.math.const['M_sun']

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.


        Returns:
            float: Gas density at a given position.
        """
        # Real zeros in density distribution?
        real_zero = True
        # INNER DISK
        inner_disk = self.math.default_disk_density(self.position,
                                                    beta=1.0, column_dens_exp=-1.0,
                                                    inner_radius=0.3 *
                                                    self.math.const['au'],
                                                    outer_radius=2.6 *
                                                    self.math.const['au'],
                                                    ref_scale_height=5. *
                                                    self.math.const['au'],
                                                    ref_radius=100. * self.math.const['au'], real_zero=real_zero)
        # RINGS
        beta = 1.26
        surf_dens_exp = -0.5
        ref_radius = 100. * self.math.const['au']
        ref_scale_height = 12. * self.math.const['au']
        # RING #1
        ring_1 = self.math.default_disk_density(self.position,
                                                beta=beta, column_dens_exp=surf_dens_exp,
                                                inner_radius=41. *
                                                self.math.const['au'],
                                                outer_radius=51. *
                                                self.math.const['au'],
                                                ref_scale_height=ref_scale_height, ref_radius=ref_radius, real_zero=real_zero)
        # RING #2
        ring_2 = self.math.default_disk_density(self.position,
                                                beta=beta, column_dens_exp=surf_dens_exp,
                                                inner_radius=155. *
                                                self.math.const['au'],
                                                outer_radius=165. *
                                                self.math.const['au'],
                                                ref_scale_height=ref_scale_height, ref_radius=ref_radius, real_zero=real_zero)
        # RING #3
        ring_3 = self.math.default_disk_density(self.position,
                                                beta=beta, column_dens_exp=surf_dens_exp, tappered_gamma=-0.0,
                                                inner_radius=269. *
                                                self.math.const['au'],
                                                outer_radius=400. *
                                                self.math.const['au'],
                                                ref_scale_height=ref_scale_height, ref_radius=ref_radius, real_zero=real_zero)
        # PAH continuum
        if self.use_cont:
            pah_cont = self.math.default_disk_density(self.position,
                                                      beta=beta, column_dens_exp=surf_dens_exp,
                                                      inner_radius=41 *
                                                      self.math.const['au'],
                                                      outer_radius=400 *
                                                      self.math.const['au'],
                                                      ref_scale_height=ref_scale_height, ref_radius=ref_radius, real_zero=real_zero)
            return [[inner_disk, ring_1, ring_2, ring_3, 0],
                    [0, 0, 0, 0, pah_cont],
                    [0, 0, 0, 0, pah_cont], ]
        return [[inner_disk, ring_1, ring_2, ring_3], [0, ring_1, ring_2, ring_3]]

    def scale_height(self, radius):
        """Calculates the scale height at a certain position.

        Args:
            radius (float) : Cylindrical radius of current position

        Returns:
            float: Scale height.
        """
        if(0.3 * self.math.const['au'] <= radius <= 2.6 * self.math.const['au']):
            scale_height = self.math.default_disk_scale_height(radius,
                                                               ref_radius=100. *
                                                               self.math.const['au'],
                                                               ref_scale_height=5. * self.math.const['au'], beta=1.)
        else:
            scale_height = self.math.default_disk_scale_height(radius,
                                                               ref_radius=100. *
                                                               self.math.const['au'],
                                                               ref_scale_height=12. * self.math.const['au'], beta=1.26)
        return scale_height


class HD169142(Model):
    """A disk model of the disk around HD169142.
    """

    def __init__(self):
        """Initialisation of the model parameters.
        """
        Model.__init__(self)

        #: Set parameters of the disk model
        # Cite: 117 +- 4 pc (Grady et al. 2007; Manoj et al. 2006; Gaia Collaboration et al. 2018)
        self.parameter['distance'] = 117.0 * self.math.const['pc']
        # Calc mass_Fraction out of themis density
        self.parameter['grid_type'] = 'cylindrical'
        self.parameter['inner_radius'] = 0.1 * self.math.const['au']
        self.parameter['outer_radius'] = 244.8 * self.math.const['au']
        # Define the used sources, dust composition and gas species
        self.parameter['stellar_source'] = 'hd169142'
        self.parameter['dust_composition'] = 'olivine_pah'
        self.parameter['detector'] = 'hd169142'
        # In the case of a cylindrical grid
        self.cylindrical_parameter['n_r'] = 100
        self.cylindrical_parameter['n_z'] = 142
        self.cylindrical_parameter['n_ph'] = 1
        self.cylindrical_parameter['sf_r'] = 1.04
        # sf_z = -1 is using scale height; sf_z = 1 is sinus;
        # sf_z > 1 is exp with step width sf_z and rest is linear
        self.cylindrical_parameter['sf_z'] = -1
        # Default disk parameter
        self.parameter['ref_radius'] = 100. * self.math.const['au']
        self.parameter['ref_scale_height'] = 10. * self.math.const['au']
        self.parameter['alpha'] = 1.625
        self.parameter['beta'] = 1.125
        # Enable multiple density distributions
        self.parameter['variable_dust'] = True
        # Init new parameter
        self.parameter['model_number'] = None
        self.parameter['r_gap_in'] = 0
        self.parameter['r_gap_out'] = 0
        # Mass fraction
        self.parameter['mass_fraction'] = 0.01

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters.
        """
        # Use extra parameter to vary the disk structure
        if extra_parameter is not None:
            if len(extra_parameter) == 1:
                self.parameter['r_gap_in'] = 5. * self.math.const['au']
                self.parameter['r_gap_out'] = 21.6442 * self.math.const['au']
                # Change mass ratios depending on the chosen model
                self.parameter['model_number'] = int(extra_parameter[0])
                if self.parameter['model_number'] == 1:
                    self.parameter['gas_mass'] = np.array([
                        [7.9099e-7, 5.8142e-3],
                        [0., 0.7733e-2 * 5.8142e-3]
                    ]) * self.math.const['M_sun']
                elif self.parameter['model_number'] == 2:
                    self.parameter['gas_mass'] = np.array([
                        [0, 0.17e-3],
                        [0.63e-3, 0.63e-3],
                        [0.255e-2, 0.255e-2],
                        [0.255e-2, 0.255e-2]
                    ])
                    self.parameter['gas_mass'][:, 0] *= 7.9099e-7 * \
                        self.math.const['M_sun'] / \
                        self.parameter['gas_mass'][:, 0].sum()
                    self.parameter['gas_mass'][:, 1] *= 5.8142e-3 * \
                        self.math.const['M_sun'] / \
                        self.parameter['gas_mass'][:, 1].sum()

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.


        Returns:
            float: Gas density at a given position.
        """
        # Init density 2D list
        density_list = np.zeros((2, 2))
        if self.parameter['model_number'] is not None:
            if self.parameter['model_number'] == 2:
                density_list = np.zeros((4, 2))
        # Calculate cylindrical radius
        radius_cy = np.sqrt(self.position[0] ** 2 + self.position[1] ** 2)
        # Set density according to region
        if radius_cy <= self.parameter['r_gap_in']:
            gas_density = self.math.default_disk_density(
                self.position,
                inner_radius=self.parameter['inner_radius'],
                outer_radius=self.parameter['outer_radius'],
                ref_radius=1. * self.math.const['au'],
                ref_scale_height=0.0346 * self.math.const['au'],
                column_dens_exp=-1.3764, beta=0.7950
            )
            density_list[:, 0] = gas_density
        elif self.parameter['r_gap_out'] <= radius_cy:
            gas_density = self.math.default_disk_density(
                self.position,
                inner_radius=self.parameter['inner_radius'],
                outer_radius=self.parameter['outer_radius'],
                ref_radius=100. * self.math.const['au'],
                ref_scale_height=9.6157 * self.math.const['au'],
                column_dens_exp=-1., beta=1.0683
            )
            density_list[:, 1] = gas_density
        return density_list

    def scale_height(self, radius):
        """Calculates the scale height at a certain position.

        Args:
            radius (float) : Cylindrical radius of current position

        Returns:
            float: Scale height.
        """
        # Set density according to region
        if radius <= self.parameter['r_gap_in']:
            scale_height = self.math.default_disk_scale_height(
                radius, ref_radius=1. * self.math.const['au'],
                ref_scale_height=0.0346 * self.math.const['au'],
                beta=0.7950
            )
        elif self.parameter['r_gap_out'] <= radius:
            scale_height = self.math.default_disk_scale_height(
                radius, ref_radius=100. * self.math.const['au'],
                ref_scale_height=9.6157 * self.math.const['au'],
                beta=1.0683
            )
        else:
            scale_height = self.math.default_disk_scale_height(
                radius, ref_radius=1. * self.math.const['au'],
                ref_scale_height=0.0346 * self.math.const['au'],
                beta=0.7950
            )
        return scale_height


class ThemisDisk(Model):
    """A standard disk model with the THEMIS dust model.
    """

    def __init__(self):
        """Initialisation of the model parameters.

        Notes:
            Shakura and Sunyaev (1973)
            Link: http://adsabs.harvard.edu/abs/1973A&A....24..337S
        """
        Model.__init__(self)

        #: Set parameters of the disk model
        self.parameter['distance'] = 109.0 * self.math.const['pc']
        self.parameter['grid_type'] = 'cylindrical'
        self.parameter['inner_radius'] = 0.2 * self.math.const['au']
        self.parameter['outer_radius'] = 350. * self.math.const['au']
        # Define the used sources, dust composition and gas species
        self.parameter['stellar_source'] = 't_tauri'
        self.parameter['dust_composition'] = 'themis'
        self.parameter['gas_species'] = 'co'
        self.parameter['detector'] = 'cartesian'
        # In the case of a cylindrical grid
        self.cylindrical_parameter['n_r'] = 100
        self.cylindrical_parameter['n_z'] = 142
        self.cylindrical_parameter['n_ph'] = 1
        self.cylindrical_parameter['sf_r'] = 1.03
        # sf_z = -1 is using scale height; sf_z = 1 is sinus;
        # sf_z > 1 is exp with step width sf_z and rest is linear
        self.cylindrical_parameter['sf_z'] = -1
        # Default disk parameter
        self.parameter['ref_radius'] = 100. * self.math.const['au']
        self.parameter['ref_scale_height'] = 10. * self.math.const['au']
        self.parameter['alpha'] = 1.625
        self.parameter['beta'] = 1.125
        # Enable multiple density distributions
        self.parameter['variable_dust'] = True
        # Init new parameter
        self.parameter['model_number'] = 0

    def update_parameter(self, extra_parameter):
        """Use this function to set model parameter with the extra parameters.
        """
        # Use extra parameter to vary the disk structure
        if extra_parameter is not None:
            if len(extra_parameter) == 4:
                self.parameter['ref_radius'] = self.math.parse(
                    extra_parameter[0], 'length')
                self.parameter['ref_scale_height'] = self.math.parse(
                    extra_parameter[1], 'length')
                self.parameter['alpha'] = float(extra_parameter[2])
                self.parameter['beta'] = float(extra_parameter[3])
            elif len(extra_parameter) == 1:
                # Change mass ratios depending on the chosen model
                self.parameter['model_number'] = int(extra_parameter[0])
                if self.parameter['model_number'] == 1:
                    self.parameter['gas_mass'] = np.array(
                        [[0.17e-3], [0.63e-4], [0.255e-3], [0.255e-3]])
                elif self.parameter['model_number'] == 2:
                    self.parameter['gas_mass'] = np.array(
                        [[0.17e-4], [0.63e-4], [0.255e-3], [0.255e-3]])
                elif self.parameter['model_number'] == 3:
                    self.parameter['gas_mass'] = np.array(
                        [[0.8e-4], [0.63e-3], [0.255e-2], [0.255e-2]])
                elif self.parameter['model_number'] in [4, 5]:
                    self.parameter['gas_mass'] = np.array(
                        [[0.17e-3], [0.63e-3], [0.255e-2], [0.255e-2]])
                    if self.parameter['model_number'] == 5:
                        self.tmp_parameter['ignored_gas_density'] = np.zeros(
                            (4, 1))
                elif self.parameter['model_number'] == 6:
                    self.parameter['gas_mass'] = np.array(
                        [[0.17e-3], [0.63e-3], [0.255e-2], [0.255e-2]])
                    self.tmp_parameter['ignored_gas_density'] = np.zeros(
                        (4, 1))
                self.parameter['mass_fraction'] = np.sum(
                    self.parameter['gas_mass'])
                print('--mass_fraction', self.parameter['mass_fraction'])
                self.parameter['gas_mass'] *= 1e-2 * self.math.const['M_sun'] / \
                    np.sum(self.parameter['gas_mass'])

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.


        Returns:
            float: Gas density at a given position.
        """
        gas_density = self.math.default_disk_density(self.position,
                                                     inner_radius=self.parameter['inner_radius'],
                                                     outer_radius=self.parameter['outer_radius'],
                                                     ref_radius=self.parameter['ref_radius'],
                                                     ref_scale_height=self.parameter['ref_scale_height'],
                                                     alpha=self.parameter['alpha'], beta=self.parameter['beta'])
        density_list = np.ones((4, 1)) * gas_density
        if self.parameter['model_number'] in [5, 6]:
            # Get limits for the models
            if self.parameter['model_number'] == 5:
                r_min = 5.
                r_max = 20.
            elif self.parameter['model_number'] == 6:
                r_min = 2.
                r_max = 20.
            # Calculate cylindrical radius
            radius_cy = np.sqrt(self.position[0] ** 2 + self.position[1] ** 2)
            # Set density of larger grains to zero to create a gap
            if r_min * self.math.const['au'] <= radius_cy <= r_max * self.math.const['au']:
                # Add negatively to take it into account for normalization
                # Same as material which was in a disk with the total disk mass but is
                # for instance accreted on a planet or star
                self.tmp_parameter['ignored_gas_density'][1:, 0] -= \
                    density_list[1:, 0] * self.volume
                density_list[1:, 0] = 0.
        return density_list

    def scale_height(self, radius):
        """Calculates the scale height at a certain position.

        Args:
            radius (float) : Cylindrical radius of current position

        Returns:
            float: Scale height.
        """
        scale_height = self.math.default_disk_scale_height(radius,
                                                           ref_radius=self.parameter['ref_radius'],
                                                           ref_scale_height=self.parameter['ref_scale_height'],
                                                           beta=self.parameter['beta'])
        return scale_height


class MultiDisk(Model):
    """The disk model with the Shakura and Sunyaev disk density profile.
    """

    def __init__(self):
        """Initialisation of the model parameters.

        Notes:
            Shakura and Sunyaev (1973)
            Link: http://adsabs.harvard.edu/abs/1973A&A....24..337S
        """
        Model.__init__(self)

        #: Set parameters of the disk model
        self.parameter['distance'] = 140.0 * self.math.const['pc']
        self.parameter['gas_mass'] = np.array(
            [[0.67], [0.33]]) * 1e-4 * self.math.const['M_sun']
        # [[1e-2], [1e-2 * 1e-3]]) * self.math.const['M_sun']
        self.parameter['grid_type'] = 'spherical'
        self.parameter['inner_radius'] = 1. * self.math.const['au']
        self.parameter['outer_radius'] = 300. * self.math.const['au']
        # Define the used sources, dust composition and gas species
        self.parameter['stellar_source'] = 't_tauri'
        self.parameter['dust_composition'] = 'silicate_pah'
        self.parameter['gas_species'] = 'co'
        self.parameter['detector'] = 'cartesian'
        # In the case of a spherical grid
        self.spherical_parameter['n_r'] = 100
        self.spherical_parameter['n_th'] = 181
        self.spherical_parameter['n_ph'] = 1
        self.spherical_parameter['sf_r'] = 1.03
        # sf_th = -1 is linear; sf_th = 1 is sinus; rest is exp with step width sf_th
        self.spherical_parameter['sf_th'] = 1.0
        # In the case of a cylindrical grid
        self.cylindrical_parameter['n_r'] = 100
        self.cylindrical_parameter['n_z'] = 181
        self.cylindrical_parameter['n_ph'] = 1
        self.cylindrical_parameter['sf_r'] = 1.03
        # sf_z = -1 is linear; sf_z = 1 is sinus; rest is exp with step width sf_z
        self.cylindrical_parameter['sf_z'] = 1.0
        # Default disk parameter
        self.ref_radius = 100. * self.math.const['au']
        self.ref_scale_height = 10. * self.math.const['au']
        self.parameter['variable_dust'] = True

    def gas_density_distribution(self):
        """Calculates the gas density at a given position.


        Returns:
            float: Gas density at a given position.
        """
        gas_density = self.math.default_disk_density(self.position,
                                                     inner_radius=self.parameter['inner_radius'],
                                                     outer_radius=self.parameter['outer_radius'])
        return [[gas_density], [gas_density]]

    def dust_temperature(self):
        """Calculates the dust temperature at a given position.

        Returns:
            float: Dust temperature at a given position.
        """
        dust_temperature = 20.
        return dust_temperature
