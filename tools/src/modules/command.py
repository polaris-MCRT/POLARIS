#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np


class CmdPolaris:
    """This class creates ".cmd" files which can be used to execute POLARIS.
    """

    def __init__(self, file_io, model, gas, dust, source,
                 bg_source, dust_source, detector, server, parse_args):
        """Initialisation of the command options.

        Args:
            file_io: Handles file input/output and all
                necessary paths.
            model: Handles the model space including various
                quantities such as the density distribution.
            gas: Handles the gas species with its parameters.
            dust: Handles the dust composition with its parameters.
            source: Handles the radiation of a radiation source
                such as stars or interstellar radiation fields.
            bg_source: Handles the radiation of a background source.
            dust_source: Handles the radiation of the dust grains.
            detector: Handles the position and wavelength
                at which a detector observes
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
        """

        self.file_io = file_io
        self.model = model
        self.gas = gas
        self.dust = dust
        self.source = source
        self.bg_source = bg_source
        self.dust_source = dust_source
        self.detector = detector
        self.server = server
        self.parse_args = parse_args

        # Get math module
        from modules.math import Math
        self.math = Math()

    def write_common_part(self, cmd_file, midplane_points=None, midplane_zoom=None, midplane_3d_param=None,
            midplane_rad_field=False, mass_fraction=None, scattering=None, mu=None, nr_threads=-1, 
            rot_axis_1=None, rot_axis_2=None):
        """Writes commands into the command file.

        Args:
            cmd_file: Instance of the command file.
            midplane_points (int): Number of pixel per axis in midplane map data.
            midplane_zoom (int): Zoom factor per axis in midplane map data.
            midplane_3d_param (List[str, ]): Write midplane fits files as 3D slices?
                (plane e.g. 1 -> xy, nr_slices, z_min, z_max)
            midplane_rad_field (bool): Write radiation field into midplane fits files?
            mass_fraction (float): Fraction of mass between the dust and gas component.
            scattering (str): Name of the scattering method, used for scattering simulations.
                (PH_ISO, PH_MIE, PH_HG)
            mu (float): Average mass of one gas particle.
            nr_threads (int): Number of cpu threads.
            rot_axis_1 (List): First rotation angle (phi) rotates around this vector.
            rot_axis_2 (List): Second rotation angle (theta) rotates around this vector.
        """

        # Overwrite default values with user input
        if self.parse_args.mass_fraction is not None:
            mass_fraction = self.parse_args.mass_fraction
        elif self.dust.parameter['abundances'] is not None:
            mass_fraction = self.dust.parameter['abundances'].sum()
        elif mass_fraction is None:
            mass_fraction = self.model.parameter['mass_fraction']
        if self.parse_args.scattering is not None:
            scattering = self.parse_args.scattering
        elif scattering is None:
            scattering = self.dust.parameter['scattering']
        if self.parse_args.midplane_points is not None:
            midplane_points = self.parse_args.midplane_points
        elif midplane_points is None:
            midplane_points = 256
        if self.parse_args.midplane_zoom is not None:
            midplane_zoom = self.parse_args.midplane_zoom
        elif midplane_zoom is None:
            midplane_zoom = 1
        if self.parse_args.midplane_3d_param is not None:
            midplane_3d_param = self.parse_args.midplane_3d_param
        if self.parse_args.midplane_rad_field:
            midplane_rad_field = self.parse_args.midplane_rad_field
        if self.parse_args.mu is not None:
            mu = self.parse_args.mu
        elif mu is None:
            mu = self.math.const['avg_gas_mass']
        if self.parse_args.nr_threads is not None:
            nr_threads = self.parse_args.nr_threads
        if self.parse_args.rot_axis_1 is not None:
            rot_axis_1 = self.parse_args.rot_axis_1
        elif rot_axis_1 is None:
            rot_axis_1 = [1, 0, 0]
        if self.parse_args.rot_axis_2 is not None:
            rot_axis_2 = self.parse_args.rot_axis_2
        elif rot_axis_2 is None:
            rot_axis_2 = [0, 1, 0]

        # Make rotation axis to unit length vectors
        if np.linalg.norm(rot_axis_1) != 1:
            rot_axis_1 /= np.linalg.norm(rot_axis_1)
        if np.linalg.norm(rot_axis_2) != 1:
            rot_axis_2 /= np.linalg.norm(rot_axis_2)

        cmd_file.write('<common>\n')
        cmd_file.write(self.dust.get_command())
        cmd_file.write('\t<axis1>\t' + str(rot_axis_1[0]) + '\t' + str(rot_axis_1[1])
                       + '\t' + str(rot_axis_1[2]) + '\n')
        cmd_file.write('\t<axis2>\t' + str(rot_axis_2[0]) + '\t' + str(rot_axis_2[1])
                       + '\t' + str(rot_axis_2[2]) + '\n')
        cmd_file.write('\n')
        cmd_file.write('\t<write_inp_midplanes>\t' + str(midplane_points) + '\n')
        cmd_file.write('\t<write_out_midplanes>\t' + str(midplane_points) + '\n')
        if midplane_3d_param is not None:
            if 0 <= len(midplane_3d_param) <= 4 and 1 <= int(midplane_3d_param[0]) <= 3:
                cmd_file.write('\t<write_3d_midplanes>')
                for i in range(len(midplane_3d_param)):
                    if 1 < i < 4:
                        cmd_file.write('\t' + str(self.math.parse(midplane_3d_param[i], 'length')))
                    else:
                        cmd_file.write('\t' + str(midplane_3d_param[i]))
                cmd_file.write('\n')
        cmd_file.write('\t<plot_inp_midplanes>\t0\n')
        cmd_file.write('\t<plot_out_midplanes>\t0\n')
        cmd_file.write('\t<midplane_zoom>\t\t' + str(midplane_zoom) + '\n')
        if midplane_rad_field:
            cmd_file.write('\t<write_radiation_field>\t1\n')
        cmd_file.write('\n')
        cmd_file.write('\t<mass_fraction>\t\t' + str(mass_fraction) + '\n')
        cmd_file.write('\t<mu>\t\t\t' + str(mu) + '\n')
        cmd_file.write('\n')
        cmd_file.write('\t<phase_function>\tPH_' + str(scattering) + '\n')
        cmd_file.write('\n')
        cmd_file.write('\t<nr_threads>\t\t' + str(nr_threads) + '\n')
        cmd_file.write('</common>\n')
        cmd_file.write('\n')

    def write_temp_part(self, cmd_file, grid_filename='grid.dat', dust_offset=False, radiation_field=False,
            temp_a_max=None, adj_tgas=None, sub_dust=True, full_dust_temp=False, add_rat=False):
        """Writes commands for temperature calculation to the command file.

        Args:
            cmd_file: Instance of the command file.
            grid_filename (str): Name of the grid file.
            dust_offset (bool): Calculate temperatures on top of
                temperatures written in the grid file.
            radiation_field (bool): Save radiation field in grid for stochastic heating
                or scattering in raytracing?
            temp_a_max (float): Calculate stochastic heating in addition to
                equilibrium temperature up to a=temp_a_max.
            adj_tgas (float): Set output gas temperature with dust
                temperature times adj_tgas.
            sub_dust (bool): Remove density in cells with temperature larger than the
                dust grain sublimation temperature.
            full_dust_temp (bool): Calculate temperature for each dust grain size individually?
            add_rat (bool): Also calculate radiative torques?
        """

        # Overwrite default values with user input
        if self.parse_args.grid_filename is not None:
            grid_filename = self.parse_args.grid_filename
        if self.parse_args.adj_tgas is not None:
            adj_tgas = self.parse_args.adj_tgas
        if self.parse_args.dust_offset:
            dust_offset = True
        if self.parse_args.full_dust_temp:
            full_dust_temp = True
        if self.parse_args.radiation_field:
            radiation_field = True
        if self.parse_args.temp_a_max is not None:
            temp_a_max = self.parse_args.temp_a_max
        if self.parse_args.conv_dens is not None:
            conv_dens = self.parse_args.conv_dens
        else:
            conv_dens = self.model.conv_parameter['conv_dens']
        if self.parse_args.conv_len is not None:
            conv_len = self.parse_args.conv_len
        else:
            conv_len = self.model.conv_parameter['conv_len']
        if self.parse_args.conv_mag is not None:
            conv_mag = self.parse_args.conv_mag
        else:
            conv_mag = self.model.conv_parameter['conv_mag']
        if self.parse_args.conv_vel is not None:
            conv_vel = self.parse_args.conv_vel
        else:
            conv_vel = self.model.conv_parameter['conv_vel']

        cmd_file.write('<task> 1\n')
        if self.source is not None and self.source.parameter['nr_photons'] > 0:
            cmd_file.write(self.source.get_command())
        cmd_file.write('\n')
        if add_rat:
            cmd_file.write('\t<cmd>\t\t\tCMD_TEMP_RAT\n')
        else:
            cmd_file.write('\t<cmd>\t\t\tCMD_TEMP\n')
        cmd_file.write('\n')
        grid_path = self.file_io.path['model'] + grid_filename
        cmd_file.write('\t<path_grid>\t\t"' + grid_path + '"\n')
        cmd_file.write('\t<path_out>\t\t"' + self.file_io.path['simulation_type'] + '"\n')
        cmd_file.write('\n')
        cmd_file.write('\t<dust_offset>\t\t' + str(int(dust_offset)) + '\n')
        cmd_file.write('\n')
        if adj_tgas is not None:
            cmd_file.write('\t<adj_tgas>\t\t' + str(adj_tgas) + '\n')
        if sub_dust:
            cmd_file.write('\t<sub_dust>\t1\n')
        if full_dust_temp:
            cmd_file.write('\t<full_dust_temp>\t1\n')
        if radiation_field:
            cmd_file.write('\t<radiation_field>\t1\n')
        if temp_a_max is not None:
            cmd_file.write('\t<stochastic_heating>\t' + str(temp_a_max) + '\n')
        cmd_file.write('\n')
        cmd_file.write('\t<conv_dens>\t\t' + str(conv_dens) + '\n')
        cmd_file.write('\t<conv_len>\t\t' + str(conv_len) + '\n')
        cmd_file.write('\t<conv_mag>\t\t' + str(conv_mag) + '\n')
        cmd_file.write('\t<conv_vel>\t\t' + str(conv_vel) + '\n')
        cmd_file.write('</task>\n')
        cmd_file.write('\n')

    def write_rat_part(self, cmd_file, grid_filename='grid.dat', conv_dens=1., 
            conv_len=1., conv_mag=1., conv_vel=1.):
        """Writes commands for radiative torque calculation to the command file.

        Args:
            cmd_file: Instance of the command file.
            grid_filename (str): Name of the grid file.
            conv_dens (float): Conversion factor for densities.
            conv_len (float): Conversion factor for spacial extents.
            conv_mag (float): Conversion factor for the magnetic field strengths.
            conv_vel (float): Conversion factor for the velocities.
        """

        # Overwrite default values with user input
        if self.parse_args.conv_dens is not None:
            conv_dens = self.parse_args.conv_dens
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_dens = self.model.conv_parameter['conv_dens']
        if self.parse_args.conv_len is not None:
            conv_len = self.parse_args.conv_len
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_len = self.model.conv_parameter['conv_len']
        if self.parse_args.conv_mag is not None:
            conv_mag = self.parse_args.conv_mag
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_mag = self.model.conv_parameter['conv_mag']
        if self.parse_args.conv_vel is not None:
            conv_vel = self.parse_args.conv_vel
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_vel = self.model.conv_parameter['conv_vel']

        cmd_file.write('<task> 1\n')
        if self.source is not None and self.source.parameter['nr_photons'] > 0:
            cmd_file.write(self.source.get_command())
        cmd_file.write('\n')
        # if self.dust_source.parameter['nr_photons'] > 0:
        #     cmd_file.write(self.dust_source.get_command())
        #     cmd_file.write('\n')
        cmd_file.write('\t<cmd>\t\t\tCMD_RAT\n')
        cmd_file.write('\n')
        if self.parse_args.grid_filename is not None:
            grid_path = self.file_io.path['model'] + self.parse_args.grid_filename
        elif os.path.isfile(self.file_io.path['simulation'] + 'temp/grid_temp.dat'):
            grid_path = self.file_io.path['simulation'] + 'temp/grid_temp.dat'
        elif os.path.isfile(self.file_io.path['simulation'] + 'temp_rat/grid_temp.dat'):
            grid_path = self.file_io.path['simulation'] + 'temp_rat/grid_temp.dat'
        else:
            grid_path = self.file_io.path['model'] + grid_filename
        cmd_file.write('\t<path_grid>\t\t"' + grid_path + '"\n')
        cmd_file.write('\t<path_out>\t\t"' + self.file_io.path['simulation_type'] + '"\n')
        cmd_file.write('\n')
        cmd_file.write('\t<conv_dens>\t\t' + str(conv_dens) + '\n')
        cmd_file.write('\t<conv_len>\t\t' + str(conv_len) + '\n')
        cmd_file.write('\t<conv_mag>\t\t' + str(conv_mag) + '\n')
        cmd_file.write('\t<conv_vel>\t\t' + str(conv_vel) + '\n')
        cmd_file.write('</task>\n')
        cmd_file.write('\n')

    def write_dust_mc_part(self, cmd_file, grid_filename='grid.dat', peel_off=True, enfsca=True,
            conv_dens=1., conv_len=1., conv_mag=1., conv_vel=1.):
        """Writes commands for Monte-Carlo radiative transfer to the command file.

        Args:
            cmd_file: Instance of the command file.
            grid_filename (str): Name of the grid file.
            peel_off (bool): Use peel_off technique?
            enfsca (bool): Enforce the first scattering of a photon packet.
            conv_dens (float): Conversion factor for densities.
            conv_len (float): Conversion factor for spacial extents.
            conv_mag (float): Conversion factor for the magnetic field strengths.
            conv_vel (float): Conversion factor for the velocities.
        """

        # Overwrite default values with user input
        if self.parse_args.peel_off is not None:
            peel_off = self.parse_args.peel_off
        if self.parse_args.enfsca is not None:
            enfsca = self.parse_args.enfsca
        if self.parse_args.conv_dens is not None:
            conv_dens = self.parse_args.conv_dens
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_dens = self.model.conv_parameter['conv_dens']
        if self.parse_args.conv_len is not None:
            conv_len = self.parse_args.conv_len
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_len = self.model.conv_parameter['conv_len']
        if self.parse_args.conv_mag is not None:
            conv_mag = self.parse_args.conv_mag
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_mag = self.model.conv_parameter['conv_mag']
        if self.parse_args.conv_vel is not None:
            conv_vel = self.parse_args.conv_vel
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_vel = self.model.conv_parameter['conv_vel']

        cmd_file.write('<task> 1\n')
        cmd_file.write(self.detector.get_dust_scattering_command())
        cmd_file.write('\n')

        if self.source is not None and self.source.parameter['nr_photons'] > 0:
            cmd_file.write(self.source.get_command())
        if self.dust_source.parameter['nr_photons'] > 0:
            cmd_file.write(self.dust_source.get_command())
            cmd_file.write('\n')
        cmd_file.write('\n')
        cmd_file.write('\t<cmd>\t\t\tCMD_DUST_SCATTERING\n')
        cmd_file.write('\n')
        if self.parse_args.grid_filename is not None:
            grid_path = self.file_io.path['model'] + self.parse_args.grid_filename
        elif os.path.isfile(self.file_io.path['simulation'] + 'temp/grid_temp.dat'):
            grid_path = self.file_io.path['simulation'] + 'temp/grid_temp.dat'
        elif os.path.isfile(self.file_io.path['simulation'] + 'temp_rat/grid_temp.dat'):
            grid_path = self.file_io.path['simulation'] + 'temp_rat/grid_temp.dat'
        else:
            grid_path = self.file_io.path['model'] + grid_filename
        cmd_file.write('\t<path_grid>\t\t"' + grid_path + '"\n')
        cmd_file.write('\t<path_out>\t\t"' + self.file_io.path['simulation_type'] + '"\n')
        cmd_file.write('\n')
        cmd_file.write('\t<peel_off>\t\t' + str(int(peel_off)) + '\n')
        cmd_file.write('\t<enfsca>\t\t' + str(int(enfsca)) + '\n')
        if not peel_off:
            cmd_file.write('\t<acceptance_angle>\t' + str(self.detector.parameter['acceptance_angle']) + '\n')
        cmd_file.write('\n')
        cmd_file.write('\t<conv_dens>\t\t' + str(conv_dens) + '\n')
        cmd_file.write('\t<conv_len>\t\t' + str(conv_len) + '\n')
        cmd_file.write('\t<conv_mag>\t\t' + str(conv_mag) + '\n')
        cmd_file.write('\t<conv_vel>\t\t' + str(conv_vel) + '\n')
        cmd_file.write('</task>\n')
        cmd_file.write('\n')

    def write_dust_part(self, cmd_file, grid_filename='grid.dat', max_subpixel_lvl=1, temp_a_max=None,
             f_c=0., f_highj=0.25, no_rt_scattering=False, conv_dens=1., conv_len=1., conv_mag=1., conv_vel=1.):
        """Writes commands for raytrace simulations to the command file.

        Args:
            cmd_file: Instance of the command file.
            grid_filename (str): Name of the grid file.
            max_subpixel_lvl (int): Maximum level of subpixel usage.
                2**(max_subpixel_lvl + 1) subpixel per pixel.
            f_c (float): Correlation factor between internal alignment and the
                other alignment mechanisms.
            f_highj (float): Amount of dust grains that rotate at high angular momentum
                and which are therefore perfectly aligned (Only for rat alignment).
            temp_a_max (float): Calculate stochastic heating in addition to
                equilibrium temperature up to a=temp_a_max.
            no_rt_scattering (bool): disable consideration of scattered light in raytracing simulations.
            conv_dens (float): Conversion factor for densities.
            conv_len (float): Conversion factor for spacial extents.
            conv_mag (float): Conversion factor for the magnetic field strengths.
            conv_vel (float): Conversion factor for the velocities.
        """

        # Overwrite default values with user input
        if self.parse_args.max_subpixel_lvl is not None:
            max_subpixel_lvl = self.parse_args.max_subpixel_lvl
        if self.parse_args.f_highj is not None:
            f_highj = self.parse_args.f_highj
        if self.parse_args.temp_a_max is not None:
            temp_a_max = self.parse_args.temp_a_max
        if self.parse_args.no_rt_scattering:
            no_rt_scattering = True
        if self.parse_args.conv_dens is not None:
            conv_dens = self.parse_args.conv_dens
        elif self.parse_args.grid_filename is not None and \
                'grid_temp.dat' not in self.parse_args.grid_filename and \
                'grid_rat.dat' not in self.parse_args.grid_filename:
            conv_dens = self.model.conv_parameter['conv_dens']
        if self.parse_args.conv_len is not None:
            conv_len = self.parse_args.conv_len
        elif self.parse_args.grid_filename is not None and \
                'grid_temp.dat' not in self.parse_args.grid_filename and \
                'grid_rat.dat' not in self.parse_args.grid_filename:
            conv_len = self.model.conv_parameter['conv_len']
        if self.parse_args.conv_mag is not None:
            conv_mag = self.parse_args.conv_mag
        elif self.parse_args.grid_filename is not None and \
                'grid_temp.dat' not in self.parse_args.grid_filename and \
                'grid_rat.dat' not in self.parse_args.grid_filename:
            conv_mag = self.model.conv_parameter['conv_mag']
        if self.parse_args.conv_vel is not None:
            conv_vel = self.parse_args.conv_vel
        elif self.parse_args.grid_filename is not None and \
                'grid_temp.dat' not in self.parse_args.grid_filename and \
                'grid_rat.dat' not in self.parse_args.grid_filename:
            conv_vel = self.model.conv_parameter['conv_vel']

        cmd_file.write('<task> 1\n')
        cmd_file.write(self.detector.get_dust_emission_command())
        cmd_file.write('\n')
        if self.source is not None and self.source.parameter['nr_photons'] > 0:
            cmd_file.write(self.source.get_command())
        cmd_file.write('\n')
        cmd_file.write(self.bg_source.get_command())
        cmd_file.write('\n')
        cmd_file.write('\t<cmd>\t\t\tCMD_DUST_EMISSION\n')
        if self.parse_args.grid_filename is not None:
            grid_path = self.file_io.path['model'] + self.parse_args.grid_filename
        elif os.path.isfile(self.file_io.path['simulation'] + 'rat/grid_rat.dat'):
            grid_path = self.file_io.path['simulation'] + 'rat/grid_rat.dat'
        elif os.path.isfile(self.file_io.path['simulation'] + 'temp/grid_temp.dat'):
            grid_path = self.file_io.path['simulation'] + 'temp/grid_temp.dat'
        elif os.path.isfile(self.file_io.path['simulation'] + 'temp_rat/grid_temp.dat'):
            grid_path = self.file_io.path['simulation'] + 'temp_rat/grid_temp.dat'
        else:
            grid_path = self.file_io.path['model'] + grid_filename
        cmd_file.write('\t<path_grid>\t\t"' + grid_path + '"\n')
        cmd_file.write('\t<path_out>\t\t"' + self.file_io.path['simulation_type'] + '"\n')
        if 'dust_' in self.parse_args.simulation_type and \
                self.parse_args.simulation_type not in ['dust_mc', 'dust_full']:
            for align in self.parse_args.simulation_type.replace('dust_','').split('_'):
                cmd_file.write('\t<align>\t\t\tALIG_' + str(align).upper() + '\n')
        cmd_file.write('\n')
        cmd_file.write('\t<max_subpixel_lvl>\t' + str(max_subpixel_lvl) + '\n')
        cmd_file.write('\t<f_highJ>\t\t' + str(f_highj) + '\n')
        cmd_file.write('\t<f_c>\t\t\t' + str(f_c) + '\n')
        cmd_file.write('\n')
        if no_rt_scattering:
            cmd_file.write('\t<rt_scattering>\t0\n')
            cmd_file.write('\n')
        if temp_a_max is not None:
            cmd_file.write('\t<stochastic_heating>\t' + str(temp_a_max) + '\n')
            cmd_file.write('\n')
        cmd_file.write('\t<conv_dens>\t\t' + str(conv_dens) + '\n')
        cmd_file.write('\t<conv_len>\t\t' + str(conv_len) + '\n')
        cmd_file.write('\t<conv_mag>\t\t' + str(conv_mag) + '\n')
        cmd_file.write('\t<conv_vel>\t\t' + str(conv_vel) + '\n')
        cmd_file.write('</task>\n')
        cmd_file.write('\n')

    def write_line_part(self, cmd_file, grid_filename='grid.dat', max_subpixel_lvl=1, turbulent_velocity=100.0,
                        do_vel_maps=True, conv_dens=1., conv_len=1., conv_mag=1., conv_vel=1.):
        """Writes commands for the radiative linetransfer to the command file.

        Args:
            cmd_file: Instance of the command file.
            grid_filename (str): Name of the grid file.
            max_subpixel_lvl (int): Maximum level of subpixel usage.
                2**(max_subpixel_lvl + 1) subpixel per pixel.
            turbulent_velocity (Quantity m/s): Turbulent velocity of the gas species [m/s].
            do_vel_maps (bool): Write out velocity channel maps?
                (large amount of pixel and channels => high memory usage!)
            conv_dens (float): Conversion factor for densities.
            conv_len (float): Conversion factor for spacial extents.
            conv_mag (float): Conversion factor for the magnetic field strengths.
            conv_vel (float): Conversion factor for the velocities.
        """

        # Overwrite default values with user input
        if self.parse_args.max_subpixel_lvl is not None:
            max_subpixel_lvl = self.parse_args.max_subpixel_lvl
        if self.parse_args.turbulent_velocity is not None:
            turbulent_velocity = self.math.parse(self.parse_args.turbulent_velocity, 'velocity')
        if self.parse_args.no_vel_maps:
            do_vel_maps = False
        if self.parse_args.conv_dens is not None:
            conv_dens = self.parse_args.conv_dens
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_dens = self.model.conv_parameter['conv_dens']
        if self.parse_args.conv_len is not None:
            conv_len = self.parse_args.conv_len
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_len = self.model.conv_parameter['conv_len']
        if self.parse_args.conv_mag is not None:
            conv_mag = self.parse_args.conv_mag
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_mag = self.model.conv_parameter['conv_mag']
        if self.parse_args.conv_vel is not None:
            conv_vel = self.parse_args.conv_vel
        elif self.parse_args.grid_filename is not None and 'grid_temp.dat' not in self.parse_args.grid_filename:
            conv_vel = self.model.conv_parameter['conv_vel']

        cmd_file.write('<task> 1\n')
        cmd_file.write(self.gas.get_command())
        cmd_file.write('\n')
        cmd_file.write(self.bg_source.get_command())
        cmd_file.write('\n')
        cmd_file.write(self.detector.get_line_command())
        cmd_file.write('\n')
        cmd_file.write('\t<cmd>\t\tCMD_LINE_EMISSION\n')
        cmd_file.write('\n')
        if self.parse_args.grid_filename is not None:
            grid_path = self.file_io.path['model'] + self.parse_args.grid_filename
        elif os.path.isfile(self.file_io.path['simulation'] + 'temp/grid_temp.dat'):
            grid_path = self.file_io.path['simulation'] + 'temp/grid_temp.dat'
        elif os.path.isfile(self.file_io.path['simulation'] + 'temp_rat/grid_temp.dat'):
            grid_path = self.file_io.path['simulation'] + 'temp_rat/grid_temp.dat'
        else:
            grid_path = self.file_io.path['model'] + grid_filename
        cmd_file.write('\t<path_grid>\t\t"' + grid_path + '"\n')
        cmd_file.write('\t<path_out>\t\t"' + self.file_io.path['simulation_type'] + '"\n')
        cmd_file.write('\n')
        if do_vel_maps:
            cmd_file.write('\t<vel_maps>\t\t1\n')
        if self.parse_args.kepler and self.source.parameter['kepler_usable']:
            cmd_file.write(
                '\t<kepler_star_mass>\t' + str(self.source.parameter['mass'] / self.math.const['M_sun']) + '\n')
        elif self.model.parameter['vel_is_speed_of_sound']:
            cmd_file.write('\t<vel_is_speed_of_sound>\t1\n')
        cmd_file.write('\t<max_subpixel_lvl>\t' + str(max_subpixel_lvl) + '\n')
        cmd_file.write(
            '\t<turbulent_velocity>\t' + str(turbulent_velocity) + '\n')
        cmd_file.write('\n')
        cmd_file.write('\t<conv_dens>\t\t' + str(conv_dens) + '\n')
        cmd_file.write('\t<conv_len>\t\t' + str(conv_len) + '\n')
        cmd_file.write('\t<conv_mag>\t\t' + str(conv_mag) + '\n')
        cmd_file.write('\t<conv_vel>\t\t' + str(conv_vel) + '\n')
        cmd_file.write('</task>\n')
        cmd_file.write('\n')

    def write_run_file(self, run_file):
        """Writes a run file to execute POLARIS directly or remotely on a server.
        Putting on a queue is also possible.

        Args:
            run_file: Instance of the run bash script file.
        """
        self.server.check_for_walltime()
        run_file.write('#!/bin/bash\n')
        run_file.write(self.server.get_command())
        # For MAC PCs with newest OS (they are not allowing the change of the LD library path in .bashrc)
        if 'darwin' in sys.platform:
            run_file.write('export LD_LIBRARY_PATH="' + self.file_io.path['polaris'] + 'CCfits/.libs/:' +
                self.file_io.path['polaris'] + 'cfitsio:${LD_LIBRARY_PATH}"\n')
            run_file.write('export DYLD_LIBRARY_PATH="' + self.file_io.path['polaris'] + 'CCfits/.libs/:' +
                self.file_io.path['polaris'] + 'cfitsio:${DYLD_LIBRARY_PATH}"\n')
        if self.parse_args.save_output:
            run_file.write(self.file_io.path['bin'] + 'polaris ' + self.file_io.path['simulation_type'] +
                'POLARIS.cmd ' + '| tee ' + self.file_io.path['simulation_type'] + '/POLARIS.out\n')
            run_file.write(r"sed -i 's/\r$//; s/\r/\r\n/g; /\r/d' " +
                self.file_io.path['simulation_type'] + '/POLARIS.out\n')
        else:
            run_file.write(self.file_io.path['bin'] + 'polaris ' + self.file_io.path['simulation_type'] +
                'POLARIS.cmd' + '\n')
        run_file.close()

    def execute_run_file(self):
        """Execute POLARIS with the command file.
        """
        os.system('chmod +x ' + self.file_io.path['simulation'] + 'run.sh')
        os.system('sh ' + self.file_io.path['simulation'] + 'run.sh')

    def queue_run_file(self):
        """Push POLARIS run file to a server queue.
        """
        os.system('chmod +x ' + self.file_io.path['simulation'] + 'run.sh')
        if self.server.parameter['queue_system'] == 'SBATCH':
            os.system('sbatch ' + self.file_io.path['simulation'] + 'run.sh')
        elif self.server.parameter['queue_system'] == 'PBS':
            os.system('qsub ' + self.file_io.path['simulation'] + 'run.sh')
        else:
            raise AttributeError('Queuing system not known!')
