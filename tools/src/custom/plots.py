#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import curve_fit

from modules.gas import GasChooser
from modules.visual import Plot


class CustomPlots:
    """This class contains modified plot routines for
    special purposes and images.
    """

    def __init__(self, basic_plot_routines, model, file_io, parse_args):
        """Initialisation of plot parameters.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
            file_io: Handles file input/output and all
                necessary paths.
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
        """
        self.basic_plots = basic_plot_routines
        self.model = model
        self.file_io = file_io
        self.parse_args = parse_args

        # Get math module
        from modules.math import Math
        self.math = Math()

        # Stokes parameter names dependent on an index
        self.stokes_parameter = {
            0: 'I',
            1: 'Q',
            2: 'U',
            3: 'V',
        }

    def plot_1(self):
        """Plot line spectrum from POLARIS simulations.
        """
        # Set data input to Jy to calculate the total flux
        self.file_io.cmap_unit = 'total'
        # Read spectrum data
        plot_data, header = self.file_io.read_spectrum(
            'line_spectrum_species_0001_line_0001')
        # Create pdf file if show_plot is not chosen
        self.file_io.init_plot_output('line_spectrum_species_0001_line_0001')
        velocity = []
        for vch in range(header['nr_channels']):
            # Get velocity of current channel
            velocity.append(1e-3 * self.math.get_velocity(vch,
                                                          header['nr_channels'], header['max_velocity']))
        for i_quantity in range(4):
            # Create Matplotlib figure
            plot = Plot(self.model, self.parse_args,
                        xlabel=r'$\mathit{v}\ [\si{\kilo\metre\per\second}]$',
                        ylabel=self.file_io.get_quantity_labels(i_quantity),
                        extent=[velocity[0], velocity[-1], None, None], with_cbar=False)
            # Plot spectrum as line
            plot.plot_line(velocity, plot_data[i_quantity, :], marker='.')
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)

    def plot_2(self):
        """Plot fits file with internal image creation.
        """
        self.model.tmp_parameter['radius_x_arcsec'] = self.math.length_conv(
            600. * self.math.const['au'], 'arcsec', self.model.parameter['distance'])
        self.model.tmp_parameter['radius_y_arcsec'] = self.math.length_conv(
            600. * self.math.const['au'], 'arcsec', self.model.parameter['distance'])
        from astropy.io import fits
        self.file_io.init_plot_output('gg_tau_miri_simulation')
        hdulist = fits.open(
            self.file_io.path['results'] + 'polaris_detector_nr0001.fits')
        header_dict = dict(
            wavelengths=[7.7e-6],
            ID=1,
            simulation_type='dust',
        )
        # Update dictionary with parameters from header
        header_dict['nr_pixel_x'] = hdulist[0].header['NAXIS1']
        header_dict['nr_pixel_y'] = hdulist[0].header['NAXIS2']
        header_dict['nr_wavelengths'] = 1
        tbldata = hdulist[0].data.T
        tbldata[np.where(tbldata < 0)] = 0
        plot = Plot(self.model, self.parse_args)
        plot.plot_imshow(tbldata, cbar_label=self.file_io.get_quantity_labels(
            0), set_bad_to_min=False)
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_5(self):
        """Plot magnetic field strength map from Zeeman splitting.
        """
        self.file_io.init_plot_output('spectrum_test')
        # Set that file_io reads data as unit/arcsec
        self.file_io.cmap_unit = 'px'
        # Read velocity channel data
        tbldata, header = self.file_io.read_vel_maps(
            'vel_channel_maps_species_0001_line_0001')
        # Calculate velocity width of a channel
        velocity_channel = np.subtract(
            np.multiply(range(header['nr_channels']), 2. *
                        header['max_velocity'] / header['nr_channels']),
            header['max_velocity'])
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, ylabel='$I\ [\mathsf{Jy}]$',
                    xlabel=r'$v\ [\si{\kilo\metre\per\second}]$')
        cut_pixel_index = int(header['nr_pixel_x'] * 1. / 2.)
        for i_x in range(0, cut_pixel_index):
            plot.plot_line(xdata=velocity_channel,
                           ydata=tbldata[0, :, i_x, cut_pixel_index], linestyle='-')

        plot.plot_line(xdata=velocity_channel,
                       ydata=tbldata[0, :, cut_pixel_index, cut_pixel_index], linestyle=':')

        for i_x in range(header['nr_pixel_x'] - 1, cut_pixel_index, -1):
            plot.plot_line(xdata=velocity_channel,
                           ydata=tbldata[0, :, i_x, cut_pixel_index], linestyle='--')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_10(self):
        """Plot integrated velocity channel map from POLARIS simulations.
        """
        self.file_io.init_plot_output(
            'comparison_subpixel', path=self.file_io.path['model'])
        # Read velocity channel data
        tbldata_i_list = []
        vmin = 1e99
        vmax = 1e-99
        for i_sb in range(3):
            self.file_io.set_path_from_str('plot', model_name='disk', simulation_type='line',
                                           simulation_name='compare_mol3d_polaris_theta_0_sp_' + str(i_sb))
            tbldata, header = self.file_io.read_int_vel_map(
                'int_channel_map_species_0001_line_0001')
            tbldata_i = tbldata[0, :, :]
            tbldata_i_list.append(tbldata_i)
            vmin = min(vmin, np.min(tbldata_i[np.nonzero(tbldata_i)]))
            vmax = max(vmax, np.max(tbldata_i[np.nonzero(tbldata_i)]))
        for i_data in range(len(tbldata_i_list)):
            # Create Matplotlib figure
            plot = Plot(self.model, self.parse_args)
            plot.plot_title(
                r'$\mathsf{Max\ level\ of\ subpixel}=' + str(i_data) + '$')
            # Plot chosen quantity to velocity map
            plot.plot_imshow(tbldata_i_list[i_data], cbar_label=self.file_io.get_quantity_labels(0, int_map=True),
                             vmin=vmin, vmax=vmax, extend='neither')
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)

    def plot_19(self):
        """Plot Raytrace and Monte-Carlo results as map from POLARIS simulations.
        """
        # Define title and amount of components per map
        simu_types = ['dust_mc', 'dust']
        n_wl = 23
        for i_wl in range(n_wl):
            # Create pdf file if show_plot is not chosen
            self.file_io.init_plot_output('polaris_total_maps_nr' + str(i_wl + 1).zfill(4),
                                          path=self.file_io.path['simulation'])
            total_map_data = np.ones(1)
            for i_simulation_type in range(len(simu_types)):
                # Set paths of each simulation
                self.file_io.set_path_from_str('plot', self.parse_args.model_name, self.parse_args.simulation_name,
                                               simu_types[i_simulation_type])
                # Read raytrace results from file
                map_data, header = self.file_io.read_emission_map(
                    'polaris_detector_nr' + str(i_wl + 1).zfill(4))
                if total_map_data.all():
                    total_map_data = map_data
                else:
                    for i in range(len(total_map_data[:, 0, 0])):
                        for i_x in range(len(total_map_data[i, :, 0])):
                            for i_y in range(len(total_map_data[i, i_x, :])):
                                total_map_data[i, i_x,
                                               i_y] += map_data[i, i_x, i_y]
            map_title = r'$\mathsf{total\ emission},\ \lambda=\SI{' + \
                        str(self.math.latex_float(
                            header['wavelength'] * 1e6)) + '}{\micro\metre}$'
            # Create one plot per component of the simulation
            for i_quantity in range(4):
                # Take colorbar label from quantity id
                cbar_label = self.file_io.get_quantity_labels(i_quantity)
                # Take data for current quantity
                tbldata = total_map_data[i_quantity, :, :]
                # skip zero content plots
                if not tbldata.any():
                    continue
                # Create Matplotlib figure
                plot = Plot(self.model, self.parse_args,
                            title=map_title)
                if i_quantity is 0:
                    # Intensity plot
                    plot.plot_imshow(tbldata, cbar_label=cbar_label)
                elif i_quantity in [1, 2, 3]:
                    # Plot polarization fluxes symmetrically around zero
                    if np.nanmax(np.abs(tbldata)) > 0.:
                        maximum_flux = np.asscalar(np.nanmax(np.abs(tbldata)))
                    else:
                        maximum_flux = 1.0
                    plot.plot_imshow(tbldata, cmap='coolwarm', vmin=-maximum_flux, vmax=maximum_flux,
                                     extend='neither', cbar_label=cbar_label)
                # Save figure to pdf file or print it on screen
                plot.save_figure(self.file_io)
            # Close pdf file
            self.file_io.close_plot_output()

    def plot_20(self):
        """Plot total SED from POLARIS simulations.
        """
        ''' Initialisation '''
        # Set data input to Jy to calculate the total flux
        if self.parse_args.cmap_unit is None:
            self.file_io.cmap_unit = 'total'
        # Create pdf file if show_plot is not chosen
        self.file_io.init_plot_output(
            'polaris_total_sed', path=self.file_io.path['simulation'])
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, xlabel=r'$\lambda\ [\si{\metre}]$',
                    ylabel=r'$\mathit{I}\ [\si{Jy}]$', with_cbar=False)

        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', self.parse_args.model_name, self.parse_args.simulation_name, 'dust')
        # Read raytrace results from file
        ray_data, ray_header = self.file_io.read_emission_sed(
            'polaris_detector_nr0001_sed')
        # Plot spectral energy distribution
        plot.plot_line(ray_header['wavelengths'], ray_data[0, :],
                       log='y', label=r'$\mathsf{thermal\ emission}$')

        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', self.parse_args.model_name, self.parse_args.simulation_name, 'dust_mc')
        # Read raytrace results from file
        mc_data, mc_header = self.file_io.read_emission_sed(
            'polaris_detector_nr0001_sed')
        # Plot spectral energy distribution (direct)
        plot.plot_line(mc_header['wavelengths'], mc_data[6, :],
                       log='y', label=r'$\mathsf{direct\ starlight}$')
        # Plot spectral energy distribution (scattered)
        plot.plot_line(mc_header['wavelengths'], mc_data[7, :],
                       log='y', label=r'$\mathsf{scattered\ starlight}$')

        wavelengths_total = []
        quantity_total = []
        offset = 0
        for i_wl in range(mc_header['nr_wavelengths'] + ray_header['nr_wavelengths']):
            if i_wl - offset >= ray_header['nr_wavelengths']:
                break
            elif i_wl >= mc_header['nr_wavelengths']:
                wavelengths_total.append(
                    ray_header['wavelengths'][i_wl - offset])
                quantity_total.append(ray_data[0, i_wl - offset])
            elif mc_header['wavelengths'][i_wl] == ray_header['wavelengths'][i_wl - offset]:
                wavelengths_total.append(ray_header['wavelengths'][i_wl])
                quantity_total.append(
                    mc_data[0, i_wl] + ray_data[0, i_wl - offset])
            elif mc_header['wavelengths'][i_wl] < ray_header['wavelengths'][i_wl - offset]:
                offset += 1
                wavelengths_total.append(mc_header['wavelengths'][i_wl])
                quantity_total.append(mc_data[0, i_wl])
        # Plot total spectral energy distribution
        plot.plot_line(wavelengths_total, quantity_total, log='xy', linestyle='--',
                       label=r'$\mathsf{total\ SED}$')

        ''' Saving and Legend '''
        # Plot the legend
        plot.plot_legend()
        # adapt limits
        plot.set_limits(
            limits=[None, None, 1e-5, np.max(quantity_total) * 1e1])
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_21(self):
        """Plot magnetic field strength map from Zeeman splitting as an animation/video.
        """
        import matplotlib.animation as animation
        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', 'mhd_run', 'animation_test', 'dust')
        # Read raytrace results from file
        raytrace_data, header = self.file_io.read_emission_map(
            'polaris_detector_nr0001')
        # Set vector size to match with 64 x 64 pixel sized image
        self.file_io.vec_field_size = 8
        # Number of angles
        n_ph = 181
        # Degree of polarization
        i_quantity = 7
        # Take colorbar label from quantity id
        cbar_label = self.file_io.get_quantity_labels(i_quantity)
        # list of tbldata
        tbldata_list = np.zeros(
            (n_ph, header['nr_pixel_x'], header['nr_pixel_y']))
        artist_list = []
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, zoom_factor=1.0,
                    ax_unit='pc', image_type='animation')
        for i_phi in range(n_ph):
            # Print progress
            print('Create animation image number:', i_phi + 1, '/', n_ph)
            # Read raytrace results from file
            raytrace_data, header = self.file_io.read_emission_map(
                'polaris_detector_nr' + str(i_phi + 1).zfill(4))
            # Take data for current quantity
            tbldata_list[i_phi, :, :] = raytrace_data[i_quantity, :, :]
            vmin = 0.
            vmax = 45.
            if i_phi < n_ph - 1:
                # Plot quantity to velocity map
                plot.plot_imshow(tbldata=tbldata_list[i_phi, :, :],
                                 vmin=vmin, vmax=vmax, extend='neither')
            else:
                plot.plot_imshow(tbldata=tbldata_list[i_phi, :, :], cbar_label=cbar_label,
                                 vmin=vmin, vmax=vmax, extend='neither')
            # For PI and P, plot polarization vectors
            vec_field_data = self.file_io.read_polarization_vectors(
                raytrace_data)
            vector_color = self.math.get_vector_color(
                tbldata_list[i_phi, :, :], log=self.parse_args.log)
            plot.plot_quiver(vec_field_data, color=vector_color,
                             vmin=vmin, vmax=vmax)
            # Make Layout tight
            plot.make_tight_layout()
            artist_list.append(plot.animation_images)
            plot.animation_images = []

        print(np.asscalar(np.nanmin(tbldata_list)),
              np.asscalar(np.nanmax(tbldata_list)))

        # Create Matplotlib figure
        animation = animation.ArtistAnimation(
            plot.fig, artist_list, interval=100, blit=True)
        # Save image to movie file
        animation.save(
            self.file_io.path['plots'] + 'continuum_polarization_animation.mp4', dpi=800)

    def plot_22(self):
        """Plot 2D graph of optical depth profile.
        """
        self.file_io.init_plot_output('tau_profile')
        data = self.file_io.read_data_file('tau_profile.dat')
        # freq = self.math.const['c'] / 7.7e-6
        N_wl = 31
        N = len(data[:, 0])
        N_pos = int(N / N_wl)
        data = data.reshape(N_pos, N_wl, 5).transpose(1, 0, 2)
        wl_list = ['1.5364e-06', '3.61733e-06', '4.98706e-06', '7.6522e-06',
                   '1.05498e-05', '1.45445e-05']

        # initial_intensity = radius.copy()
        # final_intensity = radius.copy()
        # for i in range(N):
        # radius[i] = np.sqrt(data[i, 0] ** 2 + data[i, 1] ** 2 + data[i, 2] ** 2) / self.math.const['au']
        # radius_proj[i] = np.sqrt(data[i, 0] ** 2 + (data[i, 1] * np.cos(37. / 180. * np.pi)) ** 2 + data[i, 2] ** 2) / self.math.const['au']
        # optical_depth[i] = data[i, 3]
        # mult =  self.math.const['c'] / (freq * freq)
        # mult *= 1e+26  / (140. * self.math.const['pc']) ** 2
        # initial_intensity[i] = data[i, 4] * mult * 1e3
        # final_intensity[i] = data[i, 5] * mult * 1e3

        radius = np.zeros(N_pos)
        radius_proj = np.zeros(N_pos)
        optical_depth = np.zeros((N_wl, N_pos))
        wavelength = np.zeros(N_wl)
        for i_pos in range(N_pos):
            radius[i_pos] = np.sqrt(data[0, i_pos, 1] ** 2 + data[0, i_pos, 2]
                                    ** 2 + data[0, i_pos, 3] ** 2) / self.math.const['au']
            radius_proj[i_pos] = np.sqrt(data[0, i_pos, 1] ** 2 + (data[0, i_pos, 2] * np.cos(
                37. / 180. * np.pi)) ** 2 + data[0, i_pos, 3] ** 2) / self.math.const['au']
            for i_wl in range(N_wl):
                wavelength[i_wl] = data[i_wl, i_pos, 0]
                optical_depth[i_wl, i_pos] = data[i_wl, i_pos, 4]

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, xlabel=r'$R_\textsf{midplane}\ [\si{au}]$',
                    ylabel=r'$\tau_\textsf{obs}$', with_cbar=False)
        for i_wl in range(N_wl):
            if str(wavelength[i_wl]) in wl_list:
                plot.plot_line(radius, optical_depth[i_wl, :], label=r'$\lambda=\SI{' +
                               str(round(wavelength[i_wl] * 1e6, 2)) + '}{\micro\metre}$')
                plot.plot_legend()
        plot.save_figure(self.file_io)

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, xlabel=r'$R_\textsf{plane-of-sky}\ [\si{au}]$',
                    ylabel=r'$\tau_\textsf{obs}$', with_cbar=False)
        for i_wl in range(N_wl):
            if str(wavelength[i_wl]) in wl_list:
                plot.plot_line(radius_proj, optical_depth[i_wl, :], label=r'$\lambda=\SI{' +
                               str(round(wavelength[i_wl] * 1e6, 2)) + '}{\micro\metre}$')
                plot.plot_legend()
        plot.save_figure(self.file_io)

    def plot_23(self):
        """Plot multiple spectral line profiles into one image.
        """
        self.file_io.init_plot_output('multiple_spectral_profiles')
        # Chosen quantity to plot (0: I, 1: Q, ...)
        i_quantity = 0
        # Calculation of the column density of each simulation
        col_dens = [5e-8, 1e-7, 1e-6, 1e-5]
        # Density at abundance = 1e-5
        col_dens = np.multiply(col_dens, 2.21859e+14)
        col_dens = np.multiply(
            col_dens, self.math.const['au'])  # To column density
        # Read spectrum data
        tbldata = self.file_io.read_data_file(
            'line_spectrum_species_0001_line_0001.det', skip_header=21)
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args,
                    ylabel='$F_\mathsf{' +
                    self.stokes_parameter[i_quantity] + '}\ [\mathsf{Jy}]$',
                    xlabel=r'$v\ [\si{\kilo\metre\per\second}]$', with_cbar=False,
                    extent=[tbldata[0, 0] * 1e-3, tbldata[-1, 0] * 1e-3, 0., 1.5e-4])
        for line_number in range(len(col_dens)):
            # Read spectrum data
            tbldata = self.file_io.read_data_file(
                'line_spectrum_species_0001_line_' +
                str(line_number + 1).zfill(4) + '.det',
                skip_header=21)
            # Plot spectrum as line
            plot.plot_line(tbldata[:, 0] * 1e-3, tbldata[:, i_quantity + 1], marker='.',
                           label='column density = ' + '{:.2e}'.format(col_dens[line_number], 2))
        # Plot the legend
        plot.plot_legend()
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_27(self):
        """Plot magnetic field strength weighted with the intensity from Zeeman splitting.
        """
        # Load necessary modules
        from sys import stdout
        self.file_io.init_plot_output('magnetic_field_times_intensity')
        # Set that file_io reads data as unit/arcsec
        self.file_io.cmap_unit = 'px'
        # Read velocity channel data
        tbldata_vel_map, header = self.file_io.read_vel_maps(
            'vel_channel_maps_species_0001_line_0001')
        # Read velocity channel data
        tbldata_int_map = self.file_io.read_int_vel_map(
            'int_channel_map_species_0001_line_0001')[0]
        # Get gas module
        gas_chooser = GasChooser(self.file_io, self.parse_args)
        gas = gas_chooser.get_module_from_name(header['species_name'])
        # Calculate velocity width of a channel
        vel_channel_width = 2. * \
            header['max_velocity'] / (header['nr_channels'] - 1)
        # Initialise output arrays for the derived and reference magnetic field strength
        mag_times_intensity = np.zeros(
            (header['nr_pixel_x'], header['nr_pixel_y']))
        for i_x in range(header['nr_pixel_x']):
            for i_y in range(header['nr_pixel_y']):
                stdout.write('--- Calculate magnetic field strength for each pixel: ' +
                             str(int(100. * (i_x * header['nr_pixel_y'] + i_y + 1) /
                                     (header['nr_pixel_x'] * header['nr_pixel_y']))) + ' % \r')
                stdout.flush()
                # Calculate velocity shift by comparing I and V profiles
                velocity_shift = self.math.velocity_shift_from_profile(
                    tbldata_vel_map[0, :, i_x, i_y], tbldata_vel_map[3, :, i_x, i_y], vel_channel_width)

                # Get magnetic field strength from Zeeman information of used species
                mag_times_intensity[i_x, i_y] = abs(
                    gas.shift_2_mag(velocity_shift, header['frequency'], header['i_transition']) *
                    tbldata_int_map[0, i_x, i_y])

        # Normalize the values to 1.0 as arbitrary units
        mag_times_intensity = np.divide(mag_times_intensity,
                                        np.nansum(tbldata_int_map[0, :, :]))
        print(np.nansum(mag_times_intensity))

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args)
        # Plot quantity to velocity map
        plot.plot_imshow(mag_times_intensity, cbar_label=r'$\Delta B_\mathsf{LOS}\ [\mathsf{arb.\ units}]$',
                         set_bad_to_min=True)
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_28(self):
        """Plot magnetic field strength map from Zeeman splitting.
        """
        # Load necessary modules
        from sys import stdout
        self.file_io.init_plot_output('magnetic_field_map')
        # Set that file_io reads data as unit/arcsec
        self.file_io.cmap_unit = 'px'
        # Read velocity channel data
        tbldata, header = self.file_io.read_vel_maps(
            'vel_channel_maps_species_0001_line_0001')
        # Get gas module
        gas_chooser = GasChooser(self.file_io, self.parse_args)
        gas = gas_chooser.get_module_from_name(header['species_name'])
        # Calculate velocity width of a channel
        vel_channel_width = 2. * \
            header['max_velocity'] / (header['nr_channels'] - 1)
        #: float: define global zoom factor
        zoom_factor = 1
        # Initialise output arrays for the derived and reference magnetic field strength
        derived_magnetic_field = np.zeros(
            (header['nr_pixel_x'], header['nr_pixel_y']))
        reference_magnetic_field = np.zeros(
            (header['nr_pixel_x'], header['nr_pixel_y']))
        magnetic_field_difference = np.zeros(
            (header['nr_pixel_x'], header['nr_pixel_y']))
        magnetic_field_rel_difference = np.zeros(
            (header['nr_pixel_x'], header['nr_pixel_y']))
        abs_magnetic_field_rel_difference = np.zeros(
            (header['nr_pixel_x'], header['nr_pixel_y']))
        column_density = np.zeros((header['nr_pixel_x'], header['nr_pixel_y']))

        for i_x in range(header['nr_pixel_x']):
            for i_y in range(header['nr_pixel_y']):
                stdout.write('--- Calculate magnetic field strength for each pixel: ' +
                             str(int(100. * (i_x * header['nr_pixel_y'] + i_y + 1) /
                                     (header['nr_pixel_x'] * header['nr_pixel_y']))) + ' % \r')
                stdout.flush()
                # Calculate velocity shift by comparing I and V profiles
                velocity_shift = self.math.velocity_shift_from_profile(
                    tbldata[0, :, i_x, i_y], tbldata[3, :, i_x, i_y], vel_channel_width)
                # Get magnetic field strength from Zeeman information of used species
                derived_magnetic_field[i_x, i_y] = gas.shift_2_mag(velocity_shift, header['frequency'],
                                                                   header['i_transition']) * 1e10
                # Channel number 0 contains the LOS magnetic field strength of the model
                reference_magnetic_field[i_x,
                                         i_y] = tbldata[4, 0, i_x, i_y] * 1e10
                magnetic_field_difference[i_x, i_y] = \
                    (derived_magnetic_field[i_x, i_y] -
                     reference_magnetic_field[i_x, i_y])
                magnetic_field_rel_difference[i_x, i_y] = 100.0 * (
                    (derived_magnetic_field[i_x, i_y] - reference_magnetic_field[i_x, i_y]) /
                    reference_magnetic_field[i_x, i_y])
                abs_magnetic_field_rel_difference[i_x, i_y] = 100.0 * abs(
                    derived_magnetic_field[i_x, i_y] - reference_magnetic_field[i_x, i_y]) / abs(
                    reference_magnetic_field[i_x, i_y])
                column_density[i_x, i_y] = tbldata[4, 2, i_x, i_y] * 1e10

                # Visualize the line profiles with additional information as line plots
                cut_pixel_index = int(header['nr_pixel_y'] * 2. / 4.)
                if i_y == cut_pixel_index and self.parse_args.show_plot \
                        and abs(i_x - int(header['nr_pixel_x'] * 2. / 4.)) < 20:
                    print('pixel nr: i_x =', i_y, '/ i_y =', i_x, )
                    # Create Matplotlib figure
                    plot = Plot(self.model, self.parse_args, ylabel='$I\ [\mathsf{Jy}]$', with_cbar=False,
                                xlabel=r'$v\ [\si{\kilo\metre\per\second}]$', limits=[None, None, None, None])
                    '''
                        plot.plot_line(xdata=range(header['nr_channels']),
                                       ydata=tbldata[1, :, i_x, i_y] / max(
                                           np.max(
                                               abs(tbldata[1, :, i_x, i_y])),
                                           np.max(abs(tbldata[2, :, i_x, i_y]))),
                                       marker='.', color='blue', linestyle='-',
                                       label=r'$F_\mathsf{Q}$')
                        plot.plot_line(xdata=range(header['nr_channels']),
                                       ydata=tbldata[2, :, i_x, i_y] / max(
                                           np.max(
                                               abs(tbldata[1, :, i_x, i_y])),
                                           np.max(abs(tbldata[2, :, i_x, i_y]))),
                                       marker='.', color='green', linestyle='-',
                                       label=r'$F_\mathsf{U}$')
                        '''
                    velocity_channel = np.subtract(
                        np.multiply(range(
                            header['nr_channels']), 2. * header['max_velocity'] / header['nr_channels']),
                        header['max_velocity'])
                    intensity_derivative = np.gradient(
                        tbldata[0, :, i_x, i_y], vel_channel_width)
                    plot.plot_line(xdata=velocity_channel, ydata=np.multiply(intensity_derivative[:], velocity_shift),
                                   marker='.', color='blue', label=r'$\mathsf{d}I \cdot\nu_\mathsf{measured}$')
                    # plot.plot_line(xdata=velocity_channel, ydata=intensity_derivative[:, i_x, i_y] * gas.mag_2_shift(
                    #        tbldata[4, 0, i_x, i_y], header['frequency'], header['i_transition']), marker='.',
                    #        color='green', linestyle='--', label=r'$\mathsf{d}I \cdot' + r'\nu_\mathsf{ideal}$')
                    plot.plot_line(xdata=velocity_channel, ydata=tbldata[3, :, i_x, i_y], marker='.', color='red',
                                   linestyle='-.', label=r'$V$')
                    # plot.plot_line(xdata=velocity_channel,ydata=tbldata[0, :, i_x, i_y],marker='.', color='black',
                    #                linestyle='-', label=r'$F_\mathsf{I}$')
                    # Plot the legend
                    plot.plot_legend()
                    print('derived:',
                          gas.shift_2_mag(
                              velocity_shift, header['frequency'], header['i_transition']) * 1e10,
                          'reference:', tbldata[4, 0, i_x, i_y] * 1e10)
                    plot.save_figure(self.file_io)

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args,
                    zoom_factor=zoom_factor)
        # Plot quantity to velocity map
        plot.plot_imshow(derived_magnetic_field, cbar_label=r'$B_\mathsf{LOS}\ [\si{\micro G}]$', cmap='coolwarm',
                         vmin=-
                         np.asscalar(
                             np.nanmax(np.abs(derived_magnetic_field))),
                         vmax=np.asscalar(np.nanmax(np.abs(derived_magnetic_field))), extend='neither')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args,
                    zoom_factor=zoom_factor)
        # Plot quantity to velocity map
        plot.plot_imshow(reference_magnetic_field, cbar_label=r'$B_\mathsf{LOS}\ [\si{\micro G}]$', cmap='coolwarm',
                         vmin=-
                         np.asscalar(
                             np.nanmax(np.abs(reference_magnetic_field))),
                         vmax=np.asscalar(np.nanmax(np.abs(reference_magnetic_field))), extend='neither')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args,
                    zoom_factor=zoom_factor)
        # Plot quantity to velocity map
        plot.plot_imshow(magnetic_field_difference, cbar_label=r'$\Delta B_\mathsf{LOS}\ [\si{\micro G}]$',
                         cmap='coolwarm', vmin=-np.asscalar(np.nanmax(np.abs(magnetic_field_difference))),
                         vmax=np.asscalar(np.nanmax(np.abs(magnetic_field_difference))), extend='neither')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args,
                    zoom_factor=zoom_factor)
        # Plot quantity to velocity map
        plot.plot_imshow(magnetic_field_rel_difference,
                         cbar_label=r'$\Delta B_\mathsf{LOS}/B_\mathsf{LOS}\ [\si{\percent}]$',
                         cmap='coolwarm', extend='neither',
                         vmin=-
                         np.asscalar(
                             np.nanmax(np.abs(magnetic_field_rel_difference))),
                         vmax=np.asscalar(np.nanmax(np.abs(magnetic_field_rel_difference))))
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args,
                    zoom_factor=zoom_factor)
        # Plot quantity to velocity map (optional extend='max', vmin=0., vmax=10.)
        plot.plot_imshow(abs_magnetic_field_rel_difference,
                         cbar_label=r'$|\Delta B_\mathsf{LOS}|/B_\mathsf{LOS}\ [\si{\percent}]$')
        # plot.axarr[0].set_xticks([-5000, -2500, 0, 2500, 5000])
        # plot.axarr[0].set_yticks([-5000, -2500, 0, 2500, 5000])
        plot.plot_contour(magnetic_field_rel_difference, colors='white',
                          levels=[0.1, 1.0, 3.0, 10.0], label_type='percentage_0_decimal')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_29(self):
        """Test spectrum
        """
        # Set that file_io reads data as unit/arcsec
        self.file_io.cmap_unit = 'px'
        # Read velocity channel data
        tbldata, header = self.file_io.read_vel_maps(
            'vel_channel_maps_species_0001_line_0001')
        for i_x in range(header['nr_pixel_x']):
            for i_y in range(header['nr_pixel_y']):
                # Visualize the line profiles with additional information as line plots
                x = 101
                y = 131
                if i_y in [x, y] and i_x in [x, y] and self.parse_args.show_plot:
                    print('pixel nr: i_x =', i_y, '/ i_y =', i_x, )
                    # Create Matplotlib figure
                    plot = Plot(self.model, self.parse_args, ylabel='$I\ [\mathsf{Jy}]$', with_cbar=False,
                                xlabel=r'$v\ [\si{\kilo\metre\per\second}]$', limits=[None, None, None, None])
                    velocity_channel = np.subtract(
                        np.multiply(range(
                            header['nr_channels']), 2. * header['max_velocity'] / header['nr_channels']),
                        header['max_velocity'])
                    plot.plot_line(xdata=velocity_channel, ydata=tbldata[0, :, i_x, i_y], marker='.', color='red',
                                   linestyle='-.', label=r'$I$')
                    plot.save_figure(self.file_io)

    def plot_37(self):
        """Plot schematic illustration of two Zeeman split spectral lines
        with either low or high magnetic field strengths.
        """
        self.file_io.init_plot_output(
            'derive_mag_field_spectrum_3', path=self.file_io.path['model'])
        for i_mag_field in range(2):
            # Create Matplotlib figure
            plot = Plot(self.model, self.parse_args, xlabel=r'$\mathsf{Frequency}$',
                        ylabel=r'$\mathsf{Flux}$', extent=[-3., 3., None, None], with_cbar=False)
            # Import norm function
            from scipy.stats import norm
            # Create xaxis
            tmp_xdata = np.linspace(-3, 3, 1000)

            if i_mag_field == 0:
                # Create gauss curves
                tmp_ydata_pi = norm(loc=0, scale=0.6)
                tmp_ydata_sp = norm(loc=0, scale=0.6)
                tmp_ydata_sm = norm(loc=0, scale=0.6)
                tmp_ydata = np.add(tmp_ydata_pi.pdf(tmp_xdata),
                                   np.add(0.5 * tmp_ydata_sp.pdf(tmp_xdata),
                                          0.5 * tmp_ydata_sm.pdf(tmp_xdata)))
                plot.plot_line(tmp_xdata, tmp_ydata, color='black',
                               label=r'$I$', no_ticks=True)
            elif i_mag_field == 1:
                # Create gauss curves
                tmp_ydata_pi = norm(loc=0, scale=0.6)
                tmp_ydata_sp = norm(loc=-1, scale=0.6)
                tmp_ydata_sm = norm(loc=1, scale=0.6)
                tmp_ydata = np.add(tmp_ydata_pi.pdf(tmp_xdata),
                                   np.add(0.5 * tmp_ydata_sp.pdf(tmp_xdata),
                                          0.5 * tmp_ydata_sm.pdf(tmp_xdata)))
                plot.plot_line(tmp_xdata, tmp_ydata, color='gray',
                               label=r'$I$', no_ticks=True)
                plot.plot_line(tmp_xdata, 0.5 * tmp_ydata_sp.pdf(tmp_xdata), color=plot.colorpalette[2],
                               label=r'$\sigma_+$',
                               no_ticks=True)
                plot.plot_line(tmp_xdata, tmp_ydata_pi.pdf(
                    tmp_xdata), color='#ff8000', label=r'$\pi$', no_ticks=True)
                plot.plot_line(tmp_xdata, 0.5 * tmp_ydata_sm.pdf(tmp_xdata), color=plot.colorpalette[0],
                               label=r'$\sigma_-$',
                               no_ticks=True)
                plot.plot_double_arrow(
                    arrow_origin=[-1.05, 0.38], arrow_offset=[1.1, 0.])
                plot.plot_text(
                    text_pos=[-0.5, 0.42], text=r'$\Delta\nu_\mathsf{z}$', color='black')

            # Plot the legend
            plot.plot_legend()
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)

    def plot_38(self):
        """Plot schematic illustration of two Zeeman split spectral lines
        with either low or high magnetic field strengths.
        """
        self.file_io.init_plot_output(
            'derive_mag_field_spectrum_3_german', path=self.file_io.path['model'])
        for i_mag_field in range(2):
            # Create Matplotlib figure
            plot = Plot(self.model, self.parse_args, xlabel=r'$\mathsf{Frequenz}$', with_cbar=False,
                        ylabel=r'$\mathsf{Strahlungsfluss}$', extent=[-3., 3., None, None], language='german')
            # Import norm function
            from scipy.stats import norm
            # Create xaxis
            tmp_xdata = np.linspace(-3, 3, 1000)

            if i_mag_field == 0:
                # Create gauss curves
                tmp_ydata_pi = norm(loc=0, scale=0.6)
                tmp_ydata_sp = norm(loc=0, scale=0.6)
                tmp_ydata_sm = norm(loc=0, scale=0.6)
                tmp_ydata = np.add(tmp_ydata_pi.pdf(tmp_xdata),
                                   np.add(0.5 * tmp_ydata_sp.pdf(tmp_xdata),
                                          0.5 * tmp_ydata_sm.pdf(tmp_xdata)))
                plot.plot_line(tmp_xdata, tmp_ydata, color='black',
                               label=r'$I$', no_ticks=True)
            elif i_mag_field == 1:
                # Create gauss curves
                tmp_ydata_pi = norm(loc=0, scale=0.6)
                tmp_ydata_sp = norm(loc=-1, scale=0.6)
                tmp_ydata_sm = norm(loc=1, scale=0.6)
                tmp_ydata = np.add(tmp_ydata_pi.pdf(tmp_xdata),
                                   np.add(0.5 * tmp_ydata_sp.pdf(tmp_xdata),
                                          0.5 * tmp_ydata_sm.pdf(tmp_xdata)))
                plot.plot_line(tmp_xdata, tmp_ydata, color='gray',
                               label=r'$I$', no_ticks=True)
                plot.plot_line(tmp_xdata, 0.5 * tmp_ydata_sp.pdf(tmp_xdata), color=plot.colorpalette[2],
                               label=r'$\sigma_+$',
                               no_ticks=True)
                plot.plot_line(tmp_xdata, tmp_ydata_pi.pdf(
                    tmp_xdata), color='#ff8000', label=r'$\pi$', no_ticks=True)
                plot.plot_line(tmp_xdata, 0.5 * tmp_ydata_sm.pdf(tmp_xdata), color=plot.colorpalette[0],
                               label=r'$\sigma_-$',
                               no_ticks=True)
                plot.plot_double_arrow(
                    arrow_origin=[-1.05, 0.38], arrow_offset=[1.1, 0.])
                plot.plot_text(text_pos=[-0.5, 0.42],
                               text=r'$\Delta\nu_\mathsf{z}$')

            # Plot the legend
            plot.plot_legend()
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)

    def plot_39(self):
        """Plot LOS magnetic field strength dependent on the
        column density.
        """
        # Load necessary modules
        from sys import stdout
        self.file_io.init_plot_output('mag_field_column_density')
        # Read velocity channel data
        tbldata, header = self.file_io.read_vel_maps(
            'vel_channel_maps_species_0001_line_0001')
        # Get gas module
        gas_chooser = GasChooser(self.file_io, self.parse_args)
        gas = gas_chooser.get_module_from_name(header['species_name'])
        # Calculate velocity width of a channel
        vel_channel_width = 2. * \
            header['max_velocity'] / (header['nr_channels'] - 1)
        # Initialise output arrays for the derived and reference magnetic field strength
        derived_magnetic_field = np.zeros(
            (header['nr_pixel_x'], header['nr_pixel_y']))
        reference_magnetic_field = np.zeros(
            (header['nr_pixel_x'], header['nr_pixel_y']))
        # Initialise data for column density and LOS magnetic field strength
        column_density_xdata = np.zeros(
            header['nr_pixel_x'] * header['nr_pixel_y'])
        derived_magnetic_field_ydata = np.zeros(
            header['nr_pixel_x'] * header['nr_pixel_y'])
        reference_magnetic_field_ydata = np.zeros(
            header['nr_pixel_x'] * header['nr_pixel_y'])
        for i_x in range(header['nr_pixel_x']):
            for i_y in range(header['nr_pixel_y']):
                stdout.write('--- Calculate magnetic field strength for each pixel: ' +
                             str(int(100. * (i_x * header['nr_pixel_y'] + i_y + 1) /
                                     (header['nr_pixel_x'] * header['nr_pixel_y']))) + ' % \r')
                stdout.flush()
                # Calculate velocity shift by comparing I and V profiles
                velocity_shift = self.math.velocity_shift_from_profile(
                    tbldata[0, :, i_x, i_y], tbldata[3, :, i_x, i_y], vel_channel_width)
                # Get magnetic field strength from Zeeman information of used species
                derived_magnetic_field[i_x, i_y] = gas.shift_2_mag(velocity_shift, header['frequency'],
                                                                   header['i_transition']) * 1e10
                # Channel number 0 contains the LOS magnetic field strength of the model
                reference_magnetic_field[i_x,
                                         i_y] = tbldata[4, 0, i_x, i_y] * 1e10
                # Channel number 2 contains the column density (convert to cm^-2)
                column_density = tbldata[4, 2, i_x, i_y] * 1e-4
                # mu = 2.8 to take He into account (Crutcher 2004)
                column_density_xdata[i_x * header['nr_pixel_y'] +
                                     i_y] = column_density * 2.0 / 2.8
                derived_magnetic_field_ydata[i_x * header['nr_pixel_y'] +
                                             i_y] = derived_magnetic_field[i_x, i_y]
                reference_magnetic_field_ydata[i_x * header['nr_pixel_y'] +
                                               i_y] = reference_magnetic_field[i_x, i_y]

        def plot_mass_to_flux_line(_plot):
            tmp_xdata = np.logspace(18, 26)
            m_to_f = 1 / (2 * np.pi * np.sqrt(self.math.const['G'] * 1e3))
            gas_mass = 2.8 * self.math.const['amu'] * 1e3
            factor = gas_mass / (m_to_f * 2.0 * 1e-6)
            tmp_ydata = []
            for i in range(len(tmp_xdata)):
                tmp_ydata.append(factor * tmp_xdata[i])
            _plot.plot_line(xdata=tmp_xdata, ydata=tmp_ydata,
                            label=r'$\left(\frac{B_\mathsf{LOS}}{\si{\micro G}}\right) = '
                                  + str(round(factor * 1e21, 1)) +
                            '\times 10^{-21} ' + r'\left(\frac{N_\mathsf{'
                            r'H_2}}{\mathsf{cm^{'
                            r'-2}}}\right)$',
                            color=_plot.colorpalette[2],
                            linestyle='--', log='xy')

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, with_cbar=False,
                    xlabel=r'$N_\mathsf{H}\ [\mathsf{cm^{-2}}]$', ylabel=r'$B_\mathsf{LOS}\ [\si{\micro G}]$')
        # Plot quantity to velocity map
        plot.plot_line(xdata=column_density_xdata, ydata=derived_magnetic_field_ydata, log='xy', marker='.',
                       linestyle='', label='$B_\mathsf{LOS, derived}$', color=plot.colorpalette[0])
        plot_mass_to_flux_line(plot)
        # Plot the legend
        plot.plot_legend()
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, with_cbar=False,
                    xlabel=r'$N_\mathsf{H}\ [\mathsf{cm^{-2}}]$', ylabel=r'$B_\mathsf{LOS}\ [\si{\micro G}]$')
        # Plot quantity to velocity map
        plot.plot_line(xdata=column_density_xdata, ydata=reference_magnetic_field_ydata, log='xy', marker='.',
                       linestyle='', label='$B_\mathsf{LOS, reference}$', color='green')
        plot_mass_to_flux_line(plot)
        # Plot the legend
        plot.plot_legend()
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_40(self):
        """Plot histogram of B_total to B_LOS.
        """
        # Load necessary modules
        from sys import stdout
        self.file_io.init_plot_output('mag_field_total_to_los')
        # Read velocity channel data
        tbldata, header = self.file_io.read_vel_maps(
            'vel_channel_maps_species_0001_line_0001')
        # Get gas module
        gas_chooser = GasChooser(self.file_io, self.parse_args)
        gas = gas_chooser.get_module_from_name(header['species_name'])
        # Calculate velocity width of a channel
        vel_channel_width = 2. * \
            header['max_velocity'] / (header['nr_channels'] - 1)
        # Initialise output arrays for the derived and reference magnetic field strength
        derived_magnetic_field = np.zeros(
            header['nr_pixel_x'] * header['nr_pixel_y'])
        reference_magnetic_field = np.zeros(
            header['nr_pixel_x'] * header['nr_pixel_y'])
        # Initialise data for LOS and total magnetic field strengths
        total_magnetic_field = np.zeros(
            header['nr_pixel_x'] * header['nr_pixel_y'])
        # Initialise data for LOS and total magnetic field strengths
        for i_x in range(header['nr_pixel_x']):
            for i_y in range(header['nr_pixel_y']):
                stdout.write('--- Calculate magnetic field strength for each pixel: ' +
                             str(int(100. * (i_x * header['nr_pixel_y'] + i_y + 1) /
                                     (header['nr_pixel_x'] * header['nr_pixel_y']))) + ' % \r')
                stdout.flush()
                # Calculate velocity shift by comparing I and V profiles
                velocity_shift = self.math.velocity_shift_from_profile(
                    tbldata[0, :, i_x, i_y], tbldata[3, :, i_x, i_y], vel_channel_width)
                # Get magnetic field strength from Zeeman information of used species
                derived_magnetic_field[i_x, i_y] = gas.shift_2_mag(velocity_shift, header['frequency'],
                                                                   header['i_transition']) * 1e10
                # Channel number 0 contains the LOS magnetic field strength of the model
                reference_magnetic_field[i_x * header['nr_pixel_y'] +
                                         i_y] = tbldata[4, 0, i_x, i_y] * 1e10
                # Channel number 1 contains the total magnetic field strength of the model
                total_magnetic_field[i_x * header['nr_pixel_y'] +
                                     i_y] = tbldata[4, 1, i_x, i_y] * 1e10

        mean_total_derived = np.asscalar(
            np.mean(np.divide(total_magnetic_field, derived_magnetic_field)))
        mean_total_reference = np.asscalar(
            np.mean(np.divide(total_magnetic_field, reference_magnetic_field)))

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, extent=[0, 100, None, None], with_cbar=False,
                    xlabel=r'$\frac{B_\mathsf{total}}{B_\mathsf{LOS, derived}}$', ylabel=r'$N$')
        # Plot quantity to velocity map
        plot.plot_hist(data=np.divide(total_magnetic_field, derived_magnetic_field), log=True, hist_bins=100,
                       color=plot.colorpalette[0], bin_range=[0, 100],
                       label=r'$\mathsf{average}=' + str(round(mean_total_derived, 2)) + '$')
        # Plot the legend
        plot.plot_legend()
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, extent=[0, 100, None, None], with_cbar=False,
                    xlabel=r'$\frac{B_\mathsf{total}}{B_\mathsf{LOS, reference}}$', ylabel=r'$N$')
        # Plot quantity to velocity map
        plot.plot_hist(data=np.divide(total_magnetic_field, reference_magnetic_field), log=True, hist_bins=100,
                       color='green', bin_range=[0, 100],
                       label=r'$\mathsf{average}=' + str(round(mean_total_reference, 2)) + '$')
        # Plot the legend
        plot.plot_legend()
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_41(self):
        """Plot schematic illustration of two Zeeman split spectral lines
        with either low or high magnetic field strengths.
        """
        self.file_io.init_plot_output(
            'derive_mag_field_spectrum_1', path=self.file_io.path['model'])
        simulation_names = ['illustration_low_mag_field',
                            'illustration_very_high_mag_field']
        max_value = 1.
        for i_mag_field in range(2):
            # Create Matplotlib figure
            plot = Plot(self.model, self.parse_args, xlabel=r'$\mathsf{Frequency}$', ylabel=r'$\mathsf{Flux}$',
                        extent=[-3., 3., None, None], with_cbar=False)
            # Set paths of each simulation
            self.file_io.set_path_from_str(
                'plot', 'cube_test', simulation_names[i_mag_field], 'line')
            # Read spectrum data
            tbldata = self.file_io.read_data_file(
                'line_spectrum_species_0001_line_0002.det', skip_header=21)
            tmp_xdata = np.multiply(tbldata[:, 0], 1e-3)
            tmp_ydata = tbldata[:, 1]
            if i_mag_field == 0:
                max_value = np.max(tmp_ydata)
            tmp_ydata /= max_value

            if i_mag_field == 0:
                # Plot spectrum as line
                plot.plot_line(tmp_xdata, tmp_ydata, color=plot.colorpalette[0],
                               label=r'$I$', no_ticks=True)
                # Plot the legend
                # plot.plot_legend()
                # Create Matplotlib figure

                # Plot spectrum as line
                # plot.plot_line(tmp_xdata, tmp_ydata, color=plot.colorpalette[0], label=r'$I (B=0)$')
            elif i_mag_field == 1:
                tmp_ydata_sp = tbldata[:, 4]
                tmp_ydata_sm = tmp_ydata_sp.copy()

                n = len(tmp_ydata[:])
                max_value_sigma = np.max(
                    tbldata[:, 4]) / np.max(tmp_ydata[0:int(n / 2.5)])

                tmp_ydata_sp[np.where(tmp_ydata_sp < 0.)] = 0.
                tmp_ydata_sp /= max_value_sigma
                tmp_ydata_sm[np.where(tmp_ydata_sm > 0.)] = 0.
                tmp_ydata_sm /= max_value_sigma

                # tmp_ydata_pi = np.subtract(tmp_ydata, np.add(abs(tmp_ydata_sp), abs(tmp_ydata_sm)))
                # Plot spectrum as line
                plot.plot_line(tmp_xdata, tmp_ydata, color=plot.colorpalette[2],
                               label=r'$I$', no_ticks=True)
                # plot.plot_line(tmp_xdata, abs(tmp_ydata_sp), color=plot.colorpalette[2], label=r'$\sigma_+$',
                #               no_ticks=True)
                # plot.plot_line(tmp_xdata, abs(tmp_ydata_pi), color='orange', label=r'$\pi$', no_ticks=True)
                # plot.plot_line(tmp_xdata, abs(tmp_ydata_sm), color=plot.colorpalette[0], label=r'$\sigma_-$',
                #               no_ticks=True)
                # plot.plot_double_arrow(arrow_origin=[-0.9, 0.57], arrow_offset=[1.8, 0],
                #                       color=plot.colorpalette[2])
                # plot.plot_text(text='${\propto}B$', text_pos=[0., 0.6], color=plot.colorpalette[2])
            # Plot the legend
            # plot.plot_legend()
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)

    def plot_42(self):
        """Plot schematic illustration of two Zeeman split spectral lines
        with either low or high magnetic field strengths.
        """
        self.file_io.init_plot_output(
            'derive_mag_field_spectrum_2', path=self.file_io.path['model'])
        simulation_names = ['illustration_low_mag_field',
                            'illustration_high_mag_field']

        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, with_cbar=False,
                    xlabel=[r'$\mathsf{Geschwindigkeit}$',
                            r'$\mathsf{Geschwindigkeit}$'],
                    ylabel=[r'$\mathsf{Strahlungsfluss}$', r'$\mathsf{Strahlungsfluss}$'], nr_x_images=2, nr_y_images=2,
                    extent=[-2.5, 2.5, None, None], language='german')

        for i_mag_field in range(2):
            # Set paths of each simulation
            self.file_io.set_path_from_str(
                'plot', 'cube_test', simulation_names[i_mag_field], 'line')
            # Read spectrum data
            tbldata = self.file_io.read_data_file(
                'line_spectrum_species_0001_line_0001.det', skip_header=21)

            for i_quantity in range(0, 4, 3):
                tmp_xdata = np.multiply(tbldata[:, 0], 1e-3)
                tmp_ydata = tbldata[:, i_quantity + 1]
                tmp_ydata /= np.max(tmp_ydata)
                if i_quantity == 0:
                    # Plot spectrum as line
                    plot.plot_line(ax_index=i_mag_field, xdata=tmp_xdata, ydata=tmp_ydata, color=plot.colorpalette[0],
                                   label=r'$I$')
                elif i_quantity == 3:
                    # Plot spectrum as line
                    plot.plot_line(ax_index=2 + i_mag_field, xdata=tmp_xdata, ydata=tmp_ydata, color='#ff8000',
                                   label=r'$V$')
                channel_width = tmp_xdata[1] - tmp_xdata[0]
                tmp_ydata = np.gradient(tbldata[:, 1], channel_width)

                # Define fit function to get proportionality factor
                def fit_function(vch, prop_factor):
                    return prop_factor * tmp_ydata[vch]

                # Calculate best-fit parameter
                popt, pcov = curve_fit(fit_function, range(len(tbldata[:, 0])),
                                       tbldata[:, 4])
                tmp_ydata *= popt[0]
                tmp_ydata /= np.max(tmp_ydata)
                tmp_ydata *= 0.2
                plot.plot_line(ax_index=int(i_quantity / 3) * 2 + i_mag_field, xdata=tmp_xdata, ydata=tmp_ydata,
                               linestyle='--',
                               color='black', label=r'$\mathsf{d}I/\mathsf{d}\nu$')
                # Plot proportional arrows to visualise the relation between dI/dv and V
                if i_quantity == 3:
                    if i_mag_field == 0:
                        plot.plot_annotate(ax_index=int(i_quantity / 3) * 2 + i_mag_field,
                                           text=r'${\propto}B_\mathsf{Sicht}$',
                                           text_pos=[-0.71, 0.30], pos=[-0.71, 0.85], color=plot.colorpalette[1])
                        plot.plot_annotate(ax_index=int(i_quantity / 3) * 2 + i_mag_field,
                                           text='', text_pos=[0.71, -0.26], pos=[0.71, -0.8],
                                           color=plot.colorpalette[1])
                    elif i_mag_field == 1:
                        plot.plot_annotate(ax_index=int(i_quantity / 3) * 2 + i_mag_field,
                                           text=r'${\not\propto}B_\mathsf{Sicht}$',
                                           text_pos=[-0.81, 0.30], pos=[-0.81, 0.85], color=plot.colorpalette[2])
                        plot.plot_annotate(ax_index=int(i_quantity / 3) * 2 + i_mag_field,
                                           text='', text_pos=[0.81, -0.25], pos=[0.81, -0.8],
                                           color=plot.colorpalette[2])
                # Plot the legend
                if int(i_quantity / 3) * 2 + i_mag_field in [0, 1]:
                    plot.plot_legend(ax_index=int(i_quantity / 3)
                                     * 2 + i_mag_field, loc='upper right')
                else:
                    plot.plot_legend(ax_index=int(
                        i_quantity / 3) * 2 + i_mag_field)
        plot.axarr[0].set_ylim([-0.25, 1.25])
        plot.axarr[1].set_ylim([-0.25, 1.25])
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_43(self):
        """Plot raytrace results from POLARIS simulations.
        """
        self.file_io.init_plot_output('subpixel_test')
        tbldata = []
        max_value = []
        for pixel in [256]:
            for theta in [90]:
                for zoom in [1]:
                    for subpixel in [0, 3]:
                        simulation_fits_filename = str(pixel) + '_pixel_' + str(subpixel) + '_subpixel_' \
                            + str(theta) + '_theta_' + str(zoom) + '_zoom.fits'
                        # Read raytrace results from file
                        raytrace_data, header = self.file_io.read_emission_map(
                            simulation_fits_filename)
                        # Take data for current quantity
                        tbldata.append(raytrace_data[0, :, :])
                    max_value.append(
                        max(np.max(tbldata[-1]), np.max(tbldata[-2])))
                    tbldata.append(abs(np.subtract(tbldata[-1], tbldata[-2])))
                    # tbldata.append(abs(100. * np.divide(np.subtract(tbldata[-1], tbldata[-2]),
                    #                                     np.add(tbldata[-1], tbldata[-2]))))
        for i_data in range(len(tbldata)):
            # Create Matplotlib figure
            plot = Plot(self.model, self.parse_args)
            # Take colorbar label from quantity id
            if i_data % 3 == 2:
                cbar_label = r'$\Delta$' + self.file_io.get_quantity_labels(0)
                # r'$\Delta F/F\ [\si{\percent}]$'
                plot.plot_imshow(
                    tbldata[i_data], cbar_label=cbar_label, vmin=0, extend='neither')
            else:
                cbar_label = self.file_io.get_quantity_labels(0)
                plot.plot_imshow(tbldata[i_data], cbar_label=cbar_label, vmax=max_value[int(i_data / 3)], vmin=0.,
                                 extend='neither')
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)

    def plot_66(self):
        """Plot total SED from POLARIS simulations.
        """
        # Set data input to Jy/px to calculate the total flux
        if self.parse_args.cmap_unit is None:
            self.file_io.cmap_unit = 'total'
        # Create pdf file if show_plot is not chosen
        self.file_io.init_plot_output(
            'polaris_mcfost_sed', path=self.file_io.path['simulation'])
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, xlabel=r'$\lambda\ [\si{\metre}]$',
                    ylabel=r'$\mathit{I}\ [\si{Jy}]$', with_cbar=False)

        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', self.parse_args.model_name, self.parse_args.simulation_name, 'dust')
        # Read raytrace results from file
        ray_data, ray_header = self.file_io.read_emission_sed(
            'polaris_detector_nr0001_sed')
        # Plot spectral energy distribution
        # plot.plot_line(ray_header['wavelengths'], ray_data[0, :], log='xy', label=r'$\mathsf{POLARIS\ thermal\ emission}$')

        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', self.parse_args.model_name, self.parse_args.simulation_name, 'dust_mc')
        # Read raytrace results from file
        mc_data, mc_header = self.file_io.read_emission_sed(
            'polaris_detector_nr0001_sed')
        # Plot spectral energy distribution (direct)
        # plot.plot_line(mc_header['wavelengths'], mc_data[6, :], log='xy', label=r'$\mathsf{POLARIS\ direct\ starlight}$')
        # Plot spectral energy distribution (scattered)
        # plot.plot_line(mc_header['wavelengths'], mc_data[7, :], log='xy', label=r'$\mathsf{POLARIS\ scattered\ starlight}$')

        wavelengths_total = []
        quantity_total = []
        offset = 0
        for i_wl in range(mc_header['nr_wavelengths'] + ray_header['nr_wavelengths']):
            if i_wl - offset >= ray_header['nr_wavelengths']:
                break
            elif i_wl >= mc_header['nr_wavelengths']:
                wavelengths_total.append(
                    ray_header['wavelengths'][i_wl - offset])
                quantity_total.append(ray_data[0, i_wl - offset])
            elif mc_header['wavelengths'][i_wl] == ray_header['wavelengths'][i_wl - offset]:
                wavelengths_total.append(ray_header['wavelengths'][i_wl])
                quantity_total.append(
                    mc_data[0, i_wl] + ray_data[0, i_wl - offset])
            elif mc_header['wavelengths'][i_wl] < ray_header['wavelengths'][i_wl - offset]:
                offset += 1
                wavelengths_total.append(mc_header['wavelengths'][i_wl])
                quantity_total.append(mc_data[0, i_wl])
        # Plot total spectral energy distribution
        plot.plot_line(wavelengths_total, quantity_total, log='xy', linestyle='--',
                       label=r'$\mathsf{POLARIS\ total\ SED}$')

        wavelengths = ray_header['wavelengths']
        from astropy.io import fits
        mcfost_path = '/home/rbrauer/astrophysics/mcfost/data_th/'

        def to_jansky(wl_list, data):
            res_data = np.zeros(np.shape(data))
            for i, wl in enumerate(wl_list):
                frequency = self.math.const['c'] / wl_list[i]
                res_data[i] = data[i] * 1e26 / frequency
            return res_data

        hdulist = fits.open(mcfost_path + 'sed_rt.fits.gz')
        mcfost_total_sed = to_jansky(wavelengths, hdulist[0].data[0, 0, 0, :])
        # mcfost_direct_starlight = to_jansky(wavelengths, hdulist[0].data[4, 0, 0, :])
        # mcfost_scattered_starlight = to_jansky(wavelengths, hdulist[0].data[5, 0, 0, :])
        # mcfost_thermal_reemission = to_jansky(wavelengths, hdulist[0].data[6, 0, 0, :])

        # Plot spectral energy distribution
        # plot.plot_line(wavelengths, mcfost_thermal_reemission, log='xy',
        #               label=r'$\mathsf{MCFOST\ thermal\ emission}$')
        # Plot spectral energy distribution (direct)
        # plot.plot_line(wavelengths, mcfost_direct_starlight, log='xy',
        #               label=r'$\mathsf{MCFOST\ direct\ starlight}$')
        # Plot spectral energy distribution (scattered)
        # plot.plot_line(wavelengths, mcfost_scattered_starlight, log='xy',
        #               label=r'$\mathsf{MCFOST\ scattered\ starlight}$')
        # Plot total spectral energy distribution
        plot.plot_line(wavelengths, mcfost_total_sed, log='xy', linestyle='--',
                       label=r'$\mathsf{MCFOST\ total\ SED}$')

        # hdulist = fits.open(mcfost_path + 'sed_mc.fits.gz')
        # mcfost_total_sed = to_jansky(wavelengths, hdulist[0].data[0, 0, 0, :])
        # plot.plot_line(wavelengths, mcfost_total_sed, log='xy', linestyle='--',
        #                label=r'$\mathsf{MCFOST\ total\ SED\ (mc)}$')

        # Plot the legend
        plot.plot_legend()
        # adapt limits
        plot.set_limits(limits=[None, None, 1e-10,
                                np.max(quantity_total) * 1e1])
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_67(self):
        """Plot total SED from POLARIS simulations.
        """
        # Set data input to Jy/px to calculate the total flux
        if self.parse_args.cmap_unit is None:
            self.file_io.cmap_unit = 'total'
        # Create pdf file if show_plot is not chosen
        self.file_io.init_plot_output(
            'polaris_mcfost_dust', path=self.file_io.path['simulation'])
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, xlabel=r'$\lambda\ [\si{\metre}]$',
                    ylabel=r'$\mathit{\kappa}\ [\si{m^2/kg}]$', with_cbar=False)

        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', self.parse_args.model_name, self.parse_args.simulation_name, 'dust')
        # Read raytrace results from file
        ray_data, ray_header = self.file_io.read_emission_sed(
            'polaris_detector_nr0001_sed')
        # Plot spectral energy distribution
        # plot.plot_line(ray_header['wavelengths'], ray_data[0, :], log='xy', label=r'$\mathsf{POLARIS\ thermal\ emission}$')

        wavelengths = ray_header['wavelengths']
        from astropy.io import fits

        mcfost_dust_path = '/home/rbrauer/astrophysics/mcfost/data_dust/'
        hdulist = fits.open(mcfost_dust_path + 'kappa.fits.gz')
        plot.plot_line(wavelengths, hdulist[0].data[:] * 1e-4 * 1e3, log='xy', linestyle='--',
                       label=r'$\mathsf{MCFOST\ kappa}$')

        avg_grain_mass = 4.06951e-18  # 5.55238e-20
        polaris_kappa = np.genfromtxt(
            self.file_io.path['results'] + 'dust_mixture.dat', skip_header=30)
        plot.plot_line(polaris_kappa[:, 0], polaris_kappa[:, 1] / avg_grain_mass, log='xy', linestyle='-',
                       label=r'$\mathsf{POLARIS\ kappa}$')

        # Plot the legend
        plot.plot_legend()
        # adapt limits
        # plot.set_limits(limits=[None, None, 1e-10, np.max(quantity_total) * 1e1])
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_68(self):
        """Plot Raytrace or Monte-Carlo results as map from POLARIS simulations.
        """
        # Create pdf file if show_plot is not chosen
        self.file_io.init_plot_output(
            'gg_tau_planet_emission_map', path=self.file_io.path['simulation'])
        # Read raytrace results from file
        map_data, header = self.file_io.read_emission_map(
            'polaris_detector_nr0001')
        # Create one plot per component of the simulation
        if self.parse_args.simulation_type == 'dust':
            quantity_list = [0, 6]
            contrast_factor = 1e-4
        else:
            quantity_list = [0]
            contrast_factor = 1e-5
        for i_quantity in quantity_list:
            # Take colorbar label from quantity id
            cbar_label = self.file_io.get_quantity_labels(i_quantity)
            for i_wl, wl in enumerate(header['wavelengths']):
                if wl not in [5.55047321e-06, 7.6522007e-06, 1.05497637e-05, 1.30681353e-05,
                              1.61876774e-05, 2.0051897e-05, 2.23172626e-05]:
                    continue
                map_title = r'$\lambda=\SI{' + \
                    str(self.math.latex_float(wl * 1e6)) + '}{\micro\metre}$'
                # Take data for current quantity
                tbldata = map_data[i_quantity, i_wl, :, :]
                # skip zero content plots
                if not tbldata.any():
                    continue
                # Create Matplotlib figure
                plot = Plot(self.model, self.parse_args)
                # Show title to know the wavelength
                plot.plot_title(map_title)
                if i_quantity == 0:
                    plot.plot_imshow(tbldata, cbar_label=cbar_label, set_bad_to_min=True, norm='LogNorm',
                                     vmin=contrast_factor * np.max(tbldata), vmax=np.max(tbldata), extend='min')
                elif i_quantity == 6:
                    plot.plot_imshow(
                        tbldata, cbar_label=cbar_label, set_bad_to_min=True)
                # Save figure to pdf file or print it on screen
                plot.save_figure(self.file_io)

    def plot_69(self):
        """Plot Raytrace and Monte-Carlo results as SED from POLARIS simulations.
        """
        # Set data input to Jy/px to calculate the total flux
        if self.parse_args.cmap_unit is None:
            self.file_io.cmap_unit = 'total'
        # Plot only the intensity as a full sed plot
        i_quantity = 0
        # Create pdf file if show_plot is not chosen
        self.file_io.init_plot_output('polaris_detector_nr0001_sed')
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, xlabel=r'$\lambda\ [\si{\metre}]$',
                    ylabel=self.file_io.get_quantity_labels(0), with_cbar=False)
        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', self.parse_args.model_name, self.parse_args.simulation_name, 'dust')
        # Read raytrace results from file
        ray_data, ray_header = self.file_io.read_emission_sed(
            'polaris_detector_nr0001_sed')
        # Plot spectral energy distribution
        plot.plot_line(ray_header['wavelengths'], ray_data[i_quantity,
                                                           :], label=r'$\mathsf{thermal\ emission}$')

        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', self.parse_args.model_name, self.parse_args.simulation_name, 'dust_mc')
        # Read raytrace results from file
        mc_data, mc_header = self.file_io.read_emission_sed(
            'polaris_detector_nr0001_sed')
        # Plot spectral energy distribution (direct)
        plot.plot_line(mc_header['wavelengths'], mc_data[6, :],
                       label=r'$\mathsf{direct\ starlight}$')
        # Plot spectral energy distribution (scattered)
        plot.plot_line(mc_header['wavelengths'], mc_data[7, :],
                       label=r'$\mathsf{scattered\ starlight}$')
        wavelengths_total = []
        quantity_total = []
        offset = 0
        for i_wl in range(mc_header['nr_wavelengths'] + ray_header['nr_wavelengths']):
            if i_wl - offset >= ray_header['nr_wavelengths']:
                break
            elif i_wl >= mc_header['nr_wavelengths']:
                wavelengths_total.append(
                    ray_header['wavelengths'][i_wl - offset])
                quantity_total.append(ray_data[i_quantity, i_wl - offset])
            elif mc_header['wavelengths'][i_wl] == ray_header['wavelengths'][i_wl - offset]:
                wavelengths_total.append(mc_header['wavelengths'][i_wl])
                quantity_total.append(
                    mc_data[i_quantity, i_wl] + ray_data[0, i_wl - offset])
            elif mc_header['wavelengths'][i_wl] < ray_header['wavelengths'][i_wl - offset]:
                offset += 1
                wavelengths_total.append(mc_header['wavelengths'][i_wl])
                quantity_total.append(mc_data[i_quantity, i_wl])
        # Plot total spectral energy distribution
        plot.plot_line(np.array(wavelengths_total), np.array(quantity_total), linestyle='--',
                       label=r'$\mathsf{total\ SED}$')
        # Fit SED with polfit
        from scipy.interpolate import interp1d
        f = interp1d(np.array(wavelengths_total),
                     np.array(quantity_total), kind='cubic')
        wl_array = np.multiply(np.array([3.5, 4.5, 5.6, 7.7, 10.5]), 1e-6)
        x = np.linspace(wavelengths_total[0], wavelengths_total[-1], 100)
        plot.plot_line(x, f(x), linestyle=':')
        print(wl_array)
        print(f(wl_array))
        # Plot the legend
        plot.plot_legend()
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_70(self):
        """Plot Raytrace or Monte-Carlo results as SED from POLARIS simulations.
        """
        # Set data input to Jy/px to calculate the total flux
        if self.parse_args.cmap_unit is None:
            self.file_io.cmap_unit = 'total'
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, xlabel=r'$\lambda\ [\si{\metre}]$',
                    ylabel=self.file_io.get_quantity_labels(0), with_cbar=False)
        for i in [1, 6, 8]:
            sed_data, header, _ = self.file_io.read_emission_sed('polaris_detector_nr' +
                                                                 str(i).zfill(4) + '_sed')
            # Create pdf file if show_plot is not chosen
            self.file_io.init_plot_output('polaris_detector_compare_sed')
            # Plot spectral energy distribution
            plot.plot_line(header['wavelengths'],
                           sed_data[0, :], log='xy', marker='.')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_71(self):
        """Plot P_l over tau from SED results.
        """
        # Set data input to Jy/px to calculate the total flux
        if self.parse_args.cmap_unit is None:
            self.file_io.cmap_unit = 'total'
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, xlabel=r'$\lambda\ [\si{\metre}]$',
                    ylabel=r'$\frac{\mathit{P}_\mathsf{l}}{\tau}\ [\%]$', with_cbar=False)
        sed_data, header, _ = self.file_io.read_emission_sed(
            'polaris_detector_nr0001_sed')
        # Create pdf file if show_plot is not chosen
        self.file_io.init_plot_output('polaris_detector_p_over_tau')
        # Plot spectral energy distribution
        p_over_l = np.divide(sed_data[5, :], sed_data[6, :])
        plot.plot_line(header['wavelengths'], p_over_l, log='xy', marker='.')
        plot.plot_line([0.55e-6, 0.55e-6], [p_over_l.min(),
                                            p_over_l.max()], log='xy', marker='.')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    # ------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------
    # In the following section are plotting routines related to different paper/poster/talks/proposal.
    # Index numbering is coded as follows: plot_XYYYZZZ
    # -----------------------------------------------------
    # X -> 1: paper, 2: talk, 3: poster, 4: proposal
    # Y -> Index of paper/talk/poster/proposal
    # Z -> Index of individual plot
    # e.g.  def plot_1001001(self):
    # ------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------
    # -------------------------------------- WORK IN PROGRESS ----------------------------------------
    # ------------------------------------------------------------------------------------------------

    def plot_2005001(self):
        """Plot midplane data from POLARIS simulations (from MHD simulations).
        """
        # Set vector size to match with 2048 x 2048 pixel sized image
        self.file_io.init_plot_output(
            'talk_bordeaux_midplane_cuts', path=self.file_io.path['simulation'])
        visualization_input_list = ['input_dust_number_density_xy']
        # For each item in visualization_input_list create a separate figure
        for visualization_input in visualization_input_list:
            # Read midplane data including a vector field
            tbldata, header, vec_field_data = self.file_io.read_midplane_file(
                visualization_input)
            if tbldata is not None and vec_field_data is not None:
                # Update the axes label to fit the axes used in the midplane file
                if 'xy' in visualization_input:
                    label_plane = 'xy'
                elif 'xz' in visualization_input:
                    label_plane = 'xz'
                elif 'yz' in visualization_input:
                    label_plane = 'yz'
                else:
                    raise ValueError(
                        'The chosen midplane file has no valid cut through a plane (xy, xz, yz?)!')
                # Create Matplotlib figure
                plot = Plot(self.model, self.parse_args,
                            label_plane=label_plane)
                # Limit data to reasonable values
                tbldata[np.where(tbldata <= 1e-30)] = 0
                # Plot midplane data depending on quantity derived from filename
                self.basic_plots.plot_midplane_map_base(
                    visualization_input, plot, tbldata, vec_field_data)
                # Load stellar source to get position of the binary stars
                from modules.source import SourceChooser
                stellar_source_chooser = SourceChooser(
                    self.model, self.file_io, self.parse_args)
                stellar_source = stellar_source_chooser.get_module_from_name(
                    'gg_tau_binary')
                # Plot position of binary stars (with conversion from m to au)
                stellar_source.tmp_parameter['position_star'][0] = np.divide(stellar_source.tmp_parameter['position_star'][0],
                                                                             self.math.const['au'])
                plot.plot_text(
                    stellar_source.tmp_parameter['position_star'][0], r'$+$')
                stellar_source.tmp_parameter['position_star'][1] = np.divide(stellar_source.tmp_parameter['position_star'][1],
                                                                             self.math.const['au'])
                plot.plot_text(
                    stellar_source.tmp_parameter['position_star'][1], r'$+$')
                # Save figure to pdf file or print it on screen
                plot.save_figure(self.file_io)

    def plot_2005002(self):
        """Plot Raytrace and Monte-Carlo results as SED from POLARIS simulations.
                """
        # Set data input to Jy/px to calculate the total flux
        if self.parse_args.cmap_unit is None:
            self.file_io.cmap_unit = 'total'
        # Plot only the intensity as a full sed plot
        i_quantity = 0
        # Create pdf file if show_plot is not chosen
        self.file_io.init_plot_output('sed_plus_vizir')
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, xlabel=r'$\lambda\ [\si{\metre}]$',
                    ylabel=self.file_io.get_quantity_labels(0), with_cbar=False)
        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', self.parse_args.model_name, self.parse_args.simulation_name, 'dust')
        # Read raytrace results from file
        ray_data, ray_header = self.file_io.read_emission_sed(
            'polaris_detector_nr0001_sed')
        # Plot spectral energy distribution
        plot.plot_line(ray_header['wavelengths'], ray_data[i_quantity, :], log='y',
                       label=r'$\mathsf{thermal\ emission}$')

        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', self.parse_args.model_name, self.parse_args.simulation_name, 'dust_mc')
        # Read raytrace results from file
        mc_data, mc_header = self.file_io.read_emission_sed(
            'polaris_detector_nr0001_sed')
        # Plot spectral energy distribution (direct)
        plot.plot_line(mc_header['wavelengths'], mc_data[6, :],
                       log='y', label=r'$\mathsf{direct\ starlight}$')
        # Plot spectral energy distribution (scattered)
        plot.plot_line(mc_header['wavelengths'], mc_data[7, :],
                       log='y', label=r'$\mathsf{scattered\ starlight}$')
        wavelengths_total = []
        quantity_total = []
        offset = 0
        for i_wl in range(mc_header['nr_wavelengths'] + ray_header['nr_wavelengths']):
            if i_wl - offset >= ray_header['nr_wavelengths']:
                break
            elif i_wl >= mc_header['nr_wavelengths']:
                wavelengths_total.append(
                    ray_header['wavelengths'][i_wl - offset])
                quantity_total.append(ray_data[i_quantity, i_wl - offset])
            elif mc_header['wavelengths'][i_wl] == ray_header['wavelengths'][i_wl - offset]:
                wavelengths_total.append(mc_header['wavelengths'][i_wl])
                quantity_total.append(
                    mc_data[i_quantity, i_wl] + ray_data[0, i_wl - offset])
            elif mc_header['wavelengths'][i_wl] < ray_header['wavelengths'][i_wl - offset]:
                offset += 1
                wavelengths_total.append(mc_header['wavelengths'][i_wl])
                quantity_total.append(mc_data[i_quantity, i_wl])
        # Plot total spectral energy distribution
        plot.plot_line(np.array(wavelengths_total), np.array(quantity_total), log='xy', linestyle='--',
                       label=r'$\mathsf{total\ SED}$')

        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', self.parse_args.model_name, self.parse_args.simulation_name, 'dust_mc')
        # Plot vizier data
        from astropy.io.votable import parse_single_table
        table = parse_single_table(
            self.file_io.path['results'] + 'vizier_votable.vot')
        data_flux = table.array['sed_flux']
        data_flux_error = table.array['sed_eflux'].filled(0)
        data_wl = self.math.const['c'] / (table.array['sed_freq'] * 1e9)
        # Plot total spectral energy distribution
        plot.plot_line(data_wl, data_flux, yerr=data_flux_error, log='xy', linestyle='none', marker='.', color='purple',
                       alpha=0.7, label=r'$\mathsf{VizieR\ SED}$')

        # Plot the legend
        plot.plot_legend()
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_1006001(self):
        """Plot zoom in midplane of density distribution of GG Tau disk
        """
        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', 'gg_tau_disk', 'default', 'temp')
        # Set output filename
        self.file_io.init_plot_output(
            'zoom_in_gg_tau_disk', path=self.file_io.path['model'])
        # Set midplane type
        visualization_input = 'input_dust_mass_density_xy'
        # Read midplane data (main image)
        [tbldata, vec_field_data], _, _ = self.file_io.read_midplane_file(visualization_input,
                                                                          filename='input_midplane_full.fits')
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, label_plane='xy',
                    cmap_scaling=['log'])
        # Plot midplane data depending on quantity derived from filename
        self.basic_plots.plot_midplane_map_base(visualization_input, plot, tbldata, vec_field_data,
                                                vmin=1e-14, set_bad_to_min=True, cmap='inferno')
        # Read midplane data (zoomed image, changes extent of self.model)
        [tbldata_zoom, vec_field_data_zoom], _, _ = self.file_io.read_midplane_file(
            visualization_input, filename='input_midplane_zoom.fits')
        # Create zoom plot (set model to provide extent of zoomed image)
        plot.create_zoom_axis(model=self.model, zoom_factor=3.5)
        # Plot zoomed image
        self.basic_plots.plot_midplane_map_base(visualization_input, plot, tbldata_zoom, vec_field_data_zoom,
                                                ax_index=1, plot_cbar=False, vmin=1e-14, set_bad_to_min=True, cmap='inferno')
        # Load stellar source to get position of the binary stars
        from modules.source import SourceChooser
        stellar_source_chooser = SourceChooser(
            self.model, self.file_io, self.parse_args)
        stellar_source = stellar_source_chooser.get_module_from_name(
            'gg_tau_stars')
        # Plot position of binary stars (with conversion from m to au)
        star_descr = ['Aa', 'Ab1', 'Ab2']
        for i_star, position in enumerate(stellar_source.tmp_parameter['position_star']):
            # plot.plot_marker(
            #    np.divide(position[0:2], self.math.const['au']),
            #    "*", ax_index=0, color='white', markersize=12)
            if i_star == 0:
                offset = [0, 10]
            elif i_star == 1:
                offset = [-10, 0]
            elif i_star == 2:
                offset = [10, 0]
            plot.plot_text(
                np.add(np.divide(position[0:2],
                                 self.math.const['au']), offset),
                text=r'$\text{' + star_descr[i_star] + r'}$', ax_index=1, color='white')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_1006002(self):
        """Plot GG Tau A SED comparison
        """
        # Set data input to Jy/px to calculate the total flux
        if self.parse_args.cmap_unit is None:
            self.file_io.cmap_unit = 'total'
        # Set output filename
        self.file_io.init_plot_output(
            'sed_comparison', path=self.file_io.path['model'])
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, xlabel=r'$\lambda\ [\si{\metre}]$',
                    ylabel=self.file_io.get_quantity_labels(0), with_cbar=False,
                    limits=[None, None, 1e-10, None])
        # Plot vizier data
        from astropy.io.votable import parse_single_table
        table = parse_single_table(
            self.file_io.path['model'] + 'vizier_votable.vot')
        data_flux = table.array['sed_flux']
        data_flux_error = table.array['sed_eflux'].filled(0)
        data_wl = self.math.const['c'] / (table.array['sed_freq'] * 1e9)
        # Plot total spectral energy distribution
        plot.plot_line(data_wl, data_flux, yerr=data_flux_error, log='xy',
                       linestyle='none', marker='.', color='black',
                       alpha=0.7, label=r'$\text{Photometric observations}$')
        # Set some variables
        detector_index = 202
        i_quantity = 0
        # Different model configurations
        model_list = [
            'no_circumstellar_disks',
            'only_Aa',
            'only_Ab1',
            'only_Ab2',
            'only_Aa_Ab1',
            'only_Aa_Ab2',
            'only_Ab1_Ab2',
            'default'
        ]
        model_descr = [
            'No circumstellar disks',
            'Disk around Aa',
            'Disk around Ab1',
            'Disk around Ab2',
            'Disks around Aa and Ab1',
            'Disks around Aa and Ab2',
            'Disks around Ab1 and Ab2',
            'Disks around all stars'
        ]
        # Loop over model configurations
        for i_model, model_name in enumerate(model_list):
            if i_model < 8:
                # Set paths of each simulation
                self.file_io.set_path_from_str(
                    'plot', 'gg_tau_disk', model_name, 'dust')
                # Read raytrace results from file
                sed_data, header, _ = self.file_io.read_emission_sed(
                    'polaris_detector_nr' + str(detector_index).zfill(4) + '_sed')
                # Init wavelengths and fluxes
                wavelengths = header['wavelengths']
                quantity = sed_data[i_quantity, :]
                # Plot spectral energy distribution
                plot.plot_line(wavelengths, quantity, log='xy',
                               label=r'$\text{' + model_descr[i_model] + r'}$')
            '''
            else:
                # Set paths of each simulation
                self.file_io.set_path_from_str(
                    'plot', 'gg_tau_disk', model_name, 'dust')
                # Read raytrace results from file
                plot_data, header, plot_data_type = self.file_io.read_emission_map(
                    'polaris_detector_nr' + str(detector_index).zfill(4))
                # Set paths of each simulation (for no disks)
                self.file_io.set_path_from_str(
                    'plot', 'gg_tau_disk', 'no_circumstellar_disks', 'dust')
                # Read raytrace results from file (for no disks)
                plot_data_no, header_no, plot_data_type_no = self.file_io.read_emission_map(
                    'polaris_detector_nr' + str(detector_index).zfill(4))
                # Init wavelengths and fluxes
                wavelengths = header['wavelengths']
                quantity = np.zeros((plot_data.shape[1], 3))
                # Get the image sidelength
                sidelength_x = 2. * self.model.tmp_parameter['radius_x_arcsec']
                sidelength_y = 2. * self.model.tmp_parameter['radius_y_arcsec']
                # Get number of pixel per axis
                nr_pixel_x = plot_data.shape[-2]
                nr_pixel_y = plot_data.shape[-1]
                # Find the considered pixel
                for i_x in range(nr_pixel_x):
                    for i_y in range(nr_pixel_y):
                        if np.sqrt((i_x - 256)**2 + (i_y - 256)**2) < 30:
                            quantity[:, 0] += \
                                plot_data_no[i_quantity, :, i_x, i_y]
                            quantity[:, 1] += plot_data[i_quantity, :, i_x, i_y]
                        else:
                            quantity[:, 2] += plot_data[i_quantity, :, i_x, i_y]

                # Plot spectral energy distribution
                plot.plot_line(wavelengths, quantity[:, 0], log='xy',
                               label=r'$\text{- central stars}$')
                plot.plot_line(wavelengths, np.subtract(quantity[:, 1], quantity[:, 0]), log='xy',
                               label=r'$\text{- circumstellar disks}$')
                plot.plot_line(wavelengths, quantity[:, 2], log='xy',
                               label=r'$\text{- circumbinary disk}$')
            '''
        # Plot the legend
        plot.plot_legend(loc=4)
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_1006003(self):
        """Plot two times 6 polarized intensity emission maps with differnet configurations for GG Tau A.
        """
        def cropND(img, bounding, offset=[0, 0]):
            import operator
            start = tuple(map(lambda a, da, off: a//2 - da //
                              2 + off, img.shape, bounding, offset))
            end = tuple(map(operator.add, start, bounding))
            slices = tuple(map(slice, start, end))
            return img[slices]
        # Import libraries
        from astropy.io import fits
        from scipy.ndimage.interpolation import zoom
        # Set some variables
        detector_index = 101
        i_quantity = 4
        calc = False
        model_list = [
            'no_circumstellar_disks',
            'only_Aa',
            'only_Ab1',
            'only_Ab2',
            'only_Aa_Ab1',
            'only_Aa_Ab2',
            'only_Ab1_Ab2',
            'default'
        ]
        model_descr = [
            'No circumstellar disks',
            'Disk around Aa',
            'Disk around Ab1',
            'Disk around Ab2',
            'Disks around Aa and Ab1',
            'Disks around Aa and Ab2',
            'Disks around Ab1 and Ab2',
            'Disks around all stars'
        ]
        measurement_position_list = [
            [0, -1.05],
            [0, 0.85],
            [-0.57, -0.975],
            [1.17, -0.42],
            # [-1.3, -1.1],
            # [-0.57, 0.82]
        ]
        # Set beam size (in arcsec)
        self.file_io.beam_size = 0.07
        # Take colorbar label from quantity id
        cbar_label = self.file_io.get_quantity_labels(i_quantity)
        # Define output pdf
        self.file_io.init_plot_output(
            '2x2x3_PI_emission_map', path=self.file_io.path['model'])
        # Sum up the flux inside the circle
        flux_sum = np.zeros((2, 6, len(measurement_position_list)))
        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', 'gg_tau_disk', 'default', 'dust')
        # Create pdf file if show_plot is not chosen and read map data from file
        plot_data, header, plot_data_type = self.file_io.read_emission_map(
            'polaris_detector_nr' + str(detector_index).zfill(4))
        # Create two 2x3 plots
        for i_plot in range(2):
            # Create Matplotlib figure
            plot = Plot(self.model, self.parse_args, ax_unit='arcsec',
                        nr_x_images=2, nr_y_images=3)
            for i_subplot in range(6):
                if i_subplot >= 2:
                    # Set paths of each simulation
                    self.file_io.set_path_from_str(
                        'plot', 'gg_tau_disk', model_list[i_subplot - 2 + i_plot * 4], 'dust')
                    # Create pdf file if show_plot is not chosen and read map data from file
                    plot_data, header, plot_data_type = self.file_io.read_emission_map(
                        'polaris_detector_nr' + str(detector_index).zfill(4))
                    # Take data for current quantity
                    tbldata = plot_data[i_quantity, 0, :, :]
                    # Plot map description
                    plot.plot_text(
                        text_pos=[0.03, 0.97], relative_position=True,
                        text=r'$\text{' +
                        model_descr[i_subplot -
                                    2 + i_plot * 4] + r'}$',
                        horizontalalignment='left', verticalalignment='top', ax_index=i_subplot, color='white'
                    )
                elif i_subplot in [0]:
                    hdulist = fits.open(
                        '/home/rbrauer/Documents/projects/005_gg_tau/near_infrared_imaging_paper/sub_pi.fits')
                    tbldata = cropND(
                        hdulist[0].data.T, (500, 500), offset=[0, 30]) / 3e7
                else:
                    continue
                plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=i_subplot, set_bad_to_min=True,
                                 norm='LogNorm', vmin=3e-8, vmax=1e-5, cmap='magma')
                for i_pos, center_pos in enumerate(measurement_position_list):
                    plot.plot_text(text_pos=center_pos,
                                   text=str(i_pos + 1), ax_index=i_subplot, color='black', fontsize=10,
                                   bbox=dict(boxstyle='circle, pad=0.2', facecolor='white', alpha=0.5))
                if calc:
                    # Get the image sidelength
                    sidelength_x = 2. * self.model.tmp_parameter['radius_x_arcsec']
                    sidelength_y = 2. * self.model.tmp_parameter['radius_y_arcsec']
                    # Get number of pixel per axis
                    nr_pixel_x = tbldata.shape[-2]
                    nr_pixel_y = tbldata.shape[-1]
                    # Find the considered pixel
                    for i_x in range(nr_pixel_x):
                        for i_y in range(nr_pixel_y):
                            pos = np.subtract(np.multiply(np.divide(np.add([i_x, i_y], 0.5),
                                                                    [nr_pixel_x, nr_pixel_y]), [sidelength_x, sidelength_y]),
                                            np.divide([sidelength_x, sidelength_y], 2.))
                            for i_pos, center_pos in enumerate(measurement_position_list):
                                pos_r = np.linalg.norm(
                                    np.subtract(pos, center_pos))
                                if pos_r < 0.2:
                                    flux_sum[i_plot, i_subplot,
                                            i_pos] += tbldata[i_x, i_y]
                # Hide second observation plot
                plot.remove_axes(ax_index=1)
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)
        if calc:
            print(r'\begin{tabular}{lcccc}')
            print(r'\theadstart')
            print(
                r'\thead Configuration & \thead $\boldsymbol{\textit{PI}}$ (1) & \thead $\boldsymbol{\textit{PI}}$ (2) & \thead $\boldsymbol{\textit{PI}}$ (3) & \thead $\boldsymbol{\textit{PI}}$ (4) \\')
            print(r'\tbody')
            for i_plot in range(2):
                for i_subplot in range(6):
                    if i_subplot > 1:
                        print(model_descr[i_subplot - 2 + i_plot * 4], '&', ' & '.join(
                            '$\SI{' + f'{x:1.2e}' + '}{}$' for x in flux_sum[i_plot, i_subplot, :]), r'\\')
            print('Observation &', ' & '.join(
                '$\SI{' + f'{x:1.2e}' + '}{}$' for x in flux_sum[0, 0, :]), r'\\')
            print(r'\tend')
            print(r'\end{tabular}')
            print(r'\begin{tabular}{lcccccc}')
            print(r'\theadstart')
            print(
                r'\thead Configuration & \thead $\boldsymbol{\delta \textit{PI}}$ (region 1) & \thead $\boldsymbol{\delta \textit{PI}}$ (region 2) & \thead $\boldsymbol{\delta \textit{PI}}$ (region 3) & \thead $\boldsymbol{\delta \textit{PI}}$ (region 4) & \thead \textit{PI} (region 1) / \textit{PI} (region 2) \\')
            print(r'\tbody')
            for i_plot in range(2):
                for i_subplot in range(6):
                    if i_subplot > 1:
                        print(model_descr[i_subplot - 2 + i_plot * 4], '&', ' & '.join('$\SI{-' + f'{x:1.0f}' + '}{\percent}$' for x in np.subtract(
                            100, 1e2*np.divide(flux_sum[i_plot, i_subplot, :], flux_sum[0, 2, :]))), '&', '& '.join('$\SI{' + f'{x:1.0f}' + '}{\percent}$' for x in [1e2*flux_sum[i_plot, i_subplot, 0] / flux_sum[i_plot, i_subplot, 1]]), r'\\')
            print('Observation & x & x & x & x &', ' & '.join(
                '$\SI{' + f'{x:1.0f}' + '}{\percent}$' for x in [1e2*flux_sum[0, 0, 0]/flux_sum[0, 0, 1]]), r'\\')
            print(r'\tend')
            print(r'\end{tabular}')

    def plot_1006004(self):
        """Plot two times 4 intensity emission maps with differnet configurations for GG Tau A.
        """
        # Import libraries
        from astropy.io import fits
        from scipy.ndimage.interpolation import zoom
        # Set some variables
        detector_index = 101
        i_quantity = 0
        zoom = 0.8
        # Set some lists
        model_list = [
            'no_circumstellar_disks',
            'only_Aa',
            'only_Ab1',
            'only_Ab2',
            'only_Aa_Ab1',
            'only_Aa_Ab2',
            'only_Ab1_Ab2',
            'default'
        ]
        model_descr = [
            'No circumstellar disks',
            'Disk around Aa',
            'Disk around Ab1',
            'Disk around Ab2',
            'Disks around Aa and Ab1',
            'Disks around Aa and Ab2',
            'Disks around Ab1 and Ab2',
            'Disks around all stars'
        ]
        measurement_position_list = [
            [0, -1.05],
            [0, 0.85],
            [-0.57, -0.975],
            [1.17, -0.42],
            # [-1.3, -1.1],
            # [-0.57, 0.82]
        ]
        # Set beam size (in arcsec)
        self.file_io.beam_size = 0.05
        # Take colorbar label from quantity id
        cbar_label = self.file_io.get_quantity_labels(i_quantity)
        # Define output pdf
        self.file_io.init_plot_output(
            '2x2x2_I_emission_map', path=self.file_io.path['model'])
        # Sum up the flux inside the circle
        flux_sum = np.zeros((2, 4, len(measurement_position_list)))
        # Create two 2x3 plots
        for i_plot in range(2):
            # Init plot class variable
            plot = None
            for i_subplot in range(4):
                # Set paths of each simulation
                self.file_io.set_path_from_str(
                    'plot', 'gg_tau_disk', model_list[i_subplot + i_plot * 4], 'dust')
                # Create pdf file if show_plot is not chosen and read map data from file
                plot_data, header, plot_data_type = self.file_io.read_emission_map(
                    'polaris_detector_nr' + str(detector_index).zfill(4))
                # Take data for current quantity
                tbldata = plot_data[i_quantity, 0, :, :]
                # Create Matplotlib figure
                if i_subplot == 0:
                    plot = Plot(self.model, self.parse_args, ax_unit='arcsec',
                                nr_x_images=2, nr_y_images=2)
                # Plot imshow
                plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=i_subplot, set_bad_to_min=True,
                                 norm='LogNorm', vmin=1e-7, vmax=1e-4)
                # Plot map description
                plot.plot_text(text_pos=[0.03, 0.97], relative_position=True,
                               text=r'$\text{' +
                               model_descr[i_subplot +
                                           i_plot * 4] + r'}$',
                               horizontalalignment='left', verticalalignment='top',
                               ax_index=i_subplot, color='white')
                for i_pos, center_pos in enumerate(measurement_position_list):
                    plot.plot_text(text_pos=center_pos,
                                   text=str(i_pos + 1), ax_index=i_subplot, color='black', fontsize=10,
                                   bbox=dict(boxstyle='circle, pad=0.2', facecolor='white', alpha=0.5))
                # Get the image sidelength
                sidelength_x = 2. * self.model.tmp_parameter['radius_x_arcsec']
                sidelength_y = 2. * self.model.tmp_parameter['radius_y_arcsec']
                # Get number of pixel per axis
                nr_pixel_x = tbldata.shape[-2]
                nr_pixel_y = tbldata.shape[-1]
                # Find the considered pixel
                for i_x in range(nr_pixel_x):
                    for i_y in range(nr_pixel_y):
                        pos = np.subtract(np.multiply(np.divide(np.add([i_x, i_y], 0.5),
                                                                [nr_pixel_x, nr_pixel_y]), [sidelength_x, sidelength_y]),
                                          np.divide([sidelength_x, sidelength_y], 2.))
                        for i_pos, center_pos in enumerate(measurement_position_list):
                            pos_r = np.linalg.norm(
                                np.subtract(pos, center_pos))
                            if pos_r < 0.2:
                                flux_sum[i_plot, i_subplot,
                                         i_pos] += tbldata[i_x, i_y]
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)
        print(r'\begin{tabular}{lcccc}')
        print(r'\theadstart')
        print(
            r'\thead Configuration & \thead $\boldsymbol{\textit{F}}$ (1) & \thead $\boldsymbol{\textit{F}}$ (2) & \thead $\boldsymbol{\textit{F}}$ (3) & \thead $\boldsymbol{\textit{F}}$ (4) \\')
        print(r'\tbody')
        for i_plot in range(2):
            for i_subplot in range(4): # 6
                #if i_subplot > 1:
                print(model_descr[i_subplot - 2 + i_plot * 4], '&', ' & '.join(
                    '$\SI{' + f'{x:1.2e}' + '}{}$' for x in flux_sum[i_plot, i_subplot, :]), r'\\')
        #print('Observation &', ' & '.join(
        #    '$\SI{' + f'{x:1.2e}' + '}{}$' for x in flux_sum[0, 0, :]), r'\\')
        print(r'\tend')
        print(r'\end{tabular}')
        print(r'\begin{tabular}{lcccccc}')
        print(r'\theadstart')
        print(
            r'\thead Configuration & \thead $\boldsymbol{\delta \textit{F}}$ (region 1) & \thead $\boldsymbol{\delta \textit{F}}$ (region 2) & \thead $\boldsymbol{\delta \textit{F}}$ (region 3) & \thead $\boldsymbol{\delta \textit{F}}$ (region 4) & \thead \textit{F} (region 1) / \textit{F} (region 2) \\')
        print(r'\tbody')
        for i_plot in range(2):
            for i_subplot in range(4): # 6
                #if i_subplot > 1:
                print(model_descr[i_subplot + i_plot * 4], '&', ' & '.join('$\SI{-' + f'{x:1.0f}' + '}{\percent}$' for x in np.subtract(
                    100, 1e2*np.divide(flux_sum[i_plot, i_subplot, :], flux_sum[0, 0, :]))), '&', '& '.join('$\SI{' + f'{x:1.0f}' + '}{\percent}$' for x in [1e2*flux_sum[i_plot, i_subplot, 0] / flux_sum[i_plot, i_subplot, 1]]), r'\\')
        #print('Observation & x & x & x & x &', ' & '.join(
        #    '$\SI{' + f'{x:1.0f}' + '}{\percent}$' for x in [1e2*flux_sum[0, 0, 0]/flux_sum[0, 0, 1]]), r'\\')
        print(r'\tend')
        print(r'\end{tabular}')


    def plot_1006005(self):
        """Plot GG Tau A emission maps for different inclination axis.
        """
        # Import libraries
        from astropy.io import fits
        # Set some variables
        detector_index = 1
        i_quantity = 0
        # Set beam size (in arcsec)
        self.file_io.beam_size = 0.05
        # Take colorbar label from quantity id
        cbar_label = self.file_io.get_quantity_labels(i_quantity)
        # Define output pdf
        self.file_io.init_plot_output(
            'position_angle_analysis', path=self.file_io.path['model'])
        # Create a plot for each position angle
        for PA in range(-20, 75, 5):
            # Set paths of each simulation
            self.file_io.set_path_from_str(
                'plot', 'gg_tau_disk', 'position_angle_' + str(PA), 'dust')
            # Create pdf file if show_plot is not chosen and read map data from file
            plot_data, header, plot_data_type = self.file_io.read_emission_map(
                'polaris_detector_nr' + str(detector_index).zfill(4))
            # Create Matplotlib figure
            plot = Plot(self.model, self.parse_args, ax_unit='arcsec')
            # Take data for current quantity
            tbldata = plot_data[i_quantity, 0, :, :]
            # Plot imshow
            plot.plot_imshow(tbldata, cbar_label=cbar_label, set_bad_to_min=True,
                             norm='LogNorm', vmin=1e-6, vmax=1e-1, extend='neither')
            # Plot map description
            plot.plot_text(text_pos=[0.03, 0.97], relative_position=True,
                           text=r'$\text{PA}=\SI{' + str(PA) + '}{\degree}$',
                           horizontalalignment='left', verticalalignment='top', color='white')
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)

    def plot_1006006(self):
        """Plot GG Tau A azimuthal brightness distribution.
        """
        # Import libraries
        from astropy.io import fits
        # Set some variables
        detector_index = 100
        i_quantity = 0
        # Take colorbar label from quantity id
        cbar_label = self.file_io.get_quantity_labels(i_quantity)
        # Define output pdf
        self.file_io.init_plot_output('azimuthal_brightness_distribution')
        # Create pdf file if show_plot is not chosen and read map data from file
        plot_data, header, plot_data_type = self.file_io.read_emission_map(
            'polaris_detector_nr' + str(detector_index).zfill(4))
        # Create radially averaged azimuthal profile
        inclination = -37.0 / 180. * np.pi
        inc_PA = (360. - 270. + 7.) / 180. * np.pi
        inc_offset = -0.2
        R_min = 170.  # AU
        R_max = 260.  # AU
        angles = np.linspace(0, 2*np.pi, 1000)
        data = np.zeros((len(angles), 2))
        for i in range(len(angles)):
            data[i, :] = self.math.apply_inclination(
                pos=[np.cos(angles[i]), np.sin(angles[i])],
                inclination=inclination, inc_PA=inc_PA, inc_offset=inc_offset)
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args,
                    limits=[-R_max, R_max, -R_max, R_max])
        # Take data for current quantity
        tbldata = plot_data[i_quantity, 0, :, :]
        # Plot imshow
        plot.plot_imshow(tbldata, cbar_label=cbar_label, set_bad_to_min=True,
                         norm='LogNorm', vmin=5e-4, vmax=1e-2, extend='neither')
        plot.plot_line(R_min * data[:, 0], R_min *
                       data[:, 1], no_grid=True, color='white')
        plot.plot_line(R_max * data[:, 0], R_max *
                       data[:, 1], no_grid=True, color='white')
        plot.plot_line([R_min, R_max], [0, 0], no_grid=True, color='white')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)
        # Set azimuthal parameters
        azimuthal_parameter = [0, 0, inclination, inc_PA,
                               R_min * self.math.const['au'], R_max * self.math.const['au'], inc_offset]
        # Get radially averaged azimuthal brightness profile
        position, data = self.file_io.create_azimuthal_profile(
            plot_data[i_quantity, 0, ...], azimuthal_parameter, subpixel=4, N_ph=90)
        # Angle to degree
        position *= 180. / np.pi
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, with_cbar=False,
                    xlabel=r'$\theta\ [\si{\degree}]$',
                    ylabel=cbar_label, limits=[position[0], position[-1], None, None])
        # Plot cut/radial profile
        plot.plot_line(position, data, log='y')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)
