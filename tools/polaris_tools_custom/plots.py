#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import curve_fit

from polaris_tools_modules.gas import GasChooser
from polaris_tools_modules.visual import Plot


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
        from polaris_tools_modules.math import Math
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
                    xlabel=[r'$\mathsf{Frequency}$', r'$\mathsf{Frequency}$'],
                    ylabel=[r'$\mathsf{Flux}$', r'$\mathsf{Flux}$'], nr_x_images=2, nr_y_images=2,
                    extent=[-2.5, 2.5, None, None])

        for i_mag_field in range(2):
            # Set paths of each simulation
            self.file_io.set_path_from_str(
                'plot', 'cube', simulation_names[i_mag_field], 'line')
            # Read spectrum data
            plot_data, header = self.file_io.read_spectrum('line_spectrum_species_0001_line_0001')
            velocity = []
            for vch in range(header['nr_channels']):
                # Get velocity of current channel
                velocity.append(1e-3 * self.math.get_velocity(vch, header['nr_channels'], header['max_velocity']))

            for i_quantity in range(0, 4, 3):
                tmp_xdata = velocity
                tmp_ydata = plot_data[i_quantity, :]
                tmp_ydata /= np.max(tmp_ydata)
                if i_quantity == 0:
                    # Plot spectrum as line
                    tmp_ydata *= 1.5
                    plot.plot_line(ax_index=i_mag_field, xdata=tmp_xdata, ydata=tmp_ydata, color='blue', label=r'$I$')
                elif i_quantity == 3:
                    # Plot spectrum as line
                    plot.plot_line(ax_index=2 + i_mag_field, xdata=tmp_xdata, ydata=tmp_ydata, color='#ff8000', label=r'$V$')
                channel_width = tmp_xdata[1] - tmp_xdata[0]
                tmp_ydata = np.gradient(plot_data[0, :], channel_width)

                # Define fit function to get proportionality factor
                def fit_function(vch, prop_factor):
                    return prop_factor * tmp_ydata[vch]

                # Calculate best-fit parameter
                popt, pcov = curve_fit(fit_function, range(header['nr_channels']), plot_data[3, :])
                tmp_ydata *= popt[0]
                tmp_ydata /= np.max(tmp_ydata)
                tmp_ydata *= 0.2
                plot.plot_line(ax_index=int(i_quantity / 3) * 2 + i_mag_field, xdata=tmp_xdata, ydata=tmp_ydata,
                               linestyle='--', color='black', label=r'$\mathsf{d}I/\mathsf{d}\nu$')
                # Plot proportional arrows to visualise the relation between dI/dv and V
                if i_quantity == 3:
                    if i_mag_field == 0:
                        plot.plot_annotate(ax_index=int(i_quantity / 3) * 2 + i_mag_field,
                                           text=r'${\propto}B_\mathsf{LOS}$',
                                           text_pos=[-0.71, 0.30], pos=[-0.71, 0.85], color='green')
                        plot.plot_annotate(ax_index=int(i_quantity / 3) * 2 + i_mag_field,
                                           text='', text_pos=[0.71, -0.26], pos=[0.71, -0.8],
                                           color='green')
                    elif i_mag_field == 1:
                        plot.plot_annotate(ax_index=int(i_quantity / 3) * 2 + i_mag_field,
                                           text=r'${\not\propto}B_\mathsf{LOS}$',
                                           text_pos=[-0.81, 0.30], pos=[-0.81, 0.85], color='red')
                        plot.plot_annotate(ax_index=int(i_quantity / 3) * 2 + i_mag_field,
                                           text='', text_pos=[0.81, -0.25], pos=[0.81, -0.8],
                                           color='red')
                # Plot the legend
                if int(i_quantity / 3) * 2 + i_mag_field in [0, 1]:
                    plot.plot_legend(ax_index=int(i_quantity / 3)
                                     * 2 + i_mag_field, loc='upper right')
                else:
                    plot.plot_legend(ax_index=int(
                        i_quantity / 3) * 2 + i_mag_field)
        plot.ax_list[0].set_ylim([-0.25, 1.55])
        plot.ax_list[1].set_ylim([-0.25, 1.55])
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
                # Load radiation source to get position of the binary stars
                from polaris_tools_modules.source import SourceChooser
                radiation_source_chooser = SourceChooser(
                    self.model, self.file_io, self.parse_args)
                radiation_source = radiation_source_chooser.get_module_from_name(
                    'gg_tau_binary')
                # Plot position of binary stars (with conversion from m to au)
                radiation_source.tmp_parameter['position_star'][0] = np.divide(radiation_source.tmp_parameter['position_star'][0],
                                                                               self.math.const['au'])
                plot.plot_text(
                    radiation_source.tmp_parameter['position_star'][0], r'$+$')
                radiation_source.tmp_parameter['position_star'][1] = np.divide(radiation_source.tmp_parameter['position_star'][1],
                                                                               self.math.const['au'])
                plot.plot_text(
                    radiation_source.tmp_parameter['position_star'][1], r'$+$')
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
        plot.create_zoom_axis(model=self.model, zoom_factor=3.2)
        plot.ax_list[1].spines['bottom'].set_color('0.5')
        plot.ax_list[1].spines['top'].set_color('0.5')
        plot.ax_list[1].spines['right'].set_color('0.5')
        plot.ax_list[1].spines['left'].set_color('0.5')
        # Plot zoomed image
        self.basic_plots.plot_midplane_map_base(visualization_input, plot, tbldata_zoom, vec_field_data_zoom,
                                                ax_index=1, plot_cbar=False, vmin=1e-14, set_bad_to_min=True, cmap='inferno')
        # Load radiation source to get position of the binary stars
        from polaris_tools_modules.source import SourceChooser
        radiation_source_chooser = SourceChooser(
            self.model, self.file_io, self.parse_args)
        radiation_source = radiation_source_chooser.get_module_from_name(
            'gg_tau_stars')
        # Plot position of binary stars (with conversion from m to au)
        star_descr = ['Aa', 'Ab1', 'Ab2']
        for i_star, position in enumerate(radiation_source.tmp_parameter['position_star']):
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
        from scipy.interpolate import interp1d
        data = np.genfromtxt(
            self.file_io.path['model'] + 'foreground_extinction.dat')
        ext_cross_section = interp1d(data[:, 0] * 1e-6, data[:, 1])

        def reddening(wl, value, a_v):
            C_ext_wl = ext_cross_section(wl)
            A_wl = C_ext_wl / 8.743994626531222e-17 * a_v
            tau_wl = A_wl / 1.086
            return value * np.exp(-tau_wl)

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
        for i, wl in enumerate(data_wl):
            if wl > 1e-6 and data_flux[i] < 1e-2:
                data_flux[i] = np.nan
                data_flux_error[i] = np.nan
                data_wl[i] = np.nan
        # Plot total spectral energy distribution
        plot.plot_line(data_wl, data_flux, yerr=data_flux_error, log='xy',
                       linestyle='none', marker='.', color='black',
                       alpha=0.7, label=r'$\text{Photometric observations}$')
        # Set some variables
        detector_index = 204
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
        line_cycle = [':', '-.', '-.', '-.', '--', '--', '--', '-']
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
                quantity = reddening(wavelengths, sed_data[i_quantity, :], 0.3)
                # Plot spectral energy distribution
                plot.plot_line(wavelengths, quantity, log='xy', linestyle=line_cycle[i_model],
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
        detector_index = 103
        i_quantity = 4
        calc = False
        observation = 'SCUBA'  # 'SCUBA', 'SPHERE'
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
            [0, -1.60],
            [0, 0.85],
            [-0.57, -0.975],
            [1.17, -0.42],
            [1.5, 0],
        ]
        # Set beam size (in arcsec)
        if observation == 'SCUBA':
            self.file_io.beam_size = 0.07
            vmin = 2e-7
            vmax = 8e-5
        elif observation == 'SPHERE':
            self.file_io.beam_size = 0.03
            vmin = 3e-8
            vmax = 1e-5
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
                    if observation == 'SCUBA':
                        hdulist = fits.open(
                            '/home/rbrauer/Documents/projects/005_gg_tau/near_infrared_imaging_paper/sub_pi.fits')
                        tbldata = cropND(
                            hdulist[0].data.T, (500, 500), offset=[0, 30]) / 1e7
                    elif observation == 'SPHERE':
                        hdulist = fits.open(
                            '/home/rbrauer/Documents/projects/005_gg_tau/SPHERE_observation_miriam/GG_Tau_2016-11-191_I_POL.fits')
                        tbldata = cropND(
                            hdulist[0].data.T, (390, 390), offset=[6, 14]) / 8e6
                else:
                    continue
                plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=i_subplot, set_bad_to_min=True,
                                 norm='LogNorm', vmin=vmin, vmax=vmax, cmap='magma')
                if i_plot == 1:
                    x=np.linspace(0, 2*np.pi, 100)
                    inc_pa = 7./180.*np.pi
                    pos_x = []
                    pos_y = []
                    for angle in x:
                        pos_x.append(0.05 + 0.5*2.65*np.cos(angle) * np.cos(inc_pa) - 0.5*2.1 * np.sin(angle) * np.sin(inc_pa))
                        pos_y.append(-0.3 + 0.5*2.65*np.cos(angle) * np.sin(inc_pa) + 0.5*2.1 * np.sin(angle) * np.cos(inc_pa))
                    plot.plot_line(pos_x, pos_y, ax_index=0, color='cyan', linestyle=':', no_grid=True)
                    plot.plot_line(pos_x, pos_y, ax_index=2, color='cyan', linestyle=':', no_grid=True)
                    #plot.plot_ellipse([0.05, -0.3], 2.65, 2.1, angle=7, facecolor='none',
                    #             ax_index=0, edgecolor='cyan', linestyle=':')
                    #plot.plot_ellipse([0.05, -0.3], 2.65, 2.1, angle=7, facecolor='none',
                    #             ax_index=2, edgecolor='cyan', linestyle=':')
                    x=np.linspace(np.pi*1.17, 1.97*np.pi, 100)
                    inc_pa = 0 # 7./180.*np.pi
                    pos_x = []
                    pos_y = []
                    for angle in x:
                        pos_x.append(-0.1 + 0.5*2.65*np.cos(angle) * np.cos(inc_pa) + 0.5*2.1 * np.sin(angle) * np.sin(inc_pa))
                        pos_y.append(0.3 + 0.5*2.65*np.cos(angle) * np.sin(inc_pa) + 0.5*2.1 * np.sin(angle) * np.cos(inc_pa))
                    plot.plot_line(pos_x, pos_y, ax_index=0, color='orange', linestyle=':', no_grid=True)
                    plot.plot_line(pos_x, pos_y, ax_index=2, color='orange', linestyle=':', no_grid=True)
                for i_pos, center_pos in enumerate(measurement_position_list):
                    plot.plot_text(text_pos=center_pos,
                                   text=str(i_pos + 1), ax_index=i_subplot, color='white', fontsize=10,
                                   bbox=dict(boxstyle='circle, pad=0.2', facecolor='none', edgecolor='white', alpha=0.5))
                if calc:
                    # Get the image sidelength
                    sidelength_x = 2. * \
                        self.model.tmp_parameter['radius_x_arcsec']
                    sidelength_y = 2. * \
                        self.model.tmp_parameter['radius_y_arcsec']
                    # Get number of pixel per axis
                    nr_pixel_x = tbldata.shape[-2]
                    nr_pixel_y = tbldata.shape[-1]
                    # Find the considered pixel
                    count = np.zeros(len(measurement_position_list))
                    obs_error = []
                    for i_x in range(nr_pixel_x):
                        for i_y in range(nr_pixel_y):
                            pos = np.subtract(np.multiply(np.divide(np.add([i_x, i_y], 0.5),
                                                                    [nr_pixel_x, nr_pixel_y]), [sidelength_x, sidelength_y]),
                                              np.divide([sidelength_x, sidelength_y], 2.))
                            for i_pos, center_pos in enumerate(measurement_position_list):
                                pos_r = np.linalg.norm(
                                    np.subtract(pos, center_pos))
                                if pos_r < 0.15:
                                    flux_sum[i_plot, i_subplot,
                                             i_pos] += tbldata[i_x, i_y]
                                    count[i_pos] += 1
                            pos_r = np.linalg.norm(
                                np.subtract(pos, [-1.8, -1.8]))
                            if pos_r < 0.15:
                                obs_error.append(tbldata[i_x, i_y])
                    if i_plot == 0 and i_subplot == 0:
                        flux = flux_sum[i_plot, i_subplot, :] / count
                        da = np.std(obs_error)
                        uncertainty = np.sqrt(
                            (1/flux[1:]*da)**2 + (flux[0]/flux[1:]**2*da)**2)
                        uncertainty[-1] = np.sqrt(
                            (1/flux[1]*da)**2 + (flux[-1]/flux[1]**2*da)**2)
                        bg_flux = np.divide(flux[:], np.std(obs_error))
                        print('obs. error:', uncertainty)
                        #print('BG flux:', bg_flux)
                # Hide second observation plot
                plot.remove_axes(ax_index=1)
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)
        if calc:
            print(r'\begin{tabular}{lcccccc}')
            print(
                r'Configuration & $\textit{PI}_1 [\si{Jy}]$ & $\textit{PI}_2 [\si{Jy}]$ & $\textit{PI}_3 [\si{Jy}]$ & $\textit{PI}_4 [\si{Jy}]$ & $\textit{PI}_5 [\si{Jy}]$  & $\textit{PI}_6 [\si{Jy}]$\\')
            print(r'\hline')
            for i_plot in range(2):
                for i_subplot in range(6):
                    if i_subplot > 1:
                        print(model_descr[i_subplot - 2 + i_plot * 4], '&', ' & '.join(
                            '$\SI{' + f'{x:1.2e}' + '}{}$' for x in flux_sum[i_plot, i_subplot, :]), r'\\')
            print(r'\end{tabular}')
            print(r'\begin{tabular}{lcccccc}')
            print(
                r'Configuration & $\delta \textit{PI}_1 [\si{\percent}]$ & $\delta \textit{PI}_2 [\si{\percent}]$ & $\delta \textit{PI}_3 [\si{\percent}]$ & $\delta \textit{PI}_4 [\si{\percent}]$ & $\delta \textit{PI}_5 [\si{\percent}]$ & $\delta \textit{PI}_6 [\si{\percent}]$ \\')
            print(r'\hline')
            for i_plot in range(2):
                for i_subplot in range(6):
                    if i_subplot > 1:
                        print(model_descr[i_subplot - 2 + i_plot * 4], '&', ' & '.join('$\SI{-' + f'{x:1.0f}' + '}{}$' for x in np.absolute(np.subtract(
                            100, 1e2*np.divide(flux_sum[i_plot, i_subplot, :], flux_sum[0, 2, :])))), r'\\')
            print(r'\end{tabular}')
            print(r'\begin{tabular}{lccccccc}')
            print(
                r'Configuration & $\textit{PI}_1 / \textit{PI}_2$ & $\textit{PI}_1 / \textit{PI}_3$ & $\textit{PI}_1 / \textit{PI}_4$ & $\textit{PI}_1 / \textit{PI}_5$ & $\textit{PI}_6 / \textit{PI}_3$ \\')
            print(r'\hline')
            for i_plot in range(2):
                for i_subplot in range(6):
                    if i_subplot > 1:
                        string = model_descr[i_subplot - 2 + i_plot * 4]
                        for i_x, x in enumerate(np.divide(flux_sum[i_plot, i_subplot, 0], flux_sum[i_plot, i_subplot, 1:-1])):
                            if abs(x - np.divide(flux_sum[i_plot, 0, 0], flux_sum[i_plot, 0, 1 + i_x])) <= uncertainty[i_x]:
                                string += r' & {\color{new_green}\textbf{' + \
                                    f'{x:1.2f}' + '}}'
                            else:
                                string += r' & {\color{red}' + \
                                    f'{x:1.2f}' + '}'
                        x = np.divide(
                            flux_sum[i_plot, i_subplot, -1], flux_sum[i_plot, i_subplot, 1])
                        if abs(x - np.divide(flux_sum[i_plot, 0, -1], flux_sum[i_plot, 0, 2])) <= uncertainty[i_x]:
                            string += r' & {\color{new_green}\textbf{' + \
                                f'{x:1.2f}' + '}}'
                        else:
                            string += r' & {\color{red}' + f'{x:1.2f}' + '}'
                        print(string, r'\\')
            print(r'\hline')
            if observation == 'SCUBA':
                string = 'Observation (Subaru/HiCIAO)'
            elif observation == 'SPHERE':
                string = 'Observation (SPHERE/IRDIS)'
            for i_x, x in enumerate(np.divide(flux_sum[i_plot, 0, 0], flux_sum[i_plot, 0, 1:-1])):
                string += ' & $' + f'{x:1.2f}' + '\pm' + \
                    f'{uncertainty[i_x]:1.2f}' + '$'
            x = np.divide(flux_sum[i_plot, 0, -1], flux_sum[i_plot, 0, 1])
            string += ' & $' + f'{x:1.2f}' + '\pm' + \
                f'{uncertainty[-1]:1.2f}' + '$'
            print(string, r'\\')
            print(r'\\')
            print(
                r'Configuration & $\textit{PI}_1 / \sigma_\text{BG}$    & $\textit{PI}_2 / \sigma_\text{BG}$    & $\textit{PI}_3 / \sigma_\text{BG}$ & $\textit{PI}_4 / \sigma_\text{BG}$    & $\textit{PI}_5 / \sigma_\text{BG}$ & $\textit{PI}_6 / \sigma_\text{BG}$ \\')
            if observation == 'SCUBA':
                string = 'Observation (Subaru/HiCIAO)'
            elif observation == 'SPHERE':
                string = 'Observation (SPHERE/IRDIS)'
            for i_x, x in enumerate(bg_flux):
                string += ' & $' + f'{x:1.1f}' + '$'
            print(string, r'\\')
            print(r'\end{tabular}')

    def plot_1006004(self):
        """Plot two times 4 intensity emission maps with differnet configurations for GG Tau A.
        """
        # Import libraries
        from astropy.io import fits
        from scipy.ndimage.interpolation import zoom
        # Set some variables
        detector_index = 103
        i_quantity = 0
        zoom = 0.8
        calc = False
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
            [0, -1.60],
            [0, 0.85],
            [-0.57, -0.975],
            [1.17, -0.42],
            [1.5, 0],
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
                import matplotlib
                cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                    "bbw", ["black", "royalblue", "white"])
                matplotlib.pyplot.register_cmap(cmap=cmap)
                plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=i_subplot, set_bad_to_min=True,
                                 norm='LogNorm', vmin=1e-7, vmax=5e-5, cmap='bbw')
                # Plot map description
                plot.plot_text(text_pos=[0.03, 0.97], relative_position=True,
                               text=r'$\text{' +
                               model_descr[i_subplot +
                                           i_plot * 4] + r'}$',
                               horizontalalignment='left', verticalalignment='top',
                               ax_index=i_subplot, color='white')
                for i_pos, center_pos in enumerate(measurement_position_list):
                    plot.plot_text(text_pos=center_pos,
                                   text=str(i_pos + 1), ax_index=i_subplot, color='white', fontsize=10,
                                   bbox=dict(boxstyle='circle, pad=0.2', facecolor='none', edgecolor='white', alpha=0.5))
                if calc:
                    # Get the image sidelength
                    sidelength_x = 2. * \
                        self.model.tmp_parameter['radius_x_arcsec']
                    sidelength_y = 2. * \
                        self.model.tmp_parameter['radius_y_arcsec']
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
                                if pos_r < 0.15:
                                    flux_sum[i_plot, i_subplot,
                                             i_pos] += tbldata[i_x, i_y]
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)
        if calc:
            print(r'\begin{tabular}{lcccccc}')
            print(
                r'Configuration & $F_1 [\si{Jy}]$ & $F_2 [\si{Jy}]$ & $F_3 [\si{Jy}]$ & $F_4 [\si{Jy}]$ & $F_5 [\si{Jy}]$ & $F_6 [\si{Jy}]$ \\')
            print(r'\hline')
            for i_plot in range(2):
                for i_subplot in range(4):
                    print(model_descr[i_subplot + i_plot * 4], '&', ' & '.join(
                        '$\SI{' + f'{x:1.2e}' + '}{}$' for x in flux_sum[i_plot, i_subplot, :]), r'\\')
            print(r'\end{tabular}')
            print(r'\begin{tabular}{lcccccc}')
            print(
                r'Configuration & $\delta F_1 [\si{\percent}]$ & $\delta F_2 [\si{\percent}]$ & $\delta F_3 [\si{\percent}]$ & $\delta F_4 [\si{\percent}]$ & $\delta F_5 [\si{\percent}]$ & $\delta F_6 [\si{\percent}]$ \\')
            print(r'\hline')
            for i_plot in range(2):
                for i_subplot in range(4):
                    print(model_descr[i_subplot + i_plot * 4], '&', ' & '.join('$\SI{-' + f'{x:1.0f}' + '}{}$' for x in np.absolute(np.subtract(
                        100, 1e2*np.divide(flux_sum[i_plot, i_subplot, :], flux_sum[0, 0, :])))), r'\\')
            print(r'\end{tabular}')
            print(r'\begin{tabular}{lccccc}')
            print(
                r'Configuration & $F_1 / F_2$ & $F_1 / F_3$ & $F_1 / F_4$ & $F_1 / F_5$ & $F_6 / F_3$ \\')
            print(r'\hline')
            for i_plot in range(2):
                for i_subplot in range(4):
                    print(model_descr[i_subplot + i_plot * 4], '&', ' & '.join(
                        f'{x:1.2f}' for x in np.divide(flux_sum[i_plot, i_subplot, 0], flux_sum[i_plot, i_subplot, 1:-1])), '&',
                        f'{np.divide(flux_sum[i_plot, i_subplot, -1], flux_sum[i_plot, i_subplot, 2]):1.2f}', r'\\')
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

    def plot_1006007(self):
        """Plot GG Tau A observation and one configuration
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
        detector_index = 101  # 104
        i_quantity = 0 #4
        i_subplot = 9
        #vmin = 1e-7
        #vmax = 1e-4
        vmin = 2e-7
        vmax = 8e-5
        observation = 'SCUBA'  # 'SCUBA', 'SPHERE'
        model_list = [
            'no_circumstellar_disks',
            'only_Aa',
            'only_Ab1',
            'only_Ab2',
            'only_Aa_Ab1',
            'only_Aa_Ab2',
            'only_Ab1_Ab2',
            'default',
            'vertical_Ab2',
            'vertical_Ab2_hres'
        ]
        model_descr = [
            'No circumstellar disks',
            'Disk around Aa',
            'Disk around Ab1',
            'Disk around Ab2',
            'Disks around Aa and Ab1',
            'Disks around Aa and Ab2',
            'Disks around Ab1 and Ab2',
            'Disks around all stars',
            'Vertical disk around Ab2',
            'Vertical disk around Ab2'
        ]
        measurement_position_list = [
            [0, -1.05],
            [0, -1.60],
            [0, 0.85],
            [-0.57, -0.975],
            [1.17, -0.42],
            [1.5, -0.15],
        ]
        # Set beam size (in arcsec)
        if observation == 'SCUBA':
            self.file_io.beam_size = 0.07
        elif observation == 'SPHERE':
            self.file_io.beam_size = 0.03
        # Take colorbar label from quantity id
        cbar_label = self.file_io.get_quantity_labels(i_quantity)
        # Define output pdf
        self.file_io.init_plot_output(
            'PI_emission_map_comparison_with_circles', path=self.file_io.path['model'])
        # Sum up the flux inside the circle
        flux_sum = np.zeros(len(measurement_position_list))
        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', 'gg_tau_disk', model_list[i_subplot], 'dust')
        # Create pdf file if show_plot is not chosen and read map data from file
        plot_data, header, plot_data_type = self.file_io.read_emission_map(
            'polaris_detector_nr' + str(detector_index).zfill(4))
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, ax_unit='arcsec',
                    nr_x_images=2, nr_y_images=1, size_x=6 / 1.4298620007401583)
        # Take data for current quantity
        tbldata = plot_data[i_quantity, 0, :, :]
        # plot imshow
        plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=1, set_bad_to_min=True,
                         norm='LogNorm', vmin=vmin, vmax=vmax, cmap='magma')
        # Load radiation source to get position of the binary stars
        from polaris_tools_modules.source import SourceChooser
        radiation_source_chooser = SourceChooser(
            self.model, self.file_io, self.parse_args)
        radiation_source = radiation_source_chooser.get_module_from_name(
            'gg_tau_stars')
        # Plot position of binary stars (with conversion from m to au)
        star_descr = ['Aa', 'Ab1', 'Ab2']
        for i_star, position in enumerate(radiation_source.tmp_parameter['position_star']):
            if i_star == 0:
                offset = np.array([-35, -10]) * 1.3
            elif i_star == 1:
                offset = np.array([-24, 10]) * 1.3
            elif i_star == 2:
                offset = np.array([44, -12]) * 1.2
            pos_arcsec = np.add(
                np.divide(position[0:2], self.math.const['au']), offset) / 140.
            plot.plot_text(pos_arcsec,
                           text=r'$\text{' + star_descr[i_star] + r'}$', ax_index=1, color='white')
        for i_pos, center_pos in enumerate(measurement_position_list):
            plot.plot_text(text_pos=center_pos, text=str(i_pos + 1), ax_index=1, color='white', fontsize=10, bbox=dict(boxstyle='circle, pad=0.2', facecolor='none', edgecolor='white', alpha=0.5))
        # Plot map description
        plot.plot_text(text_pos=[0.03, 0.97], relative_position=True, text=r'$\text{' +
                       model_descr[i_subplot] + r'}$', horizontalalignment='left', verticalalignment='top',
                       ax_index=1, color='white')
        # Get the image sidelength
        sidelength_x = 2. * \
            self.model.tmp_parameter['radius_x_arcsec']
        sidelength_y = 2. * \
            self.model.tmp_parameter['radius_y_arcsec']
        # Get number of pixel per axis
        nr_pixel_x = tbldata.shape[-2]
        nr_pixel_y = tbldata.shape[-1]
        # Find the considered pixel
        count = np.zeros(len(measurement_position_list))
        obs_error = []
        for i_x in range(nr_pixel_x):
            for i_y in range(nr_pixel_y):
                pos = np.subtract(np.multiply(np.divide(np.add([i_x, i_y], 0.5),
                                                        [nr_pixel_x, nr_pixel_y]), [sidelength_x, sidelength_y]),
                                  np.divide([sidelength_x, sidelength_y], 2.))
                for i_pos, center_pos in enumerate(measurement_position_list):
                    pos_r = np.linalg.norm(
                        np.subtract(pos, center_pos))
                    if pos_r < 0.15:
                        flux_sum[i_pos] += tbldata[i_x, i_y]
        # print(flux_sum[:])
        print(r'\begin{tabular}{lccccc}')
        print(
            r'Configuration & $\textit{PI}_1 [\si{Jy}]$ & $\textit{PI}_2 [\si{Jy}]$ & $\textit{PI}_3 [\si{Jy}]$ & $\textit{PI}_4 [\si{Jy}]$ & $\textit{PI}_5 [\si{Jy}]$  & $\textit{PI}_6 [\si{Jy}]$\\')
        print(r'\hline')
        print(model_descr[i_subplot], '&', ' & '.join(
            '$\SI{' + f'{x:1.2e}' + '}{}$' for x in flux_sum[:]), r'\\')
        print(r'\begin{tabular}{lccccc}')
        print(
            r'Configuration & $\textit{PI}_1 / \textit{PI}_2$ & $\textit{PI}_1 / \textit{PI}_3$ & $\textit{PI}_1 / \textit{PI}_4$ & $\textit{PI}_1 / \textit{PI}_5$ & $\textit{PI}_3 / \textit{PI}_6$ \\')
        print(r'\hline')
        string = model_descr[i_subplot]
        for i_x, x in enumerate(np.divide(flux_sum[0], flux_sum[1:-1])):
            string += r' & {\color{red}' + f'{x:1.2f}' + '}'
        string += r'& - & {\color{red}' + \
            f'{np.divide(flux_sum[-1], flux_sum[2]):1.2f}' + '}'
        print(string, r'\\')
        print(r'\end{tabular}')
        # Observation
        if observation == 'SCUBA':
            hdulist = fits.open(
                '/home/rbrauer/Documents/projects/005_gg_tau/near_infrared_imaging_paper/sub_pi.fits')
            tbldata = cropND(
                hdulist[0].data.T, (500, 500), offset=[0, 30]) / 2e7  # 1e7
        elif observation == 'SPHERE':
            hdulist = fits.open(
                '/home/rbrauer/Documents/projects/005_gg_tau/SPHERE_observation_miriam/GG_Tau_2016-11-191_I_POL.fits')
            tbldata = cropND(
                hdulist[0].data.T, (390, 390), offset=[6, 14]) / 8e6
        # plot imshow
        plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=0, set_bad_to_min=True,
                         norm='LogNorm', vmin=vmin, vmax=vmax, cmap='magma')
        for i_pos, center_pos in enumerate(measurement_position_list):
            if i_pos == 5:
                center_pos[1] += 0.15
            plot.plot_text(text_pos=center_pos, text=str(i_pos + 1), ax_index=0, color='white',
                           fontsize=10, bbox=dict(boxstyle='circle, pad=0.2', facecolor='none', edgecolor='white', alpha=0.5))
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_1006008(self):
        """Plot GG Tau A observation and one configuration
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
        detector_index = 101  # 104
        i_quantity = 0
        i_subplot = 9
        vmin = 1e-7
        vmax = 5e-5
        model_list = [
            'no_circumstellar_disks',
            'only_Aa',
            'only_Ab1',
            'only_Ab2',
            'only_Aa_Ab1',
            'only_Aa_Ab2',
            'only_Ab1_Ab2',
            'default',
            'vertical_Ab2',
            'vertical_Ab2_hres'
        ]
        model_descr = [
            'No circumstellar disks',
            'Disk around Aa',
            'Disk around Ab1',
            'Disk around Ab2',
            'Disks around Aa and Ab1',
            'Disks around Aa and Ab2',
            'Disks around Ab1 and Ab2',
            'Disks around all stars',
            'Vertical disk around Ab2',
            'Vertical disk around Ab2'
        ]
        measurement_position_list = [
            [0, -1.05],
            [0, -1.60],
            [0, 0.85],
            [-0.57, -0.975],
            [1.17, -0.42],
            [1.5, -0.15],
        ]
        # Set beam size (in arcsec)
        self.file_io.beam_size = 0.05
        # Take colorbar label from quantity id
        cbar_label = self.file_io.get_quantity_labels(i_quantity)
        # Define output pdf
        self.file_io.init_plot_output(
            'I_emission_map_comparison_with_circles', path=self.file_io.path['model'])
        # Sum up the flux inside the circle
        flux_sum = np.zeros(len(measurement_position_list))
        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', 'gg_tau_disk', model_list[i_subplot], 'dust')
        # Create pdf file if show_plot is not chosen and read map data from file
        plot_data, header, plot_data_type = self.file_io.read_emission_map(
            'polaris_detector_nr' + str(detector_index).zfill(4))
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, ax_unit='arcsec',
                    size_x=6 / 1.4298620007401583)
        # Take data for current quantity
        tbldata = plot_data[i_quantity, 0, :, :]
        # plot imshow
        import matplotlib
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "bbw", ["black", "royalblue", "white"])
        matplotlib.pyplot.register_cmap(cmap=cmap)
        plot.plot_imshow(tbldata, cbar_label=cbar_label, set_bad_to_min=True,
                         norm='LogNorm', cmap='bbw', vmin=vmin, vmax=vmax)
        # Load radiation source to get position of the binary stars
        from polaris_tools_modules.source import SourceChooser
        radiation_source_chooser = SourceChooser(
            self.model, self.file_io, self.parse_args)
        radiation_source = radiation_source_chooser.get_module_from_name(
            'gg_tau_stars')
        # Plot position of binary stars (with conversion from m to au)
        star_descr = ['Aa', 'Ab1', 'Ab2']
        for i_star, position in enumerate(radiation_source.tmp_parameter['position_star']):
            if i_star == 0:
                offset = np.array([-35, -10]) * 1.3
            elif i_star == 1:
                offset = np.array([-24, 10]) * 1.3
            elif i_star == 2:
                offset = np.array([44, -12]) * 1.2
            pos_arcsec = np.add(
                np.divide(position[0:2], self.math.const['au']), offset) / 140.
            plot.plot_text(pos_arcsec,
                           text=r'$\text{' + star_descr[i_star] + r'}$', color='white')
        for i_pos, center_pos in enumerate(measurement_position_list):
            plot.plot_text(text_pos=center_pos, text=str(i_pos + 1), color='white',
                           fontsize=10, bbox=dict(boxstyle='circle, pad=0.2', facecolor='none', edgecolor='white', alpha=0.5))
        # Plot map description
        plot.plot_text(text_pos=[0.03, 0.97], relative_position=True,
                       text=r'$\text{' + model_descr[i_subplot] + r'}$', horizontalalignment='left', verticalalignment='top', color='white')
        # Get the image sidelength
        sidelength_x = 2. * \
            self.model.tmp_parameter['radius_x_arcsec']
        sidelength_y = 2. * \
            self.model.tmp_parameter['radius_y_arcsec']
        # Get number of pixel per axis
        nr_pixel_x = tbldata.shape[-2]
        nr_pixel_y = tbldata.shape[-1]
        # Find the considered pixel
        count = np.zeros(len(measurement_position_list))
        obs_error = []
        for i_x in range(nr_pixel_x):
            for i_y in range(nr_pixel_y):
                pos = np.subtract(np.multiply(np.divide(np.add([i_x, i_y], 0.5),
                                                        [nr_pixel_x, nr_pixel_y]), [sidelength_x, sidelength_y]),
                                  np.divide([sidelength_x, sidelength_y], 2.))
                for i_pos, center_pos in enumerate(measurement_position_list):
                    pos_r = np.linalg.norm(
                        np.subtract(pos, center_pos))
                    if pos_r < 0.15:
                        flux_sum[i_pos] += tbldata[i_x, i_y]
        # print(flux_sum[:])
        print(r'\begin{tabular}{lccccc}')
        print(
            r'Configuration & $F_1 [\si{Jy}]$ & $F_2 [\si{Jy}]$ & $F_3 [\si{Jy}]$ & $F_4 [\si{Jy}]$ & $F_5 [\si{Jy}]$  & $F_6 [\si{Jy}]$\\')
        print(r'\hline')
        print(model_descr[i_subplot], '&', ' & '.join(
            '$\SI{' + f'{x:1.2e}' + '}{}$' for x in flux_sum[:]), r'\\')
        print(r'\end{tabular}')
        print(r'\begin{tabular}{lccccc}')
        print(
            r'Configuration & $F_1 / F_2$ & $F_1 / F_3$ & $F_1 / F_4$ & $F_1 / F_5$ & $F_3 / F_6$ \\')
        print(r'\hline')
        string = model_descr[i_subplot]
        for i_x, x in enumerate(np.divide(flux_sum[0], flux_sum[1:-1])):
            string += r' & ' + f'{x:1.2f}'
        string += r'& ' + f'{np.divide(flux_sum[-1], flux_sum[2]):1.2f}'
        print(string, r'\\')
        print(r'\end{tabular}')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_1006009(self):
        """Plot GG Tau A observation and two configurations
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
        detector_index = 1  # 104
        i_quantity = 4
        i_subplot = 0
        vmin = 1e-7
        vmax = 1e-4
        observation = 'SCUBA'  # 'SCUBA', 'SPHERE'
        model_list = [
            'no_circumstellar_disks',
            'only_Aa',
            'only_Ab1',
            'only_Ab2',
            'only_Aa_Ab1',
            'only_Aa_Ab2',
            'only_Ab1_Ab2',
            'default',
            'vertical_Ab2',
            'vertical_Ab2_scale_height_Aa_0_8_very_small'
        ]
        model_descr = [
            r'$\text{No circumstellar disks}$',
            r'$\text{Disk around Aa}$',
            r'$\text{Disk around Ab1}$',
            r'$\text{Disk around Ab2}$',
            r'$\text{Disks around Aa and Ab1}$',
            r'$\text{Disks around Aa and Ab2}$',
            r'$\text{Disks around Ab1 and Ab2}$',
            r'$\text{Disks around all stars}$',
            r'$\text{Vertical disk around Ab2}$',
            r'$\text{Vertical disk around Ab2}$' + '\n' +
            r'$\text{Coplanar disks around Aa+Ab1}$'
        ]
        measurement_position_list = [
            [0, -1.05],
            [0, -1.60],
            [0, 0.85],
            [-0.57, -0.975],
            [1.17, -0.42],
            [1.5, -0.15],
        ]
        # Set beam size (in arcsec)
        if observation == 'SCUBA':
            self.file_io.beam_size = 0.07
        elif observation == 'SPHERE':
            self.file_io.beam_size = 0.03
        # Take colorbar label from quantity id
        cbar_label = self.file_io.get_quantity_labels(i_quantity)
        # Define output pdf
        self.file_io.init_plot_output(
            'proposal_image', path=self.file_io.path['model'])
        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', 'gg_tau_disk', model_list[i_subplot], 'dust')
        # Create pdf file if show_plot is not chosen and read map data from file
        plot_data, header, plot_data_type = self.file_io.read_emission_map(
            'polaris_detector_nr' + str(detector_index).zfill(4))
        # Create Matplotlib figure
        plot = Plot(self.model, self.parse_args, ax_unit='arcsec',
                    nr_x_images=3, nr_y_images=1, size_x=6 / 1.4298620007401583)
        for i_plot in range(2):
            # Set paths of each simulation
            self.file_io.set_path_from_str(
                'plot', 'gg_tau_disk', model_list[i_subplot + i_plot * 9], 'dust')
            # Create pdf file if show_plot is not chosen and read map data from file
            plot_data, header, plot_data_type = self.file_io.read_emission_map(
                'polaris_detector_nr' + str(detector_index).zfill(4))
            # Take data for current quantity
            tbldata = plot_data[i_quantity, 0, :, :]
            # plot imshow
            plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=i_plot + 1, set_bad_to_min=True,
                             norm='LogNorm', vmin=vmin, vmax=vmax, cmap='magma')
            if i_plot == 1:
                # Load radiation source to get position of the binary stars
                from polaris_tools_modules.source import SourceChooser
                radiation_source_chooser = SourceChooser(
                    self.model, self.file_io, self.parse_args)
                radiation_source = radiation_source_chooser.get_module_from_name(
                    'gg_tau_stars')
                # Plot position of binary stars (with conversion from m to au)
                star_descr = ['Aa', 'Ab1', 'Ab2']
                for i_star, position in enumerate(radiation_source.tmp_parameter['position_star']):
                    if i_star == 0:
                        offset = np.array([-35, -10]) * 1.3
                    elif i_star == 1:
                        offset = np.array([-24, 10]) * 1.3
                    elif i_star == 2:
                        offset = np.array([44, -12]) * 1.2
                    pos_arcsec = np.add(
                        np.divide(position[0:2], self.math.const['au']), offset) / 140.
                    plot.plot_text(pos_arcsec,
                                   text=r'$\text{' + star_descr[i_star] + r'}$', ax_index=i_plot + 1, color='white')
            for i_pos, center_pos in enumerate(measurement_position_list):
                plot.plot_text(text_pos=center_pos, text=str(i_pos + 1), ax_index=i_plot + 1, color='white',
                               fontsize=10, bbox=dict(boxstyle='circle, pad=0.2', facecolor='none', edgecolor='white', alpha=0.5))
            # Plot map description
            plot.plot_text(text_pos=[0.03, 0.97], relative_position=True,
                           text=model_descr[i_subplot + 9 *
                                            i_plot], horizontalalignment='left',
                           verticalalignment='top', ax_index=i_plot + 1, color='white')
        # Observation
        if observation == 'SCUBA':
            hdulist = fits.open(
                '/home/rbrauer/Documents/projects/005_gg_tau/near_infrared_imaging_paper/sub_pi.fits')
            tbldata = cropND(
                hdulist[0].data.T, (500, 500), offset=[0, 30]) / 1e7
        elif observation == 'SPHERE':
            hdulist = fits.open(
                '/home/rbrauer/Documents/projects/005_gg_tau/SPHERE_observation_miriam/GG_Tau_2016-11-191_I_POL.fits')
            tbldata = cropND(
                hdulist[0].data.T, (390, 390), offset=[6, 14]) / 8e6
        # plot imshow
        plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=0, set_bad_to_min=True,
                         norm='LogNorm', vmin=vmin, vmax=vmax, cmap='magma')
        for i_pos, center_pos in enumerate(measurement_position_list):
            if i_pos == 5:
                center_pos[1] += 0.15
            plot.plot_text(text_pos=center_pos, text=str(i_pos + 1), ax_index=0, color='white',
                           fontsize=10, bbox=dict(boxstyle='circle, pad=0.2', facecolor='none', edgecolor='white', alpha=0.5))
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_3001001(self):
        """Plot GG Tau A observation and one configuration
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
        i_subplot = 0
        vmin = 2e-7
        vmax = 8e-5
        observation = 'SCUBA'  # 'SCUBA', 'SPHERE'
        model_list = [
            'no_circumstellar_disks',
            'only_Aa',
            'only_Ab1',
            'only_Ab2',
            'only_Aa_Ab1',
            'only_Aa_Ab2',
            'only_Ab1_Ab2',
            'default',
            'vertical_Ab2',
            'vertical_Ab2_hres',
        ]
        # Set beam size (in arcsec)
        if observation == 'SCUBA':
            self.file_io.beam_size = 0.07
        elif observation == 'SPHERE':
            self.file_io.beam_size = 0.03
        # Take colorbar label from quantity id
        cbar_label = self.file_io.get_quantity_labels(i_quantity)
        # Define output pdf
        self.file_io.init_plot_output(
            'PI_emission_map_comparison', path=self.file_io.path['model'])
        # Set paths of each simulation
        self.file_io.set_path_from_str(
            'plot', 'gg_tau_disk', model_list[i_subplot], 'dust')
        # Create pdf file if show_plot is not chosen and read map data from file
        plot_data, header, plot_data_type = self.file_io.read_emission_map(
            'polaris_detector_nr' + str(detector_index).zfill(4))
        # Create Matplotlib figure
        #plot = Plot(self.model, self.parse_args, ax_unit='arcsec',
        #            nr_x_images=1, nr_y_images=1)
        plot = Plot(self.model, self.parse_args, ax_unit='arcsec',
                    nr_x_images=2, nr_y_images=1)
        # Take data for current quantity
        tbldata = plot_data[i_quantity, 0, :, :]
        # plot imshow
        plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=1, set_bad_to_min=True,
                        norm='LogNorm', vmin=vmin, vmax=vmax, cmap='magma', extend='neither')
        #Load radiation source to get position of the binary stars
        from polaris_tools_modules.source import SourceChooser
        radiation_source_chooser = SourceChooser(
            self.model, self.file_io, self.parse_args)
        radiation_source = radiation_source_chooser.get_module_from_name(
            'gg_tau_stars')
        # Plot position of binary stars (with conversion from m to au)
        star_descr = ['Aa', 'Ab1', 'Ab2']
        for i_star, position in enumerate(radiation_source.tmp_parameter['position_star']):
            if i_star == 0:
                offset = np.array([-35, -10]) * 1.4
            elif i_star == 1:
                offset = np.array([-22, 10]) * 1.5
            elif i_star == 2:
                offset = np.array([40, -12]) * 1.3
            pos_arcsec = np.add(
                np.divide(position[0:2], self.math.const['au']), offset) / 140.
            # plot.plot_text(pos_arcsec,
            #               text=r'$\text{' + star_descr[i_star] + r'}$', ax_index=1, color='white')
        # Observation
        if observation == 'SCUBA':
            hdulist = fits.open(
                '/home/rbrauer/Documents/projects/005_gg_tau/near_infrared_imaging_paper/sub_pi.fits')
            tbldata = cropND(
                hdulist[0].data.T, (500, 500), offset=[0, 30]) / 2e7
        elif observation == 'SPHERE':
            hdulist = fits.open(
                '/home/rbrauer/Documents/projects/005_gg_tau/SPHERE_observation_miriam/GG_Tau_2016-11-191_I_POL.fits')
            tbldata = cropND(
                hdulist[0].data.T, (390, 390), offset=[6, 14]) / 8e6
        # plot imshow
        plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=0, set_bad_to_min=True,
                         norm='LogNorm', vmin=vmin, vmax=vmax, cmap='magma', extend='neither')
        # Save figure to pdf file or print it on screen
        plot.save_figure(self.file_io)

    def plot_3001002(self):
        """Plot GG Tau A observation and one configuration
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
        vmin = 2e-7
        vmax = 8e-5
        observation = 'SCUBA'  # 'SCUBA', 'SPHERE'
        model_list = [
            'no_circumstellar_disks',
            'only_Aa',
            'only_Ab1',
            'only_Ab2',
            'only_Aa_Ab1',
            'only_Aa_Ab2',
            'only_Ab1_Ab2',
            'default',
            'vertical_Ab2_hres',
        ]
        model_descr = [
            'No circumstellar disks',
            'Disk around Aa',
            'Disk around Ab1',
            'Disk around Ab2',
            'Disks around Aa and Ab1',
            'Disks around Aa and Ab2',
            'Disks around Ab1 and Ab2',
            'Disks around all stars',
            'Vertical disk around Ab2',
        ]
        measurement_position_list = [
            [0, -1.05],
            [0, -1.60],
            [0, 0.85],
            [-0.57, -0.975],
            [1.17, -0.42],
            [1.5, -0.15],
        ]
        # Set beam size (in arcsec)
        if observation == 'SCUBA':
            self.file_io.beam_size = 0.07
        elif observation == 'SPHERE':
            self.file_io.beam_size = 0.03
        # Take colorbar label from quantity id
        cbar_label = self.file_io.get_quantity_labels(i_quantity)
        # Define output pdf
        self.file_io.init_plot_output(
            'PI_emission_map_comparison', path=self.file_io.path['model'])
        for i_subplot in range(9):
            # Set paths of each simulation
            self.file_io.set_path_from_str(
                'plot', 'gg_tau_disk', model_list[i_subplot], 'dust')
            # Create pdf file if show_plot is not chosen and read map data from file
            plot_data, header, plot_data_type = self.file_io.read_emission_map(
                'polaris_detector_nr' + str(detector_index).zfill(4))
            # Create Matplotlib figure
            plot = Plot(self.model, self.parse_args, ax_unit='arcsec',
                        nr_x_images=2, nr_y_images=1, size_x=6 / 1.4298620007401583)
            # Take data for current quantity
            tbldata = plot_data[i_quantity, 0, :, :]
            # plot imshow
            plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=1, set_bad_to_min=True,
                             norm='LogNorm', vmin=vmin, vmax=vmax, cmap='magma', extend='neither')
            for i_pos, center_pos in enumerate(measurement_position_list):
                plot.plot_text(text_pos=center_pos,
                    text=str(i_pos + 1), ax_index=1, color='white', fontsize=10,
                    bbox=dict(boxstyle='circle, pad=0.2', facecolor='none', edgecolor='white', alpha=0.5))
            # Plot map description
            plot.plot_text(text_pos=[0.03, 0.97], relative_position=True,
                           text=model_descr[i_subplot], horizontalalignment='left',
                           verticalalignment='top', ax_index=1, color='white')
            # Observation
            if observation == 'SCUBA':
                hdulist = fits.open(
                    '/home/rbrauer/Documents/projects/005_gg_tau/near_infrared_imaging_paper/sub_pi.fits')
                tbldata = cropND(
                    hdulist[0].data.T, (500, 500), offset=[0, 30]) / 1.5e7
            elif observation == 'SPHERE':
                hdulist = fits.open(
                    '/home/rbrauer/Documents/projects/005_gg_tau/SPHERE_observation_miriam/GG_Tau_2016-11-191_I_POL.fits')
                tbldata = cropND(
                    hdulist[0].data.T, (390, 390), offset=[6, 14]) / 8e6
            # plot imshow
            plot.plot_imshow(tbldata, cbar_label=cbar_label, ax_index=0, set_bad_to_min=True,
                             norm='LogNorm', vmin=vmin, vmax=vmax, cmap='magma', extend='neither')
            for i_pos, center_pos in enumerate(measurement_position_list):
                if i_pos == 5:
                    pos = [center_pos[0], center_pos[1] + 0.1]
                else:
                    pos = center_pos
                plot.plot_text(text_pos=pos,
                    text=str(i_pos + 1), ax_index=0, color='white', fontsize=10,
                    bbox=dict(boxstyle='circle, pad=0.2', facecolor='none', edgecolor='white', alpha=0.5))
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)
