#!/usr/bin/env python
# -*- coding: utf-8 -*-

import locale

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid

from polaris_tools_modules.math import Math


class Plot:
    """This class creates various plots from results of POLARIS simulations.
    """

    def __init__(self, model, parse_args, image_type='image', nr_x_images=1, nr_y_images=1,
                 xlabel='', ylabel='', zlabel='', extent=None, limits=None, title='', ax_unit=None,
                 label_plane='xy', zoom_factor=None, zoom_x_factor=None, zoom_y_factor=None, with_cbar=True,
                 cmap_scaling=None, scale_axis_log=False, labelpad=None, language='english', size_x=6, size_y=4.5):
        """Initialisation of plot parameters.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
            parse_args: Provides all parameters chosen
                by user when executing PolarisTools.
            image_type (str): Type of plot image.
                (image, projection_3d, animation, healpix)
            nr_x_images (int): Number of images along the x-axis.
            nr_y_images (int): Number of images along the y-axis.
            ax_unit (str): Automatically format of the X- and Y-axis labels and the extent of the plot.
                (arcsec, au, pc, m, arb_units) If None -> xlabel, ylabel and extent will be used.
            label_plane (str): Set the axis label to the correct plane for midplane plots.
                'xy', 'xz', 'yz', 'rz'
            xlabel (str or list): Label of the x-axis.
            ylabel (str or list): Label of the y-axis.
            zlabel (str): Label of the z-axis. Only used if image_type == projection_3d.
            extent (List[float, float, float, float]): Extent of the plot
                [xmin, xmax, ymin, ymax].
            limits (List[float, float, float, float]): Limits of the plot
                [xmin, xmax, ymin, ymax].
            title (str): Title of the plot.
            zoom_factor (float): Zoom factor to reduce image to smaller region.
            zoom_x_factor (float): Zoom factor to reduce image to smaller region
                (Only for the X-axis).
            zoom_y_factor (float): Zoom factor to reduce image to smaller region
                (Only for the Y-axis).
            with_cbar (bool): Add space for colorbar?
            cmap_scaling (List): List with the scaling name first.
                ('symlog' + linthresh, 'power' + gamma, 'log')
            scale_axis_log (bool): Logarithmic scale of both axis (for imshow).
            labelpad (List): Padding of the axis label.
                Dimension of the list has to match the dimensions of the plot.
            language (str): Language for decimal separation.
            size_x (float): Size of the figure in x-direction [inches].
            size_y (float): Size of the figure in y-direction [inches].
        """

        # Get math module
        self.math = Math()

        ''' ##################################
        ######  Setting up parameters!  ######
        ################################## '''
        # Set plot parameter from user input
        if parse_args.ax_unit is not None:
            self.ax_unit = parse_args.ax_unit
        elif ax_unit is None:
            self.ax_unit = 'au'
        else:
            self.ax_unit = ax_unit
        if parse_args.zoom_factor is not None:
            self.zoom_x_factor = parse_args.zoom_factor
            self.zoom_y_factor = parse_args.zoom_factor
        else:
            self.zoom_x_factor = zoom_x_factor
            self.zoom_y_factor = zoom_y_factor
        if parse_args.zoom_x_factor is not None:
            self.zoom_x_factor = parse_args.zoom_x_factor
        if parse_args.zoom_y_factor is not None:
            self.zoom_y_factor = parse_args.zoom_y_factor
        if parse_args.x_scaling is not None:
            self.x_scaling = parse_args.x_scaling
        else:
            self.x_scaling = None
        if parse_args.y_scaling is not None:
            self.y_scaling = parse_args.y_scaling
        else:
            self.y_scaling = None
        if parse_args.cmap_scaling is not None:
            self.cmap_scaling = parse_args.cmap_scaling
        else:
            self.cmap_scaling = cmap_scaling
        if parse_args.extend is not None:
            self.extend = parse_args.extend
        else:
            self.extend = None
        if parse_args.vmin is not None:
            self.vmin = parse_args.vmin
        else:
            self.vmin = None
        if parse_args.vmax is not None:
            self.vmax = parse_args.vmax
        else:
            self.vmax = None
        if parse_args.bad_to_min:
            self.bad_to_min = True
        else:
            self.bad_to_min = False
        if any(i is not None for i in [parse_args.xmin, parse_args.xmax, parse_args.ymin, parse_args.ymax]):
            limits = [parse_args.xmin, parse_args.xmax,
                      parse_args.ymin, parse_args.ymax]

        # Set the axis formatter as class variable
        self.label_plane = label_plane
        # Create space in image for colorbar?
        if parse_args.no_cbar:
            self.with_cbar = False
        else:
            self.with_cbar = with_cbar

        # Set the style (Latex font)
        self.set_style(font=parse_args.font, font_size_env=parse_args.font_size_env,
                       gray_background=parse_args.gray_background)

        ''' ########################################
        ######  Creating the plots/subplots!  ######
        ######################################## '''
        # Set image type as class variable
        self.image_type = image_type
        # Create plot figure and axes for 2D or 3D plots
        if image_type == 'animation':
            # One single plot, but as an animation
            self.fig, self.ax_list = plt.subplots(figsize=(size_x, size_y))
            self.ax_list = np.array([self.ax_list])
            self.animation_images = []
        elif image_type == 'image':
            # Create figure with a size fitting to one or multiple plots
            self.fig = plt.figure(figsize=(np.power(nr_x_images, 0.7) * size_x,
                                           np.power(nr_y_images, 0.7) * size_y))
            if self.with_cbar:
                self.ax_list = AxesGrid(self.fig, 111, nrows_ncols=(nr_y_images, nr_x_images),
                                        axes_pad=0.1, cbar_mode='single', cbar_location='right', cbar_pad=0.1,
                                        cbar_size=str(4.5 / np.sqrt(nr_x_images)) + '%')
            else:
                self.ax_list = AxesGrid(self.fig, 111, nrows_ncols=(
                    nr_y_images, nr_x_images), axes_pad=0.05, aspect=False)
        elif image_type == 'projection_3d':
            # One 3D plot
            # noinspection PyUnresolvedReferences
            from mpl_toolkits.mplot3d import Axes3D
            self.fig = plt.figure(figsize=(1.8 * size_x, 1.8 * size_y))
            self.ax_list = self.fig.gca(projection='3d')
            self.ax_list = np.array([self.ax_list])
            self.ax_list[0].w_xaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))
            self.ax_list[0].w_yaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))
            self.ax_list[0].w_zaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))
        elif image_type == 'healpix':
            self.fig = None
            self.ax_list = []
        else:
            raise ValueError('Error: image_type is not known!')

        ''' ############################################
        ######  Creating plot labels and extent!  ######
        ############################################ '''
        self.update_all(extent, limits, model, label_plane,
                        nr_x_images, nr_y_images, xlabel, ylabel, zlabel, language, labelpad)

        ''' #######################################
        ######  Set various other settings!  ######
        ####################################### '''
        # Set logarithmic mode for both axis
        self.scale_axis_log = scale_axis_log
        if scale_axis_log:
            for ax in self.ax_list:
                ax.set_xscale('log')
                ax.set_yscale('log')

        # Set locale for decimal separator
        self.set_locale(language)

        # Set default colormap of imshow and similar plots
        self.cmap = parse_args.cmap

        # Add a title to the plot
        if title is not '':
            for ax in self.ax_list:
                ax.set_title(title)

        # Create linestyle cycler
        from itertools import cycle
        lines = ['-', '--', '-.', ':']
        self.linecycler = cycle(lines)

    def get_linestyle(self):
        """ Get linestyle from rotation"""
        return next(self.linecycler)

    @staticmethod
    def set_style(font_size_env=None, gray_background=False, font=None):
        """Set the style of the plot.

        Args:
            font_size_env (str): Name of the chosen font size environment.
                (paper, notebook, beamer, poster)
            gray_background (bool): Change background color if chosen to match with latex talk template.
            font (str): Define font other than default
        """
        # Set context of the plot (changes linewidth and font sizes)
        if font_size_env == 'beamer':
            mpl.rcParams.update({'font.size': 14})
        elif font_size_env == 'poster':
            mpl.rcParams.update({'font.size': 18})
        elif font_size_env == 'paper':
            mpl.rcParams.update({'font.size': 10})
        else:
            mpl.rcParams.update({'font.size': 12})
        # Change background color if chosen to match with latex talk template
        if gray_background:
            plt.rcParams['savefig.facecolor'] = '#FAFAFA'
        # Change plot style
        # plt.style.use(['seaborn-dark'])
        # Use Latex, specify the  Latex font and load siunitx
        preamble = [
            r'\usepackage[retain-zero-exponent=true, load=accepted, alsoload=astro]{siunitx}']
        if font == 'fira':
            preamble.append(r'\usepackage[sfdefault,scaled=.85]{FiraSans}')
            preamble.append(r'\usepackage{newtxsf}')
        elif font == 'bitstream':
            preamble.append(r'\usepackage[bitstream-charter]{mathdesign}')
            preamble.append(r'\usepackage[T1]{fontenc}')
        else:
            preamble.append(r'\usepackage{lmodern}')
        mpl.rc('text', usetex=True)
        mpl.rc('text.latex', preamble=preamble)
        mpl.rcParams['pgf.preamble'] = preamble
        mpl.rcParams['axes.formatter.limits'] = -3, 4

    def set_locale(self, language):
        """Set the locale of image plots (decimal separator).

        Args:
            language (str): Language for decimal separation.
        """
        preamble = mpl.rcParams['text.latex.preamble']
        if language == 'english':
            locale.setlocale(locale.LC_NUMERIC, "en_US.UTF-8")
            preamble.append(r'\sisetup{locale=US}')
        elif language == 'german':
            locale.setlocale(locale.LC_NUMERIC, "de_DE.UTF-8")
            preamble.append(r'\usepackage{icomma}')
            preamble.append(r'\sisetup{locale=DE}')
        else:
            raise ValueError('Do not know the chosen locale language!')
        plt.rcParams['axes.formatter.use_locale'] = True
        if not self.scale_axis_log:
            for i_ax in range(len(self.ax_list)):
                self.ax_list[i_ax].ticklabel_format(useLocale=True)
        mpl.rc('text.latex', preamble=preamble)

    def update_label(self, extent=None, label_plane=None, nr_x_images=1, nr_y_images=1,
                     xlabel=None, ylabel=None, zlabel=None, labelpad=None, language='english'):
        """Creates labels based on the chosen model.

        Args:
            extent (List[float, float, float, float]): Extent of the plot
                [xmin, xmax, ymin, ymax].
            label_plane (str): Set the axis label to the correct plane for midplane plots.
                'xy', 'xz', 'yz', 'rz'
            nr_x_images (int): Number of images along the x-axis.
            nr_y_images (int): Number of images along the y-axis.
            xlabel (str or list): Label of the x-axis.
            ylabel (str or list): Label of the y-axis.
            zlabel (str): Label of the z-axis. Only used if image_type == projection_3d.
            labelpad (List): Padding of the axis label.
                Dimension of the list has to match the dimensions of the plot.
            language (str): Language for decimal separation.

        Return:
            str: The axes that have to be adjusted automatically.
                ('x', 'y', 'xy', None)
        """
        # Set axis label anautomatic_axes = None
        automatic_axes = ''
        # Set axis label and extent of 2D or 3D plot
        if self.image_type == 'projection_3d':
            # Set label padding, since padding of 3D plots seems to be too low
            if labelpad is None:
                labelpad = [None, None, None]
            # Set label according to input label
            if xlabel != '' and ylabel != '' and zlabel != '':
                self.ax_list[0].set_xlabel(xlabel, labelpad=labelpad[0])
                self.ax_list[0].set_ylabel(ylabel, labelpad=labelpad[1])
                self.ax_list[0].set_zlabel(zlabel, labelpad=labelpad[2])
            else:
                if self.ax_unit == 'au':
                    if language == 'german':
                        ax_unit_str = 'AE'
                    else:
                        ax_unit_str = 'AU'
                elif self.ax_unit == 'm':
                    ax_unit_str = self.ax_unit
                elif self.ax_unit == 'arb_units':
                    if language == 'german':
                        ax_unit_str = r'willk. Einh.'
                    else:
                        ax_unit_str = r'arb.\ units'
                else:
                    raise ValueError(
                        'Error Axis label for 3D plot cannot be set correctly!')
                self.ax_list[0].set_xlabel(
                    r'$\Delta x\ [\mathsf{' + ax_unit_str + '}]$', labelpad=labelpad[0])
                self.ax_list[0].set_ylabel(
                    r'$\Delta y\ [\mathsf{' + ax_unit_str + '}]$', labelpad=labelpad[1])
                self.ax_list[0].set_zlabel(
                    r'$\Delta z\ [\mathsf{' + ax_unit_str + '}]$', labelpad=labelpad[2])
        elif self.image_type in ['image', 'animation'] and nr_x_images == 1 and nr_y_images == 1:
            if xlabel != '' and ylabel != '':
                pass
            elif self.ax_unit is not None:
                xlabel, ylabel, automatic_axes = self.create_axis_label(xlabel, ylabel,
                                                                        nr_x_images, nr_y_images, label_plane, language)
            else:
                # If not enough labels are defined, raise error
                raise ValueError(
                    'Error: ax_unit ist not set, but xlabel and ylabel neither!')
            # Set labels of the axes of the subplot
            for ax_index in range(len(self.ax_list)):
                self.ax_list[ax_index].set_xlabel(xlabel)
                self.ax_list[ax_index].set_ylabel(ylabel)
        elif self.image_type == 'image':
            if nr_x_images == 1 and nr_y_images > 1:
                if xlabel != '' and len(ylabel) == nr_y_images:
                    pass
                elif self.ax_unit is not None:
                    xlabel, ylabel, automatic_axes = self.create_axis_label(xlabel, ylabel,
                                                                            nr_x_images, nr_y_images, label_plane, language)
                else:
                    # If not enough labels are defined, raise error
                    raise ValueError(
                        'Error: The plot needs a list of ' + nr_y_images + ' ylabels and one xlabel!')
                # Set labels of the axes of the subplots
                for ax_index in range(len(self.ax_list)):
                    if ax_index == nr_y_images - 1:
                        self.ax_list[ax_index].set_xlabel(xlabel)
                    self.ax_list[ax_index].set_ylabel(ylabel[ax_index])
            elif nr_x_images > 1 and nr_y_images == 1:
                if ylabel != '' and len(xlabel) == nr_x_images:
                    pass
                elif self.ax_unit is not None:
                    xlabel, ylabel, automatic_axes = self.create_axis_label(xlabel, ylabel,
                                                                            nr_x_images, nr_y_images, label_plane, language)
                else:
                    # If not enough labels are defined, raise error
                    raise ValueError('Error: share_y plot needs a list of ' +
                                     nr_x_images + ' xlabels and one ylabel!')
                # Set labels of the axes of the subplots
                for ax_index in range(len(self.ax_list)):
                    if ax_index == 0:
                        self.ax_list[ax_index].set_ylabel(ylabel)
                    self.ax_list[ax_index].set_xlabel(xlabel[ax_index])
            elif nr_x_images > 1 and nr_y_images > 1:
                if len(xlabel) == nr_x_images and len(ylabel) == nr_y_images:
                    pass
                elif self.ax_unit is not None:
                    xlabel, ylabel, automatic_axes = self.create_axis_label(xlabel, ylabel,
                                                                            nr_x_images, nr_y_images, label_plane, language)
                else:
                    # If not enough labels are defined, raise error
                    raise ValueError('Error: share_both plot needs a list of ' +
                                     nr_x_images + ' xlabels and ' + nr_y_images + ' ylabels!')
                # Set labels of the axes of the subplots
                for ax_index in range(len(self.ax_list)):
                    if ax_index % nr_x_images == 0:
                        self.ax_list[ax_index].set_ylabel(
                            ylabel[int(ax_index / nr_x_images)])
                        if int(ax_index / nr_x_images) == nr_y_images - 1:
                            self.ax_list[ax_index].set_xlabel(xlabel[0])
                    elif ax_index > (nr_x_images * nr_y_images) - nr_x_images:
                        self.ax_list[ax_index].set_xlabel(
                            xlabel[int(ax_index - (nr_x_images * nr_y_images) + nr_x_images)])
        elif self.image_type == 'healpix':
            # Nothing needed here
            None
        else:
            raise ValueError('Error: image_type is not known!')
        # Return which axes where automatically set
        return automatic_axes

    def create_axis_label(self, xlabel, ylabel, nr_x_images, nr_y_images, label_plane, language):
        """Create axis label from auto label format.

        Args:
            xlabel (str): Input xlabel from plot class (to check if generation is neccessary).
            ylabel (str): Input ylabel from plot class (to check if generation is neccessary).
            nr_x_images (int): Dimension of xlabel output (1 -> only string).
            nr_y_images (int): Dimension of ylabel output (1 -> only string).
            label_plane (str): Set the axis label to the correct plane for midplane plots.
                'xy', 'xz', 'yz', 'rz'
            language (str): Language for decimal separation.

        Returns:
            Tuple with X-axis, Y-axis label, and a string that shows which axes are set automatically.
        """
        # Apply the correct plane for the x- and y-axis label
        if label_plane == 'xy':
            tmp_xlabel = r'$\Delta x\ '
            tmp_ylabel = r'$\Delta y\ '
        elif label_plane == 'xz':
            tmp_xlabel = r'$\Delta x\ '
            tmp_ylabel = r'$\Delta z\ '
        elif label_plane == 'yz':
            tmp_xlabel = r'$\Delta y\ '
            tmp_ylabel = r'$\Delta z\ '
        elif label_plane == 'rz':
            tmp_xlabel = r'$r\ '
            tmp_ylabel = r'$\Delta z\ '
        else:
            raise ValueError('Error: label_plane is not set correctly!')
        if self.ax_unit == 'arcsec':
            axis_label_unit = r'[\si{\arcsec}]$'
        elif self.ax_unit == 'au':
            if language == 'german':
                axis_label_unit = r'[\si{AE}]$'
            else:
                axis_label_unit = r'[\si{au}]$'
        elif self.ax_unit == 'pc':
            axis_label_unit = r'[\si{\parsec}]$'
        elif self.ax_unit == 'm':
            axis_label_unit = r'[\si{\metre}]$'
        elif self.ax_unit == 'arb_units':
            axis_label_unit = r'[\mathsf{arb.\ units}]$'
        else:
            raise ValueError('Error: ax_unit is not set correctly!')
        tmp_xlabel += axis_label_unit
        tmp_ylabel += axis_label_unit
        if xlabel != '' and ylabel == '':
            if nr_x_images > 1 and (isinstance(xlabel, str) or len(xlabel) != nr_x_images):
                raise ValueError('xlabel was not set correctly!')
            if nr_y_images > 1:
                tmp_ylabel = [tmp_ylabel for i in range(nr_y_images)]
            return xlabel, tmp_ylabel, 'y'
        elif xlabel == '' and ylabel != '':
            if nr_x_images > 1:
                tmp_xlabel = [tmp_xlabel for i in range(nr_x_images)]
            if nr_y_images > 1 and (isinstance(ylabel, str) or len(ylabel) != nr_y_images):
                raise ValueError('ylabel was not set correctly!')
            return tmp_xlabel, ylabel, 'x'
        elif xlabel == '' and ylabel == '':
            if nr_x_images > 1:
                tmp_xlabel = [tmp_xlabel] * nr_x_images
            if nr_y_images > 1:
                tmp_ylabel = [tmp_ylabel] * nr_y_images
            return tmp_xlabel, tmp_ylabel, 'xy'
        else:
            return xlabel, ylabel, ''

    def update_extent(self, extent=None, model=None, automatic_axes=None):
        """Update extent of the plot.

        Args:
            extent (List[float, float, float, float]): Extent of the plot
                [xmin, xmax, ymin, ymax].
            model: Handles the model space including various
                quantities such as the density distribution.
            automatic_axes (str): Set the axes that have to be adjusted automatically.
                ('x', 'y', 'xy', None)
        """
        if extent is not None:
            self.extent = extent
            if 'x' in automatic_axes:
                self.extent[0:2] = self.math.length_conv(
                    extent[0:2], self.ax_unit)
            if 'y' in automatic_axes:
                self.extent[2:4] = self.math.length_conv(
                    extent[2:4], self.ax_unit)
        elif self.ax_unit is not None:
            if self.ax_unit == 'arb_units':
                # Arbitrary units should go from -1 to 1.
                if self.image_type == 'projection_3d':
                    self.extent = [-1., 1., -1., 1., -1., 1.]
                elif automatic_axes == 'xy':
                    self.extent = [-1., 1., -1., 1.]
                elif automatic_axes == 'x':
                    self.extent = [-1., 1., None, None]
                elif automatic_axes == 'y':
                    self.extent = [None, None, -1., 1.]
                elif automatic_axes is None:
                    self.extent = [None, None, None, None]
            elif model is not None:
                # Get extent from model
                radius_x = model.tmp_parameter['radius_x_' + self.ax_unit]
                radius_y = model.tmp_parameter['radius_y_' + self.ax_unit]
                if self.image_type == 'projection_3d':
                    self.extent = [-radius_x, radius_x, -
                                   radius_x, radius_x, -radius_x, radius_x]
                elif automatic_axes == 'xy':
                    self.extent = [-radius_x, radius_x, -radius_y, radius_y]
                elif automatic_axes == 'x':
                    self.extent = [-radius_x, radius_x, None, None]
                elif automatic_axes == 'y':
                    self.extent = [None, None, -radius_y, radius_y]
                else:
                    self.extent = [None, None, None, None]
            else:
                raise ValueError(
                    'Without defined model, the update of the extent is not possible!')
        else:
            if self.image_type == 'projection_3d':
                self.extent = [None, None, None, None, None, None]
            else:
                self.extent = [None, None, None, None]

    def update_limits(self, limits=None):
        """If the zoom factor is set, modify the limits accordingly. If no zoom factor but limits are chosen,
        apply the limits. Else, apply not limit and use the extent.

        Args:
            limits (List[float, float, float, float]): Limit of the plot
                [xmin, xmax, ymin, ymax].
        """
        # Set limit according to zoom factor
        if (self.zoom_x_factor is not None or self.zoom_y_factor is not None) \
                and self.extent is not None and self.image_type != 'projection_3d':
            self.limits = [None, None, None, None]
            for i_limits in range(len(self.extent)):
                if self.extent[i_limits] is None:
                    raise ValueError(
                        'Without a defined extent, no zoom factor can be applied!')
                else:
                    if i_limits in [0, 1] and self.zoom_x_factor is not None:
                        self.limits[i_limits] = self.extent[i_limits] / \
                            self.zoom_x_factor
                    elif i_limits in [2, 3] and self.zoom_y_factor is not None:
                        self.limits[i_limits] = self.extent[i_limits] / \
                            self.zoom_y_factor
        elif limits is not None:
            self.limits = limits
        else:
            self.limits = None

    def update_all(self, extent=None, limits=None, model=None, label_plane=None,
                   nr_x_images=1, nr_y_images=1, xlabel=None, ylabel=None, zlabel=None, labelpad=None, language='english'):
        """Creates labels based on the chosen model.

        Args:
            extent (List[float, float, float, float]): Extent of the plot
                [xmin, xmax, ymin, ymax].
            limits (List[float, float, float, float]): Limits of the plot
                [xmin, xmax, ymin, ymax].
            model: Handles the model space including various
                quantities such as the density distribution.
            label_plane (str): Set the axis label to the correct plane for midplane plots.
                'xy', 'xz', 'yz', 'rz'
            nr_x_images (int): Number of images along the x-axis.
            nr_y_images (int): Number of images along the y-axis.
            xlabel (str or list): Label of the x-axis.
            ylabel (str or list): Label of the y-axis.
            zlabel (str): Label of the z-axis. Only used if image_type == projection_3d.
            labelpad (List): Padding of the axis label.
                Dimension of the list has to match the dimensions of the plot.
            language (str): Language for decimal separation.
        """
        # Update the labels of the image
        automatic_axes = self.update_label(extent, label_plane, nr_x_images, nr_y_images,
                                           xlabel, ylabel, zlabel, labelpad, language)
        # Update the extent of the image
        self.update_extent(extent, model, automatic_axes)
        # Update the limits of the image
        self.update_limits(limits)

    def set_xlabel(self, xlabel, ax_index=None):
        """Set xlabel after creation of plot instance.

        Args:
            xlabel (str): label of x-axis.
            ax_index (int): Index of subplot image.
        """
        if ax_index is not None:
            self.ax_list[ax_index].set_xlabel(xlabel)
        else:
            for ax in self.ax_list:
                ax.set_xlabel(xlabel)

    def set_ylabel(self, ylabel, ax_index=None):
        """Set ylabel after creation of plot instance.

        Args:
            ylabel (str): label of y-axis.
            ax_index (int): Index of subplot image.
        """
        if ax_index is not None:
            self.ax_list[ax_index].set_ylabel(ylabel)
        else:
            for ax in self.ax_list:
                ax.set_ylabel(ylabel)

    def create_zoom_axis(self, ax_index=0, zoom_factor=3, loc=1, extent=None, model=None):
        """
        Create axis to plot zoomed image on.

        Args:
            ax_index (int): Index of subplot image.
            zoom_factor (int): Size of the zoomed image.
            loc (int): Position of the zoomed image
            extent (List[float, float, float, float]): Extent zoomed image.
                [xmin, xmax, ymin, ymax].
            model: Instance of the model (to obtain the extent if not set)
        """
        # Import modules
        from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
        from mpl_toolkits.axes_grid1.inset_locator import mark_inset
        # Create sub axis
        self.ax_list = np.hstack([self.ax_list,
                                  zoomed_inset_axes(self.ax_list[ax_index], zoom_factor, loc=loc)]).ravel()
        # Update the extent of the zoomed image
        self.update_extent(extent, model, automatic_axes='xy')
        # Reomove ticks from zoom region
        self.remove_ticks(-1)
        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        mark_inset(self.ax_list[ax_index], self.ax_list[-1],
                   loc1=2, loc2=4, fc="none", ec="0.5")

    def plot_line(self, xdata, ydata, ax_index=0, log=None, step=False,
                  fill_between=False, no_ticks=False, no_grid=False, **args):
        """Plot 2D line from xdata and ydata.

        Args:
            xdata: X-axis position of the data points.
            ydata: Y-axis position of the data points.
            ax_index (int): Index of subplot image.
            log (str): If set, x- and/or y-axis is plotted logarithmically.
                (any combination of x, y)
            step (bool): Use steps instead of linearly connected lines.
            fill_between (bool): fill below curve.
            no_ticks (bool): Set true if not ticks shall be plotted.
            no_grid (bool): Set true to not show grid lines.
        """
        # Get user input for log
        if self.x_scaling == 'log' and self.y_scaling == 'log':
            log = 'xy'
        elif self.x_scaling == 'log' and self.y_scaling == 'linear':
            log = 'x'
        elif self.y_scaling == 'log' and self.x_scaling == 'linear':
            log = 'y'
        elif self.x_scaling == 'linear' and self.y_scaling == 'linear':
            log = None

        # Choose plot method depending on the log scaling
        if 'yerr' in args.keys() or 'xerr' in args.keys():
            plot_func = self.ax_list[ax_index].errorbar
        elif step:
            plot_func = self.ax_list[ax_index].step
            args['where'] = 'mid'
            if fill_between:
                self.ax_list[ax_index].fill_between(
                    xdata, ydata, step='mid', alpha=0.4)
        elif log == 'xy':
            plot_func = self.ax_list[ax_index].loglog
        elif log == 'x':
            plot_func = self.ax_list[ax_index].semilogx
        elif log == 'y':
            plot_func = self.ax_list[ax_index].semilogy
        else:
            plot_func = self.ax_list[ax_index].plot

        if fill_between and step:
            self.ax_list[ax_index].fill_between(
                xdata, ydata, alpha=0.4, step='mid', color='grey')
        elif fill_between:
            self.ax_list[ax_index].fill_between(
                xdata, ydata, alpha=0.4, color='grey')

        # Plot the 2D line data
        plot_func(xdata, ydata, **args)

        # Enable grid for better reading
        if not no_grid:
            self.ax_list[ax_index].grid(linestyle=':')

        # Remove ticks if chosen
        if no_ticks:
            self.ax_list[ax_index].set_xticklabels([])
            self.ax_list[ax_index].set_yticklabels([])

    def plot_hist(self, data, ax_index=0, label='', hist_bins=20, log=False, hist_type='bar', color=None,
                  bin_range=None):
        """Plot a 1D histogram.

        Args:
            data: 1 dimensional list of numbers.
            ax_index (int): Index of subplot image.
            label (str): Label of the histogram (for legend).
            hist_bins (int): Number of bins used to create the histogram
            log (str): If 'y', the histogram axis will be set to a log scale.
            hist_type (str): The type of histogram.
            color (str): Color of the histogram bars.
            bin_range (List): The lower and upper range of the bins.
        """
        # Change to log mode if user demands it
        if self.y_scaling == 'log':
            log = True

        # Plot the histogram
        if color is not None:
            self.ax_list[ax_index].hist(data, label=label, bins=hist_bins, log=log, histtype=hist_type, color=color,
                                        range=bin_range)
        else:
            self.ax_list[ax_index].hist(
                data, label=label, bins=hist_bins, log=log, histtype=hist_type, range=bin_range)

    def plot_marker(self, marker_pos, marker, ax_index=0, **args):
        """Plot a marker at a certain position.

        Args:
            text_pos (List[float, float]): List with the 2D marker position.
            marker (str): Marker symbol.
            ax_index (int): Index of subplot image.
            markersize (float): Size of the marker.
        """
        if len(marker_pos) != 2:
            raise ValueError(
                'Marker position of plot_marker has to consist of only 2 values!')

        self.ax_list[ax_index].plot(
            marker_pos[0], marker_pos[1], marker=marker, **args)

    def plot_bar(self, left, height, width=0.8, bottom=None, ax_index=0, orientation='vertical', align='center',
                 label='', color=None, edgecolor=None, alpha=1.0):
        """Plot single bars.

        Args:
            left: The position of the left part of the bars.
            height: The height of the bars.
            width: The width of the bars.
            bottom: The position of the bottom of the bars.
            ax_index (int): Index of subplot image.
            orientation (str): The orientation of the bar plot.
            align (str): The alignment of the bars with their x-axis position.
            label (str): Label of the histogram (for legend).
            color (str): Color of the line and markers.
            edgecolor (str): The color of the bar edges.
            alpha (float): The alpha value of the plot color.
        """
        # Plot the bars
        self.ax_list[ax_index].bar(left, height, width=width, bottom=bottom, orientation=orientation, align=align,
                                   edgecolor=edgecolor, alpha=alpha, label=label, color=color)

    def plot_pcolor(self, X, Y, tbldata, ax_index=0, cbar_label='', plot_cbar=True, extend=None, norm='Normalize',
                    cmap=None, vmin=None, vmax=None, linthresh=None, gamma=None, set_bad_to_min=False):
        """Plot 3D data in 2D colorcoded form.

        Args:
            X: Numpy array with 2 dimensions for the pixel coordinate X.
            Y: Numpy array with 2 dimensions for the pixel coordinate Y.
            tbldata: Numpy array with data for color plotting.
            ax_index (int): Index of subplot image.
            extend (str): Extend of the colorbar if clipping is used.
                If vmin is larger than the smallest value ('min').
                If vmax is smaller than the largest value ('max').
                If both is true ('both').
            cbar_label (str): Label of the colorbar.
            plot_cbar (bool): Plot the colorbar?
            norm (str): Type of color normalization (Normalize, LogNorm, SymLogNorm, ...).
            cmap: Name or instance of the colormap.
            vmin (float): Minimum value of the colorbar.
            vmax (float): Maximum value of the colorbar.
            linthresh (float): Limit under which 'SymLogNorm' is linear.
            gamma (float): Exponent for exponential  normalization.
            set_bad_to_min (bool): Set the bad color to the color of the
                minimum value of the colorbar.
        """
        # Convert X and Y from metre to ax_unit
        X = self.math.length_conv(X, self.ax_unit)
        Y = self.math.length_conv(Y, self.ax_unit)

        # Use default colormap if no one is defined
        if self.cmap is None:
            if cmap is not None:
                self.cmap = cmap
            else:
                self.cmap = 'viridis'

        # Change vmin and vmax if user demands it
        if self.vmin is not None:
            vmin = self.vmin
        if self.vmax is not None:
            vmax = self.vmax

        # Change to log mode if user demands it
        if self.cmap_scaling is not None:
            if self.cmap_scaling[0] == 'symlog' and len(self.cmap_scaling) == 2:
                norm = 'SymLogNorm'
                linthresh = float(self.cmap_scaling[1])
            elif self.cmap_scaling[0] == 'power' and len(self.cmap_scaling) == 2:
                if (vmin is None or vmin >= 0) and (vmax is None or vmax > 0) and \
                        tbldata.any():
                    norm = 'PowerNorm'
                    gamma = float(self.cmap_scaling[1])
            elif self.cmap_scaling[0] == 'log':
                if not any(i < 0. for i in tbldata.flatten()) \
                        and not all(i == 0. for i in tbldata.flatten()) \
                        and (vmin is None or vmin >= 0) \
                        and (vmax is None or vmax > 0):
                    norm = 'LogNorm'

        # Set the norm related to chosen norm
        from matplotlib.colors import Normalize, LogNorm, PowerNorm, SymLogNorm
        if norm is 'Normalize':
            norm = Normalize(vmin=vmin, vmax=vmax)
        elif norm is 'LogNorm':
            norm = LogNorm(vmin=vmin, vmax=vmax)
        elif norm is 'PowerNorm' and gamma is not None:
            norm = PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax, clip=True)
        elif norm is 'SymLogNorm' and linthresh is not None:
            norm = SymLogNorm(linthresh=linthresh, vmin=vmin, vmax=vmax)
        else:
            raise AttributeError(
                'The chosen norm for imshow plot is not found or linthresh is not set!')

        # Add colorbar extend if extend is set
        if self.extend is not None:
            extend = self.extend
        elif extend is None:
            # Set extend if vmin or vmax is set
            if vmax is not None and vmin is not None:
                extend = 'both'
            elif vmin is not None:
                extend = 'min'
            elif vmax is not None:
                extend = 'max'
            else:
                extend = 'neither'

        # Change colormap bad values to lowest values
        colormap = plt.get_cmap(self.cmap)
        if set_bad_to_min or self.bad_to_min:
            colormap.set_bad(colormap(0))

        # Plot 2D color plot
        self.image = self.ax_list[ax_index].pcolor(
            X, Y, tbldata, cmap=colormap, norm=norm)

        # Plot colorbar is chosen
        if self.with_cbar and plot_cbar:
            self.plot_colorbar(ax_index=ax_index,
                               label=cbar_label, extend=extend)

    def plot_imshow(self, tbldata, ax_index=0, cbar_label='', plot_cbar=True, extend=None, norm='Normalize',
                    cmap=None, interpolation='nearest', origin='lower', extent=None, aspect='auto',
                    vmin=None, vmax=None, linthresh=None, gamma=None, set_bad_to_min=False):
        """Plot 3D data in 2D colorcoded form.

        Args:
            tbldata: Numpy array with 2 dimensions for the pixel coordinates. The values
                are plotted as colors.
            ax_index (int): Index of subplot image.
            extend (str): Extend of the colorbar if clipping is used.
                If vmin is larger than the smallest value ('min').
                If vmax is smaller than the largest value ('max').
                If both is true ('both').
            cbar_label (str): Label of the colorbar.
            plot_cbar (bool): Plot the colorbar?
            norm (str): Type of color normalization (Normalize, LogNorm, SymLogNorm, ...).
            cmap: Name or instance of the colormap.
            interpolation (str): Type of interpolation between pixel.
            origin (str): Starting point where the data is aligned on.
            extent (List[float, float, float, float]): Extent of the plot
                [xmin, xmax, ymin, ymax].
            aspect (str): Aspect between the x- and y-axis.
            vmin (float): Minimum value of the colorbar.
            vmax (float): Maximum value of the colorbar.
            linthresh (float): Limit under which 'SymLogNorm' is linear.
            gamma (float): Exponent for exponential  normalization.
            set_bad_to_min (bool): Set the bad color to the color of the
                minimum value of the colorbar.
        """
        # Set preset extent if not set by this function
        if extent is None:
            extent = self.extent
        # Use default colormap if no one is defined
        if self.cmap is None:
            if cmap is not None:
                self.cmap = cmap
            else:
                self.cmap = 'viridis'

        # Change vmin and vmax if user demands it
        if self.vmin is not None:
            vmin = self.vmin
        if self.vmax is not None:
            vmax = self.vmax

        # Change to log mode if user demands it
        if self.cmap_scaling is not None:
            if self.cmap_scaling[0] == 'symlog' and len(self.cmap_scaling) == 2:
                norm = 'SymLogNorm'
                linthresh = float(self.cmap_scaling[1])
            elif self.cmap_scaling[0] == 'power' and len(self.cmap_scaling) == 2:
                if (vmin is None or vmin >= 0) and (vmax is None or vmax > 0) and \
                        tbldata.any():
                    norm = 'PowerNorm'
                    gamma = float(self.cmap_scaling[1])
            elif self.cmap_scaling[0] == 'log':
                if not any(i < 0. for i in tbldata.flatten()) \
                        and not all(i == 0. for i in tbldata.flatten()) \
                        and (vmin is None or vmin >= 0) \
                        and (vmax is None or vmax > 0):
                    norm = 'LogNorm'

        # Set the norm related to chosen norm
        from matplotlib.colors import Normalize, LogNorm, PowerNorm, SymLogNorm
        if norm is 'Normalize':
            norm = Normalize(vmin=vmin, vmax=vmax)
        elif norm is 'LogNorm':
            norm = LogNorm(vmin=vmin, vmax=vmax)
        elif norm is 'PowerNorm' and gamma is not None:
            norm = PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax, clip=True)
        elif norm is 'SymLogNorm' and linthresh is not None:
            norm = SymLogNorm(linthresh=linthresh, vmin=vmin, vmax=vmax)
        else:
            raise AttributeError(
                'The chosen norm for imshow plot is not found or linthresh is not set!')

        # Add colorbar extend if extend is set
        if self.extend is not None:
            extend = self.extend
        elif extend is None:
            # Set extend if vmin or vmax is set
            if vmax is not None and vmin is not None:
                extend = 'both'
            elif vmin is not None:
                extend = 'min'
            elif vmax is not None:
                extend = 'max'
            else:
                extend = 'neither'

        # Change colormap bad values to lowest values
        if '_half' in self.cmap:
            colormap = self.truncate_colormap(
                plt.get_cmap(self.cmap.replace('_half', '')), 0.0, 0.5)
        else:
            colormap = plt.get_cmap(self.cmap)
        if set_bad_to_min or self.bad_to_min:
            colormap.set_bad(colormap(0))

        # Plot 2D color plot
        self.image = self.ax_list[ax_index].imshow(tbldata.T, cmap=colormap, interpolation=interpolation, norm=norm,
                                                   origin=origin, extent=extent, aspect=aspect)
        # Plot colorbar is chosen
        if self.with_cbar and plot_cbar:
            self.plot_colorbar(ax_index=ax_index,
                               label=cbar_label, extend=extend)
        # Add imshow plot to list for animation
        if self.image_type == 'animation':
            self.animation_images.append(self.image)

    def plot_healpix(self, wmap_map, ax_index=0, cbar_label='', plot_cbar=True,
                     norm=None, cmap=None, extent=None, vmin=None, vmax=None, title='', set_bad_to_min=False):
        """Plot healpix data in various ways.

        Args:
            wmap_map: Numpy array with the quantities and a 1D data axis.
            ax_index (int): Index of subplot image.
            cbar_label (str): Label of the colorbar.
            plot_cbar (bool): Plot the colorbar?
            norm (str): Type of color normalization (Normalize, LogNorm, SymLogNorm, ...).
            cmap: Name or instance of the colormap.
            interpolation (str): Type of interpolation between pixel.
            vmin (float): Minimum value of the colorbar.
            vmax (float): Maximum value of the colorbar.
            linthresh (float): Limit under which 'SymLogNorm' is linear.
            title (str): Title of plot.
            set_bad_to_min (bool): Set the bad color to the color of the
                minimum value of the colorbar.
        """
        # Load healpy module
        import healpy as hp

        # Set preset extent if not set by this function
        if extent is None:
            extent = self.extent
        # Use default colormap if no one is defined
        if self.cmap is None:
            if cmap is not None:
                self.cmap = cmap
            else:
                self.cmap = 'viridis'

        # Change vmin and vmax if user demands it
        if self.vmin is not None:
            vmin = self.vmin
        if self.vmax is not None:
            vmax = self.vmax

        # Change to log mode if user demands it
        if self.cmap_scaling is not None:
            if self.cmap_scaling[0] == 'log':
                if not any(i < 0. for i in wmap_map.flatten()) \
                        and not all(i == 0. for i in wmap_map.flatten()) \
                        and (vmin is None or vmin >= 0) \
                        and (vmax is None or vmax > 0):
                    norm = 'log'
            elif self.cmap_scaling[0] == 'hist':
                norm = 'hist'

        # Adjust colomap for plot
        colormap = plt.get_cmap(self.cmap)
        # Healpy adjustments
        colormap.set_under('w')
        colormap.set_over(colormap(1.0))

        # Change colormap bad values to lowest values
        if set_bad_to_min or self.bad_to_min:
            colormap.set_bad(colormap(0))
        else:
            colormap.set_bad('gray')

        # Calculate size of the plot to fit the size of the healpix map
        xsize = np.ceil(np.sqrt(wmap_map.size) / 800.) * 800

        # Plot 2D color plot
        self.image = hp.mollview(wmap_map, unit=cbar_label, title=title, min=vmin, xsize=xsize,
                                 max=vmax, norm=norm, cbar=(self.with_cbar and plot_cbar), cmap=colormap, format=r'$\SI{%1.2e}{}$')

    def plot_quiver(self, vec_field_data, index=None, ax_index=0, units='width', scale_units='width', width=0.004,
                    headwidth=0., headlength=0., headaxislength=0., pivot='middle', color=None, cmap=None,
                    const_quiver=False, vmin=None, vmax=None, xmin=None, xmax=None, ymin=None, ymax=None):
        """Plot a 2D image from 3D-vector field.

        Args:
            vec_field_data: Numpy array with 2 dimensions for the pixel coordinates and
                one dimension for the vector components.
            index (list): List of two integers that define the indexes of the vec_field that will
                be used for the arrows (0 -> x, 1 -> y, 2 -> z).
            ax_index (int): Index of subplot image.
            units (str): Scaling unit for arrow width.
            scale_units (str): Scaling unit for arrow length.
            width (float): width of the arrows.
            headwidth (float): Width of the arrowhead.
            headlength (float): Length of the arrowhead.
            headaxislength (float): Length of the arrowheadaxis.
            pivot (str): Position, where the arrow is mounted to ('middle', ...).
            color (str): Color of the arrows.
            cmap: Name or Instance of the colormap to define the color of the arrows,
                if no color is set.
            const_quiver (bool): If True, all quiver vectors have the same length.
            vmin (float): Define minimum vector size.
            vmax (float): Define maximum vector size.
            xmin (float): Define minimum x value where vectors will be plotted.
            xmax (float): Define maximum x value where vectors will be plotted.
            ymin (float): Define minimum y value where vectors will be plotted.
            ymax (float): Define maximum y value where vectors will be plotted.
        """
        # Use default colormap if no one is defined
        if self.cmap is None:
            if cmap is not None:
                self.cmap = cmap
            else:
                self.cmap = 'viridis'

        # Add minimum vector size if vmin is set
        if self.vmin is not None:
            vmin = self.vmin
        #if self.vmax is not None:
        #    vmax = self.vmax
        if index is None:
            if 'xy' in self.label_plane:
                index = [0, 1]
            elif 'xz' in self.label_plane:
                index = [0, 2]
            elif 'yz' in self.label_plane:
                index = [1, 2]
            elif 'rz' in self.label_plane:
                index = [0, 2]
            else:
                raise ValueError(
                    'Error: index of the vector field components for the arrows are not set!')

        #: int: Number of pixel in one direction
        bins = len(vec_field_data[:, 0, 0])

        #: 2D array for the x-axis position
        x = np.zeros((bins, bins))
        #: 2D array for the y-axis position
        y = np.zeros((bins, bins))

        #: 2D array for the vector length in the x-direction
        vx = np.zeros((bins, bins))
        #: 2D array for the vector length in the y-direction
        vy = np.zeros((bins, bins))

        #: 2D array for the total vector length
        z = np.zeros((bins, bins))

        # If extent is not set, the quiver plot cannot be generated
        if any(self.extent[i] is None for i in range(4)):
            raise ValueError(
                'Quiver plot cannot be generated without image extent!')
        #: float: distance between two pixel on the x-axis.
        dx = self.extent[1] - self.extent[0]
        #: float: distance between two pixel on the y-axis.
        dy = self.extent[3] - self.extent[2]

        for i_x in range(bins):
            for i_y in range(bins):
                # Set the x-axis positions
                x[i_x, i_y] = float(self.extent[0]) + 1 / \
                    bins * (i_x + 0.5) * dx
                # Set the y-axis positions
                y[i_x, i_y] = float(self.extent[2]) + 1 / \
                    bins * (i_y + 0.5) * dy
                # Set the vector component in the x-direction
                vx[i_x, i_y] = vec_field_data[i_x, i_y, index[0]]
                # Set the vector component in the y-direction
                vy[i_x, i_y] = vec_field_data[i_x, i_y, index[1]]
                # Calculate the total length of the vector
                z[i_x, i_y] = np.sqrt(vx[i_x, i_y] ** 2 + vy[i_x, i_y] ** 2)

        for i_x in range(bins):
            for i_y in range(bins):
                # Plot vector only if position and length are not limited
                plot_vector = True
                if xmin is not None and xmin > x[i_x, i_y]:
                    plot_vector = False
                if xmax is not None and x[i_x, i_y] > xmax:
                    plot_vector = False
                if ymin is not None and ymin > y[i_x, i_y]:
                    plot_vector = False
                if ymax is not None and y[i_x, i_y] > ymax:
                    plot_vector = False
                if vmin is not None and z[i_x, i_y] < vmin:
                    plot_vector = False
                if vmax is not None:
                    norm = vmax  # np.max(z[i_x, i_y], vmax)
                else:
                    norm = np.nanmax(z)  
                # Plot vector if plot_vector is true
                if plot_vector:
                    if const_quiver:
                        # Set the length of each vector to 1
                        if z[i_x, i_y] > 0.:
                            vx[i_x, i_y] /= z[i_x, i_y]
                            vy[i_x, i_y] /= z[i_x, i_y]
                        else:
                            vx[i_x, i_y] = np.nan
                            vy[i_x, i_y] = np.nan
                    else:
                        # Normalize the length of the vector to 1
                        if norm > 0.:
                            vx[i_x, i_y] /= norm
                            vy[i_x, i_y] /= norm
                        else:
                            vx[i_x, i_y] = 0.
                            vy[i_x, i_y] = 0.
                else:
                    # Set the length to nan, to ignore these vectors
                    vx[i_x, i_y] = np.nan
                    vy[i_x, i_y] = np.nan

        #: float: Scale of the arrow length
        scale = 1.2 * bins
        # Increase the size of the arrows if a zoom factor is used
        if self.zoom_x_factor is not None or self.zoom_y_factor is not None:
            scale /= min(self.zoom_x_factor, self.zoom_y_factor)
        # Create colormap depending on chosen colormap (2 discrete colors)
        if '_r' in self.cmap:
            colormap = plt.get_cmap(self.cmap.replace('_r', ''), 2)
        else:
            colormap = plt.get_cmap(self.cmap + '_r', 2)

        # Plot the vector field with a constant color or with a colormap
        if color is not None:
            quiver_image = self.ax_list[ax_index].quiver(x, y, vx, vy, units=units, scale=scale, width=width,
                                                         scale_units=scale_units, headwidth=headwidth,
                                                         headlength=headlength, headaxislength=headaxislength,
                                                         pivot=pivot, color=color)
        else:
            quiver_image = self.ax_list[ax_index].quiver(x, y, vx, vy, z, units=units, scale=scale, width=width,
                                                         scale_units=scale_units, headwidth=headwidth,
                                                         headlength=headlength, headaxislength=headaxislength,
                                                         pivot=pivot, cmap=colormap)
        # Add quiver plot to list for animation
        if self.image_type == 'animation':
            self.animation_images.append(quiver_image)

    def plot_pol_vector_text(self, vec_per_width, max_pol_degree, ax_index=0, color='black', round_lvl=0):
        """Plot a line equivalent to the largest vector size with the
        corresponding degree of polarization.

        Args:
            vec_per_width (int): Amount of vectors per image width.
            max_pol_degree (float): Maximum degree of polarization
            ax_index (int): Index of subplot image.
            color (str): Color of the text and line.
            round_lvl (int): How many digits for rounding the number.
        """
        # Position of the text object in image coordinates
        if self.limits is not None:
            text_pos = [self.limits[0] + (self.limits[1] - self.limits[0]) * 0.9,
                        self.limits[2] + (self.limits[3] - self.limits[2]) * 0.05]
            line_pos = [self.limits[0] + (self.limits[1] - self.limits[0]) * 0.9,
                        self.limits[2] + (self.limits[3] - self.limits[2]) * 0.045]
            # Length of the longest polarization vector
            length = (self.limits[1] - self.limits[0]) / (1.2 * vec_per_width)
        elif self.extent is not None:
            text_pos = [self.extent[0] + (self.extent[1] - self.extent[0]) * 0.9,
                        self.extent[2] + (self.extent[3] - self.extent[2]) * 0.05]
            line_pos = [self.extent[0] + (self.extent[1] - self.extent[0]) * 0.9,
                        self.extent[2] + (self.extent[3] - self.extent[2]) * 0.045]
            # Length of the longest polarization vector
            length = (self.extent[1] - self.extent[0]) / (1.2 * vec_per_width)
        else:
            raise ValueError('Neither the extent nor the limits are set!')
        # Increase the size of the line if a zoom factor is used
        if self.zoom_x_factor is not None or self.zoom_y_factor is not None:
            length *= min(self.zoom_x_factor, self.zoom_y_factor)
        # Text including the maximum degree of polarization
        if round_lvl > 0:
            text = r'$\SI{' + str(round(max_pol_degree,
                                        round_lvl)) + r'}{\percent}$'
        else:
            text = r'$\SI{' + str(int(max_pol_degree)) + r'}{\percent}$'

        # Plot the text
        self.plot_text(text_pos, text, color=color,
                       verticalalignment='bottom', ax_index=ax_index, zorder=1)
        self.text.set_bbox(
            dict(facecolor='black', alpha=0.8, edgecolor='black'))
        # Plot the line
        self.plot_line([line_pos[0] - (length / 2.), line_pos[0] + (length / 2.)], [line_pos[1], line_pos[1]],
                       log='never', color=color, ax_index=ax_index, zorder=2, no_grid=True)

    def plot_contour(self, tbldata, ax_index=0, xaxis=None, yaxis=None, origin='lower',
                     interpolation='nearest', linestyles='-', extent=None, levels=None, colors=None, cmap=None,
                     fontsize=None, label_type='default'):
        """Plot contour plot for different levels.

        Args:
            tbldata: Numpy array with 2 dimensions for the pixel coordinates. The values
                are used to find the contour positions.
            ax_index (int): Index of subplot image.
            xaxis: Array of x-axis position.
            yaxis: Array of y-axis position.
            origin (str): Starting point where the data is aligned on.
            interpolation (str): Type of interpolation between pixel.
            linestyles (str): Style of the lines ('-', ':', ...).
            extent (List[float, float, float, float]): Extent of the plot
                [xmin, xmax, ymin, ymax].
            levels (List): List of values at which the contour line
                has to be plotted.
            colors (str): Color of the contour lines.
            cmap: Name or Instance of the colormap.
            fontsize (float): Change the fontsize of the ticks.
            label_type (str): Type of contour label.
                ('percentage', 'percentage_1_decimal', 'tau', 'tau_1_decimal', '0_decimal', 'default', None)
        """
        # Set preset extent if not set by this function
        if extent is None:
            extent = self.extent
        # Use default colormap if no one is defined
        if self.cmap is None:
            if cmap is not None:
                self.cmap = cmap
            else:
                self.cmap = 'viridis'

        # Plot the contour lines with a constant color or with a colormap
        if xaxis is not None and yaxis is not None:
            if colors is not None:
                contour = self.ax_list[ax_index].contour(xaxis, yaxis, tbldata.T, origin=origin,
                                                         interpolation=interpolation, linestyles=linestyles,
                                                         levels=levels, colors=colors)
            else:
                contour = self.ax_list[ax_index].contour(xaxis, yaxis, tbldata.T, origin=origin,
                                                         interpolation=interpolation, linestyles=linestyles,
                                                         levels=levels, cmap=self.cmap)
        else:
            if colors is not None:
                contour = self.ax_list[ax_index].contour(tbldata.T, origin=origin, interpolation=interpolation,
                                                         linestyles=linestyles, extent=extent, levels=levels,
                                                         colors=colors)
            else:
                contour = self.ax_list[ax_index].contour(tbldata.T, origin=origin, interpolation=interpolation,
                                                         linestyles=linestyles, extent=extent, levels=levels, cmap=self.cmap)

        if label_type == 'percentage':
            plt.clabel(contour, inline=1,
                       fmt=r'$\SI{%1.0f}{\percent}$', fontsize=fontsize)
        elif label_type == 'percentage_1_decimal':
            plt.clabel(contour, inline=1,
                       fmt=r'$\SI{%1.1f}{\percent}$', fontsize=fontsize)
        elif label_type == '0_decimal':
            plt.clabel(contour, inline=1,
                       fmt=r'$\SI{%1.0f}{}$', fontsize=fontsize)
        elif label_type == 'tau':
            plt.clabel(contour, inline=1,
                       fmt=r'$\tau=%1.0f$', fontsize=fontsize)
        elif label_type == 'tau_1_decimal':
            plt.clabel(contour, inline=1,
                       fmt=r'$\tau=%1.1f$', fontsize=fontsize)
        elif label_type == 'default':
            plt.clabel(contour, inline=1, fontsize=fontsize)

    def plot_surface(self, xmesh, ymesh, zmesh, ax_index=0, rstride=8, cstride=8, alpha=0.8, cmap=None,
                     color=None):
        """Plot 3D surface.

        Args:
            xmesh: Mesh grid of x-axis coordinates.
            ymesh: Mesh grid of y-axis coordinates.
            zmesh: Mesh grid of z-axis coordinates.
            ax_index (int): Index of subplot image.
            rstride (int): Array row stride (step size).
            cstride (int): Array column stride (step size).
            alpha (float): Alpha value of the colored surface.
            cmap (str): Colormap to color the surface.
            color (str): Color of the surface.
        """
        # Use default colormap if no one is defined
        if self.cmap is None:
            if cmap is not None:
                self.cmap = cmap
            else:
                self.cmap = 'viridis'

        if self.image_type != 'projection_3d':
            raise AttributeError('3D plot environment is needed!')
        if color is not None:
            self.ax_list[ax_index].plot_surface(xmesh, ymesh, zmesh, rstride=rstride, cstride=cstride, alpha=alpha,
                                                color=color)
        else:
            self.ax_list[ax_index].plot_surface(xmesh, ymesh, zmesh, rstride=rstride, cstride=cstride, alpha=alpha,
                                                cmap=self.cmap)

    def plot_quiver_3d(self, xmesh, ymesh, zmesh, u, v, w, ax_index=0, length=10, pivot='middle', color='white',
                       arrow_length_ratio=0.6):
        """Plot quiver in 3D mode.

        Args:
            xmesh: Mesh grid of x-axis coordinates.
            ymesh: Mesh grid of y-axis coordinates.
            zmesh: Mesh grid of z-axis coordinates.
            u: Vector x-axis component.
            v: Vector y-axis component.
            w: Vector z-axis component.
            ax_index (int): Index of subplot image.
            length (float): The length of each vector.
            pivot (str): Anchor of the vectors.
            color (str): Color of the vectors.
            arrow_length_ratio (float): The ratio of the arrow head with respect to the quiver.
        """
        if self.image_type != 'projection_3d':
            raise AttributeError('3D plot environment is needed!')
        self.ax_list[ax_index].quiver3D(xmesh, ymesh, zmesh, u, v, w, length=length, pivot=pivot, color=color,
                                        arrow_length_ratio=arrow_length_ratio)

    def plot_text(self, text_pos, text, ax_index=0, color='k', relative_position=False,
                  horizontalalignment='center', verticalalignment='center', zdir=(0, 0, 1), **args):
        """Plot text in the image.

        Args:
            text_pos (List[float, float]): List with the 2D text positions.
            text (str): Text that will be plotted.
            ax_index (int): Index of subplot image.
            color (str): Color of the text.
            relative_position (bool): Position of text is in relation to the axes extent.
                (x in [0, 1], y in [0, 1])
            horizontalalignment (str): Horizontal alignment of the text.
                (left, center, right)
            verticalalignment (str): Vertical alignment of the text.
                (top, center, bottom)
            zdir (List): Direction of the text object in 3D plots.
        """
        # Calculate relative positions inside the image if chosen
        if relative_position:
            text_pos[0] = text_pos[0] * \
                (self.extent[1] - self.extent[0]) + self.extent[0]
            text_pos[1] = text_pos[1] * \
                (self.extent[3] - self.extent[2]) + self.extent[2]
        # Plot text
        if self.image_type == 'projection_3d':
            self.text = self.ax_list[ax_index].text(text_pos[0], text_pos[1], text_pos[2], text, zdir, color=color,
                                                    horizontalalignment=horizontalalignment, verticalalignment=verticalalignment, **args)
        else:
            # Plot the text object
            self.text = self.ax_list[ax_index].text(text_pos[0], text_pos[1], text, color=color, horizontalalignment=horizontalalignment,
                                                    verticalalignment=verticalalignment, **args)

    def plot_rectangle(self, pos, width, height, ax_index=0, facecolor=None, edgecolor='none', alpha=0.1, hatch=None):
        """Plot rectangle in the image.

        Args:
            pos (List): Position of the rectangle (2D).
            width (float): Width of the rectangle.
            height (float): Height of the rectangle.
            ax_index (int): Index of subplot image.
            facecolor (str): Color of the rectangle surface.
            edgecolor (str): Color of the rectangle edges.
            alpha (float): Transparency of the rectangle.
            hatch (str): symbols to fill rectangle
        """
        self.ax_list[ax_index].add_patch(patches.Rectangle(pos, width, height, facecolor=facecolor,
                                                           edgecolor=edgecolor, alpha=alpha, hatch=hatch))

    def plot_circle(self, pos, size, ax_index=0, **args):
        """Plot circle in the image.

        Args:
            pos (List): Position of the circle (2D).
            size (float): Size of the circle.
            ax_index (int): Index of subplot image.
            args: Additional arguments.
        """
        self.ax_list[ax_index].add_patch(patches.Circle(pos, size, **args))

    def plot_ellipse(self, xy, width, height, angle=0, ax_index=0, **args):
        """Plot circle in the image.

        Args:
            pos (List): Position of the circle (2D).
            size (float): Size of the circle.
            ax_index (int): Index of subplot image.
            args: Additional arguments.
        """
        self.ax_list[ax_index].add_patch(patches.Ellipse(xy, width, height, angle, **args))

    def plot_wedge(self, pos, radius, theta1, theta2, width=None, ax_index=0, **args):
        """Plot wedge in the image.

        Args:
            pos (List): Position of the wedge (2D).
            radius (float): Radius of the wedge.
            theta1 (float): First angle limit of the wedge.
            theta2 (float): Second angle limit of the wedge.
            width (float): Width of the wedge.
            ax_index (int): Index of subplot image.
            args: Additional arguments.
        """
        self.ax_list[ax_index].add_patch(patches.Wedge(
            pos, radius, theta1, theta2, width, **args))

    def plot_title(self, text, ax_index=0, **args):
        """Plot title on top of an image.

        Args:
            text (str): Text that will be plotted.
            ax_index (int): Index of subplot image.
        """
        self.ax_list[ax_index].set_title(text, **args)

    def plot_arrow(self, arrow_origin, arrow_offset, ax_index=0, head_width=0.1, head_length=0.1,
                   width=0.03, color='black'):
        """Plot an arrow in the image.

        Args:
            arrow_origin (List): Origin of the arrow (x, y).
            arrow_offset (List): Offset of the arrow (dx, dy).
            ax_index (int): Index of subplot image.
            head_width (float): Width of the arrow head.
            head_length (float): Length of the arrow head.
            width (float): Width of arrow.
            color (str): Color of the arrow.
        """
        self.ax_list[ax_index].arrow(arrow_origin[0], arrow_origin[1], arrow_offset[0], arrow_offset[1],
                                     head_width=head_width, head_length=head_length, width=width, fc=color, ec=color)

    def plot_annotate(self, text, pos, text_pos, ax_index=0, color=None, ha='center', va='center'):
        """Plot an arrow in the image.

        Args:
            text (str): Text to show.
            pos (List): Position to point at (x, y).
            text_pos (List): Position of the text (dx, dy).
            ax_index (int): Index of subplot image.
            color (str): Color of the annotation.
            ha (str): Horizontal alignment.
            va (str): Vertical alignment.
        """
        self.ax_list[ax_index].annotate(text, xy=pos, xytext=text_pos, color=color,
                                        arrowprops=dict(linewidth=mpl.rcParams['lines.linewidth'],
                                                        arrowstyle="-|>", color=color, shrinkA=4),
                                        horizontalalignment=ha, verticalalignment=va)

    def plot_double_arrow(self, arrow_origin, arrow_offset, ax_index=0, color='black'):
        """Plot an arrow in the image.

        Args:
            arrow_origin (List): Origin of the arrow (x, y).
            arrow_offset (List): Offset of the arrow (dx, dy).
            ax_index (int): Index of subplot image.
            color (str): Color of the double arrow.
        """
        self.ax_list[ax_index].annotate('', xy=arrow_origin, xytext=(np.add(arrow_origin, arrow_offset)),
                                        arrowprops=dict(linewidth=mpl.rcParams['lines.linewidth'],
                                                        facecolor=color, edgecolor=color, arrowstyle='<->'))

    def plot_colorbar(self, ax_index=0, label='', extend='neither'):
        """Plot colorbar that includes the colors of an image.

        Args:
            ax_index (int): Index of subplot image.
            label: Label of the colorbar.
            extend (str): Extend of the colorbar if clipping is used.
                If vmin is larger than the smallest value ('min').
                If vmax is smaller than the largest value ('max').
                If both is true ('both').
        """
        # Plot the colorbar
        cbar = plt.colorbar(
            self.image, cax=self.ax_list[ax_index].cax, extend=extend)
        # Set label of the colorbar
        cbar.set_label(label)

    @staticmethod
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        """Use only part of a colormap.

        Args:
            cmap: Colormap instance.
            minval (float): Minimum value to extract.
            maxval (float): Maximum value to extract.
            n (int): Number of values.
        """
        new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(
                n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    def plot_legend(self, loc=0, ncol=1, fancybox=True, ax_index=0, bbox_to_anchor=None):
        """Plot a legend.

        Args:
            loc: Location index (0 -> find best position)
            ncol (int): Number of columns.
            ax_index (int): Index of subplot image.
            fancybox (bool): Use fancybox.
            bbox_to_anchor: Coordinates to achor legend to.
        """
        # Plot the legend
        if bbox_to_anchor is not None:
            self.ax_list[ax_index].legend(
                loc=loc, ncol=ncol, fancybox=fancybox, bbox_to_anchor=bbox_to_anchor)
        else:
            self.ax_list[ax_index].legend(
                loc=loc, ncol=ncol, fancybox=fancybox)

    def check_ticks(self):
        """Check if ticks are overlapping for vel_channel plot.
        """
        for ax in self.ax_list:
            xticks = ax.get_xticks()
            if xticks[0] == self.extent[0] or xticks[-1] == self.extent[1]:
                ax.set_xticks(xticks[1:-1])
            yticks = ax.get_yticks()
            if yticks[0] == self.extent[2] or yticks[-1] == self.extent[3]:
                ax.set_yticks(yticks[1:-1])

    def remove_ticks(self, ax_index=None, axis='xy'):
        """ Remove the ticks from a plot.
        """
        if ax_index is None:
            for ax in self.ax_list:
                if 'x' in axis.lower():
                    ax.set_xticks([])
                if 'y' in axis.lower():
                    ax.set_yticks([])
        else:
            if 'x' in axis.lower():
                self.ax_list[ax_index].set_xticks([])
            if 'y' in axis.lower():
                self.ax_list[ax_index].set_yticks([])

    def remove_axes(self, ax_index=None):
        """ Remove the axes of a plot.
        """
        if ax_index is None:
            for ax in self.ax_list:
                ax.set_visible(False)
        else:
            self.ax_list[ax_index].set_axis_off()

    @staticmethod
    def make_tight_layout():
        """Make plot layout tight.
        """
        plt.tight_layout()

    def set_limits(self, limits=None):
        """Limit axes of the plot.
        """
        if limits is not None:
            if len(limits) == 4:
                self.limits = limits
        if self.limits is not None:
            if self.image_type == 'projection_3d':
                for ax in self.ax_list:
                    if self.limits[0] is not None or self.limits[1] is not None:
                        ax.set_xlim3d(self.limits[0], self.limits[1])
                    if self.limits[2] is not None or self.limits[3] is not None:
                        ax.set_ylim3d(self.limits[2], self.limits[3])
                    if self.limits[4] is not None or self.limits[5] is not None:
                        ax.set_zlim3d(self.limits[4], self.limits[5])
            else:
                for ax in self.ax_list:
                    if self.limits[0] is not None or self.limits[1] is not None:
                        ax.set_xlim(self.limits[0], self.limits[1])
                    if self.limits[2] is not None or self.limits[3] is not None:
                        ax.set_ylim(self.limits[2], self.limits[3])

    def save_figure(self, file_io, crop_image=True, subplots_adjust=None):
        """Saves the figure either to pdf file or display it on the screen.

        Args:
            file_io: Instance of the file input/output class.
            crop_image (bool): Crop the image to use as less white space as possible.
            subplots_adjust (dict): Override default subplot_adjust commands.
        """
        # Apply limits to plots before saving the output file or printing the image
        self.set_limits()
        # Adjust space between subplots or custom choice
        if self.image_type == 'image':
            if subplots_adjust is not None:
                self.fig.subplots_adjust(**subplots_adjust)
            else:
                self.fig.subplots_adjust(wspace=0.05, hspace=0.05)
        # Reduce the unnecessary white space in the image
        if crop_image:
            bbox_inches = 'tight'
        else:
            # 3D plot sometimes are not properly shown if crop is True
            bbox_inches = None
        # Decide which kind of saving is used
        if file_io.plot_output == 'show':
            # Make Layout tight
            # plt.tight_layout()
            # Display image on the screen
            plt.show()
        elif file_io.plot_output == 'pdf':
            if file_io.pdf_file_instance is None:
                raise ValueError(
                    'Output file was not initialized! (custom plot?)')
            # Save image to pdf file with tight bboxes to save additional white space
            # plt.tight_layout()
            plt.savefig(file_io.pdf_file_instance, format='pdf',
                        bbox_inches=bbox_inches, dpi=300)
        elif file_io.plot_output == 'tex':
            # Save image as latex tikz command file
            file_io.add_output_image()
            plt.savefig(file_io.plot_output_filename, format='pgf')
        else:
            raise ValueError(
                'Output type for plots not known! (' + file_io.plot_output + ')')
        # Close Matplotlib figure
        plt.close()
