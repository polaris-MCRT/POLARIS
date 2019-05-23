#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
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
            plot = Plot(self.model, self.parse_args, xlabel=r'$\mathit{v}\ [\si{\kilo\metre\per\second}]$',
                        ylabel=self.file_io.get_quantity_labels(i_quantity),
                        extent=[velocity[0], velocity[-1], None, None], with_cbar=False)
            # Plot spectrum as line
            plot.plot_line(velocity, plot_data[i_quantity, :], marker='.')
            # Save figure to pdf file or print it on screen
            plot.save_figure(self.file_io)
