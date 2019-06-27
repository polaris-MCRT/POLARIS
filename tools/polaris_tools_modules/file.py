#!/usr/bin/env python
# -*- coding: utf-8 -*-

import array
import os

import numpy as np
from astropy.io import fits


class FileIO:
    """This class handles many different kinds of file input/output
    like write a plot to pdf or read data written out by POLARIS.
    """

    def __init__(self, model, server, parse_args, polaris_dir, tool_type, beam_size=None,
                 cmap_unit='arcsec', vec_field_size=None, ax_unit=None):
        """Initialisation of the input/output parameters.

        Args:
            model: Handles the model space including various
                quantities such as the density distribution.
            server: Handles the usage on a server/cluster system
            parse_args : Provides all parameters chosen
                by user when executing polaris package.
            polaris_dir (str): Path to polaris package directory.
            tool_type (str): The toolkit that executes this function.
                'plot'  : polaris-plot
                'run'   : polaris-run
                'remote': polaris-remote
                'grid'  : polaris-gen
            cmap_unit (str): Convert unit of maps into unit/cmap_unit.
                ('arcsec', 'absolute', 'px', 'nuF')
            vec_field_size (int): Number of pixel in one dimension for vector plot.
                (vector size needs to fit an integer amount into an image)
        """
        # Get model module
        self.model = model
        self.server = server
        self.parse_args = parse_args

        # Get math module
        from polaris_tools_modules.math import Math
        self.math = Math()

        # Set variables used for plotting and file read
        if tool_type is 'plot':
            if parse_args.show_plot:
                #: bool: show plot or save it to file
                self.plot_output = 'show'
            elif parse_args.tex_output:
                #: bool: create tikz latex files for output
                self.plot_output = 'tex'
            else:
                self.plot_output = 'pdf'
            #: str: Filename of the output plot file
            self.plot_output_filename = ''
            #: Instance of the pdf file to plot multiple pdf pages into one pdf file
            self.pdf_file_instance = None
            #: int: Index of the current plot image
            self.image_index = 0
            #: int: Number of quantities
            self.n_quantities_map = 8
            self.n_quantities_sed = 7
            self.n_quantities_gas = 6
            #: float: Size of the beam, if used [arcsec].
            if parse_args.beam_size is not None:
                self.beam_size = self.math.parse(
                    parse_args.beam_size, 'angle', 'arcsec')
            else:
                self.beam_size = beam_size
            #: int: Number of pixel in one dimension for vector plot
            if parse_args.vec_field_size is not None:
                self.vec_field_size = parse_args.vec_field_size
            elif vec_field_size is not None:
                self.vec_field_size = vec_field_size
            else:
                if self.parse_args.visualization_type == 'midplane':
                    self.vec_field_size = 8
                else:
                    self.vec_field_size = 10
            #: bool: Convert unit of maps into unit/arcsec instead of unit/pixel
            if parse_args.cmap_unit is not None:
                self.cmap_unit = parse_args.cmap_unit
            elif self.beam_size is not None:
                self.cmap_unit = 'beam'
            else:
                self.cmap_unit = 'arcsec'
            if parse_args.ax_unit is not None:
                self.ax_unit = parse_args.ax_unit
            elif ax_unit is None:
                self.ax_unit = 'au'
            else:
                self.ax_unit = ax_unit

        #: dict: Dictionary with all paths for POLARIS
        self.path = {}

        #: List: Main simulation types to choose from
        self.main_simu_types = ['temp', 'rat', 'dust_mc',
                                'dust_full', 'dust', 'line', 'zeeman']

        #: List: alignment methods to choose from
        self.alignment_mechanisms = ['pa', 'rat', 'gold', 'idg', 'internal']

        #: str: Path to polaris home directory.
        self.polaris_dir = polaris_dir

        if tool_type in ['plot', 'run', 'remote']:
            simulation_name = parse_args.simulation_name
        else:
            simulation_name = None
        if tool_type in ['plot', 'run']:
            simulation_type = parse_args.simulation_type
        else:
            simulation_type = None

        # Set path from user input.
        self.set_path_from_str(tool_type, model_name=parse_args.model_name, simulation_name=simulation_name,
                               simulation_type=simulation_type)

    def set_path_from_str(self, tool_type, model_name, simulation_name, simulation_type):
        """Sets all paths used by a given toolkit depending on input strings.

        Args:
            tool_type (str): The toolkit that executes this function.
                'grid'   : polaris-gen
                'run'    : polaris-run
                'plot'   : polaris-plot
                'remote' : polaris-remote
                'extra'  : polaris-extra
                'test'   : polaris-test
            model_name (str): Name of the model (see model.py).
            simulation_name (str): Name of the simulation (see polaris-run.in).
            simulation_type (str): Type of the simulation (see polaris-run.in).
        """
        # Path to directory with the polaris package
        self.path['polaris'] = self.polaris_dir
        # Path to directory with the polaris and PolarisTools binaries
        self.path['bin'] = self.path['polaris'] + 'bin/'
        # Path to directory with input files besides the command file
        self.path['input'] = self.path['polaris'] + 'input/'
        # Path to directory with input files besides the command file
        self.path['dust'] = self.path['input'] + 'dust/'
        # Path to projects directory
        self.path['projects'] = self.path['polaris'] + 'projects/'
        # Path to directory with the PolarisTools source files and test cases
        self.path['tools'] = self.path['polaris'] + 'tools/'
        # Path to chosen model directory
        self.path['model'] = self.path['projects'] + str(model_name) + '/'
        if tool_type in ['run', 'gen']:
            # Create model directory, if it does not exist
            if not os.path.lexists(self.path['model']):
                os.mkdir(self.path['model'])
        if tool_type in ['plot', 'run', 'remote']:
            # Path to the directory of the chosen simulation name
            self.path['simulation'] = self.path['model'] + \
                str(simulation_name) + '/'
            if tool_type == 'run':
                # Create simulation directory, if it does not exist
                if not os.path.lexists(self.path['simulation']):
                    os.mkdir(self.path['simulation'])
                # Update simulation_directory relative to the server path
                self.server.parameter['simulation_directory'] = self.path['simulation']
        if tool_type in ['plot', 'run', 'remote', 'extra']:
            # Path to the directory with the gas files (LEIDEN database)
            self.path['gas'] = self.path['input'] + 'gas/'
            # Path to the directory with the gas files (JPL database)
            self.path['jpl'] = self.path['gas'] + 'jpl_catalog/'
            # Path to the directory with the gas files (CDMS database)
            self.path['cdms'] = self.path['gas'] + 'cdms_catalog/'
            # Path to the directory with the dustem input files
            self.path['dustem'] = self.path['dust'] + 'dustem/'
            # Path to the directory with the dustem input files
            self.path['trust'] = self.path['dust'] + 'trust/DustModel/'
        if tool_type in ['plot', 'run']:
            # Path to the directory of the chosen simulation type (heat, rat, ...)
            self.path['simulation_type'] = self.path['simulation'] + \
                self.get_updated_simulation_type(simulation_type) + '/'
            # Create simulation directory, if it does not exist
            if not os.path.lexists(self.path['simulation_type']):
                os.mkdir(self.path['simulation_type'])
            # Path to the directory in which polaris package saves the plots
            self.path['plots'] = self.path['simulation_type'] + 'plots/'
            # Create simulation directory, if it does not exist
            if not os.path.lexists(self.path['plots']):
                os.mkdir(self.path['plots'])
            # Path to the directory with the POLARIS results
            self.path['results'] = self.path['simulation_type'] + 'data/'
            # Create simulation directory, if it does not exist
            if not os.path.lexists(self.path['results']):
                os.mkdir(self.path['results'])
        elif tool_type in ['remote']:
            self.server.check_for_remote()
            # Server address to push/pull files to/from
            self.path['server'] = self.server.parameter['user_id'] + \
                '@' + self.server.parameter['address']
            # Path to polaris package on the server
            self.path['server_polaris'] = self.path['server'] + \
                self.server.parameter['server_polaris_dir']
            # Path to the POLARIS model on the server
            self.path['server_model'] = self.path['model'].replace(
                self.path['polaris'], self.path['server_polaris'])
            # Path to the POLARIS results on the server
            self.path['server_simulation'] = self.path['simulation'].replace(
                self.path['polaris'], self.path['server_polaris'])
        elif tool_type in ['test']:
            # Path to the test cases used by polaris-test
            self.path['testing'] = self.path['polaris'] + 'testing/'

    # -----------------------------------------------------
    # -----------------------------------------------------
    # The following functions are reading or writing data.
    # -----------------------------------------------------
    # -----------------------------------------------------

    def init_plot_output(self, filename, path=None):
        """Creates a pdf file object to save plots into.

        Args:
            filename (str): Filename of the PDF file that will be created.
            path (str): Path to pdf file. Default plots path is used,
                if this is not set.
        """
        if self.plot_output == 'tex':
            # If path is not set, use default value
            if path is None:
                path = self.path['plots']
            self.plot_output_filename = path + filename + \
                '_image_' + str(self.image_index) + '.tex'
        elif self.plot_output == 'pdf':
            # Load necessary modules
            from matplotlib.backends.backend_pdf import PdfPages
            # If path is not set, use default value
            if path is None:
                path = self.path['plots']
            self.plot_output_filename = path + filename + '.pdf'
            # Create PDF file. Images are saved to this file
            self.pdf_file_instance = PdfPages(self.plot_output_filename)
        elif self.plot_output != 'show':
            raise ValueError('Output option ' +
                             self.plot_output + ' is not known!')

    def close_plot_output(self):
        if self.plot_output == 'pdf' and self.pdf_file_instance is not None:
            self.pdf_file_instance.close()

    def add_output_image(self):
        """Rename output_filename to allow saving multiple pdf pages instead as multiple .tex files.
        """
        self.image_index += 1
        file_type = self.plot_output_filename[-4:]
        self.plot_output_filename = self.plot_output_filename.replace(
            '_image_' + str(self.image_index - 1) + file_type, '_image_' + str(self.image_index) + file_type)

    def read_dust_file(self, dust_parameter_dict):
        """Reads dust grain sizes and wavelengths from POLARIS dust catalog files.

        Args:
            dust_parameter_dict (dict): Dictionary of the dust parameters.

        Returns:
            Lists with the dust grain sizes and wavelengths.
        """
        if dust_parameter_dict['dust_cat_file'] is None or '.nk' in dust_parameter_dict['dust_cat_file']:
            raise ValueError('The chosen dust component has no related dust catalog '
                             'or is a refractive index file!')
        size_list = []
        wavelength_list = []
        # Get dust data from catalog file
        dust_file = open(self.path['dust'] +
                         dust_parameter_dict['dust_cat_file'], 'r')
        i_line = 0
        i_print_dust_size = -1
        i_print_wavelength = -1
        for line in dust_file.readlines():
            if '#a_eff' == line.strip():
                i_print_dust_size = i_line + 2
            elif '#wavelength' == line.strip():
                i_print_wavelength = i_line + 2
            elif i_line == i_print_dust_size:
                for data in line.split():
                    size_list.append(float(data))
            elif i_line == i_print_wavelength:
                for data in line.split():
                    wavelength_list.append(float(data))
            i_line += 1
        dust_file.close()
        return size_list, wavelength_list

    def get_quantity_labels(self, i_quantity, int_map=False):
        """Creates automatically colorbar labels.

        Args:
            i_quantity (int): Index of the quantity which will be shown on the colorbar.
            int_map (bool): Is the colorbar for an integrated velocity channel map?

        Returns:
            str: Predefined labels for the colorbar.
        """

        # If beam_label is True or a beam convolution is performed, change to Jy/beam
        quantity = r'\mathit{F}'
        if self.cmap_unit == 'arcsec':
            unit = 'Jy/as^2'
        elif self.cmap_unit == 'px':
            unit = 'Jy/px'
        elif self.cmap_unit == 'total':
            unit = 'Jy'
        elif self.cmap_unit == 'nuF':
            unit = r'\watt\per\metre\squared'
            quantity = r'\mathit{\nu F}'
        elif self.cmap_unit == 'beam':
            unit = 'Jy/beam'
        else:
            self.cmap_unit = 'arcsec'
            unit = 'Jy/as^2'
        # Integrated velocity maps have their unit multiplied by velocity
        if int_map:
            unit += r'\cdot\kilo\metre\per\second'
        if self.parse_args.simulation_type == 'dust_mc':
            #: list: List with all available colorbar labels for monte_carlo plots
            quantity_labels = [r'$' + quantity + r'_\mathsf{I}\ [\si{' + unit + '}]$',
                               r'$' + quantity +
                               r'_\mathsf{Q}\ [\si{' + unit + '}]$',
                               r'$' + quantity +
                               r'_\mathsf{U}\ [\si{' + unit + '}]$',
                               r'$' + quantity +
                               r'_\mathsf{V}\ [\si{' + unit + '}]$',
                               r'$\mathit{PI}\ [\mathsf{' + unit + '}]$',
                               r'$\mathit{P}_\mathsf{l}\ [\%]$',
                               r'$\mathit{I}_\mathsf{direct}\ [\si{' + unit + '}]$',
                               r'$\mathit{I}_\mathsf{scattered}\ [\si{' + unit + '}]$']
        elif 'dust' in self.parse_args.simulation_type:
            #: list: List with all available colorbar labels for raytrace plots
            quantity_labels = [r'$' + quantity + r'_\mathsf{I}\ [\si{' + unit + '}]$',
                               r'$' + quantity +
                               r'_\mathsf{Q}\ [\si{' + unit + '}]$',
                               r'$' + quantity +
                               r'_\mathsf{U}\ [\si{' + unit + '}]$',
                               r'$' + quantity +
                               r'_\mathsf{V}\ [\si{' + unit + '}]$',
                               r'$\mathit{PI}\ [\mathsf{' + unit + '}]$',
                               r'$\mathit{P}_\mathsf{l}\ [\%]$',
                               r'$\tau$',
                               r'$\mathit{n}_\mathsf{H}\ [\si{m^{-2}}]$']
        elif self.parse_args.simulation_type in ['line', 'zeeman']:
            #: list: List with all available colorbar labels for spectral line plots
            quantity_labels = [r'$' + quantity + r'_\mathsf{I}\ [\si{' + unit + '}]$',
                               r'$' + quantity +
                               r'_\mathsf{Q}\ [\si{' + unit + '}]$',
                               r'$' + quantity +
                               r'_\mathsf{U}\ [\si{' + unit + '}]$',
                               r'$' + quantity +
                               r'_\mathsf{V}\ [\si{' + unit + '}]$',
                               r'$\tau$',
                               r'$\mathit{n}_\mathsf{H}\ [\si{m^{-2}}]$']
        else:
            raise ValueError('Error: For the simulation_type: ' + str(self.parse_args.simulation_type) +
                             ' is no automatic label generation for the colorbar available!')
        if len(quantity_labels) <= i_quantity or i_quantity < 0:
            raise ValueError('For the quantity number ' +
                             str(i_quantity) + ', no label is providable!')
        return quantity_labels[i_quantity]

    def read_emission_map_header(self, hdulist):
        """Reads the header of a fits file with POLARIS Raytrace or Monte-Carlo map results.

        Args:
            hdulist: Fits HDU list object.

        Returns:
            dict: Dictionary with the information in the header.
            str: Description of the data contained in the fits (healpix, map) or
                chosen postprocessing of the data (cut, radial)
        """
        # Check type of Polaris fits file data
        plot_data_type = self.check_fits_data(hdulist)
        if plot_data_type == 'healpix':
            hdu_index = 1
        else:
            hdu_index = 0
        #: dict: Dictionary with parameters from header
        header_dict = dict(
            wavelengths=[],
            ID=hdulist[hdu_index].header['ID'],
            emission_type=hdulist[hdu_index].header['ETYPE'],
            simulation_type='dust' if 'thermal' in hdulist[hdu_index].header['ETYPE'] else 'dust_mc',
        )
        if plot_data_type == 'healpix':
            header_dict['nsides'] = int(hdulist[hdu_index].header['NSIDE'])
            header_dict['nr_pixel_x'] = int(
                hdulist[hdu_index].header['LASTPIX'] - hdulist[hdu_index].header['FIRSTPIX'] + 1)
            header_dict['nr_wavelengths'] = int(
                hdulist[hdu_index].header['TFIELDS'] / 6)
            # Put the wavelengths in list
            for i_wl in range(header_dict['nr_wavelengths']):
                header_dict['wavelengths'].append(float(
                    hdulist[hdu_index].header['TTYPE' + str(6 * i_wl + 1)].split()[-2]))
        else:
            # Update dictionary with parameters from header
            header_dict['angle1'] = hdulist[hdu_index].header['RANGLE1']
            header_dict['angle2'] = hdulist[hdu_index].header['RANGLE2']
            header_dict['sidelength_x'] = hdulist[hdu_index].header['CDELT1'] * \
                hdulist[hdu_index].header['NAXIS1']
            header_dict['sidelength_y'] = hdulist[hdu_index].header['CDELT2'] * \
                hdulist[hdu_index].header['NAXIS2']
            header_dict['nr_pixel_x'] = hdulist[hdu_index].header['NAXIS1']
            header_dict['nr_pixel_y'] = hdulist[hdu_index].header['NAXIS2']
            header_dict['nr_wavelengths'] = hdulist[hdu_index].header['NAXIS3']
            header_dict['distance'] = hdulist[hdu_index].header['DISTANCE']
            # Adjust the unit of the sidelengths
            header_dict['sidelength_x'] *= self.math.conv_length_factor(hdulist[0].header['CUNIT1'].lower(),
                                                                        header_dict['distance'])
            header_dict['sidelength_y'] *= self.math.conv_length_factor(hdulist[0].header['CUNIT2'].lower(),
                                                                        header_dict['distance'])
            # Put the wavelengths in list
            for i_wl in range(header_dict['nr_wavelengths']):
                header_dict['wavelengths'].append(
                    hdulist[0].header['WAVELENGTH' + str(i_wl + 1)])
            # Check if the sidelength of the image ist the same as of the model
            self.check_sidelengths(header_dict)
            # Check if the distance in the fits file is the same as in the model
            self.check_distance(header_dict)
        return header_dict, plot_data_type

    def read_emission_sed_header(self, hdulist):
        """Reads the header of a fits file with POLARIS Raytrace or Monte-Carlo SED results.

        Args:
            hdulist: Fits HDU list object.

        Returns:
            dict: Dictionary with the information in the header.
            str: Description of the data contained in the fits (healpix, map) or
                chosen postprocessing of the data (cut, radial)
        """
        # Check type of Polaris fits file data
        plot_data_type = self.check_fits_data(hdulist)
        #: dict: Dictionary with parameters from header
        header_dict = dict(
            angle1=hdulist[0].header['RANGLE1'],
            angle2=hdulist[0].header['RANGLE2'],
            nr_wavelengths=hdulist[0].header['NAXIS1'],
            wavelengths=[],
            distance=hdulist[0].header['DISTANCE'],
            ID=hdulist[0].header['ID'],
            simulation_type='dust' if 'thermal' in hdulist[0].header['ETYPE'] else 'dust_mc',
        )
        # Put the wavelengths in list
        for i_wl in range(header_dict['nr_wavelengths']):
            header_dict['wavelengths'].append(
                hdulist[0].header['WAVELENGTH' + str(i_wl + 1)])
        # Check if the distance in the fits file is the same as in the model
        self.check_distance(header_dict)
        return header_dict, plot_data_type

    def read_midplane_file_header(self, hdulist):
        """Reads the header of a file with POLARIS midplane results.

        Args:
            hdulist: Fits HDU list object.

        Returns:
            dict: Dictionary with the information in the header.
            str: Description of the data contained in the fits (healpix, map) or
                chosen postprocessing of the data (cut, radial)
        """
        # Check type of Polaris fits file data
        plot_data_type = self.check_fits_data(hdulist)
        #: dict: Dictionary with parameters from header
        header_dict = dict(
            nr_pixel_x=hdulist[0].header['NAXIS1'],
            nr_pixel_y=hdulist[0].header['NAXIS2'],
            sidelength_x=hdulist[0].header['CDELT1'] *
            hdulist[0].header['NAXIS1'],
            sidelength_y=hdulist[0].header['CDELT2'] *
            hdulist[0].header['NAXIS2'],
            midplane_quantities=[],
        )
        # Put the midplane quantities in list
        for i in range(hdulist[0].header['NAXIS4']):
            header_dict['midplane_quantities'].append(
                hdulist[0].header['MIDPLANE' + str(i + 1)])
        # Adjust the unit of the sidelengths
        header_dict['sidelength_x'] *= self.math.conv_length_factor(
            hdulist[0].header['CUNIT1'].lower())
        header_dict['sidelength_y'] *= self.math.conv_length_factor(
            hdulist[0].header['CUNIT2'].lower())
        # Check if the sidelength of the image ist the same as of the model
        self.check_sidelengths(header_dict)
        return header_dict, plot_data_type

    def read_vel_maps_header(self, hdulist):
        """Reads the header from velocity channel map fits files.

        Args:
            hdulist: Fits HDU list object.

         Returns:
            dict: Dictionary with the information in the header.
            str: Description of the data contained in the fits (healpix, map) or
                chosen postprocessing of the data (cut, radial)
        """
        # Check type of Polaris fits file data
        plot_data_type = self.check_fits_data(hdulist)
        if plot_data_type == 'healpix':
            hdu_index = 1
        else:
            hdu_index = 0
        #: dict: Dictionary with parameters from header
        header_dict = dict(
            species_name=hdulist[hdu_index].header['GAS_SPECIES'],
            i_transition=hdulist[hdu_index].header['TRANS'],
            frequency=hdulist[hdu_index].header['FREQ'],
            nr_channels=hdulist[hdu_index].header['CHANNELS'],
            max_velocity=hdulist[hdu_index].header['MAXVEL'],
        )
        if plot_data_type == 'healpix':
            header_dict['nsides'] = int(hdulist[hdu_index].header['NSIDE'])
            header_dict['nr_pixel_x'] = int(
                hdulist[hdu_index].header['LASTPIX'] - hdulist[hdu_index].header['FIRSTPIX'] + 1)
            header_dict['nr_columns'] = int(
                hdulist[hdu_index].header['TFIELDS'])
        else:
            # Update dictionary with parameters from header
            header_dict['angle1'] = hdulist[hdu_index].header['RANGLE1']
            header_dict['angle2'] = hdulist[hdu_index].header['RANGLE2']
            header_dict['sidelength_x'] = hdulist[hdu_index].header['CDELT1'] * \
                hdulist[hdu_index].header['NAXIS1']
            header_dict['sidelength_y'] = hdulist[hdu_index].header['CDELT2'] * \
                hdulist[hdu_index].header['NAXIS2']
            header_dict['nr_pixel_x'] = hdulist[hdu_index].header['NAXIS1']
            header_dict['nr_pixel_y'] = hdulist[hdu_index].header['NAXIS2']
            header_dict['distance'] = hdulist[hdu_index].header['DISTANCE']
            # Adjust the unit of the sidelengths
            header_dict['sidelength_x'] *= self.math.conv_length_factor(hdulist[0].header['CUNIT1'].lower(),
                                                                        header_dict['distance'])
            header_dict['sidelength_y'] *= self.math.conv_length_factor(hdulist[0].header['CUNIT2'].lower(),
                                                                        header_dict['distance'])
            # Check if the sidelength of the image ist the same as of the model
            self.check_sidelengths(header_dict)
            # Check if the distance in the fits file is the same as in the model
            self.check_distance(header_dict)
        return header_dict, plot_data_type

    def read_spectrum_header(self, hdulist):
        """Reads the header from spectrum fits files.

        Args:
            hdulist: Fits HDU list object.

         Returns:
            dict: Dictionary with the information in the header.
            str: Description of the data contained in the fits (healpix, map) or
                chosen postprocessing of the data (cut, radial)
        """
        #: dict: Dictionary with parameters from header
        header_dict = dict(
            species_name=hdulist[0].header['GAS_SPECIES'],
            i_transition=hdulist[0].header['TRANS'],
            frequency=hdulist[0].header['FREQ'],
            nr_channels=hdulist[0].header['CHANNELS'],
            max_velocity=hdulist[0].header['MAXVEL'],
        )
        return header_dict

    def read_emission_map(self, filename):
        """Reads the data of a fits file with POLARIS Raytrace or Monte-Carlo map results.

        Args:
            filename (str): Filename to the text file with Raytrace or Monte-Carlo results.

        Returns:
            Tuple with:
            1) Numpy array with 2 dimensions for the pixel coordinates
            and one dimension for the Stokes and polarization data.
            2) Dictionary with information of the file header.
        """
        # Load data and header from fits if available
        if os.path.isfile(self.path['results'] + filename + '.fits'):
            hdulist = fits.open(self.path['results'] + filename + '.fits')
        else:
            raise FileExistsError('--- Hint: No emission map fits file exists (' +
                                  self.path['results'] + filename + '.fits' + ')')
        #: dict: Dictionary with the information in the header.
        header, plot_data_type = self.read_emission_map_header(hdulist)
        if plot_data_type == 'healpix':
            if header['simulation_type'] == 'dust_full':
                raise ValueError('Emission maps with the healpix background raytracing grid cannot be '
                                 'combined with Monte-Carlo results!')
            import healpy as hp
            data = hp.read_map(
                self.path['results'] + filename + '.fits', field=None, verbose=False)
            wmap_map = np.zeros(
                (self.n_quantities_map, header['nr_wavelengths'], data.shape[1]))
            #: Amount of arcseconds per pixel squared to convert flux from Jy/pixel into Jy/arcsec^2
            arcsec_squared_per_pixel = 4 * np.pi * \
                (180 / np.pi * 60 * 60) ** 2 / header['nr_pixel_x']
            for i_wl in range(header['nr_wavelengths']):
                i_col = i_wl * int(self.n_quantities_map - 3)
                # Set the data array from fits input
                wmap_map[0:4, i_wl, :] = data[i_col:i_col + 4, :]
                wmap_map[6:7, i_wl, :] = data[i_col + 4: i_col + 5, :]
                wmap_map[4, i_wl, :] = np.sqrt(np.add(np.power(data[i_col + 1, :], 2),
                                                      np.power(data[i_col + 2, :], 2)))
            # Apply intensity unit conversion
            if self.cmap_unit == 'total' or self.cmap_unit == 'px':
                pass
            elif self.cmap_unit == 'nuF':
                for i_wl in range(header['nr_wavelengths']):
                    wmap_map[0:5, i_wl, :] *= 1e-26 * \
                        self.math.const['c'] / header['wavelengths'][i_wl]
            else:
                wmap_map[0:5, :, :] /= arcsec_squared_per_pixel
            # Add column density at the end
            wmap_map[7, :, :] = data[int(
                self.n_quantities_map - 3) * header['nr_wavelengths'], :]
            # Apply beam convolution if chosen
            wmap_map = self.beam_conv(
                wmap_map, plot_data_type, np.sqrt(arcsec_squared_per_pixel))
            # Calculate degree of polarization
            wmap_map[5, i_wl, :] = 1e2 * np.divide(wmap_map[4, i_wl, :], wmap_map[0, i_wl, :],
                                                   out=np.zeros_like(wmap_map[4, i_wl, :]), where=wmap_map[0, i_wl, :] != 0)
            return wmap_map, header, plot_data_type
        else:
            #: Data from the file with the monte carlo results
            data = np.transpose(hdulist[0].data, (0, 1, 3, 2))
            #: 3 dimensional numpy array to save the results into
            tbldata = np.zeros((self.n_quantities_map, header['nr_wavelengths'],
                                header['nr_pixel_x'], header['nr_pixel_y']))
            #: Amount of arcseconds per pixel to convert flux from Jy/pixel into Jy/arcsec^2
            arcsec_squared_per_pixel = (2. * self.model.tmp_parameter['radius_x_arcsec'] / header['nr_pixel_x']) * \
                (2. *
                 self.model.tmp_parameter['radius_y_arcsec'] / header['nr_pixel_y'])
            # Set the data array from fits input
            tbldata[0:4, :, :, :] = data[0:4, :, :, :]
            tbldata[6:8, :, :, :] = data[4:6, :, :, :]
            tbldata[4, :, :, :] = np.sqrt(
                np.add(np.power(data[1, :, :, :], 2), np.power(data[2, :, :, :], 2)))
            # Apply intensity unit conversion
            if self.cmap_unit == 'total' or self.cmap_unit == 'px':
                pass
            elif self.cmap_unit == 'nuF':
                for i_wl in range(header['nr_wavelengths']):
                    tbldata[0:5, i_wl, :, :] *= 1e-26 * \
                        self.math.const['c'] / header['wavelengths'][i_wl]
                if header['simulation_type'] in ['dust_mc', 'dust_full']:
                    for i_wl in range(header['nr_wavelengths']):
                        tbldata[6:8, i_wl, :, :] *= 1e-26 * \
                            self.math.const['c'] / header['wavelengths'][i_wl]
            else:
                tbldata[0:5, :, :, :] /= arcsec_squared_per_pixel
                if header['simulation_type'] in ['dust_mc', 'dust_full']:
                    tbldata[6:8, :, :, :] /= arcsec_squared_per_pixel
            # Apply beam convolution if chosen
            tbldata = self.beam_conv(
                tbldata, plot_data_type, np.sqrt(arcsec_squared_per_pixel))
            # Calculate degree of polarization
            tbldata[5, :, :, :] = 1e2 * np.divide(tbldata[4, :, :, :], tbldata[0, :, :, :],
                                                  out=np.zeros_like(tbldata[4, :, :, :]), where=tbldata[0, :, :, :] != 0)
            if plot_data_type == 'map':
                return tbldata, header, plot_data_type
            elif plot_data_type == 'cut':
                # Create cuts instead of returning the map data
                position, data = self.create_cut(tbldata, cut_parameter=self.parse_args.cut_parameter,
                                                 N_r=2 * (header['nr_pixel_x'] + header['nr_pixel_y']))
                return [position, data], header, plot_data_type
            elif plot_data_type == 'radial':
                # Create cuts instead of returning the map data
                position, data = self.create_radial_profile(tbldata,
                                                            radial_parameter=self.parse_args.radial_parameter)
                return [position, data], header, plot_data_type
            else:
                raise ValueError(
                    'The emission map data cannot be delivered in the chosen format!')

    def read_emission_sed(self, filename):
        """Reads the data of a fits file with POLARIS Raytrace or Monte-Carlo SED results.

        Args:
            filename (str): Filename to the text file with Raytrace or Monte-Carlo results.

        Returns:
            Tuple with:
            1) Numpy array with 2 dimensions for the pixel coordinates
            and one dimension for the Stokes and polarization data.
            2) Dictionary with information of the file header.
        """
        # Load data and header from fits if available
        if os.path.isfile(self.path['results'] + filename + '.fits'):
            hdulist = fits.open(self.path['results'] + filename + '.fits')
        else:
            raise FileExistsError('--- Hint: No emission sed fits file exists (' +
                                  self.path['results'] + filename + '.fits)')
        #: dict: Dictionary with the information in the header.
        header, plot_data_type = self.read_emission_sed_header(hdulist)
        #: Data from the file with the monte carlo results
        data = hdulist[0].data[:, 0, :]
        #: 3 dimensional numpy array to save the results into
        tbldata = np.zeros((self.n_quantities_sed, header['nr_wavelengths']))
        # Set the data array from fits input
        tbldata[0:4, :] = data[0:4, :]
        tbldata[6, :] = data[4, :]
        tbldata[4, :] = np.sqrt(
            np.add(np.power(data[1, :], 2), np.power(data[2, :], 2)))
        # Ignore divide by zero
        with np.errstate(divide='ignore', invalid='ignore'):
            tbldata[5, :] = 1e2 * np.divide(tbldata[4, :], data[0, :])
        # Apply intensity unit conversion
        if self.cmap_unit == 'nuF':
            for i_wl in range(header['nr_wavelengths']):
                tbldata[0:5, i_wl] *= 1e-26 * \
                    self.math.const['c'] / header['wavelengths'][i_wl]
        return tbldata, header, plot_data_type

    def read_midplane_file(self, visualization_input, skip_not_known=False, filename=None):
        """Reads the data of a file with midplane results.

        Args:
            visualization_input (str): Identifier of the midplane data.
                (e.g. input_gas_density_xy)
            skip_not_known (bool): Ignore if midplane does not exists (if 'all' was chosen).
            filename (str): Alternative file in the results directory to use.

        Returns:
            Tuple with:
            1) numpy array with 2 dimensions for the pixel coordinates
            and one dimension for the vector data.
            2) numpy array with 2 dimensions for the pixel coordinates.
            The values are plotted as colors.
            3) Dictionary with information of the file header.
        """
        # Load data and header from fits if available
        from astropy.io import fits
        if filename is not None:
            if os.path.isfile(self.path['results'] + filename):
                hdulist = fits.open(self.path['results'] + filename)
            else:
                raise FileExistsError('Error: No midplane file exists!')
        elif 'input' in visualization_input:
            if os.path.isfile(self.path['results'] + 'input_midplane.fits'):
                hdulist = fits.open(
                    self.path['results'] + 'input_midplane.fits')
            else:
                raise FileExistsError('Error: No midplane file exists!')
        elif 'output' in visualization_input:
            if os.path.isfile(self.path['results'] + 'output_midplane.fits'):
                hdulist = fits.open(
                    self.path['results'] + 'output_midplane.fits')
            else:
                raise FileExistsError('Error: No midplane file exists!')
        else:
            raise ValueError(
                'Error: visualization_input is not set correctly (try \'all\')!')
        #: dict: Dictionary with the information in the header.
        header, plot_data_type = self.read_midplane_file_header(hdulist)
        #: int: Number of bins per dimension per polarization vector
        if header['nr_pixel_x'] % 2 != 0 and self.vec_field_size % 2 == 0:
            offset_x = (header['nr_pixel_x'] - 1) % self.vec_field_size
        else:
            offset_x = header['nr_pixel_x'] % self.vec_field_size
        while offset_x % 2 != 0:
            offset_x += self.vec_field_size
        vector_bins_x = int(
            (header['nr_pixel_x'] - offset_x) / self.vec_field_size)
        offset_x = int(offset_x / 2)

        if header['nr_pixel_y'] % 2 != 0 and self.vec_field_size % 2 == 0:
            offset_y = (header['nr_pixel_y'] - 1) % self.vec_field_size
        else:
            offset_y = header['nr_pixel_y'] % self.vec_field_size
        while offset_y % 2 != 0:
            offset_y += self.vec_field_size
        vector_bins_y = int(
            (header['nr_pixel_y'] - offset_y) / self.vec_field_size)
        offset_y = int(offset_y / 2)
        #: Numpy array for the vector data
        vec_field_data = np.zeros((vector_bins_x, vector_bins_y, 5))
        # Get index for different cuts through the model
        if 'xy' in visualization_input:
            cut_index = 0
        elif 'xz' in visualization_input:
            cut_index = 1
        elif 'yz' in visualization_input:
            cut_index = 2
        else:
            raise ValueError('Midplane cut is not one of xy, xz, yz!')
        #: float: potential multiplier to modify midplane data
        midplane_mod = 1.
        # Get index for different cuts through the model
        try:
            midlane_index = -1
            if 'gas_number_density' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'gas_number_density ' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'gas_mass_density' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'gas_mass_density ' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'dust_number_density' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'dust_number_density ' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
                    elif 'gas_number_density ' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        midplane_mod = self.model.parameter['mass_fraction']
                        break
            elif 'dust_mass_density' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'dust_mass_density ' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
                    elif 'gas_mass_density ' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        midplane_mod = self.model.parameter['mass_fraction']
                        break
            elif 'gas_temperature' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'gas_temperature ' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'dust_temperature' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'dust_temperature ' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'rat_aalig' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'rat_aalig' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'delta' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'delta' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'mag' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'mag_total' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'vel' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'vel_total' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'mach' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'mach' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'larm' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'larm' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'dust_choice' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'dust_choice' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'a_min' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'a_min' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            elif 'a_max' in visualization_input:
                for quantity in header['midplane_quantities']:
                    if 'a_max' in quantity:
                        midlane_index = header['midplane_quantities'].index(
                            quantity)
                        break
            if midlane_index == -1:
                raise ValueError('visualization_input: ' +
                                 visualization_input + ' not found!')
        except:
            if skip_not_known:
                return None, None, plot_data_type
            else:
                raise ValueError('Midplane file does not include the chosen quantity! '
                                 'Choose another quantity or try \'all\'.')
        data = np.transpose(hdulist[0].data, (0, 1, 3, 2))
        # Init steps to jump over the middle row/column is necessary
        middle_step_x = 0
        middle_step_y = 0
        if 'vel' in visualization_input or 'mag' in visualization_input:
            for i_x in range(offset_x, header['nr_pixel_x'] - offset_x):
                # Skip middle cell if vector size is even and pixel size is odd
                if header['nr_pixel_x'] % 2 != 0 and self.vec_field_size % 2 == 0 \
                        and i_x == int((header['nr_pixel_x'] - 1) / 2):
                    if middle_step_x == 0:
                        middle_step_x = 1
                    continue
                # Next axis loop
                for i_y in range(offset_y, header['nr_pixel_y'] - offset_y):
                    # Skip middle cell if vector size is even and pixel size is odd
                    if header['nr_pixel_y'] % 2 != 0 and self.vec_field_size % 2 == 0 \
                            and i_y == int((header['nr_pixel_y'] - 1) / 2):
                        if middle_step_y == 0:
                            middle_step_y = 1
                        continue
                    #: int: X-axis index to the vector
                    i_x_vec = int(
                        (i_x - offset_x - middle_step_x) / self.vec_field_size)
                    #: int: Y-axis index to the vector
                    i_y_vec = int(
                        (i_y - offset_y - middle_step_y) / self.vec_field_size)
                    # Calculating the vectors and average over vec_field_size^2
                    vec_field_data[i_x_vec, i_y_vec, 0] += midplane_mod * \
                        data[midlane_index + 1, cut_index, i_x, i_y] / \
                        self.vec_field_size ** 2
                    vec_field_data[i_x_vec, i_y_vec, 1] += midplane_mod * \
                        data[midlane_index + 2, cut_index, i_x, i_y] / \
                        self.vec_field_size ** 2
                    vec_field_data[i_x_vec, i_y_vec, 2] += midplane_mod * \
                        data[midlane_index + 3, cut_index, i_x, i_y] / \
                        self.vec_field_size ** 2
        # Calculating the x and y positions
        for i_x_vec in range(vector_bins_x):
            for i_y_vec in range(vector_bins_y):
                # X axis
                vec_field_data[i_x_vec, i_y_vec, 3] = ((
                    i_x_vec + 0.5) * self.vec_field_size + offset_x) / header['nr_pixel_x']
                if header['nr_pixel_x'] % 2 != 0 and self.vec_field_size % 2 == 0:
                    if i_x_vec == int(vector_bins_x / 2):
                        vec_field_data[i_x_vec, i_y_vec,
                                       3] += 0.5 / header['nr_pixel_x']
                    elif i_x_vec > int(vector_bins_x / 2):
                        vec_field_data[i_x_vec, i_y_vec,
                                       3] += 1 / header['nr_pixel_x']
                # Y axis
                vec_field_data[i_x_vec, i_y_vec, 4] = ((
                    i_y_vec + 0.5) * self.vec_field_size + offset_y) / header['nr_pixel_y']
                if header['nr_pixel_y'] % 2 != 0 and self.vec_field_size % 2 == 0:
                    if i_y_vec == int(vector_bins_y / 2):
                        vec_field_data[i_x_vec, i_y_vec,
                                       4] += 1 / header['nr_pixel_y']
                    if i_y_vec > int(vector_bins_y / 2):
                        vec_field_data[i_x_vec, i_y_vec,
                                       4] += 1 / header['nr_pixel_y']
        #: Numpy array for the midplane data
        tbldata = hdulist[0].data[midlane_index,
                                  cut_index, :, :].T * midplane_mod
        if plot_data_type == 'map':
            return [tbldata, vec_field_data], header, plot_data_type
        elif plot_data_type == 'cut':
            # Create cuts instead of returning the map data
            position, data = self.create_cut(tbldata, cut_parameter=self.parse_args.cut_parameter,
                                             N_r=2 * (header['nr_pixel_x'] + header['nr_pixel_y']))
            return [position, data], header, plot_data_type
        elif plot_data_type == 'radial':
            # Create cuts instead of returning the map data
            position, data = self.create_radial_profile(tbldata,
                                                        radial_parameter=self.parse_args.radial_parameter)
            return [position, data], header, plot_data_type
        else:
            raise ValueError(
                'The midplane map data cannot be delivered in the chosen format!')

    def read_vel_maps(self, filename):
        """Reads the data of a file with velocity channel maps.

        Args:
            filename (str): Filename to the fits file with velocity channel maps.

        Returns:
            Tuple with:
            1) Numpy array with one dimension for the channel, one dimension for
            the data (tokes I, Q, ...) and 2 dimensions for the pixel coordinates.
            2) Dictionary with information of the file header.
        """
        if os.path.isfile(self.path['results'] + filename + '_vel_' + str(1).zfill(4) + '.fits'):
            # Load data and header from extra data fits if available
            hdulist_extra = fits.open(
                self.path['results'] + filename + '_extra.fits')
        else:
            raise ValueError('Error: The chosen fits file ' +
                             filename + ' does not exist!')
        #: dict: Dictionary with the information in the header.
        header, plot_data_type = self.read_vel_maps_header(hdulist_extra)
        if plot_data_type == 'healpix':
            import healpy as hp
            #: Amount of arcseconds per pixel squared to convert flux from Jy/pixel into Jy/arcsec^2
            arcsec_squared_per_pixel = 4 * np.pi * \
                (180 / np.pi * 60 * 60) ** 2 / header['nr_pixel_x']
            #: List: Load the velocity channels into a list
            data_list = [hp.read_map(self.path['results'] + filename + '_vel_' + str(vch + 1).zfill(4) +
                                     '.fits', field=None, verbose=False) for vch in range(header['nr_channels'])]
            # Load the extra data (zeeman, column dens)
            data_extra = hp.read_map(
                self.path['results'] + filename + '_extra.fits', field=None, verbose=False)
            #: Numpy array for the vel_map data
            wmap_map = np.zeros(
                (self.n_quantities_map, header['nr_channels'], header['nr_pixel_x']))
            # Fill the tbldata with the input data
            for vch in range(header['nr_channels']):
                wmap_map[0:5, vch, ...] = data_list[vch][0:5, ...]
                # Add the extra data
                if len(data_extra.shape) == 2:
                    if vch < data_extra.shape[0]:
                        wmap_map[5, vch, ...] = data_extra[vch, ...]
                else:
                    if vch == 0:
                        wmap_map[5, 0, ...] = data_extra[...]
            # Change unit if per arcseconds was chosen
            if self.cmap_unit == 'total' or self.cmap_unit == 'px':
                pass
            else:
                wmap_map[0:4, :, :] /= arcsec_squared_per_pixel
            wmap_map = self.beam_conv(
                wmap_map, plot_data_type, np.sqrt(arcsec_squared_per_pixel))
            return wmap_map, header, plot_data_type
        else:
            #: Amount of arcseconds per pixel to convert flux from Jy/pixel into Jy/arcsec^2
            arcsec_squared_per_pixel = (2. * self.model.tmp_parameter['radius_x_arcsec'] / header['nr_pixel_x']) * \
                (2. *
                 self.model.tmp_parameter['radius_y_arcsec'] / header['nr_pixel_y'])
            #: List: Load the velocity channels into a list
            hdulist_list = [fits.open(self.path['results'] + filename + '_vel_' + str(vch + 1).zfill(4) + '.fits')
                            for vch in range(0, header['nr_channels'])]
            #: Numpy array for the vel_map data
            tbldata = np.zeros(
                (6, header['nr_channels'], header['nr_pixel_x'], header['nr_pixel_y']))
            # Fill the tbldata with the input data
            for vch in range(header['nr_channels']):
                tbldata[0:5, vch, :, :] = np.transpose(
                    hdulist_list[vch][0].data, (0, 2, 1))
                # Add the extra data
                if vch < hdulist_extra[0].data.shape[0]:
                    tbldata[5, vch, :, :] = np.transpose(
                        hdulist_extra[0].data, (0, 2, 1))[vch, ...]
            # Change unit if per arcseconds was chosen
            if self.cmap_unit == 'total' or self.cmap_unit == 'px':
                pass
            else:
                tbldata[0:4, :, :, :] /= arcsec_squared_per_pixel
            # Convolve with beam if necessary
            tbldata = self.beam_conv(
                tbldata, plot_data_type, np.sqrt(arcsec_squared_per_pixel))
            if plot_data_type == 'map':
                return tbldata, header, plot_data_type
            elif plot_data_type == 'cut':
                # Create cuts instead of returning the map data
                position, data = self.create_cut(tbldata, cut_parameter=self.parse_args.cut_parameter,
                                                 N_r=2 * (header['nr_pixel_x'] + header['nr_pixel_y']))
                return [position, data], header, plot_data_type
            elif plot_data_type == 'radial':
                # Create cuts instead of returning the map data
                position, data = self.create_radial_profile(tbldata,
                                                            radial_parameter=self.parse_args.radial_parameter)
                return [position, data], header, plot_data_type
            else:
                raise ValueError(
                    'The velocity map data cannot be delivered in the chosen format!')

    def read_int_vel_map(self, filename):
        """Reads the data of a file with integrated velocity channel map.

        Args:
            filename (str): Filename to the fits file with integrated velocity channel map.

        Returns:
            Tuple with:
            1) Numpy array with one dimension for the data (tokes I, Q, ...)
            and 2 dimensions for the pixel coordinates.
            2) Dictionary with information of the file header.
        """
        if os.path.isfile(self.path['results'] + filename + '.fits'):
            # Load data and header from fits if available
            hdulist = fits.open(self.path['results'] + filename + '.fits')
        else:
            raise ValueError(
                'Error: hdulist cannot be returned if no fits file exists!')
        #: dict: Dictionary with the information in the header.
        header, plot_data_type = self.read_vel_maps_header(hdulist)
        if plot_data_type == 'healpix':
            import healpy as hp
            wmap_map = hp.read_map(
                self.path['results'] + filename + '.fits', field=None, verbose=False)
            #: Amount of arcseconds per pixel squared to convert flux from Jy/pixel into Jy/arcsec^2
            arcsec_squared_per_pixel = 4 * np.pi * \
                (180 / np.pi * 60 * 60) ** 2 / header['nr_pixel_x']
            if self.cmap_unit == 'total' or self.cmap_unit == 'px':
                pass
            else:
                wmap_map[0:4, :] /= arcsec_squared_per_pixel
            wmap_map = self.beam_conv(
                wmap_map, plot_data_type, np.sqrt(arcsec_squared_per_pixel))
            return wmap_map, header, plot_data_type
        else:
            #: Numpy array for the midplane data
            tbldata = np.transpose(hdulist[0].data, (0, 2, 1))
            #: Amount of arcseconds per pixel to convert flux from Jy/pixel into Jy/arcsec^2
            arcsec_squared_per_pixel = (2. * self.model.tmp_parameter['radius_x_arcsec'] / header['nr_pixel_x']) * \
                (2. *
                 self.model.tmp_parameter['radius_y_arcsec'] / header['nr_pixel_y'])
            if self.cmap_unit == 'total' or self.cmap_unit == 'px':
                pass
            else:
                tbldata[0:4, :, :] /= arcsec_squared_per_pixel
            tbldata = self.beam_conv(
                tbldata, plot_data_type, np.sqrt(arcsec_squared_per_pixel))
            if plot_data_type == 'map':
                return tbldata, header, plot_data_type
            elif plot_data_type == 'cut':
                # Create cuts instead of returning the map data
                position, data = self.create_cut(tbldata, cut_parameter=self.parse_args.cut_parameter,
                                                 N_r=2 * (header['nr_pixel_x'] + header['nr_pixel_y']))
                return [position, data], header, plot_data_type
            elif plot_data_type == 'radial':
                # Create cuts instead of returning the map data
                position, data = self.create_radial_profile(tbldata,
                                                            radial_parameter=self.parse_args.radial_parameter)
                return [position, data], header, plot_data_type
            else:
                raise ValueError(
                    'The integrated velocity map data cannot be delivered in the chosen format!')

    def read_spectrum(self, filename):
        """Reads data from spectrum file into a numpy array.

        Args:
            filename (str): Filename to the spectrum file.

        Returns:
            Numpy array including the data from the spectrum.
        """
        # Load data and header from fits if available
        if os.path.isfile(self.path['results'] + filename + '.fits'):
            hdulist = fits.open(self.path['results'] + filename + '.fits')
        else:
            raise FileExistsError('--- Hint: No emission sed fits file exists (' +
                                  self.path['results'] + filename + '.fits' + ')')
        #: dict: Dictionary with the information in the header.
        header = self.read_spectrum_header(hdulist)
        #: Data from the file with the monte carlo results
        tbldata = hdulist[0].data[:, 0, :]
        # Apply intensity unit conversion
        if self.cmap_unit == 'nuF':
            for i_vch in range(header['nr_channels']):
                tbldata[0:5, i_vch] *= 1e-26 * header['frequency']
        return tbldata, header

    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    # The following functions are changing data read in by the file_io class.
    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------

    def beam_conv(self, tbldata, plot_data_type, arcsec_per_pixel):
        """Convolves tbldata with beam.

        Args:
            tbldata: N dimensional numpy array with results
                from Raytrace or Monte-Carlo simulations.
            plot_data_type (str): Description of the data contained in the fits
                (healpix, map) or chosen postprocessing of the data (cut, radial)
            arcsec_per_pixel (float): Amount of arcseconds per pixel

        Returns:
            3 dimensional numpy array with results from raytrace simulations
            convolved with the beam.
        """
        # If not chosen, skip convolution
        if self.beam_size is None:
            return tbldata
        #: float: Full width at half maximum of the beam size [pixel]
        fwhm_px = self.beam_size / arcsec_per_pixel
        if plot_data_type == 'healpix':
            # Load necessary modules
            from healpy import smoothing
        else:
            # Load necessary modules
            from astropy.convolution import Gaussian2DKernel, convolve
            #: float: Standard deviation of the beam size [pixel]
            stddev_px = fwhm_px / (2. * np.sqrt(2 * np.log(2)))
            #: 2D gaussian kernel to convolve the raytrace results
            gauss = Gaussian2DKernel(stddev_px, x_size=int(
                20 * stddev_px), y_size=int(20 * stddev_px))
        # Only convolve quantities related to a flux
        if 'dust' in self.parse_args.simulation_type:
            quantity_list = [0, 1, 2, 3, 4]
            if self.parse_args.simulation_type in ['dust_mc', 'dust_full']:
                quantity_list.extend([6, 7])
            for i_quantity in quantity_list:
                if plot_data_type == 'healpix':
                    for i_wl in range(np.size(tbldata, 1)):
                        tbldata[i_quantity, i_wl, ...] = smoothing(tbldata[i_quantity, i_wl, ...], verbose=False,
                                                                   fwhm=self.math.angle_conv(self.beam_size, 'arcsec'))
                        # P_l is not flux / beam but should be convolved
                        if self.cmap_unit == 'beam':
                            tbldata[i_quantity, i_wl, ...] *= np.pi * \
                                (self.beam_size / 2.) ** 2
                else:
                    for i_wl in range(np.size(tbldata, 1)):
                        tbldata[i_quantity, i_wl, ...] = convolve(
                            tbldata[i_quantity, i_wl, ...], gauss)
                        # P_l is not flux / beam but should be convolved
                        if self.cmap_unit == 'beam':
                            tbldata[i_quantity, i_wl, ...] *= np.pi * \
                                (self.beam_size / 2.) ** 2
                        # Due to FFT, some values are negativ but very small. Make them positive again
                        if i_quantity in [0, 4, 6, 7]:
                            tbldata[i_quantity, i_wl, ...] = np.absolute(
                                tbldata[i_quantity, i_wl, ...])
            # Set the vector field size to match nearly with the beam size
            self.vec_field_size = int(fwhm_px)
        elif self.parse_args.visualization_type == 'int_map':
            for i_quantity in [0, 1, 2, 3, 4]:
                if plot_data_type == 'healpix':
                    tbldata[i_quantity, ...] = smoothing(tbldata[i_quantity, ...], verbose=False,
                                                         fwhm=self.math.angle_conv(self.beam_size, 'arcsec'))
                    # P_l is not flux / beam but should be convolved
                    if i_quantity != 4:
                        tbldata[i_quantity, ...] *= np.pi * \
                            (self.beam_size / 2.) ** 2
                else:
                    tbldata[i_quantity, ...] = convolve(
                        tbldata[i_quantity, ...], gauss)
                    # P_l is not flux / beam but should be convolved
                    if i_quantity != 4 and self.cmap_unit == 'beam':
                        tbldata[i_quantity, ...] *= np.pi * \
                            (self.beam_size / 2.) ** 2
                    # Due to FFT, some values are negativ but very small. Make them positive again
                    if i_quantity in [0, 4]:
                        tbldata[i_quantity, ...] = np.absolute(
                            tbldata[i_quantity, ...])
        elif self.parse_args.visualization_type == 'vel_map':
            for i_quantity in [0, 1, 2, 3, 5]:
                if plot_data_type == 'healpix':
                    for i_channel in range(np.size(tbldata, 1)):
                        tbldata[i_quantity, i_channel, ...] = smoothing(tbldata[i_quantity, i_channel, ...],
                                                                        verbose=False, fwhm=self.math.angle_conv(self.beam_size, 'arcsec'))
                        # P_l is not flux / beam but should be convolved
                        if i_quantity != 4 and self.cmap_unit == 'beam':
                            tbldata[i_quantity, i_channel, ...] *= np.pi * \
                                (self.beam_size / 2.) ** 2
                else:
                    for i_channel in range(np.size(tbldata, 1)):
                        tbldata[i_quantity, i_channel, ...] = convolve(
                            tbldata[i_quantity, i_channel, ...], gauss)
                        # P_l is not flux / beam but should be convolved
                        if i_quantity != 4 and self.cmap_unit == 'beam':
                            tbldata[i_quantity, ...] *= np.pi * \
                                (self.beam_size / 2.) ** 2
                        # Due to FFT, some values are negativ but very small. Make them positive again
                        if i_quantity in [0, 4]:
                            tbldata[i_quantity, i_channel, ...] = np.absolute(
                                tbldata[i_quantity, i_channel, ...])
        return tbldata

    def read_polarization_vectors(self, tbldata, min_intensity=0.):
        """Creates polarization vectors from raytrace results.

        Args:
            tbldata: 3 dimensional numpy array with results
                from raytrace simulations.
            min_intensity (float): Minimum flux (in Jy) that is required 
                to obtain the corresponding vector.

        Returns:
            Numpy array with 2 dimensions for the pixel coordinates
            and one dimension for the polarization angle.
        """
        #: int: Number of bins of the numpy array
        bins_x = len(tbldata[0, 0, :])
        bins_y = len(tbldata[0, :, 0])
        #: int: Number of bins per dimension per polarization vector
        if bins_x % 2 != 0 and self.vec_field_size % 2 == 0:
            offset_x = (bins_x - 1) % self.vec_field_size
        else:
            offset_x = bins_x % self.vec_field_size
        while offset_x % 2 != 0:
            offset_x += self.vec_field_size
        vector_bins_x = int(
            (bins_x - offset_x) / self.vec_field_size)
        offset_x = int(offset_x / 2)

        if bins_y % 2 != 0 and self.vec_field_size % 2 == 0:
            offset_y = (bins_y - 1) % self.vec_field_size
        else:
            offset_y = bins_y % self.vec_field_size
        while offset_y % 2 != 0:
            offset_y += self.vec_field_size
        vector_bins_y = int(
            (bins_y - offset_y) / self.vec_field_size)
        offset_y = int(offset_y / 2)
        #: Stokes I flux averaged over the polarization vector area
        i_vec_avg = np.zeros((vector_bins_x, vector_bins_y))
        #: Stokes Q flux averaged over the polarization vector area
        q_vec_avg = i_vec_avg.copy()
        #: Stokes U flux averaged over the polarization vector area
        u_vec_avg = i_vec_avg.copy()
        #: Degree of polarization averaged over the polarization vector area
        p_vec_avg = i_vec_avg.copy()
        #: Count how many pixel with intensity are added to the vector
        p_vec_count = i_vec_avg.copy()
        #: x center position of the averaged vector area
        x_vec_avg = i_vec_avg.copy()
        #: y center position of the averaged vector area
        y_vec_avg = i_vec_avg.copy()
        #: Polarization angle calculated from the averaged Q and U components
        pol_angle = i_vec_avg.copy()
        #: Numpy array that includes the polarization vectors for the whole map
        vec_field_data = np.zeros((vector_bins_x, vector_bins_y, 5))
        # Init steps to jump over the middle row/column is necessary
        middle_step_x = 0
        middle_step_y = 0
        for i_x in range(offset_x, bins_x - offset_x):
            # Skip middle cell if vector size is even and pixel size is odd
            if bins_x % 2 != 0 and self.vec_field_size % 2 == 0 \
                    and i_x == int((bins_x - 1) / 2):
                if middle_step_x == 0:
                    middle_step_x = 1
                continue
            for i_y in range(offset_y, bins_y - offset_y):
                # Skip middle cell if vector size is even and pixel size is odd
                if bins_y % 2 != 0 and self.vec_field_size % 2 == 0 \
                        and i_y == int((bins_y - 1) / 2):
                    if middle_step_y == 0:
                        middle_step_y = 1
                    continue
                #: int: X-axis index to the vector
                i_x_vec = int((i_x - offset_x - middle_step_x) /
                              self.vec_field_size)
                #: int: Y-axis index to the vector
                i_y_vec = int((i_y - offset_y - middle_step_y) /
                              self.vec_field_size)
                # Averaging Stokes I
                i_vec_avg[i_x_vec, i_y_vec] += tbldata[0,
                                                       i_x, i_y] / self.vec_field_size ** 2
                # Averaging Stokes Q
                q_vec_avg[i_x_vec, i_y_vec] += tbldata[1,
                                                       i_x, i_y] / self.vec_field_size ** 2
                # Averaging Stokes U
                u_vec_avg[i_x_vec, i_y_vec] += tbldata[2,
                                                       i_x, i_y] / self.vec_field_size ** 2
                # Calculating pol vector length
                if tbldata[0, i_x, i_y] > 1e-200:
                    p_vec_avg[i_x_vec, i_y_vec] += 1e2 * np.sqrt(
                        tbldata[1, i_x, i_y] ** 2 + tbldata[2, i_x, i_y] ** 2) / (tbldata[0, i_x, i_y])
                    p_vec_count[i_x_vec,
                                i_y_vec] += 1
        for i_x_vec in range(vector_bins_x):
            for i_y_vec in range(vector_bins_y):
                if p_vec_count[i_x_vec, i_y_vec] <= 0:
                    p_vec_avg[i_x_vec, i_y_vec] = np.nan
                else:
                    i_vec_avg[i_x_vec, i_y_vec] /= p_vec_count[i_x_vec, i_y_vec]
                    if i_vec_avg[i_x_vec, i_y_vec] < min_intensity:
                        p_vec_avg[i_x_vec, i_y_vec] = np.nan
                    else:
                        p_vec_avg[i_x_vec,
                                  i_y_vec] /= p_vec_count[i_x_vec, i_y_vec]
                # Calculating the x and y positions
                # X axis
                x_vec_avg[i_x_vec, i_y_vec] = ((
                    i_x_vec + 0.5) * self.vec_field_size + offset_x) / bins_x
                if bins_x % 2 != 0 and self.vec_field_size % 2 == 0:
                    if i_x_vec == int(vector_bins_x / 2):
                        x_vec_avg[i_x_vec, i_y_vec] += 0.5 / bins_x
                    elif i_x_vec > int(vector_bins_x / 2):
                        x_vec_avg[i_x_vec, i_y_vec] += 1 / bins_x
                # Y axis
                y_vec_avg[i_x_vec, i_y_vec] = ((
                    i_y_vec + 0.5) * self.vec_field_size + offset_y) / bins_y
                if bins_y % 2 != 0 and self.vec_field_size % 2 == 0:
                    if i_y_vec == int(vector_bins_y / 2):
                        y_vec_avg[i_x_vec, i_y_vec] += 0.5 / bins_y
                    elif i_y_vec > int(vector_bins_y / 2):
                        y_vec_avg[i_x_vec, i_y_vec] += 1 / bins_y
                # Calculating the polarization angle from averaged Q and U components
                pol_angle[i_x_vec, i_y_vec] = self.math.angle_from_stokes(
                    q_vec_avg[i_x_vec, i_y_vec], u_vec_avg[i_x_vec, i_y_vec])
                if np.isnan(pol_angle[i_x_vec, i_y_vec]):
                    return None
        # Put the 2D polarization vector data into numpy array
        vec_field_data[:, :, 0] = np.cos(pol_angle) * p_vec_avg
        vec_field_data[:, :, 1] = np.sin(pol_angle) * p_vec_avg
        vec_field_data[:, :, 3] = x_vec_avg
        vec_field_data[:, :, 4] = y_vec_avg
        return vec_field_data

    def read_polarization_vectors_healpix(self, tbldata, header_dict):
        """Creates polarization vectors from raytrace healpix results.

        Args:
            tbldata: 2 dimensional numpy array with results
                from raytrace healpix simulations.
            header_dict (dict): Dictionary with header information.
        """
        from scipy.optimize import fsolve
        import matplotlib.pyplot as plt
        import healpy as hp

        def aux_angle(phi):
            # x = np.linspace(-np.pi, np.pi, 101)
            initial_guess = 0.0

            def func(x): return np.pi * np.sin(x) - 2.0 * phi - np.sin(2 * phi)

            solution = fsolve(func, initial_guess, xtol=0.001)
            Phi = solution[0]
            return Phi

        def mollweide2plane(lon, lat, lon0):
            theta = aux_angle(lat)
            x = 2.0 * (lon - lon0) * np.cos(theta) / np.pi
            y = np.sin(theta)
            return x, y

        def plane2mollweide(x, y, lon0):
            theta = np.arcsin(y)
            lon = lon0 + np.pi * x / (2 * np.cos(theta))
            lat = np.arcsin((2 * theta + np.sin(2 * theta)) / np.pi)
            return lon, lat

        def in_ellipse(x, y, dx, dy):
            if abs(x) > dx or abs(y) > dy:
                return False
            elif (x ** 2 / dx ** 2 + y ** 2 / dy ** 2) > 1.0:
                return False
            else:
                return True

        #: int: Number of bins per dimension per polarization vector
        vector_bins_x = 60
        vector_bins_y = 30
        #: Stokes I flux averaged over the polarization vector area
        x = np.linspace(-2, 2, vector_bins_x)
        #: Stokes Q flux averaged over the polarization vector area
        y = np.linspace(-1, 1, vector_bins_y)
        #: Polarization angle calculated from the averaged Q and U components
        pol_angle = np.zeros((vector_bins_x, vector_bins_y))
        #: Numpy array that includes the polarization vectors for the whole map
        # vec_field_data = np.zeros((vector_bins_x, vector_bins_y, 3))

        pol_x1 = []
        pol_x2 = []
        pol_y1 = []
        pol_y2 = []

        pol_a = np.zeros(header_dict['nsides'] * header_dict['nsides'] * 12)

        for i_x in range(vector_bins_x):
            for i_y in range(vector_bins_y):
                if in_ellipse(x[i_x], y[i_y], 2, 1):
                    lon1, lat1 = plane2mollweide(-x[i_x], -y[i_y], np.pi)
                    hp_id = hp.ang2pix(header_dict['nsides'],
                                       lat1 + np.pi / 2., lon1 + np.pi)

                    # Calculating the polarization angle from averaged Q and U components
                    pol_angle[i_x, i_y] = self.math.angle_from_stokes(
                        tbldata[1, hp_id], tbldata[2, hp_id])

                    if tbldata[0, hp_id] > 1e-200:
                        p_vec_avg = np.sqrt(
                            tbldata[1, hp_id] ** 2 + tbldata[2, hp_id] ** 2) / tbldata[0, hp_id] * 0.04
                    else:
                        p_vec_avg = 0

                    pol_a[hp_id] = pol_angle[i_x, i_y]
                    pol_a[hp_id] += 2 * np.pi

                    ddx = 0.36 * np.cos(pol_angle[i_x, i_y])
                    ddy = 0.36 * np.sin(pol_angle[i_x, i_y])

                    lon2 = lon1 - 0.5 * np.pi * ddx
                    lat2 = lat1 - 0.5 * np.pi * ddy

                    dx1, dy1 = mollweide2plane(lon2, lat2, 1.0 * np.pi)

                    lon3 = lon1 + 0.5 * np.pi * ddx
                    lat3 = lat1 + 0.5 * np.pi * ddy

                    dx2, dy2 = mollweide2plane(lon3, lat3, 1.0 * np.pi)

                    dx = dx2 - dx1
                    dy = dy2 - dy1

                    l = np.sqrt(dx ** 2 + dy ** 2)
                    dx *= p_vec_avg / l
                    dy *= p_vec_avg / l

                    pol_x1.append(x[i_x] - dx)
                    pol_x2.append(x[i_x] + dx)

                    pol_y1.append(y[i_y] - dy)
                    pol_y2.append(y[i_y] + dy)

        for i in range(0, len(pol_x1)):
            dx1 = pol_x1[i]
            dx2 = pol_x2[i]

            dy1 = pol_y1[i]
            dy2 = pol_y2[i]

            p_x = [dx1, dx2]
            p_y = [dy1, dy2]

            plt.plot(p_x, p_y, '-', color='white')

    def create_cut(self, tbldata, cut_parameter, N_r=1000):
        """"Calculates cut through map data depending on cut_parameters

        Args:
            tbldata: A 3D array with a size of n x nr_of_pixel x nr_of_pixel.
                (n is number of quantities/wavlelengths/frequencies...)
            cut_parameter (List[float, float, float]): The cut parameter.
                (cut angle, cut center position x and y)
            N_r (int): Number of positions along the cut.

        Returns:
            List, List: Position and Radial profile data as 1D array and nD array.
        """
        # Get cut parameter
        cut_angle = cut_parameter[0] / 180. * np.pi
        if len(cut_parameter) == 1:
            center_pos = [0, 0]
            width = 0
            n_width = 1
        elif len(cut_parameter) == 3:
            center_pos = cut_parameter[1:3]
            width = 0
            n_width = 1
        elif len(cut_parameter) == 5:
            center_pos = cut_parameter[1:3]
            width = float(cut_parameter[3])
            n_width = int(cut_parameter[4])
        else:
            raise ValueError(
                'Cut parameter have to be the cut angle and optionally the cut center positions!')
        sidelength_x = 2. * \
            self.model.tmp_parameter['radius_x_' + self.ax_unit]
        sidelength_y = 2. * \
            self.model.tmp_parameter['radius_y_' + self.ax_unit]
        # Calc maximum length in map
        max_len = np.sqrt(abs(center_pos[0] + sidelength_x / 2) ** 2 +
                          abs(center_pos[1] + sidelength_y / 2) ** 2)
        # nr_offset_index = len(tbldata.shape) - 2
        nr_pixel_x = tbldata.shape[-2]
        nr_pixel_y = tbldata.shape[-1]
        cut_position = np.linspace(-max_len, max_len, N_r)
        cut_data = np.zeros(np.shape(tbldata)[0:-2] + (N_r,))
        for i_r in range(N_r):
            for i_width in range(n_width):
                d_width = 0.5 * width * (i_width + 0.5) / (n_width)
                pos = [cut_position[i_r] * np.cos(cut_angle) + d_width * np.sin(cut_angle) + center_pos[0],
                       cut_position[i_r] * np.sin(cut_angle) + d_width * np.cos(cut_angle) + center_pos[1]]
                pos_id_x = np.multiply(
                    np.divide(np.add(pos[0], sidelength_x / 2.), sidelength_x), nr_pixel_x)
                pos_id_y = np.multiply(
                    np.divide(np.add(pos[1], sidelength_y / 2.), sidelength_y), nr_pixel_y)
                if pos_id_x < nr_pixel_x and pos_id_x >= 0 and pos_id_y < nr_pixel_y and pos_id_y >= 0:
                    cut_data[..., i_r] += tbldata[...,
                                                  int(pos_id_x), int(pos_id_y)] / n_width
                else:
                    cut_position[i_r] = np.nan
        # Remove nan values
        not_nan_list = ~np.isnan(cut_position)
        cut_position = cut_position[not_nan_list]
        cut_data = cut_data[..., not_nan_list]
        # Combine values of the same pixel to one value
        # tmp_var = np.var(cut_data, axis=tuple(x for x in range(nr_offset_index)))
        # i_r = 1
        # while i_r <= len(cut_position) - 1:
        #     if np.all(cut_data[..., i_r] == cut_data[..., i_r - 1]):
        #         if i_r == len(cut_position) - 1:
        #             cut_position[i_r - 1] = (cut_position[i_r] + cut_position[i_r - 1]) / 2.
        #         elif np.all(cut_data[..., i_r] != cut_data[..., i_r + 1]):
        #             cut_position[i_r - 1] = (cut_position[i_r] + cut_position[i_r - 1]) / 2.
        #         cut_data = np.delete(cut_data, i_r, axis=-1)
        #         cut_position = np.delete(cut_position, i_r)
        #     else:
        #         i_r += 1
        return cut_position, cut_data

    def create_radial_profile(self, tbldata, radial_parameter, N_r=100, subpixel=4):
        """"Calculates azimuthally averaged radial profiles through map data 
        depending on radial_parameters.

        Args:
            tbldata: A 3D array with a size of n x nr_of_pixel x nr_of_pixel.
                (n is number of quantities/wavlelengths/frequencies...)
            radial_parameter (List[float, float]): The radial parameter.
                (center position x and y)
            N_r (int): Number of positions along the radial profile.
            subpixel (int): define the mount of subpixel per axis.

        Returns:
            List, List: Position and Radial profile data as 1D array and nD array.
        """
        # Get cut parameter
        center_pos = radial_parameter[0:2]
        sidelength_x = 2. * \
            self.model.tmp_parameter['radius_x_' + self.ax_unit]
        sidelength_y = 2. * \
            self.model.tmp_parameter['radius_y_' + self.ax_unit]
        # Calc maximum length in map
        max_len = np.sqrt(abs(center_pos[0] + sidelength_x / 2) ** 2 +
                          abs(center_pos[1] + sidelength_y / 2) ** 2)
        # Get number of axes not used for the pixel
        # nr_offset_index = len(tbldata.shape) - 2
        # Get number of pixel per axis
        nr_pixel_x = tbldata.shape[-2]
        nr_pixel_y = tbldata.shape[-1]
        # Update subpixel, if too small
        while subpixel * nr_pixel_x < 2 * N_r * np.sqrt(2):
            subpixel += 1
        radial_position = np.linspace(0, max_len, N_r)
        radial_number = np.zeros(N_r)
        radial_data = np.zeros(np.shape(tbldata)[0:-2] + (N_r,))
        for i_x in range(subpixel * nr_pixel_x):
            for i_y in range(subpixel * nr_pixel_y):
                index = [i_x, i_y]
                pos = np.subtract(np.multiply(np.divide(np.add(index, 0.5),
                                                        [subpixel * nr_pixel_x, subpixel * nr_pixel_y]),
                                              [sidelength_x, sidelength_y]),
                                  np.divide([sidelength_x, sidelength_y], 2.))
                pos_r = np.linalg.norm(np.subtract(pos, center_pos))
                i_r = int(pos_r / max_len * (N_r - 1))
                radial_data[..., i_r] += tbldata[...,
                                                 int(i_x / float(subpixel)), int(i_y / float(subpixel))]
                radial_number[i_r] += 1
        # Normalization
        for i_r in range(N_r):
            if radial_number[i_r] > 0:
                radial_data[..., i_r] /= radial_number[i_r]
            else:
                radial_data[..., i_r] = 0.
        return radial_position, radial_data

    def create_azimuthal_profile(self, tbldata, azimuthal_parameter, N_ph=180, subpixel=4):
        """"Calculates radially averaged azimuthal profiles through map data 
        depending on azimuthal_parameter.

        Args:
            tbldata: A 3D array with a size of n x nr_of_pixel x nr_of_pixel.
                (n is number of quantities/wavlelengths/frequencies...)
            azimuthal_parameter (List[float, float, float, float]): The azimuthal parameter.
                (center position x and y, inclination, inclination position angle, min and max radius of ring)
            N_ph (int): Number of positions along the azimuthal profile.
            subpixel (int): define the mount of subpixel per axis.

        Returns:
            List, List: Position and azimuthal profile data as 1D array and nD array.
        """
        # Get cut parameter
        center_pos = azimuthal_parameter[0:2]
        sidelength_x = 2. * \
            self.model.tmp_parameter['radius_x_' + self.ax_unit]
        sidelength_y = 2. * \
            self.model.tmp_parameter['radius_y_' + self.ax_unit]
        # Convert input limits to used unit
        R_min = self.math.length_conv(azimuthal_parameter[4], self.ax_unit)
        R_max = self.math.length_conv(azimuthal_parameter[5], self.ax_unit)
        # Get number of pixel per axis
        nr_pixel_x = tbldata.shape[-2]
        nr_pixel_y = tbldata.shape[-1]
        azimuthal_position = np.linspace(0, 2. * np.pi, N_ph)
        azimuthal_number = np.zeros(N_ph)
        azimuthal_data = np.zeros(np.shape(tbldata)[0:-2] + (N_ph,))
        for i_x in range(subpixel * nr_pixel_x):
            for i_y in range(subpixel * nr_pixel_y):
                index = [i_x, i_y]
                pos = np.subtract(np.multiply(np.divide(np.add(index, 0.5),
                                                        [subpixel * nr_pixel_x, subpixel * nr_pixel_y]),
                                              [sidelength_x, sidelength_y]),
                                  np.divide([sidelength_x, sidelength_y], 2.))
                real_pos = self.math.apply_inclination(
                    pos=np.subtract(pos, center_pos),
                    inclination=azimuthal_parameter[2],
                    inc_PA=azimuthal_parameter[3],
                    inc_offset=azimuthal_parameter[6], inv=True)
                if(R_min <= np.linalg.norm(real_pos) <= R_max):
                    pos_ph = np.arctan2(-real_pos[1], -real_pos[0])
                    i_ph = int((pos_ph + np.pi) / (2. * np.pi) * N_ph)
                    azimuthal_data[..., i_ph] += tbldata[..., int(i_x / float(subpixel)),
                                                         int(i_y / float(subpixel))]
                    azimuthal_number[i_ph] += 1

        # Normalization
        for i_ph in range(N_ph):
            if azimuthal_number[i_ph] > 0:
                azimuthal_data[..., i_ph] /= azimuthal_number[i_ph]
            else:
                azimuthal_data[..., i_ph] = 0.
        return azimuthal_position, azimuthal_data

    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------
    # The following functions checking integrity of the data or updating variables.
    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------

    def get_updated_simulation_type(self, simulation_type):
        """Check and update the chosen simulation type.

        Args:
            simulation_type (str): Type of the simulation (see polaris-run.in).

        Returns:
            simulation_type (str): Updated type of simulation.
        """
        # Start with simulation type not found among the valid choices
        found = False
        # Is ONLY one of the main simulation types chosen?
        for types in self.main_simu_types:
            if types == simulation_type or types in simulation_type.split('_')[0]:
                found = True
        # In the case of ray, multiple additions are valid to choose aligment
        if 'dust_' in simulation_type and simulation_type not in ['dust_mc', 'dust_full']:
            # Sort the list of aligment mechanisms to have a unique name for each combination
            simu_list = simulation_type.replace('dust_', '').split('_')
            simu_list.sort()
            simulation_type = 'dust_' + '_'.join(simu_list)
            # Each chosen alignment mechanisms needs to be valid
            if any(align not in self.alignment_mechanisms for align in simu_list):
                found = False
        # Error message if chosen simulation_types are not valid
        if not found:
            raise ValueError(
                'Wrong simulation_type chosen (run_plaris -h for info)')
        return simulation_type

    def check_fits_data(self, hdulist):
        """Check the way data is saved in the fits file.

        Args:
            hdulist: Fits HDU list object.

        Returns:
            str: Identifier of the containing data.
        """
        # Get type of plot created from map results
        try:
            if hdulist[1].header['PIXTYPE'].lower() == 'healpix':
                return 'healpix'
        except:
            if self.parse_args.visualization_type in ['midplane', 'map', 'int_map', 'vel_map', 'mag_field', 'velocity', 'custom']:
                if self.parse_args.cut_parameter is not None and self.parse_args.radial_parameter is not None:
                    raise ValueError('Choose only one of --cut and --radial!')
                elif self.parse_args.cut_parameter is not None:
                    return 'cut'
                elif self.parse_args.radial_parameter is not None:
                    return 'radial'
                else:
                    return 'map'
            elif self.parse_args.visualization_type in ['sed', 'spectrum']:
                return 'spectrum'
        return None

    def check_sidelengths(self, header_dict):
        """Check if sidelengths in the fits header are the same as set in the model.

        Args:
            header_dict (dict): Dictionary with header information.
        """
        # Check if the sidelength of the image ist the same as of the model
        if abs(header_dict['sidelength_x'] - 2. * self.model.tmp_parameter['radius_x_m']) > \
                (header_dict['sidelength_x'] + 2. * self.model.tmp_parameter['radius_x_m']) * 1e-10 or \
                abs(header_dict['sidelength_y'] - 2. * self.model.tmp_parameter['radius_y_m']) > \
                (header_dict['sidelength_y'] + 2. * self.model.tmp_parameter['radius_y_m']) * 1e-10:
            print('HINT: The map sidelengths are not the same as the model extent defined in model.py!\n'
                  'This occurs usually if a zoom_factor was chosen, but it might be a problem if not.')
            self.model.adjust_extent(
                header_dict['sidelength_x'], header_dict['sidelength_y'])

    def check_distance(self, header_dict):
        """Check if distance in the fits header is the same as set in the model.

        Args:
            header_dict (dict): Dictionary with header information.
        """
        # Check if the distance in the fits file is the same as in the model
        if abs(header_dict['distance'] - self.model.parameter['distance']) > header_dict['distance'] * 1e-10:
            print(
                'HINT: The distance set in model.py is not the same as in the Polaris fits output file!')
            self.model.adjust(header_dict['distance'], 'distance')
