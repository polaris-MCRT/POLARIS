#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np


class DustCreator:
    """The DustCreator class is able to create various dust components for POLARIS.
    """

    def __init__(self, file_io, parse_args, parameter):
        """Initialisation of the dust creator parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
            parameter (dict): Dictionary with parameter to create dust grain catalog.
        """
        self.file_io = file_io
        self.parse_args = parse_args
        self.parameter = parameter

        # Get math module
        from modules.math import Math
        self.math = Math()

        # Elements to create scattering matrix
        self.elements = np.zeros(16)
        self.elements[0] = 1  # S11
        self.elements[1] = 2  # S12
        self.elements[2] = 0  # S13
        self.elements[3] = 0  # S14
        self.elements[4] = 2  # S21
        self.elements[5] = 1  # S22
        self.elements[6] = 0  # S23
        self.elements[7] = 0  # S24
        self.elements[8] = 0  # S31
        self.elements[9] = 0  # S32
        self.elements[10] = 3  # S33
        self.elements[11] = 4  # S34
        self.elements[12] = 0  # S41
        self.elements[13] = 0  # S42
        self.elements[14] = -4  # S43
        self.elements[15] = 3  # S44

        self.nr_matrix_elements = 4

    def create_polaris_dust_from_dustem(self):
        """Create dust catalog from dustem input files.
        """
        # Import method to check if two values are close to each other
        from math import isclose
        print('--- Reading dust component: ' + self.parameter['input_file'])
        # Read the list of usable wavelengths from dustem
        self.read_dustem_wavelengths()
        # Read the cross sections from dustem
        self.read_dustem_cross_sections()
        # Read the Henyey-Greenstein parameter g
        self.read_dustem_g_factors()
        # Read the heat capacity
        self.read_dustem_calorimetry()
        print('------ Finished reading of dustem data!')
        # Write dust component data as POLARIS file
        self.write_dust_to_file(' (from DUSTEM code)')
        # Write dust component calorimetry so POLARIS can read it
        self.write_dust_calorimetry_to_file()
        print('------ Finished creation of POLARIS dust component!')

    def create_polaris_dust_from_trust(self):
        """Create dust catalog from trust benchmark input files.
        """
        # Import method to check if two values are close to each other
        from math import isclose
        print('--- Reading dust component: ' + self.parameter['input_file'])
        # Read optical properties from trus dust grain file
        self.read_trust_grain_inputs()
        # Read the heat capacity
        self.read_trust_calorimetry()
        # Read scattering matrix
        #self.read_trust_scattering()
        # Add one size
        if 'Sil' in self.parameter['input_file']:
            self.add_size(3.50e-10)
        print('------ Finished reading of TRUST benchmark data!')
        # Write dust component data as POLARIS file
        self.write_dust_to_file(' (from TRUST benchmark)')
        # Write dust component calorimetry so POLARIS can read it
        self.write_dust_calorimetry_to_file()
        print('------ Finished creation of POLARIS dust component!')
        #self.write_dust_scat_to_file_init()
        #self.write_dust_scat_to_file()
        print('------ Finished creation of polaris scattering data!')

    def create_polaris_dust_from_mie(self):
        """Create dust catalog with Mie theory.

        Args:
            parameter (dict): Dictionary with parameter to create dust grain catalog.
        """
        self.calc_mie()
        print('------ Finished creation of MIEX dust data!                                 ')
        self.write_dust_to_file(' (made with MIE)')
        print('------ Finished creation of polaris dust component!')
        self.write_dust_scat_to_file_init()
        self.write_dust_scat_to_file()
        print('------ Finished creation of polaris scattering data!')

    def calc_mie(self):
        """Calculate the Mie scattering with PyMieScatt.
        """
        from PyMieScatt import MatrixElements, MieQ
        nk_data = np.genfromtxt(self.file_io.path['dust_nk'] + self.parameter['input_file'])

        # Adjust unit
        if self.parameter['nk_file_microns']:
            nk_data[:, 0] *= 1e-6

        # Set wavelength list
        self.parameter['wavelength_list'] = nk_data[:, 0]

        # Init arrays for dust properties
        self.init_optical_properties()

        import warnings
        warnings.filterwarnings('ignore')
        for i_a, a_eff in enumerate(self.parameter['size_list']):
            for i_wl, (wavelength, n, k) in enumerate(nk_data):
                print('Calculating the Mie optical properties for size ID:',
                    str(i_a + 1) + '/' + str(self.parameter['size_list'].size), 'and wl ID:',
                        str(i_wl + 1) + '/' + str(nk_data.shape[0]) + '      ', end='\r')
                efficiencies_dict = MieQ(n + 1j * k, wavelength * 1e9, 2 * a_eff * 1e9, asDict=True)
                for efficiency in ['Qext', 'Qsca', 'Qabs', 'g']:
                    self.parameter[efficiency][i_a, i_wl] = efficiencies_dict[efficiency]
                for i_theta, theta in enumerate(self.parameter['theta_angle_list']):
                    cos_theta = np.cos(theta / 180. * np.pi)
                    self.parameter['scat_matrix'][i_a, i_wl, 0, 0, i_theta, :] = MatrixElements(n + 1j * k,
                        wavelength * 1e9, 2 * a_eff * 1e9, cos_theta)

    def read_dustem_wavelengths(self):
        #Open wavelength file
        wavelengths_file = open(self.file_io.path['dustem'] + 'oprop/' + self.parameter['dustem_wl_file'] + '.DAT')
        # Skip header
        [wavelengths_file.readline() for i in range(3)]
        # Set number of wavlengths
        self.parameter['wavelength_list'] = np.zeros(int(wavelengths_file.readline().split()[0]))
        # Set the individual wavelength
        for i_wl in range(self.parameter['wavelength_list'].size):
            data = wavelengths_file.readline().split()
            # Convert the wavelength from microns to m
            self.parameter['wavelength_list'][i_wl] = float(data[0]) * 1e-6

    def read_dustem_cross_sections(self):
        # Open cross sections file
        dustem_cross_sections_file = open(self.file_io.path['dustem'] + 'oprop/Q_' +
            self.parameter['input_file'] + '.DAT')
        # Skip header
        [dustem_cross_sections_file.readline() for i in range(4)]
        # Set number of dust grain sizes
        self.parameter['size_list'] = np.zeros(int(dustem_cross_sections_file.readline().split()[0]))
        # Set dust grain sizes
        aeff_list = dustem_cross_sections_file.readline().split()
        for i, size in enumerate(self.parameter['size_list']):
            # Convert the grain sizes from microns to m
            self.parameter['size_list'][i] = float(aeff_list[i]) * 1e-6
        # Skip space between
        [dustem_cross_sections_file.readline() for i in range(2)]
        # Update the cross sections to fit with the sizes and wavelengths
        self.init_optical_properties()
        # Set cross sections qabs
        for i_wl, wl in enumerate(self.parameter['wavelength_list']):
            data_list = dustem_cross_sections_file.readline().split()
            for i_dust, size in enumerate(self.parameter['size_list']):
                self.parameter['Qabs'][i_dust, i_wl] = float(data_list[i_dust])
        # Skip space between
        dustem_cross_sections_file.readline()
        # Set cross sections qsca and calculate qext
        for i_wl, wl in enumerate(self.parameter['wavelength_list']):
            data_list = dustem_cross_sections_file.readline().split()
            for i_dust, size in enumerate(self.parameter['size_list']):
                self.parameter['Qsca'][i_dust, i_wl] = float(data_list[i_dust])
                self.parameter['Qext'][i_dust, i_wl] = self.parameter['Qabs'][i_dust, i_wl] + \
                    self.parameter['Qsca'][i_dust, i_wl]

    def read_dustem_g_factors(self):
        # Open g factor file
        dustem_g_factor_file = open(self.file_io.path['dustem'] + 'oprop/G_' + self.parameter['input_file'] + '.DAT')
        # Skip header
        [dustem_g_factor_file.readline() for i in range(6)]
        # Set number of dust grain sizes
        nr_of_dust_species = int(dustem_g_factor_file.readline().split()[0])
        # Check if number of dust grain sizes agree with previous value
        if len(self.parameter['size_list']) != nr_of_dust_species:
            raise ValueError('grain sizes do not fit!')
        # Skip space between
        [dustem_g_factor_file.readline() for i in range(2)]
        # Set the g factor
        for i_wl, wl in enumerate(self.parameter['wavelength_list']):
            data_list = dustem_g_factor_file.readline().split()
            for i_dust, size in enumerate(self.parameter['size_list']):
                self.parameter['g'][i_dust, i_wl, 0] = float(data_list[i_dust])

    def read_dustem_calorimetry(self):
        # Open heat capacity file
        dustem_hcap_file = open(self.file_io.path['dustem'] + 'hcap/C_' + self.parameter['input_file'] + '.DAT')
        # Skip header
        [dustem_hcap_file.readline() for i in range(8)]
        # Set number of dust grain sizes
        nr_of_dust_species = int(dustem_hcap_file.readline().split()[0])
        # Check if number of dust grain sizes agree with previous value
        if len(self.parameter['size_list']) != nr_of_dust_species:
            raise ValueError('grain sizes do not fit!')
        # Set the dust grain sizes
        aeff = [float(i) for i in dustem_hcap_file.readline().split()]
        for i_dust, size in enumerate(self.parameter['size_list']):
            # Check if dust grain sizes agree with previous values
            if aeff[i_dust] * 1e-6 != self.parameter['size_list'][i_dust]:
                raise ValueError('grain sizes do not fit!')
        # Set number of temperature for heat capacity
        self.parameter['temp_list'] = np.zeros(int(dustem_hcap_file.readline().split()[0]))
        # Update the calorimetry array to fit with the temperatures
        self.init_calorimetries(self.parameter['temp_list'].size)
        # Set the temperatures and heat capacities
        for i_temp, temp in enumerate(self.parameter['temp_list']):
            data_list = dustem_hcap_file.readline().split()
            self.parameter['temp_list'][i_temp] = 10 ** float(data_list[0])
            for i_dust, size in enumerate(self.parameter['size_list']):
                # from erg/K/cm^3 = 1e-7 * 1e6 * J/m^3
                self.parameter['calorimetry_list'][i_dust, i_temp] = 1e-7 * 1e6 * 10 ** float(data_list[i_dust + 1])

    def read_trust_grain_inputs(self):
        #Open wavelength file
        trust_grain_file = open(self.file_io.path['trust'] + 'GrainInputs/' + self.parameter['input_file'] + '.dat')
        # Skip header
        [trust_grain_file.readline() for i in range(7)]
        # Set dust grain sizes
        self.parameter['size_list'] = np.zeros(int(trust_grain_file.readline().split()[1]))
        # Set number of wavlengths
        self.parameter['wavelength_list'] = np.zeros(int(trust_grain_file.readline().split()[1]))

        # Skip variable header part
        line = None
        while line != '\n':
            line = trust_grain_file.readline()

        # Update the cross sections to fit with the sizes and wavelengths
        self.init_optical_properties()
        # Set optical properties
        for i_dust, size in enumerate(self.parameter['size_list']):
            self.parameter['size_list'][i_dust] = float(trust_grain_file.readline().split()[0]) * 1e-6
            trust_grain_file.readline()
            for i_wl in range(self.parameter['wavelength_list'].size):
                data = trust_grain_file.readline().split()
                self.parameter['wavelength_list'][i_wl] = float(data[1]) * 1e-6
                self.parameter['Qabs'][i_dust, i_wl] = float(data[2])
                self.parameter['Qsca'][i_dust, i_wl] = float(data[3])
                self.parameter['Qext'][i_dust, i_wl] = float(data[4])
                self.parameter['g'][i_dust, i_wl, 0] = float(data[5])
            trust_grain_file.readline()
        trust_grain_file.close()

    def read_trust_calorimetry(self):
        if 'PAH' in self.parameter['input_file'] or 'Gra' in self.parameter['input_file']:
            filename = 'Graphitic_Calorimetry_1000.dat'
        elif 'Sil' in self.parameter['input_file']:
            filename = 'Silicate_Calorimetry_1000.dat'
        else:
            raise ValueError('Dust component name: ' + self.parameter['input_file'] + ' is not known!')
        trust_calorimetry_data = np.genfromtxt(self.file_io.path['trust'] + 'GrainInputs/' + filename, skip_header=3)
        # Set number of temperature for heat capacity
        self.parameter['temp_list'] = np.zeros(trust_calorimetry_data.shape[0])
        # Update the calorimetry array to fit with the temperatures
        self.init_calorimetries(self.parameter['temp_list'].size, 1)
        # Set type of the calorimetry
        self.parameter['calorimetry_type'] = 'enthalpy'
        # Set the temperatures and heat capacities
        for i_temp, temp in enumerate(self.parameter['temp_list']):
            self.parameter['temp_list'][i_temp] = trust_calorimetry_data[i_temp, 0]
            # from [erg/gm] to [J/m^3]: "1e-7 * 1e3 * material_density"
            self.parameter['calorimetry_list'][0, i_temp] = trust_calorimetry_data[i_temp, 1] * 1e-7 * 1e3 * \
                self.parameter['material_density']

    def read_trust_scattering(self):
        for i_theta in range(0, 181):
            filename = 'ZDA_BARE_GR_S_ESM_' + str(i_theta).zfill(3) + 'deg.dat'
            trust_scattering_matrix = np.genfromtxt(self.file_io.path['trust'] + 'ScatMatrix/' + filename,
                skip_header=20)
            for i_a, a_eff in enumerate(self.parameter['size_list']):
                for i_wl, wavelength in enumerate(self.parameter['wavelength_list']):
                    self.parameter['scat_matrix'][i_a, i_wl, 0, 0, i_theta, :] = trust_scattering_matrix[i_wl, 1:]

    def write_dust_to_file(self, description):
        dust_file = open(self.file_io.path['dust'] + self.parameter['dust_cat_file'] + '.dat', 'w')
        dust_file.write('#Robert Brauer (email: robert.brauer@cea.fr)\n')
        dust_file.write('\n')
        dust_file.write('#Post-Doc at CEA Saclay\n')
        dust_file.write('#DRF / IRFU / Service d\'Astrophysique\n')
        dust_file.write('#Orme des Merisiers, Bât 709\n')
        dust_file.write('#91191 Gif sur Yvette\n')
        dust_file.write('#France\n')
        import datetime
        dust_file.write('#History:   ' + str(datetime.date.today()) + '\n')
        dust_file.write('\n')
        dust_file.write('#string ID\n')
        dust_file.write(str(self.parameter['dust_cat_file'] + description) + '\n')
        dust_file.write('\n')
        dust_file.write(
            '#nr. of dust species #wavelength #inc. angles #aspect ratio #density [kg/m^3] #sub.temp #delta  #align\n')
        dust_file.write(
            str(self.parameter['size_list'].size) + '\t' +
            str(self.parameter['wavelength_list'].size) + '\t' +
            str(self.parameter['inc_angle_list'].size) + '\t' +
            str(self.parameter['aspect_ratio']) + '\t' +
            str(self.parameter['material_density']) + '\t' +
            str(self.parameter['subl_temp']) + '\t' +
            str(self.parameter['rat_delta']) + '\t' +
            str(self.parameter['align']) + '\n')
        dust_file.write('\n')
        dust_file.write('#a_eff\n')
        dust_file.write('#')
        for i_dust, size in enumerate(self.parameter['size_list']):
            if i_dust < self.parameter['size_list'].size - 1:
                dust_file.write(str(i_dust) + '\t')
            else:
                dust_file.write(str(i_dust) + '\n')
        for i_dust, size in enumerate(self.parameter['size_list']):
            if i_dust < self.parameter['size_list'].size - 1:
                dust_file.write(str(self.parameter['size_list'][i_dust]) + '\t')
            else:
                dust_file.write(str(self.parameter['size_list'][i_dust]) + '\n')
        dust_file.write('\n')
        dust_file.write('#wavelength\n')
        dust_file.write('#')
        for i_wl in range(self.parameter['wavelength_list'].size):
            if i_wl < self.parameter['wavelength_list'].size - 1:
                dust_file.write(str(i_wl) + '\t')
            else:
                dust_file.write(str(i_wl) + '\n')
        for i_wl in range(self.parameter['wavelength_list'].size):
            if i_wl < self.parameter['wavelength_list'].size - 1:
                dust_file.write(str(self.parameter['wavelength_list'][i_wl]) + '\t')
            else:
                dust_file.write(str(self.parameter['wavelength_list'][i_wl]) + '\n')
        dust_file.write('\n')
        dust_file.write('#')
        nr_quantities = 7 + self.parameter['inc_angle_list'].size * 2
        for i_quantity in range(nr_quantities):
            if i_quantity < nr_quantities - 1:
                dust_file.write(str(i_quantity) + '\t')
            else:
                dust_file.write(str(i_quantity) + '\n')
        dust_file.write('#Qext1\tQext2\tQabs1\tQabs2\tQsca1\tQsca2\tdQphas\t')
        for i_angle in range(self.parameter['inc_angle_list'].size):
            dust_file.write('Qtrq' + str(int(i_angle * 90 / self.parameter['inc_angle_list'].size)) + '\t')
        for i_angle in range(self.parameter['inc_angle_list'].size):
            if i_angle < self.parameter['inc_angle_list'].size - 1:
                dust_file.write('g' + str(int(i_angle * 90 / self.parameter['inc_angle_list'].size)) + '\t')
            else:
                dust_file.write('g' + str(int(i_angle * 90 / self.parameter['inc_angle_list'].size)) + '\n')
        for i_wl in range(self.parameter['wavelength_list'].size):
            for i_dust, size in enumerate(self.parameter['size_list']):
                dust_file.write(str(self.parameter['Qext'][i_dust, i_wl]) + '\t' +
                                str(self.parameter['Qext'][i_dust, i_wl]) + '\t' +
                                str(self.parameter['Qabs'][i_dust, i_wl]) + '\t' +
                                str(self.parameter['Qabs'][i_dust, i_wl]) + '\t' +
                                str(self.parameter['Qsca'][i_dust, i_wl]) + '\t' +
                                str(self.parameter['Qsca'][i_dust, i_wl]) + '\t' +
                                str(self.parameter['d_Qphas'][i_dust, i_wl]) + '\t')
                for i_angle in range(self.parameter['inc_angle_list'].size):
                    dust_file.write(str(self.parameter['Qtrq'][i_dust, i_wl, i_angle]) + '\t')
                for i_angle in range(self.parameter['inc_angle_list'].size):
                    if i_angle < self.parameter['inc_angle_list'].size - 1:
                        dust_file.write(str(self.parameter['g'][i_dust, i_wl, i_angle]) + '\t')
                    else:
                        dust_file.write(str(self.parameter['g'][i_dust, i_wl, i_angle]) + '\n')
        dust_file.close()

    def write_dust_scat_to_file_init(self):
        os.chdir(self.file_io.path['dust'])
        if not os.path.isdir(self.parameter['dust_cat_file'] + '/'):
            os.mkdir(self.parameter['dust_cat_file'] + '/')
        scat_info_file = open(self.parameter['dust_cat_file'] + '/scat.inf', 'w')
        scat_info_file.write(str(self.parameter['size_list'].size) + '\t' +
                             str(self.parameter['wavelength_list'].size) + '\t' +
                             str(self.parameter['inc_angle_list'].size) + '\t' +
                             str(self.parameter['phi_angle_list'].size) + '\t' +
                             str(self.parameter['theta_angle_list'].size) + '\n')
        scat_info_file.write(str(self.nr_matrix_elements) + '\n')
        elements_string = str()
        for i, element in enumerate(self.elements):
            if i == self.elements.size - 1:
                elements_string += str(element) + '\n'
            else:
                elements_string += str(element) + '\t'
        scat_info_file.write(elements_string)
        scat_info_file.close()

    def write_dust_scat_to_file(self):
        import struct
        os.chdir(self.file_io.path['dust'])
        for i_wl, wl in enumerate(self.parameter['wavelength_list']):
            scat_bin_file = open(self.parameter['dust_cat_file'] + '/wID' + str(i_wl + 1).zfill(3) + '.sca', 'wb')
            for i_dust, size in enumerate(self.parameter['size_list']):
                for i_inc, inc in enumerate(self.parameter['inc_angle_list']):
                    for i_ph, phi in enumerate(self.parameter['phi_angle_list']):
                        for i_th, theta in enumerate(self.parameter['theta_angle_list']):
                            for i_mat in range(self.nr_matrix_elements):
                                scat_bin_file.write(struct.pack('f', self.parameter['scat_matrix'][i_dust, i_wl,
                                    i_inc, i_ph, i_th, i_mat]))
            scat_bin_file.close()

    def write_dust_calorimetry_to_file(self):
        os.chdir(self.file_io.path['dust'])
        if not os.path.isdir(self.parameter['dust_cat_file'] + '/'):
            os.mkdir(self.parameter['dust_cat_file'] + '/')
        calorimetry_data_file = open(self.parameter['dust_cat_file'] + '/calorimetry.dat', 'w')
        calorimetry_data_file.write('# nr. of temperatures\n')
        calorimetry_data_file.write(str(self.parameter['temp_list'].size) + '\n')
        line_string = str()
        for i_temp, temp in enumerate(self.parameter['temp_list']):
            if i_temp < self.parameter['temp_list'].size - 1:
                line_string += str(self.parameter['temp_list'][i_temp]) + '\t'
            else:
                line_string += str(self.parameter['temp_list'][i_temp]) + '\n'
        calorimetry_data_file.write('# temperature: T [K]\n')
        calorimetry_data_file.write(line_string)
        calorimetry_data_file.write('# type of calorimetry\n')
        if self.parameter['calorimetry_type'] == 'heat_capacity':
            calorimetry_data_file.write('0\n')
            calorimetry_data_file.write('# heat capacity C [J/K/m^3]\n')
            sign = 'C'
        elif self.parameter['calorimetry_type'] == 'enthalpy':
            calorimetry_data_file.write('1\n')
            calorimetry_data_file.write('# enthalpy H [J/m^3]\n')
            sign = 'H'
        if self.parameter['calorimetry_list'].shape[0] == 1:
            calorimetry_data_file.write('# ' + sign + '(T_0)\n')
            calorimetry_data_file.write('# ' + sign + '(T_1)\n')
            calorimetry_data_file.write('# ...\n')
            for i_temp, temp in enumerate(self.parameter['temp_list']):
                calorimetry_data_file.write(str(self.parameter['calorimetry_list'][0, i_temp]) + '\n')
        else:
            calorimetry_data_file.write('# ' + sign + '(T_0, a_0), ' + sign +
                '(T_0, a_1), ' + sign + '(T_0, a_2), ... \n')
            calorimetry_data_file.write('# ' + sign + '(T_1, a_0), ' + sign +
                '(T_1, a_1), ' + sign + '(T_1, a_2), ... \n')
            for i_temp, temp in enumerate(self.parameter['temp_list']):
                line_string = str()
                for i_dust, size in enumerate(self.parameter['size_list']):
                    if i_dust < self.parameter['size_list'].size - 1:
                        line_string += str(self.parameter['calorimetry_list'][i_dust, i_temp]) + '\t'
                    else:
                        line_string += str(self.parameter['calorimetry_list'][i_dust, i_temp]) + '\n'
                calorimetry_data_file.write(line_string)
        calorimetry_data_file.close()

    def init_optical_properties(self):
        self.parameter['Qsca'] = np.zeros((self.parameter['size_list'].size,
            self.parameter['wavelength_list'].size))
        self.parameter['Qabs'] = np.zeros((self.parameter['size_list'].size,
            self.parameter['wavelength_list'].size))
        self.parameter['Qext'] = np.zeros((self.parameter['size_list'].size,
            self.parameter['wavelength_list'].size))
        self.parameter['d_Qphas'] = np.zeros((self.parameter['size_list'].size,
            self.parameter['wavelength_list'].size))
        self.parameter['Qtrq'] = np.zeros((self.parameter['size_list'].size,
            self.parameter['wavelength_list'].size, 1))
        self.parameter['g'] = np.zeros((self.parameter['size_list'].size,
            self.parameter['wavelength_list'].size, 1))
        self.parameter['scat_matrix'] = np.zeros((self.parameter['size_list'].size,
            self.parameter['wavelength_list'].size,
            self.parameter['inc_angle_list'].size,
            self.parameter['phi_angle_list'].size,
            self.parameter['theta_angle_list'].size,
            self.nr_matrix_elements))

    def init_calorimetries(self, nr_of_temperatures, nr_of_dust_species=None):
        self.parameter['temp_list'] = np.zeros(nr_of_temperatures)
        if nr_of_dust_species is not None:
            self.parameter['calorimetry_list'] = np.zeros((nr_of_dust_species, self.parameter['temp_list'].size))
        else:
            self.parameter['calorimetry_list'] = np.zeros((self.parameter['size_list'].size,
                self.parameter['temp_list'].size))

    def add_size(self, size):
        from scipy.interpolate import InterpolatedUnivariateSpline

        if size < self.parameter['size_list'][0]:
            index = 0
        elif size >= sel000f.parameter['size_list'][-1]:
            index = self.parameter['size_list'].size
        else:
            for i_dust, size in enumerate(self.parameter['size_list']):
                if self.parameter['size_list'][i_dust] >= size:
                    index = i_dust
                    break

        parameter_tmp = self.parameter.copy()
        for quantity in ['Qabs', 'Qsca', 'Qext', 'd_Qphas', 'Qtrq', 'g', 'scat_matrix']:
            self.parameter[quantity] = np.insert(self.parameter[quantity], index, 0, axis=0)

        for i_wl in range(self.parameter['wavelength_list'].size):
            for quantity in ['Qabs', 'Qsca', 'Qext']:
                s = InterpolatedUnivariateSpline(self.parameter['size_list'], parameter_tmp[quantity][:, i_wl], k=2)
                self.parameter[quantity][index, i_wl] = s(size)
            s = InterpolatedUnivariateSpline(self.parameter['size_list'], parameter_tmp['g'][:, i_wl, 0], k=2)
            self.parameter['g'][index, i_wl, 0] = s(size)
        if not np.all(self.parameter['scat_matrix'] == 0):
            for i_wl, wl in enumerate(self.parameter['wavelength_list']):
                for i_inc, inc in enumerate(self.parameter['inc_angle_list']):
                    for i_ph, phi in enumerate(self.parameter['phi_angle_list']):
                        for i_th, theta in enumerate(self.parameter['theta_angle_list']):
                            for i_mat in range(self.nr_matrix_elements):
                                s = InterpolatedUnivariateSpline(self.parameter['size_list'],
                                    parameter_tmp['scat_matrix'][:, i_wl, i_inc, i_ph, i_th, i_mat], k=2)
                                self.parameter['scat_matrix'][index, i_wl, i_inc, i_ph, i_th, i_mat] = s(size)

        self.parameter['size_list'] = np.insert(self.parameter['size_list'], index, size)
