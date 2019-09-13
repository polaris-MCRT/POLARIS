#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os


class GasCreator:
    """The GasCreator class is able to create various gas species for POLARIS.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        self.file_io = file_io
        self.parse_args = parse_args

        # Get math module
        from polaris_tools_modules.math import Math
        self.math = Math()

    def convert_database(self, database, database_code):
        """Converts files fom the cdms/jpl database into the LAMBDA database format.

        Args:
            database (str): Name of the database from which we take the data to convert it.
                (Possible options: jpl, cdms)
            database_code (int): Unique code of the gas species/atom database file.

        Notes:
            JPL: http://spec.jpl.nasa.gov/
            CDMS: http://www.astro.uni-koeln.de/cdms/catalog
            LAMBDA: http://home.strw.leidenuniv.nl/~moldata/
        """
        # Init name of species and partition function with default values
        species_name = ''
        rot_spin_partition_func = 0
        if database == 'cdms':
            # Load partition functions from cdms database
            doc_cdms = open(self.file_io.path[database] + 'c' + str(database_code).zfill(6) + '.out', 'r')
            # Search for the entry with the partition function at T = 300 K
            i_line = 0
            for line in doc_cdms.readlines():
                if '300.000' in line:
                    rot_spin_partition_func = float(line.split()[1])
                if 'ID=' in line:
                    species_name = line.partition(', ID=')[0]
                if i_line == 0:
                    species_name = line.split()[0]
                i_line += 1
            doc_cdms.close()
        elif database == 'jpl':
            # Load partition functions from documentation of the chosen species from jpl database
            doc_jpl = open(self.file_io.path[database] + 'd' + str(database_code).zfill(6) + '.cat', 'r')
            # Search for the entry with the partition function at T = 300 K
            for line in doc_jpl.readlines():
                if 'Q(300.0)=' in line:
                    rot_spin_partition_func = float(line.split()[-1])
                if 'Name:' in line:
                    species_name = line.split()[-1]
            doc_jpl.close()
        else:
            raise ValueError('Database: ' + database + ' not known!')

        # If the partition function was not loaded correctly
        if species_name == '' or rot_spin_partition_func == 0:
            raise EnvironmentError('Name of the species and partition function are not sufficiently loaded!')

        # Init dictionary where each entry is one energy level
        energy_level = dict(
            # Transition index [1/cm]
            ENERGIES=[],
            # Degeneracy of the level
            WEIGHT=[],
            # Related quantum numbers
            QN_1=[],
            QN_2=[],
            QN_3=[],
            QN_4=[],
            QN_5=[],
            QN_6=[],
        )

        if os.path.isfile(self.file_io.path[database] + 'c' + str(database_code).zfill(6) + '.egy'):
            # Load database file
            data = open(self.file_io.path[database] + 'c' + str(database_code).zfill(6) + '.egy', 'r')
            # Get the energy levels from database
            for line in data.readlines():
                formatted_line = line.split()
                # Add the energy level to the dictionary
                energy_level['ENERGIES'].append(abs(float(formatted_line[2])))
                energy_level['WEIGHT'].append(int(formatted_line[5].replace(':', '')))
                if len(formatted_line) >= 7:
                    energy_level['QN_1'].append(int(formatted_line[6]))
                else:
                    energy_level['QN_1'].append(-1)
                if len(formatted_line) >= 8:
                    energy_level['QN_2'].append(int(formatted_line[7]))
                else:
                    energy_level['QN_2'].append(-1)
                if len(formatted_line) >= 9:
                    energy_level['QN_3'].append(int(formatted_line[8]))
                else:
                    energy_level['QN_3'].append(-1)
                if len(formatted_line) >= 10:
                    energy_level['QN_4'].append(int(formatted_line[9]))
                else:
                    energy_level['QN_4'].append(-1)
                if len(formatted_line) >= 11:
                    energy_level['QN_5'].append(int(formatted_line[10]))
                else:
                    energy_level['QN_5'].append(-1)
                if len(formatted_line) >= 12:
                    energy_level['QN_6'].append(int(formatted_line[11]))
                else:
                    energy_level['QN_6'].append(-1)

        else:
            # Load database file
            data = open(self.file_io.path[database] + 'c' + str(database_code).zfill(6) + '.cat', 'r')
            # Per default, the quantum number are integers
            half_value_qn = False
            # Get the upper energy level from each transition
            for line in data.readlines():
                #: float: Energy of lower transition [1/cm]
                energy_lower = float(line[31:41])
                #: float: Energy of upper transition [1/cm]
                energy_upper = energy_lower + float(line[0:13]) * 1e-2 * 1e6 / self.math.const['c']
                # Obtain each available quantum number and calculate the weight from it
                if line[55:57].rstrip('\n') == '  ' or not line[55:57].rstrip('\n'):
                    qn_1 = -1
                else:
                    qn_1 = int(line[55:57].rstrip('\n'))
                # Obtain each available quantum number and calculate the weight from it
                if line[57:59].rstrip('\n') == '  ' or not line[57:59].rstrip('\n'):
                    qn_2 = -1
                else:
                    qn_2 = int(line[57:59].rstrip('\n'))
                # Obtain each available quantum number and calculate the weight from it
                if line[59:61].rstrip('\n') == '  ' or not line[59:61].rstrip('\n'):
                    qn_3 = -1
                else:
                    qn_3 = int(line[59:61].rstrip('\n'))
                # Obtain each available quantum number and calculate the weight from it
                if line[61:63].rstrip('\n') == '  ' or not line[61:63].rstrip('\n'):
                    qn_4 = -1
                else:
                    qn_4 = int(line[61:63].rstrip('\n'))
                if line[63:65].rstrip('\n') == '  ' or not line[63:65].rstrip('\n'):
                    qn_5 = -1
                else:
                    qn_5 = int(line[63:65].rstrip('\n'))
                if line[65:67].rstrip('\n') == '  ' or not line[65:67].rstrip('\n'):
                    qn_6 = -1
                else:
                    qn_6 = int(line[65:67].rstrip('\n'))
                weight = self.math.get_weight_from_qn(database_code, qn_1, qn_2, qn_3, qn_4, qn_5, qn_6)
                # For the upper energy level, if the calculated weight is not same as the weight of the transition,
                # the quantum numbers are not integers but .5 values.
                if weight != int(line[41:44]):
                    if weight - int(line[41:44]) == 1:
                        half_value_qn = True
                    else:
                        print(weight, int(line[41:44]))
                        raise ValueError('weight from quantum numbers does not fit with weight of spectral transition!')
                # Add the energy level to the dictionary
                energy_level['ENERGIES'].append(abs(energy_upper))
                energy_level['QN_1'].append(int(qn_1))
                energy_level['QN_2'].append(int(qn_2))
                energy_level['QN_3'].append(int(qn_3))
                energy_level['QN_4'].append(int(qn_4))
                energy_level['QN_5'].append(int(qn_5))
                energy_level['QN_6'].append(int(qn_6))
                # Adjust the weight according to the format of the quantum numbers
                if half_value_qn:
                    energy_level['WEIGHT'].append(weight - 1)
                else:
                    energy_level['WEIGHT'].append(weight)
            data.close()

            # Load database file
            data = open(self.file_io.path[database] + 'c' + str(database_code).zfill(6) + '.cat', 'r')
            # Get the lower energy level from each transition
            for line in data.readlines():
                #: float: Energy of lower transition [1/cm]
                energy_lower = float(line[31:41])
                # Obtain each available quantum number and calculate the weight from it
                if line[67:69].rstrip('\n') == '  ' or not line[67:69].rstrip('\n'):
                    qn_1 = -1
                else:
                    qn_1 = int(line[67:69].rstrip('\n'))
                # Obtain each available quantum number and calculate the weight from it
                if line[69:71].rstrip('\n') == '  ' or not line[69:71].rstrip('\n'):
                    qn_2 = -1
                else:
                    qn_2 = int(line[69:71].rstrip('\n'))
                # Obtain each available quantum number and calculate the weight from it
                if line[71:73].rstrip('\n') == '  ' or not line[71:73].rstrip('\n'):
                    qn_3 = -1
                else:
                    qn_3 = int(line[71:73].rstrip('\n'))
                # Obtain each available quantum number and calculate the weight from it
                if line[73:75].rstrip('\n') == '  ' or not line[73:75].rstrip('\n'):
                    qn_4 = -1
                else:
                    qn_4 = int(line[73:75].rstrip('\n'))
                if line[73:75].rstrip('\n') == '  ' or not line[73:75].rstrip('\n'):
                    qn_5 = -1
                else:
                    qn_5 = int(line[73:75].rstrip('\n'))
                if line[75:77].rstrip('\n') == '  ' or not line[75:77].rstrip('\n'):
                    qn_6 = -1
                else:
                    qn_6 = int(line[75:77].rstrip('\n'))
                weight = self.math.get_weight_from_qn(database_code, qn_1, qn_2, qn_3, qn_4, qn_5, qn_6)
                # Add the energy level to the dictionary
                energy_level['ENERGIES'].append(abs(energy_lower))
                energy_level['QN_1'].append(int(qn_1))
                energy_level['QN_2'].append(int(qn_2))
                energy_level['QN_3'].append(int(qn_3))
                energy_level['QN_4'].append(int(qn_4))
                energy_level['QN_5'].append(int(qn_5))
                energy_level['QN_6'].append(int(qn_6))
                # Adjust the weight according to the format of the quantum numbers
                if half_value_qn:
                    energy_level['WEIGHT'].append(weight - 1)
                else:
                    energy_level['WEIGHT'].append(weight)
            data.close()

            # Algorithm to sort the energy levels with their energy
            for i_level in range(len(energy_level['ENERGIES'])):
                print(i_level,"/",len(energy_level['ENERGIES']))
                for j_level in range(i_level + 1, len(energy_level['ENERGIES'])):
                    # If the next energy level has a higher energy than the current one, switch them
                    if energy_level['ENERGIES'][i_level] > energy_level['ENERGIES'][j_level]:
                        # Switch each entry of the dictionary
                        for keys in energy_level.keys():
                            tmp = energy_level[keys][j_level]
                            energy_level[keys][j_level] = energy_level[keys][i_level]
                            energy_level[keys][i_level] = tmp

            # Delete duplicates of the energy levels
            index = len(energy_level['ENERGIES']) - 1
            while index != 0:
                index -= 1
                # If all quantum numbers are the same, the energy level is the same
                # However, the energy value might vary slightly
                if energy_level['QN_1'][index] == energy_level['QN_1'][index + 1] \
                        and energy_level['QN_2'][index] == energy_level['QN_2'][index + 1] \
                        and energy_level['QN_3'][index] == energy_level['QN_3'][index + 1] \
                        and energy_level['QN_4'][index] == energy_level['QN_4'][index + 1] \
                        and energy_level['QN_5'][index] == energy_level['QN_5'][index + 1] \
                        and energy_level['QN_6'][index] == energy_level['QN_6'][index + 1]:
                    del energy_level['ENERGIES'][index + 1]
                    del energy_level['WEIGHT'][index + 1]
                    del energy_level['QN_1'][index + 1]
                    del energy_level['QN_2'][index + 1]
                    del energy_level['QN_3'][index + 1]
                    del energy_level['QN_4'][index + 1]
                    del energy_level['QN_5'][index + 1]
                    del energy_level['QN_6'][index + 1]

        # Create a new entry for the index of the energy levels
        energy_level['LEVEL'] = []
        for i_level in range(len(energy_level['ENERGIES'])):
            energy_level['LEVEL'].append(i_level + 1)

        # Init dictionary where each entry is one spectral line transition
        radiative_transitions = dict(
            # Transition index
            TRANS=[],
            # Index of upper transition
            UP=[],
            # Index of lower transition
            LOW=[],
            # Einstein coefficient A [1/s]
            EINSTEINA=[],
            # Frequency of transition [GHz]
            FREQ=[],
            # Energy of transition [K]
            E=[],
        )

        # Load database file
        data = open(self.file_io.path[database] + 'c' + str(database_code).zfill(6) + '.cat', 'r')
        # Go through each transition again
        i_trans = 0
        for line in data.readlines():
            i_trans += 1
            # Set index of transition
            radiative_transitions['TRANS'].append(i_trans)
            #: float: Energy of lower transition [1/cm]
            energy_lower = float(line[31:41])
            #: float: Energy of upper transition [1/cm]
            energy_upper = energy_lower + float(line[0:13]) * 1e-2 * 1e6 / self.math.const['c']
            # Find index of lower energy level
            if energy_lower == 0.:
                # If the lower energy level is the zero level, it has the index 1
                radiative_transitions['LOW'].append(1)
            else:
                # Init comparison values
                min_value_lower = 999.9
                min_index_lower = 0
                # The lower energy level is found by taking the one with the closest energy
                for i_level in range(len(energy_level['LEVEL'])):
                    if abs(energy_level['ENERGIES'][i_level] - energy_lower) < min_value_lower:
                        min_value_lower = abs(energy_level['ENERGIES'][i_level] - energy_lower)
                        min_index_lower = i_level
                # Set index of lower transition
                radiative_transitions['LOW'].append(energy_level['LEVEL'][min_index_lower])
            # Init comparison values
            min_value_upper = 999.9
            min_index_upper = 0
            # The upper energy level is found by taking the one with the closest energy
            for i_level in range(len(energy_level['LEVEL'])):
                if abs(energy_level['ENERGIES'][i_level] - energy_upper) < min_value_upper:
                    min_value_upper = abs(energy_level['ENERGIES'][i_level] - energy_upper)
                    min_index_upper = i_level
            # Set index of upper transition
            radiative_transitions['UP'].append(energy_level['LEVEL'][min_index_upper])
            #: float: Frequency of lower transition [Hz]
            freq_lower = energy_lower * 1e2 * self.math.const['c']
            #: float: Frequency of upper transition [Hz]
            freq_upper = freq_lower + float(line[0:13]) * 1e6
            #: float: Energy of lower transition
            energy_lower = freq_lower * self.math.const['h']
            exp_val_lower = np.exp(-energy_lower / (self.math.const['k_B'] * 300))
            #: float: Energy of upper transition
            energy_upper = freq_upper * self.math.const['h']
            exp_val_upper = np.exp(-energy_upper / (self.math.const['k_B'] * 300))
            #: float: Line intensity at 300K
            line_int = 10 ** float(line[21:29])
            #: float: Squared value of the transition frequency [MHz]
            trans_freq_sq = float(line[0:13]) ** 2
            #: float: Upper state degeneracy
            upper_deg = float(line[41:44])
            #: float: Einstein coefficient A
            einstein_a = line_int * trans_freq_sq * (rot_spin_partition_func / upper_deg) * 1 / (
                    exp_val_lower - exp_val_upper) * 2.7964e-16
            # Set einstein coefficient A
            radiative_transitions['EINSTEINA'].append(einstein_a)
            # Set frequency of transition [GHz]
            radiative_transitions['FREQ'].append(float(line[0:13]) * 1e-3)
            #: float: Transition frequency [Hz]
            trans_freq = float(line[0:13]) * 1e6
            #: float: Brightness temperature
            bright_temp = self.math.const['h'] * trans_freq / self.math.const['k_B']
            # Set Energy of transition in Kelvin
            radiative_transitions['E'].append(bright_temp)
        data.close()

        # Write the data as Leiden database file
        leiden_file = open(self.file_io.path['gas'] + species_name + '_' + database + '.dat', 'w')
        leiden_file.write('!GAS_SPECIES\n')
        leiden_file.write(species_name + '\n')
        leiden_file.write('!MOLECULAR WEIGHT\n')
        leiden_file.write(str(database_code).zfill(6)[0:3].upper() + '\n')
        leiden_file.write('!NUMBER OF ENERGY LEVELS\n')
        leiden_file.write(str(len(energy_level['LEVEL'])) + '\n')
        leiden_file.write('!LEVEL + ENERGIES[cm^-1] + WEIGHT + QNUM\n')
        for i_level in range(len(energy_level['LEVEL'])):
            qn_string = ''
            if energy_level['QN_1'][i_level] != -1:
                qn_string += str(energy_level['QN_1'][i_level]).rjust(2)
            if energy_level['QN_2'][i_level] != -1:
                qn_string += str(energy_level['QN_2'][i_level]).rjust(2)
            if energy_level['QN_3'][i_level] != -1:
                qn_string += str(energy_level['QN_3'][i_level]).rjust(2)
            if energy_level['QN_4'][i_level] != -1:
                qn_string += str(energy_level['QN_4'][i_level]).rjust(2) + '    '
            if energy_level['QN_5'][i_level] != -1:
                qn_string += str(energy_level['QN_5'][i_level]).rjust(2) + '    '
            if energy_level['QN_6'][i_level] != -1:
                qn_string += str(energy_level['QN_6'][i_level]).rjust(2) + '    '

            leiden_file.write('{:5d}'.format(energy_level['LEVEL'][i_level])
                              + '{:16.8f}'.format(energy_level['ENERGIES'][i_level])
                              + '{:9.1f}'.format(energy_level['WEIGHT'][i_level])
                              + '     ' + qn_string + '\n')
        leiden_file.write('!NUMBER OF RADIATIVE TRANSITIONS\n')
        leiden_file.write(str(len(radiative_transitions['TRANS'])) + '\n')
        leiden_file.write('!TRANS + UP + LOW + EINSTEINA[s^-1] + FREQ[GHz] + E_u[K]\n')
        for i_trans in range(len(radiative_transitions['TRANS'])):
            leiden_file.write('{:5d}'.format(radiative_transitions['TRANS'][i_trans])
                              + '{:6d}'.format(radiative_transitions['UP'][i_trans])
                              + '{:6d}'.format(radiative_transitions['LOW'][i_trans])
                              + '{:12.3e}'.format(radiative_transitions['EINSTEINA'][i_trans])
                              + '{:16.7f}'.format(radiative_transitions['FREQ'][i_trans])
                              + '{:10.2f}'.format(radiative_transitions['E'][i_trans]) + '\n')
        leiden_file.close()

    def create_zeeman_file(self):
        """Create Zeeman database file which can be used by POLARIS.
        """
        if self.parse_args.gas_species.upper() == 'CN':
            j_upper_list = [0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 1.5]
            f_upper_list = [0.5, 1.5, 1.5, 1.5, 2.5, 0.5, 1.5]
            j_lower_list = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
            f_lower_list = [1.5, 0.5, 1.5, 0.5, 1.5, 0.5, 1.5]
            transition_index = [11, 12, 13, 14, 15, 16, 17]

            zeeman_file = open(self.file_io.path['gas'] + 'cn_zeeman.dat', 'w')
            zeeman_file.write('!Gas species name\n')
            zeeman_file.write('CN\n')
            zeeman_file.write('!Gas species radius for collision calculations\n')
            # The gas_species radius is calculated by adding the covalent
            # radii of all atoms together.
            zeeman_file.write(str(self.math.covalent_radii['C'] +
                                  self.math.covalent_radii['N']) + '\n')
            zeeman_file.write('!Number of transitions with Zeeman effect\n')
            zeeman_file.write(str(len(j_upper_list)) + '\n')

            zeeman_shifts = np.zeros(len(j_upper_list))
            for i_trans in range(len(j_upper_list)):
                lande_upper = self.math.lande_g_cn(1, j_upper_list[i_trans], f_upper_list[i_trans])
                lande_lower = self.math.lande_g_cn(0, j_lower_list[i_trans], f_lower_list[i_trans])
                zeeman_file.write('!Transition index in Leiden database\n')
                zeeman_file.write(str(transition_index[i_trans]) + '\n')
                zeeman_file.write('!Lande factor of upper level\n')
                zeeman_file.write(str(lande_upper) + '\n')
                zeeman_file.write('!Lande factor of lower level\n')
                zeeman_file.write(str(lande_lower) + '\n')
                zeeman_file.write('!Number of Zeeman sublevels in the upper level\n')
                zeeman_file.write(str((f_upper_list[i_trans] * 2) + 1) + '\n')
                zeeman_file.write('!Number of Zeeman sublevels in the lower level\n')
                zeeman_file.write(str((f_lower_list[i_trans] * 2) + 1) + '\n')

                if f_lower_list[i_trans] == f_upper_list[i_trans]:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_zero(f, m_f, transition)
                elif (f_lower_list[i_trans] - f_upper_list[i_trans]) == +1:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_plus(f, m_f, transition)
                elif (f_lower_list[i_trans] - f_upper_list[i_trans]) == -1:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_minus(f, m_f, transition)

                for i_upper_trans in np.arange(-float(f_upper_list[i_trans]), float(f_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(f_lower_list[i_trans]), float(f_lower_list[i_trans] + 1), 1):
                        if i_lower_trans == i_upper_trans:
                            zeeman_file.write('!Line strength of pi transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans) + ')\n')
                            rel_str_pi = 2.0 * relative_strength(f_upper_list[i_trans], i_upper_trans, 0)
                            zeeman_file.write(str(rel_str_pi) + '\n')
                for i_upper_trans in np.arange(-float(f_upper_list[i_trans]), float(f_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(f_lower_list[i_trans]), float(f_lower_list[i_trans] + 1), 1):
                        if i_lower_trans - i_upper_trans == +1:
                            zeeman_file.write('!Line strength of sigma_+ transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans + 1) + ')\n')
                            rel_str_sigma_p = relative_strength(f_upper_list[i_trans], i_upper_trans, +1)
                            zeeman_file.write(str(rel_str_sigma_p) + '\n')
                for i_upper_trans in np.arange(-float(f_upper_list[i_trans]), float(f_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(f_lower_list[i_trans]), float(f_lower_list[i_trans] + 1), 1):
                        if i_lower_trans - i_upper_trans == -1:
                            zeeman_file.write('!Line strength of sigma_- transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans - 1) + ')\n')
                            rel_str_sigma_m = relative_strength(f_upper_list[i_trans], i_upper_trans, -1)
                            zeeman_file.write(str(rel_str_sigma_m) + '\n')
                            zeeman_shifts[i_trans] += 2.0 * rel_str_sigma_m * (
                                    (lande_upper * i_upper_trans) - (lande_lower * i_lower_trans))
            zeeman_file.close()
            zeeman_shifts = np.multiply(zeeman_shifts, self.math.const['mu_B'] / self.math.const['h'] * 1e-10 * 2.0)
            print('2*nu/B =', zeeman_shifts[0], 'Hz/muG')
        elif self.parse_args.gas_species.upper() == 'SO':
            n_upper_list = [2, 2, 3, 4, 5]
            j_upper_list = [3, 1, 4, 3, 6]
            n_lower_list = [1, 1, 2, 3, 4]
            j_lower_list = [2, 2, 3, 2, 5]
            transition_index = [3, 8, 10, 28, 34]

            zeeman_file = open(self.file_io.path['gas'] + 'so_zeeman.dat', 'w')
            zeeman_file.write('!Gas species name\n')
            zeeman_file.write('SO\n')
            zeeman_file.write('!Gas species radius for collision calculations\n')
            # The gas_species radius is calculated by adding the covalent
            # radii of all atoms together.
            zeeman_file.write(str(self.math.covalent_radii['S'] +
                                  self.math.covalent_radii['O']) + '\n')
            zeeman_file.write('!Number of transitions with Zeeman effect\n')
            zeeman_file.write(str(len(j_upper_list)) + '\n')

            zeeman_shifts = np.zeros(len(j_upper_list))
            for i_trans in range(len(j_upper_list)):
                lande_upper = self.math.lande_g_so(n_upper_list[i_trans], j_upper_list[i_trans])
                lande_lower = self.math.lande_g_so(n_lower_list[i_trans], j_lower_list[i_trans])
                zeeman_file.write('!Transition index in Leiden database\n')
                zeeman_file.write(str(transition_index[i_trans]) + '\n')
                zeeman_file.write('!Lande factor of upper level\n')
                zeeman_file.write(str(lande_upper) + '\n')
                zeeman_file.write('!Lande factor of lower level\n')
                zeeman_file.write(str(lande_lower) + '\n')
                zeeman_file.write('!Number of Zeeman sublevels in the upper level\n')
                zeeman_file.write(str((j_upper_list[i_trans] * 2) + 1) + '\n')
                zeeman_file.write('!Number of Zeeman sublevels in the lower level\n')
                zeeman_file.write(str((j_lower_list[i_trans] * 2) + 1) + '\n')

                if j_lower_list[i_trans] == j_upper_list[i_trans]:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_zero(f, m_f, transition)
                elif (j_lower_list[i_trans] - j_upper_list[i_trans]) == +1:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_plus(f, m_f, transition)
                elif (j_lower_list[i_trans] - j_upper_list[i_trans]) == -1:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_minus(f, m_f, transition)

                for i_upper_trans in np.arange(-float(j_upper_list[i_trans]), float(j_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(j_lower_list[i_trans]), float(j_lower_list[i_trans] + 1), 1):
                        if i_lower_trans == i_upper_trans:
                            zeeman_file.write('!Line strength of pi transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans) + ')\n')
                            rel_str_pi = 2.0 * relative_strength(j_upper_list[i_trans], i_upper_trans, 0)
                            zeeman_file.write(str(rel_str_pi) + '\n')
                for i_upper_trans in np.arange(-float(j_upper_list[i_trans]), float(j_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(j_lower_list[i_trans]), float(j_lower_list[i_trans] + 1), 1):
                        if i_lower_trans - i_upper_trans == +1:
                            zeeman_file.write('!Line strength of sigma_+ transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans + 1) + ')\n')
                            rel_str_sigma_p = relative_strength(j_upper_list[i_trans], i_upper_trans, +1)
                            zeeman_file.write(str(rel_str_sigma_p) + '\n')
                for i_upper_trans in np.arange(-float(j_upper_list[i_trans]), float(j_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(j_lower_list[i_trans]), float(j_lower_list[i_trans] + 1), 1):
                        if i_lower_trans - i_upper_trans == -1:
                            zeeman_file.write('!Line strength of sigma_- transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans - 1) + ')\n')
                            rel_str_sigma_m = relative_strength(j_upper_list[i_trans], i_upper_trans, -1)
                            zeeman_file.write(str(rel_str_sigma_m) + '\n')
                            zeeman_shifts[i_trans] += 2.0 * rel_str_sigma_m * (
                                    (lande_upper * i_upper_trans) - (lande_lower * i_lower_trans))
            zeeman_file.close()
            zeeman_shifts = np.multiply(zeeman_shifts, self.math.const['mu_B'] / self.math.const['h'] * 1e-10 * 2.0)
            print('2*nu/B =', zeeman_shifts[0], 'Hz/muG')
        elif self.parse_args.gas_species.upper() == 'CCS':
            n_upper_list = [0, 1, 2, 3]
            j_upper_list = [1, 2, 3, 4]
            n_lower_list = [1, 0, 1, 2]
            j_lower_list = [0, 1, 2, 3]
            transition_index = [4, 6, 11, 15]

            zeeman_file = open(self.file_io.path['gas'] + 'ccs_zeeman.dat', 'w')
            zeeman_file.write('!Gas species name\n')
            zeeman_file.write('CCS\n')
            zeeman_file.write('!Gas species radius for collision calculations\n')
            # The gas_species radius is calculated by adding the covalent
            # radii of all atoms together.
            zeeman_file.write(str(2 * self.math.covalent_radii['C'] + self.math.covalent_radii['S']) + '\n')
            zeeman_file.write('!Number of transitions with Zeeman effect\n')
            zeeman_file.write(str(len(j_upper_list)) + '\n')

            zeeman_shifts = np.zeros(len(j_upper_list))
            for i_trans in range(len(j_upper_list)):
                lande_upper = self.math.lande_g_ccs(n_upper_list[i_trans], j_upper_list[i_trans])
                lande_lower = self.math.lande_g_ccs(n_lower_list[i_trans], j_lower_list[i_trans])
                zeeman_file.write('!Transition index in Leiden database\n')
                zeeman_file.write(str(transition_index[i_trans]) + '\n')
                zeeman_file.write('!Lande factor of upper level\n')
                zeeman_file.write(str(lande_upper) + '\n')
                zeeman_file.write('!Lande factor of lower level\n')
                zeeman_file.write(str(lande_lower) + '\n')
                zeeman_file.write('!Number of Zeeman sublevels in the upper level\n')
                zeeman_file.write(str((j_upper_list[i_trans] * 2) + 1) + '\n')
                zeeman_file.write('!Number of Zeeman sublevels in the lower level\n')
                zeeman_file.write(str((j_lower_list[i_trans] * 2) + 1) + '\n')

                if j_lower_list[i_trans] == j_upper_list[i_trans]:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_zero(f, m_f, transition)
                elif (j_lower_list[i_trans] - j_upper_list[i_trans]) == +1:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_plus(f, m_f, transition)
                elif (j_lower_list[i_trans] - j_upper_list[i_trans]) == -1:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_minus(f, m_f, transition)

                for i_upper_trans in np.arange(-float(j_upper_list[i_trans]), float(j_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(j_lower_list[i_trans]), float(j_lower_list[i_trans] + 1), 1):
                        if i_lower_trans == i_upper_trans:
                            zeeman_file.write('!Line strength of pi transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans) + ')\n')
                            rel_str_pi = 2.0 * relative_strength(j_upper_list[i_trans], i_upper_trans, 0)
                            zeeman_file.write(str(rel_str_pi) + '\n')
                for i_upper_trans in np.arange(-float(j_upper_list[i_trans]), float(j_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(j_lower_list[i_trans]), float(j_lower_list[i_trans] + 1), 1):
                        if i_lower_trans - i_upper_trans == +1:
                            zeeman_file.write('!Line strength of sigma_+ transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans + 1) + ')\n')
                            rel_str_sigma_p = relative_strength(j_upper_list[i_trans], i_upper_trans, +1)
                            zeeman_file.write(str(rel_str_sigma_p) + '\n')
                for i_upper_trans in np.arange(-float(j_upper_list[i_trans]), float(j_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(j_lower_list[i_trans]), float(j_lower_list[i_trans] + 1), 1):
                        if i_lower_trans - i_upper_trans == -1:
                            zeeman_file.write('!Line strength of sigma_- transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans - 1) + ')\n')
                            rel_str_sigma_m = relative_strength(j_upper_list[i_trans], i_upper_trans, -1)
                            zeeman_file.write(str(rel_str_sigma_m) + '\n')
                            zeeman_shifts[i_trans] += 2.0 * rel_str_sigma_m * (
                                    (lande_upper * i_upper_trans) - (lande_lower * i_lower_trans))
            zeeman_file.close()
            zeeman_shifts = np.multiply(zeeman_shifts, self.math.const['mu_B'] / self.math.const['h'] * 1e-10 * 2.0)
            print('2*nu/B =', zeeman_shifts[0], 'Hz/muG')
        elif self.parse_args.gas_species.upper() == 'H1':
            l_upper_list = [0]
            f_upper_list = [1]
            l_lower_list = [0]
            f_lower_list = [0]
            transition_index = [1]

            zeeman_file = open(self.file_io.path['gas'] + 'h1_zeeman.dat', 'w')
            zeeman_file.write('!Gas species name\n')
            zeeman_file.write('H1\n')
            zeeman_file.write('!Gas species radius for collision calculations\n')
            # The gas_species radius is calculated by adding the covalent
            # radii of all atoms together.
            zeeman_file.write(str(self.math.covalent_radii['H']) + '\n')
            zeeman_file.write('!Number of transitions with Zeeman effect\n')
            zeeman_file.write(str(len(l_upper_list)) + '\n')

            zeeman_shifts = np.zeros(len(l_upper_list))
            for i_trans in range(len(l_upper_list)):
                lande_upper = self.math.lande_g_h1(l_upper_list[i_trans], f_upper_list[i_trans])
                lande_lower = self.math.lande_g_h1(l_lower_list[i_trans], f_lower_list[i_trans])
                zeeman_file.write('!Transition index in Leiden database\n')
                zeeman_file.write(str(transition_index[i_trans]) + '\n')
                zeeman_file.write('!Lande factor of upper level\n')
                zeeman_file.write(str(lande_upper) + '\n')
                zeeman_file.write('!Lande factor of lower level\n')
                zeeman_file.write(str(lande_lower) + '\n')
                zeeman_file.write('!Number of Zeeman sublevels in the upper level\n')
                zeeman_file.write(str((f_upper_list[i_trans] * 2) + 1) + '\n')
                zeeman_file.write('!Number of Zeeman sublevels in the lower level\n')
                zeeman_file.write(str((f_lower_list[i_trans] * 2) + 1) + '\n')

                if f_lower_list[i_trans] == f_upper_list[i_trans]:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_zero(f, m_f, transition)
                elif (f_lower_list[i_trans] - f_upper_list[i_trans]) == +1:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_plus(f, m_f, transition)
                elif (f_lower_list[i_trans] - f_upper_list[i_trans]) == -1:
                    def relative_strength(f, m_f, transition):
                        return self.math.relative_line_strength_minus(f, m_f, transition)

                for i_upper_trans in np.arange(-float(f_upper_list[i_trans]), float(f_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(f_lower_list[i_trans]), float(f_lower_list[i_trans] + 1), 1):
                        if i_lower_trans == i_upper_trans:
                            zeeman_file.write('!Line strength of pi transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans) + ')\n')
                            rel_str_pi = 2.0 * relative_strength(f_upper_list[i_trans], i_upper_trans, 0)
                            zeeman_file.write(str(rel_str_pi) + '\n')
                for i_upper_trans in np.arange(-float(f_upper_list[i_trans]), float(f_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(f_lower_list[i_trans]), float(f_lower_list[i_trans] + 1), 1):
                        if i_lower_trans - i_upper_trans == +1:
                            zeeman_file.write('!Line strength of sigma_+ transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans + 1) + ')\n')
                            rel_str_sigma_p = relative_strength(f_upper_list[i_trans], i_upper_trans, +1)
                            zeeman_file.write(str(rel_str_sigma_p) + '\n')
                for i_upper_trans in np.arange(-float(f_upper_list[i_trans]), float(f_upper_list[i_trans] + 1), 1):
                    for i_lower_trans in np.arange(-float(f_lower_list[i_trans]), float(f_lower_list[i_trans] + 1), 1):
                        if i_lower_trans - i_upper_trans == -1:
                            zeeman_file.write('!Line strength of sigma_- transition (M\'='
                                              + str(i_upper_trans) + ' -> M\'\'=' + str(i_upper_trans - 1) + ')\n')
                            rel_str_sigma_m = relative_strength(f_upper_list[i_trans], i_upper_trans, -1)
                            zeeman_file.write(str(rel_str_sigma_m) + '\n')
                            zeeman_shifts[i_trans] += 2.0 * rel_str_sigma_m * (
                                    (lande_upper * i_upper_trans) - (lande_lower * i_lower_trans))
            zeeman_file.close()
            zeeman_shifts = np.multiply(zeeman_shifts, self.math.const['mu_B'] / self.math.const['h'] * 1e-10 * 2.0)
            print('2*nu/B =', zeeman_shifts[0], 'Hz/muG')

    def create_hydrogen_file(self):
        """Create LAMBDA database file for H1.
        """
        n_max = 20
        alpha = self.math.const['e'] ** 2 / (
                4 * np.pi * self.math.const['epsilon_0'] * self.math.const['hbar'] * self.math.const['c'])
        energy_level = dict(
            # ID of the energy level
            LEVEL=[],
            # Transition index [1/cm]
            ENERGIES=[],
            # Degeneracy of the level
            WEIGHT=[],
            # Related quantum numbers
            QN=[],
        )
        energy_level_zero = \
            -self.math.const['m_e'] * self.math.const['c'] ** 2 * (
                    1 - (1 + (alpha / (1 - 0.5 - 0.5 + np.sqrt((0.5 + 0.5) ** 2 - alpha ** 2))) ** 2) ** (-0.5))
        gamma = 2.7928
        energy_level_zero += (self.math.const['m_e'] / self.math.const['m_p']) * alpha ** 4 * self.math.const['m_e'] * \
                             self.math.const['c'] ** 2 * 4 * gamma / 3 * (- 3 / 2.)
        i_level = 0
        for n in range(1, n_max + 1):
            j_max = (n - 1) + 0.5
            for j in np.arange(0.5, j_max + 1):
                energy_level_tmp = -self.math.const['m_e'] * self.math.const['c'] ** 2 * (
                        1 - (1 + (alpha / (n - j - 0.5 + np.sqrt((j + 0.5) ** 2 - alpha ** 2))) ** 2) ** (-0.5))
                if n == 1 and j == 0.5:
                    n_extra = 2
                else:
                    n_extra = 1
                for f in range(n_extra):
                    i_level += 1
                    if n_extra == 2:
                        gamma = 2.7928
                        hf_energy = (self.math.const['m_e'] / self.math.const['m_p']) * alpha ** 4 * \
                                    self.math.const['m_e'] * self.math.const['c'] ** 2 * 4 * gamma / \
                                    (3 * n ** 3) * (f * (f + 1) - 3 / 2.)
                        energy_level['QN'].append('{:2d}'.format(n) + '{:2d}'.format(int(j)) + '{:2d}'.format(int(f)))
                    else:
                        hf_energy = 0.
                        energy_level['QN'].append('{:2d}'.format(n) + '{:2d}'.format(int(j)) + ' -')
                    wave_number = ((energy_level_tmp + hf_energy - energy_level_zero) * 1e-2 / (
                            self.math.const['h'] * self.math.const['c']))
                    energy_level['ENERGIES'].append(wave_number)
                    energy_level['WEIGHT'].append(((j * 2) + 1) * 2)

                    energy_level['LEVEL'].append(i_level)

        radiative_transitions = dict(
            # Transition index
            TRANS=[],
            # Index of upper transition
            UP=[],
            # Index of lower transition
            LOW=[],
            # Einstein coefficient A [1/s]
            EINSTEINA=[],
            # Frequency of transition [GHz]
            FREQ=[],
            # Energy of transition [K]
            E=[],
        )
        radiative_transitions['TRANS'].append(1)
        radiative_transitions['LOW'].append(1)
        radiative_transitions['UP'].append(2)
        radiative_transitions['EINSTEINA'].append(2.876e-15)
        radiative_transitions['FREQ'].append(1.4204058)
        radiative_transitions['E'].append(0.07)

        leiden_file = open(self.file_io.path['gas'] + 'H1_converted.dat', 'w')
        leiden_file.write('!GAS SPECIES\n')
        leiden_file.write('H1' + '\n')
        leiden_file.write('!MOLECULAR WEIGHT\n')
        leiden_file.write('1\n')
        leiden_file.write('!NUMBER OF ENERGY LEVELS\n')
        leiden_file.write(str(len(energy_level['LEVEL'])) + '\n')
        leiden_file.write('!LEVEL + ENERGIES[cm^-1] + WEIGHT + QNUM\n')
        for i_level in range(len(energy_level['LEVEL'])):
            leiden_file.write('{:5d}'.format(energy_level['LEVEL'][i_level])
                              + '{:16.8f}'.format(energy_level['ENERGIES'][i_level])
                              + '{:9.1f}'.format(energy_level['WEIGHT'][i_level])
                              + '     ' + str(energy_level['QN'][i_level]) + '\n')
        leiden_file.write('!NUMBER OF RADIATIVE TRANSITIONS\n')
        leiden_file.write(str(len(radiative_transitions['TRANS'])) + '\n')
        leiden_file.write('!TRANS + UP + LOW + EINSTEINA[s^-1] + FREQ[GHz] + E_u[K]\n')
        for i_trans in range(len(radiative_transitions['TRANS'])):
            leiden_file.write('{:5d}'.format(radiative_transitions['TRANS'][i_trans])
                              + '{:5d}'.format(radiative_transitions['UP'][i_trans])
                              + '{:5d}'.format(radiative_transitions['LOW'][i_trans])
                              + '{:12.3e}'.format(radiative_transitions['EINSTEINA'][i_trans])
                              + '{:17.7f}'.format(radiative_transitions['FREQ'][i_trans])
                              + '{:10.2f}'.format(radiative_transitions['E'][i_trans]) + '\n')
        leiden_file.close()
