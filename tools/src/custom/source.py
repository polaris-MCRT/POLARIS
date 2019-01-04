#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from modules.base import StellarSource

"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_sources_dict(dictionary):
    sources_dict = {
        'f_type': FType,
        'gg_tau_stars': GGTauStars,
        'hd97048': HD97048,
        'hd169142': HD169142,
        'custom': CustomStar,
    }
    dictionary.update(sources_dict)


class CustomStar(StellarSource):
    """Change this to the star you want to use.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the stellar source parameters.

        Args:
            file_io : Handles file input/output and all necessary paths.
        """
        StellarSource.__init__(self, file_io, parse_args)

        # Position of the star [m, m, m]
        self.parameter['position'] = [0, 0, 0]
        # Effective temperature of the star [K]
        self.parameter['temperature'] = 4000
        # Radius of the star [R_sun] or luminosity [L_sun]
        self.parameter['radius'] = 2.0 * self.math.const['R_sun']
        # Number of photons if no number is chosen via --photons
        self.parameter['nr_photons'] = 1e6
        # Can the velocity field be calculated by only this star in the center?
        self.parameter['kepler_usable'] = True
        # Mass of the star [M_sun] (for Keplerian rotation)
        self.parameter['mass'] = 0.7 * self.math.const['M_sun']

    def get_command(self):
        """Provides stellar source command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the stellar source.
        """
        '''To add multiple stars, use the following:
        new_command_line = str()
        self.parameter['temperature'] = 8000
        self.parameter['radius'] = 4.0 * self.math.const['R_sun']
        new_command_line += self.get_command_line()
        self.parameter['temperature'] = 5000
        self.parameter['radius'] = 3.0 * self.math.const['R_sun']
        new_command_line += self.get_command_line()
        return new_command_line
        '''
        return self.get_command_line()


class FType(StellarSource):
    """Change this to the star you want to use.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the stellar source parameters.

        Args:
            file_io : Handles file input/output and all necessary paths.
        """
        StellarSource.__init__(self, file_io, parse_args)

        # Position of the star [m, m, m]
        self.parameter['position'] = [0, 0, 0]
        # Effective temperature of the star [K]
        self.parameter['temperature'] = 6500
        # Radius of the star [R_sun]
        self.parameter['radius'] = 1.3 * self.math.const['R_sun']
        # Mass of the star [M_sun] (for Keplerian rotation)
        self.parameter['mass'] = 0.7 * self.math.const['M_sun']
        # Number of photons if no number is chosen via --photons
        self.parameter['nr_photons'] = 1e6


class GGTauStars(StellarSource):
    """The BinaryStar class is three pre-main sequence stars in the center
    and a planet in the disk.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the stellar source parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.

        Notes:
            - White et. al 1999 (separation of components, luminosity)
                "http://iopscience.iop.org/article/10.1086/307494/pdf"
            - Correia et. al 2008 (stellar temperature/spectral type)
                "https://arxiv.org/pdf/astro-ph/0608674.pdf"

        """
        StellarSource.__init__(self, file_io, parse_args)
        # ------ Default nr. of photons -----
        self.parameter['nr_photons'] = 1e8

        # ------ Effective temperatures -----
        # Cite: temperatures and spectral types (Di Folco et al. 2014)
        # GG Tau Aa: M0 spectral type
        self.T_Aa = 3700.
        # GG Tau Aa: M2 spectral type
        self.T_Ab1 = 3300.
        # GG Tau Aa: M3 spectral type
        self.T_Ab2 = 3100.
        # Cite:
        self.T_planet = 839.9

        # ------ Luminosities -----
        # Cite: luminosity of Aa (White et al. 1999)
        self.L_Aa = 0.84 * self.math.const['L_sun']
        # Cite: luminosity of Aa (Hartigan et al. 2003)
        #self.L_Aa = (0.38 + 0.122) * self.math.const['L_sun']
        # Cite: luminosity of Aa (White et al. 1999 and Di Folco et al. 2014)
        self.L_Ab1 = 0.89 * 0.71 * self.math.const['L_sun']
        self.L_Ab2 = (1 - 0.89) * 0.71 * self.math.const['L_sun']
        # Cite: luminosity of Aa (Hartigan et al. 2003 and Di Folco et al. 2014)
        #self.L_Ab1 = 0.89 * (0.2 + 0.079) * self.math.const['L_sun']
        #self.L_Ab2 = (1 - 0.89) * (0.2 + 0.079) * self.math.const['L_sun']
        # Cite:
        self.L_planet = 1e-1 * \
            (1.4e-5 + 1.863234318727217e-3) * self.math.const['L_sun']

        # ------ Half-major axis of the stars -----
        # Cite: separations (White et al. 1999)
        rot_angle_2 = 25. + 15.
        self.a_Aab = 36. / 2. * np.sqrt(
            (np.cos(rot_angle_2 / 180 * np.pi) / np.cos(37 / 180 * np.pi))**2 + 
            np.sin(rot_angle_2 / 180 * np.pi)**2) 
        self.a_Ab12 = 5. / 2.
        self.a_planet = 260. + 20.

        # Cite: position angle (Di Folco et al. 1999)
        self.angle_Aa = 3. / 2. * np.pi
        self.angle_Ab = self.angle_Aa + np.pi
        # Cite: position of planet (Dutrey et al. 2014)
        self.angle_planet = np.pi * (360. - 127.) / 180.

        # Add planet to sources?
        self.add_planet = False

        '''
        From Robert Brunngräber!
        1 million years
        1 M_jup:  T = 839.9 K, L = 1.409e-5 L_sun
        10 M_jup: T = 2423 K,  L = 1.950e-3 L_sun

        10 million years
        1 M_jup:  T = 528.7 K, L = 1.423e-6 L_sun
        10 M_jup: T = 1591 K,  L = 1.262e-4 L_sun

        Accreting circumplanetary disks: observational signatures (Zhu 2015)
        # print(self.math.planet_lum(5.683e26, 1e-8, 58232e3))
        # Saturn like star -> L = 5.37e-4 L_sun
        #                  -> T = 4000 K

        '''
        #: dict: Parameters for the binary components
        self.tmp_parameter = {
            # New: M0, M2, M3 (http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt)
            'temperature': [self.T_Aa, self.T_Ab1, self.T_Ab2],
            'luminosity': [self.L_Aa, self.L_Ab1, self.L_Ab2],
            'position_star': [[-1.0, self.a_Aab * self.math.const['au'] * np.sin(self.angle_Aa), 0.],
                              [self.a_Aab * self.math.const['au'] * np.cos(self.angle_Ab) -
                               self.a_Ab12 * self.math.const['au'],
                               self.a_Aab * self.math.const['au'] * np.sin(self.angle_Ab), 0.],
                              [self.a_Aab * self.math.const['au'] * np.cos(self.angle_Ab) +
                               self.a_Ab12 * self.math.const['au'],
                               self.a_Aab * self.math.const['au'] * np.sin(self.angle_Ab), 0.]]
        }

    def get_command(self):
        new_command_line = str()
        for i_comp in range(len(self.tmp_parameter['temperature'])):
            self.parameter['temperature'] = self.tmp_parameter['temperature'][i_comp]
            self.parameter['luminosity'] = self.tmp_parameter['luminosity'][i_comp]
            self.parameter['position'] = self.tmp_parameter['position_star'][i_comp]
            new_command_line += self.get_command_line()
        if self.add_planet:
            self.parameter['temperature'] = self.T_planet
            self.parameter['luminosity'] = self.L_planet
            self.parameter['position'] = [self.a_planet * self.math.const['au'] * np.sin(self.angle_planet),
                                          self.a_planet * self.math.const['au'] * np.cos(self.angle_planet), 0.]
            new_command_line += self.get_command_line()
        return new_command_line


class HD97048(StellarSource):
    """Change this to the star you want to use.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the stellar source parameters.

        Args:
            file_io : Handles file input/output and all necessary paths.
        """
        StellarSource.__init__(self, file_io, parse_args)

        # Position of the star [m, m, m]
        self.parameter['position'] = [0, 0, 0]
        # Effective temperature of the star [K]
        self.parameter['temperature'] = 1e4
        # Radius of the star [R_sun]
        self.parameter['radius'] = 2.0 * self.math.const['R_sun']
        # Mass of the star [M_sun] (for Keplerian rotation)
        self.parameter['mass'] = 1.0 * self.math.const['M_sun']
        # Number of photons if no number is chosen via --photons
        self.parameter['nr_photons'] = 1e6


class HD169142(StellarSource):
    """Change this to the star you want to use.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the stellar source parameters.

        Args:
            file_io : Handles file input/output and all necessary paths.
        """
        StellarSource.__init__(self, file_io, parse_args)

        # Position of the star [m, m, m]
        self.parameter['position'] = [0, 0, 0]
        # Effective temperature of the star [K]
        self.parameter['temperature'] = 7800.
        # Radius of the star [R_sun]
        self.parameter['luminosity'] = 9.8 * self.math.const['L_sun']
        # Number of photons if no number is chosen via --photons
        self.parameter['nr_photons'] = 1e6