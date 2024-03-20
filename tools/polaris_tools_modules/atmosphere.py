# =====================================================================================

import numpy as np
from polaris_tools_modules.math import Math

# =====================================================================================
class AtmosphereRoutines:
    def __init__(self):
        # Temperature at standard conditions [K]
        self.con_temp_std = 273.15
        # Pressure at standard conditions [Pa = kg m^(-1) s^(-3)]
        self.con_pres_std = 101325.0
        # Density at standard conditions assuming ideal gas [m^(-3)]
        self.const_dens_std = self.con_pres_std / self.con_temp_std / Math().const["k_B"]

        # Molar mass [kg mol^(-1)]
        # Cox, 2000, Allen's astrophysical quantities (Springer)
        # https://ui.adsabs.harvard.edu/abs/2000asqu.book.....C
        self.const_molar_mass = {
            "1H2": 2.016 * 1.0e-3,
            "1H2-16O": 18.015 * 1.0e-3,
            "2He" : 4.003 * 1.0e-3,
            "12C-1H4": 16.043 * 1.0e-3,
            "12C-16O": 28.010 * 1.0e-3,
            "12C-16O2": 44.01 * 1.0e-3,
            "14N2": 28.013 * 1.0e-3,
            "14N-1H3": 17.031 * 1.0e-3,
            "16O2": 31.999 * 1.0e-3,
            "40Ar": 39.948 * 1.0e-3,
            "31P-1H3": 33.998 * 1.0e-3,
            "air": 28.964 * 1.0e-3,
            "rayleigh": 1.0e-3,
            }
        # "rayleigh" is an artificial rayleigh scattering particle with depolarization 0

        # Depolarization factor
        # Hansen & Travis, 1974, Space Sci. Rev. 16, 527
        # https://ui.adsabs.harvard.edu/abs/1974SSRv...16..527H
        # Penndorf, 1957, J. Opt. Soc. Am. 47, 176
        # https://ui.adsabs.harvard.edu/abs/1957JOSA...47..176P
        # Bates, 1984, Planet. Space Sci. 32, 785
        # https://ui.adsabs.harvard.edu/abs/1984P&SS...32..785B
        # Sneep & Ubachs, 2005, J. Quant. Spectr. Rad. Transf. 92, 293
        # https://ui.adsabs.harvard.edu/abs/2005JQSRT..92..293S
        self.const_depolarization = {
            "1H2": 0.02,
            "2He": 0.0,
            "12C-1H4": 0.0,
            "12C-16O": 0.01,
            "12C-16O2": 0.09,
            "14N2": 0.03,
            "16O2": 0.06,
            "40Ar": 0.0,
            "air": 0.03,
            "rayleigh": 0.0,
            }
        # "rayleigh" is an artificial rayleigh scattering particle with depolarization 0

        # Particle dependent constants A, B for refractive index
        # Cox, 2000, Allen's astrophysical quantities (Springer)
        # https://ui.adsabs.harvard.edu/abs/2000asqu.book.....C
        # Sneep & Ubachs, 2005, J. Quant. Spectr. Rad. Transf. 92, 293
        # https://ui.adsabs.harvard.edu/abs/2005JQSRT..92..293S
        self.refractive_const = {
            "1H2": [13.58 * 1.0e-5, 7.52 * 1.0e-3],
            "2He": [3.48 * 1.0e-5, 2.3 * 1.0e-3],
            "12C-1H4": [46.66 * 1.0e-5, 8.62 * 1.0e-3],
            "12C-16O": [32.70 * 1.0e-5, 8.10 * 1.0e-3],
            "12C-16O2": [43.9 * 1.0e-5, 6.40 * 1.0e-3],
            "14N2": [29.06 * 1.0e-5, 7.70 * 1.0e-3],
            "14N-1H3": [37.0 * 1.0e-5, 12.00 * 1.0e-3],
            "16O2": [26.63 * 1.0e-5, 5.07 * 1.0e-3],
            "40Ar":  [27.92 * 1.0e-5, 5.60 * 1.0e-3],
            "air": [28.71 * 1.0e-5, 5.67 * 1.0e-3],
            "rayleigh": [1.0e-5, 1.0e-3],
            }
        # "rayleigh" is an artificial rayleigh scattering particle with depolarization 0

    # =================================================================================
    def getMolarMass(self, particle):
        """ Get the molar mass for a given particle

            Parameters
            ----------
            particle : str
                Particle defined in self.refractive_const{}, e.g. 1H2

            Returns
            -------
            float
                Molar mass [kg mol^(-1)]

            Example
            -------
            >>> getMolarMass('1H2')
        """
        mol_weight = self.const_molar_mass.get(particle)

        if isinstance(mol_weight, float):
            return mol_weight
        else:
            ValueError(f"Molar mass is not defined for {particle}")

    # =================================================================================
    def getDepolarization(self, particle):
        """ Get the depolarization factor for a given particle

            Parameters
            ----------
            particle : str
                Particle defined in self.refractive_const{}, e.g. 1H2

            Returns
            -------
            float
                Depolarization factor [dimensionless]

            Example
            -------
            >>> getDepolarization('1H2')
        """
        depol = self.const_depolarization.get(particle)

        if isinstance(depol, float):
            return depol
        else:
            ValueError(f"Depolarization factor is not defined for {particle}")

    # =================================================================================
    def getRefractiveIndex(self, wavelength, particle):
        """ Calculate the real refractive index at a given wavelength for a given particle
        Cox, 2000, Allen's astrophysical quantities (Springer)
        https://ui.adsabs.harvard.edu/abs/2000asqu.book.....C

            Parameters
            ----------
            wavelength : float or array-like
                Wavelength [m]
            particle : str
                Particle defined in self.refractive_const{}, e.g. 1H2

            Returns
            -------
            float
                Real part of refractive index [dimensionless]

            Example
            -------
            >>> getRefractiveIndex(5.5e-7, '1H2')
            >>> getRefractiveIndex([5.0e-7, 5.5e-7, 6.0e-7], '1H2')
        """
        A, B = self.refractive_const.get(particle)
        if isinstance(A, float) and isinstance(B, float):
            return A * (1.0 + B / (wavelength * 1.0e6)**2) + 1.0
        else:
            raise ValueError(f"A or B is not defined for {particle}")

    # =================================================================================
    def getRayleighCrossSection(self, wavelength, particle):
        """ Calculate the rayleigh scattering cross section at a given wavelength for a given particle
        Sneep & Ubachs, 2005, J. Quant. Spectr. Rad. Transf. 92, 293
        https://ui.adsabs.harvard.edu/abs/2005JQSRT..92..293S

            Parameters
            ----------
            wavelength : float
                Wavelength [m]
            particle : str
                Particle defined in self.refractive_const{}, e.g. 1H2

            Returns
            -------
            float
                Rayleigh scattering cross section [m^2]

            Example
            -------
            >>> getRayleighCrossSection(5.5e-7)
            >>> getRayleighCrossSection([5.0e-7, 5.5e-7, 6.0e-7], '1H2')
        """
        n_r = self.getRefractiveIndex(wavelength, particle)
        depolarization = self.getDepolarization(particle)

        A = 24.0 * np.pi**3 / wavelength**4 / self.const_dens_std**2
        B = (n_r**2 - 1.0) / (n_r**2 + 2.0)
        C = (6.0 + 3.0 * depolarization) / (6.0 - 7.0 * depolarization)
        return A * B * B * C

    # =================================================================================
    def getRayleighCrossSection2(self, wavelength, refractive_index, depolarization=0.0):
        """ Calculate the rayleigh scattering cross section at a given wavelength for a refractive index formula
        Sneep & Ubachs, 2005, J. Quant. Spectr. Rad. Transf. 92, 293
        https://ui.adsabs.harvard.edu/abs/2005JQSRT..92..293S

            Parameters
            ----------
            wavelength : float
                Wavelength [m]
            refractive_index : str
                User defined formula for the wavelength-dependent refractive index, wavelength has to expressed as 'x'
            depolarization : float, default = 0
                Depolarization factor

            Returns
            -------
            float
                Rayleigh scattering cross section [m^2]

            Example
            -------
            >>> getRayleighCrossSection2(5.5e-7, '13.58e-5 * (1.0 + 7.52e-3 / (x * 1.0e6)**2) + 1.0', 0.02)
        """
        try:
            n_r = lambda x: eval(refractive_index)
        except:
            ValueError(f"Could not evaluate {refractive_index} into a formula. Is the formula a string and the wavelength expressed as 'x'?")

        A = 24.0 * np.pi**3 / wavelength**4 / self.const_dens_std**2
        B = (n_r**2 - 1.0) / (n_r**2 + 2.0)
        C = (6.0 + 3.0 * depolarization) / (6.0 - 7.0 * depolarization)
        return A * B * B * C

    # =================================================================================
    def getRayleighScatMatrix(self, depolarization=0.0, nang=91):
        """ Calculate the Rayleigh scattering matrix for the given particle type
        Hansen & Travis, 1974, Space Sci. Rev. 16, 527
        https://ui.adsabs.harvard.edu/abs/1974SSRv...16..527H

            Paramters
            ---------
            depolarization : float, default = 0
                Depolarization factor
            nang : float, default = 91
                Number of scattering angles between 0 and pi/2 (equally distributed)

            Returns
            -------
            float
                Scattering angles [rad]
            float
                Matrix element F11 [dimensionless]
            float
                Matrix element F12 [dimensionless]
            float
                Matrix element F22 [dimensionless]
            flaot
                Matrix element F33 [dimensionless]
            float
                Matrix element F44 [dimensionless]
        """
        theta = np.linspace(0.0, np.pi, 2*nang-1)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)

        delta_1 = (1.0 - depolarization) / (1.0 + 0.5 * depolarization)
        delta_2 = (1.0 - 2.0 * depolarization) / (1.0 - depolarization)

        f11 = 0.75 * delta_1 * (1.0 + cos_theta**2) + (1.0 - delta_1)
        f12 = -0.75 * delta_1 * sin_theta**2
        f22 = 0.75 * delta_1 * (1.0 + cos_theta**2)
        f33 = 1.5 * delta_1 * cos_theta
        f44 = 1.5 * delta_1 * delta_2 * cos_theta

        return theta, f11, f12, f22, f33, f44

    # =================================================================================
    def getAltitudeFromPressure(self, pressures, temperatures, molar_mass, gravity, surface_pressure=1e5):
        """ Calculate the altitude assuming hydrostatic equilibrium and an ideal gas

            Parameters
            ----------
            pressures : list of float
                atmospheric pressure profile [Pa] (Pressure at the bottom, ..., Pressure at the top)
            temperatures : list of float
                atmospheric temperature profile [K] (Temperature at the bottom, ..., Temperature at the top)
            molar_mass : float
                Molar mass of the gas particles [kg mol^(-1)]
            gravity : float
                Gravitational acceleration of the planet [m/s^2]
            surface_pressure : float, default = 1e5
                Pressure where the planetary radius is located [Pa]

            Returns
            -------
            list of float
                Altitude [m]

            Example
            -------
            >>> getAltitudeFromPressure([1e3, 1e4, 1e5], [280, 290, 300], getMolarMass('1H2'), 25)
        """
        boundaries = np.zeros_like(pressures, dtype=float)

        scale_height = Math().const["R"] * temperatures / (molar_mass * gravity)
        for p, pressure in enumerate(pressures[0:-1]):
            boundaries[p+1] = boundaries[p] - scale_height[p+1] * np.log(pressures[p+1] / pressure)

        # shift planet_radius to 1 bar (default) boundary
        # np.interp assumes increasing x data points -> reverse arrays
        boundaries -= np.interp(surface_pressure, pressures[::-1], boundaries[::-1])

        return boundaries

    # =================================================================================
    def getNumberDensityFromOpticalDepth(self, optical_depth, length, c_ext):
        """ Calculate the number density of each atmospheric layer for a given optical depth

            Parameters
            -----------
            optical_depth : float
                Optical depth of atmospheric layer [dimensionless]
            length : float
                Length of atmospheric layer [m]
            c_ext : float
                Extinction cross section of the particle [m^2]

            Returns
            -------
            list of float
                Number density
        """
        return optical_depth / c_ext / length

    # =================================================================================
    def getNumberDensityFromPressure(self, pressure, temperature):
        """ Calculate the number density at a given pressure assuming ideal gas law

            Parameters
            ----------
            pressure : float
                Atmospheric pressure [Pa]
            temperature : float
                Atmospheric temperature [K]

            Returns
            -------
            float
                Number density [m^-3]
        """
        return pressure / (Math().const["k_B"] * temperature)

    # =================================================================================
    def getNumberDensityFromPressureAltitude(self, diff_pressure, diff_altitude, molar_mass, gravity):
        """ Calculate the number density from pressure and altitude difference assuming ideal gas law

            Parameters
            ----------
            diff_pressure : float
                Atmospheric pressure difference [Pa]
            diff_altitude : float
                Atmospheric altitude difference [m]
            molar_mass : float
                Molar mass of the gas particles [kg mol^(-1)]
            gravity : float
                Gravitational acceleration of the planet [m/s^2]

            Returns
            -------
            float
                Number density [m^-3]
        """
        return diff_pressure / (molar_mass / Math().const["N_A"] * gravity * diff_altitude)

    # =================================================================================
    def getTemperatureProfile(self, pressures, T_int, gravity, T_irr, kappa_th=1e-3, kappa_v=1e-3):
        """ Calculate the temperature of each atmospheric layer for a given pressure profile
        Guillot, 2010, A&A 520, A27
        https://ui.adsabs.harvard.edu/abs/2010A&A...520A..27G

            Parameters
            ----------
            pressures : list of float
                Atmospheric pressure profile [Pa] (Pressure at the bottom, ..., Pressure at the top)
            T_int : float
                Temperature of the internal flux [K]
                (T_eff^4 = T_eq^4 + T_int^4,
                where T_eq = T_star * np.sqrt(0.5 * R_star / distance) is the equilibrium temperature of the planet)
            gravity : float
                Gravitational acceleration of the planet [m/s^2]
            T_irr : float
                Tempertature of the irradiation [K]
                (T_irr = T_star * np.sqrt(R_star / distance))
            kappa_th : float, default = 1e-3
                mean opacity at thermal wavelengths (infrared) [m^2/kg]
            kappa_v : float, default = 1e-3
                mean opacity at the wavelengths that characterize the incoming stellar irradiation (visual) [m^2/kg]

            Returns
            -------
            list of float
                Temperatures [K]

            Example
            -------
            >>> getTemperatureProfile([1e3, 1e4, 1e5], 100, 25, 1300)
        """
        optical_depths = kappa_th * pressures / gravity
        gamma = kappa_v / kappa_th
        f = 0.25
        sqrt3 = np.sqrt(3.0)
        res = 0.75 * T_int**4 * (2.0 / 3.0 + optical_depths)
        res += 0.75 * T_irr**4 * f * (2.0 / 3.0 + 1.0 / (gamma * sqrt3) + (gamma / sqrt3 - 1.0 / (gamma * sqrt3)) * np.exp(-gamma * optical_depths * sqrt3))
        return res**(0.25)
    
    # =================================================================================
    def getTemperatureProfile2(self, pressures, T_eff, gravity, kappa_th=1e-3):
        """ Calculate the temperature of each atmospheric layer for a given pressure profile
            Simplified model to self.getTemperatureProfile()

            Parameters
            ----------
            pressures : list of float
                Atmospheric pressure profile [Pa] (Pressure at the bottom, ..., Pressure at the top)
            T_eff : float
                Effective Temperature [K]
                (T_eff^4 = T_eq^4 + T_int^4,
                where T_eq = T_star * np.sqrt(0.5 * R_star / distance) is the equilibrium temperature of the planet)
            gravity : float
                Gravitational acceleration of the planet [m/s^2]
            kappa_th : float, default = 1e-3
                mean opacity at thermal wavelengths (infrared) [m^2/kg]

            Returns
            -------
            list of float
                Temperatures [K]

            Example
            -------
            >>> getSimpleTemperatureProfile([1e3, 1e4, 1e5], 100, 25, 1300)
        """
        optical_depths = kappa_th * pressures / gravity
        res = 0.75 * T_eff**4 * (2.0 / 3.0 + optical_depths)
        return res**(0.25)

    # =================================================================================
    def writeInfScatFile(self, particle_name, wavelengths, scattering_angles, scattering_matrix):
        """ Write .inf and .sca files for user defined scattering matrices

            Parameters
            ----------
            particle_name : str
                name of particle (will be the directory name)
            wavelengths : list of float
                Wavelengths [m]
            scattering_angles : list of float
                Scattering angles [degree]
            scattering matrix : array((wavelengths, scattering angles, scattering matrix element))
                Scattering matrix elements [dimensionless]
        """
        import struct
        import os

        if (np.shape(scattering_matrix)[0] != len(wavelengths) or
            np.shape(scattering_matrix)[1] != len(scattering_angles) or
            np.shape(scattering_matrix)[2] != 16):
            raise ValueError("Shape of scattering matrix array does not fit to size of wavelength or scattering angle array")

        path = os.path.join("..", "input", "dust")
        os.makedirs(os.path.join(path, particle_name), exist_ok=True)

        # write info file
        with open(os.path.join(path, particle_name, "scat.inf", "w")) as file:
            file.write("#nr. of dust species #wav. #inc. angles #phi angle #theta angle\n")
            file.write(f"1 {len(wavelengths)} 1 1 {len(scattering_angles)}\n")
            file.write("#data length\n")
            file.write("16\n")
            file.write("#position of matrix elements\n")
            file.write("#M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44\n")
            file.write("1 2 3 4 ")
            file.write("5 6 7 8 ")
            file.write("9 10 11 12 ")
            file.write("13 14 15 16\n")

        # write binary files
        for w, wave in enumerate(wavelengths):
            with open(os.path.join(path, particle_name, f"wID{w+1:03d}.sca"), "wb") as file:
                for t, scat in enumerate(scattering_angles):
                    for e in range(16):
                        data = struct.pack("f", scattering_matrix[w][t][e])
                        file.write(data)

    # =================================================================================
    def writeCrossSectionFile(self, filename, wavelengths, molar_mass, c_abs, c_sca, depol, stringID, comment=None):
        """ Write cross section (.dat) file with wavelength dependent cross sections

            Parameters
            ----------
            filename : str
                Path and name of the .dat file
            wavelengths list of float:
                Wavelengths [m]
            molar_mass weight : float
                Molar mass of particle
            c_abs : list of float
                Absorption Cross sections [m^2]
            c_sca : list of float
                Scattering Cross sections [m^2]
            depol : float
                Depolarization factor of particles [dimensionless]
            stringID : str
                String ID of the file
            comment : str, default = None
                additional comments
        """
        with open(filename, "w") as file:
            file.write("#string ID\n")
            file.write(f"{stringID}\n")
            if comment:
                file.write(f"\n{comment}\n")

            file.write("\n#wavelength [m] #molar mass [kg/mol]\n")
            file.write(f"{len(wavelengths)} {molar_mass:.9f}\n")

            file.write("\n#wavelength [m]\n")

            for wave in wavelengths:
                file.write(f"{wave:1.9e} ")

            file.write("\n\n#Cext [m^2] Cabs [m^2] Csca [m^2] depolarization\n")
            c_ext = c_abs + c_sca
            for w, wave in enumerate(wavelengths):
                file.write(f"{c_ext[w]:.9e} {c_abs[w]:.9e} {c_sca[w]:.9e} {depol:.3f}\n")

# =====================================================================================
# End of class "GasRoutines"
# =====================================================================================
