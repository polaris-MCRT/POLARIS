# The command file consists of an arbitrary number of tasks.
# Commands that are common to all tasks can be defined in an extra block in the beginning as the first block of the command file.
# Lines marked with a # or a ! are ignored and can be used as comments.

# common commands to all tasks
<common>

    # Dust

    # considered dust grain compositions (ALL) "path_to_file" "size_keyword" mass_fraction mass_density radius_min radius_max exponent
    <dust_component>	"/PATH/TO/POLARIS/input/dust/silicate_d03.nk" "plaw" 1.0 3500.0 5e-09 2.5e-07 -3.5

    # considered alignment mechanism of non-spherical dust grains (DUST_EMISSION)
    <align> ALIG_PA or ALIG_RAT or ALIG_IDG or ALIG_GOLD or ALIG_INTERNAL or ALIG_NONPA

    # enables/disables sublimation of dust grains, sublimation temperature is set in the dust catalog (TEMP, TEMP_RAT)
    <sub_dust> 1 (yes) or 0 (no) 

    # ratio of grains at high-J attractor point for the radiative torque alignment (DUST_EMISSION)
    <f_highJ> 0.25

    # correlation factor for the radiative torque alignment (DUST_EMISSION)
    <f_c> 0.6

    # reference value for the radiative torque efficiency (RAT, TEMP_RAT)
    <Q_ref> 0.4

    # exponent for the radiative torque efficiency power law (RAT, TEMP_RAT)
    <alpha_Q> 3.0

    # Rayleigh reduction factor for not-perfect alignment or RAT alignment if imperfect internal alignment is not used (TEMP, RAT, TEMP_RAT, DUST_EMISSION)
    <R_rayleigh> 1.0

    # overwrite the gas temperature with the dust temperature times TEMP_FACTOR (TEMP, TEMP_RAT)
    <adj_tgas> 1.0

    # gas-to-dust mass ratio (ALL)
    <mass_fraction> 0.01

    # enables/disables the use of the dust temperature of the input grid as an offset (TEMP, TEMP_RAT)
    <dust_offset> 1 (yes) or 0 (no)

    # enables/disables the saving of the radiation field, required by <stochastic_heating>
    <radiation_field> 1 (yes) or 0 (no)

    # disables the usage of the radiation field to include scattering in the ray-tracing (DUST_EMISSION)
    <rt_scattering> 0

    # enables/disables the calculation of the dust temperature for each dust grain size (TEMP, TEMP_RAT)
    <full_dust_temp> 1 (yes) or 0 (no)

    # consider the emission of stochastically heated dust grains with a size < SIZE_LIMIT (DUST_EMISSION)
    <stochastic_heating> 1e-8

    # Set the prefactor for alignment limitation of the magnetic field
    <larm_f> 4.1e-19

    # phase function used for the dust grains, Henyey-Greenstein, Mie or isotropic (TEMP, RAT, DUST_SCATTERING)
    <phase_function> PH_HG or PH_MIE or PH_ISO


    # Gas

    # gas species for line emission (LINE_EMISSION) "path" POP_index (1: LTE, 2: FEP, 3: LVG) abundance
    <gas_species> "/PATH/TO/POLARIS/input/gas/co.dat" 2 1e-7

    # average atomic mass unit per gas particle (ALL)
    <mu> 2.0

    # enables/disables the creation of velocity channel maps (LINE_EMISSION)
    <vel_maps> 1 (yes) or 0 (no)

    # enables/disables the interpretation of the velocity of the grid as multiple of the speed of sound
    <vel_is_speed_of_sound> 1 (yes) or 0 (no)

    # mass of the central star and enable the use of Keplerian rotation as the velocity field (LINE_EMISSION)
    <kepler_star_mass> 1.0

    # turbulent velocity of the gas phase to add extra Doppler broadening (LINE_EMISSION)
    <turbulent_velocity> 0


    # Visualization

    # enables the creation of midplane '.fits' files, input and output grid (ALL)
    <write_inp_midplanes> 256
    <write_out_midplanes> 256

    # enables the creation of 3D midplane files, requires <write_xyz_midplanes> command (ALL)
    <write_3d_midplanes> PLANE or PLANE NR_SLICES or PLANE NR_SLICES Z_MIN, Z_MAX

    # enables the inclusion of the radiation field into the midplane '.fits' files (ALL)
    <write_radiation_field> 1 (urad) or 2 (J) or 3 (J field)

    # enables/Disables the inclusion of the Mathis field G0 value into the midplane '.fits' files (ALL) see Mathis et al., 1983; Camps et al., 2015 
    <write_g_zero> 1 (yes) or 0 (no)

    # zoom onto the midplane cut images (ALL)
    <midplane_zoom> 1

    # enables/disables the creation of Gnuplot files to visualize the grid, how many scalar points (ALL)
    <nr_gnu_points> 4000

    # enables/disables the creation of Gnuplot files to visualize the grid, how many vectors (ALL)
    <nr_gnu_vectors> 4000

    # enables/disables the creation of Gnuplot files to visualize the grid, how many lines points (ALL)
    <max_lines> 300

    # enables/disables the creation of AMIRA files, input and output grid (ALL)
    <amira_inp_points> 100
    <amira_out_points> 500


    # Numerical options

    # Set the number of processor cores on which POLARIS can run (ALL)
    <nr_threads> 1

</common>

# task block
# The number 1 is optional. Task blocks with a 0 will be skipped completely.
<task> 1

    # Detectors

    # detector (polar, sliced) for dust emission (DUST_EMISSION) nr_pixel = "Npixel"> λmin λmax Nλ IDsource α1 [°] α2 [°] D [m] dx [m] dy [m] ∆x [m] ∆y [m]  
    <detector_dust nr_pixel = "256*256">	0.00085	0.00085	1	1	0.0	0.0	4.319948614054068e+18
    <detector_dust_polar nr_pixel = "256"> 1e-6 2e-6 2 1 35.1 0 4.31998E+18
    <detector_dust_slice nr_pixel = "256"> 1e-6 2e-6 2 1 35.1 0 4.31998E+18

    # healpix detector for dust emission (DUST_EMISSION) nr_sides = "Nsides"> λmin λmax Nλ IDsource rx [m] ry [m] rz [m] lmin [°] lmax [°] bmin [°] bmax [°] dx [m] dy [m] ∆x [m] ∆y [m]
    <detector_dust_healpix nr_sides = "32"> 1e-6 1e-5 2 1 1.65e+20 0 0 -80 80 -45 45

    # detector for dust scattering (DUST_SCATTERING) nr_pixel = "Npixel"> λmin λmax Nλ α1 [°] α2 [°] D [m] dx [m] dy [m]
    <detector_dust_mc nr_pixel = "256"> 1e-6 2e-6 2 15.1 0 4.31998E+18

    # detector (polar, sliced) for spectral line emission (DUST_EMISSION) nr_pixel = "Npixel" vel_channels = "Nvel"> IDgas species IDtransition IDsource vmin,max [m/s] α1 [°] α2 [°] D [m]
    <detector_line nr_pixel = "256" vel_channels = "35"> 1 2 1 1e5 15.1 0 4.31998E+18
    <detector_line_polar nr_pixel = "256" vel_channels = "35"> 1 5 1 1e5 35.1 0 4.31998E+18
    <detector_line_slice nr_pixel = "256" vel_channels = "101"> 1 5 1 1e5 35.1 0 4.31998E+18

    # healpix detector for spectral line emission (DUST_EMISSION) nr_sides = "Nsides" vel_channels = "Nvel"> IDgas species IDtransition IDsource vmin,max [m/s] rx [m] ry [m] rz [m] lmin [°] lmax [°] bmin [°] bmax [°] vx [m/s] vy [m/s] vz [m/s]
    <detector_line_healpix nr_sides = "32" vel_channels = "35"> 1 2 1 1e5 1.96e+20 0 0 -180 180 -3.5 3.5 1.5e3 -1e3 0

    # detector for synchrotron emission (SYNCHROTRON) nr_pixel = "Npixel"> λmin λmax Nλ IDsource α1 [°] α2 [°] D [m] dx [m] dy [m] ∆x [m] ∆y [m]
    <detector_sync nr_pixel = "256"> 0.71 0.81 2 1 15.1 0 4.31998E+18

    # healpix detector for synchrotron emission (SYNCHROTRON) nr_sides = "Nsides"> λmin λmax Nλ IDsource rx [m] ry [m] rz [m] lmin [°] lmax [°] bmin [°] bmax [°]
    <detector_sync_healpix nr_sides = "512"> 0.03 0.7 11 1.65e+20 0 0 -80 80 -45 45

    # maximum level of subpixeling (DUST_EMISSION, LINE_EMISSION)
    <max_subpixel_lvl> 1

    # rotation axis for first and second rotation angle (ALL) axis_x axis_y axis_z
    <axis1>	1	0	0
    <axis2>	0	1	0


    # Radiation sources

    # stellar radiation source (TEMP, RAT, DUST_SCATTERING, DUST_EMISSION) nr_photons = "Nph"> x [m] y [m] z [m] R [Rsun] T [K] q u
    # or nr_photons = "Nph"> x [m] y [m] z [m] "path"
    <source_star nr_photons = "5e6"> 1.15e2 4.16e2 1.8e1 12 17000

    # starfield radiation source (TEMP, RAT, DUST_SCATTERING) nr_photons = "Nph"> x [m] y [m] z [m] R [Rsun] T [K] μf [m]
    # or nr_photons = "Nph"> x [m] y [m] z [m] μf [m] "path"
    <source_starfield nr_photons = "5e6"> 1.15e2 4.16e2 1.8e1 12 23e4 1e18

    # interstellar radiation field as radiation source (TEMP, RAT) nr_photons = "Nph"> "path"
    # or nr_photons = "Nph"> G0
    # or nr_photons = "Nph"> G0 fradius
    <source_isrf nr_photons = "1e7"> "/PATH/TO/POLARIS/input/interstellar_radiation_field.dat"

    # background radiation source (DUST_EMISSION, LINE_EMISSION, SYNCHROTRON) nr_photons = "Nph"> f [a.u.] Teff [K] q u v α1 [°] α2 [°]
    # or nr_photons = "Nph"> "path" α1 [°] α2 [°]
    <source_background nr_photons = "1e5"> 1e18 30000 1 0 0 90 45

    # use dust grains as radiation source, activates self-scattering calculations (DUST_SCATTERING) nr_photons = "Nph">
    <source_dust nr_photons = "1e6">

    # laser radiation source (DUST_EMISSION, DUST_SCATTERING) nr_photons = "Nph"> x [m] y [m] z [m] dx [m] dy [m] dz [m] P [W] λ0 [m] FWHM [m] q u
    <source_laser nr_photons = "1e6"> 0 0 0 0 0 0 1e-3 1e-6 1e-9 0 0


    # Simulation types

    # calculate temperature distribution 
    <cmd>          CMD_TEMP

    # calculate temperature distribution and RAT alignment radius
    <cmd>          CMD_TEMP_RAT

    # calculate RAT alignment radius
    <cmd>          CMD_RAT

    # calculate Monte Carlo dust scattering
    <cmd>          CMD_DUST_SCATTERING

    # ray-tracing for thermal emission of dust grains
    <cmd>          CMD_DUST_EMISSION

    # ray-tracing of the spectral line emission of the gas
    <cmd>          CMD_LINE_EMISSION

    # ray-tracing of synchrotron emission and Faraday rotation
    <cmd>          CMD_SYNCHROTRON


    # General

    # Path to the used grid (ALL)
    <path_grid>		"/PATH/TO/GRID"

    # Path to the POLARIS results (ALL)
    <path_out>		"/PATH/TO/RESULTS"

    # first considered detector index as occurrence in command file (DUST_EMISSION, LINE_EMISSION, SYNCHROTRON)
    <start> 1

    # last considered detector index as occurrence in command file (DUST_EMISSION, LINE_EMISSION, SYNCHROTRON)
    <stop> 1

    # conversion factor for the density (ALL)
    <conv_dens> 1

    # conversion factor for the length (ALL)
    <conv_len> 1

    # conversion factor for the magnetic field strength (ALL)
    <conv_mag> 1

    # conversion factor for the velocity (ALL)
    <conv_vel> 1

    # enables/disables the use of enforced first scattering (TEMP, RAT, DUST_SCATTERING)
    <enfsca> 1 (yes) or 0 (no)

    # enables/disables the use of the peel-off technique (DUST_SCATTERING)
    <peel_off> 1 (yes) or 0 (no)

    # acceptance angle, if peel-off is not used (DUST_SCATTERING)
    <acceptance_angle> 1

</task>

