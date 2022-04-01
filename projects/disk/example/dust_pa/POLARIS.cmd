# common commands to all tasks
<common>

	# set the considered dust grain compositions
	# "path_to_file" "size_keyword" mass_fraction mass_density radius_min radius_max exponent
	<dust_component>	"/YOUR/POLARIS/PATH/input/dust/silicate_oblate.dat" "plaw" 0.625 3800 5e-09 2e-06 -3.5
	<dust_component>	"/YOUR/POLARIS/PATH/input/dust/graphite_oblate.dat" "plaw" 0.375 2250 5e-09 2e-06 -3.5

	# phase function for dust scattering
	<phase_function>	PH_ISO

	# Enables the creation of midplane '.fits' files and define the number of pixel (input and output)
	<write_inp_midplanes>	256
	<write_out_midplanes>	256

	# dust to gas mass ratio
	<mass_fraction>		0.01

	# number of CPU cores that will be used by POLARIS
	<nr_threads>		-1

</common>

# first task
<task> 1
	
	# detector for dust emission = "N_pixel x N_pixel"> wavelength_min wavelength_max wavelength_nr ID_source rot_1 rot_2 distance
	<detector_dust nr_pixel = "256*256">	0.00085	0.00085	1	1	90.0	0.0	4.319948614054068e+18

	# background radiation source
	# background = "N_photons"> scaling_factor temperature q u v rot_1 rot_2
	<source_background nr_photons = "1e5">	0.0	0.0	0.0	0.0	0.0	0.0	0.0

	# rotation axis for first and second rotation angle
	<axis1>	1	0	0
	<axis2>	0	1	0

	# ray-tracing for thermal emission of dust grains
	<cmd>			CMD_DUST_EMISSION

	# path of the input grid file
	<path_grid>		"/YOUR/POLARIS/PATH/projects/disk/example/temp/grid.dat"

	# path for all output data
	<path_out>		"/YOUR/POLARIS/PATH/projects/disk/example/dust_pa/"
	
	# Considered alignment mechanism of non-spherical dust grains
	<align>			ALIG_PA

	# Set the f_{high,j} value for the radiative torque mechanism
	<f_highJ>		0.25

	# Set the f_{c} value for the internal alignment correlation
	<f_c>			0.0

</task>

