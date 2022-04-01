# common commands to all tasks
<common>

	# Enables the creation of midplane '.fits' files and define the number of pixel (input and output)
	<write_inp_midplanes>	256
	<write_out_midplanes>	256

	Average atomic mass unit per gas particle
	<mu>			2.0

	# number of CPU cores that will be used by POLARIS
	<nr_threads>		-1

</common>

# first task
<task> 1

	# define the gas species that is used by POLARIS
	# "path_to_file" index_lvl_pop abundance
	<gas_species>	"/YOUR/POLARIS/PATH/input/gas/c18o.dat"	1	1e-05

	# background radiation source
	# background = "N_photons"> scaling_factor temperature q u v rot_1 rot_2
	<source_background nr_photons = "1e5">	0.0	0.0	0.0	0.0	0.0	0.0	0.0

	# detector for line emission = "N_pixel x N_pixel" vel_channeles = "N_velocity_bins"> ID_gas ID_transition ID_source nu_min_max rot_1 rot_2 distance
	<detector_line nr_pixel = "256*256" vel_channels = "35">	1	1	1	3000	90.0	0.0	4.319948614054068e+18

	# rotation axis for first and second rotation angle
	<axis1>	1	0	0
	<axis2>	0	1	0

	# ray-tracing of the spectral line emission
	<cmd>		CMD_LINE_EMISSION

	# path of the input grid file
	<path_grid>		"/YOUR/POLARIS/PATH/projects/disk/example/temp/grid.dat"

	# path for all output data
	<path_out>		"/YOUR/POLARIS/PATH/projects/disk/example/line/"

	# Enables the creation of velocity channel maps
	<vel_maps>		1

	# Set the mass of the central star and enable the use of Keplerian rotation as the velocity field
	<kepler_star_mass>	0.7

	# Set the turbulent velocity of the gas phase to add extra Doppler broadening
	<turbulent_velocity>	100.0

</task>

