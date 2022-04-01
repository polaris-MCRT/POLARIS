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
	# "path_to_file" index_lvl_pop abundance "path_to_zeeman_file"
	<gas_species>	"/YOUR/POLARIS/PATH/input/gas/oh@hfs.dat"	1	1e-07	"/YOUR/POLARIS/PATH/input/gas/oh@hfs_zeeman.dat"

	# background radiation source
	# background = "N_photons"> scaling_factor temperature q u v rot_1 rot_2
	<source_background nr_photons = "1e5">	0.0	0.0	0.0	0.0	0.0	0.0	0.0

	# detector for line emission = "N_pixel x N_pixel" vel_channeles = "N_velocity_bins"> ID_gas ID_transition ID_source nu_min_max rot_1 rot_2 distance
	<detector_line nr_pixel = "256*256" vel_channels = "35">	1	2	1	1000	0.0	0.0	4.319948614054068e+18

	# ray-tracing of the spectral line emission
	<cmd>		CMD_LINE_EMISSION

	# path of the input grid file
	<path_grid>		"/YOUR/POLARIS/PATH/projects/disk/example/temp/grid.dat"

	# path for all output data
	<path_out>		"/YOUR/POLARIS/PATH/projects/disk/example/zeeman/"

	# Enables the creation of velocity channel maps
	<vel_maps>		1

	# Set the turbulent velocity of the gas phase to add extra Doppler broadening
	<turbulent_velocity>	100.0

</task>

