# common commands to all tasks
<common>

	# set the considered dust grain compositions
	# "path_to_file" "size_keyword" mass_fraction mass_density radius_min radius_max exponent
	<dust_component>	"input/dust/silicate_d03.nk" "plaw" 0.625 3500.0 5e-09 2.5e-07 -3.5
	<dust_component>	"input/dust/graphite_perpend_d03.nk" "plaw" 0.25 2250 5e-09 2.5e-07 -3.5
	<dust_component>	"input/dust/graphite_parallel_d03.nk" "plaw" 0.125 2250 5e-09 2.5e-07 -3.5

	# phase function for dust scattering
	<phase_function>	PH_MIE

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

	# star in the center as radiation source
	# star = "N_photons"> pos_x pos_y pos_z radius temperature
	<source_star nr_photons = "1e6">	0	0	0	0.9	4000

	# dust heating
	<cmd>			CMD_TEMP

	# path of the input grid file
	<path_grid>		"projects/disk/grid.dat"

	# path for all output data
	<path_out>		"projects/disk/example/temp/"

</task>

