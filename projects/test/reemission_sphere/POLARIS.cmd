<common>

	<dust_component>	"input/refractive_indices/silicate_d03.nk" "plaw" 1.0 3800.0 1e-06 1e-06 -3.5
    <phase_function>	PH_MIE

	<mass_fraction>		0.01

	<nr_threads>		-1

</common>

<task> 1

	<cmd>			CMD_DUST_EMISSION

	<detector_dust_polar nr_pixel = "255*255">	1e-6	1e-2	10	1	0.0	0.0	4.32e+18
	<detector_dust nr_pixel = "255*255">		1e-6	1e-2	10	1	0.0	0.0	4.32e+18

	<max_subpixel_lvl>	2

	<path_grid>		"projects/test/reemission_sphere/grid_3D_sphere_const_T_m1e-5.dat"
	<path_out>		"projects/test/reemission_sphere/dust/"

</task>
