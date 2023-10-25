<common>

	<dust_component>	"input/dust/silicate_d03.nk" "plaw" 1.0 3800.0 5e-09 250e-9 -3.5
    <phase_function>	PH_MIE

	<mass_fraction>		0.01

	<nr_threads>		-1

</common>

<task> 1

    <cmd>			CMD_DUST_SCATTERING

	<detector_dust_mc nr_pixel = "255*255">	1e-6	1e-3	4	0.00	0.00	4.32e+18
	<detector_dust_mc nr_pixel = "255*255">	1e-6	1e-3	4	26.0	30.0	4.32e+18
	<detector_dust_mc nr_pixel = "255*255">	1e-6	1e-3	4	74.0	67.0	4.32e+18
	<detector_dust_mc nr_pixel = "255*255">	1e-6	1e-3	4	111.0	99.0	4.32e+18
	<detector_dust_mc nr_pixel = "255*255">	1e-6	1e-3	4	130.0	170.0	4.32e+18

	<source_star nr_photons = "1e6">	0	0	0	2	4500

	<path_grid>		"projects/test/stellar_scattering_sphere/grid_3Dsphere_m1e-3.dat"
	<path_out>		"projects/test/stellar_scattering_sphere/dust_mc/"

	<peel_off>		1
	<enfsca>		1

</task>
