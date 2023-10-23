<common>
	<dust_component>	"input/dust/silicate_d03.nk" "plaw" 1.0 3800.0 5e-09 250e-9 -3.5
    <phase_function>	PH_MIE

	<mass_fraction>		0.01

	<nr_threads>		-1
</common>

<task> 1
	<cmd>			CMD_DUST_SCATTERING

	<detector_dust_mc nr_pixel = "255*255">	50e-9	1e-3	20	0.0	0.0	4.32e+18

	<source_star nr_photons = "1">	0	0	0	2	4500

	<path_grid>		"projects/test/stellar_sed/grid_empty_1Dsphere.dat"
	<path_out>		"projects/test/stellar_sed/dust_mc/"

	<peel_off>		1
	<enfsca>		0
</task>
