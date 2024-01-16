<common>
    # considered dust grain composition "path_to_file" "size_distribution_keyword" mass_fraction mass_density a_min a_max a_0 sigma
	<dust_component>	"/zfshome/supas357/progs/POLARIS_DustyPlasma/polaris/input/dust/2_0p02_660.nk" "logn" 1.0 1.0 2e-07 4e-07 3e-07 0.1

	# phase function used for the dust grains, Henyey-Greenstein, Mie or isotropic (TEMP, RAT, DUST_SCATTERING)
	<phase_function>	PH_MIE

    # enables the creation of midplane '.fits' files displaying the input density distribution  with the given pixel resolution
	<write_inp_midplanes>	256
	<midplane_zoom>		1

    # Set the number of processor cores on which POLARIS can run
	<nr_threads>		-1
</common>

<task> 1
    # detectors for dust scattering (DUST_SCATTERING) nr_pixel = "Npixel"> λmin λmax Nλ α1 [°] α2 [°] D [m]
	<detector_dust_mc nr_pixel = "51*51">	6.6e-07	6.6e-07	1	0.0	90.0	0.5
	<detector_dust_mc nr_pixel = "51*51">	6.6e-07	6.6e-07	1	270.0	0.0	0.5

	# rotation axis for first and second rotation angle (ALL) axis_x axis_y axis_z, defualt detector position at +z
	<axis1>	1	0	0
	<axis2>	0	1	0

    # laser radiation source nr_photons = "Nph"> x [m] y [m] z [m] dx [m] dy [m] dz [m] P [W] λ0 [m] FWHM [m] q u
	<source_laser nr_photons = "1e7">	-1   0	0	1	0	0	1.0	6.6e-07	 10e-09	0.0	1

    # calculate Monte Carlo dust scattering
	<cmd>			CMD_DUST_SCATTERING

    # Path to the used grid
	<path_grid>		"/zfshome/supas357/progs/POLARIS_DustyPlasma/polaris/projects/constantCylinder/example2.grid"

	# Path to the POLARIS results
	<path_out>		"/zfshome/supas357/progs/POLARIS_DustyPlasma/polaris/projects/constantCylinder/example2/"

    # enables/disables the use of the peel-off technique
	<peel_off>		1

	# enables/disables the use of enforced first scattering
	<enfsca>		1

</task>