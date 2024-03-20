# common commands to all tasks
<common>

    # set the considered gaseous compositions
    # "path_to_file"
    <gas_component id = "0">    "input/cross_sections/molecular_hydrogen.dat"
    # phase function for scattering
    <phase_function id = "0">    PH_RAYLEIGH

    # set the considered dust/cloud grain compositions
    # "path_to_file" "size_keyword" mass_fraction mass_density radius_min radius_max parameter
    <dust_component id = "1">    "input/refractive_indices/water_s81.nk" "plaw-ed" 1.0 1000.0 1.0e-09 1.0e-5 7.0 0.0 2.0e-7 1.0

    #    size_keyword = "plaw": 1 parameter
    #        n(r) = constant * r^(p0)
    #    size_keyword = "plaw-ed": 4 parameter
    #        n(r) = constant * r^(p0) * exp( -((r - p1) / p2)^(p3) )
    #        p0 = (1 - 3 * v_eff) / v_eff
    #        p3 = r_eff * v_eff
    
    # phase function for dust scattering
    <phase_function id = "1">    PH_MIE

    # number of CPU cores that will be used by POLARIS
    <nr_threads>    -1

</common>

<task> 1

    # rotation axis for first and second rotation angle
    <axis1>    1    0    0
    <axis2>    0    0    1

    # detector for dust scattering = "N_pixel"> wavelength_min wavelength_max wavelength_nr rot_1 rot_2 distance (sidelength_x sidelength_y)
    # use one detector for each phase angle
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 0.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 5.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 10.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 15.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 20.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 25.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 30.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 35.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 40.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 45.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 50.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 55.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 60.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 65.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 70.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 75.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 80.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 85.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 90.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 95.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 100.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 105.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 110.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 115.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 120.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 125.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 130.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 135.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 140.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 145.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 150.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 155.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 160.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 165.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 170.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 175.0 3.0856775814671917e+17

    # spatially extended star as radiation source
    # star = "N_photons"> pos_x pos_y pos_z radius temperature biased_emission
    <source_extended_star nr_photons = "1e6">    0 -14959787070 0 1.0 6000 1

    # Monte Carlo dust scattering in planetary atmosphere
    <cmd>    CMD_PLANET_SCATTERING

    # path of the input grid file
    <path_grid>    "projects/cloudy/cloudy.dat"
    # path for all output data
    <path_out>    "projects/cloudy/"

    # surface reflection
    # "type" parameter
    <surface_reflection>    "lambertian" 0.0

    # enable peel off
    <peel_off>    1
    # enable enforced first scattering
    <enfsca>    1

</task>
