# common commands to all tasks
<common>

    # set the considered gaseous compositions
    # "path_to_file"
    <gas_component id = "0">    "input/cross_sections/molecular_hydrogen.dat"
    # phase function for scattering
    <phase_function id = "0">    PH_RAYLEIGH

    # set the considered dust/cloud grain compositions
    # "path_to_file" "size_keyword" mass_fraction mass_density radius_min radius_max parameter
    <dust_component id = "1">    "input/refractive_indices/silicate_d03.nk" "plaw" 1.0 3800.0 1.0e-07 1.0e-5 -3.0
    
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
    # use one detector for each phase angle, here the phase angle is 30 deg, ring inclination is 25 deg
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 65.0 30.0 3.0856775814671917e+17

    # spatially extended star as radiation source
    # star = "N_photons"> pos_x pos_y pos_z radius temperature biased_emission
    <source_extended_star nr_photons = "1e8">    0 -13558171514 6322279208 1.0 6000 1

    # Monte Carlo dust scattering in planetary atmosphere
    <cmd>    CMD_PLANET_SCATTERING

    # path of the input grid file
    <path_grid>    "projects/ringed/ringed.dat"
    # path for all output data
    <path_out>    "projects/ringed/"

    # surface reflection
    # "type" parameter
    <surface_reflection>    "lambertian" 0.0

    # enable peel off
    <peel_off>    1
    # enable enforced first scattering
    <enfsca>    1

</task>
