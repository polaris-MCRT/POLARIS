# common commands to all tasks
<common>

    # set the considered gaseous compositions
    # "path_to_file" fraction
    <gas_component>    "input/cross_sections/molecular_hydrogen.dat" 1.0
    <gas_component>    "input/cross_sections/methane_k94.dat" 0.01
    
    # phase function for scattering
    <phase_function>    PH_RAYLEIGH

    # number of CPU cores that will be used by POLARIS
    <nr_threads>    -1

</common>

<task> 1

    # rotation axis for first and second rotation angle
    <axis1>    1    0    0
    <axis2>    0    0    1

    # detector for dust scattering = "N_pixel"> wavelength_min wavelength_max wavelength_nr rot_1 rot_2 distance (sidelength_x sidelength_y)
    # use one detector for each phase angle
    <detector_dust_mc nr_pixel = "64">    4.0e-07 1.0e-06 61 90.0 45.0 3.0856775814671917e+17

    # spatially extended star as radiation source
    # star = "N_photons"> pos_x pos_y pos_z radius temperature biased_emission
    <source_extended_star nr_photons = "1e7">    0 -14959787070 0 1.0 6000 1

    # Monte Carlo dust scattering in planetary atmosphere
    <cmd>    CMD_PLANET_SCATTERING

    # path of the input grid file
    <path_grid>    "projects/rayleigh/rayleigh.dat"
    # path for all output data
    <path_out>    "projects/methane/"

    # surface reflection
    # "type" parameter
    <surface_reflection>    "lambertian" 0.0

    # enable peel off
    <peel_off>    1
    # enable enforced first scattering
    <enfsca>    1

</task>
