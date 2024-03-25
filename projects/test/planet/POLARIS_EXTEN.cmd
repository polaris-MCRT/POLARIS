<common>

    <phase_function>    PH_RAYLEIGH

    <nr_threads>    -1

    <source_extended_star nr_photons = "1e6">    0 -14959787070 0 1.0 6000 1

    <cmd>    CMD_PLANET_SCATTERING

    <peel_off>    1
    <enfsca>    1

    # rotation axis for first and second rotation angle
    <axis1>    1    0    0
    <axis2>    0    0    1

</common>

<task> 1

    <gas_component>    "input/cross_sections/rayleigh_1.0.dat"

    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 0.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 10.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 20.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 30.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 40.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 50.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 60.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 70.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 80.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 90.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 100.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 110.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 120.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 130.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 140.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 150.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 160.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 170.0 3.0856775814671917e+17

    <path_grid>    "projects/test/planet/no_atmosphere.dat"
    <path_out>    "projects/test/planet/lambertian_exten/"

    <surface_reflection>    "lambertian" 1.0

</task>

<task> 1

    <gas_component>    "input/cross_sections/rayleigh_1.0.dat"

    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 0.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 10.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 20.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 30.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 40.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 50.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 60.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 70.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 80.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 90.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 100.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 110.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 120.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 130.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 140.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 150.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 160.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 170.0 3.0856775814671917e+17

    <path_grid>    "projects/test/planet/no_atmosphere.dat"
    <path_out>    "projects/test/planet/lommelseeliger_exten/"

    <surface_reflection>    "lommelseeliger" 1.0

</task>

<task> 1

    <gas_component>    "input/cross_sections/rayleigh_1.0.dat"

    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 0.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 7.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 17.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 27.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 37.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 47.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 57.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 67.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 77.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 87.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 97.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 107.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 117.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 127.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 137.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 147.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 157.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 167.5 3.0856775814671917e+17

    <path_grid>    "projects/test/planet/rayleigh_tau0.3.dat"
    <path_out>    "projects/test/planet/rayleigh_tau0.3_surface1.0_exten/"

    <surface_reflection>    "lambertian" 1.0

</task>

<task> 1

    <gas_component>    "input/cross_sections/rayleigh_1.0.dat"

    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 0.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 7.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 17.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 27.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 37.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 47.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 57.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 67.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 77.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 87.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 97.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 107.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 117.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 127.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 137.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 147.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 157.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 167.5 3.0856775814671917e+17

    <path_grid>    "projects/test/planet/rayleigh_tau5.0.dat"
    <path_out>    "projects/test/planet/rayleigh_tau5.0_surface0.0_exten/"

    <surface_reflection>    "lambertian" 0.0

</task>

<task> 1

    <gas_component>    "input/cross_sections/rayleigh_0.8.dat"

    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 0.0 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 7.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 17.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 27.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 37.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 47.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 57.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 67.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 77.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 87.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 97.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 107.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 117.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 127.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 137.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 147.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 157.5 3.0856775814671917e+17
    <detector_dust_mc nr_pixel = "64">    5.5e-07 5.5e-07 1 90.0 167.5 3.0856775814671917e+17

    <path_grid>    "projects/test/planet/rayleigh_tau10.0.dat"
    <path_out>    "projects/test/planet/rayleigh_tau10.0_surface0.3_exten/"

    <surface_reflection>    "lambertian" 0.3

</task>
