#include "Dust.h"
#include "CommandParser.h"
#include "Grid.h"
#include "MathFunctions.h"
#include "typedefs.h"
#include <complex>

void CDustComponent::initDustProperties()
{
    // Init pointer arrays for grain sizes and size distributions
    a_eff = new double[nr_of_dust_species];
    a_eff_1_5 = new double[nr_of_dust_species];
    a_eff_3_5 = new double[nr_of_dust_species];
    a_eff_2 = new double[nr_of_dust_species];
    mass = new double[nr_of_dust_species];

    // Init pointer arrays for dust optical properties
    Qext1 = new double *[nr_of_dust_species];
    Qext2 = new double *[nr_of_dust_species];
    Qabs1 = new double *[nr_of_dust_species];
    Qabs2 = new double *[nr_of_dust_species];
    Qsca1 = new double *[nr_of_dust_species];
    Qsca2 = new double *[nr_of_dust_species];
    Qcirc = new double *[nr_of_dust_species];
    HGg = new double *[nr_of_dust_species];

    // Set pointer array values to zero or add second dimension
    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        a_eff[a] = 0;
        a_eff_1_5[a] = 0;
        a_eff_3_5[a] = 0;
        a_eff_2[a] = 0;
        mass[a] = 0;

        Qext1[a] = new double[nr_of_wavelength];
        Qext2[a] = new double[nr_of_wavelength];
        Qabs1[a] = new double[nr_of_wavelength];
        Qabs2[a] = new double[nr_of_wavelength];
        Qsca1[a] = new double[nr_of_wavelength];
        Qsca2[a] = new double[nr_of_wavelength];
        Qcirc[a] = new double[nr_of_wavelength];
        HGg[a] = new double[nr_of_wavelength];

        // Set pointer array values to zero
        for(uint w = 0; w < nr_of_wavelength; w++)
        {
            Qext1[a][w] = 0;
            Qext2[a][w] = 0;
            Qabs1[a][w] = 0;
            Qabs2[a][w] = 0;
            Qsca1[a][w] = 0;
            Qsca2[a][w] = 0;
            Qcirc[a][w] = 0;
            HGg[a][w] = 0;
        }
    }

    // Qtrq and HG_g_Factor are splines over the incident angle
    Qtrq = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g_factor = new spline[nr_of_dust_species * nr_of_wavelength];
}

void CDustComponent::initScatteringMatrixArray()
{
    // Init counter and percentage to show progress
    ullong per_counter = 0;
    float last_percentage = 0;

    // Init maximum counter value
    uint max_counter = nr_of_wavelength * nr_of_dust_species;

    // First init of the scattering matrix spline
    sca_mat = new Matrix2D ****[nr_of_dust_species];

#pragma omp parallel for
    for(int a = 0; a < int(nr_of_dust_species); a++)
    {
        // Second init of the scattering matrix
        sca_mat[a] = new Matrix2D ***[nr_of_wavelength];
        for(int w = 0; w < nr_of_wavelength; w++)
        {
            // Increase counter used to show progress
            per_counter++;

            // Calculate percentage of total progress per source
            float percentage = 100.0 * float(per_counter) / float(max_counter);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
#pragma omp critical
                {
                    printIDs();
                    cout << "- allocating memory: " << percentage << " [%]                      \r";
                    last_percentage = percentage;
                }
            }

            // Third init of the scattering matrix
            sca_mat[a][w] = new Matrix2D **[nr_of_incident_angles];
            for(int inc = 0; inc < nr_of_incident_angles; inc++)
            {
                // Fourth init of the scattering matrix
                sca_mat[a][w][inc] = new Matrix2D *[nr_of_scat_phi];
                for(int sph = 0; sph < nr_of_scat_phi; sph++)
                {
                    // Fifth init of the scattering matrix
                    sca_mat[a][w][inc][sph] = new Matrix2D[nr_of_scat_theta];
                    for(int sth = 0; sth < nr_of_scat_theta; sth++)
                    {
                        // sixth init of the scattering matrix
                        sca_mat[a][w][inc][sph][sth].resize(4, 4);
                    }
                }
            }
        }
    }
}

void CDustComponent::initCalorimetry()
{
    // Get calorimetric temperatures from first component
    calorimetry_temperatures = new double[nr_of_calorimetry_temperatures];
    for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
        calorimetry_temperatures[t] = 0;

    // Init 2D array for the enthalpy
    enthalpy = new double *[nr_of_dust_species];

    // Add second dimension
    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        enthalpy[a] = new double[nr_of_calorimetry_temperatures];
        for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
            enthalpy[a][t] = 0;
    }
}

bool CDustComponent::readDustParameterFile(parameters & param, uint dust_component_choice)
{
    // Get Path to dust parameters file
    string path = param.getDustPath(dust_component_choice);

    // Init variables
    CCommandParser ps;
    unsigned char ru[4] = { '|', '/', '-', '\\' };
    string line;
    dlist values, wavelength_list_dustcat;

    // temporary variables for wavelength interpolation
    spline *eff_wl, *Qtrq_wl, *HG_g_factor_wl;
    uint nr_of_wavelength_dustcat;

    // Get min and max dust grain size
    double a_min = param.getSizeMin(dust_component_choice);
    double a_max = param.getSizeMax(dust_component_choice);

    // Init text file reader for dust cat file
    fstream reader(path.c_str());

    // Error message if the read does not work
    if(reader.fail())
    {
        cout << "\nERROR: Cannot open dust parameters file:" << endl;
        cout << path << endl;
        return false;
    }

    // Init progress counter
    uint line_counter = 0;
    uint char_counter = 0;
    uint cmd_counter = 0;
    uint eff_counter = 0;

    while(getline(reader, line))
    {
        // Show progress
        if(line_counter % 500 == 0)
        {
            char_counter++;
            printIDs();
            cout << "- reading dust parameters file: " << ru[(uint)char_counter % 4] << "             \r";
        }

        // Format the text file line
        ps.formatLine(line);

        // Increase line counter
        line_counter++;

        // If the line is empty -> skip
        if(line.size() == 0)
            continue;

        if(cmd_counter > 0)
        {
            // Parse the values of the current line
            values = ps.parseValues(line);

            // If no values found -> skip
            if(values.size() == 0)
                continue;
        }

        // Increase the command counter
        cmd_counter++;

        switch(cmd_counter)
        {
            case 1:
                // The first line contains the name of the dust component
                stringID = line;
                break;

            case 2:
                // The second line needs 8 values
                if(values.size() != 8)
                {
                    cout << "\nERROR: Line " << line_counter << " wrong amount of numbers!" << endl;
                    return false;
                }

                // The number of dust grain sizes
                nr_of_dust_species = (uint)values[0];

                // The number of wavelength used by the dust catalog
                nr_of_wavelength_dustcat = (uint)values[1];

                // The number of incident angles
                nr_of_incident_angles = (uint)values[2];

                // The aspect ratio between the longer and shorter axis
                aspect_ratio = values[3];

                // The material density (only used if no one was set in the command file)
                if(material_density == 0)
                {
                    if(values[4] == 0)
                    {
                        printIDs();
                        cout << "\nERROR: dust bulk mass is zero!" << endl;
                        return false;
                    }
                    material_density = values[4];
                }

                // The sublimation temperature
                sub_temp = values[5];

                // The delta value fot the RAT alignment theory
                delta_rat = values[6];

                // The boolean value if the dust catalog is aligned or not
                if(values[7] == 1)
                    is_align = true;
                else
                    is_align = false;

                // Calculate the GOLD alignment g factor
                gold_g_factor = 0.5 * (aspect_ratio * aspect_ratio - 1);

                // Init splines for wavelength interpolation of the dust optical
                // properties
                eff_wl = new spline[nr_of_dust_species * NR_OF_EFF];
                Qtrq_wl = new spline[nr_of_dust_species * nr_of_incident_angles];
                HG_g_factor_wl = new spline[nr_of_dust_species * nr_of_incident_angles];

                break;

            case 3:
                // The second line needs a values per dust grain size
                if(values.size() != nr_of_dust_species)
                {
                    cout << "\nERROR: Line " << line_counter << " wrong amount of effective radii!" << endl;
                    return false;
                }

                // Init arrays for grain size, size distribution, and mass
                a_eff = new double[nr_of_dust_species];
                a_eff_2 = new double[nr_of_dust_species];
                a_eff_1_5 = new double[nr_of_dust_species];
                a_eff_3_5 = new double[nr_of_dust_species];
                mass = new double[nr_of_dust_species];

                // Calculate the grain size distribution
                calcSizeDistribution(values, mass);

                // Check if size limits are inside grain sizes and set global ones
                if(!checkGrainSizeLimits(a_min, a_max))
                    return false;
                break;

            case 4:
                // The fourth line needs a values per wavelength of the dust grain catalog
                if(values.size() != nr_of_wavelength_dustcat)
                {
                    cout << "\nERROR: Line " << line_counter << " wrong amount of wavelength!" << endl;
                    return false;
                }

                // Init an array for the wavelengths of the dust grain catalog
                wavelength_list_dustcat.resize(nr_of_wavelength_dustcat);

                // Fill the wavelength array
                for(uint w = 0; w < nr_of_wavelength_dustcat; w++)
                    wavelength_list_dustcat[w] = values[w];

                if(wavelength_list[0] < wavelength_list_dustcat[0] * 0.1 ||
                   wavelength_list[nr_of_wavelength - 1] * 0.1 >
                       wavelength_list_dustcat[nr_of_wavelength_dustcat - 1])
                {
                    cout << "\nHINT: The wavelength range is out of the limits of the "
                            "catalog. "
                         << "This may cause problems!" << endl;
                    cout << "      wavelength range          : " << wavelength_list[0] << " [m] to "
                         << wavelength_list[nr_of_wavelength - 1] << " [m]" << endl;
                    cout << "      wavelength range (catalog): " << wavelength_list_dustcat[0] << " [m] to "
                         << wavelength_list_dustcat[nr_of_wavelength_dustcat - 1] << " [m]" << endl;
                }

                break;

            default:
                // Init boolean value to check if the lines contain the right values
                bool rec = false;

                // Get the wavelength and grain size indizes
                uint w = int(eff_counter / nr_of_dust_species);
                uint a = eff_counter % nr_of_dust_species;

                // At the first wavelength, resize the splines for the wavelengths
                if(w == 0)
                    for(uint i = 0; i < NR_OF_EFF - 1; i++)
                        eff_wl[a * NR_OF_EFF + i].resize(nr_of_wavelength_dustcat);

                // Each line contains NR_OF_EFF - 1 plus 2 times nr_of_incident_angles
                // values
                if(values.size() == NR_OF_EFF - 1 + 2 * nr_of_incident_angles)
                {
                    // he current line contains enough values
                    rec = true;

                    // Set the dust grain optical properties
                    for(uint i = 0; i < NR_OF_EFF - 1; i++)
                        eff_wl[a * NR_OF_EFF + i].setValue(w, wavelength_list_dustcat[w], values[i]);

                    // For each incident angles, get Qtrq and HG g factor
                    for(uint i_inc = 0; i_inc < nr_of_incident_angles; i_inc++)
                    {
                        // At the first wavelength, resize the splines for the wavelengths
                        if(w == 0)
                        {
                            Qtrq_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g_factor_wl[a * nr_of_incident_angles + i_inc].resize(
                                nr_of_wavelength_dustcat);
                        }

                        // Set the radiative torque efficiency
                        Qtrq_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w, wavelength_list_dustcat[w], values[NR_OF_EFF - 1 + i_inc]);

                        // Set the Henyey-Greenstein factor
                        HG_g_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            values[NR_OF_EFF - 1 + i_inc + nr_of_incident_angles]);
                    }

                    // Increace the counter for the dust grain optical properties
                    eff_counter++;
                }

                // At the last wavelength, activate the splines
                if(w == nr_of_wavelength_dustcat - 1)
                {
                    // For the dust grain optical properties
                    for(uint i = 0; i < NR_OF_EFF - 1; i++)
                        eff_wl[a * NR_OF_EFF + i].createSpline();

                    // For the Qtrq and HG g factor for each incident angle
                    for(uint i_inc = 0; i_inc < nr_of_incident_angles; i_inc++)
                    {
                        Qtrq_wl[a * nr_of_incident_angles + i_inc].createSpline();
                        HG_g_factor_wl[a * nr_of_incident_angles + i_inc].createSpline();
                    }
                }

                //  If a line was not correct, show error
                if(!rec)
                {
                    cout << "WARING: Wrong amount of values in line " << line_counter << "!" << endl;
                    return false;
                }
                break;
        }
    }

    // If not a line per combination of grain size and wavelength was found in the
    // catalog, show error
    if(eff_counter != nr_of_dust_species * nr_of_wavelength_dustcat)
    {
        cout << stringID << endl;
        cout << "\nERROR: Wrong amount of efficiencies in file!" << endl;
        return false;
    }

    // Close the text file reader for the dust catalog
    reader.close();

    // Init pointer arrays for dust optical properties
    Qext1 = new double *[nr_of_dust_species];
    Qext2 = new double *[nr_of_dust_species];
    Qabs1 = new double *[nr_of_dust_species];
    Qabs2 = new double *[nr_of_dust_species];
    Qsca1 = new double *[nr_of_dust_species];
    Qsca2 = new double *[nr_of_dust_species];
    Qcirc = new double *[nr_of_dust_species];
    HGg = new double *[nr_of_dust_species];

    // Init splines for incident angle interpolation of Qtrq and HG g factor
    Qtrq = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g_factor = new spline[nr_of_dust_species * nr_of_wavelength];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        // Resize the splines for each grain size
        Qext1[a] = new double[nr_of_wavelength];
        Qext2[a] = new double[nr_of_wavelength];
        Qabs1[a] = new double[nr_of_wavelength];
        Qabs2[a] = new double[nr_of_wavelength];
        Qsca1[a] = new double[nr_of_wavelength];
        Qsca2[a] = new double[nr_of_wavelength];
        Qcirc[a] = new double[nr_of_wavelength];
        HGg[a] = new double[nr_of_wavelength];

        for(uint w = 0; w < nr_of_wavelength; w++)
        {
            // Resize the splines of Qtrq and HG g factor for each wavelength
            Qtrq[w * nr_of_dust_species + a].resize(nr_of_incident_angles);
            HG_g_factor[w * nr_of_dust_species + a].resize(nr_of_incident_angles);

            if(sizeIndexUsed(a))
            {
                // Calculate the difference between two incident angles
                double d_ang;
                if(nr_of_incident_angles > 1)
                    d_ang = PI / double(nr_of_incident_angles - 1);
                else
                    d_ang = 1;

                // Set the splines of Qtrq and HG g factor for each incident angle
                for(uint i_inc = 0; i_inc < nr_of_incident_angles; i_inc++)
                {
                    Qtrq[w * nr_of_dust_species + a].setValue(
                        i_inc,
                        i_inc * d_ang,
                        Qtrq_wl[a * nr_of_incident_angles + i_inc].getValue(wavelength_list[w]));
                    HG_g_factor[w * nr_of_dust_species + a].setValue(
                        i_inc,
                        i_inc * d_ang,
                        HG_g_factor_wl[a * nr_of_incident_angles + i_inc].getValue(wavelength_list[w],
                                                                                   CONST));
                }

                // Calculate the average HG g factor over all angles
                double avg_HG_g_factor = HG_g_factor[w * nr_of_dust_species + a].getAverageY();

                // Show error if the g factor is lower than -1 and larger than 1
                if(avg_HG_g_factor > 1 || avg_HG_g_factor < -1)
                {
                    cout << "\nERROR: Henyey-Greenstein g factor is smaller than -1 or "
                            "larger than 1!";
                    return false;
                }

                // Set the splines of the dust grain optical properties
                if(isAligned())
                {
                    Qext1[a][w] = eff_wl[a * NR_OF_EFF + 0].getValue(wavelength_list[w]);
                    Qext2[a][w] = eff_wl[a * NR_OF_EFF + 1].getValue(wavelength_list[w]);
                    Qabs1[a][w] = eff_wl[a * NR_OF_EFF + 2].getValue(wavelength_list[w]);
                    Qabs2[a][w] = eff_wl[a * NR_OF_EFF + 3].getValue(wavelength_list[w]);
                    Qsca1[a][w] = eff_wl[a * NR_OF_EFF + 4].getValue(wavelength_list[w]);
                    Qsca2[a][w] = eff_wl[a * NR_OF_EFF + 5].getValue(wavelength_list[w]);
                    Qcirc[a][w] = eff_wl[a * NR_OF_EFF + 6].getValue(wavelength_list[w]);
                }
                else
                {
                    double tmpQext = 1.0 / 3.0 *
                                     (2.0 * eff_wl[a * NR_OF_EFF + 0].getValue(wavelength_list[w]) +
                                      eff_wl[a * NR_OF_EFF + 1].getValue(wavelength_list[w]));
                    double tmpQabs = 1.0 / 3.0 *
                                     (2.0 * eff_wl[a * NR_OF_EFF + 2].getValue(wavelength_list[w]) +
                                      eff_wl[a * NR_OF_EFF + 3].getValue(wavelength_list[w]));
                    double tmpQsca = 1.0 / 3.0 *
                                     (2.0 * eff_wl[a * NR_OF_EFF + 4].getValue(wavelength_list[w]) +
                                      eff_wl[a * NR_OF_EFF + 5].getValue(wavelength_list[w]));

                    Qext1[a][w] = tmpQext;
                    Qext2[a][w] = tmpQext;
                    Qabs1[a][w] = tmpQabs;
                    Qabs2[a][w] = tmpQabs;
                    Qsca1[a][w] = tmpQsca;
                    Qsca2[a][w] = tmpQsca;
                    Qcirc[a][w] = 0;
                }
                HGg[a][w] = avg_HG_g_factor;
            }
            else
            {
                Qext1[a][w] = 0;
                Qext2[a][w] = 0;
                Qabs1[a][w] = 0;
                Qabs2[a][w] = 0;
                Qsca1[a][w] = 0;
                Qsca2[a][w] = 0;
                Qcirc[a][w] = 0;
                HGg[a][w] = 0;
            }

            // Activate the splines of Qtrq and HG g factor
            Qtrq[w * nr_of_dust_species + a].createSpline();
            HG_g_factor[w * nr_of_dust_species + a].createSpline();
        }
    }

    // Read the scattering matrix if MIE scattering should be used
    // With the same grid of wavelengths and grain sizes
    if(param.getPhaseFunctionID() == PH_MIE)
        if(!readScatteringMatrices(path, nr_of_wavelength_dustcat, wavelength_list_dustcat))
            return false;

    // Remove temporary pointer arrays
    delete[] eff_wl;
    delete[] Qtrq_wl;
    delete[] HG_g_factor_wl;

    return true;
}

bool CDustComponent::readDustRefractiveIndexFile(parameters & param,
                                                 uint dust_component_choice,
                                                 double a_min_mixture,
                                                 double a_max_mixture)
{
    // Init variables
    CCommandParser ps;
    unsigned char ru[4] = { '|', '/', '-', '\\' };
    string line;
    dlist values, values_aeff;

    // temporary variables for wavelength interpolation
    spline refractive_index_real, refractive_index_imag;
    uint nr_of_wavelength_dustcat;

    // Get min and max dust grain size
    double a_min = param.getSizeMin(dust_component_choice);
    double a_max = param.getSizeMax(dust_component_choice);

    // Set number of grain sizes for Mie-theory (1 if only one grain size is used)
    if(a_min_mixture == a_max_mixture)
        nr_of_dust_species = 1;
    else
        nr_of_dust_species = MIE_SIZE_STEPS;
    values_aeff.resize(nr_of_dust_species);

    // Init dust grain sizes
    CMathFunctions::LogList(a_min_mixture, a_max_mixture, values_aeff, 10);

    // Get Path to dust parameters file
    string path = param.getDustPath(dust_component_choice);

    // Init text file reader for dust cat file
    fstream reader(path.c_str());

    // Error message if the read does not work
    if(reader.fail())
    {
        cout << "\nERROR: Cannot open dust refractive index file:" << endl;
        cout << path << endl;
        return false;
    }

    // Init progress counter
    uint line_counter = 0;
    uint char_counter = 0;
    uint cmd_counter = 0;
    uint wl_counter = 0;

    while(getline(reader, line))
    {
        // Show progress
        if(line_counter % 500 == 0)
        {
            char_counter++;
            printIDs();
            cout << "- reading dust parameters file: " << ru[(uint)char_counter % 4] << "             \r";
        }

        // Format the text file line
        ps.formatLine(line);

        // Increase line counter
        line_counter++;

        // If the line is empty -> skip
        if(line.size() == 0)
            continue;

        if(cmd_counter > 0)
        {
            // Parse the values of the current line
            values = ps.parseValues(line);

            // If no values found -> skip
            if(values.size() == 0)
                continue;
        }

        // Increase the command counter
        cmd_counter++;

        switch(cmd_counter)
        {
            case 1:
                // The first line contains the name of the dust component
                stringID = line;
                break;

            case 2:
                // The second line needs 7 values
                if(values.size() != 7)
                {
                    cout << "\nERROR: Line " << line_counter << " wrong amount of numbers!" << endl;
                    return false;
                }

                // The number of wavelength used by the dust catalog
                nr_of_wavelength_dustcat = (uint)values[0];

                // The number of incident angles
                nr_of_incident_angles = 1; // For non-spherical: (uint) values[1];

                // The aspect ratio between the longer and shorter axis
                aspect_ratio = 1; // For non-spherical: values[2];

                // The material density (only used if no one was set in the command file)
                if(material_density == 0)
                    material_density = values[3];

                // The sublimation temperature
                sub_temp = values[4];

                // The delta value fot the RAT alignment theory
                delta_rat = 0; // For non-spherical: values[5];

                // The boolean value if the dust catalog is aligned or not
                if(values[6] == 1)
                    is_align = true;
                else
                    is_align = false;

                // Calculate the GOLD alignment g factor
                gold_g_factor = 0.5 * (aspect_ratio * aspect_ratio - 1);

                // Init splines for wavelength interpolation of the dust optical
                // properties
                refractive_index_real.resize(nr_of_wavelength_dustcat);
                refractive_index_imag.resize(nr_of_wavelength_dustcat);

                // Set size parameters
                a_eff = new double[nr_of_dust_species];
                a_eff_2 = new double[nr_of_dust_species];
                a_eff_1_5 = new double[nr_of_dust_species];
                a_eff_3_5 = new double[nr_of_dust_species];
                mass = new double[nr_of_dust_species];

                // Calculate the grain size distribution
                calcSizeDistribution(values_aeff, mass);

                // Check if size limits are inside grain sizes and set global ones
                if(!checkGrainSizeLimits(a_min, a_max))
                    return false;
                break;

            default:
                // Init boolean value to check if the lines contain the right values
                bool rec = false;

                // Each line contains NR_OF_EFF - 1 plus 2 times nr_of_incident_angles
                // values
                if(values.size() == 3)
                {
                    // he current line contains enough values
                    rec = true;

                    // Set the dust grain optical properties
                    refractive_index_real.setValue(wl_counter, values[0], values[1]);
                    refractive_index_imag.setValue(wl_counter, values[0], values[2]);

                    // Increace the counter for the dust grain optical properties
                    wl_counter++;
                }

                //  If a line was not correct, show error
                if(!rec)
                {
                    cout << "WARING: Wrong amount of values in line " << line_counter << "!" << endl;
                    return false;
                }
                break;
        }
    }

    // At the last wavelength, activate the splines
    refractive_index_real.createSpline();
    refractive_index_imag.createSpline();

    // If not a line per combination of grain size and wavelength was found in the
    // catalog, show error
    if(wl_counter != nr_of_wavelength_dustcat)
    {
        cout << stringID << endl;
        cout << "\nERROR: Wrong amount of efficiencies in file!" << endl;
        return false;
    }

    // Close the text file reader for the dust catalog
    reader.close();

    // Init pointer arrays for dust optical properties
    Qext1 = new double *[nr_of_dust_species];
    Qext2 = new double *[nr_of_dust_species];
    Qabs1 = new double *[nr_of_dust_species];
    Qabs2 = new double *[nr_of_dust_species];
    Qsca1 = new double *[nr_of_dust_species];
    Qsca2 = new double *[nr_of_dust_species];
    Qcirc = new double *[nr_of_dust_species];
    HGg = new double *[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        Qext1[a] = new double[nr_of_wavelength];
        Qext2[a] = new double[nr_of_wavelength];
        Qabs1[a] = new double[nr_of_wavelength];
        Qabs2[a] = new double[nr_of_wavelength];
        Qsca1[a] = new double[nr_of_wavelength];
        Qsca2[a] = new double[nr_of_wavelength];
        Qcirc[a] = new double[nr_of_wavelength];
        HGg[a] = new double[nr_of_wavelength];
    }

    // Init splines for incident angle interpolation of Qtrq and HG g factor
    Qtrq = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g_factor = new spline[nr_of_dust_species * nr_of_wavelength];

    // Set variables for scattering via Mie theory
    nr_of_scat_mat_elements = 4;

    // --- Phi angle
    nr_of_scat_phi = 1;

    // --- Theta angle
    nr_of_scat_theta = 2 * NANG - 1;

    // Init normal scattering matrix array
    initScatteringMatrixArray();

    // Init counter and percentage to show progress
    ullong per_counter = 0;
    float last_percentage = 0;

    // Init error check
    bool error = false;

    // Init maximum counter value
    uint max_counter = nr_of_dust_species * nr_of_wavelength;

#pragma omp parallel for schedule(dynamic) collapse(2)
    for(int a = 0; a < int(nr_of_dust_species); a++)
    {
        for(int w = 0; w < int(nr_of_wavelength); w++)
        {
            // Skip everything else if error was found
            if(error)
                continue;

            // Increase counter used to show progress
            per_counter++;

            // Calculate percentage of total progress per source
            float percentage = 100.0 * float(per_counter) / float(max_counter);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
#pragma omp critical
                {
                    printIDs();
                    cout << "- calculating optical properties: " << percentage
                         << " [%]                      \r";
                    last_percentage = percentage;
                }
            }

            // Resize the splines of Qtrq and HG g factor for each wavelength
            Qtrq[w * nr_of_dust_species + a].resize(nr_of_incident_angles);
            HG_g_factor[w * nr_of_dust_species + a].resize(nr_of_incident_angles);

            if(sizeIndexUsed(a))
            {
                // Init variables and pointer arrays
                double *S11, *S12, *S33, *S34;
                S11 = new double[nr_of_scat_theta];
                S12 = new double[nr_of_scat_theta];
                S33 = new double[nr_of_scat_theta];
                S34 = new double[nr_of_scat_theta];

                // Set size index and refractive index as complex number
                double x = 2.0 * PI * a_eff[a] / wavelength_list[w];
                dcomplex refractive_index = dcomplex(refractive_index_real.getValue(wavelength_list[w]),
                                                     refractive_index_imag.getValue(wavelength_list[w]));

                // Calculate Mie-scattering
                if(!CMathFunctions::calcWVMie(x,
                                              refractive_index,
                                              Qext1[a][w],
                                              Qabs1[a][w],
                                              Qsca1[a][w],
                                              HGg[a][w],
                                              S11,
                                              S12,
                                              S33,
                                              S34))
                    error = true;

                // Set missing Efficiencies for other axis
                Qext2[a][w] = Qext1[a][w];
                Qabs2[a][w] = Qabs1[a][w];
                Qsca2[a][w] = Qsca1[a][w];
                Qcirc[a][w] = 0;

                for(int inc = 0; inc < nr_of_incident_angles; inc++)
                {
                    for(int sph = 0; sph < nr_of_scat_phi; sph++)
                    {
                        for(int sth = 0; sth < nr_of_scat_theta; sth++)
                        {
                            sca_mat[a][w][inc][sph][sth](0, 0) = S11[sth]; // S11
                            sca_mat[a][w][inc][sph][sth](1, 1) = S11[sth]; // S22

                            sca_mat[a][w][inc][sph][sth](0, 1) = S12[sth]; // S12
                            sca_mat[a][w][inc][sph][sth](1, 0) = S12[sth]; // S21

                            sca_mat[a][w][inc][sph][sth](2, 2) = S33[sth]; // S33
                            sca_mat[a][w][inc][sph][sth](3, 3) = S33[sth]; // S44

                            sca_mat[a][w][inc][sph][sth](2, 3) = S34[sth];  // S34
                            sca_mat[a][w][inc][sph][sth](3, 2) = -S34[sth]; // S43
                        }
                    }
                }
                delete[] S11;
                delete[] S12;
                delete[] S33;
                delete[] S34;
            }
            else
            {
                Qext1[a][w] = 0;
                Qext2[a][w] = 0;
                Qabs1[a][w] = 0;
                Qabs2[a][w] = 0;
                Qsca1[a][w] = 0;
                Qsca2[a][w] = 0;
                Qcirc[a][w] = 0;
                HGg[a][w] = 0;
            }

            // Activate the splines of Qtrq and HG g factor
            Qtrq[w * nr_of_dust_species + a].createSpline();
            HG_g_factor[w * nr_of_dust_species + a].createSpline();
        }
    }

    // Set that the scattering matrix was successfully read
    scat_loaded = true;

    if(error)
    {
        cout << "ERROR: Problem with optical properties calculation" << endl;
        return false;
    }
    return true;
}

bool CDustComponent::readScatteringMatrices(string path,
                                            uint nr_of_wavelength_dustcat,
                                            dlist wavelength_list_dustcat)
{
    // Init variables
    string::size_type pos = 0;
    bool disable_mie_scattering = true;
    CCommandParser ps;
    dlist values;
    string line;

    // Erase the ".dat" from the path
    if(path.find(".dat") != string::npos)
    {
        pos = path.find(".dat");
        path.erase(pos, 4);
    }

    // Create the filename of the scattering info file
    string inf_filename = path;
    inf_filename += SEP;
    inf_filename += "scat.inf";

    // Init text file reader for scattering info file
    ifstream inf_reader(inf_filename.c_str());

    // Error message if the read does not work
    if(inf_reader.fail())
    {
        cout << "\nERROR: Cannot open scattering matrix info file:" << endl;
        cout << inf_filename << endl;
        return false;
    }

    // Init progress counter
    uint line_counter = 0;
    uint cmd_counter = 0;

    // Go through each line of the info file
    while(getline(inf_reader, line))
    {
        // Format the text file line
        ps.formatLine(line);

        // Increase line counter
        line_counter++;

        // If the line is empty -> skip
        if(line.size() == 0)
            continue;

        // Parse the values of the current line
        values = ps.parseValues(line);

        // If no values found -> skip
        if(values.size() == 0)
            continue;

        // Increase the command counter
        cmd_counter++;

        switch(cmd_counter)
        {
            case 1:
                // The first line needs 5 values
                if(values.size() != 5)
                {
                    cout << "\nERROR: Wrong amount of dust component parameters in:" << endl;
                    cout << inf_filename.c_str() << " line " << line_counter << "!" << endl;
                    return false;
                }

                // The number of dust grain sizes
                if(values[0] != nr_of_dust_species)
                {
                    cout << "\nERROR: " << inf_filename.c_str() << " line " << line_counter << "!" << endl;
                    cout << "Number of dust species does not match the number in the "
                            "dust parameters file!"
                         << endl;
                    return false;
                }

                // The number of wavelength used by the dust catalog
                if(values[1] != nr_of_wavelength_dustcat)
                {
                    cout << "\nERROR: " << inf_filename.c_str() << " line " << line_counter << "!" << endl;
                    cout << "Number of wavelength does not match the number in the dust "
                            "parameters file!"
                         << endl;
                    return false;
                }

                // The number of incident angles
                if(values[2] != nr_of_incident_angles)
                {
                    cout << "\nERROR: " << inf_filename.c_str() << " line " << line_counter << "!" << endl;
                    cout << "Number of incident angles does not match the number in the "
                            "dust parameters file!"
                         << endl;
                    return false;
                }

                // The number of phi angles (outgoing radiation)
                nr_of_scat_phi = uint(values[3]);

                // The number of theta angles (outgoing radiation)
                nr_of_scat_theta = uint(values[4]);

                break;

            case 2:
                // The second line needs one value
                if(values.size() != 1)
                {
                    cout << "\nERROR: Wrong amount of scattering matrix elements in:" << endl;
                    cout << inf_filename.c_str() << " line " << line_counter << "!" << endl;
                    return false;
                }

                // The number of used scattering elements
                nr_of_scat_mat_elements = uint(values[0]);
                break;

            case 3:
                // The third line needs 16 values
                if(values.size() != 16)
                {
                    cout << "\nERROR: Wrong amount of matrix elements in:" << endl;
                    cout << inf_filename.c_str() << " line 3!" << endl;
                    return false;
                }

                // The relation which scattering matrix entry is used at which position in
                // the 4x4 matrix
                for(uint i = 0; i < 16; i++)
                {
                    elements[i] = int(values[i]);
                    if(elements[i] != 0)
                        disable_mie_scattering = false;
                }
                break;
        }
    }

    // Close the file reader
    inf_reader.close();

    // If scattering matrix is empty, disable mie scattering for the corresponding dust
    // component
    if(disable_mie_scattering)
    {
        phID = PH_HG;
        return true;
    }

    // Init counter and percentage to show progress
    ullong per_counter = 0;
    float last_percentage = 0;

    // Init maximum counter value
    uint max_counter = nr_of_incident_angles * nr_of_dust_species * nr_of_scat_phi;

    // First init of the scattering matrix interp
    interp ***** sca_mat_wl = new interp ****[nr_of_dust_species];

#pragma omp parallel for
    for(int a = 0; a < int(nr_of_dust_species); a++)
    {
        // Second init of the scattering matrix interp
        sca_mat_wl[a] = new interp ***[nr_of_incident_angles];
        for(int inc = 0; inc < nr_of_incident_angles; inc++)
        {
            // Third init of the scattering matrix interp
            sca_mat_wl[a][inc] = new interp **[nr_of_scat_phi];
            for(int sph = 0; sph < nr_of_scat_phi; sph++)
            {
                // Increase counter used to show progress
                per_counter++;

                // Calculate percentage of total progress per source
                float percentage = 100.0 * float(per_counter) / float(max_counter);

                // Show only new percentage number if it changed
                if((percentage - last_percentage) > PERCENTAGE_STEP)
                {
#pragma omp critical
                    {
                        printIDs();
                        cout << "- allocating memory: " << percentage << " [%]                      \r";
                        last_percentage = percentage;
                    }
                }

                // Fourth init of the scattering matrix interp
                sca_mat_wl[a][inc][sph] = new interp *[nr_of_scat_theta];
                for(int sth = 0; sth < nr_of_scat_theta; sth++)
                {
                    // Fifth init of the scattering matrix interp
                    sca_mat_wl[a][inc][sph][sth] = new interp[nr_of_scat_mat_elements];

                    // Resize the linear for each scattering matrix element
                    for(int mat = 0; mat < nr_of_scat_mat_elements; mat++)
                        sca_mat_wl[a][inc][sph][sth][mat].resize(nr_of_wavelength_dustcat);
                }
            }
        }
    }

    // Set counter and percentage to show progress
    per_counter = 0;
    last_percentage = 0;

    // Init error bool
    bool error = false;

    // Set maximum counter value
    max_counter = nr_of_wavelength_dustcat * nr_of_dust_species;

#pragma omp parallel for
    for(int w = 0; w < int(nr_of_wavelength_dustcat); w++)
    {
        if(error)
            continue;

        // Init characters to create the filename for each wavelength
        char str_ID_tmp[32];
        char str_ID_end[32];
        float tmp_val = 0;

        // Set the characters with the current indizes
#ifdef WINDOWS
        strcpy_s(str_ID_tmp, "wID%03d.sca");
        sprintf_s(str_ID_end, str_ID_tmp, w + 1);
#else
        strcpy(str_ID_tmp, "wID%03d.sca");
        sprintf(str_ID_end, str_ID_tmp, w + 1);
#endif

        // Create the filename
        string bin_filename = path;
        bin_filename += SEP;
        bin_filename += str_ID_end;

        // Init text file reader for scattering binary file
        ifstream bin_reader(bin_filename.c_str(), ios::in | ios::binary);

        // Error message if the read does not work
        if(bin_reader.fail())
        {
            cout << "\nERROR: Cannot open scattering matrix file:" << endl;
            cout << bin_filename.c_str() << endl;
            bin_reader.close();
            error = true;
        }

        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            // Increase counter used to show progress
            per_counter++;

            // Calculate percentage of total progress per source
            float percentage = 100.0 * float(per_counter) / float(max_counter);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
#pragma omp critical
                {
                    printIDs();
                    cout << "- loading matrices: " << percentage << " [%]                      \r";
                    last_percentage = percentage;
                }
            }

            for(uint inc = 0; inc < nr_of_incident_angles; inc++)
            {
                for(uint sph = 0; sph < nr_of_scat_phi; sph++)
                {
                    for(uint sth = 0; sth < nr_of_scat_theta; sth++)
                    {
                        for(uint mat = 0; mat < nr_of_scat_mat_elements; mat++)
                        {
                            // Read value from binary file
                            bin_reader.read((char *)&tmp_val, 4);

                            // Set value of interp to the read value
                            sca_mat_wl[a][inc][sph][sth][mat].setValue(
                                w, wavelength_list_dustcat[w], tmp_val);
                        }
                    }
                }
            }
        }
        // Close scattering binary file
        bin_reader.close();
    }

    if(error)
        return false;

    // Init normal scattering matrix array
    initScatteringMatrixArray();

    // Fill values of the scattering matrix via interpolation
#pragma omp parallel for
    for(int a = 0; a < int(nr_of_dust_species); a++)
    {
        if(sizeIndexUsed(a))
        {
            for(int inc = 0; inc < nr_of_incident_angles; inc++)
            {
                for(int sph = 0; sph < nr_of_scat_phi; sph++)
                {
                    for(int sth = 0; sth < nr_of_scat_theta; sth++)
                    {
                        for(uint e = 0; e < 16; e++)
                        {
                            double sign = CMathFunctions::sgn(elements[e]);
                            int pos = abs(elements[e]);
                            if(pos > 0)
                            {
                                for(int w = 0; w < nr_of_wavelength; w++)
                                    sca_mat[a][w][inc][sph][sth](e) =
                                        sign * sca_mat_wl[a][inc][sph][sth][pos].getValue(wavelength_list[w]);
                            }
                        }
                        delete[] sca_mat_wl[a][inc][sph][sth];
                    }
                    delete[] sca_mat_wl[a][inc][sph];
                }
                delete[] sca_mat_wl[a][inc];
            }
        }
        delete[] sca_mat_wl[a];
    }
    delete[] sca_mat_wl;

    // Set that the scattering matrix was successfully read
    scat_loaded = true;

    return true;
}

bool CDustComponent::readCalorimetryFile(parameters & param, uint dust_component_choice)
{
    // Init variables
    string::size_type pos = 0;
    CCommandParser ps;
    dlist values;
    string line;

    // Get Path to dust parameters file
    string path = param.getDustPath(dust_component_choice);

    // Erase the ".dat" from the path
    if(path.find(".dat") != string::npos)
    {
        pos = path.find(".dat");
        path.erase(pos, 4);
    }
    else if(path.find(".nk") != string::npos)
    {
        pos = path.find(".nk");
        path.erase(pos, 4);
    }

    // Create the filename of the calorimetry file
    string calo_filename = path;
    calo_filename += SEP;
    calo_filename += "calorimetry.dat";

    // Init text file reader for scattering info file
    ifstream calo_reader(calo_filename.c_str());

    // Error message if the read does not work
    if(calo_reader.fail())
    {
        cout << "\nERROR: Cannot open calorimetry file:" << endl;
        cout << calo_filename << endl;
        return false;
    }

    // Init progress counter
    uint line_counter = 0;
    uint cmd_counter = 0;

    // Go through each line of the info file
    while(getline(calo_reader, line))
    {
        // Format the text file line
        ps.formatLine(line);

        // Increase line counter
        line_counter++;

        // If the line is empty -> skip
        if(line.size() == 0)
            continue;

        // Parse the values of the current line
        values = ps.parseValues(line);

        // If no values found -> skip
        if(values.size() == 0)
            continue;

        // Increase the command counter
        cmd_counter++;

        switch(cmd_counter)
        {
            case 1:
                // The first line needs one value
                if(values.size() != 1)
                {
                    cout << "\nERROR: Wrong amount of calorimetry temperatures in:" << endl;
                    cout << calo_filename.c_str() << " line " << line_counter << "!" << endl;
                    return false;
                }

                // The amount of calorimetry temperatures
                nr_of_calorimetry_temperatures = values[0];

                // Init array for the calorimetry temperatures
                calorimetry_temperatures = new double[nr_of_calorimetry_temperatures];

                // Init 2D array for the enthalpy
                enthalpy = new double *[nr_of_dust_species];

                // Add second dimension
                for(uint a = 0; a < nr_of_dust_species; a++)
                    enthalpy[a] = new double[nr_of_calorimetry_temperatures];
                break;

            case 2:
                // The second line needs a value per calorimetric temperature
                if(values.size() != nr_of_calorimetry_temperatures)
                {
                    cout << "\nERROR: Wrong calorimetry temperatures in:" << endl;
                    cout << calo_filename.c_str() << " line " << line_counter << "!" << endl;
                    return false;
                }

                // Set the calorimetric temperatures
                for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
                    calorimetry_temperatures[t] = double(values[t]);
                break;

            case 3:
                // The second line needs one value
                if(values.size() != 1)
                {
                    cout << "\nERROR: Wrong calorimetry type in:" << endl;
                    cout << calo_filename.c_str() << " line " << line_counter << "!" << endl;
                    return false;
                }

                // The unit of the calorimetry data
                calorimetry_type = uint(values[0]);

                // Only heat capacity or enthalpy are possible
                if(calorimetry_type != CALO_HEAT_CAP && calorimetry_type != CALO_ENTHALPY)
                {
                    cout << "\nERROR: Wrong calorimetry type in:" << endl;
                    cout << calo_filename.c_str() << " line " << line_counter << "!" << endl;
                    return false;
                }
                break;

            default:
                // The other lines have either one or a value per dust grain size
                if(values.size() != 1 && values.size() != nr_of_dust_species)
                {
                    cout << "\nERROR: Wrong amount of dust species in:" << endl;
                    cout << calo_filename.c_str() << " line " << line_counter << "!" << endl;
                    return false;
                }

                // Get temperature index
                uint t = cmd_counter - 4;

                for(uint a = 0; a < nr_of_dust_species; a++)
                {
                    // Use either a value per grain size or one value for all
                    double tmp_value;
                    if(values.size() == 1)
                        tmp_value = double(values[0]);
                    else
                        tmp_value = double(values[a]);

                    // If heat capacity, perform integration
                    if(calorimetry_type == CALO_HEAT_CAP)
                    {
                        if(t == 0)
                            enthalpy[a][t] = tmp_value * calorimetry_temperatures[t];
                        else
                            enthalpy[a][t] =
                                enthalpy[a][t - 1] +
                                tmp_value * (calorimetry_temperatures[t] - calorimetry_temperatures[t - 1]);
                    }
                    else if(calorimetry_type == CALO_ENTHALPY)
                    {
                        // Enthalpy is already in the right unit
                        enthalpy[a][t] = tmp_value;
                    }
                    else
                    {
                        // Reset enthalpy if wrong type
                        enthalpy[a][t] = 0;
                    }
                }
                break;
        }
    }
    // Close calorimetry file reader
    calo_reader.close();

    // Multiply the specific enthalpy with the grain size to get the enthalpy
    for(uint a = 0; a < nr_of_dust_species; a++)
        for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
            enthalpy[a][t] *= 4.0 / 3.0 * PI * a_eff[a] * a_eff[a] * a_eff[a];

    // Set that the calorimetry file was successfully loaded
    calorimetry_loaded = true;

    return true;
}

bool CDustComponent::writeComponent(string path_data, string path_plot)
{
    // Do not write dust data if no dust component was chosen
    if(nr_of_wavelength == 0)
        return true;

    // Init character variables to store filenames
    char str_comp_ID_tmp[1024];
    char str_comp_ID_end[1024];
    char str_mix_ID_tmp[1024];
    char str_mix_ID_end[1024];
    char str_frac_tmp[1024];
    char str_frac_end[1024];

    // Init strings for various filenames/titles
    string path_cross, path_eff, path_kappa, path_diff, path_g;
    string path_scat, str_title, gnu_title, plot_sign = "points";

    // Check if enough points to draw lines
    if(nr_of_wavelength > 1)
        plot_sign = "lines";

        // Set the characters with the current indizes
#ifdef WINDOWS
    strcpy_s(str_comp_ID_tmp, "%03d");
    sprintf_s(str_comp_ID_end, str_comp_ID_tmp, i_component + 1);

    strcpy_s(str_mix_ID_tmp, "%03d");
    sprintf_s(str_mix_ID_end, str_mix_ID_tmp, i_mixture + 1);

    strcpy_s(str_frac_tmp, "%.05f");
    sprintf_s(str_frac_end, str_frac_tmp, fraction);
#else
    strcpy(str_comp_ID_tmp, "%03d");
    sprintf(str_comp_ID_end, str_comp_ID_tmp, i_component + 1);

    strcpy(str_mix_ID_tmp, "%03d");
    sprintf(str_mix_ID_end, str_mix_ID_tmp, i_mixture + 1);

    strcpy(str_frac_tmp, "%.05f");
    sprintf(str_frac_end, str_frac_tmp, fraction);
#endif

    if(is_mixture)
    {
        // For the dust mixture, set the different strings
        str_title = "#Dust mixture\n " + stringID;
        gnu_title = "#Dust mixture ";
        gnu_title += str_mix_ID_end;
        gnu_title += "\n";
        path_cross = path_plot + "dust_mixture_" + str_mix_ID_end + "_cross.plt";
        path_eff = path_plot + "dust_mixture_" + str_mix_ID_end + "_eff.plt";
        path_kappa = path_plot + "dust_mixture_" + str_mix_ID_end + "_kappa.plt";
        path_diff = path_plot + "dust_mixture_" + str_mix_ID_end + "_diff.plt";
        path_data += "dust_mixture_";
        path_data += str_mix_ID_end;
        path_data += ".dat";
        path_g = path_plot + "dust_mixture_" + str_mix_ID_end + "_g.plt";
        path_scat = path_plot + "dust_mixture_" + str_mix_ID_end + "_scat.plt";

        // Format the strings
        uint pos = 0;
        while(str_title.find("\n ") != string::npos)
        {
            pos = uint(str_title.find("\n "));
            str_title.replace(pos, 2, "\n#");
        }

        pos = 0;
        while(str_title.find("- ") != string::npos)
        {
            pos = uint(str_title.find("- "));
            str_title.replace(pos, 2, "#-");
        }

        pos = 0;
        while(gnu_title.find("\n") != string::npos)
        {
            pos = uint(gnu_title.find("\n"));
            gnu_title.replace(pos, 1, "; ");
        }

        pos = 0;
        while(gnu_title.find("#") != string::npos)
        {
            pos = uint(gnu_title.find("#"));
            gnu_title.replace(pos, 1, " ");
        }
    }
    else
    {
        // For the dust components, set the different strings
        str_title = "Dust mixture ";
        str_title += str_mix_ID_end;
        str_title += ", component nr.: ";
        str_title += str_comp_ID_end;
        str_title += " mat.: " + stringID;
        str_title += " frac.: ";
        str_title += str_frac_end;
        gnu_title = str_title;
        path_cross = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_cross.plt";
        path_eff = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_eff.plt";
        path_kappa = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_kappa.plt";
        path_diff = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_diff.plt";
        path_g = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_g.plt";
        path_scat = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_scat.plt";
        path_data += "dust_comp_";
        path_data += str_comp_ID_end;
        path_data += ".dat";
    }

    // Init text file writer for cross-sections
    ofstream cross_writer(path_cross.c_str());

    // Error message if the write does not work
    if(cross_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n" << path_cross << endl;
        return false;
    }

    // Init plot limits
    double Cmin = 1e100, Cmax = -1e100;

    // Find min and max values over the wavelength (checks only Cext1 and Cabs1)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
    {
        if(getCext1(i) < Cmin && getCext1(i) > 0)
            Cmin = getCext1(i);
        if(getCabs1(i) < Cmin && getCabs1(i) > 0)
            Cmin = getCabs1(i);
        if(getCsca1(i) < Cmin && getCsca1(i) > 0)
            Cmin = getCsca1(i);
        if(getCext2(i) < Cmin && getCext2(i) > 0)
            Cmin = getCext2(i);
        if(getCabs2(i) < Cmin && getCabs2(i) > 0)
            Cmin = getCabs2(i);
        if(getCsca2(i) < Cmin && getCsca2(i) > 0)
            Cmin = getCsca2(i);

        if(getCext1(i) > Cmax)
            Cmax = getCext1(i);
        if(getCabs1(i) > Cmax)
            Cmax = getCabs1(i);
        if(getCsca1(i) > Cmax)
            Cmax = getCsca1(i);
        if(getCext2(i) > Cmax)
            Cmax = getCext2(i);
        if(getCabs2(i) > Cmax)
            Cmax = getCabs2(i);
        if(getCsca2(i) > Cmax)
            Cmax = getCsca2(i);
    }

    // Add a bit more space for good visualization
    Cmin *= 0.9;
    Cmax *= 1.10;

    // Add Gnuplot commands to file
    cross_writer << "reset" << endl;
    if(nr_of_wavelength > 1)
        cross_writer << "set log x" << endl;
    cross_writer << "set log y" << endl;
    cross_writer << "set grid" << endl;
    if(nr_of_wavelength > 1)
        cross_writer << "set xrange[" << wavelength_list[wavelength_offset] << ":" << wavelength_list.back()
                     << "]" << endl;
    cross_writer << "set yrange[" << Cmin << ":" << Cmax << "]" << endl;
    cross_writer << "set format x \"%.1te%02T\"" << endl;
    cross_writer << "set format y \"%.1te%02T\"" << endl;
    cross_writer << "set ylabel \'C_{avg} [m^{2}]\'" << endl;
    cross_writer << "set xlabel \'{/Symbol l} [m]\'" << endl;
    cross_writer << "set title \"" << gnu_title << "\"" << endl;
    cross_writer << "plot \'-\' with " << plot_sign << " title \'C_{ext,x}\' lc rgb \"#0000F0\","
                 << "\'-\' with " << plot_sign << " title \'C_{ext,y}\' lc rgb \"#000090\","
                 << "\'-\' with " << plot_sign << " title \'C_{abs,x}\' lc rgb \"#FF0000\","
                 << "\'-\' with " << plot_sign << " title \'C_{abs,y}\' lc rgb \"#900000\","
                 << "\'-\' with " << plot_sign << " title \'C_{sca,x}\' lc rgb \"#FFFF00\","
                 << "\'-\' with " << plot_sign << " title \'C_{sca,y}\' lc rgb \"#909000\"" << endl;

    // Add Cext1 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getCext1(i) > 0)
            cross_writer << wavelength_list[i] << "\t" << getCext1(i) << endl;
    cross_writer << "e" << endl;

    // Add Cext2 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getCext2(i) > 0)
            cross_writer << wavelength_list[i] << "\t" << getCext2(i) << endl;
    cross_writer << "e" << endl;

    // Add Cabs1 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getCabs1(i) > 0)
            cross_writer << wavelength_list[i] << "\t" << getCabs1(i) << endl;
    cross_writer << "e" << endl;

    // Add Cabs2 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getCabs2(i) > 0)
            cross_writer << wavelength_list[i] << "\t" << getCabs2(i) << endl;
    cross_writer << "e" << endl;

    // Add Csca1 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getCsca1(i) > 0)
            cross_writer << wavelength_list[i] << "\t" << getCsca1(i) << endl;
    cross_writer << "e" << endl;

    // Add Csca2 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getCsca2(i) > 0)
            cross_writer << wavelength_list[i] << "\t" << getCsca2(i) << endl;
    cross_writer << "e" << endl;

    // Close text file writer
    cross_writer.close();

    // ------------------------------------------------------

    // Init text file writer for efficiencies
    ofstream eff_writer(path_eff.c_str());

    // Error message if the write does not work
    if(eff_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n" << path_eff << endl;
        return false;
    }

    // Init plot limits
    double Qmin = 1e100, Qmax = -1e100;

    // Find min and max values over the wavelength (checks only Cext1 and Cabs1)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
    {
        if(getQext1(i) < Qmin && getQext1(i) > 0)
            Qmin = getQext1(i);
        if(getQabs1(i) < Qmin && getQabs1(i) > 0)
            Qmin = getQabs1(i);
        if(getQsca1(i) < Qmin && getQsca1(i) > 0)
            Qmin = getQsca1(i);
        if(getQext2(i) < Qmin && getQext2(i) > 0)
            Qmin = getQext2(i);
        if(getQabs2(i) < Qmin && getQabs2(i) > 0)
            Qmin = getQabs2(i);
        if(getQsca2(i) < Qmin && getQsca2(i) > 0)
            Qmin = getQsca2(i);

        if(getQext1(i) > Qmax)
            Qmax = getQext1(i);
        if(getQabs1(i) > Qmax)
            Qmax = getQabs1(i);
        if(getQsca1(i) > Qmax)
            Qmax = getQsca1(i);
        if(getQext2(i) > Qmax)
            Qmax = getQext2(i);
        if(getQabs2(i) > Qmax)
            Qmax = getQabs2(i);
        if(getQsca2(i) > Qmax)
            Qmax = getQsca2(i);
    }

    // Add a bit more space for good visualization
    Qmin *= 0.9;
    Qmax *= 1.10;

    // Add Gnuplot commands to file
    eff_writer << "reset" << endl;
    if(nr_of_wavelength > 1)
        eff_writer << "set log x" << endl;
    eff_writer << "set log y" << endl;
    eff_writer << "set grid" << endl;
    if(nr_of_wavelength > 1)
        eff_writer << "set xrange[" << wavelength_list[wavelength_offset] << ":" << wavelength_list.back()
                   << "]" << endl;
    eff_writer << "set yrange[" << Qmin << ":" << Qmax << "]" << endl;
    eff_writer << "set format x \"%.1te%02T\"" << endl;
    eff_writer << "set format y \"%.1te%02T\"" << endl;
    eff_writer << "set ylabel \'Q_{avg}\'" << endl;
    eff_writer << "set xlabel \'{/Symbol l} [m]\'" << endl;
    eff_writer << "set title \"" << gnu_title << "\"" << endl;
    eff_writer << "plot \'-\' with " << plot_sign << " title \'Q_{ext,x}\' lc rgb \"#0000F0\","
               << "\'-\' with " << plot_sign << " title \'Q_{ext,y}\' lc rgb \"#000090\","
               << "\'-\' with " << plot_sign << " title \'Q_{abs,x}\' lc rgb \"#FF0000\","
               << "\'-\' with " << plot_sign << " title \'Q_{abs,y}\' lc rgb \"#900000\","
               << "\'-\' with " << plot_sign << " title \'Q_{sca,x}\' lc rgb \"#FFFF00\","
               << "\'-\' with " << plot_sign << " title \'Q_{sca,y}\' lc rgb \"#909000\"" << endl;

    // Add Cext1 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getQext1(i) > 0)
            eff_writer << wavelength_list[i] << "\t" << getQext1(i) << endl;
    eff_writer << "e" << endl;

    // Add Cext2 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getQext2(i) > 0)
            eff_writer << wavelength_list[i] << "\t" << getQext2(i) << endl;
    eff_writer << "e" << endl;

    // Add Cabs1 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getQabs1(i) > 0)
            eff_writer << wavelength_list[i] << "\t" << getQabs1(i) << endl;
    eff_writer << "e" << endl;

    // Add Cabs2 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getQabs2(i) > 0)
            eff_writer << wavelength_list[i] << "\t" << getQabs2(i) << endl;
    eff_writer << "e" << endl;

    // Add Csca1 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getQsca1(i) > 0)
            eff_writer << wavelength_list[i] << "\t" << getQsca1(i) << endl;
    eff_writer << "e" << endl;

    // Add Csca2 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getQsca2(i) > 0)
            eff_writer << wavelength_list[i] << "\t" << getQsca2(i) << endl;
    eff_writer << "e" << endl;

    // Close text file writer
    eff_writer.close();

    // ------------------------------------------------------

    // Init text file writer for mass cross-sections
    ofstream kappa_writer(path_kappa.c_str());

    // Error message if the write does not work
    if(kappa_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n" << path_eff << endl;
        return false;
    }

    // Init plot limits
    double Kappa_min = 1e100, Kappa_max = -1e100;

    // Find min and max values over the wavelength (checks only Cext1 and Cabs1)
    // Find min and max values over the wavelength (checks only Cext1 and Cabs1)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
    {
        if(getKappaExt1(i) < Kappa_min && getKappaExt1(i) > 0)
            Kappa_min = getKappaExt1(i);
        if(getKappaAbs1(i) < Kappa_min && getKappaAbs1(i) > 0)
            Kappa_min = getKappaAbs1(i);
        if(getKappaSca1(i) < Kappa_min && getKappaSca1(i) > 0)
            Kappa_min = getKappaSca1(i);
        if(getKappaExt2(i) < Kappa_min && getKappaExt2(i) > 0)
            Kappa_min = getKappaExt2(i);
        if(getKappaAbs2(i) < Kappa_min && getKappaAbs2(i) > 0)
            Kappa_min = getKappaAbs2(i);
        if(getKappaSca2(i) < Kappa_min && getKappaSca2(i) > 0)
            Kappa_min = getKappaSca2(i);

        if(getKappaExt1(i) > Kappa_max)
            Kappa_max = getKappaExt1(i);
        if(getKappaAbs1(i) > Kappa_max)
            Kappa_max = getKappaAbs1(i);
        if(getKappaSca1(i) > Kappa_max)
            Kappa_max = getKappaSca1(i);
        if(getKappaExt2(i) > Kappa_max)
            Kappa_max = getKappaExt2(i);
        if(getKappaAbs2(i) > Kappa_max)
            Kappa_max = getKappaAbs2(i);
        if(getKappaSca2(i) > Kappa_max)
            Kappa_max = getKappaSca2(i);
    }

    // Add a bit more space for good visualization
    Kappa_min *= 0.9;
    Kappa_max *= 1.10;

    // Add Gnuplot commands to file
    kappa_writer << "reset" << endl;
    if(nr_of_wavelength > 1)
        kappa_writer << "set log x" << endl;
    kappa_writer << "set log y" << endl;
    kappa_writer << "set grid" << endl;
    if(nr_of_wavelength > 1)
        kappa_writer << "set xrange[" << wavelength_list[wavelength_offset] << ":" << wavelength_list.back()
                     << "]" << endl;
    kappa_writer << "set yrange[" << Kappa_min << ":" << Kappa_max << "]" << endl;
    kappa_writer << "set format x \"%.1te%02T\"" << endl;
    kappa_writer << "set format y \"%.1te%02T\"" << endl;
    kappa_writer << "set ylabel \'{/Symbol k}_{avg}  [m^2/kg]\'" << endl;
    kappa_writer << "set xlabel \'{/Symbol l} [m]\'" << endl;
    kappa_writer << "set title \"" << gnu_title << "\"" << endl;
    kappa_writer << "plot \'-\' with " << plot_sign << " title \'{/Symbol k}_{ext,x}\' lc rgb \"#0000F0\","
                 << "\'-\' with " << plot_sign << " title \'{/Symbol k}_{ext,y}\' lc rgb \"#000090\","
                 << "\'-\' with " << plot_sign << " title \'{/Symbol k}_{abs,x}\' lc rgb \"#FF0000\","
                 << "\'-\' with " << plot_sign << " title \'{/Symbol k}_{abs,y}\' lc rgb \"#900000\","
                 << "\'-\' with " << plot_sign << " title \'{/Symbol k}_{sca,x}\' lc rgb \"#FFFF00\","
                 << "\'-\' with " << plot_sign << " title \'{/Symbol k}_{sca,y}\' lc rgb \"#909000\"" << endl;

    // Add Cext1 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getKappaExt1(i) > 0)
            kappa_writer << wavelength_list[i] << "\t" << getKappaExt1(i) << endl;
    kappa_writer << "e" << endl;

    // Add Cext2 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getKappaExt2(i) > 0)
            kappa_writer << wavelength_list[i] << "\t" << getKappaExt2(i) << endl;
    kappa_writer << "e" << endl;

    // Add Cabs1 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getKappaAbs1(i) > 0)
            kappa_writer << wavelength_list[i] << "\t" << getKappaAbs1(i) << endl;
    kappa_writer << "e" << endl;

    // Add Cabs2 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getKappaAbs2(i) > 0)
            kappa_writer << wavelength_list[i] << "\t" << getKappaAbs2(i) << endl;
    kappa_writer << "e" << endl;

    // Add Csca1 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getKappaSca1(i) > 0)
            kappa_writer << wavelength_list[i] << "\t" << getKappaSca1(i) << endl;
    kappa_writer << "e" << endl;

    // Add Csca2 data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getKappaSca2(i) > 0)
            kappa_writer << wavelength_list[i] << "\t" << getKappaSca2(i) << endl;
    kappa_writer << "e" << endl;

    // Close text file writer
    kappa_writer.close();

    // ------------------------------------------------------

    // Init text file writer for cross-section differences
    ofstream diff_writer(path_diff.c_str());

    // Error message if the write does not work
    if(diff_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n" << path_cross << endl;
        return false;
    }

    // Init plot limits
    Cmin = 1e100, Cmax = -1e100;

    // Find min and max values over the wavelength
    // Checks only differences between Cext and Cabs as well as Ccirc
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
    {
        if(abs(getCext1(i) - getCext2(i)) < Cmin && abs(getCext1(i) - getCext2(i)) > 0)
            Cmin = abs(getCext1(i) - getCext2(i));
        if(abs(getCabs1(i) - getCabs2(i)) < Cmin && abs(getCabs1(i) - getCabs2(i)) > 0)
            Cmin = abs(getCabs1(i) - getCabs2(i));
        if(abs(getCcirc(i)) < Cmin && abs(getCcirc(i)) > 0)
            Cmin = abs(getCcirc(i));

        if(abs(getCext1(i) - getCext2(i)) > Cmax)
            Cmax = abs(getCext1(i) - getCext2(i));
        if(abs(getCabs1(i) - getCabs2(i)) > Cmax)
            Cmax = abs(getCabs1(i) - getCabs2(i));
        if(abs(getCcirc(i)) > Cmax)
            Cmax = abs(getCcirc(i));
    }

    // Add a bit more space for good visualization
    Cmin *= 0.9;
    Cmax *= 1.10;

    // Add Gnuplot commands to file
    diff_writer << "reset" << endl;
    if(nr_of_wavelength > 1)
        diff_writer << "set log x" << endl;
    diff_writer << "set log y" << endl;
    diff_writer << "set grid" << endl;
    if(nr_of_wavelength > 1)
        diff_writer << "set xrange[" << wavelength_list[wavelength_offset] << ":" << wavelength_list.back()
                    << "]" << endl;
    diff_writer << "set yrange[" << Cmin << ":" << Cmax << "]" << endl;
    diff_writer << "set format x \"%.1te%02T\"" << endl;
    diff_writer << "set format y \"%.1te%02T\"" << endl;
    diff_writer << "set ylabel \'C_{avg} [m^{2}]\'" << endl;
    diff_writer << "set xlabel \'{/Symbol l} [m]\'" << endl;
    diff_writer << "set title \"" << gnu_title << "\"" << endl;
    diff_writer << "plot \'-\' with " << plot_sign << " title \'|dC_{ext}| (Cpol)\' lc rgb \"#0000FF\","
                << "\'-\' with " << plot_sign << " title \'|dC_{abd}|\'  lc rgb \"#FF0000\","
                << "\'-\' with " << plot_sign << " title \'|dC_{phas}| (C_{circ})\' lc rgb \"#800080\""
                << endl;

    // Add Cext difference to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(abs(getCext1(i) - getCext2(i)) > 0)
            diff_writer << wavelength_list[i] << "\t" << abs(getCext1(i) - getCext2(i)) << endl;
    diff_writer << "e" << endl;

    // Add Cabs difference to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(abs(getCabs1(i) - getCabs2(i)) > 0)
            diff_writer << wavelength_list[i] << "\t" << abs(getCabs1(i) - getCabs2(i)) << endl;
    diff_writer << "e" << endl;

    // Add Ccirc data to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(abs(getCcirc(i)) > 0)
            diff_writer << wavelength_list[i] << "\t" << abs(getCcirc(i)) << endl;
    diff_writer << "e" << endl;

    // Close text file writer
    diff_writer.close();

    // ------------------------------------------------------

    // Init text file writer for dust grain parameters
    ofstream data_writer(path_data.c_str());

    // Error message if the write does not work
    if(data_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n" << path_data << endl;
        return false;
    }

    // Add Gnuplot commands to file
    data_writer << "#material" << endl;
    data_writer << str_title << endl;
    data_writer << "#weight" << endl;
    data_writer << getWeight() << endl;
    data_writer << "#aspect ratio" << endl;
    data_writer << aspect_ratio << endl;
    data_writer << "#sublimation temp [K]" << endl;
    data_writer << sub_temp << endl;
    if(!is_mixture)
    {
        data_writer << "#size distribution keyword" << endl;
        data_writer << size_keyword << endl;
        data_writer << "#size distribution parameters [1-7]" << endl;
        for(uint i = 0; i < NR_OF_SIZE_DIST_PARAM; i++)
            data_writer << size_parameter[i] << endl;
    }
    data_writer << "#can align" << endl;
    data_writer << is_align << endl;
    data_writer << "#density [kg]" << endl;
    data_writer << material_density << endl;
    data_writer << "#gold g factor" << endl;
    data_writer << gold_g_factor << endl;
    data_writer << "#f_highJ" << endl;
    data_writer << f_highJ << endl;
    data_writer << "#f_correlation" << endl;
    data_writer << f_cor << endl;
    data_writer << "#min. ID\tmax. ID" << endl;
    data_writer << a_min_global + 1 << "\t" << a_max_global + 1 << endl;
    data_writer << "#min. radius\tmax. radius [m]" << endl;
    data_writer << getSizeMin() << "\t" << getSizeMax() << endl;

    // Add data header to file
    data_writer << "#wavelength\tavgCext1\tavgCext2\tavgCabs1\tavgCabs2\tavgCsca1\tavgCsc"
                   "a2\tCcirc\tavgHGg"
                << endl;

    // Add cross-sections to file
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        data_writer << wavelength_list[i] << "\t" << getCext1(i) << "\t" << getCext2(i) << "\t" << getCabs1(i)
                    << "\t" << getCabs2(i) << "\t" << getCsca1(i) << "\t" << getCsca2(i) << "\t"
                    << getCcirc(i) << "\t" << getHGg(i) << endl;

    // Close text file writer
    data_writer.close();

    // ------------------------------------------------------

    // Init text file writer for Henyey-Greenstein g factor
    ofstream g_writer(path_g.c_str());

    // Error message if the write does not work
    if(g_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n" << path_g << endl;
        return false;
    }

    // Init plot limits
    Cmin = 1e100, Cmax = -1e100;

    // Find min and max values over the wavelength
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
    {
        if(abs(getHGg(i)) < Cmin && abs(getHGg(i)) > 0)
            Cmin = abs(getHGg(i));
        if(abs(getHGg(i)) > Cmax)
            Cmax = abs(getHGg(i));
    }

    // Add a bit more space for good visualization
    Cmin *= 0.9;
    Cmax *= 1.10;

    // Add Gnuplot commands to file
    g_writer << "reset" << endl;
    if(nr_of_wavelength > 1)
        g_writer << "set log x" << endl;
    g_writer << "set log y" << endl;
    g_writer << "set grid" << endl;
    g_writer << "unset key" << endl;
    if(nr_of_wavelength > 1)
        g_writer << "set xrange[" << wavelength_list[wavelength_offset] << ":" << wavelength_list.back()
                 << "]" << endl;
    g_writer << "set yrange[" << Cmin << ":" << Cmax << "]" << endl;
    g_writer << "set format x \"%.1te%02T\"" << endl;
    g_writer << "set format y \"%.1te%02T\"" << endl;
    g_writer << "set ylabel \'|g_{avg}|'" << endl;
    g_writer << "set xlabel \'{/Symbol l} [m]\'" << endl;
    g_writer << "set title \"" << gnu_title << "\"" << endl;
    g_writer << "plot \'-\' with " << plot_sign << " lc rgb \"#0000F0\"" << endl;

    // Add HG g factor to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getHGg(i) > 0)
            g_writer << wavelength_list[i] << "\t" << getHGg(i) << endl;
    g_writer << "e" << endl;

    // Close text file writer
    g_writer.close();

    // ------------------------------------------------------

    if(phID == PH_MIE && is_mixture)
    {
        // Init text file writer for scattering matrix
        ofstream scat_writer(path_scat.c_str());

        // Error message if the write does not work
        if(scat_writer.fail())
        {
            cout << "\nERROR: Cannot write to " << path_scat << endl;
            return false;
        }

        // Init plot limits
        double S11min = 1e100, S11max = -1e100;
        double S12min = 1e100, S12max = -1e100;

        // Init pointer arrays for scattering matrix elements
        double **S11, **S12;

        // Get weight from grain size distribution
        double weight = getWeight();

        // Init pointer dimension
        S11 = new double *[nr_of_wavelength - wavelength_offset];
        S12 = new double *[nr_of_wavelength - wavelength_offset];
        for(uint w = wavelength_offset; w < nr_of_wavelength; w++)
        {
            uint wID = w - wavelength_offset;
            S11[wID] = new double[nr_of_scat_theta];
            S12[wID] = new double[nr_of_scat_theta];
            for(uint sth = 0; sth < nr_of_scat_theta; sth++)
            {
                // Init and reset variables
                double sum = 0;

                double * S11_tmp = new double[nr_of_dust_species];
                double * S12_tmp = new double[nr_of_dust_species];
                for(uint a = 0; a < nr_of_dust_species; a++)
                {
                    if(sizeIndexUsed(a))
                    {
                        double Csca_tmp = getCscaMean(a, w);
                        sum += Csca_tmp;
                        double rel_amount = a_eff_3_5[a] / weight;
                        S11_tmp[a] = Csca_tmp * rel_amount * // getScatteredFractionMie(a, w, sth) *
                                     getScatteringMatrixElement(a, w, 0, 0, sth, 0, 0);
                        S12_tmp[a] = Csca_tmp * rel_amount * // getScatteredFractionMie(a, w, sth) *
                                     getScatteringMatrixElement(a, w, 0, 0, sth, 0, 1);
                    }
                    else
                    {
                        S11_tmp[a] = 0;
                        S12_tmp[a] = 0;
                    }
                }
                S11[wID][sth] = 1 / sum *
                                CMathFunctions::integ_dust_size(
                                    a_eff, S11_tmp, nr_of_dust_species, a_min_global, a_max_global);
                S12[wID][sth] = 1 / sum *
                                CMathFunctions::integ_dust_size(
                                    a_eff, S12_tmp, nr_of_dust_species, a_min_global, a_max_global);

                if(S11[wID][sth] < S11min)
                    S11min = S11[wID][sth];
                if(S12[wID][sth] < S12min)
                    S12min = S12[wID][sth];

                if(S11[wID][sth] > S11max)
                    S11max = S11[wID][sth];
                if(S12[wID][sth] > S12max)
                    S12max = S12[wID][sth];
            }
        }

        // Add Gnuplot commands to file
        scat_writer << "reset" << endl;
        scat_writer << "set grid" << endl;
        scat_writer << "set multiplot layout 2,1 rowsfirst" << endl;

        if(nr_of_wavelength > 1)
            scat_writer << "set xrange[" << 0 << ":" << nr_of_scat_theta << "]" << endl;
        scat_writer << "set yrange[" << S11min << ":" << S11max << "]" << endl;
        scat_writer << "set format x \"%.1f\"" << endl;
        scat_writer << "set format y \"%.1te%02T\"" << endl;
        scat_writer << "set ylabel \'S11\'" << endl;
        scat_writer << "set xlabel \'{\u03B8} [°]\'" << endl;
        scat_writer << "set title \"" << gnu_title << "\"" << endl;
        scat_writer << "plot ";
        for(uint w = wavelength_offset; w < nr_of_wavelength; w++)
        {
            scat_writer << "\'-\' with " << plot_sign << " title \'" << wavelength_list[w] << " [m]\'";
            if(w != nr_of_wavelength - 1)
                scat_writer << ",";
            else
                scat_writer << endl;
        }

        for(uint w = wavelength_offset; w < nr_of_wavelength; w++)
        {
            uint wID = w - wavelength_offset;

            for(uint sth = 0; sth < nr_of_scat_theta; sth++)
                scat_writer << sth << "\t" << S11[wID][sth] << endl;
            scat_writer << "e" << endl;
        }

        if(nr_of_wavelength > 1)
            scat_writer << "set xrange[" << 0 << ":" << nr_of_scat_theta << "]" << endl;
        scat_writer << "set yrange[" << S12min << ":" << S12max << "]" << endl;
        scat_writer << "set format x \"%.1f\"" << endl;
        scat_writer << "set format y \"%.1te%02T\"" << endl;
        scat_writer << "set ylabel \'S12\'" << endl;
        scat_writer << "set xlabel \'{\u03B8} [°]\'" << endl;
        scat_writer << "set title \"" << gnu_title << "\"" << endl;
        scat_writer << "plot ";
        for(uint w = wavelength_offset; w < nr_of_wavelength; w++)
        {
            scat_writer << "\'-\' with " << plot_sign << " title \'" << wavelength_list[w] << " [m]\'";
            if(w != nr_of_wavelength - 1)
                scat_writer << ",";
            else
                scat_writer << endl;
        }

        for(uint w = wavelength_offset; w < nr_of_wavelength; w++)
        {
            uint wID = w - wavelength_offset;

            for(uint sth = 0; sth < nr_of_scat_theta; sth++)
                scat_writer << sth << "\t" << S12[wID][sth] << endl;
            scat_writer << "e" << endl;
        }

        for(uint w = wavelength_offset; w < nr_of_wavelength; w++)
        {
            uint wID = w - wavelength_offset;
            delete[] S11[wID];
            delete[] S12[wID];
        }
        delete[] S11;
        delete[] S12;

        // Close text file writer
        scat_writer.close();
    }

    return true;
}

void CDustComponent::preCalcEffProperties(parameters & param)
{
    // -------------- Calculate average-mass of effective grain size --------------
    avg_mass = getAvgMass();

    // -------------- Calculate cross-sections of effective grain size --------------
    tCext1 = new double[nr_of_wavelength];
    tCext2 = new double[nr_of_wavelength];
    tCabs1 = new double[nr_of_wavelength];
    tCabs2 = new double[nr_of_wavelength];
    tCsca1 = new double[nr_of_wavelength];
    tCsca2 = new double[nr_of_wavelength];
    tCcirc = new double[nr_of_wavelength];
    tHGg = new double[nr_of_wavelength];

    for(uint w = 0; w < nr_of_wavelength; w++)
    {
        tCext1[w] = getCext1(w);
        tCext2[w] = getCext2(w);
        tCabs1[w] = getCabs1(w);
        tCabs2[w] = getCabs2(w);
        tCsca1[w] = getCsca1(w);
        tCsca2[w] = getCsca2(w);
        tCcirc[w] = getCcirc(w);
        tHGg[w] = getHGg(w);
    }

    // -------------- Calculate emission of effective grain size --------------
    if(param.isTemperatureSimulation())
    {
        // Get number of temperatures from tab_temp spline
        uint nr_of_temperatures = tab_temp.size();

        // Init spline for absorption/emission energy interpolation
        tab_em_eff.resize(nr_of_temperatures);

        for(uint t = 0; t < nr_of_temperatures; t++)
        {
            // Get temperature from tab_temp spline
            double tmp_temp = tab_temp.getValue(t);

            // Init a temporary array for QB values
            double * tmpQB = new double[nr_of_wavelength];

            // Calculate absorption cross-section times Planck function for each
            // wavelength
            for(uint w = 0; w < nr_of_wavelength; w++)
            {
                // Calculate mean absorption cross-section
                double meanCabs = (2.0 * tCabs1[w] + tCabs2[w]) / 3.0;

                // Set absroption/emission energy at current wavelength and temperature
                tmpQB[w] = meanCabs * tab_planck[w].getValue(t);
            }

            // Calculate QB integrated over all wavelengths
            tab_em_eff.setValue(
                t, CMathFunctions::integ(wavelength_list, tmpQB, 0, nr_of_wavelength - 1), tmp_temp);

            // Delete pointer array
            delete[] tmpQB;
        }

        // Create spline for interpolation
        tab_em_eff.createSpline();
    }

    // -------------- Calculate probability lists --------------
    // Init pointer array of prob_lists
    dust_prob = new prob_list[nr_of_wavelength];
    abs_prob = new prob_list[nr_of_wavelength];
    sca_prob = new prob_list[nr_of_wavelength];

    for(uint w = 0; w < nr_of_wavelength; w++)
    {
        // Init pointer array for integration
        double * amount = new double[nr_of_dust_species];
        double * amount_abs = new double[nr_of_dust_species];
        double * amount_sca = new double[nr_of_dust_species];

        // Add relative amount of dust grains in each size bin
        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            if(sizeIndexUsed(a))
            {
                amount[a] = a_eff_3_5[a];
                amount_abs[a] = a_eff_3_5[a] * getCabsMean(a, w);
                amount_sca[a] = a_eff_3_5[a] * getCscaMean(a, w);
            }
            else
            {
                amount[a] = 0;
                amount_abs[a] = 0;
                amount_sca[a] = 0;
            }
        }

        // Init propability list to pick random grain size
        dust_prob[w].resize(nr_of_dust_species);
        abs_prob[w].resize(nr_of_dust_species);
        sca_prob[w].resize(nr_of_dust_species);

        // Set each grain size bin with the relative amount
        CMathFunctions::probListInteg(a_eff, amount, dust_prob[w]);
        CMathFunctions::probListInteg(a_eff, amount_abs, abs_prob[w]);
        CMathFunctions::probListInteg(a_eff, amount_sca, sca_prob[w]);

        delete[] amount;
        delete[] amount_abs;
        delete[] amount_sca;
    }
}

void CDustComponent::preCalcAbsorptionRates()
{
    // Get number of temperatures from tab_temp spline
    uint const nr_of_temperatures = tab_temp.size();

    // Create tabulated absorption/emission rates spline
    tab_em = new spline[nr_of_dust_species];
    if(calorimetry_loaded)
        tab_em_inv = new spline[nr_of_dust_species];

    // Init counter for progress
    ullong per_counter = 0;

#pragma omp parallel for
    for(int a = 0; a < int(nr_of_dust_species); a++)
    {
        // Resize tabulated absorption/emission rates spline
        tab_em[a].resize(nr_of_temperatures);
        if(calorimetry_loaded)
            tab_em_inv[a].resize(nr_of_temperatures);

        // Init a temporary array for QB values
        double * tmpQB = new double[nr_of_wavelength];

        // Set each entry of tab_em with integrated value of the Planck function times
        // the absorption cross-section
        for(uint t = 0; t < nr_of_temperatures; t++)
        {
            // Get temperature from tab_temp spline
            double tmp_temp = tab_temp.getValue(t);

            // Calculate absorption cross-section times Planck function for each
            // wavelength
            for(uint w = 0; w < WL_STEPS; w++)
                tmpQB[w] = getCabsMean(a, w) * tab_planck[w].getValue(t);

            // Calculate QB integrated over all wavelengths
            double tt = CMathFunctions::integ(wavelength_list, tmpQB, 0, WL_STEPS - 1);

            // Set tab_em spline with the integrated value
            tab_em[a].setValue(t, tt, tmp_temp);
            if(calorimetry_loaded)
                tab_em_inv[a].setValue(t, tmp_temp, tt);
        }

        // Create spline for absorption/emission rates
        tab_em[a].createSpline();
        if(calorimetry_loaded)
            tab_em_inv[a].createSpline();

#pragma omp critical
        {
            // Show progress
            if(per_counter % 2 == 0)
            {
                printIDs();
                cout << "- pre-calculation of absorption rates: "
                     << 100.0 * float(per_counter) / float(nr_of_dust_species - 1)
                     << " [%]                         \r";
            }
        }

        // Increase progress counter
        per_counter++;

        // Delete pointer array
        delete[] tmpQB;
    }
}

void CDustComponent::preCalcMieScatteringProb()
{
    // Init arrays of splines/interp
    avg_scattering_frac = new spline *[nr_of_dust_species];
    phase_pdf = new interp *[nr_of_dust_species];

    // Calculate differences between two angles
    double d_ang = PI / double(nr_of_scat_theta - 1);

    // Init counter for progress
    ullong per_counter = 0;

#pragma omp parallel for
    for(int a = 0; a < int(nr_of_dust_species); a++)
    {
        // Init pointer arrays
        double * S11_tmp = new double[nr_of_scat_theta];
        double * S11_solid_angle = new double[nr_of_scat_theta];
        double * tmp_scat_frac = new double[nr_of_scat_theta];

        // Init arrays of interp
        avg_scattering_frac[a] = new spline[nr_of_wavelength];
        phase_pdf[a] = new interp[nr_of_wavelength];
        for(uint w = 0; w < nr_of_wavelength; w++)
        {
            // Resize splines for scattering angles
            avg_scattering_frac[a][w].resize(nr_of_scat_theta);
            phase_pdf[a][w].resize(nr_of_scat_theta);

            for(uint sth = 0; sth < nr_of_scat_theta; sth++)
            {
                // Calculate the scattering propability in a certain direction
                S11_tmp[sth] = double(sca_mat[a][w][0][0][sth](0, 0));

                // Calculate the modified angle for integration
                S11_solid_angle[sth] = PIx2 * cos(double(sth) * d_ang);

                // Set the current value for the integration array tmp_scat_frac
                if(sth == 0)
                    tmp_scat_frac[sth] = 0;
                else
                    tmp_scat_frac[sth] = -CMathFunctions::integ(S11_solid_angle, S11_tmp, 0, sth);
            }

            // Integral of the scattering S11 value over the full sphere
            double int_scat_frac = tmp_scat_frac[nr_of_scat_theta - 1];

            for(uint sth = 0; sth < nr_of_scat_theta; sth++)
            {
                if(int_scat_frac > 0)
                {
                    // Set the cumulative distribution function of being scattered at a
                    // certain angle
                    avg_scattering_frac[a][w].setValue(
                        sth, tmp_scat_frac[sth] / int_scat_frac, double(sth) * d_ang);

                    // Set the phase function of how much is scattered at a certain angle
                    phase_pdf[a][w].setValue(sth, double(sth) * d_ang, S11_tmp[sth] / int_scat_frac);
                }
            }

            // Activate spline for the cumulative distribution function
            avg_scattering_frac[a][w].createSpline();
        }

        delete[] S11_tmp;
        delete[] S11_solid_angle;
        delete[] tmp_scat_frac;

#pragma omp critical
        {
            // Show progress
            if(per_counter % 2 == 0)
            {
                printIDs();
                cout << "- pre-calculation of Mie probabilities: "
                     << 100.0 * float(per_counter) / float(nr_of_dust_species - 1)
                     << " [%]                         \r";
            }
        }

        // Increase progress counter
        per_counter++;
    }
}

void CDustComponent::preCalcWaveProb()
{
    // Get number of temperatures from tab_temp spline
    uint nr_of_temperatures = tab_temp.size();

    // Init spline for tabulated average Planck function fractions
    avg_planck_frac = new prob_list[nr_of_dust_species * nr_of_temperatures];

    // Init counter for progress
    ullong per_counter = 0;

#pragma omp parallel for
    // Set each entry of avg_planck_frac with integrated value of pl_mean
    for(int t = 0; t < int(nr_of_temperatures); t++)
    {
        // Init array of mean absorption cross-section times dPlanck/dT for each
        // wavelength and grain size
        double ** pl_mean;
        pl_mean = new double *[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            pl_mean[a] = new double[WL_STEPS];

        // Get temperature from tab_temp spline
        double temp = tab_temp.getValue(uint(t));

        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            // Resize avg_planck_frac to number of wavelengths
            avg_planck_frac[t * nr_of_dust_species + a].resize(WL_STEPS);
        }

        // Calculate dPlanck/dT times mean absorption cross-section for each wavelength
        for(uint w = 0; w < WL_STEPS; w++)
        {
            // Calculate dPlanck/dT
            double t_pl = CMathFunctions::dplanck_dT(wavelength_list[w], temp);

            for(uint a = 0; a < nr_of_dust_species; a++)
            {
                // Calculate dPlanck/dT times mean absorption cross-section
                pl_mean[a][w] = getCabsMean(a, w) * t_pl;
            }
        }

        for(uint a = 0; a < nr_of_dust_species; a++)
            CMathFunctions::probListInteg(
                wavelength_list, pl_mean[a], avg_planck_frac[t * nr_of_dust_species + a]);

        // Delete pointer array
        for(uint a = 0; a < nr_of_dust_species; a++)
            delete[] pl_mean[a];
        delete[] pl_mean;

#pragma omp critical
        {
            // Show progress
            if(per_counter % 10 == 0)
            {
                printIDs();
                cout << "- precalculation of wavelength-probabilities: "
                     << 100.0 * float(per_counter) / float(nr_of_temperatures - 1) << " [%]           \r";
            }
        }

        // Increase progress counter
        per_counter++;
    }
}

void CDustComponent::preCalcTemperatureLists(double minTemp, double maxTemp, uint nr_of_temperatures)
{
    // Init variables
    double tmp_temp, tmp_lambda;

    // Init spline for tabulated Planck function values
    tab_planck = new spline[nr_of_wavelength];

    // Resize tabulated temperature spline
    tab_temp.resize(nr_of_temperatures);

    // Set each entry of tab_temp with the corresponding temperature (exponential curve)
    for(uint t = 0; t < nr_of_temperatures; t++)
    {
        // Calculate the temperature of a certain index
        tmp_temp = minTemp * pow(double(maxTemp) / minTemp, double(t) / (nr_of_temperatures - 1));

        // Add the temperature to the spline
        tab_temp.setValue(t, double(t), tmp_temp);
    }

    // Create spline for tabulated temperatures
    tab_temp.createSpline();

    // Set each entry of tab_planck with a spline for each wavelength
    for(uint w = 0; w < nr_of_wavelength; w++)
    {
        // Resize tabulated Planck spline
        tab_planck[w].resize(nr_of_temperatures);

        // Show progress
        printIDs();
        cout << "- pre-calculation of Planck functions: " << 100.0 * float(w) / float(nr_of_wavelength)
             << "                                \r";

        // Set each entry of tab_planck with the corresponding Planck function values
        // that depend on the temperature
        for(uint t = 0; t < nr_of_temperatures; t++)
        {
            // Get temperature from tab_temp spline
            double temp = tab_temp.getValue(t);

            // Calculate Planck function value for temperature and wavelength
            double pl = CMathFunctions::planck(wavelength_list[w], temp);

            // Set tab_planck spline with the Planck function value
            tab_planck[w].setValue(t, temp, pl);
        }

        // Create spline for tabulated Planck function values
        tab_planck[w].createSpline();
    }
}

bool CDustComponent::calcSizeDistribution(dlist values, double * mass)
{
    // Calculates various grain size distributions.
    // From DustEM Code "https://www.ias.u-psud.fr/DUSTEM/"
    // Caution by using dustem parameters without conversion (its in cgs and POLARIS needs
    // the parameters in SI)

    // Set a_eff and a_eff^2 from input
    for(int a = 0; a < nr_of_dust_species; a++)
    {
        a_eff[a] = values[a];
        a_eff_2[a] = a_eff[a] * a_eff[a];
    }

    // Calculate the dust grain size distribution depending on user input
    if(size_keyword.find("plaw") != std::string::npos)
    {
        // Power-law distribution
        for(int a = 0; a < nr_of_dust_species; a++)
            a_eff_3_5[a] = pow(a_eff[a], size_parameter[0]);

        // Add exponential decay if demanded
        if(size_keyword.find("-ed") != std::string::npos)
        {
            for(int a = 0; a < nr_of_dust_species; a++)
                if(a_eff[a] > size_parameter[1])
                    a_eff_3_5[a] *=
                        exp(-pow((a_eff[a] - size_parameter[1]) / size_parameter[2], size_parameter[3]));
        }
        // Add curvature term if demanded
        if(size_keyword.find("-cv") != std::string::npos)
        {
            double au, zeta, zxp, gama;
            if(size_keyword.find("-ed") != std::string::npos)
            {
                au = size_parameter[4];
                zeta = abs(size_parameter[5]);
                zxp = sign(size_parameter[5]);
                gama = size_parameter[6];
            }
            else
            {
                au = size_parameter[1];
                zeta = abs(size_parameter[2]);
                zxp = sign(size_parameter[2]);
                gama = size_parameter[3];
            }
            for(int a = 0; a < nr_of_dust_species; a++)
                a_eff_3_5[a] *= pow(1.0 + zeta * pow(a_eff[a] / au, gama), zxp);
        }
    }
    else if(size_keyword.find("logn") != std::string::npos)
    {
        // Log-normal distribution
        if(size_parameter[0] == 0 || size_parameter[1] == 0)
        {
            cout << "\nERROR: Centroid or sigma of log-normal cannot be 0!" << endl;
            return false;
        }
        for(int a = 0; a < nr_of_dust_species; a++)
        {
            double aux = log(a_eff[a]);
            double argu = -0.5 * pow((aux - log(size_parameter[0])) / size_parameter[1], 2);
            a_eff_3_5[a] = exp(argu);
        }
    }
    else if(size_keyword.find("zda") != std::string::npos)
    {
        // ZDA distribution from Camps/Trust benchmark
        for(int a = 0; a < nr_of_dust_species; a++)
        {
            double a_eff_micron = a_eff[a] * 1e6;
            double log_g = size_parameter[0] + size_parameter[1] * log10(a_eff_micron);
            if(size_parameter[2] != 0 || size_parameter[3] != 0 || size_parameter[4] != 0)
                log_g -=
                    size_parameter[3] * pow(abs(log10(a_eff_micron / size_parameter[4])), size_parameter[2]);
            if(size_parameter[5] != 0 || size_parameter[6] != 0 || size_parameter[7] != 0)
                log_g -=
                    size_parameter[6] * pow(abs(log10(a_eff_micron / size_parameter[7])), size_parameter[5]);
            if(size_parameter[8] != 0 || size_parameter[9] != 0 || size_parameter[10] != 0)
                log_g -= size_parameter[9] * pow(abs(a_eff_micron - size_parameter[10]), size_parameter[8]);
            if(size_parameter[11] != 0 || size_parameter[12] != 0 || size_parameter[13] != 0)
                log_g -= size_parameter[12] * pow(abs(a_eff_micron - size_parameter[13]), size_parameter[11]);
            a_eff_3_5[a] = pow(10, log_g);
        }
    }

    // Set size distribution times a_eff^2 and mass of grains a certain size
    for(int a = 0; a < nr_of_dust_species; a++)
    {
        // Set relative abundance of dust grains times their squared radius
        a_eff_1_5[a] = a_eff_3_5[a] * a_eff_2[a];

        // For PAHs, their mass can be calculated as standard grains or specifically for
        // PAHs
        if(size_keyword.find("pah") != std::string::npos)
            mass[a] = getPahMass(a);
        else
            mass[a] = 4.0 / 3.0 * PI * a_eff[a] * a_eff[a] * a_eff[a] * material_density;
    }
    return true;
}

bool CDustComponent::add(double ** size_fraction, CDustComponent * comp)
{
    // Get global min and max grain sizes
    double a_min = comp->getSizeMin();
    double a_max = comp->getSizeMax();

    // Use the highes a_max and lowest a_min of the components
    if(a_min_global > a_min)
        a_min_global = a_min;
    if(a_max_global < a_max)
        a_max_global = a_max;

    // The first component to add to the mixture initializes the variables
    if(comp->getComponentId() == 0)
    {
        // set that the current dust component was made via mixing
        is_mixture = true;

        // Get common parameters from all dust components
        nr_of_dust_species = comp->getNrOfDustSpecies();
        nr_of_incident_angles = comp->getNrOfIncidentAngles();
        nr_of_scat_theta = comp->getNrOfScatTheta();
        nr_of_scat_phi = comp->getNrOfScatPhi();
        nr_of_scat_mat_elements = comp->getNrOfScatMatElements();
        nr_of_calorimetry_temperatures = comp->getNrOfCalorimetryTemperatures();
        f_cor = comp->getCorrelationFactor();
        f_highJ = comp->getFHighJ();

        // Init dust properties to be filled with grain properties
        initDustProperties();

        // Set grain size grid, which needs to be done only once
        for(int a = 0; a < nr_of_dust_species; a++)
        {
            a_eff[a] = comp->getEffectiveRadius(a);
            a_eff_2[a] = comp->getEffectiveRadius_2(a);
        }

        // Resize Qtrq and Henyey-Greenstein g factor
        for(uint i = 0; i < nr_of_dust_species * nr_of_wavelength; i++)
        {
            Qtrq[i].resize(nr_of_incident_angles);
            HG_g_factor[i].resize(nr_of_incident_angles);
        }

        // If the scattering matrix is read in, add them together as well
        if(comp->getScatLoaded())
            initScatteringMatrixArray();

        // Set the phase function initially
        phID = comp->getPhaseFunctionID();

        // If the colarimetry data is read in, add them together as well
        if(comp->getCalorimetryLoaded())
        {
            // Init calorimetry data
            initCalorimetry();

            // Set calorimetric temperatures, which needs to be done only once
            for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
                calorimetry_temperatures[t] = comp->getCalorimetricTemperature(t);
        }
    }

    for(int a = 0; a < nr_of_dust_species; a++)
    {
        // Check if the dust grain size is really the same
        if(a_eff[a] != comp->getEffectiveRadius(a))
        {
            cout << "\nERROR: Wrong grain size at position " << a + 1 << " ( " << comp->getEffectiveRadius(a)
                 << " )" << endl;
            return false;
        }

        if(comp->sizeIndexUsed(a, a_min, a_max))
        {
            // Mix size distribution with their relative fraction
            a_eff_3_5[a] += size_fraction[a][0];
            a_eff_1_5[a] += size_fraction[a][0] * comp->getEffectiveRadius_2(a);
            mass[a] += size_fraction[a][1] * comp->getMass(a);
        }
    }

    if(comp->getCalorimetryLoaded())
        for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
        {
            // Check if the calorimetry data is really the same
            if(calorimetry_temperatures[t] != comp->getCalorimetricTemperature(t))
            {
                cout << "\nERROR: Wrong calorimetric temperature at position " << t + 1 << " ( "
                     << comp->getCalorimetricTemperature(t) << " )" << endl;
                return false;
            }
        }

    if(comp->getScatLoaded())
    {
        // Mix scattering matrix for mixture
        for(int a = 0; a < nr_of_dust_species; a++)
            for(uint w = 0; w < nr_of_wavelength; w++)
                for(uint inc = 0; inc < nr_of_incident_angles; inc++)
                    for(uint sph = 0; sph < nr_of_scat_phi; sph++)
                        for(uint sth = 0; sth < nr_of_scat_theta; sth++)
                            for(uint i_mat = 0; i_mat < nr_of_scat_mat_elements; i_mat++)
                                for(uint j_mat = 0; j_mat < nr_of_scat_mat_elements; j_mat++)
                                    sca_mat[a][w][inc][sph][sth](i_mat, j_mat) +=
                                        size_fraction[a][1] *
                                        comp->getScatteringMatrixElement(a, w, inc, sph, sth, i_mat, j_mat);
    }

    if(comp->getCalorimetryLoaded())
    {
        // Mix enthalpy for mixture
        for(int a = 0; a < nr_of_dust_species; a++)
            for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
                enthalpy[a][t] += size_fraction[a][1] * comp->getEnthalpy(a, t);
    }

    // Show progress
    printIDs();
    cout << "- mixing average cross sections\r";

    // Mix optical properties of the dust grains
    for(uint w = 0; w < nr_of_wavelength; w++)
    {
        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            if(comp->sizeIndexUsed(a, a_min, a_max))
            {
                // Add optical properties on top of the mixture ones
                addQext1(a, w, size_fraction[a][1] * comp->getQext1(a, w));
                addQext2(a, w, size_fraction[a][1] * comp->getQext2(a, w));
                addQabs1(a, w, size_fraction[a][1] * comp->getQabs1(a, w));
                addQabs2(a, w, size_fraction[a][1] * comp->getQabs2(a, w));
                addQsca1(a, w, size_fraction[a][1] * comp->getQsca1(a, w));
                addQsca2(a, w, size_fraction[a][1] * comp->getQsca2(a, w));
                addQcirc(a, w, size_fraction[a][1] * comp->getQcirc(a, w));
                addHGg(a, w, size_fraction[a][1] * comp->getHGg(a, w));
            }
        }
    }

    if(comp->isAligned())
    {
        // Show progress
        printIDs();
        cout << "- mixing Qtrq and HG g                         \r";

        // Init variables
        double tmpHGgX, tmpQtrqX;
        double tmpHGgY, tmpQtrqY;

        // Calculate the difference between two incident angles
        double d_ang;
        if(nr_of_incident_angles > 1)
            d_ang = PI / double(nr_of_incident_angles - 1);
        else
            d_ang = 1;

        for(uint w = 0; w < nr_of_wavelength; w++)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
            {
                uint i = w * nr_of_dust_species + a;
                for(uint i_inc = 0; i_inc < nr_of_incident_angles; i_inc++)
                {
                    // Get incident angle and value of Qtrq and Henyey-Greenstein g factor
                    comp->getQtrq(i, i_inc, tmpQtrqX, tmpQtrqY);
                    comp->getHG_g_factor(i, i_inc, tmpHGgX, tmpHGgY);

                    // Add the values on top of the mixture Qtrq and Henyey-Greenstein g
                    // factor
                    Qtrq[i].addValue(i_inc, i_inc * d_ang, size_fraction[a][1] * tmpQtrqY);
                    HG_g_factor[i].addValue(i_inc, i_inc * d_ang, size_fraction[a][1] * tmpHGgY);
                }
            }
        }
    }

    // Mix various parameters
    // Have to be mixed for each grain size in the future!
    aspect_ratio += comp->getFraction() * comp->getAspectRatio();
    gold_g_factor += comp->getFraction() * comp->getGoldFactor();
    delta_rat += comp->getFraction() * comp->getDeltaRat();

    // Add all dust-to-gas mass ratios together
    dust_mass_fraction += comp->getDustMassFraction();
    fraction += comp->getFraction();

    // Check for scattering phase function (use HG if one or more components use HG)
    if(comp->getPhaseFunctionID() < phID)
        phID = comp->getPhaseFunctionID();

    // Use the lowest sublimation temperature (use multiple mixtures for higher accuracy
    // of the sublimation)
    if(sub_temp > comp->getSublimationTemperature())
        sub_temp = comp->getSublimationTemperature();

    if(comp->getComponentId() == comp->getNrOfComponents() - 1)
    {
        // Activate splines of Qtrq and Henyey-Greenstein g factor
        for(uint i = 0; i < nr_of_dust_species * nr_of_wavelength; i++)
        {
            Qtrq[i].createSpline();
            HG_g_factor[i].createSpline();
        }
    }

    // Create StringID for print parameter
    createStringID(comp);

    return true;
}

uint CDustComponent::getInteractingDust(CGridBasic * grid, photon_package * pp, uint cross_section)
{
    // Get wavelength index from photon package
    uint w = pp->getWavelengthID();

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // If only one grain size is used -> return this size (interpolation instead?)
    if(a_min == a_max)
        return CMathFunctions::biListIndexSearch(a_min, a_eff, nr_of_dust_species);

    switch(cross_section)
    {
        case CROSS_ABS:
            if(abs_prob != 0)
                return findSizeID(pp, abs_prob[w], a_min, a_max);
            break;

        case CROSS_SCA:
            if(sca_prob != 0)
                return findSizeID(pp, sca_prob[w], a_min, a_max);
            break;

        default:
            if(dust_prob != 0)
                return findSizeID(pp, dust_prob[w], a_min, a_max);
            break;
    }

    // Init pointer array for integration
    double * amount = new double[nr_of_dust_species];

    // Add relative amount of dust grains in each size bin
    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a, a_min, a_max))
        {
            switch(cross_section)
            {
                case CROSS_ABS:
                    amount[a] = a_eff_3_5[a] * getCabsMean(a, w);
                    break;

                case CROSS_SCA:
                    amount[a] = a_eff_3_5[a] * getCscaMean(a, w);
                    break;

                default:
                    amount[a] = a_eff_3_5[a];
                    break;
            }
        }
        else
            amount[a] = 0;
    }

    // Init propability list to pick random grain size
    prob_list prob(nr_of_dust_species);

    // Set each grain size bin with the relative amount
    CMathFunctions::probListInteg(a_eff, amount, prob);

    delete[] amount;

    uint a = findSizeID(pp, prob, a_min, a_max);

    return a;
}

void CDustComponent::calcPACrossSections(uint a, uint w, cross_sections & cs, double theta)
{
    // Init variables
    double sCext, dCext, sCsca, dCsca, sCabs, dCabs, sCcirc;
    double sinsq_th, cossq_th;

    // Get addition and subtraction of efficiencies (times PI * a^2 later)
    sCext = (getQext1(a, w) + getQext2(a, w)) / 2.0;
    dCext = (getQext1(a, w) - getQext2(a, w)) / 2.0;
    sCsca = (getQsca1(a, w) + getQsca2(a, w)) / 2.0;
    dCsca = (getQsca1(a, w) - getQsca2(a, w)) / 2.0;
    sCabs = (getQabs1(a, w) + getQabs2(a, w)) / 2.0;
    dCabs = (getQabs1(a, w) - getQabs2(a, w)) / 2.0;
    sCcirc = getQcirc(a, w) / 2.0;

    // Calculate sin(theta)^2
    sinsq_th = sin(theta);
    sinsq_th *= sinsq_th;

    // Calculate cos(theta)^2
    cossq_th = 1 - sinsq_th;

    // Calculate cross-sections of perfectly aligned dust grains
    cs.Cext = sCext + dCext * cossq_th;
    cs.Cpol = dCext * sinsq_th;
    cs.Csca = sCsca + dCsca * cossq_th;
    cs.Cabs = sCabs + dCabs * cossq_th;
    cs.Cpabs = dCext * sinsq_th;
    cs.Ccirc = sCcirc * sinsq_th;

    // Convert from efficiencies to cross-sections
    cs *= PI * a_eff_2[a];
}

void CDustComponent::calcCrossSections(CGridBasic * grid,
                                       photon_package * pp,
                                       double theta,
                                       uint i_density,
                                       uint a,
                                       cross_sections & cs)
{
    // Get wavelength index
    uint w = pp->getWavelengthID();

    // For random alignment use average cross-sections
    if(!is_align || alignment == ALIG_RND)
    {
        cs.Cext = getCextMean(a, w);
        cs.Cabs = getCabsMean(a, w);
        cs.Csca = getCscaMean(a, w);
        return;
    }

    // Perfect alignment can be calculated efficiently
    if((alignment & ALIG_PA) == ALIG_PA)
    {
        calcPACrossSections(a, w, cs, theta);
        return;
    }

    // Init variables
    double Rrat = 0, Rgold = 0, Ridg = 0, Rent = 1;
    double delta = 1;
    double a_alig = 1;

    // Get dust temperature from grid
    double Td;
    if(grid->getTemperatureFieldInformation() == TEMP_FULL)
        Td = grid->getDustTemperature(pp, i_density, a);
    else
        Td = grid->getDustTemperature(pp, i_density);

    // Get information from grid
    Vector3D B = grid->getMagField(pp);
    double Blen = B.length();
    double Tg = grid->getGasTemperature(pp);
    double ng = grid->getGasNumberDensity(pp);
    double a_limit = CMathFunctions::calc_larm_limit(Blen, Td, Tg, ng, aspect_ratio, larm_f);

    // Init cross-sections
    double Cext, Cpol, Cabs, Cpabs, Csca, Ccirc;

    // Calculate the parameters for radiative torque alignment
    if((alignment & ALIG_RAT) == ALIG_RAT)
    {
        a_alig = grid->getAlignedRadius(pp, i_density);
        if(a_eff[a] > a_alig)
        {
            if((alignment & ALIG_INTERNAL) == ALIG_INTERNAL)
                Rrat = f_highJ + (1 - f_highJ) * getInternalRAT();
            else
                Rrat = 1;
        }
    }

    // Calculate the parameters for GOLD mechanical alignment
    if((alignment & ALIG_GOLD) == ALIG_GOLD)
    {
        Vector3D v = grid->getVelocityField(pp);
        double vlength = v.length();
        double mach = vlength / sqrt(con_kB * Tg / (mu * m_H));

        if(mach > MACH_LIMIT)
        {
            Rgold = calcGoldReductionFactor(v, B);
            if((alignment & ALIG_INTERNAL) == ALIG_INTERNAL)
            {
                double cossq_zeta = getInternalGOLD(Td, Tg, vlength);
                Rgold *= 1.5 * (cossq_zeta - 1.0 / 3.0) * (1 + f_cor);
            }
        }
    }

    // Calculate the parameters for imperfect Davis-Greenstein alignment
    if((alignment & ALIG_IDG) == ALIG_IDG)
    {
        delta = delta0 * CMathFunctions::calc_delta(Blen, Td, Tg, ng);

        double zeta_sq = (a_eff[a] + delta * (Td / Tg)) / (a_eff[a] + delta);

        if(zeta_sq >= 1)
            zeta_sq = 0.999999;
        if(Td >= Tg)
            zeta_sq = 0.999999;

        double cossq_beta = (1 - sqrt(zeta_sq / (1 - zeta_sq)) * asin(sqrt(1 - zeta_sq))) / (1 - zeta_sq);

        Ridg = 1.5 * (cossq_beta - 1.0 / 3.0);

        if((alignment & ALIG_INTERNAL) == ALIG_INTERNAL)
        {
            double cossq_zeta = getInternalIDG(Td, Tg);
            Ridg *= 1.5 * (cossq_zeta - 1.0 / 3.0) * (1 + f_cor);
        }
    }

    if(a_eff[a] < a_limit)
        Rent = combinedRFactor(Ridg, Rrat, Rgold);
    else
        Rent = Ridg;

    // Calculate sin(theta)^2
    double sinsq_theta = sin(theta);
    sinsq_theta *= sinsq_theta;

    double c_s_ext = getCext1(a, w);
    double c_p_ext = getCext2(a, w);
    double c_av_ext = getCextMean(a, w);

    double cx_ext = c_av_ext + Rent / 3.0 * (c_s_ext - c_p_ext) * (1 - 3.0 * sinsq_theta);
    double cy_ext = c_av_ext + Rent / 3.0 * (c_s_ext - c_p_ext);

    cs.Cext = 0.5 * (cx_ext + cy_ext);
    cs.Cpol = 0.5 * Rent * (c_s_ext - c_p_ext) * sinsq_theta;

    double c_s_abs = getCabs1(a, w);
    double c_p_abs = getCabs2(a, w);
    double c_av_abs = getCabsMean(a, w);

    double cx_abs = c_av_abs + Rent / 3.0 * (c_s_abs - c_p_abs) * (1 - 3.0 * sinsq_theta);
    double cy_abs = c_av_abs + Rent / 3.0 * (c_s_abs - c_p_abs);

    cs.Cabs = 0.5 * (cx_abs + cy_abs);
    cs.Cpabs = 0.5 * Rent * (c_s_abs - c_p_abs) * sinsq_theta;

    double c_s_sca = getCsca1(a, w);
    double c_p_sca = getCsca2(a, w);
    double c_av_sca = getCscaMean(a, w);

    double cx_sca = c_av_sca + Rent / 3.0 * (c_s_sca - c_p_sca) * (1 - 3.0 * sinsq_theta);
    double cy_sca = c_av_sca + Rent / 3.0 * (c_s_sca - c_p_sca);

    cs.Csca = 0.5 * (cx_sca + cy_sca);
    cs.Ccirc = 0.5 * getCcirc(a, w) * Rent * sinsq_theta;
}

void CDustComponent::convertTempInQB(CGridBasic * grid,
                                     cell_basic * cell,
                                     uint i_density,
                                     double min_gas_density,
                                     bool use_gas_temp)
{
    if(grid->getGasNumberDensity(cell) < min_gas_density)
        return;

    // Set use dust offset later
    dust_offset = true;

    // Init variables
    double temp_offset;
    double * rel_weight = 0;

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, cell);
    double a_max = getSizeMax(grid, cell);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, cell);

    if(grid->getTemperatureFieldInformation() != TEMP_FULL)
    {
        // Get integration over the dust size distribution
        rel_weight = getRelWeight(a_min, a_max, size_param);
    }

    // Get temperature from grid as offset (assuming only ONE dust/gas temperature in
    // input grid)
    if(use_gas_temp)
        temp_offset = grid->getGasTemperature(cell);
    else
        temp_offset = grid->getDustTemperature(cell);

    // Check if temp_offset is larger zero
    if(temp_offset <= 0)
        return;

    // Find temperature index
    uint tID = findTemperatureID(temp_offset);

    // Init temporary pointer arrays for temperature and absorption rate
    double * qb_offset = new double[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a, a_min, a_max))
        {
            // Get offset energy of dust grain with tID
            qb_offset[a] = getQB(a, tID);

            if(grid->getTemperatureFieldInformation() == TEMP_FULL)
            {
                // Save offset energy to grid
                grid->setQBOffset(cell, i_density, a, qb_offset[a]);
            }
            else
            {
                // Multiply it with the amount of dust grains in the current bin
                qb_offset[a] *= rel_weight[a];
            }
        }
        else
        {
            // Set offset energy to zero if grain size is not used
            qb_offset[a] = 0;
        }
    }

    if(grid->getTemperatureFieldInformation() != TEMP_FULL)
    {
        // Get average absorption rate via interpolation
        double avg_qb_offset =
            CMathFunctions::integ_dust_size(a_eff, qb_offset, nr_of_dust_species, a_min, a_max);

        // Save average offset energy to grid
        grid->setQBOffset(cell, i_density, avg_qb_offset);
    }

    // Delete pointer array
    delete[] qb_offset;
    if(grid->getTemperatureFieldInformation() != TEMP_FULL)
        delete[] rel_weight;
}

bool CDustComponent::adjustTempAndWavelengthBW(CGridBasic * grid,
                                               photon_package * pp,
                                               uint i_density,
                                               bool use_energy_density)
{
    // Init variables
    double t_dust_new = 0;
    uint wIDnew = 0;

    // Get random number for wavelength of reemitting photon
    double rnd1 = pp->getRND();

    // Get the interacting dust grain size
    uint a = getInteractingDust(grid, pp);

    // Get dust temperature of the emitting grains in current cell
    t_dust_new = updateDustTemperature(grid, pp, i_density, a, use_energy_density);

    // If new dust grain temperature is larger zero -> start reemission
    if(t_dust_new > 0)
    {
        // Find temperature ID for the obtained temperature
        uint tIDnew = findTemperatureID(t_dust_new);

        // Pick a random wavelength from plank distribution related to the temperature
        wIDnew = findWavelengthID(a, tIDnew, rnd1);
    }
    else
        return false;

    // Update the energy of the photon package, in case of energy density
    if(use_energy_density)
        updateStokesVector(pp, wIDnew + 1);

    // Mol3D uses the upper value of the wavelength interval,
    // used for the selection of the emitting wavelengths from source!
    pp->setWavelengthID(wIDnew + 1);

    return true;
}

double CDustComponent::updateDustTemperature(CGridBasic * grid,
                                             photon_package * pp,
                                             uint i_density,
                                             uint a,
                                             bool use_energy_density)
{
    // Init variables
    double temp = 0, avg_temp = 0;

    // Get absorpion rate from grid
    double abs_rate = getAbsRate(grid, pp, a, use_energy_density);

    // If dust offset is true, add dust temperature from grid cell
    if(dust_offset)
    {
        if(grid->getTemperatureFieldInformation() == TEMP_FULL)
            abs_rate += grid->getQBOffset(pp, i_density, a);
        else if(grid->getTemperatureFieldInformation() == TEMP_EFF ||
                grid->getTemperatureFieldInformation() == TEMP_SINGLE)
            abs_rate += grid->getQBOffset(pp, i_density);
    }

    // Get temperature depending on the absorption rate (Minimum temperature is TEMP_MIN)
    temp = max(double(TEMP_MIN), findTemperature(a, abs_rate));

    if(sublimate)
        if(temp >= sub_temp)
            temp = 0;

    // Update min and max temperatures for visualization
    if(temp > max_temp)
        max_temp = temp;
    if(temp < min_temp)
        min_temp = temp;

    return temp;
}

void CDustComponent::calcTemperature(CGridBasic * grid,
                                     cell_basic * cell,
                                     uint i_density,
                                     bool use_energy_density)
{
    // Calculate the temperature only for cells with a density not zero
    if(getNumberDensity(grid, cell, i_density) == 0)
        return;

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, cell);
    double a_max = getSizeMax(grid, cell);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, cell);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    // Init temporary pointer arrays for absorption rate and temperature
    double * abs_rate = new double[nr_of_dust_species];
    double * temp = new double[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a, a_min, a_max))
        {
            // Check if dust grains should have been stochastically heated
            if(a_eff[a] <= getStochasticHeatingMaxSize())
            {
                // Init and resize spline for absorbed energy per wavelength
                spline abs_rate_per_wl;
                abs_rate_per_wl.resize(WL_STEPS);

                // Get radiation field and calculate absorbed energy for each wavelength
                for(uint w = 0; w < WL_STEPS; w++)
                {
                    double abs_rate_wl_tmp = grid->getRadiationField(cell, w) * getCabsMean(a, w);
                    abs_rate_per_wl.setValue(w, wavelength_list[w], abs_rate_wl_tmp);
                }

                // Activate spline of absorbed energy for each wavelength
                abs_rate_per_wl.createSpline();

                // Get pointer array of the temperature propabilities
                long double * temp_probability = getStochasticProbability(a, abs_rate_per_wl);

                // Reset absorpion rate
                abs_rate[a] = 0;

                // Set the temperature propabilities in the grid
                for(uint t = 0; t < getNrOfCalorimetryTemperatures(); t++)
                {
                    uint tID = findTemperatureID(calorimetry_temperatures[t]);
                    abs_rate[a] += temp_probability[t] * getQB(a, tID);
                }

                // Delete pointer array
                delete[] temp_probability;
            }
            else
            {
                // Get absorpion rate from grid
                abs_rate[a] = getAbsRate(grid, cell, a, use_energy_density);
            }

            // Add offset on absorption rate
            if(dust_offset)
            {
                if(grid->getTemperatureFieldInformation() == TEMP_FULL)
                    abs_rate[a] += grid->getQBOffset(cell, i_density, a);
                else if(grid->getTemperatureFieldInformation() == TEMP_EFF ||
                        grid->getTemperatureFieldInformation() == TEMP_SINGLE)
                    abs_rate[a] += grid->getQBOffset(cell, i_density);
            }

            // Calculate temperature from absorption rate
            temp[a] = max(double(TEMP_MIN), findTemperature(a, abs_rate[a]));

            // Consider sublimation temperature
            if(sublimate && grid->getTemperatureFieldInformation() == TEMP_FULL)
                if(temp[a] >= getSublimationTemperature())
                    temp[a] = 0;

            if(grid->getTemperatureFieldInformation() == TEMP_EFF ||
               grid->getTemperatureFieldInformation() == TEMP_SINGLE)
            {
                // Multiply with the amount of dust grains in the current bin for
                // integration
                abs_rate[a] *= rel_weight[a];
            }
            else if(grid->getTemperatureFieldInformation() == TEMP_FULL ||
                    (grid->getTemperatureFieldInformation() == TEMP_STOCH &&
                     a_eff[a] <= getStochasticHeatingMaxSize()))
            {
                // Set dust temperature in grid
                grid->setDustTemperature(cell, i_density, a, temp[a]);

                // Update min and max temperatures for visualization
                if(temp[a] > max_temp)
                    max_temp = temp[a];
                if(temp[a] < min_temp)
                    min_temp = temp[a];
            }

            // Multiply with the amount of dust grains in the current bin for integration
            temp[a] *= rel_weight[a];
        }
        else
        {
            // Set absorption rate to zero
            abs_rate[a] = 0;
            temp[a] = 0;
        }
    }

    double avg_temp;
    if(grid->getTemperatureFieldInformation() == TEMP_EFF ||
       grid->getTemperatureFieldInformation() == TEMP_SINGLE)
    {
        // Get average absorption rate via interpolation
        double avg_abs_rate =
            CMathFunctions::integ_dust_size(a_eff, abs_rate, nr_of_dust_species, a_min, a_max);

        // Calculate average temperature from absorption rate
        avg_temp = findTemperature(grid, cell, avg_abs_rate);
    }
    else
        avg_temp = CMathFunctions::integ_dust_size(a_eff, temp, nr_of_dust_species, a_min, a_max);

    // Delete pointer array
    delete[] rel_weight;
    delete[] abs_rate;
    delete[] temp;

    if(sublimate)
        if(avg_temp >= getSublimationTemperature())
        {
            // Set temperature to zero
            avg_temp = 0;

            // Remove sublimated dust from grid
            grid->adjustDustDensity(cell, i_density, 0);
        }

    // Set average dust temperature in grid
    grid->setDustTemperature(cell, i_density, max(double(TEMP_MIN), avg_temp));

    // Update min and max temperatures for visualization
    if(avg_temp > max_temp)
        max_temp = avg_temp;
    if(avg_temp < min_temp)
        min_temp = avg_temp;
}

void CDustComponent::calcAlignedRadii(CGridBasic * grid, cell_basic * cell, uint i_density)
{
    // Calculate the aligned radii only for cells with a density not zero
    if(getNumberDensity(grid, cell, i_density) == 0)
    {
        grid->setAlignedRadius(cell, i_density, a_eff[nr_of_dust_species - 1]);
        return;
    }

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, cell);
    double a_max = getSizeMax(grid, cell);

    // default value of the alignment radius
    double a_alig = getSizeMax(grid, cell);
    double th = 0;
    double dir = 0;

    // Aspect ratio of the grain
    double s = getAspectRatio();

    // alpha_1 ~ delta
    double alpha_1 = 1; // getDeltaRat();

    // Get grid values
    double T_gas = grid->getGasTemperature(cell);
    double n_g = grid->getGasNumberDensity(cell);
    double vol = grid->getVolume(cell);

    // Get average molecular weight
    double mu = grid->getMu();

    // Get thermal velocity
    double v_th = sqrt(2.0 * con_kB * T_gas / (mu * m_H));

    // Loop over all considered grain sizes
    double omega_old = 0;

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a, a_min, a_max))
        {
            // Get dust temperature from grid
            double T_dust;
            if(grid->getTemperatureFieldInformation() == TEMP_FULL ||
               (grid->getTemperatureFieldInformation() == TEMP_STOCH &&
                a_eff[a] <= getStochasticHeatingMaxSize()))
                T_dust = grid->getDustTemperature(cell, i_density, a);
            else
                T_dust = grid->getDustTemperature(cell, i_density);

            // Minor and major axis
            double a_minor = a_eff[a] * pow(s, 2. / 3.);
            double a_major = a_eff[a] * pow(s, -1. / 3.);

            // Moment of inertia along a_1
            double I_p = 8. * PI / 15. * getMaterialDensity(a) * a_minor * pow(a_major, 4);

            // Thermal angular momentum
            double J_th = sqrt(I_p * con_kB * T_gas);

            // Init. pointer arrays
            double * arr_product = new double[nr_of_wavelength];
            double * du = new double[nr_of_wavelength];
            double * ddir = new double[nr_of_wavelength];
            double * dth = new double[nr_of_wavelength];

            // Drag by gas
            double tau_gas = 3. / (4 * PIsq) * I_p / (mu * n_g * m_H * v_th * alpha_1 * pow(a_eff[a], 4));

            for(uint w = 0; w < nr_of_wavelength; w++)
            {
                // Init variables
                Vector3D en_dir;
                double arr_en_dens = 0;

                // Get radiation field (4 * PI * vol * J)
                grid->getSpecLength(cell, w, arr_en_dens, en_dir);

                // If the radiation field is zero -> set arrays to zero and move on
                if(arr_en_dens == 0)
                {
                    arr_product[w] = 0;
                    du[w] = 0;
                    ddir[w] = 0;
                    dth[w] = 0;
                    continue;
                }

                // Get angle between magnetic field and radiation field
                double theta = grid->getTheta(cell, en_dir);

                // Anisotropy parameter
                double gamma = en_dir.length() / arr_en_dens;

                // arr_en_dens = 4 * PI * vol * J -> 4 * PI / c * J
                arr_en_dens /= double(vol * con_c);

                du[w] = wavelength_list[w] * arr_en_dens;

                // Radiative torque efficiency as a power-law
                double Qr = 0.14;

                if(wavelength_list[w] > 1.8 * a_eff[a])
                    Qr = 0.14 / pow(wavelength_list[w] / (1.8 * a_eff[a]), 3);

                double cos_theta = abs(cos(theta));

                Qr *= cos_theta;

                // Qr=getQrat(a, w, 0.0);
                arr_product[w] =
                    arr_en_dens * (wavelength_list[w] / PIx2) * Qr * gamma * PI * pow(a_eff[a], 2);

                ddir[w] = wavelength_list[w] * arr_en_dens * gamma;
                dth[w] = wavelength_list[w] * arr_en_dens * cos_theta;
            }

            // Perform integration for total radiation field
            double u = CMathFunctions::integ(wavelength_list, du, 0, nr_of_wavelength - 1);
            dir = CMathFunctions::integ(wavelength_list, ddir, 0, nr_of_wavelength - 1);
            th = CMathFunctions::integ(wavelength_list, dth, 0, nr_of_wavelength - 1);

            dir /= u;
            th /= u;

            // drag by thermal emission
            double FIR = 1.40e10 * pow(u, 2. / 3.) / (a_eff[a] * n_g * sqrt(T_gas));

            // double FIR = CMathFunctions::integ(wavelength_list, dFIR, 0,
            // nr_of_wavelength - 1);
            double omega_frac = CMathFunctions::integ(wavelength_list, arr_product, 0, nr_of_wavelength - 1);

            double tau_drag = tau_gas / (1. + FIR);
            omega_frac *= tau_drag / J_th;

            // Delete pointer array
            delete[] arr_product;
            delete[] du;
            delete[] ddir;
            delete[] dth;

            if(omega_frac >= SUPERTHERMAL_LIMIT)
            {
                // linear interpolation
                if(a > 1)
                {
                    double a1 = a_eff[a - 1];
                    double a2 = a_eff[a];

                    double o1 = omega_old - SUPERTHERMAL_LIMIT;
                    double o2 = omega_frac - SUPERTHERMAL_LIMIT;

                    a_alig = a1 - o1 * (a2 - a1) / (o2 - o1);
                }
                else
                    a_alig = a_min;

                break;
            }

            // keep the prev. omega fraction for interpolation
            omega_old = omega_frac;
        }
    }

    // Check for proper size range
    if(a_alig < a_min)
        a_alig = a_min;

    if(a_alig > a_max)
        a_alig = a_max;

    // Set aligned grain size in grid
    grid->setAlignedRadius(cell, i_density, a_alig);
    grid->setAvgDir(cell, dir);
    grid->setAvgTheta(cell, th);

    // Update aligned grain size limits
    if(a_alig < min_a_alig)
        min_a_alig = a_alig;
    if(a_alig > max_a_alig)
        max_a_alig = a_alig;
}

double CDustComponent::calcGoldReductionFactor(Vector3D & v, Vector3D & B)
{
    // Init variables
    double s;
    double g = gold_g_factor;
    Vector3D v_proj = (v * B) / B.sq_length() * B;
    double len_vz = v_proj.sq_length();
    double len_vxyz = v.length();
    double R, cossq_beta;

    len_vxyz *= len_vxyz;

    if(len_vxyz == len_vz)
        return -0.499999999999;

    if(len_vxyz == 3 * len_vz)
        return 0;

    s = -0.5 * (len_vxyz - 3 * len_vz) / (len_vxyz - len_vz);

    if(s <= -0.5)
        s = -0.4999999999;

    if(s < 0)
    {
        if(g < 0)
        {
            cossq_beta =
                (sqrt(-g) * asin(sqrt(-s / (1 + g)))) / (s * atan(sqrt((s * g) / (1 + s + g)))) - 1 / s;
        }
        else
        {
            cossq_beta =
                (sqrt(g) * asin(sqrt(-s / (1 + g)))) / (s * atan(sqrt((-s * g) / (1 + s + g)))) - 1 / s;
        }
    }
    else
    {
        if(g < 0)
        {
            cossq_beta =
                (sqrt(-g) * asinh(sqrt(s / (1 + g)))) / (s * atanh(sqrt(-(s * g) / (1 + s + g)))) - 1 / s;
        }
        else
        {
            cossq_beta =
                (sqrt(g) * asinh(sqrt(s / (1 + g)))) / (s * atanh(sqrt((s * g) / (1 + s + g)))) - 1 / s;
        }
    }

    R = 1.5 * (cossq_beta - 1.0 / 3.0);

    if(R <= -0.5)
        R = -0.499999999999;

    return R;
}

void CDustComponent::calcStochasticHeatingPropabilities(CGridBasic * grid,
                                                        cell_basic * cell,
                                                        uint i_density,
                                                        dlist & wavelength_list_full)
{
    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, cell);
    double a_max = getSizeMax(grid, cell);

    if(getNumberDensity(grid, cell, i_density) == 0)
    {
        // If density is zero, set propability to 1 for lowest temperature
        for(uint a = 0; a < nr_of_dust_species; a++)
            if(a_eff[a] <= getStochasticHeatingMaxSize() && sizeIndexUsed(a, a_min, a_max))
                grid->setDustTempProbability(cell, i_density, a, 0, 1.0);
        return;
    }

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        // Check if dust grains should have been stochastically heated
        if(a_eff[a] <= getStochasticHeatingMaxSize() && sizeIndexUsed(a, a_min, a_max))
        {
            // Init and resize spline for absorbed energy per wavelength
            spline abs_rate_per_wl;
            abs_rate_per_wl.resize(WL_STEPS);

            // Get radiation field and calculate absorbed energy for each wavelength
            for(uint w = 0; w < WL_STEPS; w++)
                abs_rate_per_wl.setValue(
                    w, wavelength_list_full[w], grid->getRadiationField(cell, w) * getCabsMean(a, w));

            // Activate spline of absorbed energy for each wavelength
            abs_rate_per_wl.createSpline();

            // Get pointer array of the temperature propabilities
            long double * temp_probability = getStochasticProbability(a, abs_rate_per_wl);

            // Set the temperature propabilities in the grid
            for(uint t = 0; t < getNrOfCalorimetryTemperatures(); t++)
                grid->setDustTempProbability(cell, i_density, a, t, temp_probability[t]);

            // Delete pointer array
            delete[] temp_probability;
        }
    }
}

double CDustComponent::getCalorimetryA(uint a, uint f, uint i, spline & abs_rate_per_wl)
{
    // Calculation of A from Eq. (23) in Camps et al. (2015)

    // Init varibales
    double res = 0;

    if(f > i)
    {
        // Calculate difference between enthalpy of two temperatures
        double enthalpy_diff = getEnthalpy(a, f) - getEnthalpy(a, i);

        if(enthalpy_diff > 0)
        {
            // Get wavelength according to the enthalpy difference
            double wavelength = con_h * con_c / (enthalpy_diff);

            // Calculate A
            res = abs_rate_per_wl.getValue(wavelength) * con_h * con_c * getEnthalpyBinWidth(a, f) /
                  pow(enthalpy_diff, 3);
        }
    }
    else if(f == i - 1)
    {
        // Calculate difference between enthalpy of two temperatures
        double enthalpy_diff = getEnthalpy(a, i) - getEnthalpy(a, f);

        // Get calorimetric temperature to index i
        double T_i = getCalorimetricTemperature(i);

        // Calculate A
        res = tab_em_inv[a].getValue(T_i) * PIx4 / enthalpy_diff;
    }
    else
        cout << "\nHINT: Error at the getCalorimetryA fuction (stochastic heating)";

    // Returning A needs to be at least zero
    return max(double(0), res);
}

long double * CDustComponent::getStochasticProbability(uint a, spline & abs_rate_per_wl)
{
    // Calculate the propability of a dust grain to have a certain temperature
    // See Sect. 4.3. in Camps et al. (2015)

    // Init variables
    Matrix2D B_mat;

    // Resize 2D matrix
    B_mat.resize(nr_of_calorimetry_temperatures, nr_of_calorimetry_temperatures);

    // Go through each possible transition between temperature states and set B matrix
    // values
    for(uint f = nr_of_calorimetry_temperatures - 1; f > 0; f--)
        if(f == nr_of_calorimetry_temperatures - 1)
            for(uint i = 0; i < nr_of_calorimetry_temperatures - 1; i++)
                B_mat(f, i) = getCalorimetryA(a, f, i, abs_rate_per_wl);
        else
            for(uint i = 0; i < f; i++)
                B_mat(f, i) = (B_mat(f + 1, i) + getCalorimetryA(a, f, i, abs_rate_per_wl));

    // Init propability array
    long double * X_vec = new long double[nr_of_calorimetry_temperatures];

    // Set first value to 1
    X_vec[0] = numeric_limits<long double>::min();
    for(uint i = 1; i < nr_of_calorimetry_temperatures; i++)
    {
        // Set the other values to zero
        X_vec[i] = 0;

        // Calculate propability from lower temperatures
        long double caloA = getCalorimetryA(a, i - 1, i, abs_rate_per_wl);
        if(caloA > 0)
        {
            for(uint j = 0; j < i; j++)
                X_vec[i] += (long double)B_mat(i, j) * X_vec[j] / caloA;
        }
        else
            X_vec[i] = 0;
    }

    // Init sum for normalization
    long double X_sum = 0;

    // Calculate the sum of the propability
    for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
        X_sum += X_vec[t];

    // Perform normalization or  reset propability
    for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
        if(isinf(X_sum) || isnan(X_sum))
        {
            X_vec[t] = 0;
            cout << "\nERROR: Wrong values in stochastic heating calculation!" << endl;
        }
        else if(X_sum > 0)
            X_vec[t] = X_vec[t] / X_sum;
        else
            X_vec[t] = 0;

    return X_vec;
}

StokesVector CDustComponent::calcEmissivitiesHz(CGridBasic * grid, photon_package * pp, uint i_density)
{
    // Get extinction and absorption cross-sections
    double Cext = getCextMean(grid, pp);
    double Cabs = getCabsMean(grid, pp);

    // Get wavelength index
    uint w = pp->getWavelengthID();

    // Calculate frequency
    double frequency = con_c / wavelength_list[w];

    // Get Planck emission at this frequency
    double pl_hz = CMathFunctions::planck_hz(frequency, grid->getDustTemperature(pp, i_density));

    // Get number density of dust grains
    double dens_dust = getNumberDensity(grid, pp, i_density);

    // Fill Stokes vector including optical depth
    StokesVector tmp_stokes(Cabs * pl_hz * dens_dust, 0, 0, 0, Cext * dens_dust);

    return tmp_stokes;
}

double CDustComponent::calcEmissivities(CGridBasic * grid, photon_package * pp, uint i_density)
{
    // Init variables to calculate the emission/extinction
    double temp_dust, pl_abs;

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    // Get wavelength index of photon package
    uint w = pp->getWavelengthID();

    // Init temporary array for integration
    double * pl_abs_tmp = new double[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a, a_min, a_max))
        {
            // Calculate emission/extinction according to the information inside of the grid
            if(a_eff[a] <= getStochasticHeatingMaxSize())
            {
                // Consider stochastic heating for the emission if chosen
                for(uint t = 0; t < getNrOfCalorimetryTemperatures(); t++)
                {
                    // Get current calorimetric temperature
                    temp_dust = getCalorimetricTemperature(t);

                    // Get propability that the dust grains have the current calorimetric
                    // temperature
                    double temp_probability = grid->getDustTempProbability(pp, i_density, a, t);

                    // Add relative emissivity from this temperature
                    pl_abs_tmp[a] +=
                        getCabsMean(a, w) * rel_weight[a] * temp_probability * getTabPlanck(w, temp_dust);
                }
            }
            else
            {
                // Consider the temperature of every dust grain size or an average
                // temperature
                if(grid->getTemperatureFieldInformation() == TEMP_FULL)
                    temp_dust = grid->getDustTemperature(pp, i_density, a);
                else
                    temp_dust = grid->getDustTemperature(pp, i_density);

                // Calculate the emission of the dust grains
                pl_abs_tmp[a] = getCabsMean(a, w) * rel_weight[a] * getTabPlanck(w, temp_dust);
            }
        }
        else
        {
            // Init planck lists
            pl_abs_tmp[a] = 0;
        }
    }

    // Perform integration to obtain a more precise result for the efficient emission and
    // cross section
    pl_abs = CMathFunctions::integ_dust_size(a_eff, pl_abs_tmp, nr_of_dust_species, a_min, a_max);

    // Delete pointer array
    delete[] rel_weight;
    delete[] pl_abs_tmp;

    // Multiply with number density
    pl_abs *= getNumberDensity(grid, pp, i_density);

    return pl_abs;
}

StokesVector CDustComponent::getRadFieldScatteredFraction(CGridBasic * grid,
                                                          photon_package * pp,
                                                          uint i_density,
                                                          Vector3D en_dir,
                                                          double energy)
{
    if(energy > 1e-200)
    {
        // Get local min and max grain sizes
        double a_min = getSizeMin(grid, pp);
        double a_max = getSizeMax(grid, pp);

        // Get local size parameter for size distribution
        double size_param = getSizeParam(grid, pp);

        // Init  and calculate the cross-sections
        cross_sections cs;

        // Get wavelength index of photon package
        uint w = pp->getWavelengthID();

        // Get theta of scattering
        double scattering_theta = acos(en_dir * pp->getDirection());

        // Get angle between the magnetic field and the photon direction
        double mag_field_theta = grid->getThetaMag(pp);

        // Get integration over the dust size distribution
        double * rel_weight = getRelWeight(a_min, a_max, size_param);

        uint thID;
        if(phID == PH_MIE)
        {
            // Get index of theta scattering
            thID = getScatThetaID(scattering_theta);
        }

        // Init temporary Stokes array for integration
        StokesVector * scatter_stokes = new StokesVector[nr_of_dust_species];

        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            if(sizeIndexUsed(a, a_min, a_max))
            {
                // Get cross sections and relative weight of the current dust grain size
                calcCrossSections(grid, pp, mag_field_theta, i_density, a, cs);

                // Multiply energy with scattered fraction in the theta angle and the
                // scattering cross-section
                scatter_stokes[a].setI(energy * rel_weight[a] * cs.Csca *
                                       getScatteredFraction(a, w, scattering_theta));

                if(phID == PH_MIE)
                {
                    // Get scattering matrix
                    const Matrix2D & mat_sca = getScatteringMatrix(a, w, 0, 0, thID);

                    // Multiply Stokes vector with scattering matrix
                    double i_1 = scatter_stokes[a].I();
                    if(i_1 > 1e-200)
                    {
                        scatter_stokes[a] = mat_sca * scatter_stokes[a];
                        scatter_stokes[a] *= i_1 / scatter_stokes[a].I();
                    }
                    else
                        scatter_stokes[a].clear();
                }
            }
        }

        // Perform integration for the emission
        StokesVector final_stokes(
            CMathFunctions::integ_dust_size(a_eff, scatter_stokes, nr_of_dust_species, a_min, a_max));

        // Delete pointer arrays
        delete[] rel_weight;
        delete[] scatter_stokes;

        // Get rotation angle to rotate back into the map/detector frame
        double phi_map = CMathFunctions::getRotationAngleObserver(en_dir, pp->getEY(), pp->getEX());

        // Rotate Stokes Vector to be in agreement with the detector plane
        final_stokes.rot(phi_map);

        // Multiply with number density
        final_stokes *= getNumberDensity(grid, pp, i_density);

        return final_stokes;
    }
    else
    {
        return StokesVector();
    }
}

StokesVector CDustComponent::calcEmissivitiesEmi(CGridBasic * grid,
                                                 photon_package * pp,
                                                 uint i_density,
                                                 double phi,
                                                 double energy,
                                                 Vector3D en_dir)
{
    // Init variables to calculate the emission/extinction
    double temp_dust = 0;
    double scattering_theta = 0, mag_field_theta = 0, phi_map = 0;
    uint temp_info = grid->getTemperatureFieldInformation();
    uint thID = 0;

    // Precalculate values for scattering
    if(energy > 1e-200)
    {
        scattering_theta = acos(en_dir * pp->getDirection());
        phi_map = CMathFunctions::getRotationAngleObserver(en_dir, pp->getEY(), pp->getEX());
        if(phID == PH_MIE)
        {
            // Get index of theta scattering
            thID = getScatThetaID(scattering_theta);
        }
    }

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Init  and calculate the cross-sections
    cross_sections cs;

    // Get wavelength index of photon package
    uint w = pp->getWavelengthID();

    // Get angle between the magnetic field and the photon direction
    mag_field_theta = grid->getThetaMag(pp);

    // Calculate orientation of the Stokes vector in relation to the magnetic field
    double sin_2ph = sin(2.0 * phi);
    double cos_2ph = cos(2.0 * phi);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    // Init temporary Stokes array for integration
    StokesVector * tmp_stokes = new StokesVector[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a, a_min, a_max))
        {
            // Get cross sections and relative weight of the current dust grain size
            calcCrossSections(grid, pp, mag_field_theta, i_density, a, cs);

            // Calculate emission/extinction according to the information inside of the
            // grid
            if(getEffectiveRadius(a) <= getStochasticHeatingMaxSize())
            {
                // Consider stochastic heating for the emission if chosen
                for(uint t = 0; t < getNrOfCalorimetryTemperatures(); t++)
                {
                    // Get current calorimetric temperature
                    temp_dust = getCalorimetricTemperature(t);

                    // Get propability that the dust grains have the current calorimetric
                    // temperature
                    double pl = grid->getDustTempProbability(pp, i_density, a, t);

                    // Get relative Planck emission
                    pl *= rel_weight[a] * getTabPlanck(w, temp_dust);

#ifdef CAMPS_BENCHMARK
                    // To perform Camps et. al (2015) benchmark.
                    tmp_stokes[a].addI(cs.Cabs * pl);
#else
                    // Add relative emissivity from this temperature
                    tmp_stokes[a].addI(cs.Cabs * pl);
                    tmp_stokes[a].addQ(cs.Cpabs * pl * cos_2ph);
                    tmp_stokes[a].addU(cs.Cpabs * pl * sin_2ph);
#endif
                }
            }
            else
            {
                // Consider the temperature of every dust grain size or an average
                // temperature
                if(temp_info == TEMP_FULL)
                    temp_dust = grid->getDustTemperature(pp, i_density, a);
                else
                    temp_dust = grid->getDustTemperature(pp, i_density);

                double pl = rel_weight[a] * getTabPlanck(w, temp_dust);

#ifdef CAMPS_BENCHMARK
                // To perform Camps et. al (2015) benchmark.
                tmp_stokes[a].addQ(cs.Cabs * pl;
#else
                // Add relative emissivity from this temperature
                tmp_stokes[a].addI(cs.Cabs * pl);
                tmp_stokes[a].addQ(cs.Cpabs * pl * cos_2ph);
                tmp_stokes[a].addU(cs.Cpabs * pl * sin_2ph);
#endif
            }

            // Add scattering component if radiation field is stored in grid
            if(energy > 1e-200)
            {
                // Init variables
                StokesVector scatter_stokes;

                // Multiply energy with scattered fraction in the theta angle and the
                // scattering cross-section
                scatter_stokes.setI(energy * rel_weight[a] * cs.Csca *
                                    getScatteredFraction(a, w, scattering_theta));

                if(phID == PH_MIE)
                {
                    // Get scattering matrix
                    const Matrix2D & mat_sca = getScatteringMatrix(a, w, 0, 0, thID);

                    // Multiply Stokes vector with scattering matrix
                    double i_1 = scatter_stokes.I();
                    if(i_1 > 1e-200)
                    {
                        scatter_stokes = mat_sca * scatter_stokes;
                        scatter_stokes *= i_1 / scatter_stokes.I();
                    }
                    else
                        scatter_stokes.clear();
                }

                // Rotate Stokes Vector to be in agreement with the detector plane
                scatter_stokes.rot(phi_map);

#ifndef CAMPS_BENCHMARK
                // Add scattered light to the Stokes vector
                tmp_stokes[a].addS(scatter_stokes);
#endif
            }
        }
    }

    // Perform integration for the emission
    StokesVector final_stokes(
        CMathFunctions::integ_dust_size(a_eff, tmp_stokes, nr_of_dust_species, a_min, a_max));

    // Delete pointer arrays
    delete[] rel_weight;
    delete[] tmp_stokes;

    // Multiply with number density
    final_stokes *= getNumberDensity(grid, pp, i_density);

    return final_stokes;
}

void CDustComponent::calcExtCrossSections(CGridBasic * grid,
                                          photon_package * pp,
                                          uint i_density,
                                          double & avg_Cext,
                                          double & avg_Cpol,
                                          double & avg_Ccirc)
{
    // Init  and calculate the cross-sections
    cross_sections cs;

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    // Get angle between the magnetic field and the photon direction
    double mag_field_theta = grid->getThetaMag(pp);

    // Init temporary cross-section array for integration
    double * Cext = new double[nr_of_dust_species];
    double * Cpol = new double[nr_of_dust_species];
    double * Ccirc = new double[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a, a_min, a_max))
        {
            // Get cross sections and relative weight of the current dust grain size
            calcCrossSections(grid, pp, mag_field_theta, i_density, a, cs);

            // Add relative cross-sections for integration
            Cext[a] = cs.Cext * rel_weight[a];
            Cpol[a] = cs.Cpol * rel_weight[a];
            Ccirc[a] = cs.Ccirc * rel_weight[a];
        }
        else
        {
            // Set cross-sections to zero for the unused grain size
            Cext[a] = 0;
            Cpol[a] = 0;
            Ccirc[a] = 0;
        }
    }

    // Perform integration for the cross-sections
    avg_Cext = CMathFunctions::integ_dust_size(a_eff, Cext, nr_of_dust_species, a_min, a_max);
    avg_Cpol = CMathFunctions::integ_dust_size(a_eff, Cpol, nr_of_dust_species, a_min, a_max);
    avg_Ccirc = CMathFunctions::integ_dust_size(a_eff, Ccirc, nr_of_dust_species, a_min, a_max);

    // Delete pointer arrays
    delete[] rel_weight;
    delete[] Cext;
    delete[] Cpol;
    delete[] Ccirc;
}

photon_package CDustComponent::getEscapePhoton(CGridBasic * grid,
                                               photon_package * pp,
                                               uint a,
                                               Vector3D obs_ex,
                                               Vector3D dir_obs)
{
    switch(phID)
    {
        case PH_MIE:
            return getEscapePhotonMie(grid, pp, a, obs_ex, dir_obs);

        default:
        {
            // Init variables
            double len, dens, Cext, tau_obs = 0;

            // Get wavelength index of the photon package
            uint w = pp->getWavelengthID();

            // Determination of the scattering angle (phi, theta) towards the observing
            // map in the photon frame. Get the rotation matrix of the photon (photon
            // space to lab space)
            Matrix2D D_photon = pp->getD();
            D_photon.transpose();
            Vector3D dir_rlp = D_photon * dir_obs;

            // Calculate the theta angle to the observer
            double theta_photon_to_obs = acos(dir_rlp.Z());

            // Calculate the fraction that is scattered into this theta direction
            double scattered_fraction = getScatteredFraction(a, w, theta_photon_to_obs);

            // Get the Stokes vector of the current photon package
            StokesVector tmp_stokes = pp->getStokesVector();

            // Reduce the photon package Stokes vector by albedo and scattering fraction
            tmp_stokes *= scattered_fraction * getCscaMean(a, w) / getCextMean(a, w);

            // Init temporary photon package
            photon_package pp_res;

            // Set the photon package at the position of the current photon
            pp_res.setPosition(pp->getPosition());
            pp_res.setPositionCell(pp->getPositionCell());

            // Synchronize the direction and wavelength as well
            pp_res.setDirection(dir_obs);
            pp_res.setWavelengthID(w);

            // Transport the photon package through the grid
            while(grid->next(&pp_res))
            {
                // Get the traveled distance
                len = pp_res.getTmpPathLength();

                // Get the current density
                dens = getNumberDensity(grid, &pp_res);

                // Get the current mean extinction cross-section
                Cext = getCextMean(grid, &pp_res);

                // Add the optical depth of the current path to the total optical depth
                tau_obs += Cext * len * dens;
            }

            // Reduce the Stokes vector by the optical depth
            tmp_stokes *= exp(-tau_obs);

            // Set the new Stokes vector to the photon package
            pp_res.setStokesVector(tmp_stokes);

            return pp_res;
        }
    }
}

photon_package CDustComponent::getEscapePhotonMie(CGridBasic * grid,
                                                  photon_package * pp,
                                                  uint a,
                                                  Vector3D obs_ex,
                                                  Vector3D dir_obs)
{
    // Init variables
    double len, dens, Cext, phi_fraction = 1, tau_obs = 0;

    // Get wavelength index of the photon package
    uint w = pp->getWavelengthID();

    // Get the Stokes vector of the current photon package
    StokesVector tmp_stokes = pp->getStokesVector();

    // Init temporary photon package
    photon_package pp_res;

    // Set the photon package at the position of the current photon
    pp_res.setPosition(pp->getPosition());
    pp_res.setPositionCell(pp->getPositionCell());

    // Synchronize the direction and wavelength as well
    pp_res.setD(pp->getD());
    pp_res.setWavelengthID(w);

    // Determination of the scattering angle (phi, theta) towards the observing map in the
    // photon frame. Get the rotation matrix of the photon (photon space to lab space)
    Matrix2D D_photon = pp->getD();
    D_photon.transpose();
    Vector3D dir_rlp = D_photon * dir_obs;

    // Calculate the theta and phi angle to the observer
    double phi_photon_to_obs = atan3(dir_rlp.Y(), -dir_rlp.X());
    double theta_photon_to_obs = acos(dir_rlp.Z());

    // Update the coordinate space of the photon
    pp_res.updateCoordSystem(phi_photon_to_obs, theta_photon_to_obs);

    // Get the theta angle index to obtain the scattering matrix
    uint thID = getScatThetaID(theta_photon_to_obs);

    // Create the scattering matrix with the local parameters
    const Matrix2D & mat_sca = getScatteringMatrix(a, w, 0, 0, thID);

    if(phID == PH_MIE)
    {
        // Get PHIPAR to take non equal distribution of phi angles into account
        double PHIPAR = (sqrt(pow(tmp_stokes.Q(), 2) + pow(tmp_stokes.U(), 2)) / tmp_stokes.I()) *
                        (-mat_sca(0, 1) / mat_sca(0, 0));

        // Get cos(2 * phi)
        double cos_2_phi = 1.0 - 2.0 * pow(sin(phi_photon_to_obs), 2);

        // Calculate the fraction that is scattered into this phi direction
        phi_fraction = (1.0 - PHIPAR * cos_2_phi);
    }

    // Calculate the fraction that is scattered into this theta direction
    double theta_fraction = getScatteredFractionMie(a, w, theta_photon_to_obs);

    // Reduce Stokes vector by albedo and scattering propability into theta and phi
    tmp_stokes *= theta_fraction * phi_fraction * getCscaMean(a, w) / getCextMean(a, w);

    // Backup Stokes vector
    double stokes_1_bak = tmp_stokes.I();

    // Rotate Stokes vector to new photon direction
    tmp_stokes.rot(phi_photon_to_obs);

    // Multiply Stokes vector with scattering matrix
    tmp_stokes = mat_sca * tmp_stokes;

    // Normalize Stokes vector to preserve total intensity
    tmp_stokes *= stokes_1_bak / tmp_stokes.I();

    // Transport the photon package through the grid
    while(grid->next(&pp_res))
    {
        // Get the traveled distance
        len = pp_res.getTmpPathLength();

        // Get the current density
        dens = getNumberDensity(grid, &pp_res);

        // Get the current mean extinction cross-section
        Cext = getCextMean(grid, &pp_res);

        // Add the optical depth of the current path to the total optical depth
        tau_obs += Cext * len * dens;
    }

    // Reduce the Stokes vector by the optical depth
    tmp_stokes *= exp(-tau_obs);

    // Rotate photon package into the coordinate space of the detector
    double rot_angle_phot_obs =
        CMathFunctions::getRotationAngleObserver(obs_ex, pp_res.getEX(), pp_res.getEY());
    tmp_stokes.rot(rot_angle_phot_obs);

    // Set the new Stokes vector to the photon package
    pp_res.setStokesVector(tmp_stokes);

    return pp_res;
}

double CDustComponent::getCellEmission(CGridBasic * grid, photon_package * pp, uint i_density)
{
    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    // Get current cell
    cell_basic * cell = pp->getPositionCell();

    // Get Volume of current cell
    double vol = grid->getVolume(cell);

    // Get dust number density of current cell
    double dens = getNumberDensity(grid, cell, i_density);

    // Get wavelength of photon package
    uint w = pp->getWavelengthID();

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        // Get dust temperature from grid
        double temp;
        if(grid->getTemperatureFieldInformation() == TEMP_FULL)
            temp = grid->getDustTemperature(pp, i_density, a);
        else
            temp = grid->getDustTemperature(pp, i_density);

        // Calculate energy of current grain size
        rel_weight[a] *= getCabsMean(a, w) * getTabPlanck(w, temp);
    }

    // Calclate the total energy via interpolation
    double total_energy =
        dens * vol * PIx4 *
        CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);

    // Delete pointer array
    delete[] rel_weight;

    return total_energy;
}

void CDustComponent::henyeygreen(photon_package * pp, uint a, bool adjust_stokes)
{
    // Init variables
    double cos_theta, theta, phi;
    double g = 0;

    // Get two random numbers for the new direction
    double z1 = pp->getRND();
    double z2 = pp->getRND();

    // Get the current wavelength
    double w = pp->getWavelengthID();

    // Get the Henyey-Greenstein g factor
    g = getHGg(a, w);

    // If g is very close to zero, use random direction
    if(abs(g) < 0.5e-5)
    {
        pp->calcRandomDirection();
        pp->updateCoordSystem();
        return;
    }

    // Set g factor to the boundaries if larger/smaller
    if(g < -1)
        g = -0.99999;
    if(g > 1)
        g = 0.99999;

    // Get cosine theta from g factor
    cos_theta = (1.0 - g * g) / (1.0 - g + 2.0 * g * z1);
    cos_theta = 1 + g * g - cos_theta * cos_theta;
    cos_theta /= (2.0 * g);

    // Get theta from cosine theta
    theta = acos(cos_theta);

    // equal distribution of the phi angles
    phi = PIx2 * z2;

    // Update the photon package with the new direction
    pp->updateCoordSystem(phi, theta);

    if(adjust_stokes)
    {
        StokesVector tmp_stokes = pp->getStokesVector();
        double theta_fraction = getScatteredFraction(w, a, theta);
        tmp_stokes *= theta_fraction * getCscaMean(a, w) / getCextMean(a, w);
        pp->setStokesVector(tmp_stokes);
    }
}

void CDustComponent::miesca(photon_package * pp, uint a, bool adjust_stokes)
{
    // Init variables
    double HELP, phi, phi1, PHIPAR = 0, GAMMA = 1, hd1, hd2;

    // Get wavelength of photon package
    uint w = pp->getWavelengthID();

    // Get theta angle from distribution
    double theta = findTheta(a, w, pp->getRND());

    // Get theta index from theta angle
    uint thID = getScatThetaID(theta);

    // Get Stokes vector from photon package
    StokesVector tmp_stokes = pp->getStokesVector();

    // Get scattering matrix
    const Matrix2D & mat_sca = getScatteringMatrix(a, w, 0, 0, thID);

    // Get PHIPAR to take non equal distribution of phi angles into account
    PHIPAR = (sqrt(tmp_stokes.Q() * tmp_stokes.Q() + tmp_stokes.U() * tmp_stokes.U()) / tmp_stokes.I()) *
             (-mat_sca(0, 1) / mat_sca(0, 0));

    // Obtain phi angle
    bool hl1 = false;
    double rndx = pp->getRND();
    if(abs(PHIPAR) < 0.1)
        phi = rndx * PIx2;
    else
    {
        uint run_counter = 0;
        do
        {
            HELP = rndx * PIx4;
            phi = HELP + PHIPAR * sin(HELP);
            phi1 = 0.0;

            if(run_counter > 1000)
            {
                phi = pp->getRND() * PIx2;
                break;
            }

            do
            {
                hd1 = abs(phi - phi1);
                phi1 = phi;
                phi = HELP + PHIPAR * sin(phi1);
                hd2 = abs(phi - phi1);

                run_counter++;

                if(run_counter > 1000)
                    break;

                if(abs(phi - phi1) <= 0.0175)
                {
                    break;
                }
                else
                {
                    if(abs(hd1 - hd2) < 1.0e-15)
                    {
                        hl1 = true;
                        break;
                    }
                }

            } while(true);

        } while(hl1);

        phi = phi / 2.0;

        if(tmp_stokes.Q() != 0.0)
        {
            GAMMA = 0.5 * atan2(tmp_stokes.U(), tmp_stokes.Q());
            if(tmp_stokes.U() < 0.0)
                GAMMA = PI + GAMMA;
        }
        else
        {
            if(tmp_stokes.U() < 0.0)
                GAMMA = PI4x3;
            else if(tmp_stokes.U() > 0.0)
                GAMMA = PI4;
        }

        phi = PI - GAMMA + phi;
    }

    // Update the photon package with the new direction
    pp->updateCoordSystem(phi, theta);

    if(adjust_stokes)
    {
        tmp_stokes *= getCscaMean(a, w) / getCextMean(a, w);

        double i_1 = tmp_stokes.I();
        tmp_stokes.rot(phi);
        tmp_stokes = mat_sca * tmp_stokes;
        tmp_stokes *= i_1 / tmp_stokes.I();

        pp->setStokesVector(tmp_stokes);
    }
}

bool CDustMixture::createDustMixtures(parameters & param, string path_data, string path_plot)
{
    // Do not load dust component if not required
    uint nr_of_total_components = param.getTotalNrOfDustComponents();
    if(nr_of_total_components == 0)
    {
        if(param.getCommand() == CMD_LINE_EMISSION || param.getCommand() == CMD_OPIATE ||
           param.getCommand() == CMD_SYNCHROTRON)
            return true;
        else
            return false;
    }

    // dust_choices_to_index takes the different dust_choices and gives the
    // index of the dust mixture (from 0, 3, 5 -> 0, 1, 2)
    dust_choices_to_index.resize(param.getMaxDustComponentChoice() + 1);

    // dust_choices includes all set dust_i_mixture values (example: 1, 3, 8, 66)
    dust_choices = param.getDustComponentChoices();

    // Connect the dust_choices with the id of the final dust mixture
    uint nr_of_dust_mixtures = dust_choices.size();
    for(uint i_mixture = 0; i_mixture < nr_of_dust_mixtures; i_mixture++)
        for(uint dust_choice = 0; dust_choice <= param.getMaxDustComponentChoice(); dust_choice++)
            if(dust_choices[i_mixture] == dust_choice)
                dust_choices_to_index[dust_choice] = i_mixture;

    // Get number of dust mixtures and read the parameters files of their components
    mixed_component = new CDustComponent[nr_of_dust_mixtures];
    for(uint i_mixture = 0; i_mixture < nr_of_dust_mixtures; i_mixture++)
    {
        // Set up the relation between dust_i_mixture and the "real" dust index used in
        // this code.
        uilist unique_components;

        // Get the "dust choice" of the current dust component
        uint current_dust_choice = dust_choices[i_mixture];

        // Combine in this loop only dust components which have the same "dust_choice"
        nr_of_components = 0;
        for(uint i_comp = 0; i_comp < nr_of_total_components; i_comp++)
            if(param.getDustChoiceFromComponentId(i_comp) == current_dust_choice)
            {
                nr_of_components++;
                unique_components.push_back(i_comp);
            }

        // Init the minimum and maximum grain size limits of the dust mixture
        double a_min_mixture = 1e200, a_max_mixture = 0;

        // Find the minimum and maximum grain size limits of the dust mixture
        for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
        {
            // Get the global id of the current dust component
            uint dust_component_choice = unique_components[i_comp];

            // Get min and max grain sizes of the dust component
            double a_min = param.getSizeMin(dust_component_choice);
            double a_max = param.getSizeMax(dust_component_choice);

            // Use the highest a_max and lowest a_min of the components
            if(a_min_mixture > a_min)
                a_min_mixture = a_min;
            if(a_max_mixture < a_max)
                a_max_mixture = a_max;
        }

        // Set the minimum and maximum grain size limits of the dust mixture
        mixed_component[i_mixture].setSizeMin(a_min_mixture);
        mixed_component[i_mixture].setSizeMax(a_max_mixture);

        // Init the sum of the current dust composition
        double fraction_sum = 0;

        // Calculate the sum of all fractions
        for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
        {
            // Get the global id of the current dust component
            uint dust_component_choice = unique_components[i_comp];

            double fraction = param.getDustFraction(dust_component_choice);
            fraction_sum += fraction;
        }

        /* Not necessary
        if(!param.getIndividualDustMassFractions() && fraction_sum != 1.0)
        {
            cout << "\nERROR: Fractions of dust materials do not add up to 100 % "
                 << "(" << sum << "%)!" << endl;
            return false;
        }
        */

        // Init single components pointer array
        single_component = new CDustComponent[nr_of_components];
        for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
        {
            // Get the global id of the current dust component
            uint dust_component_choice = unique_components[i_comp];

            // Set the ID and number of the components
            setIDs(single_component[i_comp], i_comp, nr_of_components, i_mixture, nr_of_dust_mixtures);

            // Get dust component fraction of mixture
            double fraction = param.getDustFraction(dust_component_choice);

            // Set mass fractions of each
            if(param.getIndividualDustMassFractions())
            {
                single_component[i_comp].setIndividualDustMassFractions(true);
                single_component[i_comp].setDustMassFraction(fraction);
                single_component[i_comp].setFraction(fraction / fraction_sum);
            }
            else
            {
                single_component[i_comp].setDustMassFraction(fraction * param.getDustMassFraction() /
                                                             fraction_sum);
                single_component[i_comp].setFraction(fraction / fraction_sum);
            }

            // Get size distribution parameters
            string size_keyword = param.getDustSizeKeyword(dust_component_choice);
            single_component[i_comp].setSizeParameter(size_keyword,
                                                      param.getDustSizeParameter(dust_component_choice));

            // Get material density and similar user defined parameters
            single_component[i_comp].setMaterialDensity(param.getMaterialDensity(dust_component_choice));
            single_component[i_comp].setFHighJ(param.getFHighJ());
            single_component[i_comp].setDelta0(param.getDelta0());
            single_component[i_comp].setLarmF(param.getLarmF());
            single_component[i_comp].setMu(param.getMu());
            single_component[i_comp].setPhaseFunctionID(param.getPhaseFunctionID());

            // Get global wavelength grid
            single_component[i_comp].setWavelengthList(wavelength_list, wavelength_offset);

            // Get Path to dust parameters file
            string path = param.getDustPath(dust_component_choice);

            // Read dust input data (refractive index or optical properties)
            if(path.find(".nk") != std::string::npos)
            {
                // Read refractive index and use Mie theory to get optical properties
                if(!single_component[i_comp].readDustRefractiveIndexFile(
                       param, dust_component_choice, a_min_mixture, a_max_mixture))
                    return false;
            }
            else
            {
                // Read the dust catalog file for each dust component (including scatter
                // matrix)
                if(!single_component[i_comp].readDustParameterFile(param, dust_component_choice))
                    return false;
            }

            // Read the calorimetric file for each dust component
            if(param.getStochasticHeatingMaxSize() > single_component[i_comp].getSizeMin())
                if(!single_component[i_comp].readCalorimetryFile(param, dust_component_choice))
                {
                    cout << "\nERROR: Cannot open calorimetry file, which is required "
                            "for stochastic heating!"
                         << endl;
                    return false;
                }

            // Write dust component files, if multiple components will be mixed together
            if(nr_of_components > 1)
                single_component[i_comp].writeComponent(path_data, path_plot);
        }

        // Set the ID and number of the mixtures
        setIDs(mixed_component[i_mixture], 0, nr_of_components, i_mixture, nr_of_dust_mixtures);

        // Mix components together
        if(!mixComponents(param, i_mixture))
            return false;

        // Delete single components
        killSingleComponents();
    }

    // Init foreground extinction
    extinction_magnitude = param.getForegroundExtinctionMagnitude();
    extinction_magnitude_wavelength = param.getForegroundExtinctionWavelength();
    extinction_dust_mixture = param.getForegroundExtinctionDustMixture();

    return true;
}

void CDustMixture::printParameter(parameters & param, CGridBasic * grid)
{
    // If no mixture was defined, show basic information
    if(getNrOfMixtures() == 0)
    {
        cout << CLR_LINE;
        cout << "Dust parameters (No dust component defined!)                            "
                "            "
             << endl;
        cout << SEP_LINE;
        return;
    }

    // Show title
    cout << CLR_LINE;
    cout << "Dust parameters                                                             "
            "            "
         << endl;
    cout << SEP_LINE;

    // Show the full wavelength grid for temo and RAT calculation
    if(param.isMonteCarloSimulation() ||
       (param.getCommand() == CMD_DUST_EMISSION && param.getStochasticHeatingMaxSize() > 0))
        cout << "- Number of wavelengths   : " << WL_STEPS << "  (" << WL_MIN << " [m] - " << WL_MAX
             << " [m])" << endl;

    // Monte-Carlo scattering is only used for temp, rat and scatter maps
    if(param.isMonteCarloSimulation() || param.getCommand() == CMD_DUST_SCATTERING ||
       param.getScatteringToRay())
        cout << "- Phase function          : " << getPhaseFunctionStr() << endl;

    // Enforced first scattering method is only used for Monte-Carlo scattering maps
    if(param.getCommand() == CMD_DUST_SCATTERING)
    {
        cout << "- Enforced first scat.    : ";
        if(param.getEnfScattering())
            cout << "enabled          " << endl;
        else
            cout << "disabled         " << endl;
    }

    // Show information about dust emission simulations
    if(param.getCommand() == CMD_DUST_EMISSION)
    {
        cout << "- Temperature distr.      : ";
        if(grid->getTemperatureFieldInformation() == TEMP_FULL)
        {
            cout << "temperatures for all grain sizes found in grid" << endl;
            if(param.getStochasticHeatingMaxSize() > 0)
                cout << "                            calculate stochastic heating up to "
                        "grain size of "
                     << param.getStochasticHeatingMaxSize() << " [m]" << endl;
        }
        else if(grid->getTemperatureFieldInformation() == TEMP_STOCH)
        {
            cout << "temperatures for effective grain size found in grid" << endl
                 << "                            including stochastically heated "
                    "temperatures"
                 << endl;
            if(param.getStochasticHeatingMaxSize() > 0)
                cout << "    HINT: Stochastic heating was already calculated. This "
                        "should not happen!"
                     << endl;
        }
        else if(grid->getTemperatureFieldInformation() == TEMP_EFF)
        {
            cout << "temperature for effective grain size found in grid" << endl;
            if(param.getStochasticHeatingMaxSize() > 0)
                cout << "                            calculate stochastic heating up to "
                        "grain size of "
                     << param.getStochasticHeatingMaxSize() << " [m]" << endl;
        }
        else if(grid->getTemperatureFieldInformation() == TEMP_SINGLE)
        {
            cout << "temperature for effective grain size found in grid" << endl;
            cout << "                            single temperature for all density "
                    "distributions!"
                 << endl;
        }
        else
            cout << "not available (This should not happen!)" << endl;

        if(param.getAligRAT())
        {
            cout << "- Alignment radii         : ";
            if(grid->useDustChoice() && grid->getNrAlignedRadii() == 1)
                cout << "found a common radius for all dust mixtures" << endl;
            else if(grid->useDustChoice() && grid->getNrAlignedRadii() == 1)
                cout << "found a common radius for all dust mixtures and density dist." << endl;
            else if(!grid->useDustChoice() && grid->getNrAlignedRadii() == getNrOfMixtures())
                cout << "found a separate radius for each density distribution" << endl;
            else
                cout << "ERROR: This should not happen!" << endl;
        }

        cout << "- Include scattered light : ";
        if(grid->getRadiationFieldAvailable())
        {
            if(param.getScatteringToRay())
            {
                cout << "yes, based on the radiation field" << endl
                     << "    HINT: Only one dominant radiation source and mostly single "
                        "scattering?"
                     << endl
                     << "          -> If not, use <rt_scattering> 0                      "
                        "          "
                     << endl;
            }
            else
                cout << "no, disabled via <rt_scattering> 0" << endl;
        }
        else
        {
            if(param.getScatteringToRay())
                cout << "yes, radiation field will be calculated before raytracing" << endl;
            else
            {
                cout << "no, radiation field not found in grid and radiation sources "
                        "missing"
                     << endl
                     << "                            "
                     << "try to define point source(s) and/or the dust source" << endl;
            }
        }

        cout << "- Foreground Extinction   : ";
        if(extinction_magnitude > 0)
        {
            cout << "A_lambda = " << extinction_magnitude
                 << " at wavelength = " << extinction_magnitude_wavelength << " [m] " << endl;
            if(extinction_dust_mixture == MAX_UINT)
                cout << "                            based on ISM grains with MRN-size distribution" << endl;
            else
                cout << "                            based on dust mixture id = " << extinction_dust_mixture
                     << endl;
        }
        else
            cout << "no, enable via <foreground_extinction> A_lambda wavelength" << endl;

        cout << "Observed wavelengths:" << endl;

        dlist dust_ray_detectors = param.getDustRayDetectors();
        for(uint i = 0; i < dust_ray_detectors.size(); i += NR_OF_RAY_DET)
        {
            // Index of current detector
            uint pos = i / NR_OF_RAY_DET;

            if(uint(dust_ray_detectors[i + 2]) > 1)
                cout << "- Emission detetector " << (pos + 1) << "   : from " << dust_ray_detectors[i + 0]
                     << " [m] to " << dust_ray_detectors[i + 1] << " [m] with "
                     << uint(dust_ray_detectors[i + 2]) << " logarithmic values" << endl;
            else if(uint(dust_ray_detectors[i + 2]) == 1)
                cout << "- Emission detetector " << (pos + 1) << "   : " << dust_ray_detectors[i + 0]
                     << " [m]" << endl;
        }
    }

    // Show information about dust scattering simulations
    if(param.getCommand() == CMD_DUST_SCATTERING)
    {
        cout << "- Scattering method       : ";
        if(param.getPeelOff())
            cout << "use peel-off technique" << endl;
        else
            cout << "use acceptance angle (" << param.getAcceptanceAngle() << "°)" << endl;
        cout << "Observed wavelengths:" << endl;

        dlist dust_mc_detectors = param.getDustMCDetectors();
        for(uint i = 0; i < dust_mc_detectors.size(); i += NR_OF_MC_DET)
        {
            uint pos = i / NR_OF_MC_DET;

            if(uint(dust_mc_detectors[i + 2]) > 1)
                cout << "- Scattering detetector " << (pos + 1) << " : from " << dust_mc_detectors[i + 0]
                     << " [m]) to " << dust_mc_detectors[i + 1] << " [m]) with "
                     << uint(dust_mc_detectors[i + 2]) << " logarithmic values" << endl;
            else if(uint(dust_mc_detectors[i + 2]) == 1)
                cout << "- Scattering detetector " << (pos + 1) << "   : " << dust_mc_detectors[i + 0]
                     << " [m]" << endl;
        }
    }

    // Show information about dust temperature distribution simulations
    if(param.isTemperatureSimulation())
    {
        cout << "- Temperature calculation : ";
        if(param.getDustTempMulti())
            cout << "for all grain sizes" << endl;
        else
            cout << "for effective grain size" << endl;

        if(param.getStochasticHeatingMaxSize() > 0 && !param.getSaveRadiationField())
            cout << "                          : including stochastic heating up to "
                 << param.getStochasticHeatingMaxSize() << " [m]" << endl;
        else if(param.getStochasticHeatingMaxSize() > 0 && param.getSaveRadiationField())
        {
            cout << "\nHINT: Stochastic heating and saving the radiation field is chosen." << endl
                 << "      The radiation field will be saved and stochastic heating "
                    "should be set"
                 << endl
                 << "      with CMD_DUST_EMISSION simulation." << endl;
        }
    }

    // Show information about saving the radiation field
    if(param.isMonteCarloSimulation())
    {
        cout << "- Save radiation field    : ";
        if(param.getSaveRadiationField())
            cout << "yes (stochastic heating and full raytracing possible)" << endl;
        else
            cout << "no" << endl;
    }

    // Show information about each dust mixture
    for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
    {
        cout << SEP_LINE;
        cout << "Dust mixture " << (i_mixture + 1) << "/" << getNrOfMixtures() << " (Choice ID "
             << param.getDustChoiceFromMixtureId(i_mixture) << ")" << endl;

        cout << "- Avg. grain mass         : " << getAvgMass(i_mixture) << " [kg]" << endl;

        if(param.getCommand() == CMD_DUST_EMISSION && !param.getAligRANDOM())
        {
            if(mixed_component[i_mixture].isAligned())
                cout << "- Affected by alignment   : Yes" << endl;
            else
                cout << "- Affected by alignment   : No" << endl;
        }

        double total_dust_mass = 0;
        for(long i_cell = 0; i_cell < grid->getMaxDataCells(); i_cell++)
        {
            cell_basic * cell = grid->getCellFromIndex(i_cell);
            total_dust_mass += getMassDensity(grid, cell, i_mixture) * grid->getVolume(cell);
        }
        cout << "- Total mass              : " << total_dust_mass / M_sun << " [M_sun], " << total_dust_mass
             << " [kg]" << endl;
        cout << mixed_component[i_mixture].getStringID();
    }
    cout << SEP_LINE;

    // If dust optical properties are constant -> pre calc some values
    if(grid->useConstantGrainSizes())
        preCalcEffProperties(param);
}

bool CDustMixture::mixComponents(parameters & param, uint i_mixture)
{
    // If only one component is defined for this mixture -> change only pointer
    if(nr_of_components == 1)
    {
        // Point mixture pointer on single component
        mixed_component[i_mixture] = single_component[0];

        // Set single component to mixture
        mixed_component[i_mixture].setIsMixture(true);

        // Reset pointer of single component
        single_component = 0;

        // Create StringID for print parameter
        mixed_component[i_mixture].createStringID(&mixed_component[i_mixture]);

        // Pre-calculate various quantities
        if(!preCalcDustProperties(param, i_mixture))
            return false;

        return true;
    }

    // Get global wavelength grid (including bin width)
    mixed_component[i_mixture].setWavelengthList(wavelength_list, wavelength_offset);

    // Set Scat matrix loaded and calorimetry loaded to true
    // as long as no component has not loaded them
    mixed_component[i_mixture].setScatLoaded(true);
    mixed_component[i_mixture].setCalorimetryLoaded(true);

    // Get Relative fractions for each size bin
    double *** size_fraction = getSizeFractions();

    // Get common parameters grid
    uint nr_of_dust_species = single_component[0].getNrOfDustSpecies();
    uint nr_of_incident_angles = single_component[0].getNrOfIncidentAngles();

    for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
    {
        // Check if the components have the same amount of grain sizes
        if(nr_of_dust_species != single_component[i_comp].getNrOfDustSpecies())
        {
            cout << "\nERROR: Component Nr. " << i_comp + 1 << " has a different amount of dust species!"
                 << endl;
            return false;
        }

        // Check if the components have the same amount of incident angles
        if(nr_of_incident_angles != single_component[i_comp].getNrOfIncidentAngles())
        {
            cout << "\nERROR: Component Nr. " << i_comp + 1 << " has a different amount of incident angles!"
                 << endl;
            return false;
        }

        // If scattering matrix was not loaded for one component -> not loaded for the
        // mixture
        if(!single_component[i_comp].getScatLoaded())
            mixed_component[i_mixture].setScatLoaded(false);

        // If calorimetry data was not loaded for one component .> not loaded for the
        // mixture
        if(!single_component[i_comp].getCalorimetryLoaded())
            mixed_component[i_mixture].setCalorimetryLoaded(false);

        // Only if no component can be aligned, do not use alignment of mixture
        if(single_component[i_comp].isAligned())
            mixed_component[i_mixture].setIsAligned(true);

        if(single_component[i_comp].getIndividualDustMassFractions())
            mixed_component[i_mixture].setIndividualDustMassFractions(true);

        // Add parameters of each component together to the mixture
        if(!mixed_component[i_mixture].add(size_fraction[i_comp], &single_component[i_comp]))
            return false;
    }

    // Pre-calculate various quantities
    if(!preCalcDustProperties(param, i_mixture))
        return false;

    // Delete pointer array
    delete[] size_fraction;

    return true;
}

bool CDustMixture::preCalcDustProperties(parameters & param, uint i_mixture)
{
    // Set various parameters for the mixture
    mixed_component[i_mixture].setSublimate(param.isSublimate());

    // Calculate wavelength differences for temperature (reemission)
    if(param.isMonteCarloSimulation())
        if(!mixed_component[i_mixture].calcWavelengthDiff())
        {
            cout << "\nERROR: The wavelength grid only has one wavelength which is not "
                    "enough for temp calculation!\n"
                 << "       Change WL_STEPS to more than one in the typedefs.h file!" << endl;
            return false;
        }

    // Create the temperature grid
    if(mixed_component[i_mixture].getCalorimetryLoaded() &&
       param.getStochasticHeatingMaxSize() > mixed_component[i_mixture].getSizeMin())
    {
        // Adjust the temperature grid limits to fit with calorimetry data
        double temp_min = mixed_component[i_mixture].getMinCalorimetricTemperature();
        double temp_max = mixed_component[i_mixture].getMaxCalorimetricTemperature();
        mixed_component[i_mixture].preCalcTemperatureLists(temp_min, temp_max, TEMP_STEP);

        // If calorimetry data was loaded, allow stochastic heating
        mixed_component[i_mixture].setStochasticHeatingMaxSize(param.getStochasticHeatingMaxSize());
    }
    else
        mixed_component[i_mixture].preCalcTemperatureLists(TEMP_MIN, TEMP_MAX, TEMP_STEP);

    // Pre-calculate reemission probability of the dust mixture at different temperatures
    if(param.isMonteCarloSimulation())
        mixed_component[i_mixture].preCalcWaveProb();

    // Pre-calculate temperature to total emission relation (either for temperature or
    // stochastic heating calculation)
    if(param.isMonteCarloSimulation() ||
       (param.getCommand() == CMD_DUST_EMISSION &&
        param.getStochasticHeatingMaxSize() > mixed_component[i_mixture].getSizeMin()))
        mixed_component[i_mixture].preCalcAbsorptionRates();

    // If scattering matrix was loaded, pre-calculate phase function and scattering
    // distribution
    if(mixed_component[i_mixture].getScatLoaded())
        mixed_component[i_mixture].preCalcMieScatteringProb();
    else if(mixed_component[i_mixture].getPhaseFunctionID() == PH_MIE)
    {
        cout << "Error Mie scattering is chosen, but no scattering matrix was read!" << endl;
        return false;
    }

    // Use random alignment for Monte-Carlo simulations
    if(param.isMonteCarloSimulation() || param.getCommand() == CMD_DUST_SCATTERING)
        mixed_component[i_mixture].setAlignmentMechanism(ALIG_RND);
    else
        mixed_component[i_mixture].setAlignmentMechanism(param.getAlignmentMechanism());

    return true;
}
