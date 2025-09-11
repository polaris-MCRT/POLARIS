/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include <string.h>
#include <cmath>
#include "DustComponent.hpp"
#include "CommandParser.hpp"
#include "GridBasic.hpp"
#include "MathFunctions.hpp"
#include "Typedefs.hpp"
#include "Parameters.hpp"

void CDustComponent::initDustProperties()
{
    // Init pointer arrays for grain sizes and size distributions
    a_eff = new double[nr_of_dust_species];
    grain_distribution_x_aeff_sq = new double[nr_of_dust_species];
    grain_size_distribution = new double[nr_of_dust_species];
    a_eff_squared = new double[nr_of_dust_species];
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
    HGg2 = new double *[nr_of_dust_species];
    HGg3 = new double *[nr_of_dust_species];

    CextMean = new double *[nr_of_dust_species];
    CabsMean = new double *[nr_of_dust_species];
    CscaMean = new double *[nr_of_dust_species];

    // Set pointer array values to zero or add second dimension
    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        a_eff[a] = 0;
        grain_distribution_x_aeff_sq[a] = 0;
        grain_size_distribution[a] = 0;
        a_eff_squared[a] = 0;
        mass[a] = 0;

        Qext1[a] = new double[nr_of_wavelength];
        Qext2[a] = new double[nr_of_wavelength];
        Qabs1[a] = new double[nr_of_wavelength];
        Qabs2[a] = new double[nr_of_wavelength];
        Qsca1[a] = new double[nr_of_wavelength];
        Qsca2[a] = new double[nr_of_wavelength];
        Qcirc[a] = new double[nr_of_wavelength];
        HGg[a] = new double[nr_of_wavelength];
        HGg2[a] = new double[nr_of_wavelength];
        HGg3[a] = new double[nr_of_wavelength];

        CextMean[a] = new double [nr_of_wavelength];
        fill(CextMean[a], CextMean[a] + nr_of_wavelength, 0);
        CabsMean[a] = new double [nr_of_wavelength];
        fill(CabsMean[a], CabsMean[a] + nr_of_wavelength, 0);
        CscaMean[a] = new double [nr_of_wavelength];
        fill(CscaMean[a], CscaMean[a] + nr_of_wavelength, 0);

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
            HGg2[a][w] = 0;
            HGg3[a][w] = 1;
        }
    }

    // Qtrq and parameters for Henyey-Greenstein phase function are splines over the incident angle
    Qtrq = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g_factor = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g2_factor = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g3_factor = new spline[nr_of_dust_species * nr_of_wavelength];
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

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        // Second init of the scattering matrix
        sca_mat[a] = new Matrix2D ***[nr_of_wavelength];
        for(uint w = 0; w < nr_of_wavelength; w++)
        {
            // Increase counter used to show progress
            per_counter++;

            // Calculate percentage of total progress per source
            float percentage = 100.0 * float(per_counter) / float(max_counter);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
                printIDs();
                cout << "- allocating memory: " << percentage << " [%]                      \r";
                last_percentage = percentage;
            }

            // Third init of the scattering matrix
            sca_mat[a][w] = new Matrix2D **[nr_of_incident_angles];
            if(nr_of_scat_theta[a][w]!=0)
            {
                for(uint inc = 0; inc < nr_of_incident_angles; inc++)
                {
                    // Fourth init of the scattering matrix
                    sca_mat[a][w][inc] = new Matrix2D *[nr_of_scat_phi];
                    for(uint sph = 0; sph < nr_of_scat_phi; sph++)
                    {
                        // Fifth init of the scattering matrix
                        sca_mat[a][w][inc][sph] = new Matrix2D[nr_of_scat_theta[a][w]];
                        for(uint sth = 0; sth < nr_of_scat_theta[a][w]; sth++)
                        {
                            // sixth init of the scattering matrix
                            sca_mat[a][w][inc][sph][sth].resize(4, 4);
                        }
                    }
                }
            }
        }
    }
}

void CDustComponent::initScatteringMatrixArray(uint nr_of_scat_theta_tmp)
{
    // Init counter and percentage to show progress
    ullong per_counter = 0;
    float last_percentage = 0;

    // Init maximum counter value
    uint max_counter = nr_of_wavelength * nr_of_dust_species;

    // First init of the scattering matrix and scattering angle spline
    sca_mat = new Matrix2D ****[nr_of_dust_species];
    scat_theta = new double **[nr_of_dust_species];
    nr_of_scat_theta = new uint *[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        // Second init of the scattering matrix and scattering angle
        sca_mat[a] = new Matrix2D ***[nr_of_wavelength];
        scat_theta[a] = new double *[nr_of_wavelength];
        nr_of_scat_theta[a] = new uint [nr_of_wavelength];

        for(uint w = 0; w < nr_of_wavelength; w++)
        {
            // Increase counter used to show progress
            per_counter++;

            // Calculate percentage of total progress per source
            float percentage = 100.0 * float(per_counter) / float(max_counter);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
                printIDs();
                cout << "- allocating memory: " << percentage << " [%]                      \r";
                last_percentage = percentage;
            }

            // Third init of the scattering matrix
            sca_mat[a][w] = new Matrix2D **[nr_of_incident_angles];
            scat_theta[a][w] = new double[nr_of_scat_theta_tmp];
            nr_of_scat_theta[a][w] = nr_of_scat_theta_tmp;

            for(uint inc = 0; inc < nr_of_incident_angles; inc++)
            {
                // Fourth init of the scattering matrix
                sca_mat[a][w][inc] = new Matrix2D *[nr_of_scat_phi];
                for(uint sph = 0; sph < nr_of_scat_phi; sph++)
                {
                    // Fifth init of the scattering matrix
                    sca_mat[a][w][inc][sph] = new Matrix2D[nr_of_scat_theta_tmp];
                    for(uint sth = 0; sth < nr_of_scat_theta_tmp; sth++)
                    {
                        // sixth init of the scattering matrix
                        sca_mat[a][w][inc][sph][sth].resize(4, 4);
                    }
                }
            }
        }
    }
}

void CDustComponent::initNrOfScatThetaArray()
{
    nr_of_scat_theta = new uint *[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
        nr_of_scat_theta[a] = new uint [nr_of_wavelength]();
}

void CDustComponent::initScatThetaArray()
{
    scat_theta = new double **[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        scat_theta[a] = new double *[nr_of_wavelength];
        for(uint w = 0; w < nr_of_wavelength; w++)
            scat_theta[a][w] = new double [nr_of_scat_theta[a][w]]();
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
    spline *eff_wl, *Qtrq_wl, *HG_g_factor_wl, *HG_g2_factor_wl, *HG_g3_factor_wl;
    uint nr_of_wavelength_dustcat;

    // Get min and max dust grain size
    double a_min = param.getSizeMin(dust_component_choice);
    double a_max = param.getSizeMax(dust_component_choice);

    // Init text file reader for dust cat file
    fstream reader(path.c_str());

    // Error message if the read does not work
    if(reader.fail())
    {
        cout << ERROR_LINE << "Cannot open dust parameters file:" << endl;
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
                    cout << ERROR_LINE << "Line " << line_counter << " wrong amount of numbers!" << endl;
                    return false;
                }

                // The number of dust grain sizes
                nr_of_dust_species = (uint)values[0];

                // The number of wavelength used by the dust catalog
                nr_of_wavelength_dustcat = (uint)values[1];

                // The number of incident angles
                nr_of_incident_angles = (uint)values[2];

                // The aspect ratio of minor to major dust grain axes
                aspect_ratio = values[3];

                // The material density (only used if no one was set in the command file)
                if(material_density == 0)
                {
                    if(values[4] == 0)
                    {
                        printIDs();
                        cout << ERROR_LINE << "dust bulk mass is zero!" << endl;
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

                // Init splines for wavelength interpolation of the dust optical properties
                eff_wl = new spline[nr_of_dust_species * NR_OF_EFF];
                Qtrq_wl = new spline[nr_of_dust_species * nr_of_incident_angles];
                HG_g_factor_wl = new spline[nr_of_dust_species * nr_of_incident_angles];
                HG_g2_factor_wl = new spline[nr_of_dust_species * nr_of_incident_angles];
                HG_g3_factor_wl = new spline[nr_of_dust_species * nr_of_incident_angles];

                break;

            case 3:
                // The second line needs a values per dust grain size
                if(values.size() != nr_of_dust_species)
                {
                    cout << ERROR_LINE << "Line " << line_counter << " wrong amount of effective radii!" << endl;
                    return false;
                }

                // Init arrays for grain size, size distribution, and mass
                a_eff = new double[nr_of_dust_species];
                a_eff_squared = new double[nr_of_dust_species];
                grain_distribution_x_aeff_sq = new double[nr_of_dust_species];
                grain_size_distribution = new double[nr_of_dust_species];
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
                    cout << ERROR_LINE << "Line " << line_counter << " wrong amount of wavelength!" << endl;
                    return false;
                }

                // Init an array for the wavelengths of the dust grain catalog
                wavelength_list_dustcat.resize(nr_of_wavelength_dustcat);

                // Fill the wavelength array
                for(uint w = 0; w < nr_of_wavelength_dustcat; w++)
                    wavelength_list_dustcat[w] = values[w];

                if(wavelength_list[0] < wavelength_list_dustcat[0] ||
                    wavelength_list[nr_of_wavelength - 1] > wavelength_list_dustcat[nr_of_wavelength_dustcat - 1])
                {
                    cout << WARNING_LINE << "The wavelength range is out of the limits of the catalog. This may cause problems!\n"
                        << "         wavelength range          : " << wavelength_list[0] << " [m] to "
                        << wavelength_list[nr_of_wavelength - 1] << " [m]\n"
                        << "         wavelength range (catalog): " << wavelength_list_dustcat[0] << " [m] to "
                        << wavelength_list_dustcat[nr_of_wavelength_dustcat - 1] << " [m]" << endl;
                    if(!IGNORE_WAVELENGTH_RANGE)
                    {
                        cout << "         To continue, set 'IGNORE_WAVELENGTH_RANGE' to 'true' in src/Typedefs.h and recompile!" << endl;
                        return false;
                    }
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

                // Each line contains NR_OF_EFF - 1 plus nr_of_incident_angles values
                if(values.size() == NR_OF_EFF - 1 + nr_of_incident_angles)
                {
                    // the current line contains enough values
                    rec = true;

                    // Set the dust grain optical properties
                    for(uint i = 0; i < NR_OF_EFF - 1; i++)
                        eff_wl[a * NR_OF_EFF + i].setValue(w, wavelength_list_dustcat[w], values[i]);

                    // For each incident angles, get Qtrq
                    for(uint i_inc = 0; i_inc < nr_of_incident_angles; i_inc++)
                    {
                        // At the first wavelength, resize the splines for the wavelengths
                        if(w == 0)
                        {
                            Qtrq_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g2_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g3_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                        }

                        // Set the radiative torque efficiency
                        Qtrq_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w, wavelength_list_dustcat[w], values[NR_OF_EFF - 1 + i_inc]);

                        // Set the Henyey-Greenstein g factor
                        HG_g_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            0.0);
                        
                        // Set the second Henyey-Greenstein factor to 0
                        HG_g2_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            0.0);
                        
                        // Set the third Henyey-Greenstein factor to 1
                        HG_g3_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            1.0);
                    }

                    // Increace the counter for the dust grain optical properties
                    eff_counter++;
                }

                // Alternatively, each line contains additional column with g for HG
                if(values.size() == NR_OF_EFF - 1 + 2 * nr_of_incident_angles)
                {
                    // the current line contains enough values
                    rec = true;

                    // Set the dust grain optical properties
                    for(uint i = 0; i < NR_OF_EFF - 1; i++)
                        eff_wl[a * NR_OF_EFF + i].setValue(w, wavelength_list_dustcat[w], values[i]);

                    // For each incident angles, get Qtrq and HG parameter
                    for(uint i_inc = 0; i_inc < nr_of_incident_angles; i_inc++)
                    {
                        // At the first wavelength, resize the splines for the wavelengths
                        if(w == 0)
                        {
                            Qtrq_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g2_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g3_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                        }

                        // Set the radiative torque efficiency
                        Qtrq_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w, wavelength_list_dustcat[w], values[NR_OF_EFF - 1 + i_inc]);

                        // Set the Henyey-Greenstein g factor
                        HG_g_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            values[NR_OF_EFF - 1 + i_inc + nr_of_incident_angles]);

                        // Set the second Henyey-Greenstein factor
                        HG_g2_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            0.0);

                        // Set the third Henyey-Greenstein factor to 0
                        HG_g3_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            1.0);
                    }

                    // Increace the counter for the dust grain optical properties
                    eff_counter++;
                }

                // Alternatively, each line contains additional two columns with g and alpha for Draine HG
                if(values.size() == NR_OF_EFF - 1 + 3 * nr_of_incident_angles)
                {
                    // the current line contains enough values
                    rec = true;

                    // Set the dust grain optical properties
                    for(uint i = 0; i < NR_OF_EFF - 1; i++)
                        eff_wl[a * NR_OF_EFF + i].setValue(w, wavelength_list_dustcat[w], values[i]);

                    // For each incident angles, get Qtrq and HG parameters
                    for(uint i_inc = 0; i_inc < nr_of_incident_angles; i_inc++)
                    {
                        // At the first wavelength, resize the splines for the wavelengths
                        if(w == 0)
                        {
                            Qtrq_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g2_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g3_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                        }

                        // Set the radiative torque efficiency
                        Qtrq_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w, wavelength_list_dustcat[w], values[NR_OF_EFF - 1 + i_inc]);

                        // Set the Henyey-Greenstein g factor
                        HG_g_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            values[NR_OF_EFF - 1 + i_inc + nr_of_incident_angles]);

                        // Set the second Henyey-Greenstein factor
                        HG_g2_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            values[NR_OF_EFF - 1 + i_inc + 2 * nr_of_incident_angles]);

                        // Set the third Henyey-Greenstein factor to 1
                        HG_g3_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            1.0);
                    }

                    // Increace the counter for the dust grain optical properties
                    eff_counter++;
                }

                // Alternatively, each line contains additional three columns with g1, g2 and weight for TTHG
                if(values.size() == NR_OF_EFF - 1 + 4 * nr_of_incident_angles)
                {
                    // the current line contains enough values
                    rec = true;

                    // Set the dust grain optical properties
                    for(uint i = 0; i < NR_OF_EFF - 1; i++)
                        eff_wl[a * NR_OF_EFF + i].setValue(w, wavelength_list_dustcat[w], values[i]);

                    // For each incident angles, get Qtrq and HG parameters
                    for(uint i_inc = 0; i_inc < nr_of_incident_angles; i_inc++)
                    {
                        // At the first wavelength, resize the splines for the wavelengths
                        if(w == 0)
                        {
                            Qtrq_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g2_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                            HG_g3_factor_wl[a * nr_of_incident_angles + i_inc].resize(nr_of_wavelength_dustcat);
                        }

                        // Set the radiative torque efficiency
                        Qtrq_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w, wavelength_list_dustcat[w], values[NR_OF_EFF - 1 + i_inc]);

                        // Set the Henyey-Greenstein g factor
                        HG_g_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            values[NR_OF_EFF - 1 + i_inc + nr_of_incident_angles]);

                        // Set the second Henyey-Greenstein factor
                        HG_g2_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            values[NR_OF_EFF - 1 + i_inc + 2 * nr_of_incident_angles]);

                        // Set the third Henyey-Greenstein factor
                        HG_g3_factor_wl[a * nr_of_incident_angles + i_inc].setValue(
                            w,
                            wavelength_list_dustcat[w],
                            values[NR_OF_EFF - 1 + i_inc + 3 * nr_of_incident_angles]);
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
                        HG_g2_factor_wl[a * nr_of_incident_angles + i_inc].createSpline();
                        HG_g3_factor_wl[a * nr_of_incident_angles + i_inc].createSpline();
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
        cout << ERROR_LINE << "Wrong amount of efficiencies in file!" << endl;
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
    HGg2 = new double *[nr_of_dust_species];
    HGg3 = new double *[nr_of_dust_species];

    CextMean = new double *[nr_of_dust_species];
    CabsMean = new double *[nr_of_dust_species];
    CscaMean = new double *[nr_of_dust_species];

    // Init splines for incident angle interpolation of Qtrq, HG parameters
    Qtrq = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g_factor = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g2_factor = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g3_factor = new spline[nr_of_dust_species * nr_of_wavelength];

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
        HGg2[a] = new double[nr_of_wavelength];
        HGg3[a] = new double[nr_of_wavelength];

        CextMean[a] = new double [nr_of_wavelength];
        fill(CextMean[a], CextMean[a] + nr_of_wavelength, 0);
        CabsMean[a] = new double [nr_of_wavelength];
        fill(CabsMean[a], CabsMean[a] + nr_of_wavelength, 0);
        CscaMean[a] = new double [nr_of_wavelength];
        fill(CscaMean[a], CscaMean[a] + nr_of_wavelength, 0);

        for(uint w = 0; w < nr_of_wavelength; w++)
        {
            // Resize the splines of Qtrq and HG g factor for each wavelength
            Qtrq[w * nr_of_dust_species + a].resize(nr_of_incident_angles);
            HG_g_factor[w * nr_of_dust_species + a].resize(nr_of_incident_angles);
            HG_g2_factor[w * nr_of_dust_species + a].resize(nr_of_incident_angles);
            HG_g3_factor[w * nr_of_dust_species + a].resize(nr_of_incident_angles);

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
                        HG_g_factor_wl[a * nr_of_incident_angles + i_inc].getValue(wavelength_list[w], CONST));
                    HG_g2_factor[w * nr_of_dust_species + a].setValue(
                        i_inc,
                        i_inc * d_ang,
                        HG_g2_factor_wl[a * nr_of_incident_angles + i_inc].getValue(wavelength_list[w], CONST));
                    HG_g3_factor[w * nr_of_dust_species + a].setValue(
                        i_inc,
                        i_inc * d_ang,
                        HG_g3_factor_wl[a * nr_of_incident_angles + i_inc].getValue(wavelength_list[w], CONST));
                }

                // Calculate the average parameters for Henyey-Greenstein phase function over all angles
                double avg_HG_g_factor = HG_g_factor[w * nr_of_dust_species + a].getAverageY();
                double avg_HG_g2_factor = HG_g2_factor[w * nr_of_dust_species + a].getAverageY();
                double avg_HG_g3_factor = HG_g3_factor[w * nr_of_dust_species + a].getAverageY();

                if(avg_HG_g_factor <= -1.0 || avg_HG_g_factor >= 1.0) {
                    cout << ERROR_LINE << "Henyey-Greenstein g factor is invalid: " << avg_HG_g_factor << endl;
                    return false;
                }

                if(param.getPhaseFunctionID(dust_component_choice) == PH_DHG) {
                    if(avg_HG_g2_factor < 0.0 || avg_HG_g2_factor > 1.0) {
                        cout << ERROR_LINE << "Henyey-Greenstein factor alpha is invalid: " << avg_HG_g2_factor << endl;
                        return false;
                    }
                } else if(param.getPhaseFunctionID(dust_component_choice) == PH_TTHG) {
                    if(avg_HG_g2_factor <= -1.0 || avg_HG_g2_factor >= 1.0) {
                        cout << ERROR_LINE << "Henyey-Greenstein factor g2 is invalid: " << avg_HG_g2_factor << endl;
                        return false;
                    }
                    if(avg_HG_g2_factor * avg_HG_g3_factor > 0.0) {
                        cout << ERROR_LINE << "Henyey-Greenstein g1 and g2 must have different signs." << endl;
                        return false;
                    }
                }

                if(avg_HG_g3_factor < 0.0 || avg_HG_g3_factor > 1.0) {
                    cout << ERROR_LINE << "Henyey-Greenstein weight factor is invalid: " << avg_HG_g3_factor << endl;
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
                HGg2[a][w] = avg_HG_g2_factor;
                HGg3[a][w] = avg_HG_g3_factor;
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
                HGg2[a][w] = 0;
                HGg3[a][w] = 1;
            }

            // Activate the splines of Qtrq and HG g factor
            Qtrq[w * nr_of_dust_species + a].createSpline();
            HG_g_factor[w * nr_of_dust_species + a].createSpline();
            HG_g2_factor[w * nr_of_dust_species + a].createSpline();
            HG_g3_factor[w * nr_of_dust_species + a].createSpline();

            CextMean[a][w] = PI * a_eff_squared[a] * (2.0 * Qext1[a][w] + Qext2[a][w]) / 3.0;
            CabsMean[a][w] = PI * a_eff_squared[a] * (2.0 * Qabs1[a][w] + Qabs2[a][w]) / 3.0;
            CscaMean[a][w] = PI * a_eff_squared[a] * (2.0 * Qsca1[a][w] + Qsca2[a][w]) / 3.0;
        }
    }

    // Read the scattering matrix if MIE scattering should be used
    // With the same grid of wavelengths and grain sizes
    if(param.getPhaseFunctionID(dust_component_choice) == PH_MIE)
        if(!readScatteringMatrices(path, nr_of_wavelength_dustcat, wavelength_list_dustcat))
            return false;

    // Remove temporary pointer arrays
    delete[] eff_wl;
    delete[] Qtrq_wl;
    delete[] HG_g_factor_wl;
    delete[] HG_g2_factor_wl;
    delete[] HG_g3_factor_wl;

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

    // Set number of grain sizes for Mie theory (1 if only one grain size is used)
    if(a_min_mixture == a_max_mixture)
        nr_of_dust_species = 1;
    else
        nr_of_dust_species = MIE_NR_DUST_SIZE;
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
        cout << ERROR_LINE << "Cannot open dust refractive index file:" << endl;
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
                    cout << ERROR_LINE << "Line " << line_counter << " wrong amount of numbers!" << endl;
                    return false;
                }

                // The number of wavelength used by the dust catalog
                nr_of_wavelength_dustcat = (uint)values[0];

                // The number of incident angles
                nr_of_incident_angles = 1; // For non-spherical: (uint) values[1];

                // The aspect ratio of minor to major dust grain axes
                aspect_ratio = 1; // For non-spherical: values[2];

                // The material density (only used if no one was set in the command file)
                if(material_density == 0)
                    material_density = values[3];

                // The sublimation temperature
                sub_temp = values[4];

                // The delta value fot the RAT alignment theory
                delta_rat = 1; // a sphere has delta = 1; for non-spherical: values[5];
                // see e.g. Draine & Weingartner, 1996 ApJ 470:551, 1997 ApJ 480:663

                // The boolean value if the dust catalog is aligned or not
                if(values[6] == 1)
                    is_align = true;
                else
                    is_align = false;

                // Calculate the GOLD alignment g factor
                gold_g_factor = 0.5 * (aspect_ratio * aspect_ratio - 1);

                // Init splines for wavelength interpolation of the dust optical properties
                refractive_index_real.resize(nr_of_wavelength_dustcat);
                refractive_index_imag.resize(nr_of_wavelength_dustcat);

                // Set size parameters
                a_eff = new double[nr_of_dust_species];
                a_eff_squared = new double[nr_of_dust_species];
                grain_distribution_x_aeff_sq = new double[nr_of_dust_species];
                grain_size_distribution = new double[nr_of_dust_species];
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

    if(wavelength_list[0] < refractive_index_real.getX(0) ||
        wavelength_list[nr_of_wavelength - 1] > refractive_index_real.getX(nr_of_wavelength_dustcat - 1))
    {
        cout << WARNING_LINE << "The wavelength range is out of the limits of the catalog. This may cause problems!\n"
            << "         wavelength range          : " << wavelength_list[0] << " [m] to "
            << wavelength_list[nr_of_wavelength - 1] << " [m]\n"
            << "         wavelength range (catalog): " << refractive_index_real.getX(0) << " [m] to "
            << refractive_index_real.getX(nr_of_wavelength_dustcat - 1) << " [m]" << endl;
        if(!IGNORE_WAVELENGTH_RANGE)
        {
            cout << "         To continue, set 'IGNORE_WAVELENGTH_RANGE' to 'true' in src/Typedefs.h and recompile!" << endl;
            return false;
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
        cout << ERROR_LINE << "Wrong amount of efficiencies in file!" << endl;
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
    HGg2 = new double *[nr_of_dust_species];
    HGg3 = new double *[nr_of_dust_species];

    CextMean = new double *[nr_of_dust_species];
    CabsMean = new double *[nr_of_dust_species];
    CscaMean = new double *[nr_of_dust_species];

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
        HGg2[a] = new double[nr_of_wavelength];
        HGg3[a] = new double[nr_of_wavelength];

        CextMean[a] = new double [nr_of_wavelength];
        fill(CextMean[a], CextMean[a] + nr_of_wavelength, 0);
        CabsMean[a] = new double [nr_of_wavelength];
        fill(CabsMean[a], CabsMean[a] + nr_of_wavelength, 0);
        CscaMean[a] = new double [nr_of_wavelength];
        fill(CscaMean[a], CscaMean[a] + nr_of_wavelength, 0);
    }

    // Init splines for incident angle interpolation of Qtrq and parameters for Henyey-Greenstein phase function
    Qtrq = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g_factor = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g2_factor = new spline[nr_of_dust_species * nr_of_wavelength];
    HG_g3_factor = new spline[nr_of_dust_species * nr_of_wavelength];

    // Set variables for scattering via Mie theory
    nr_of_scat_mat_elements = 4;

    // --- Phi angle
    nr_of_scat_phi = 1;

    // --- Theta angle
    uint nr_of_scat_theta_start = 2 * NANG - 1;

    // Init normal scattering matrix array
    initNrOfScatThetaArray();
    initScatThetaArray();
    initScatteringMatrixArray();

    // Init counter and percentage to show progress
    ullong per_counter = 0;
    float last_percentage = 0;

    // Init error check
    bool error = false;

    double max_rel_diff = 0.0;

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
            HG_g2_factor[w * nr_of_dust_species + a].resize(nr_of_incident_angles);
            HG_g3_factor[w * nr_of_dust_species + a].resize(nr_of_incident_angles);

            if(sizeIndexUsed(a))
            {
                // Init variables and pointer arrays
                double *S11_start, *S12_start, *S33_start, *S34_start;
                S11_start = new double[nr_of_scat_theta_start];
                S12_start = new double[nr_of_scat_theta_start];
                S33_start = new double[nr_of_scat_theta_start];
                S34_start = new double[nr_of_scat_theta_start];

                dlist scat_angle_start(nr_of_scat_theta_start);
                for(uint i_scat_ang=0; i_scat_ang < nr_of_scat_theta_start; i_scat_ang++)
                {
                    scat_angle_start[i_scat_ang] = i_scat_ang * PI/(nr_of_scat_theta_start-1);
                }

                // Set size index and refractive index as complex number
                double x = 2.0 * PI * a_eff[a] / wavelength_list[w];
#if BENCHMARK == PINTE
                dcomplex refractive_index = dcomplex(refractive_index_real.getValue(wavelength_list[w], LOGLINEAR),
                                                     refractive_index_imag.getValue(wavelength_list[w], LOGLINEAR));
#else
                dcomplex refractive_index = dcomplex(refractive_index_real.getValue(wavelength_list[w], LOGLINEAR),
                                                     refractive_index_imag.getValue(wavelength_list[w], LOGLINEAR));
#endif

                // Calculate Mie-scattering
                if(!CMathFunctions::calcWVMie(x,
                                              scat_angle_start,
                                              refractive_index,
                                              Qext1[a][w],
                                              Qabs1[a][w],
                                              Qsca1[a][w],
                                              HGg[a][w],
                                              S11_start,
                                              S12_start,
                                              S33_start,
                                              S34_start))
                    error = true;

                dlist S11_final(1), S12_final(1), S33_final(1), S34_final(1), scat_angle_final(1);
                S11_final[0] = S11_start[0];
                S12_final[0] = S12_start[0];
                S33_final[0] = S33_start[0];
                S34_final[0] = S34_start[0];
                scat_angle_final[0] = scat_angle_start[0];

                // set default HG values if optical properties are calcualted with Mie-scattering
                HGg2[a][w] = 0.0;
                HGg3[a][w] = 1.0;

                double current_S11_rel_diff;

                for(uint i_scat_ang=0; i_scat_ang < nr_of_scat_theta_start -1; i_scat_ang++)
                {
                    dlist S11_tmp(1), S12_tmp(1), S33_tmp(1), S34_tmp(1);
                    dlist scat_angle_tmp(2);

                    S11_tmp[0] = S11_start[i_scat_ang+1];
                    S12_tmp[0] = S12_start[i_scat_ang+1];
                    S33_tmp[0] = S33_start[i_scat_ang+1];
                    S34_tmp[0] = S34_start[i_scat_ang+1];
                    scat_angle_tmp[0] = scat_angle_start[i_scat_ang+1];
                    scat_angle_tmp[1] = 0.5 * (scat_angle_final.back() + scat_angle_start[i_scat_ang+1]);

                    while(true)
                    {
                        current_S11_rel_diff = abs( S11_tmp.back() - S11_final.back() ) / max( S11_tmp.back(), S11_final.back() );

                        // Subdivide scattering angles only if size parameter x is not too large.
                        // Large x lead to extrem forward scattering and the subdivision criterion will
                        // get almost impossible to achive for small scattering angles.
                        // The limit of x=100 is somewhat arbitrary.
                        while(x < 100.0 && current_S11_rel_diff > MAX_MIE_SCA_REL_DIFF)
                        {
                            double *pointer_s11_tmp, *pointer_s12_tmp, *pointer_s33_tmp, *pointer_s34_tmp;
                            pointer_s11_tmp = new double[1];
                            pointer_s12_tmp = new double[1];
                            pointer_s33_tmp = new double[1];
                            pointer_s34_tmp = new double[1];

                            dlist scat_angle_calc(1);
                            scat_angle_calc[0] = scat_angle_tmp.back();

                            if(!CMathFunctions::calcWVMie(x,
                                                          scat_angle_calc,
                                                          refractive_index,
                                                          Qext1[a][w],
                                                          Qabs1[a][w],
                                                          Qsca1[a][w],
                                                          HGg[a][w],
                                                          pointer_s11_tmp,
                                                          pointer_s12_tmp,
                                                          pointer_s33_tmp,
                                                          pointer_s34_tmp))
                                error = true;

                            S11_tmp.push_back(pointer_s11_tmp[0]);
                            S12_tmp.push_back(pointer_s12_tmp[0]);
                            S33_tmp.push_back(pointer_s33_tmp[0]);
                            S34_tmp.push_back(pointer_s34_tmp[0]);
                            scat_angle_tmp.push_back( 0.5 * (scat_angle_final.back() + scat_angle_tmp.back()) );

                            delete[] pointer_s11_tmp;
                            delete[] pointer_s12_tmp;
                            delete[] pointer_s33_tmp;
                            delete[] pointer_s34_tmp;

                            current_S11_rel_diff = abs( S11_tmp.back() - S11_final.back() ) / max( S11_tmp.back(), S11_final.back() );
                        }
                        scat_angle_tmp.pop_back();

                        S11_final.push_back(S11_tmp.back());
                        S12_final.push_back(S12_tmp.back());
                        S33_final.push_back(S33_tmp.back());
                        S34_final.push_back(S34_tmp.back());
                        scat_angle_final.push_back(scat_angle_tmp.back());

                        S11_tmp.pop_back();
                        S12_tmp.pop_back();
                        S33_tmp.pop_back();
                        S34_tmp.pop_back();
                        scat_angle_tmp.pop_back();

                        if(S11_tmp.size() < 1)
                            break;
                        else
                            scat_angle_tmp.push_back( 0.5 * (scat_angle_final.back() + scat_angle_tmp.back()) );
                    }
                    S11_tmp.clear();
                    S12_tmp.clear();
                    S33_tmp.clear();
                    S34_tmp.clear();
                    scat_angle_tmp.clear();
                }
                scat_angle_start.clear();
                delete[] S11_start;
                delete[] S12_start;
                delete[] S33_start;
                delete[] S34_start;

                uint nr_of_scat_theta_final = scat_angle_final.size();

                // Set missing Efficiencies for other axis
                Qext2[a][w] = Qext1[a][w];
                Qabs2[a][w] = Qabs1[a][w];
                Qsca2[a][w] = Qsca1[a][w];
                Qcirc[a][w] = 0;

                if(scat_theta[a][w] != 0){
                    delete[] scat_theta[a][w];
                }

                double diff_tmp;
                for(uint sth = 1; sth < nr_of_scat_theta_final; sth++)
                {
                    diff_tmp = abs(S11_final[sth-1] - S11_final[sth]) / max(S11_final[sth-1], S11_final[sth]);
                    max_rel_diff = max(diff_tmp, max_rel_diff);
                }

                scat_theta[a][w] = new double[nr_of_scat_theta_final];

                for(uint inc = 0; inc < nr_of_incident_angles; inc++)
                {
                    sca_mat[a][w][inc] = new Matrix2D*[nr_of_scat_phi];
                    for(uint sph = 0; sph < nr_of_scat_phi; sph++)
                    {
                        sca_mat[a][w][inc][sph] = new Matrix2D[nr_of_scat_theta_final];
                        for(uint sth = 0; sth < nr_of_scat_theta_final; sth++)
                        {
                            sca_mat[a][w][inc][sph][sth].resize(4, 4);

                            sca_mat[a][w][inc][sph][sth](0, 0) = S11_final[sth]; // S11
                            sca_mat[a][w][inc][sph][sth](1, 1) = S11_final[sth]; // S22

                            sca_mat[a][w][inc][sph][sth](0, 1) = S12_final[sth]; // S12
                            sca_mat[a][w][inc][sph][sth](1, 0) = S12_final[sth]; // S21

                            sca_mat[a][w][inc][sph][sth](2, 2) = S33_final[sth]; // S33
                            sca_mat[a][w][inc][sph][sth](3, 3) = S33_final[sth]; // S44

                            sca_mat[a][w][inc][sph][sth](2, 3) = S34_final[sth]; // S34
                            sca_mat[a][w][inc][sph][sth](3, 2) = -S34_final[sth]; // S43

                            scat_theta[a][w][sth] = scat_angle_final[sth];
                        }
                    }
                }

                S11_final.clear();
                S12_final.clear();
                S33_final.clear();
                S34_final.clear();

                nr_of_scat_theta[a][w] = nr_of_scat_theta_final;
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
                HGg2[a][w] = 0;
                HGg3[a][w] = 1;
            }

            // Activate the splines of Qtrq and HG g factor
            Qtrq[w * nr_of_dust_species + a].createSpline();
            HG_g_factor[w * nr_of_dust_species + a].createSpline();
            HG_g2_factor[w * nr_of_dust_species + a].createSpline();
            HG_g3_factor[w * nr_of_dust_species + a].createSpline();

            CextMean[a][w] = PI * a_eff_squared[a] * (2.0 * Qext1[a][w] + Qext2[a][w]) / 3.0;
            CabsMean[a][w] = PI * a_eff_squared[a] * (2.0 * Qabs1[a][w] + Qabs2[a][w]) / 3.0;
            CscaMean[a][w] = PI * a_eff_squared[a] * (2.0 * Qsca1[a][w] + Qsca2[a][w]) / 3.0;
        } // end of wavelength loop
    } // end of grain size loop

    // Set that the scattering matrix was successfully read
    scat_loaded = true;

    if(error)
    {
        cout << ERROR_LINE << "Problem with optical properties calculation" << endl;
        return false;
    }

    printIDs();
    cout << "- calculating optical properties: done          " << endl;

    if(max_rel_diff > 0.5) // arbitrary limit
    {
        cout << WARNING_LINE << "Number of scattering angles might be too low (max rel diff = " << max_rel_diff << ")." << endl;
        cout << "  If required, increase 'NANG' or decrease 'MAX_MIE_SCA_REL_DIFF' (for x < 100) in src/Typedefs.h and recompile!" << endl;
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
        cout << ERROR_LINE << "Cannot open scattering matrix info file:" << endl;
        cout << inf_filename << endl;
        return false;
    }

    // Init progress counter
    uint line_counter = 0;
    uint cmd_counter = 0;

    uint nr_of_scat_theta_tmp = 2*NANG - 1;

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
                    cout << ERROR_LINE << "Wrong amount of dust component parameters in:" << endl;
                    cout << inf_filename.c_str() << " line " << line_counter << "!" << endl;
                    return false;
                }

                // The number of dust grain sizes
                if(values[0] != nr_of_dust_species)
                {
                    cout << ERROR_LINE << inf_filename.c_str() << " line " << line_counter << "!" << endl;
                    cout << "Number of dust species does not match the number in the "
                            "dust parameters file!"
                         << endl;
                    return false;
                }

                // The number of wavelength used by the dust catalog
                if(values[1] != nr_of_wavelength_dustcat)
                {
                    cout << ERROR_LINE << inf_filename.c_str() << " line " << line_counter << "!" << endl;
                    cout << "Number of wavelength does not match the number in the dust "
                            "parameters file!"
                         << endl;
                    return false;
                }

                // The number of incident angles
                if(values[2] != nr_of_incident_angles)
                {
                    cout << ERROR_LINE << inf_filename.c_str() << " line " << line_counter << "!" << endl;
                    cout << "Number of incident angles does not match the number in the "
                            "dust parameters file!"
                         << endl;
                    return false;
                }

                // The number of phi angles (outgoing radiation)
                nr_of_scat_phi = uint(values[3]);

                // The number of theta angles (outgoing radiation)
                nr_of_scat_theta_tmp = uint(values[4]);

                break;

            case 2:
                // The second line needs one value
                if(values.size() != 1)
                {
                    cout << ERROR_LINE << "Wrong amount of scattering matrix elements in:" << endl;
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
                    cout << ERROR_LINE << "Wrong amount of matrix elements in:" << endl;
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
        for(uint inc = 0; inc < nr_of_incident_angles; inc++)
        {
            // Third init of the scattering matrix interp
            sca_mat_wl[a][inc] = new interp **[nr_of_scat_phi];
            for(uint sph = 0; sph < nr_of_scat_phi; sph++)
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
                sca_mat_wl[a][inc][sph] = new interp *[nr_of_scat_theta_tmp];
                for(uint sth = 0; sth < nr_of_scat_theta_tmp; sth++)
                {
                    // Fifth init of the scattering matrix interp
                    sca_mat_wl[a][inc][sph][sth] = new interp[nr_of_scat_mat_elements];

                    // Resize the linear for each scattering matrix element
                    for(uint mat = 0; mat < nr_of_scat_mat_elements; mat++)
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
            cout << ERROR_LINE << "Cannot open scattering matrix file:" << endl;
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
                    for(uint sth = 0; sth < nr_of_scat_theta_tmp; sth++)
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
    initScatteringMatrixArray(nr_of_scat_theta_tmp);

    // Fill values of the scattering matrix via interpolation
    #pragma omp parallel for
    for(int a = 0; a < int(nr_of_dust_species); a++)
    {
        if(sizeIndexUsed(a))
        {
            for(uint inc = 0; inc < nr_of_incident_angles; inc++)
            {
                for(uint sph = 0; sph < nr_of_scat_phi; sph++)
                {
                    for(uint sth = 0; sth < nr_of_scat_theta_tmp; sth++)
                    {
                        for(uint e = 0; e < 16; e++)
                        {
                            double sign = CMathFunctions::sgn(elements[e]);
                            int pos = abs(elements[e]);
                            if(pos > 0)
                            {
                                for(uint w = 0; w < nr_of_wavelength; w++)
                                {
                                    sca_mat[a][w][inc][sph][sth](e) =
                                        sign *
                                        sca_mat_wl[a][inc][sph][sth][pos - 1].getValue(wavelength_list[w]);
                                    scat_theta[a][w][sth] = PI * double(sth) / double(nr_of_scat_theta_tmp - 1);
                                }
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
        cout << ERROR_LINE << "Cannot open calorimetry file:" << endl;
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
                    cout << ERROR_LINE << "Wrong amount of calorimetry temperatures in:" << endl;
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
                    cout << ERROR_LINE << "Wrong calorimetry temperatures in:" << endl;
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
                    cout << ERROR_LINE << "Wrong calorimetry type in:" << endl;
                    cout << calo_filename.c_str() << " line " << line_counter << "!" << endl;
                    return false;
                }

                // The unit of the calorimetry data
                calorimetry_type = uint(values[0]);

                // Only heat capacity or enthalpy are possible
                if(calorimetry_type != CALO_HEAT_CAP && calorimetry_type != CALO_ENTHALPY)
                {
                    cout << ERROR_LINE << "Wrong calorimetry type in:" << endl;
                    cout << calo_filename.c_str() << " line " << line_counter << "!" << endl;
                    return false;
                }
                break;

            default:
                // The other lines have either one or a value per dust grain size
                if(values.size() != 1 && values.size() != nr_of_dust_species)
                {
                    cout << ERROR_LINE << "Wrong amount of dust species in:" << endl;
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

bool CDustComponent::writeComponentData(string path_data)
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
    string str_title, plot_title;

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
#if BENCHMARK == PINTE
        string path_mueller = path_data + "dust_mixture_" + str_mix_ID_end + "_mueller.dat";

        ofstream mueller_matrix_file(path_mueller.c_str());

        // Error message if the write does not work
        if(mueller_matrix_file.fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n" << path_cross << endl;
            return false;
        }

        for(uint sth = 0; sth < nr_of_scat_theta[0][0]; sth++)
        {
            mueller_matrix_file << scat_theta[0][0][sth] << TAB;
            mueller_matrix_file << sca_mat[0][0][0][0][sth](0, 0) << TAB;
            mueller_matrix_file << sca_mat[0][0][0][0][sth](0, 1) << TAB;
            mueller_matrix_file << sca_mat[0][0][0][0][sth](2, 2) << TAB;
            mueller_matrix_file << sca_mat[0][0][0][0][sth](2, 3) << endl;
        }
        mueller_matrix_file.close();
#endif

        // For the dust mixture, set the different strings
        str_title = "#Dust mixture\n " + stringID;
        plot_title = "#Dust mixture ";
        plot_title += str_mix_ID_end;
        plot_title += "\n";
        path_data = path_data + "dust_mixture_" + str_mix_ID_end + ".dat";

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
        while(plot_title.find("\n") != string::npos)
        {
            pos = uint(plot_title.find("\n"));
            plot_title.replace(pos, 1, "; ");
        }

        pos = 0;
        while(plot_title.find("#") != string::npos)
        {
            pos = uint(plot_title.find("#"));
            plot_title.replace(pos, 1, " ");
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
        plot_title = str_title;
        path_data = path_data + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + ".dat";
    }

    // Init text file writer for dust grain parameters
    ofstream data_writer(path_data.c_str());

    // Error message if the write does not work
    if(data_writer.fail())
    {
        cout << ERROR_LINE << "Cannot write to:\n" << path_data << endl;
        return false;
    }

    // Add plot commands to file
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
    data_writer << "#density [kg/m^3]" << endl;
    data_writer << getMaterialDensity() << endl;
    data_writer << "#gold g factor" << endl;
    data_writer << gold_g_factor << endl;
    data_writer << "#f_highJ" << endl;
    data_writer << f_highJ << endl;
    data_writer << "#f_correlation" << endl;
    data_writer << f_cor << endl;
    data_writer << "#min. grain size [m]\tmax. grain size [m]" << endl;
    data_writer << a_min_global << "\t" << a_max_global << endl;
    data_writer << "#min. radius [m]\tmax. radius [m]" << endl;
    data_writer << getSizeMin() << "\t" << getSizeMax() << endl;

    // Add data header to file
    data_writer << "#wavelength [m]\tavgCext1 [m^2]\tavgCext2 [m^2]\tavgCabs1 [m^2]\tavgCabs2 [m^2]\tavgCsca1 [m^2]\tavgCsca2 [m^2]\tCcirc [m^2]\tavgHGg\t"
                   "avgQext1\tavgQext2\tavgQabs1\tavgQabs2\tavgQsca1\tavgQsca2\t"
                   "avgKext1 [m^2/kg]\tavgKext2 [m^2/kg]\tavgKabs1 [m^2/kg]\tavgKabs2 [m^2/kg]\tavgKsca1 [m^2/kg]\tavgKsca2 [m^2/kg]"
                << endl;

    // Add cross-sections to file
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        data_writer << wavelength_list[i] << "\t" << getCext1(i) << "\t" << getCext2(i) << "\t" << getCabs1(i)
                    << "\t" << getCabs2(i) << "\t" << getCsca1(i) << "\t" << getCsca2(i) << "\t"
                    << getCcirc(i) << "\t" << getHGg(i) << "\t" << getQext1(i) << "\t" << getQext2(i) << "\t"
                    << getQabs1(i) << "\t" << getQabs2(i) << "\t" << getQsca1(i) << "\t" << getQsca2(i)
                    << "\t" << getKappaExt1(i) << "\t" << getKappaExt2(i) << "\t" << getKappaAbs1(i) << "\t"
                    << getKappaAbs2(i) << "\t" << getKappaSca1(i) << "\t" << getKappaSca2(i) << endl;

    // Close text file writer
    data_writer.close();

    return true;
}

bool CDustComponent::writeComponentPlot(string path_plot)
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
    string path_scat, path_size_dist, str_title, plot_title;
    string plot_sign = "points";

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
#if BENCHMARK == PINTE
        string path_mueller = path_data + "dust_mixture_" + str_mix_ID_end + "_mueller.dat";

        ofstream mueller_matrix_file(path_mueller.c_str());

        // Error message if the write does not work
        if(mueller_matrix_file.fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n" << path_cross << endl;
            return false;
        }

        for(uint sth = 0; sth < nr_of_scat_theta[0][0]; sth++)
        {
            mueller_matrix_file << scat_theta[0][0][sth] << TAB;
            mueller_matrix_file << sca_mat[0][0][0][0][sth](0, 0) << TAB;
            mueller_matrix_file << sca_mat[0][0][0][0][sth](0, 1) << TAB;
            mueller_matrix_file << sca_mat[0][0][0][0][sth](2, 2) << TAB;
            mueller_matrix_file << sca_mat[0][0][0][0][sth](2, 3) << endl;
        }
        mueller_matrix_file.close();
#endif

        // For the dust mixture, set the different strings
        str_title = "#Dust mixture\n " + stringID;
        plot_title = "#Dust mixture ";
        plot_title += str_mix_ID_end;
        plot_title += "\n";
        path_cross = path_plot + "dust_mixture_" + str_mix_ID_end + "_cross.plt";
        path_eff = path_plot + "dust_mixture_" + str_mix_ID_end + "_eff.plt";
        path_kappa = path_plot + "dust_mixture_" + str_mix_ID_end + "_kappa.plt";
        path_diff = path_plot + "dust_mixture_" + str_mix_ID_end + "_diff.plt";
        path_g = path_plot + "dust_mixture_" + str_mix_ID_end + "_g.plt";
        path_scat = path_plot + "dust_mixture_" + str_mix_ID_end + "_scat.plt";
        path_size_dist = path_plot + "dust_mixture_" + str_mix_ID_end + "_size_dist.plt";

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
        while(plot_title.find("\n") != string::npos)
        {
            pos = uint(plot_title.find("\n"));
            plot_title.replace(pos, 1, "; ");
        }

        pos = 0;
        while(plot_title.find("#") != string::npos)
        {
            pos = uint(plot_title.find("#"));
            plot_title.replace(pos, 1, " ");
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
        plot_title = str_title;
        path_cross = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_cross.plt";
        path_eff = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_eff.plt";
        path_kappa = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_kappa.plt";
        path_diff = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_diff.plt";
        path_g = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_g.plt";
        path_scat = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_scat.plt";
        path_size_dist = path_plot + "dust_mixture_" + str_mix_ID_end + "_comp_" + str_comp_ID_end + "_size_dist.plt";
    }

    // Init text file writer for cross-sections
    ofstream cross_writer(path_cross.c_str());

    // Error message if the write does not work
    if(cross_writer.fail())
    {
        cout << ERROR_LINE << "Cannot write to:\n" << path_cross << endl;
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

    // Add plot commands to file
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
    cross_writer << "set title \"" << plot_title << "\"" << endl;
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
        cout << ERROR_LINE << "Cannot write to:\n" << path_eff << endl;
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

    // Add plot commands to file
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
    eff_writer << "set title \"" << plot_title << "\"" << endl;
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
        cout << ERROR_LINE << "Cannot write to:\n" << path_eff << endl;
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

    // Add plot commands to file
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
    kappa_writer << "set title \"" << plot_title << "\"" << endl;
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
        cout << ERROR_LINE << "Cannot write to:\n" << path_cross << endl;
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

    // Add plot commands to file
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
    diff_writer << "set title \"" << plot_title << "\"" << endl;
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

    // Init text file writer for Henyey-Greenstein g factor
    ofstream g_writer(path_g.c_str());

    // Error message if the write does not work
    if(g_writer.fail())
    {
        cout << ERROR_LINE << "Cannot write to:\n" << path_g << endl;
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

    // Add plot commands to file
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
    g_writer << "set title \"" << plot_title << "\"" << endl;
    g_writer << "plot \'-\' with " << plot_sign << " lc rgb \"#0000F0\"" << endl;

    // Add HG g factor to file (if larger than 0)
    for(uint i = wavelength_offset; i < nr_of_wavelength; i++)
        if(getHGg(i) > 0)
            g_writer << wavelength_list[i] << "\t" << getHGg(i) << endl;
    g_writer << "e" << endl;

    // Close text file writer
    g_writer.close();

    // ------------------------------------------------------

    if(MAX_MIE_SCA_REL_DIFF >= 1 && phID == PH_MIE && is_mixture)
    {
        // Init text file writer for scattering matrix
        ofstream scat_writer(path_scat.c_str());

        uint nr_of_scat_theta_tmp = 2 * NANG - 1;

        // Error message if the write does not work
        if(scat_writer.fail())
        {
            cout << ERROR_LINE << "Cannot write to " << path_scat << endl;
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
            S11[wID] = new double[nr_of_scat_theta_tmp];
            S12[wID] = new double[nr_of_scat_theta_tmp];
            for(uint sth = 0; sth < nr_of_scat_theta_tmp; sth++)
            {
                // Init and reset variables
                double sum = 0;

                double * S11_tmp = new double[nr_of_dust_species];
                double * S12_tmp = new double[nr_of_dust_species];
                for(uint a = 0; a < nr_of_dust_species; a++)
                {
                    if(sizeIndexUsed(a) && nr_of_scat_theta[a][w] != 0)
                    {
                        double Csca_tmp = getCscaMean(a, w);
                        sum += Csca_tmp;
                        double rel_amount = grain_size_distribution[a] / weight;
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

                delete[] S11_tmp;
                delete[] S12_tmp;
            }
        }

        // Add plot commands to file
        scat_writer << "reset" << endl;
        scat_writer << "set grid" << endl;
        scat_writer << "set multiplot layout 2,1 rowsfirst" << endl;

        if(nr_of_wavelength > 1)
            scat_writer << "set xrange[" << 0 << ":" << 180 << "]" << endl;
        scat_writer << "set yrange[" << S11min << ":" << S11max << "]" << endl;
        scat_writer << "set format x \"%.1f\"" << endl;
        scat_writer << "set format y \"%.1te%02T\"" << endl;
        scat_writer << "set ylabel \'S11\'" << endl;
        scat_writer << "set xlabel \'{\u03B8} [°]\'" << endl;
        scat_writer << "set title \"" << plot_title << "\"" << endl;
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

            for(uint sth = 0; sth < nr_of_scat_theta_tmp; sth++)
                scat_writer << 180 * sth / (nr_of_scat_theta_tmp - 1) << "\t" << S11[wID][sth] << endl;
            scat_writer << "e" << endl;
        }

        if(nr_of_wavelength > 1)
            scat_writer << "set xrange[" << 0 << ":" << 180 << "]" << endl;
        scat_writer << "set yrange[" << S12min << ":" << S12max << "]" << endl;
        scat_writer << "set format x \"%.1f\"" << endl;
        scat_writer << "set format y \"%.1te%02T\"" << endl;
        scat_writer << "set ylabel \'S12\'" << endl;
        scat_writer << "set xlabel \'{\u03B8} [°]\'" << endl;
        scat_writer << "set title \"" << plot_title << "\"" << endl;
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

            for(uint sth = 0; sth < nr_of_scat_theta_tmp; sth++)
                scat_writer << 180 * sth / (nr_of_scat_theta_tmp - 1) << "\t" << S12[wID][sth] << endl;
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

    // ------------------------------------------------------

    // Init text file writer for size distribution
    ofstream size_dist_writer(path_size_dist.c_str());

    // Error message if the write does not work
    if(size_dist_writer.fail())
    {
        cout << ERROR_LINE << "Cannot write to:\n" << path_size_dist << endl;
        return false;
    }

    // Get weight
    double weight = getWeight();

    // Init plot limits
    double SizeDistMin = 1e100, SizeDistMax = -1e100;

    // Find min and max values over the wavelength
    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a))
        {
            if(getGrainSizeDistribution(a) / weight < SizeDistMin && getGrainSizeDistribution(a) / weight > 0)
                SizeDistMin = getGrainSizeDistribution(a) / weight;
            if(getGrainSizeDistribution(a) / weight > SizeDistMax)
                SizeDistMax = getGrainSizeDistribution(a) / weight;
        }
    }

    // Add a bit more space for good visualization
    SizeDistMin *= 0.9;
    SizeDistMax *= 1.10;

    // Add plot commands to file
    size_dist_writer << "reset" << endl;
    if(nr_of_dust_species > 1)
        size_dist_writer << "set log x" << endl;
    size_dist_writer << "set log y" << endl;
    size_dist_writer << "set grid" << endl;
    size_dist_writer << "unset key" << endl;
    if(nr_of_dust_species > 1)
        size_dist_writer << "set xrange[" << getSizeMin() << ":" << getSizeMax() << "]" << endl;
    size_dist_writer << "set yrange[" << SizeDistMin << ":" << SizeDistMax << "]" << endl;
    size_dist_writer << "set format x \"%.1te%02T\"" << endl;
    size_dist_writer << "set format y \"%.1te%02T\"" << endl;
    size_dist_writer << "set ylabel \'{/Symbol d}N(a)'" << endl;
    size_dist_writer << "set xlabel \'a [m]\'" << endl;
    size_dist_writer << "set title \"" << plot_title << "\"" << endl;
    size_dist_writer << "plot \'-\' with " << plot_sign << " lc rgb \"#0000F0\"" << endl;

    // Add HG g factor to file (if larger than 0)
    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a))
            size_dist_writer << getEffectiveRadius(a) << "\t" << getGrainSizeDistribution(a) / weight << endl;
        else
            size_dist_writer << getEffectiveRadius(a) << "\t0" << endl;
    }
    size_dist_writer << "e" << endl;

    // Close text file writer
    size_dist_writer.close();

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
    tHGg2 = new double[nr_of_wavelength];
    tHGg3 = new double[nr_of_wavelength];

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
        tHGg2[w] = getHGg2(w);
        tHGg3[w] = getHGg3(w);
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
                tmpQB[w] = meanCabs * CMathFunctions::planck(wavelength_list[w], tmp_temp);
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
                amount[a] = grain_size_distribution[a];
                amount_abs[a] = grain_size_distribution[a] * getCabsMean(a, w);
                amount_sca[a] = grain_size_distribution[a] * getCscaMean(a, w);
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

    // Resize wavelength list for stochastic heating propabilities
    dlist wl_list(WL_STEPS);
    CMathFunctions::LogList(WL_MIN, WL_MAX, wl_list, 10);

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
        double * tmpQB = new double[WL_STEPS];

        // Set each entry of tab_em with integrated value of the Planck function times
        // the absorption cross-section
        for(uint t = 0; t < nr_of_temperatures; t++)
        {
            // Get temperature from tab_temp spline
            double tmp_temp = tab_temp.getValue(t);

            // Calculate absorption cross-section times Planck function for each
            // wavelength
            for(uint w = 0; w < WL_STEPS; w++)
            {
                uint wID = getWavelengthID(wl_list[w]);
                // the wavelength axis of Cabs may be greater than WL_STEPS if detectors are defined
                tmpQB[w] = getCabsMean(a, wID) * CMathFunctions::planck(wl_list[w], tmp_temp);
            }

            // Calculate QB integrated over all wavelengths
            double tt = CMathFunctions::integ(wl_list, tmpQB, 0, WL_STEPS - 1);

            // Set tab_em spline with the integrated value
            tab_em[a].setValue(t, tt, tmp_temp);
            if(calorimetry_loaded)
                tab_em_inv[a].setValue(t, tmp_temp, tt);
        }

        // Create spline for absorption/emission rates
        tab_em[a].createSpline();
        if(calorimetry_loaded)
            tab_em_inv[a].createSpline();

        // Increase progress counter
        per_counter++;

        // Delete pointer array
        delete[] tmpQB;
    }
}

void CDustComponent::preCalcMieScatteringProb()
{
    // Init arrays of splines/interp
    avg_scattering_frac = new interp *[nr_of_dust_species];
    phase_pdf = new interp *[nr_of_dust_species];

    // Init counter for progress
    ullong per_counter = 0;

    #pragma omp parallel for
    for(int a = 0; a < int(nr_of_dust_species); a++)
    {
        // Init arrays of interp
        avg_scattering_frac[a] = new interp[nr_of_wavelength];
        phase_pdf[a] = new interp[nr_of_wavelength];

        if(sizeIndexUsed(a))
        {
            for(uint w = 0; w < nr_of_wavelength; w++)
            {
                // Init pointer arrays
                double * S11_tmp = new double[nr_of_scat_theta[a][w]];
                double * S11_solid_angle = new double[nr_of_scat_theta[a][w]];
                double * tmp_scat_frac = new double[nr_of_scat_theta[a][w]];

                // Resize splines for scattering angles
                avg_scattering_frac[a][w].resize(nr_of_scat_theta[a][w]);
                phase_pdf[a][w].resize(nr_of_scat_theta[a][w]);

                for(uint sth = 0; sth < nr_of_scat_theta[a][w]; sth++)
                {
                    // Calculate the scattering propability in a certain direction
                    S11_tmp[sth] = double(sca_mat[a][w][0][0][sth](0, 0));

                    // Calculate the modified angle for integration
                    S11_solid_angle[sth] = PIx2 * cos(scat_theta[a][w][sth]);

                    // Set the current value for the integration array tmp_scat_frac
                    if(sth == 0)
                        tmp_scat_frac[sth] = 0;
                    else
                        tmp_scat_frac[sth] = -CMathFunctions::integ(S11_solid_angle, S11_tmp, 0, sth);
                }

                // Integral of the scattering S11 value over the full sphere
                double int_scat_frac = tmp_scat_frac[nr_of_scat_theta[a][w] - 1];

                for(uint sth = 0; sth < nr_of_scat_theta[a][w]; sth++)
                {
                    if(int_scat_frac > 0)
                    {
                        // Set the cumulative distribution function of being scattered at a
                        // certain angle
                        avg_scattering_frac[a][w].setValue(
                            sth, tmp_scat_frac[sth] / int_scat_frac, scat_theta[a][w][sth]);

                        // Set the phase function of how much is scattered at a certain angle
                        phase_pdf[a][w].setValue(sth, scat_theta[a][w][sth], S11_tmp[sth] / int_scat_frac);
                    }
                }

                // Activate spline for the cumulative distribution function
                //avg_scattering_frac[a][w].createSpline();

                delete[] S11_tmp;
                delete[] S11_solid_angle;
                delete[] tmp_scat_frac;
            }

            // Increase progress counter
            per_counter++;
        }
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
            pl_mean[a] = new double[nr_of_wavelength];

        // Get temperature from tab_temp spline
        double temp = tab_temp.getValue(uint(t));

        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            // Resize avg_planck_frac to number of wavelengths
            avg_planck_frac[t * nr_of_dust_species + a].resize(nr_of_wavelength);
        }

        // Calculate dPlanck/dT times mean absorption cross-section for each wavelength
        for(uint w = 0; w < nr_of_wavelength; w++)
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

        // Increase progress counter
        per_counter++;
    }
}

void CDustComponent::preCalcTemperatureLists(double minTemp, double maxTemp, uint nr_of_temperatures)
{
    // Init variables
    double tmp_temp;

    // Init spline for tabulated Planck function values
    tab_planck = new spline[nr_of_wavelength];

    // Resize tabulated temperature spline
    tab_temp.resize(nr_of_temperatures);

    // Set each entry of tab_temp with the corresponding temperature (log distribution)
    for(uint t = 0; t < nr_of_temperatures; t++)
    {
        // Calculate the temperature of a certain index
        tmp_temp = minTemp * pow(maxTemp / minTemp, double(t) / double(nr_of_temperatures - 1));

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
        // printIDs();
        // cout << "- pre-calculation of Planck functions: " << 100.0 * float(w) / float(nr_of_wavelength)
        //      << "                                \r";

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
    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        a_eff[a] = values[a];
        a_eff_squared[a] = a_eff[a] * a_eff[a];
    }

    // Calculate the dust grain size distribution depending on user input
    if(size_keyword.find("plaw") != std::string::npos)
    {
        // Power-law distribution
        for(uint a = 0; a < nr_of_dust_species; a++)
            grain_size_distribution[a] = pow(a_eff[a], size_parameter[0]);

        // Add exponential decay if demanded
        if(size_keyword.find("-ed") != std::string::npos)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                if(a_eff[a] > size_parameter[1])
                    grain_size_distribution[a] *=
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
                zxp = Vector3D::sign(size_parameter[5]);
                gama = size_parameter[6];
            }
            else
            {
                au = size_parameter[1];
                zeta = abs(size_parameter[2]);
                zxp = Vector3D::sign(size_parameter[2]);
                gama = size_parameter[3];
            }
            for(uint a = 0; a < nr_of_dust_species; a++)
                grain_size_distribution[a] *= pow(1.0 + zeta * pow(a_eff[a] / au, gama), zxp);
        }
    }
    else if(size_keyword.find("logn") != std::string::npos)
    {
        // Log-normal distribution
        if(size_parameter[0] == 0 || size_parameter[1] == 0)
        {
            cout << ERROR_LINE << "Centroid or sigma of log-normal cannot be 0!" << endl;
            return false;
        }
        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            double argu = -0.5 * pow(log(a_eff[a] / size_parameter[0]) / size_parameter[1], 2);
            grain_size_distribution[a] = exp(argu) / a_eff[a];
        }
    }
    else if(size_keyword.find("zda") != std::string::npos)
    {
        // ZDA distribution from Camps/Trust benchmark
        for(uint a = 0; a < nr_of_dust_species; a++)
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
            grain_size_distribution[a] = pow(10, log_g);
        }
    }

    // Set size distribution times a_eff^2 and mass of grains a certain size
    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        // Set relative abundance of dust grains times their squared radius
        grain_distribution_x_aeff_sq[a] = grain_size_distribution[a] * a_eff_squared[a];

        // For PAHs, their mass can be calculated as standard grains or specifically for
        // PAHs
        if(size_keyword.find("pah") != std::string::npos)
            mass[a] = getPahMass(a);
        else
            mass[a] = 4.0 / 3.0 * PI * a_eff[a] * a_eff[a] * a_eff[a] * material_density;
    }
    return true;
}

bool CDustComponent::add(double ** size_fraction, CDustComponent * comp, uint ** nr_of_scat_theta_tmp, double *** scat_theta_tmp)
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
        nr_of_scat_phi = comp->getNrOfScatPhi();
        nr_of_scat_mat_elements = comp->getNrOfScatMatElements();
        nr_of_calorimetry_temperatures = comp->getNrOfCalorimetryTemperatures();
        f_cor = comp->getFcorr();
        f_highJ = comp->getFHighJ();
        Q_ref = comp->getQref();
        alpha_Q = comp->getAlphaQ();
        R_rayleigh = comp->getRayleighReductionFactor();
        larm_f = comp->getLarmF();

        // Init dust properties to be filled with grain properties
        initDustProperties();

        // Set grain size grid, which needs to be done only once
        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            a_eff[a] = comp->getEffectiveRadius(a);
            a_eff_squared[a] = comp->getEffectiveRadiusSquared(a);
        }

        // Resize Qtrq and parameters for Henyey-Greenstein phase function
        for(uint i = 0; i < nr_of_dust_species * nr_of_wavelength; i++)
        {
            Qtrq[i].resize(nr_of_incident_angles);
            HG_g_factor[i].resize(nr_of_incident_angles);
            HG_g2_factor[i].resize(nr_of_incident_angles);
            HG_g3_factor[i].resize(nr_of_incident_angles);
        }

        // If the scattering matrix is read in, add them together as well
        if(comp->getScatLoaded())
        {
            initNrOfScatThetaArray();
            SetNrOfScatTheta(nr_of_scat_theta_tmp);
            initScatThetaArray();
            SetScatTheta(scat_theta_tmp);
            initScatteringMatrixArray();
        }

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

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        // Check if the dust grain size is really the same
        if(a_eff[a] != comp->getEffectiveRadius(a))
        {
            cout << ERROR_LINE << "Wrong grain size at position " << a + 1 << " ( " << comp->getEffectiveRadius(a)
                 << " )" << endl;
            return false;
        }

        if(comp->sizeIndexUsed(a, a_min, a_max))
        {
            // Mix size distribution with their relative fraction
            grain_size_distribution[a] += size_fraction[a][0];
            grain_distribution_x_aeff_sq[a] += size_fraction[a][0] * comp->getEffectiveRadiusSquared(a);
            mass[a] += size_fraction[a][1] * comp->getMass(a);
        }
    }

    if(comp->getCalorimetryLoaded())
        for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
        {
            // Check if the calorimetry data is really the same
            if(calorimetry_temperatures[t] != comp->getCalorimetricTemperature(t))
            {
                cout << ERROR_LINE << "Wrong calorimetric temperature at position " << t + 1 << " ( "
                     << comp->getCalorimetricTemperature(t) << " )" << endl;
                return false;
            }
        }

    if(comp->getScatLoaded())
    {
        // Mix scattering matrix for mixture
        for(uint a = 0; a < nr_of_dust_species; a++)
            if(comp->sizeIndexUsed(a))
                for(uint w = 0; w < nr_of_wavelength; w++)
                    for(uint inc = 0; inc < nr_of_incident_angles; inc++)
                        for(uint sph = 0; sph < nr_of_scat_phi; sph++)
                            for(uint sth_mix = 0; sth_mix < nr_of_scat_theta[a][w]; sth_mix++)
                            {
                                uint sth_comp = lower_bound(comp->getScatTheta(a,w),
                                                            comp->getScatTheta(a,w) + comp->getNrOfScatTheta(a,w),
                                                            scat_theta[a][w][sth_mix]) - comp->getScatTheta(a,w);

                                if(scat_theta[a][w][sth_mix] == comp->getScatTheta(a,w,sth_comp))
                                    for(uint i_mat = 0; i_mat < nr_of_scat_mat_elements; i_mat++)
                                        for(uint j_mat = 0; j_mat < nr_of_scat_mat_elements; j_mat++)
                                            sca_mat[a][w][inc][sph][sth_mix](i_mat, j_mat) +=
                                                size_fraction[a][1] *
                                                comp->getScatteringMatrixElement(a, w, inc, sph, sth_comp, i_mat, j_mat);
                                else
                                    for(uint i_mat = 0; i_mat < nr_of_scat_mat_elements; i_mat++)
                                        for(uint j_mat = 0; j_mat < nr_of_scat_mat_elements; j_mat++)
                                            sca_mat[a][w][inc][sph][sth_mix](i_mat, j_mat) =
                                                sca_mat[a][w][inc][sph][sth_mix-1](i_mat, j_mat);
                            }
    }

    if(comp->getCalorimetryLoaded())
    {
        // Mix enthalpy for mixture
        for(uint a = 0; a < nr_of_dust_species; a++)
            for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
                enthalpy[a][t] += size_fraction[a][1] * comp->getEnthalpy(a, t);
    }

    // Show progress
    printIDs();
    cout << "- mixing average cross sections\r";

    // Mix optical properties of the dust grains
    for(uint w = 0; w < nr_of_wavelength; w++)
        for(uint a = 0; a < nr_of_dust_species; a++)
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
                addHGg2(a, w, size_fraction[a][1] * comp->getHGg2(a, w));
                addHGg3(a, w, size_fraction[a][1] * comp->getHGg3(a, w));

                CextMean[a][w] = PI * a_eff_squared[a] * (2.0 * getQext1(a, w) + getQext2(a, w)) / 3.0;
                CabsMean[a][w] = PI * a_eff_squared[a] * (2.0 * getQabs1(a, w) + getQabs2(a, w)) / 3.0;
                CscaMean[a][w] = PI * a_eff_squared[a] * (2.0 * getQsca1(a, w) + getQsca2(a, w)) / 3.0;
            }

    if(comp->isAligned())
    {
        // Show progress
        printIDs();
        cout << "- mixing Qtrq and HG g                         \r";

        // Init variables
        double tmpHGgX, tmpHGg2X, tmpHGg3X, tmpQtrqX;
        double tmpHGgY, tmpHGg2Y, tmpHGg3Y, tmpQtrqY;

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
                    // Get incident angle and value of Qtrq and parameters for Henyey-Greenstein phase function
                    comp->getQtrq(i, i_inc, tmpQtrqX, tmpQtrqY);
                    comp->getHG_g_factor(i, i_inc, tmpHGgX, tmpHGgY);
                    comp->getHG_g2_factor(i, i_inc, tmpHGg2X, tmpHGg2Y);
                    comp->getHG_g3_factor(i, i_inc, tmpHGg3X, tmpHGg3Y);

                    // Add the values on top of the mixture Qtrq and Henyey-Greenstein g
                    // factor
                    Qtrq[i].addValue(i_inc, i_inc * d_ang, size_fraction[a][1] * tmpQtrqY);
                    HG_g_factor[i].addValue(i_inc, i_inc * d_ang, size_fraction[a][1] * tmpHGgY);
                    HG_g2_factor[i].addValue(i_inc, i_inc * d_ang, size_fraction[a][1] * tmpHGg2Y);
                    HG_g3_factor[i].addValue(i_inc, i_inc * d_ang, size_fraction[a][1] * tmpHGg3Y);
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
    // if(comp->getPhaseFunctionID() < phID)
    //     phID = comp->getPhaseFunctionID();

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
            HG_g2_factor[i].createSpline();
            HG_g3_factor[i].createSpline();
        }
    }

    // Create StringID for print parameter
    createStringID(comp);

    return true;
}

uint CDustComponent::getInteractingDust(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen, uint cross_section) const
{
    // Get wavelength index from photon package
    uint w = pp->getDustWavelengthID();

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, *pp);
    double a_max = getSizeMax(grid, *pp);

    // If only one grain size is used -> return this size (interpolation instead?)
    if(a_min == a_max)
        return CMathFunctions::biListIndexSearch(a_min, a_eff, nr_of_dust_species);

    switch(cross_section)
    {
        case CROSS_ABS:
            if(abs_prob != 0)
                return findSizeID(pp, abs_prob[w], a_min, a_max, rand_gen);
            break;

        case CROSS_SCA:
            if(sca_prob != 0)
                return findSizeID(pp, sca_prob[w], a_min, a_max, rand_gen);
            break;

        default:
            if(dust_prob != 0)
                return findSizeID(pp, dust_prob[w], a_min, a_max, rand_gen);
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
                    amount[a] = grain_size_distribution[a] * getCabsMean(a, w);
                    break;

                case CROSS_SCA:
                    amount[a] = grain_size_distribution[a] * getCscaMean(a, w);
                    break;

                default:
                    amount[a] = grain_size_distribution[a];
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

    uint a = findSizeID(pp, prob, a_min, a_max, rand_gen);

    return a;
}

void CDustComponent::calcPACrossSections(uint a, uint w, cross_sections & cs, double theta) const
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
    cs.Cpabs = dCabs * sinsq_th;
    cs.Ccirc = sCcirc * sinsq_th;

    // Convert from efficiencies to cross-sections
    cs *= PI * a_eff_squared[a];
}

void CDustComponent::calcNONPACrossSections(uint a, uint w, cross_sections & cs, double theta, double Rent) const
{
    double sinsq_th = sin(theta);
    sinsq_th *= sinsq_th;

    double c_s_ext = getCext1(a, w);
    double c_p_ext = getCext2(a, w);
    double c_av_ext = getCextMean(a, w);

    double cx_ext = c_av_ext + Rent / 3.0 * (c_s_ext - c_p_ext) * (1 - 3.0 * sinsq_th);
    double cy_ext = c_av_ext + Rent / 3.0 * (c_s_ext - c_p_ext);

    cs.Cext = 0.5 * (cx_ext + cy_ext);
    cs.Cpol = 0.5 * Rent * (c_s_ext - c_p_ext) * sinsq_th;

    double c_s_abs = getCabs1(a, w);
    double c_p_abs = getCabs2(a, w);
    double c_av_abs = getCabsMean(a, w);

    double cx_abs = c_av_abs + Rent / 3.0 * (c_s_abs - c_p_abs) * (1 - 3.0 * sinsq_th);
    double cy_abs = c_av_abs + Rent / 3.0 * (c_s_abs - c_p_abs);

    cs.Cabs = 0.5 * (cx_abs + cy_abs);
    cs.Cpabs = 0.5 * Rent * (c_s_abs - c_p_abs) * sinsq_th;

    double c_s_sca = getCsca1(a, w);
    double c_p_sca = getCsca2(a, w);
    double c_av_sca = getCscaMean(a, w);

    double cx_sca = c_av_sca + Rent / 3.0 * (c_s_sca - c_p_sca) * (1 - 3.0 * sinsq_th);
    double cy_sca = c_av_sca + Rent / 3.0 * (c_s_sca - c_p_sca);

    cs.Csca = 0.5 * (cx_sca + cy_sca);
    cs.Ccirc = 0.5 * getCcirc(a, w) * Rent * sinsq_th;
}

void CDustComponent::calcCrossSections(CGridBasic * grid,
                                       const photon_package & pp,
                                       uint i_density,
                                       uint a,
                                       double mag_field_theta,
                                       cross_sections & cs) const
{
    // Get wavelength index
    uint w = pp.getDustWavelengthID();

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
        calcPACrossSections(a, w, cs, mag_field_theta);
        return;
    }

    // Not perfect alignment
    if((alignment & ALIG_NONPA) == ALIG_NONPA)
    {
        if(R_rayleigh == -1)
        {
            // if R_rayleigh == -1, use alignment efficiency from grid
            double R_tmp = grid->getMagField(pp).length();
            if(R_tmp > 1.0)
                R_tmp = 1.0;
            if(R_tmp < 0.0)
                R_tmp = 0.0;
            calcNONPACrossSections(a, w, cs, mag_field_theta, R_tmp);
        }
        else
        {
            calcNONPACrossSections(a, w, cs, mag_field_theta, R_rayleigh);
        }
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

    // Calculate the parameters for radiative torque alignment
    if((alignment & ALIG_RAT) == ALIG_RAT)
    {
        a_alig = grid->getAlignedRadius(pp, i_density);
        if(a_eff[a] > a_alig)
        {
            if((alignment & ALIG_INTERNAL) == ALIG_INTERNAL)
                Rrat = f_highJ + (1 - f_highJ) * getInternalRAT();
            else
                Rrat = R_rayleigh;
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
                Rgold *= (1.5 * cossq_zeta - 0.5) * (1 + f_cor);
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

        Ridg = 1.5 * cossq_beta - 0.5;

        if((alignment & ALIG_INTERNAL) == ALIG_INTERNAL)
        {
            double cossq_zeta = getInternalIDG(Td, Tg);
            Ridg *= (1.5 * cossq_zeta - 0.5) * (1 + f_cor);
        }
    }

    if(a_eff[a] < a_limit)
        Rent = combinedRFactor(Ridg, Rrat, Rgold);
    else
        Rent = Ridg;

    // Calculate sin(theta)^2
    double sinsq_theta = sin(mag_field_theta);
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
    if(grid->getGasNumberDensity(*cell) < min_gas_density)
        return;

    // Set use dust offset later
    dust_offset = true;

    // Init variables
    double temp_offset;
    double * rel_weight = 0;

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, *cell);
    double a_max = getSizeMax(grid, *cell);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, *cell);

    if(grid->getTemperatureFieldInformation() != TEMP_FULL)
    {
        // Get integration over the dust size distribution
        rel_weight = getRelWeight(a_min, a_max, size_param);
    }

    // Get temperature from grid as offset (assuming only ONE dust/gas temperature in
    // input grid)
    if(use_gas_temp)
        temp_offset = grid->getGasTemperature(*cell);
    else
        temp_offset = grid->getDustTemperature(*cell);

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
                                               bool use_energy_density,
                                               CRandomGenerator * rand_gen)
{
    // Init variables
    double t_dust_new = 0;
    uint wIDnew = 0;

    // Get random number for wavelength of re-emitted photon
    double rnd1 = rand_gen->getRND();

    // Get the interacting dust grain size
    uint a = getInteractingDust(grid, pp, rand_gen);

    // Get dust temperature of the emitting grains in current cell
    t_dust_new = updateDustTemperature(grid, *pp, i_density, a, use_energy_density);

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
    pp->setWavelength(wavelength_list[wIDnew + 1], wIDnew + 1);

    return true;
}

double CDustComponent::updateDustTemperature(CGridBasic * grid,
                                             const photon_package & pp,
                                             uint i_density,
                                             uint a,
                                             bool use_energy_density)
{
    // Init variables
    double temp = 0;

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
    temp = findTemperature(a, abs_rate);

    if(sublimate)
        if(temp >= sub_temp)
            temp = 0;

    // Update min and max temperatures for visualization
    max_temp = max(max_temp,temp);
    // min_temp = min(min_temp,temp);

    return temp;
}

void CDustComponent::calcTemperature(CGridBasic * grid,
                                     cell_basic * cell,
                                     uint i_density,
                                     bool use_energy_density)
{
    // Calculate the temperature only for cells with a density not zero
    if(getNumberDensity(grid, *cell, i_density) == 0)
        return;

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, *cell);
    double a_max = getSizeMax(grid, *cell);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, *cell);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    // Init temporary pointer arrays for absorption rate and temperature
    double * abs_rate = new double[nr_of_dust_species];
    double temp;

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a, a_min, a_max))
        {
            // Check if dust grains should have been stochastically heated
            if(a_eff[a] <= getStochasticHeatingMaxSize())
            {
                // Init and resize spline for absorbed energy per wavelength
                spline abs_rate_per_wl;
                abs_rate_per_wl.resize(nr_of_wavelength);

                // Get radiation field and calculate absorbed energy for each wavelength
                for(uint w = 0; w < nr_of_wavelength; w++)
                {
                    double abs_rate_wl_tmp = grid->getRadiationField(*cell, w) * getCabsMean(a, w);
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
                abs_rate[a] = getAbsRate(grid, *cell, a, use_energy_density);
            }

            // Add offset on absorption rate
            if(dust_offset)
            {
                if(grid->getTemperatureFieldInformation() == TEMP_FULL)
                    abs_rate[a] += grid->getQBOffset(*cell, i_density, a);
                else if(grid->getTemperatureFieldInformation() == TEMP_EFF ||
                        grid->getTemperatureFieldInformation() == TEMP_SINGLE)
                    abs_rate[a] += grid->getQBOffset(*cell, i_density);
            }

            // Calculate temperature from absorption rate
            temp = max(double(TEMP_MIN), findTemperature(a, abs_rate[a]));

            // Consider sublimation temperature
            if(sublimate && grid->getTemperatureFieldInformation() == TEMP_FULL)
                if(temp >= sub_temp)
                    temp = TEMP_MIN;

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
                grid->setDustTemperature(cell, i_density, a, temp);

                // Update min and max temperatures for visualization
                max_temp = max(max_temp,temp);
                // if(temp < min_temp)
                //     min_temp = temp;
            }
        }
        else
        {
            // Set absorption rate to zero
            abs_rate[a] = 0;
        }
    }

    // Get average absorption rate via interpolation
    double avg_abs_rate =
        CMathFunctions::integ_dust_size(a_eff, abs_rate, nr_of_dust_species, a_min, a_max);

    // Calculate average temperature from absorption rate
    double avg_temp = findTemperature(grid, cell, avg_abs_rate);

    // Delete pointer array
    delete[] rel_weight;
    delete[] abs_rate;

    if(sublimate)
        if(avg_temp >= sub_temp)
        {
            // Set temperature to zero
            avg_temp = TEMP_MIN;

            // Remove sublimated dust from grid
            // Not if rad field can be used for stochastic heating later
            if(grid->specLengthIsVector())
                grid->adjustDustDensity(cell, i_density, 0);
        }

    // Set average dust temperature in grid
    grid->setDustTemperature(cell, i_density, max(double(TEMP_MIN), avg_temp));

    // Update min and max temperatures for visualization
    max_temp = max(max_temp,avg_temp);
    // if(avg_temp < min_temp)
    //     min_temp = avg_temp;
}

void CDustComponent::calcAlignedRadii(CGridBasic * grid, cell_basic * cell, uint i_density)
{
    // For details, see Reissl et al. 2020, A&A 640:A118, doi:10.1051/0004-6361/201937177

    // Calculate the aligned radii only for cells with a non-zero density
    if(getNumberDensity(grid, *cell, i_density) == 0)
    {
        grid->setAlignedRadius(cell, i_density, a_eff[nr_of_dust_species - 1]);
        return;
    }

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, *cell);
    double a_max = getSizeMax(grid, *cell);

    // default value of the alignment radius
    double a_alig = getSizeMax(grid, *cell);
    double th = 0;
    double dir = 0;

    // Aspect ratio of the grain
    double s = getAspectRatio();

    // alpha_1 ~ delta
    double alpha_1 = 1; // getDeltaRat();

    // Get grid values
    double T_gas = grid->getGasTemperature(*cell);
    double n_g = grid->getGasNumberDensity(*cell);
    double vol = grid->getVolume(*cell);

    // Get average molecular weight
    double mu = grid->getMu();

    // Get thermal velocity
    double v_th = sqrt(2.0 * con_kB * T_gas / (mu * m_H));

    // Loop over all considered grain sizes
    double omega_old = 0;

    for(uint a = 0; a < nr_of_dust_species; a++)
        if(sizeIndexUsed(a, a_min, a_max))
        {
            // Minor and major axis
            double a_minor = a_eff[a] * pow(s, 2. / 3.);
            double a_major = a_eff[a] * pow(s, -1. / 3.);

            // Moment of inertia along a_1
            double I_p = 8. * PI / 15. * getMaterialDensity(a) * a_minor * pow(a_major, 4);

            // Thermal angular momentum
            double J_th = sqrt(I_p * con_kB * T_gas);

            // Init. pointer arrays
            double * Gamma_rad = new double[nr_of_wavelength];
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
                grid->getSpecLength(*cell, w, &arr_en_dens, &en_dir);

                // If the radiation field is zero -> set arrays to zero and move on
                if(arr_en_dens == 0)
                {
                    Gamma_rad[w] = 0;
                    du[w] = 0;
                    ddir[w] = 0;
                    dth[w] = 0;
                    continue;
                }

                // Get angle between magnetic field and radiation field
                double theta = grid->getTheta(*cell, en_dir);

                // Anisotropy parameter
                double gamma = en_dir.length() / arr_en_dens;

                // arr_en_dens = 4 * PI * vol * J -> 4 * PI / c * J
                arr_en_dens /= double(vol * con_c);

                du[w] = wavelength_list[w] * arr_en_dens;

                // Radiative torque efficiency as a power-law
                double Qr = Q_ref;

                if(wavelength_list[w] > 1.8 * a_eff[a])
                    Qr = Q_ref / pow(wavelength_list[w] / (1.8 * a_eff[a]), alpha_Q);

                double cos_theta = abs(cos(theta));

                Qr *= cos_theta;

                // Qr=getQrat(a, w, 0.0);
                Gamma_rad[w] =
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
            double omega_frac = CMathFunctions::integ(wavelength_list, Gamma_rad, 0, nr_of_wavelength - 1);

            double tau_drag = tau_gas / (1. + FIR);
            omega_frac *= tau_drag / J_th;

            // Delete pointer array
            delete[] Gamma_rad;
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

double CDustComponent::calcGoldReductionFactor(const Vector3D & v, const Vector3D & B) const
{
    // Init variables
    double s;
    double g = gold_g_factor;
    Vector3D v_proj = (v * B) / B.sq_length() * B;
    double len_vz = v_proj.sq_length();
    double len_vxyz = v.sq_length();
    double R, cossq_beta;

    if(len_vxyz == 3 * len_vz)
        return 0.0;

    if(len_vxyz == len_vz)
        return 0.0;

    if(gold_g_factor == 0)
        return 0.0;

    s = -0.5 * (len_vxyz - 3 * len_vz) / (len_vxyz - len_vz);

    // if(s <= -0.5)
    //     s = -0.4999999999;

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
                (sqrt(g) * asin(sqrt(-s / (1 + g)))) / (s * atanh(sqrt((-s * g) / (1 + s + g)))) - 1 / s;
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
                (sqrt(g) * asinh(sqrt(s / (1 + g)))) / (s * atan(sqrt((s * g) / (1 + s + g)))) - 1 / s;
        }
    }

    R = 1.5 * cossq_beta - 0.5;

    // if(R <= -0.5)
    //     R = -0.499999999999;

    return R;
}

void CDustComponent::calcStochasticHeatingPropabilities(CGridBasic * grid,
                                                        cell_basic * cell,
                                                        uint i_density,
                                                        dlist & wl_list)
{
    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, *cell);
    double a_max = getSizeMax(grid, *cell);

    if(getNumberDensity(grid, *cell, i_density) == 0)
    {
        // If density is zero, set propability to 1 for lowest temperature
        for(uint a = 0; a < nr_of_dust_species; a++)
            if(a_eff[a] <= getStochasticHeatingMaxSize() && sizeIndexUsed(a, a_min, a_max))
                grid->setDustTempProbability(cell, i_density, a, 0, 1.0);
        return;
    }

    // size of wl_list must equal WL_STEPS
    if(wl_list.size() != WL_STEPS)
    {
        cout << ERROR_LINE << "Size of the wavelength array does not match (stochastic heating)" << endl;
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
            {
                uint wID = getWavelengthID(wl_list[w]);
                // the wavelength axis of Cabs may be greater than WL_STEPS if detectors are defined
                abs_rate_per_wl.setValue(w, wl_list[w], grid->getRadiationField(*cell, w) * getCabsMean(a, wID));
            }

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

double CDustComponent::getCalorimetryA(uint a, uint f, uint i, const spline & abs_rate_per_wl) const
{
    // Calculation of A from Eq. (23) in Camps et al. 2015, A&A 580, A87

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

        if(enthalpy_diff > 0)
        {
            // Get calorimetric temperature to index i
            double T_i = getCalorimetricTemperature(i);

            // Calculate A
            res = tab_em_inv[a].getValue(T_i) * PIx4 / enthalpy_diff;
        }
    }
    else
        cout << ERROR_LINE << "Error at the getCalorimetryA fuction (stochastic heating)";

    // Returning A needs to be at least zero
    return max(double(0), res);
}

long double * CDustComponent::getStochasticProbability(uint a, const spline & abs_rate_per_wl) const
{
    // Calculate the propability of a dust grain to have a certain temperature
    // See Sect. 4.3. in Camps et al. 2015, A&A 580, A87

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
    // X_vec[0] = numeric_limits<long double>::min();
    X_vec[0] = 1.0;
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

        // rescale to prevent infinities if propability increases drastically
        if(X_vec[i] > 1e10)
        {
            for(uint j = 0; j <= i; j++)
                X_vec[j] /= X_vec[i];
        }
    }

    // Init sum for normalization
    long double X_sum = 0;

    // Calculate the sum of the propability
    for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
        X_sum += X_vec[t];

    // Perform normalization or reset propability
    if(isinf(X_sum) || isnan(X_sum) || X_sum <= 0)
    {
        for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
            X_vec[t] = 0;
        cout << ERROR_LINE << "Wrong values in stochastic heating calculation!" << endl;
    }
    else
    {
        for(uint t = 0; t < nr_of_calorimetry_temperatures; t++)
            X_vec[t] /= X_sum;
    }

    return X_vec;
}

void CDustComponent::calcEmissivityHz(CGridBasic * grid,
                                      const photon_package & pp,
                                      uint i_density,
                                      StokesVector * dust_emissivity) const
{
    // Get extinction and absorption cross-sections
    double Cext = getCextMean(grid, pp);
    double Cabs = getCabsMean(grid, pp);

    // Get wavelength index
    uint w = pp.getDustWavelengthID();

    // Calculate frequency
    double frequency = con_c / wavelength_list[w];

    // Get Planck emission at this frequency
    double pl_hz = CMathFunctions::planck_hz(frequency, grid->getDustTemperature(pp, i_density));

    // Get number density of dust grains
    double dens_dust = getNumberDensity(grid, pp, i_density);

    // Fill Stokes vector including optical depth
    dust_emissivity->addI(Cabs * pl_hz * dens_dust);
    dust_emissivity->addT(Cext * dens_dust);
}

double CDustComponent::calcEmissivity(CGridBasic * grid, const photon_package & pp, uint i_density) const
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
    uint w = pp.getDustWavelengthID();

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
                        getCabsMean(a, w) * rel_weight[a] * temp_probability * getPlanck(w, temp_dust);
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
                pl_abs_tmp[a] = getCabsMean(a, w) * rel_weight[a] * getPlanck(w, temp_dust);
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
                                                          const photon_package & pp,
                                                          uint i_density,
                                                          const Vector3D & en_dir,
                                                          double energy) const
{
    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Init  and calculate the cross-sections
    cross_sections cs;

    // Get wavelength index of photon package
    uint w = pp.getDustWavelengthID();

    // Get angle between the magnetic field and the photon direction
    double mag_field_theta = !is_align || alignment == ALIG_RND ? 0 : grid->getThetaMag(pp);

    // Get theta of scattering
    double cos_scattering_theta = en_dir * pp.getDirection();
    double scattering_theta;

    if(cos_scattering_theta < -1.0)
    {
        scattering_theta = PI;
    }
    else if(cos_scattering_theta > 1.0)
    {
        scattering_theta = 0.0;
    }
    else
    {
        scattering_theta = acos(cos_scattering_theta);
    }

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    // Init temporary Stokes array for integration
    StokesVector * scatter_stokes = new StokesVector[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a, a_min, a_max))
        {
            // Get index of theta scattering
            uint thID = phID == PH_MIE ? getScatThetaID(scattering_theta,a,w) : 0;

            // Get cross sections and relative weight of the current dust grain size
            calcCrossSections(grid, pp, i_density, a, mag_field_theta, cs);

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
    double phi_map = CMathFunctions::getRotationAngleObserver(en_dir, pp.getEY(), pp.getEX());

    // Rotate Stokes Vector to be in agreement with the detector plane
    final_stokes.rot(phi_map);

    // Multiply with number density
    final_stokes *= getNumberDensity(grid, pp, i_density);

    return final_stokes;
}

StokesVector CDustComponent::calcEmissivityEmi(CGridBasic * grid,
                                               const photon_package & pp,
                                               uint i_density,
                                               uint emission_component,
                                               double phi,
                                               double energy,
                                               Vector3D en_dir) const
{
    // Init variables to calculate the emission/extinction
    double temp_dust = 0;
    double tmp_planck = 0;
    double scattering_theta = 0, phi_map = 0;
    uint temp_info = grid->getTemperatureFieldInformation();

    // Precalculate values for scattering
    if(energy > 1e-200)
    {
        scattering_theta = acos(en_dir * pp.getDirection());
        phi_map = CMathFunctions::getRotationAngleObserver(en_dir, pp.getEY(), pp.getEX());
    }

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Init  and calculate the cross-sections
    cross_sections cs;

    // Get wavelength index of photon package
    uint w = pp.getDustWavelengthID();

    // Get angle between the magnetic field and the photon direction
    double mag_field_theta = !is_align || alignment == ALIG_RND ? 0 : grid->getThetaMag(pp);

    // Calculate orientation of the Stokes vector in relation to the magnetic field
    double sin_2ph = sin(2.0 * phi);
    double cos_2ph = cos(2.0 * phi);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    // If dust temperature is const for all grains, calc planck only once
    if(temp_info != TEMP_FULL)
    {
        temp_dust = grid->getDustTemperature(pp, i_density);
        tmp_planck = getPlanck(w, temp_dust);
    }

    // Init temporary Stokes array for integration
    StokesVector * tmp_stokes = new StokesVector[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a, a_min, a_max))
        {
            // Get cross sections and relative weight of the current dust grain size
            calcCrossSections(grid, pp, i_density, a, mag_field_theta, cs);

            // Calculate emission/extinction according to the information inside of the
            // grid
            if(getEffectiveRadius(a) <= getStochasticHeatingMaxSize())
            {
                if(emission_component == DUST_EMI_FULL || emission_component == DUST_EMI_STOCH)
                {
                    // Get current calorimetric temperature
                    for(uint t = 0; t < getNrOfCalorimetryTemperatures(); t++)
                    {
                        // Get current calorimetric temperature
                        temp_dust = getCalorimetricTemperature(t);

                        // Get propability that the dust grains have the current calorimetric
                        // temperature
                        double pl = grid->getDustTempProbability(pp, i_density, a, t);

                        // Get relative Planck emission
                        pl *= rel_weight[a] * getPlanck(w, temp_dust);

#if BENCHMARK == CAMPS
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
            }
            else if(emission_component == DUST_EMI_FULL || emission_component == DUST_EMI_TEMP)
            {
                // Consider the temperature of every dust grain size or an average
                // temperature
                if(temp_info == TEMP_FULL)
                {
                    temp_dust = grid->getDustTemperature(pp, i_density, a);
                    tmp_planck = getPlanck(w, temp_dust);
                }

                double pl = rel_weight[a] * tmp_planck;

#if BENCHMARK == CAMPS
                // To perform Camps et. al (2015) benchmark.
                tmp_stokes[a].addI(cs.Cabs * pl);
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
                    // Get index of theta scattering
                    uint thID = getScatThetaID(scattering_theta, a, w);

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

#if BENCHMARK == CAMPS
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
                                          const photon_package & pp,
                                          uint i_density,
                                          double * avg_Cext,
                                          double * avg_Cpol,
                                          double * avg_Ccirc) const
{
    // Init  and calculate the cross-sections
    cross_sections cs;

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get angle between the magnetic field and the photon direction
    double mag_field_theta = !is_align || alignment == ALIG_RND ? 0 : grid->getThetaMag(pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    // Init temporary cross-section array for integration
    double * Cext = new double[nr_of_dust_species];
    double * Cpol = new double[nr_of_dust_species];
    double * Ccirc = new double[nr_of_dust_species];

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        if(sizeIndexUsed(a, a_min, a_max))
        {
            // Get cross sections and relative weight of the current dust grain size
            calcCrossSections(grid, pp, i_density, a, mag_field_theta, cs);

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
    *avg_Cext = CMathFunctions::integ_dust_size(a_eff, Cext, nr_of_dust_species, a_min, a_max);
    *avg_Cpol = CMathFunctions::integ_dust_size(a_eff, Cpol, nr_of_dust_species, a_min, a_max);
    *avg_Ccirc = CMathFunctions::integ_dust_size(a_eff, Ccirc, nr_of_dust_species, a_min, a_max);

    // Delete pointer arrays
    delete[] rel_weight;
    delete[] Cext;
    delete[] Cpol;
    delete[] Ccirc;
}

void CDustComponent::getEscapePhoton(CGridBasic * grid,
                                     photon_package * pp,
                                     uint a,
                                     Vector3D obs_ex,
                                     Vector3D dir_obs,
                                     photon_package * pp_escape) const
{
    switch(phID)
    {
        case PH_MIE:
            getEscapePhotonMie(grid, pp, a, obs_ex, dir_obs, pp_escape);
            break;
        default:
        {
            // Get wavelength index of the photon package
            uint w = pp->getDustWavelengthID();

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
            StokesVector tmp_stokes = *pp->getStokesVector();

            // Reduce the photon package Stokes vector by scattering fraction
            tmp_stokes *= scattered_fraction;

            // Init temporary photon package
            pp_escape->setPosition(pp->getPosition());
            pp_escape->setPositionCell(pp->getPositionCell());

            // Set the photon package at the position of the current photon
            pp_escape->setDirection(dir_obs);
            pp_escape->setWavelength(wavelength_list[w], w);

            // Synchronize the direction and wavelength as well
            pp_escape->setStokesVector(tmp_stokes);
            break;
        }
    }
}

void CDustComponent::getEscapePhotonMie(CGridBasic * grid,
                                        photon_package * pp,
                                        uint a,
                                        Vector3D obs_ex,
                                        Vector3D dir_obs,
                                        photon_package * pp_escape) const
{
    // Get wavelength index of the photon package
    uint w = pp->getDustWavelengthID();

    // Get the Stokes vector of the current photon package
    StokesVector tmp_stokes = *pp->getStokesVector();

    // Init temporary photon package
    pp_escape->setPosition(pp->getPosition());
    pp_escape->setPositionCell(pp->getPositionCell());

    // Set the photon package at the position of the current photon
    pp_escape->setD(pp->getD());
    pp_escape->setWavelength(wavelength_list[w], w);

    // Determination of the scattering angle (phi, theta) towards the observing map in the
    // photon frame. Get the rotation matrix of the photon (photon space to lab space)
    Matrix2D D_photon = pp_escape->getD();
    D_photon.transpose();
    Vector3D dir_rlp = D_photon * dir_obs;

    // Calculate the theta and phi angle to the observer
    double phi_photon_to_obs = Vector3D::atan3(dir_rlp.Y(), -dir_rlp.X());
    double theta_photon_to_obs = acos(dir_rlp.Z());

    // Update the coordinate space of the photon
    pp_escape->updateCoordSystem(phi_photon_to_obs, theta_photon_to_obs);

    // Get the theta angle index to obtain the scattering matrix
    uint thID = getScatThetaID(theta_photon_to_obs,a,w);

    // Create the scattering matrix with the local parameters
    const Matrix2D & mat_sca = getScatteringMatrix(a, w, 0, 0, thID);

    // Default phi distribution is isotropic
    double phi_fraction = 1;
    // Get PHIPAR to take non equal distribution of phi angles into account
    if(tmp_stokes.I() > 0.0 && mat_sca(0, 0) > 0.0)
    {
        double q_i = tmp_stokes.Q() / tmp_stokes.I();
        double u_i = tmp_stokes.U() / tmp_stokes.I();
        double phipar = sqrt(q_i * q_i + u_i * u_i) * (-mat_sca(0, 1) / mat_sca(0, 0));

        double gamma = 0.5 * Vector3D::atan3(tmp_stokes.Q(), tmp_stokes.U());
        double cos_2_phi = cos(2.0 * (phi_photon_to_obs + gamma - PI));

        // Calculate the fraction that is scattered into this phi direction
        phi_fraction = (1.0 - phipar * cos_2_phi);
    }
    else
    {
        if(tmp_stokes.I() <= 0.0)
            cout << ERROR_LINE << "Photon package intensity is zero or negative!\n" << endl;
        if(mat_sca(0, 0) <= 0.0)
            cout << ERROR_LINE << "First scattering matrix element is zero or negative!\n" << endl;
        
        tmp_stokes.clear();
        pp_escape->setStokesVector(tmp_stokes);
        return;
    }

    // Calculate the fraction that is scattered into this theta direction
    double theta_fraction = getScatteredFractionMie(a, w, theta_photon_to_obs);

    // Reduce Stokes vector by scattering propability into theta and phi
    tmp_stokes *= theta_fraction * phi_fraction;

    // Backup Stokes vector
    double stokes_1_bak = tmp_stokes.I();
    if(stokes_1_bak > 1e-200)
    {
        // Rotate Stokes vector to new photon direction
        tmp_stokes.rot(phi_photon_to_obs);
        // Multiply Stokes vector with scattering matrix
        tmp_stokes *= mat_sca;
        // Normalize Stokes vector to preserve total intensity
        tmp_stokes *= stokes_1_bak / tmp_stokes.I();
    }
    else
    {
        tmp_stokes.clear();
        pp_escape->setStokesVector(tmp_stokes);
        return;
    }

    // Rotate photon package into the coordinate space of the detector
    // (x or r)-axis of photon is y-axis of detector
    // (y or l)-axis of photon is negative x-axis of detector
    // see Figure 12 in O. Fischer (1993)
    double rot_angle_phot_obs =
        CMathFunctions::getRotationAngleObserver(obs_ex, pp_escape->getEX(), -1*pp_escape->getEY());
    tmp_stokes.rot(rot_angle_phot_obs);

    // The scattering part is based on O. Fischer (1993)
    // But on our detectors, U is defined the other way round
    tmp_stokes.multU(-1);

    // Set the new Stokes vector to the photon package
    pp_escape->setStokesVector(tmp_stokes);
}

double CDustComponent::getCellEmission(CGridBasic * grid, const photon_package & pp, uint i_density) const
{
    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    // Get Volume of current cell
    double vol = grid->getVolume(pp);

    // Get dust number density of current cell
    double dens = getNumberDensity(grid, pp, i_density);

    // Get wavelength of photon package
    uint w = pp.getDustWavelengthID();

    if(grid->getTemperatureFieldInformation() == TEMP_FULL)
    {
        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            // Get dust temperature from grid
            double temp = grid->getDustTemperature(pp, i_density, a);

            // Calculate energy of current grain size
            rel_weight[a] *= getCabsMean(a, w) * getPlanck(w, temp);
        }
    }
    else
    {
        // Get dust temperature from grid
        double temp = grid->getDustTemperature(pp, i_density);
        double planck_tmp = getPlanck(w, temp);

        // Calculate energy of current grain size
        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] *= getCabsMean(a, w) * planck_tmp;
    }

    // Calculate the total energy via integration
    double total_energy =
        dens * vol * PIx4 *
        CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);

    // Delete pointer array
    delete[] rel_weight;

    return total_energy;
}

void CDustComponent::henyeygreen(photon_package * pp, uint a, CRandomGenerator * rand_gen)
{
    // Calculate a new random direction based on the Henyey-Greenstein function
    // Henyey & Greenstein 1941, ApJ 93, 70

    // Init variables
    double cos_theta, theta, phi;

    // Get the current wavelength
    double w = pp->getDustWavelengthID();

    // Get the Henyey-Greenstein g factor
    double g = getHGg(a, w);

    if(g <= -1.0 || g >= 1.0) {
        cout << ERROR_LINE << "Henyey-Greenstein g factor is invalid: " << g << endl;
        return;
    }

    double rnd = rand_gen->getRND();
    if(abs(g) < 2.0 * EPS_DOUBLE) { // If g is close to zero, use random isotropic direction
        cos_theta = 2.0 * rnd - 1.0;
    } else { // get cosine theta, see e.g. Eq. (19) in Witt 1977, ApJS 35, 1
        cos_theta = (1.0 - g * g) / (1.0 - g + 2.0 * g * rnd);
        cos_theta = 1.0 + g * g - cos_theta * cos_theta;
        cos_theta /= (2.0 * g);
    }

    if(cos_theta < -1.0) {
        cos_theta = -1.0;
    }
    if(cos_theta > 1.0) {
        cos_theta = 1.0;
    }

    // Get theta from cosine theta
    theta = acos(cos_theta);

    // equal distribution of the phi angles
    phi = PIx2 * rand_gen->getRND();

    // Update the photon package with the new direction
    pp->updateCoordSystem(phi, theta);
}

void CDustComponent::drainehenyeygreen(photon_package * pp, uint a, CRandomGenerator * rand_gen)
{
    // Calculate a new random direction based on the Draine Henyey-Greenstein function
    // Draine 2003, ApJ 598, 1017

    // Init variables
    double cos_theta, theta, phi;

    // Get the current wavelength
    double w = pp->getDustWavelengthID();

    // Get the parameters for Henyey-Greenstein phase function
    double g = getHGg(a, w);
    double alpha = getHGg2(a, w);

    if(g <= -1.0 || g >= 1.0) {
        cout << ERROR_LINE << "Henyey-Greenstein g factor is invalid: " << g << endl;
        return;
    }

    if(alpha < 0.0 || alpha > 1.0) {
        cout << ERROR_LINE << "Henyey-Greenstein alpha factor is invalid: " << alpha << endl;
        return;
    }

    if(alpha < 2.0 * EPS_DOUBLE) { // If alpha is close to zero, use Henyey-Greenstein
        henyeygreen(pp, a, rand_gen);
        return;
    }

    double rnd = rand_gen->getRND();
    if(abs(g) < 2.0 * EPS_DOUBLE) {
        // If g is close to zero, use Cardano's formula
        // the integral of the phase function leads to a cubic equation
        // x^3 + px + q = 0
        // where x = cos(theta), p = 3/alpha, q = (1 + 3/alpha - 2*rnd*(1 + 3/alpha))
        // solve for x with Cardan's formula
        // we only get one real solution, since discr = (q/2)^2 + (p/3)^3 > 0
        double coeff_p = 3.0 / alpha;
        double coeff_q = (1.0 + coeff_p - 2.0 * rnd * (1.0 + coeff_p));
        double discr = 0.25 * coeff_q * coeff_q + coeff_p * coeff_p * coeff_p / 27.0;

        double u_plus = cbrt(-0.5 * coeff_q + sqrt(discr));
        // double u_minus = cbrt(-0.5 * coeff_q - sqrt(discr));
        // cos_theta = u_plus + u_minus;

        // since u_plus * u_minus = -p/3, we can also write cos_theta as
        cos_theta = u_plus - coeff_p / (3.0 * u_plus);
        // or
        // cos_theta = u_minus - coeff_p / (3.0 * u_minus);
    } else {
        // find cos_theta with Brent's method
        cos_theta = CMathFunctions::findRootBrent(-1.0, 1.0, &CMathFunctions::getDHGIntegral, {g, alpha, rnd});
    }

    if(cos_theta < -1.0) {
        cos_theta = -1.0;
    }
    if(cos_theta > 1.0) {
        cos_theta = 1.0;
    }

    // Get theta from cosine theta
    theta = acos(cos_theta);

    // equal distribution of the phi angles
    phi = PIx2 * rand_gen->getRND();

    // Update the photon package with the new direction
    pp->updateCoordSystem(phi, theta);
}

void CDustComponent::threeparamhenyeygreen(photon_package * pp, uint a, CRandomGenerator * rand_gen)
{
    // Calculate a new random direction based on the three parameter Henyey-Greenstein function
    // Kattawar 1975, JQSRT 15, 839
    // Witt 1977, ApJS 31, 1

    // Init variables
    double cos_theta, theta, phi;

    // Get the current wavelength
    double w = pp->getDustWavelengthID();

    // Get the parameters for Henyey-Greenstein phase function
    double g1 = getHGg(a, w);
    double g2 = getHGg2(a, w);
    double weight = getHGg3(a, w);

    if(g1 <= -1.0 || g1 >= 1.0) {
        cout << ERROR_LINE << "Henyey-Greenstein g factor is invalid: " << g1 << endl;
        return;
    }

    if(g2 <= -1.0 || g2 >= 1.0) {
        cout << ERROR_LINE << "Henyey-Greenstein g2 factor is invalid: " << g2 << endl;
        return;
    }

    if(g1 * g2 > 0.0) {
        cout << ERROR_LINE << "Henyey-Greenstein g1 and g2 must have different signs." << endl;
        return;
    }
    
    if(weight < 0.0 || weight > 1.0) {
        cout << ERROR_LINE << "Henyey-Greenstein weight factor is invalid: " << weight << endl;
        return;
    }

    if(1.0 - weight < 2.0 * EPS_DOUBLE) { // If weight is close to one, use Henyey-Greenstein
        henyeygreen(pp, a, rand_gen);
        return;
    }

    double rnd = rand_gen->getRND();
    if(abs(g1) < 2.0 * EPS_DOUBLE && abs(g2) < 2.0 * EPS_DOUBLE){ // if both g1 and g2 are close to zero, use random isotropic direction
        cos_theta = 2.0 * rnd - 1.0;
    } else { // find cos_theta with Brent's method
        cos_theta = CMathFunctions::findRootBrent(-1.0, 1.0, &CMathFunctions::getTTHGIntegral, {g1, g2, weight, rnd});
    }

    if(cos_theta < -1.0) {
        cos_theta = -1.0;
    }
    if(cos_theta > 1.0) {
        cos_theta = 1.0;
    }

    // Get theta from cosine theta
    theta = acos(cos_theta);

    // equal distribution of the phi angles
    phi = PIx2 * rand_gen->getRND();

    // Update the photon package with the new direction
    pp->updateCoordSystem(phi, theta);
}

void CDustComponent::miesca(photon_package * pp, uint a, CRandomGenerator * rand_gen)
{
    // Get wavelength of photon package
    uint w = pp->getDustWavelengthID();

    // Get theta angle from distribution
    double theta = findTheta(a, w, rand_gen->getRND());
    double phi;

    // Get theta index from theta angle
    uint thID = getScatThetaID(theta,a,w);

    // Get Stokes vector from photon package
    StokesVector tmp_stokes = *pp->getStokesVector();

    // Get scattering matrix
    const Matrix2D & mat_sca = getScatteringMatrix(a, w, 0, 0, thID);

    double phipar;
    // Get PHIPAR to take non equal distribution of phi angles into account
    if(tmp_stokes.I() > 0.0 && mat_sca(0, 0) > 0.0) {
        double q_i = tmp_stokes.Q() / tmp_stokes.I();
        double u_i = tmp_stokes.U() / tmp_stokes.I();
        phipar = sqrt(q_i * q_i + u_i * u_i) * (-mat_sca(0, 1) / mat_sca(0, 0));
    } else {
        if(tmp_stokes.I() <= 0.0) {
            cout << ERROR_LINE << "Photon package intensity is zero or negative!\n" << endl;
        }
        if(mat_sca(0, 0) <= 0.0) {
            cout << ERROR_LINE << "First scattering matrix element is zero or negative!\n" << endl;
        }

        tmp_stokes.clear();
        pp->setStokesVector(tmp_stokes);
        return;
    }

    double gamma = 0.5 * Vector3D::atan3(tmp_stokes.Q(), tmp_stokes.U());

    // find phi (Kepler's equation) with Brent's method
    phi = CMathFunctions::findRootBrent(0.0, PIx2, &CMathFunctions::getPhiIntegral, {phipar, rand_gen->getRND()});
    phi = PI - gamma + phi;

    // Update the photon package with the new direction
    pp->updateCoordSystem(phi, theta);

    double i_1 = tmp_stokes.I();
    if(i_1 > 1e-200) {
        tmp_stokes.rot(phi);
        tmp_stokes *= mat_sca;
        tmp_stokes *= i_1 / tmp_stokes.I();
    } else {
        tmp_stokes.clear();
    }

    pp->setStokesVector(tmp_stokes);
}

// ----------------------------------------------------------------------
// ----------- Efficiencies for grain size and wavelength ---------------
// ----------------------------------------------------------------------
inline double CDustComponent::getQext1(uint a, uint w) const
{
    return Qext1[a][w];
}

inline double CDustComponent::getQext2(uint a, uint w) const
{
    return Qext2[a][w];
}

inline double CDustComponent::getQabs1(uint a, uint w) const
{
    return Qabs1[a][w];
}

inline double CDustComponent::getQabs2(uint a, uint w) const
{
    return Qabs2[a][w];
}

inline double CDustComponent::getQsca1(uint a, uint w) const
{
    return Qsca1[a][w];
}

inline double CDustComponent::getQsca2(uint a, uint w) const
{
    return Qsca2[a][w];
}

inline double CDustComponent::getQcirc(uint a, uint w) const
{
    return Qcirc[a][w];
}

inline double CDustComponent::getHGg(uint a, uint w) const
{
    return HGg[a][w];
}

inline double CDustComponent::getHGg2(uint a, uint w) const
{
    return HGg2[a][w];
}

inline double CDustComponent::getHGg3(uint a, uint w) const
{
    return HGg3[a][w];
}

void CDustComponent::setHGg(uint a, uint w, double val)
{
    HGg[a][w] = val;
}

void CDustComponent::setHGg2(uint a, uint w, double val)
{
    HGg2[a][w] = val;
}

void CDustComponent::setHGg3(uint a, uint w, double val)
{
    HGg3[a][w] = val;
}

// ------------------------------------------------------------------------------------
// ----------- Add values to efficiencies for grain size and wavelength ---------------
// ------------------------------------------------------------------------------------
void CDustComponent::addQext1(uint a, uint w, double val)
{
    Qext1[a][w] += val;
}

void CDustComponent::addQext2(uint a, uint w, double val)
{
    Qext2[a][w] += val;
}

void CDustComponent::addQabs1(uint a, uint w, double val)
{
    Qabs1[a][w] += val;
}

void CDustComponent::addQabs2(uint a, uint w, double val)
{
    Qabs2[a][w] += val;
}

void CDustComponent::addQsca1(uint a, uint w, double val)
{
    Qsca1[a][w] += val;
}

void CDustComponent::addQsca2(uint a, uint w, double val)
{
    Qsca2[a][w] += val;
}

void CDustComponent::addQcirc(uint a, uint w, double val)
{
    Qcirc[a][w] += val;
}

void CDustComponent::addHGg(uint a, uint w, double val)
{
    HGg[a][w] += val;
}

void CDustComponent::addHGg2(uint a, uint w, double val)
{
    HGg2[a][w] += val;
}

void CDustComponent::addHGg3(uint a, uint w, double val)
{
    HGg3[a][w] += val;
}

// ------------------------------------------------------------------------
// ----------- Cross-sections for grain size and wavelength ---------------
// ------------------------------------------------------------------------
double CDustComponent::getCext1(uint a, uint w) const
{
    return PI * a_eff_squared[a] * getQext1(a, w);
}

double CDustComponent::getCext2(uint a, uint w) const
{
    return PI * a_eff_squared[a] * getQext2(a, w);
}

double CDustComponent::getCabs1(uint a, uint w) const
{
    return PI * a_eff_squared[a] * getQabs1(a, w);
}

double CDustComponent::getCabs2(uint a, uint w) const
{
    return PI * a_eff_squared[a] * getQabs2(a, w);
}

double CDustComponent::getCsca1(uint a, uint w) const
{
    return PI * a_eff_squared[a] * getQsca1(a, w);
}

double CDustComponent::getCsca2(uint a, uint w) const
{
    return PI * a_eff_squared[a] * getQsca2(a, w);
}

double CDustComponent::getCcirc(uint a, uint w) const
{
    return PI * a_eff_squared[a] * getQcirc(a, w);
}

// -------------------------------------------------------------------------------
// ----------- Cross-sections for wavelength mixed in current cell ---------------
// -------------------------------------------------------------------------------
double CDustComponent::getCext1(CGridBasic * grid, const photon_package & pp) const
{
    // Get wavelength of photon package
    uint w = pp.getDustWavelengthID();

    // Return precalculated value if available
    if(tCext1 != 0)
        return tCext1[w];

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_weight[a] *= a_eff_squared[a] * getQext1(a, w);
    double res =
        PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
    delete[] rel_weight;
    return res;
}

double CDustComponent::getCext2(CGridBasic * grid, const photon_package & pp) const
{
    // Get wavelength of photon package
    uint w = pp.getDustWavelengthID();

    // Return precalculated value if available
    if(tCext2 != 0)
        return tCext2[w];

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_weight[a] *= a_eff_squared[a] * getQext2(a, w);
    double res =
        PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
    delete[] rel_weight;
    return res;
}

double CDustComponent::getCabs1(CGridBasic * grid, const photon_package & pp) const
{
    // Get wavelength of photon package
    uint w = pp.getDustWavelengthID();

    // Return precalculated value if available
    if(tCabs1 != 0)
        return tCabs1[w];

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_weight[a] *= a_eff_squared[a] * getQabs1(a, w);
    double res =
        PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
    delete[] rel_weight;
    return res;
}

double CDustComponent::getCabs2(CGridBasic * grid, const photon_package & pp) const
{
    // Get wavelength of photon package
    uint w = pp.getDustWavelengthID();

    // Return precalculated value if available
    if(tCabs2 != 0)
        return tCabs2[w];

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_weight[a] *= a_eff_squared[a] * getQabs2(a, w);
    double res =
        PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
    delete[] rel_weight;
    return res;
}

double CDustComponent::getCsca1(CGridBasic * grid, const photon_package & pp) const
{
    // Get wavelength of photon package
    uint w = pp.getDustWavelengthID();

    // Return precalculated value if available
    if(tCsca1 != 0)
        return tCsca1[w];

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_weight[a] *= a_eff_squared[a] * getQsca1(a, w);
    double res =
        PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
    delete[] rel_weight;
    return res;
}

double CDustComponent::getCsca2(CGridBasic * grid, const photon_package & pp) const
{
    // Get wavelength of photon package
    uint w = pp.getDustWavelengthID();

    // Return precalculated value if available
    if(tCsca2 != 0)
        return tCsca2[w];

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_weight[a] *= a_eff_squared[a] * getQsca2(a, w);
    double res =
        PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
    delete[] rel_weight;
    return res;
}

double CDustComponent::getCcirc(CGridBasic * grid, const photon_package & pp) const
{
    // Get wavelength of photon package
    uint w = pp.getDustWavelengthID();

    // Return precalculated value if available
    if(tCcirc != 0)
        return tCcirc[w];

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_weight[a] *= a_eff_squared[a] * getQcirc(a, w);
    double res =
        PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
    delete[] rel_weight;
    return res;
}

double CDustComponent::getHGg(CGridBasic * grid, const photon_package & pp) const
{
    // Get wavelength of photon package
    uint w = pp.getDustWavelengthID();

    // Return precalculated value if available
    if(tHGg != 0)
        return tHGg[w];

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_weight[a] *= getHGg(a, w);
    double res = CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
    delete[] rel_weight;
    return res;
}

double CDustComponent::getHGg2(CGridBasic * grid, const photon_package & pp) const
{
    // Get wavelength of photon package
    uint w = pp.getDustWavelengthID();

    // Return precalculated value if available
    if(tHGg2 != 0)
        return tHGg2[w];

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_weight[a] *= getHGg2(a, w);
    double res = CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
    delete[] rel_weight;
    return res;
}

double CDustComponent::getHGg3(CGridBasic * grid, const photon_package & pp) const
{
    // Get wavelength of photon package
    uint w = pp.getDustWavelengthID();

    // Return precalculated value if available
    if(tHGg3 != 0)
        return tHGg3[w];

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, pp);
    double a_max = getSizeMax(grid, pp);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, pp);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_weight[a] *= getHGg3(a, w);
    double res = CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
    delete[] rel_weight;
    return res;
}

// -----------------------------------------------------------------------------
// ----------- Average cross-sections for grain size and wavelength ------------
// -----------------------------------------------------------------------------

inline double CDustComponent::getCextMean(uint a, uint w) const
{
    return CextMean[a][w];
}

inline double CDustComponent::getCabsMean(uint a, uint w) const
{
    return CabsMean[a][w];
}

inline double CDustComponent::getCscaMean(uint a, uint w) const
{
    return CscaMean[a][w];
}

// -----------------------------------------------------------------------------
// ----------- Average cross-sections for wavelength mixed in current cell -----
// -----------------------------------------------------------------------------
double CDustComponent::getCextMean(CGridBasic * grid, const photon_package & pp) const
{
    return (2.0 * getCext1(grid, pp) + getCext2(grid, pp)) / 3.0;
}

double CDustComponent::getCabsMean(CGridBasic * grid, const photon_package & pp) const
{
    return (2.0 * getCabs1(grid, pp) + getCabs2(grid, pp)) / 3.0;
}

double CDustComponent::getCscaMean(CGridBasic * grid, const photon_package & pp) const
{
    return (2.0 * getCsca1(grid, pp) + getCsca2(grid, pp)) / 3.0;
}

// ------------------------------------------------------------------------------------
// ----------- Average cross-sections for wavelength mixed with global conditions -----
// ------------------------------------------------------------------------------------
double CDustComponent::getCextMean(double w) const
{
    return (2.0 * getCext1(w) + getCext2(w)) / 3.0;
}

double CDustComponent::getCabsMean(double w) const
{
    return (2.0 * getCabs1(w) + getCabs2(w)) / 3.0;
}

double CDustComponent::getCscaMean(double w) const
{
    return (2.0 * getCsca1(w) + getCsca2(w)) / 3.0;
}

// ------------------------------------------------------------------------------
// ----------- Cross-sections for wavelength mixed with global limits -----------
// ------------------------------------------------------------------------------
double CDustComponent::getCext1(uint w) const
{
    double * Cext1 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Cext1[a] = grain_distribution_x_aeff_sq[a] * getQext1(a, w);
    double res =
        PI / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Cext1, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Cext1;
    return res;
}

double CDustComponent::getCext2(uint w) const
{
    double * Cext2 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Cext2[a] = grain_distribution_x_aeff_sq[a] * getQext2(a, w);
    double res =
        PI / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Cext2, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Cext2;
    return res;
}

double CDustComponent::getCabs1(uint w) const
{
    double * Cabs1 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Cabs1[a] = grain_distribution_x_aeff_sq[a] * getQabs1(a, w);
    double res =
        PI / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Cabs1, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Cabs1;
    return res;
}

double CDustComponent::getCabs2(uint w) const
{
    double * Cabs2 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Cabs2[a] = grain_distribution_x_aeff_sq[a] * getQabs2(a, w);
    double res =
        PI / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Cabs2, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Cabs2;
    return res;
}

double CDustComponent::getCsca1(uint w) const
{
    double * Csca1 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Csca1[a] = grain_distribution_x_aeff_sq[a] * getQsca1(a, w);
    double res =
        PI / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Csca1, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Csca1;
    return res;
}

double CDustComponent::getCsca2(uint w) const
{
    double * Csca2 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Csca2[a] = grain_distribution_x_aeff_sq[a] * getQsca2(a, w);
    double res =
        PI / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Csca2, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Csca2;
    return res;
}

double CDustComponent::getCcirc(uint w) const
{
    double * Ccirc = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Ccirc[a] = grain_distribution_x_aeff_sq[a] * getQcirc(a, w);
    double res =
        PI / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Ccirc, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Ccirc;
    return res;
}

double CDustComponent::getHGg(uint w) const
{
    double * HGg = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        HGg[a] = grain_size_distribution[a] * getHGg(a, w);
    double res =
        1.0 / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, HGg, nr_of_dust_species, a_min_global, a_max_global);
    delete[] HGg;
    return res;
}

double CDustComponent::getHGg2(uint w) const
{
    double * HGg2 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        HGg2[a] = grain_size_distribution[a] * getHGg2(a, w);
    double res =
        1.0 / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, HGg2, nr_of_dust_species, a_min_global, a_max_global);
    delete[] HGg2;
    return res;
}

double CDustComponent::getHGg3(uint w) const
{
    double * HGg3 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        HGg3[a] = grain_size_distribution[a] * getHGg3(a, w);
    double res =
        1.0 / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, HGg3, nr_of_dust_species, a_min_global, a_max_global);
    delete[] HGg3;
    return res;
}

// ------------------------------------------------------------------------------
// ----------- Efficiencies for wavelength mixed with global limits -----------
// ------------------------------------------------------------------------------
double CDustComponent::getQext1(uint w) const
{
    double * Qext1 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Qext1[a] = grain_size_distribution[a] * getQext1(a, w);
    double res =
        1.0 / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Qext1, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Qext1;
    return res;
}

double CDustComponent::getQext2(uint w) const
{
    double * Qext2 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Qext2[a] = grain_size_distribution[a] * getQext2(a, w);
    double res =
        1.0 / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Qext2, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Qext2;
    return res;
}

double CDustComponent::getQabs1(uint w) const
{
    double * Qabs1 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Qabs1[a] = grain_size_distribution[a] * getQabs1(a, w);
    double res =
        1.0 / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Qabs1, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Qabs1;
    return res;
}

double CDustComponent::getQabs2(uint w) const
{
    double * Qabs2 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Qabs2[a] = grain_size_distribution[a] * getQabs2(a, w);
    double res =
        1.0 / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Qabs2, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Qabs2;
    return res;
}

double CDustComponent::getQsca1(uint w) const
{
    double * Qsca1 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Qsca1[a] = grain_size_distribution[a] * getQsca1(a, w);
    double res =
        1.0 / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Qsca1, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Qsca1;
    return res;
}

double CDustComponent::getQsca2(uint w) const
{
    double * Qsca2 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Qsca2[a] = grain_size_distribution[a] * getQsca2(a, w);
    double res =
        1.0 / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Qsca2, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Qsca2;
    return res;
}

double CDustComponent::getQcirc(uint w) const
{
    double * Qcirc = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Qcirc[a] = grain_size_distribution[a] * getQcirc(a, w);
    double res =
        1.0 / getWeight() *
        CMathFunctions::integ_dust_size(a_eff, Qcirc, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Qcirc;
    return res;
}

// -----------------------------------------------------------------------------------
// ----------- Mass cross-sections for wavelength mixed with global limits -----------
// -----------------------------------------------------------------------------------
double CDustComponent::getKappaExt1(uint w) const
{
    double * Cext1 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Cext1[a] = grain_distribution_x_aeff_sq[a] * getQext1(a, w);
    double res =
        PI / getWeight() / getAvgMass() *
        CMathFunctions::integ_dust_size(a_eff, Cext1, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Cext1;
    return res;
}

double CDustComponent::getKappaExt2(uint w) const
{
    double * Cext2 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Cext2[a] = grain_distribution_x_aeff_sq[a] * getQext2(a, w);
    double res =
        PI / getWeight() / getAvgMass() *
        CMathFunctions::integ_dust_size(a_eff, Cext2, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Cext2;
    return res;
}

double CDustComponent::getKappaAbs1(uint w) const
{
    double * Cabs1 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Cabs1[a] = grain_distribution_x_aeff_sq[a] * getQabs1(a, w);
    double res =
        PI / getWeight() / getAvgMass() *
        CMathFunctions::integ_dust_size(a_eff, Cabs1, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Cabs1;
    return res;
}

double CDustComponent::getKappaAbs2(uint w) const
{
    double * Cabs2 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Cabs2[a] = grain_distribution_x_aeff_sq[a] * getQabs2(a, w);
    double res =
        PI / getWeight() / getAvgMass() *
        CMathFunctions::integ_dust_size(a_eff, Cabs2, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Cabs2;
    return res;
}

double CDustComponent::getKappaSca1(uint w) const
{
    double * Csca1 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Csca1[a] = grain_distribution_x_aeff_sq[a] * getQsca1(a, w);
    double res =
        PI / getWeight() / getAvgMass() *
        CMathFunctions::integ_dust_size(a_eff, Csca1, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Csca1;
    return res;
}

double CDustComponent::getKappaSca2(uint w) const
{
    double * Csca2 = new double[nr_of_dust_species];
    for(uint a = 0; a < nr_of_dust_species; a++)
        Csca2[a] = grain_distribution_x_aeff_sq[a] * getQsca2(a, w);
    double res =
        PI / getWeight() / getAvgMass() *
        CMathFunctions::integ_dust_size(a_eff, Csca2, nr_of_dust_species, a_min_global, a_max_global);
    delete[] Csca2;
    return res;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double CDustComponent::getDeltaRat() const
{
    return delta_rat;
}

double * CDustComponent::getRelWeight(double a_min, double a_max, double size_param) const
{
    double * rel_weight = new double[nr_of_dust_species];

    if(getRelWeightTab(0) != MAX_DOUBLE)
        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] = getRelWeightTab(a);
    else
    {
        // Create normalization factor
        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] = grain_size_distribution[a] * pow(a_eff[a], size_param);

        double weight = CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);

        // Create the final relative mass distribution
        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            rel_weight[a] /= weight;
        }
    }

    return rel_weight;
}

double CDustComponent::getRelWeightTab(uint a) const
{
    if(relWeightTab != 0)
        return relWeightTab[a];
    return MAX_DOUBLE;
}

void CDustComponent::preCalcRelWeight()
{
    relWeightTab = new double[nr_of_dust_species];
    fill(relWeightTab,relWeightTab+nr_of_dust_species,MAX_DOUBLE);
    relWeightTab = getRelWeight(a_min_global, a_max_global);
}

double CDustComponent::getWeight() const
{
    double weight =
        CMathFunctions::integ_dust_size(a_eff, grain_size_distribution, nr_of_dust_species, a_min_global, a_max_global);
    return weight;
}

double CDustComponent::getMassWeight()
{
    // Init pointer array of relative mass of each dust grain size bin
    double * rel_mass = new double[nr_of_dust_species];

    // Set the relative mass of each dust grain size bin
    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_mass[a] = grain_size_distribution[a] * mass[a];

    // Calculate the average mass of the dust grains in the current cell
    double mass_weight =
        CMathFunctions::integ_dust_size(a_eff, rel_mass, nr_of_dust_species, a_min_global, a_max_global);

    // Delete pointer array
    delete[] rel_mass;

    return mass_weight;
}

double CDustComponent::getGoldFactor()
{
    return gold_g_factor;
}

double CDustComponent::getLarmF()
{
    return larm_f;
}

double CDustComponent::getQref()
{
    return Q_ref;
}

double CDustComponent::getAlphaQ()
{
    return alpha_Q;
}

double CDustComponent::getRayleighReductionFactor()
{
    return R_rayleigh;
}

void CDustComponent::setAlignmentMechanism(uint al)
{
    alignment = al;
}

void CDustComponent::setPhaseFunctionID(uint ph)
{
    phID = ph;
}

uint CDustComponent::getPhaseFunctionID()
{
    return phID;
}

uint CDustComponent::getNrOfStochasticSizes()
{
    uint nr_stochastic_sizes = 0;
    if(stochastic_heating_max_size > 0)
        for(uint a = 0; a < nr_of_dust_species; a++)
            if(a_eff[a] <= stochastic_heating_max_size)
                nr_stochastic_sizes++;
    return nr_stochastic_sizes;
}

// double getMinDustTemp()
// {
//     return min_temp;
// }

double CDustComponent::getMaxDustTemp()
{
    return max_temp;
}

double CDustComponent::getMinAlignedRadius()
{
    return min_a_alig;
}

double CDustComponent::getMaxAlignedRadius()
{
    return max_a_alig;
}

double CDustComponent::getScatteringMatrixElement(uint a,
                                    uint w,
                                    uint incID,
                                    uint sphID,
                                    uint sthID,
                                    uint i_mat,
                                    uint j_mat) const
{
    if(sca_mat == 0)
        return 0;
        
    return sca_mat[a][w][incID][sphID][sthID](i_mat, j_mat);
}

const Matrix2D & CDustComponent::getScatteringMatrix(uint a, uint w, uint incID, uint sphID, uint sthID) const
{
    return sca_mat[a][w][incID][sphID][sthID];
}

void CDustComponent::cleanScatteringData()
{
    if(sca_mat != 0)
    {
        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            if(sizeIndexUsed(a))
            {
                for(uint w = 0; w < nr_of_wavelength; w++)
                {
                    if(nr_of_scat_theta[a][w] != 0)
                    {
                        for(uint inc = 0; inc < nr_of_incident_angles; inc++)
                        {
                            for(uint sph = 0; sph < nr_of_scat_phi; sph++)
                                delete[] sca_mat[a][w][inc][sph];
                            delete[] sca_mat[a][w][inc];
                        }
                        delete[] sca_mat[a][w];
                    }
                }
                delete[] sca_mat[a];
            }
        }
        delete[] sca_mat;
        sca_mat = 0;
    }
}

double CDustComponent::getScatteredFraction(uint a, uint w, double theta) const
{
    double res = 1;

    switch(phID)
    {
        case PH_ISO:
            res = 1.0 / PIx4;
            break;

        case PH_HG:
            res = CMathFunctions::phaseFunctionHG(getHGg(a, w), theta);
            break;
        
        case PH_DHG:
            res = CMathFunctions::phaseFunctionDHG(getHGg(a, w), getHGg2(a, w), theta);
            break;
        
        case PH_TTHG:
            res = CMathFunctions::phaseFunctionTTHG(getHGg(a, w), getHGg2(a, w), getHGg3(a, w), theta);
            break;

        case PH_MIE:
            res = getScatteredFractionMie(a, w, theta);
            break;
    }
    return res;
}

double CDustComponent::getScatteredFractionMie(uint a, uint w, double theta) const
{
    return phase_pdf[a][w].getValue(theta);
}

double CDustComponent::getScatteredFractionMie(uint a, uint w, uint sth) const
{
    return phase_pdf[a][w].getValue(sth);
}

void CDustComponent::scatter(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen)
{
    switch(phID)
    {
        case PH_HG:
        {
            uint a = getInteractingDust(grid, pp, rand_gen, CROSS_SCA);
            henyeygreen(pp, a, rand_gen);
            break;
        }

        case PH_DHG:
        {
            uint a = getInteractingDust(grid, pp, rand_gen, CROSS_SCA);
            drainehenyeygreen(pp, a, rand_gen);
            break;
        }

        case PH_TTHG:
        {
            uint a = getInteractingDust(grid, pp, rand_gen, CROSS_SCA);
            threeparamhenyeygreen(pp, a, rand_gen);
            break;
        }

        case PH_MIE:
        {
            uint a = getInteractingDust(grid, pp, rand_gen, CROSS_SCA);
            miesca(pp, a, rand_gen);
            break;
        }

        default:
            pp->setRandomDirection(rand_gen->getRND(), rand_gen->getRND());
            pp->updateCoordSystem();
            break;
    }
}

double CDustComponent::getStochasticHeatingMaxSize() const
{
    return stochastic_heating_max_size;
}

void CDustComponent::setStochasticHeatingMaxSize(double val)
{
    stochastic_heating_max_size = val;
}

double CDustComponent::getFraction() const
{
    return fraction;
}

void CDustComponent::setFraction(double val)
{
    fraction = val;
}

double CDustComponent::getDustMassFraction() const
{
    return dust_mass_fraction;
}

void CDustComponent::setDustMassFraction(double val)
{
    dust_mass_fraction = val;
}

void CDustComponent::setSizeParameter(string size_key, dlist size_parameter_list)
{
    size_keyword = size_key;
    size_parameter = size_parameter_list;
}

string CDustComponent::getDustSizeKeyword()
{
    return size_keyword;
}

string CDustComponent::getDustSizeParameterString()
{
    stringstream str_stream;
    for(uint i = 0; i < NR_OF_SIZE_DIST_PARAM; i++)
        if(size_parameter[i] != 0)
        {
            if(i > 0)
                str_stream << ", ";
            str_stream << size_parameter[i];
        }
        else
        {
            if(i == 0)
                return "flat";
        }
    return str_stream.str();
}

double CDustComponent::convDensityToNumber(CGridBasic * grid, const cell_basic & cell, bool from_gas) const
{
    double avg_mass = getAvgMass(grid, cell);

    if(avg_mass == 0)
        return 0;

    double conversion_factor = 1.0;
    if(from_gas)
    {
        if(!grid->getGasIsMassDensity())
            conversion_factor *= grid->getMu() * m_H;
        conversion_factor *= getDustMassFraction();
        conversion_factor /= avg_mass;
    }
    else if(grid->getDustIsMassDensity())
        conversion_factor = 1 / avg_mass;

    return conversion_factor;
}

double CDustComponent::convDensityToMass(CGridBasic * grid, const cell_basic & cell, bool from_gas) const
{
    double conversion_factor = 1.0;
    if(from_gas)
    {
        if(!grid->getGasIsMassDensity())
            conversion_factor *= grid->getMu() * m_H;
        conversion_factor *= getDustMassFraction();
    }
    else if(!grid->getDustIsMassDensity())
        conversion_factor *= getAvgMass(grid, cell);

    return conversion_factor;
}

double CDustComponent::getNumberDensity(CGridBasic * grid, const cell_basic & cell) const
{
    if(grid->useDustDensities())
        return convDensityToNumber(grid, cell) * grid->getDustDensity(cell);
    else
        return convDensityToNumber(grid, cell, true) * grid->getGasDensity(cell);
}

double CDustComponent::getNumberDensity(CGridBasic * grid, const photon_package & pp) const
{
    return getNumberDensity(grid, *pp.getPositionCell());
}

double CDustComponent::getNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const
{
    if(grid->useDustDensities())
        return convDensityToNumber(grid, cell) * grid->getDustDensity(cell, i_density);
    else
        return convDensityToNumber(grid, cell, true) * grid->getGasDensity(cell, i_density);
}

double CDustComponent::getNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_density) const
{
    return getNumberDensity(grid, *pp.getPositionCell(), i_density);
}

double CDustComponent::getMassDensity(CGridBasic * grid, const cell_basic & cell) const
{
    if(grid->useDustDensities())
        return convDensityToMass(grid, cell) * grid->getDustDensity(cell);
    else
        return getDustMassFraction() * grid->getGasMassDensity(cell);
}

uint CDustComponent::getMassDensity(CGridBasic * grid, const photon_package & pp) const
{
    return getMassDensity(grid, *pp.getPositionCell());
}

double CDustComponent::getMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const
{
    if(grid->useDustDensities())
        return convDensityToMass(grid, cell) * grid->getDustDensity(cell, i_density);
    else
        return getDustMassFraction() * grid->getGasMassDensity(cell, i_density);
}

void CDustComponent::setMu(double mu_)
{
    mu = mu_;
}

void CDustComponent::setMaterialDensity(double dens)
{
    material_density = dens;
}

bool CDustComponent::checkGrainSizeLimits(double a_min, double a_max)
{
    if(a_min == 0)
    {
        // If no minimum grain size is set, start with the first one
        a_min_global = a_eff[0];
    }
    else if(a_min <= a_eff[nr_of_dust_species - 1])
    {
        // Check if inside boundaries
        a_min_global = a_min;
    }
    else
    {
        cout << ERROR_LINE << "Minimum grain size must be smaller than " << a_eff[nr_of_dust_species - 1]
                << "!        " << endl;
        return false;
    }

    if(a_max == 0)
    {
        // If no maximum grain size is set, end with the last one
        a_max_global = a_eff[nr_of_dust_species - 1];
    }
    else if(a_max >= a_eff[0])
    {
        // Check if inside boundaries
        a_max_global = a_max;
    }
    else
    {
        cout << ERROR_LINE << "Maximum grain size must be larger than " << a_eff[0] << "!        " << endl;
        return false;
    }

    // If the minimum grain size is larger than the maximum, no calculation possible
    if(a_min_global > a_max_global)
    {
        cout << ERROR_LINE << "Minimum grain size (" << a_min_global
                << ") must be smaller than maximum grain size (" << a_max_global << ")!" << endl;
        return false;
    }

    return true;
}

double CDustComponent::combinedRFactor(double R1, double R2, double R3) const
{
    double res = (R1 + R2 + R3 + R1 * R2 + R1 * R3 + R2 * R3 + 3 * R1 * R2 * R3) /
                    (1 + 2 * R1 * R2 + 2 * R1 * R3 + 2 * R2 * R3 + 2 * R1 * R2 * R3);

    if(res > 1)
        res = 1;
    if(res < -0.5)
        res = -0.5;

    return res;
}

void CDustComponent::setDelta0(double val)
{
    delta0 = val;
}

void CDustComponent::setLarmF(double val)
{
    larm_f = val;
}

double CDustComponent::getAspectRatio()
{
    return aspect_ratio;
}

double CDustComponent::getSizeParam(CGridBasic * grid, const cell_basic & cell) const
{
    double size_param = grid->getGrainSizeParam(cell, component_id);
    if(size_param != 0)
        return size_param;
    return 0;
}

double CDustComponent::getSizeParam(CGridBasic * grid, const photon_package & pp) const
{
    return getSizeParam(grid, *pp.getPositionCell());
}

double CDustComponent::getSizeParam() const
{
    return 0;
}

double CDustComponent::getSizeMin(CGridBasic * grid, const cell_basic & cell) const
{
    double a_min = grid->getMinGrainRadius(cell, component_id);
    if(a_min > 0)
        return a_min;
    return a_min_global;
}

double CDustComponent::getSizeMin(CGridBasic * grid, const photon_package & pp) const
{
    return getSizeMin(grid, *pp.getPositionCell());
}

double CDustComponent::getSizeMin() const
{
    return a_min_global;
}

void CDustComponent::setSizeMin(double val)
{
    a_min_global = val;
}

double CDustComponent::getSizeMax(CGridBasic * grid, const cell_basic & cell) const
{
    double a_max = grid->getMaxGrainRadius(cell, component_id);
    if(a_max != 0 && a_max < nr_of_dust_species)
        return a_max;
    return a_max_global;
}

double CDustComponent::getSizeMax(CGridBasic * grid, const photon_package & pp) const
{
    return getSizeMax(grid, *pp.getPositionCell());
}

double CDustComponent::getSizeMax() const
{
    return a_max_global;
}

void CDustComponent::setSizeMax(double val)
{
    a_max_global = val;
}

bool CDustComponent::sizeIndexUsed(uint a, double a_min, double a_max) const
{
    if((a_eff[a] - a_min) > -a_min * 1e-5 && (a_max - a_eff[a]) > -a_max * 1e-5)
        return true;
    else
    {
        if(a < nr_of_dust_species - 1)
            if(a_eff[a] < a_min && a_eff[a + 1] >= a_min)
                return true;
        if(a > 0)
            if(a_eff[a - 1] < a_max && a_eff[a] >= a_max)
                return true;
    }
    return false;
}

bool CDustComponent::sizeIndexUsed(uint a) const
{
    if((a_eff[a] - a_min_global) > -a_min_global * 1e-5 &&
        (a_max_global - a_eff[a]) > -a_max_global * 1e-5)
        return true;
    else
    {
        if(a < nr_of_dust_species - 1)
            if(a_eff[a] < a_min_global && a_eff[a + 1] >= a_min_global)
                return true;
        if(a > 0)
            if(a_eff[a - 1] < a_max_global && a_eff[a] >= a_max_global)
                return true;
    }
    return false;
}

double CDustComponent::getInternalIDG(double Td, double Tg) const
{
    double alpha = sqrt((1 + 0.5 / (aspect_ratio * aspect_ratio)) * (1 + Td / Tg));
    double delta = 2 / (1 + aspect_ratio * aspect_ratio) - 1;

    double x = alpha * delta;
    if(x <= 1e-4)
        x = 1e-4;
    double sq_x = sqrt(x);
    double erfi = CMathFunctions::getErfi(sq_x);
    if(erfi <= 1e-4)
        erfi = 1e-4;

    double res = exp(x) / (PIsq * sq_x * erfi) - 1. / (2. * x);
    if(res != res)
        res = 1;
    if(res > 1)
        res = 1;
    if(res < 1e-5)
        res = 1e-5;

    return res;
}

double CDustComponent::getInternalGOLD(double Td, double Tg, double vg) const
{
    double h = 2 / (1 + aspect_ratio * aspect_ratio);
    double alpha = 0.5 / Td * (2 / h + 1) * (0.5 * (Td + Tg) + mu * m_H * vg * vg / (6 * con_kB));
    double delta = h - 1;

    double x = alpha * delta;
    if(x <= 1e-4)
        x = 1e-4;
    double sq_x = sqrt(x);
    double erfi = CMathFunctions::getErfi(sq_x);
    if(erfi <= 1e-4)
        erfi = 1e-4;

    double res = exp(x) / (PIsq * sq_x * erfi) - 1. / (2. * x);

    if(res != res)
        res = 1;
    if(res > 1)
        res = 1;
    if(res < 1e-5)
        res = 1e-5;

    return res;
}

double CDustComponent::getInternalRAT() const
{
    // x = delta, alpha = 1;
    double x = 2. / (1. + aspect_ratio * aspect_ratio) - 1.;
    if(x <= 1e-4)
        x = 1e-4;
    double sq_x = sqrt(x);
    double erfi = CMathFunctions::getErfi(sq_x);
    if(erfi <= 1e-4)
        erfi = 1e-4;

    double res = exp(x) / (PIsq * sq_x * erfi) - 1. / (2. * x);

    if(res != res)
        res = 1;
    if(res > 1)
        res = 1;
    if(res < 1e-5)
        res = 1e-5;

    return res;
}

double CDustComponent::findTemperature(uint a, double qb) const
{
    // Return temperature for energy (qb) value (grain size a)
    double temp = tab_em[a].getValue(qb);

    if(temp < TEMP_MIN)
        temp = double(TEMP_MIN);

    return temp;
}

double CDustComponent::findTemperature(CGridBasic * grid, cell_basic * cell, double qb) const
{
    if(tCabs1 != 0 && tCabs2 != 0)
    {
        // Get temperature from current absorbed energy
        return max(double(TEMP_MIN), tab_em_eff.getValue(qb));
    }

    // Get number of temperatures from tab_temp spline
    uint nr_of_temperatures = tab_temp.size();

    // Init spline for absorption/emission energy interpolation
    spline tmp_tab_em(nr_of_temperatures);

    // Init a temporary array for QB values
    double * tmpQB = new double[nr_of_wavelength];

    // Init a temporary array for planck values
    double * tmpCabs = new double[nr_of_wavelength];

    // init photon package for position and wavelength
    photon_package pp = photon_package();

    // Set position of the photon into current cell
    pp.setPositionCell(cell);

    for(uint w = 0; w < nr_of_wavelength; w++)
    {
        // Set wavelength of photon package
        pp.setWavelength(wavelength_list[w], w);

        // Pre calculate absorption cross-sections
        tmpCabs[w] = getCabsMean(grid, pp);
    }

    for(uint t = 0; t < nr_of_temperatures; t++)
    {
        // Get temperature from tab_temp spline
        double tmp_temp = tab_temp.getValue(t);

        // Calculate absorption cross-section times Planck function for each
        // wavelength
        for(uint w = 0; w < nr_of_wavelength; w++)
            tmpQB[w] = tmpCabs[w] * CMathFunctions::planck(wavelength_list[w], tmp_temp);

        // Calculate QB integrated over all wavelengths
        tmp_tab_em.setValue(
            t, CMathFunctions::integ(wavelength_list, tmpQB, 0, nr_of_wavelength - 1), tmp_temp);
    }

    // Delete multiple pointer
    delete[] tmpCabs;
    delete[] tmpQB;

    // Create spline for interpolation
    tmp_tab_em.createSpline();

    // Return temperature from current absorbed energy
    return max(double(TEMP_MIN), tmp_tab_em.getValue(qb));
}

uint CDustComponent::findTemperatureID(double t) const
{
    return tab_temp.getYIndex(t);
}

double CDustComponent::getPlanck(uint w, double temp) const
{
    double pl = CMathFunctions::planck(wavelength_list[w], temp);
    return max(1e-200, pl);
}

double CDustComponent::getAbsRate(CGridBasic * grid, const cell_basic & cell, uint a, bool use_energy_density) const
{
    double abs_rate = 0;
    if(use_energy_density)
    {
        double * sp_energy = new double[nr_of_wavelength];
        for(uint w = 0; w < nr_of_wavelength; w++)
            sp_energy[w] = (grid->getSpecLength(cell, w) * getCabsMean(a, w));
        abs_rate = CMathFunctions::integ(wavelength_list, sp_energy, 0, nr_of_wavelength - 1);

        delete[] sp_energy;
    }
    else
        for(uint w = 0; w < nr_of_wavelength; w++)
            abs_rate += grid->getSpecLength(cell, w) * getCabsMean(a, w);

    return abs(abs_rate) / (PIx4 * grid->getVolume(cell));
}

double CDustComponent::getAbsRate(CGridBasic * grid, const photon_package & pp, uint a, bool use_energy_density) const
{
    return getAbsRate(grid, *pp.getPositionCell(), a, use_energy_density);
}

void CDustComponent::setFHighJ(double val)
{
    f_highJ = val;
}

void CDustComponent::setFcorr(double val)
{
    f_cor = val;
}

void CDustComponent::setQref(double val)
{
    Q_ref = val;
}

void CDustComponent::setAlphaQ(double val)
{
    alpha_Q = val;
}

void CDustComponent::setRayleighReductionFactor(double val)
{
    R_rayleigh = val;
}

double CDustComponent::getQrat(uint a, uint w, double theta) const
{
    if(theta < 0)
        theta = 0;
    if(theta > PI2)
        theta = PI - theta;
    if(theta > PI2)
        theta = PI2;

    theta += PI2;

    double res = Qtrq[w * nr_of_dust_species + a].getValue(theta);

    if(res <= 0)
        res = getQratApproximation(a, w);
    if(res > 3)
        res = 3;
    return res;
}

double CDustComponent::getQratApproximation(uint a, uint w) const
{
    double ww = wavelength_list[w];
    double aa = getEffectiveRadius(a);

    if(ww / aa < 1)
        return 1;

    double res = pow(ww / aa, -3.0);
    return res;
}

double CDustComponent::getQB(uint a, uint tID) const
{
    return tab_em[a].getX(tID);
}

uint CDustComponent::findWavelengthID(uint a, uint tIDnew, double rnd) const
{
    return avg_planck_frac[tIDnew * nr_of_dust_species + a].getIndex(rnd);
}

double CDustComponent::findTheta(uint a, uint w, double rnd) const
{
    return avg_scattering_frac[a][w].getValue(rnd);
}

uint CDustComponent::findSizeID(photon_package * pp, prob_list & prob, double a_min, double a_max, CRandomGenerator * rand_gen) const
{
    uint a;

    // Ensure that the grain size index is one of the used ones
    while(true)
    {
        // Find the related grain size index
        a = prob.getIndex(rand_gen->getRND());

        // If in case the size index is outside of the used grain sizes, pick another
        // one
        if(sizeIndexUsed(a, a_min, a_max))
            break;
    }

    return a;
}

double CDustComponent::getAvgMass(CGridBasic * grid, const cell_basic & cell) const
{
    if(avg_mass != 0)
        return avg_mass;

    // Get local min and max grain sizes
    double a_min = getSizeMin(grid, cell);
    double a_max = getSizeMax(grid, cell);

    // Get local size parameter for size distribution
    double size_param = getSizeParam(grid, cell);

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    // Create the final mass distribution
    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        rel_weight[a] *= mass[a];
    }

    // Calculate the average mass of the dust grains in the current cell
    double avg_mass =
        CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);

    // Delete pointer array
    delete[] rel_weight;

    return avg_mass;
}

double CDustComponent::getAvgMass() const
{
    // Get integration over the dust size distribution
    double weight = getWeight();

    // Init pointer array of relative mass of each dust grain size bin
    double * rel_mass = new double[nr_of_dust_species];

    // Set the relative mass of each dust grain size bin
    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_mass[a] = grain_size_distribution[a] * mass[a];

    // Calculate the average mass of the dust grains in the current cell
    double avg_mass =
        1 / weight *
        CMathFunctions::integ_dust_size(a_eff, rel_mass, nr_of_dust_species, a_min_global, a_max_global);

    // Delete pointer array
    delete[] rel_mass;

    return avg_mass;
}

double CDustComponent::getPahMass(uint a)
{
    // See Draine and Lee (2000)
    // Infrared Emission from Interstellar Dust. I. Stochastic Heating of Small Grains

    // Calculate the amount of carbon atoms of the PAH molecule
    double nr_carbon_atoms = round(pow(a_eff[a] * 1e9, 3) * 468.0);

    // Calculate the amount of hydrogen atoms of the PAH molecule
    double nr_hydrogen_atoms;
    if(nr_carbon_atoms <= 25)
        nr_hydrogen_atoms = floor(0.5 * nr_carbon_atoms + 0.5);
    else if(nr_carbon_atoms <= 100)
        nr_hydrogen_atoms = floor(2.5 * sqrt(nr_carbon_atoms) + 0.5);
    else
        nr_hydrogen_atoms = floor(0.25 * nr_carbon_atoms + 0.5);

    return (12 * nr_carbon_atoms + nr_hydrogen_atoms) * m_H;
}

double CDustComponent::getEnthalpy(uint a, uint t) const
{
    return enthalpy[a][t];
}

double CDustComponent::getEnthalpyBinWidth(uint a, uint t) const
{
    double diff = 0;
    if(t == nr_of_calorimetry_temperatures - 1)
        diff = enthalpy[a][t] - enthalpy[a][t - 1];
    else if(t == 0)
        diff = enthalpy[a][t + 1] - enthalpy[a][t];
    else
        diff = (enthalpy[a][t + 1] - enthalpy[a][t - 1]) / 2.0;
    return diff;
}

double CDustComponent::getCalorimetricTemperature(uint i) const
{
    return calorimetry_temperatures[i];
}

double CDustComponent::getMinCalorimetricTemperature() const
{
    return calorimetry_temperatures[0];
}

double CDustComponent::getMaxCalorimetricTemperature() const
{
    return calorimetry_temperatures[nr_of_calorimetry_temperatures - 1];
}

uint CDustComponent::getNrOfCalorimetryTemperatures() const
{
    return nr_of_calorimetry_temperatures;
}

uint CDustComponent::getNrOfDustSpecies() const
{
    return nr_of_dust_species;
}

void CDustComponent::setNrOfDustSpecies(uint val)
{
    nr_of_dust_species = val;
}

void CDustComponent::setNrOfWavelength(uint val)
{
    nr_of_wavelength = val;
}

uint CDustComponent::getNrOfWavelength() const
{
    return nr_of_wavelength;
}

uint CDustComponent::getWavelengthID(double wavelength)
{
    dlist::iterator it = find(wavelength_list.begin(), wavelength_list.end(), wavelength);
    if(it != wavelength_list.end())
        return distance(wavelength_list.begin(), it);

    cout << WARNING_LINE << "Wavelength not found!" << endl;
    return 0;
}

uint CDustComponent::getNrOfIncidentAngles() const
{
    return nr_of_incident_angles;
}

uint CDustComponent::getNrOfScatTheta(uint a, uint w) const
{
    return nr_of_scat_theta[a][w];
}

uint ** CDustComponent::getNrOfScatTheta()
{
    return nr_of_scat_theta;
}

uint CDustComponent::getNrOfScatPhi() const
{
    return nr_of_scat_phi;
}

uint CDustComponent::getNrOfScatMatElements() const
{
    return nr_of_scat_mat_elements;
}

string CDustComponent::getStringID() const
{
    return stringID;
}

void CDustComponent::createStringID(CDustComponent * comp)
{
    // Init string stream
    stringstream str_stream;
    str_stream.str("");

    // Start with printing the dust components (if mixture or only single component)
    if(stringID.length() == 0 || comp->getNrOfComponents() == 1)
        str_stream << "Dust components:" << endl;

    // Use the fraction as dust to gas mass ratio (if user chosen it)
    double fraction = comp->getFraction();
    string fraction_string = "mass ratio";
    if(individual_dust_fractions)
    {
        fraction *= dust_mass_fraction;
        fraction_string = "dust-to-gas mass ratio";
    }

    // Fill the string with various parameters and format it
    char tmp_str[1024];
#ifdef WINDOWS
    sprintf_s(tmp_str,
                "- %s\n    %s: %g, size distr. : \"%s\" (%s), size: %g [m] - %g [m]\n",
                comp->getStringID().c_str(),
                fraction_string.c_str(),
                fraction,
                comp->getDustSizeKeyword().c_str(),
                comp->getDustSizeParameterString().c_str(),
                comp->getSizeMin(),
                comp->getSizeMax());
#else
    sprintf(tmp_str,
            "- %s\n    %s: %g, size distr. : \"%s\" (%s), size: %g [m] - %g [m]\n",
            comp->getStringID().c_str(),
            fraction_string.c_str(),
            fraction,
            comp->getDustSizeKeyword().c_str(),
            comp->getDustSizeParameterString().c_str(),
            comp->getSizeMin(),
            comp->getSizeMax());
#endif

    // Add formatted string to stream
    str_stream << tmp_str;

    // Add stream to stringID or replace it if only one component
    if(comp->getNrOfComponents() == 1)
        stringID = str_stream.str();
    else
        stringID += str_stream.str();
}

void CDustComponent::setIndividualDustMassFractions(bool val)
{
    individual_dust_fractions = val;
}

bool CDustComponent::getIndividualDustMassFractions()
{
    return individual_dust_fractions;
}

void CDustComponent::setWavelengthList(dlist _wavelength_list, uint _wavelength_offset)
{
    wavelength_list = _wavelength_list;
    nr_of_wavelength = wavelength_list.size();
    wavelength_offset = _wavelength_offset;
}

bool CDustComponent::calcWavelengthDiff()
{
    // Set width of each wavelength bin
    wavelength_diff.resize(nr_of_wavelength);

    // Differences need at least two wavelengths
    if(nr_of_wavelength == 1)
        return false;

    for(uint w = 0; w < nr_of_wavelength; w++)
    {
        if(w == nr_of_wavelength - 1)
            wavelength_diff[w] = wavelength_list[w] - wavelength_list[w - 1];
        else if(w == 0)
            wavelength_diff[w] = wavelength_list[w + 1] - wavelength_list[w];
        else
            wavelength_diff[w] = (wavelength_list[w + 1] - wavelength_list[w - 1]) / 2.0;
    }
    return true;
}

void CDustComponent::updateStokesVector(photon_package * pp, uint wnew) const
{
    // Get wavelength of photon package
    uint w = pp->getDustWavelengthID();
    *pp->getStokesVector() *= wavelength_diff[w] / wavelength_diff[wnew];
}

double CDustComponent::getEffectiveRadius(uint a) const
{
    return a_eff[a];
}

double * CDustComponent::getEffectiveRadii()
{
    return a_eff;
}

double CDustComponent::getGrainDistributionXRadiusSq(uint a) const
{
    return grain_distribution_x_aeff_sq[a];
}

double CDustComponent::getGrainSizeDistribution(uint a) const
{
    return grain_size_distribution[a];
}

double * CDustComponent::getGrainSizeDistribution()
{
    return grain_size_distribution;
}

double CDustComponent::getEffectiveRadiusSquared(uint a) const
{
    return a_eff_squared[a];
}

double CDustComponent::getMass(uint a) const
{
    return mass[a];
}

double CDustComponent::getVolume(uint a) const
{
    return 4.0 / 3.0 * PI * a_eff[a] * a_eff[a] * a_eff[a];
}

double CDustComponent::getMaterialDensity(uint a) const
{
    return mass[a] / getVolume(a);
}

double CDustComponent::getMaterialDensity() const
{
    // Get local min and max grain sizes
    double a_min = getSizeMin();
    double a_max = getSizeMax();

    // Get local size parameter for size distribution
    double size_param = getSizeParam();

    // Get integration over the dust size distribution
    double * rel_weight = getRelWeight(a_min, a_max, size_param);

    for(uint a = 0; a < nr_of_dust_species; a++)
        rel_weight[a] *= mass[a] / getVolume(a);

    double res = CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);

    delete[] rel_weight;

    return res;
}

double CDustComponent::getFHighJ() const
{
    return f_highJ;
}

double CDustComponent::getFcorr() const
{
    return f_cor;
}

bool CDustComponent::isAligned() const
{
    return is_align;
}

void CDustComponent::setIsAligned(bool val)
{
    is_align = val;
}

void CDustComponent::setIsMixture(bool val)
{
    is_mixture = val;
}

double CDustComponent::getSublimationTemperature()
{
    return sub_temp;
}

void CDustComponent::getQtrq(uint i, uint j, double & x, double & y)
{
    x = Qtrq[i].getX(j);
    y = Qtrq[i].getY(j);
}

void CDustComponent::getHG_g_factor(uint i, uint j, double & x, double & y)
{
    x = HG_g_factor[i].getX(j);
    y = HG_g_factor[i].getY(j);
}

void CDustComponent::getHG_g2_factor(uint i, uint j, double & x, double & y)
{
    x = HG_g2_factor[i].getX(j);
    y = HG_g2_factor[i].getY(j);
}

void CDustComponent::getHG_g3_factor(uint i, uint j, double & x, double & y)
{
    x = HG_g3_factor[i].getX(j);
    y = HG_g3_factor[i].getY(j);
}

bool CDustComponent::getCalorimetryLoaded()
{
    return calorimetry_loaded;
}

void CDustComponent::setCalorimetryLoaded(bool val)
{
    calorimetry_loaded = val;
}

void CDustComponent::setIDs(uint i_comp, uint nr_components, uint i_mix, uint nr_mixtures)
{
    i_component = i_comp;
    nr_of_components = nr_components;

    i_mixture = i_mix;
    nr_of_mixtures = nr_mixtures;
}

void CDustComponent::printIDs()
{
    cout << "-> Dust mixture " << i_mixture + 1 << "/" << nr_of_mixtures << ", ";
    if(!is_mixture)
        cout << "component " << i_component + 1 << "/" << nr_of_components << " ";
}

bool CDustComponent::getScatLoaded()
{
    return scat_loaded;
}

void CDustComponent::setScatLoaded(bool val)
{
    scat_loaded = val;
}

void CDustComponent::setSublimate(bool val)
{
    sublimate = val;
}

uint CDustComponent::getComponentId()
{
    return i_component;
}

uint CDustComponent::getNrOfComponents()
{
    return nr_of_components;
}

void CDustComponent::SetNrOfScatTheta(uint ** nr_of_scat_theta_tmp)
{
    if(nr_of_scat_theta != 0)
    {
        for(uint a = 0; a < nr_of_dust_species; a++)
            delete[] nr_of_scat_theta[a];
        delete[] nr_of_scat_theta;
    }
    nr_of_scat_theta = nr_of_scat_theta_tmp;
}

void CDustComponent::SetScatTheta(double *** scat_theta_tmp)
{
    if(scat_theta != 0)
    {
        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            for(uint w = 0; w < nr_of_wavelength; w++)
                delete[] scat_theta[a][w];
            delete[] scat_theta[a];
        }
        delete[] scat_theta;
    }
    scat_theta = scat_theta_tmp;
}

double * CDustComponent::getScatTheta(uint a, uint w)
{
    return scat_theta[a][w];
}

double CDustComponent::getScatTheta(uint a, uint w, uint sth) const
{
    return scat_theta[a][w][sth];
}

uint CDustComponent::getScatThetaID(double theta, uint a, uint w) const
{
    // Returns the index of the value in scat_theta[a][w] closest to theta
    uint sth;

    // If no refinement was done, don't use binary
    // search to save (lots of) time
    if( (2*NANG-1) == nr_of_scat_theta[a][w])
        sth = uint(theta/PI * 2*(NANG-1) + 0.5);
    else
    {
        if(theta == scat_theta[a][w][0])
            return 0;

        sth = CMathFunctions::biListIndexSearch(theta,scat_theta[a][w],nr_of_scat_theta[a][w]);
        if(sth != MAX_UINT && (scat_theta[a][w][sth+1] + scat_theta[a][w][sth]) / 2. <= theta )
            sth++;
    }

    return sth;
}
