/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include <cmath>
#include "DustMixture.hpp"
#include "CommandParser.hpp"
#include "GridBasic.hpp"
#include "MathFunctions.hpp"
#include "Typedefs.hpp"
#include "Parameters.hpp"

bool CDustMixture::createDustMixtures(parameters & param, string path_data, string path_plot)
{
    // Do not load dust component if not required
    uint nr_of_total_components = param.getTotalNrOfDustComponents();
    if(nr_of_total_components == 0)
    {
        if(param.getCommand() == CMD_LINE_EMISSION || param.getCommand() == CMD_OPIATE ||
           param.getCommand() == CMD_SYNCHROTRON || param.getCommand() == CMD_FREE_FREE)
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
            single_component[i_comp].setFcorr(param.getFcorr());
            single_component[i_comp].setQref(param.getQref());
            single_component[i_comp].setAlphaQ(param.getAlphaQ());
            single_component[i_comp].setRayleighReductionFactor(param.getRayleighReductionFactor());

            single_component[i_comp].setDelta0(param.getDelta0());
            single_component[i_comp].setLarmF(param.getLarmF());
            single_component[i_comp].setMu(param.getMu());
            // single_component[i_comp].setPhaseFunctionID(param.getPhaseFunctionID());
            single_component[i_comp].setPhaseFunctionID(param.getPhaseFunctionID(current_dust_choice));

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
                    cout << ERROR_LINE << "Cannot open calorimetry file, which is required "
                            "for stochastic heating!"
                         << endl;
                    return false;
                }

            // Write dust component files, if multiple components will be mixed together
            if(nr_of_components > 1)
            {
                if(param.getWriteDustFiles())
                    single_component[i_comp].writeComponentPlot(path_plot);
                single_component[i_comp].writeComponentData(path_data);
            }
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

void CDustMixture::printParameters(parameters & param, CGridBasic * grid)
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
    cout << SEP_LINE;
    cout << "Dust parameters                                                             "
            "            "
         << endl;
    cout << SEP_LINE;

    // Show the full wavelength grid for TEMP and RAT calculation
    if(param.isTemperatureSimulation() || param.isRatSimulation())
        cout << "- Number of wavelengths   : " << WL_STEPS << "  (" << WL_MIN << " [m] - " << WL_MAX
             << " [m])" << endl;

    // Monte-Carlo scattering is only used for temp, rat and scatter maps
    // if(param.isMonteCarloSimulation() || scattering_to_raytracing)
    //     cout << "- Phase function          : " << getPhaseFunctionStr() << endl;

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
                cout << INFO_LINE << "Stochastic heating was already calculated. This "
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
                cout << ERROR_LINE << "This should not happen!" << endl;
        }

        cout << "- Include scattered light : ";
        if(grid->isRadiationFieldAvailable())
        {
            if(scattering_to_raytracing)
                cout << "yes, based on the radiation field" << endl;
            else
                cout << "no, disabled via <rt_scattering> 0" << endl;
        }
        else
        {
            if(scattering_to_raytracing)
                cout << "yes, radiation field will be calculated before raytracing" << endl;
            else if(!param.getScatteringToRay())
                cout << "no, disabled via <rt_scattering> 0" << endl;
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

        if(param.splitDustEmission())
        {
            cout << "- Split dust emission     : "
                 << "yes, as additional entries in fits file" << endl;
        }

        cout << "Observed wavelengths:" << endl;
        dlist dust_ray_detectors = param.getDustRayDetectors();
        for(uint i = 0; i < dust_ray_detectors.size(); i += NR_OF_RAY_DET)
        {
            // Index of current detector
            uint pos = i / NR_OF_RAY_DET;

            if(uint(dust_ray_detectors[i + 2]) > 1)
                {if(USE_LOG_SPACING)
                    {cout << "- Emission detetector " << (pos + 1) << "   : from " << dust_ray_detectors[i + 0]
                     << " [m] to " << dust_ray_detectors[i + 1] << " [m] with "
                     << uint(dust_ray_detectors[i + 2]) << " logarithmic values" << endl;}
                else
                    {cout << "- Emission detetector " << (pos + 1) << "   : from " << dust_ray_detectors[i + 0]
                     << " [m] to " << dust_ray_detectors[i + 1] << " [m] with "
                     << uint(dust_ray_detectors[i + 2]) << " linear values" << endl;}
                }
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
                {if(USE_LOG_SPACING)
                    {cout << "- Scattering detetector " << (pos + 1) << " : from " << dust_mc_detectors[i + 0]
                     << " [m]) to " << dust_mc_detectors[i + 1] << " [m]) with "
                     << uint(dust_mc_detectors[i + 2]) << " logarithmic values" << endl;}
                else
                    {cout << "- Scattering detetector " << (pos + 1) << " : from " << dust_mc_detectors[i + 0]
                     << " [m]) to " << dust_mc_detectors[i + 1] << " [m]) with "
                     << uint(dust_mc_detectors[i + 2]) << " linear values" << endl;}
                }
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
            cout << INFO_LINE << "Stochastic heating and saving the radiation field is chosen." << endl;
            cout << "  The radiation field will be saved and stochastic heating should be set" << endl;
            cout << "  with CMD_DUST_EMISSION simulation." << endl;
        }
    }

    // Show information about saving the radiation field
    if(param.isTemperatureSimulation() || param.isRatSimulation())
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

        cout << "- Phase function          : " << getPhaseFunctionStr(i_mixture) << endl;
        cout << "- Avg. grain mass         : " << getAvgMass(i_mixture) << " [kg]" << endl;

        if(param.getCommand() == CMD_DUST_EMISSION && !param.getAligRANDOM())
        {
            if(mixed_component[i_mixture].isAligned())
                cout << "- Affected by alignment   : Yes" << endl;
            else
                cout << "- Affected by alignment   : No" << endl;
        }

        double total_dust_mass = 0;
        for(ulong i_cell = 0; i_cell < grid->getMaxDataCells(); i_cell++)
        {
            cell_basic * cell = grid->getCellFromIndex(i_cell);
            total_dust_mass += getMassDensity(grid, *cell, i_mixture) * grid->getVolume(*cell);
        }
        
        ulong max_cells = grid->getMaxDataCells();
        uint marked_cells=mixed_component[i_mixture].getNrMarked(max_cells);
        
        if(marked_cells>0)
            cout << "- Nr. of sub. marker      : " << marked_cells << endl;
        
        cout << "- Total mass              : " << total_dust_mass / M_sun << " [M_sun], " << total_dust_mass
             << " [kg]" << endl;
        cout << mixed_component[i_mixture].getStringID();
    }
    cout << SEP_LINE;

    // If dust optical properties are constant -> pre calc some values
    if(grid->useConstantGrainSizes())
    {
        preCalcEffProperties(param);
        #if (USE_PRECALC_TABLE)
            preCalcRelWeight();
        #endif
    }
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
    nr_of_dust_species = single_component[0].getNrOfDustSpecies();
    uint nr_of_incident_angles = single_component[0].getNrOfIncidentAngles();

    for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
    {
        // Check if the components have the same amount of grain sizes
        if(nr_of_dust_species != single_component[i_comp].getNrOfDustSpecies())
        {
            cout << ERROR_LINE << "Component Nr. " << i_comp + 1 << " has a different amount of dust species!"
                 << endl;
            return false;
        }

        // Check if the components have the same amount of incident angles
        if(nr_of_incident_angles != single_component[i_comp].getNrOfIncidentAngles())
        {
            cout << ERROR_LINE << "Component Nr. " << i_comp + 1 << " has a different amount of incident angles!"
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
    }

    // Get all scattering thetas of the individual components and count all unique values for the mixed_component
    uint ** nr_of_scat_theta;
    double *** scat_theta;

    if(mixed_component[i_mixture].getScatLoaded())
        getNrOfUniqueScatTheta(nr_of_scat_theta, scat_theta);

    for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
    {
        // Add parameters of each component together to the mixture
        if(!mixed_component[i_mixture].add(size_fraction[i_comp], &single_component[i_comp], nr_of_scat_theta, scat_theta))
            return false;
    }

    // Pre-calculate various quantities
    if(!preCalcDustProperties(param, i_mixture))
        return false;

    // Delete pointer array
    for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
    {
        for(uint a = 0; a < nr_of_dust_species; a++)
            if(mixed_component[i_mixture].sizeIndexUsed(a))
                delete[] size_fraction[i_comp][a];
        delete[] size_fraction[i_comp];
    }
    delete[] size_fraction;

    return true;
}

void CDustMixture::getNrOfUniqueScatTheta(uint ** & nr_of_scat_theta, double *** & scat_theta)
{
    nr_of_scat_theta = new uint *[nr_of_dust_species];
    scat_theta = new double **[nr_of_dust_species];
    for(uint a=0; a < nr_of_dust_species; a++)
    {
        nr_of_scat_theta[a] = new uint[nr_of_wavelength];
        scat_theta[a] = new double *[nr_of_wavelength];
        for(uint w=0; w < nr_of_wavelength; w++)
        {
            dlist scat_theta_tmp;
            for(uint i_comp=0; i_comp < nr_of_components; i_comp++)
                scat_theta_tmp.insert(scat_theta_tmp.end(),
                                      single_component[i_comp].getScatTheta(a,w),
                                      single_component[i_comp].getScatTheta(a,w) + single_component[i_comp].getNrOfScatTheta(a,w));

            sort(scat_theta_tmp.begin(),scat_theta_tmp.end());
            scat_theta_tmp.erase( unique(scat_theta_tmp.begin(),scat_theta_tmp.end()), scat_theta_tmp.end() );

            nr_of_scat_theta[a][w] = scat_theta_tmp.size();

            scat_theta[a][w] = new double[nr_of_scat_theta[a][w]];
            copy(scat_theta_tmp.begin(),scat_theta_tmp.end(),scat_theta[a][w]);

            scat_theta_tmp.clear();
        }
    }
}

bool CDustMixture::preCalcDustProperties(parameters & param, uint i_mixture)
{
    // Set various parameters for the mixture
    mixed_component[i_mixture].setSubStatus(param.getSubStatus());

    // Calculate wavelength differences for temperature (reemission)
    if(param.isTemperatureSimulation() || param.isRatSimulation())
        if(!mixed_component[i_mixture].calcWavelengthDiff())
        {
            cout << ERROR_LINE << "The wavelength grid only has one wavelength which is not enough for temp calculation!\n"
                 << "       Change WL_STEPS to more than one in the src/Typedefs.h and recompile!" << endl;
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
    if(param.isTemperatureSimulation() || param.isRatSimulation())
        mixed_component[i_mixture].preCalcWaveProb();

    // Pre-calculate temperature to total emission relation (either for temperature or
    // stochastic heating calculation)
    if((param.isTemperatureSimulation() || param.isRatSimulation()) ||
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
    if(param.isMonteCarloSimulation())
        mixed_component[i_mixture].setAlignmentMechanism(ALIG_RND);
    else
        mixed_component[i_mixture].setAlignmentMechanism(param.getAlignmentMechanism());

    return true;
}

uint CDustMixture::getMixtureID(CGridBasic * grid, const cell_basic & cell) const
{
    uint dust_choice = grid->getDustChoiceID(cell);
    return dust_choices_to_index[dust_choice];
}

uint CDustMixture::getMixtureID(CGridBasic * grid, const photon_package & pp) const
{
    uint dust_choice = grid->getDustChoiceID(pp);
    return dust_choices_to_index[dust_choice];
}

double CDustMixture::getAvgMass(uint i_mixture)
{
    return mixed_component[i_mixture].getAvgMass();
}

bool CDustMixture::writeComponentData(string path_data)
{
    if(mixed_component != 0)
        for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            if(!mixed_component[i_mixture].writeComponentData(path_data))
                return false;
    return true;
}

bool CDustMixture::writeComponentPlot(string path_plot)
{
    if(mixed_component != 0)
        for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            if(!mixed_component[i_mixture].writeComponentPlot(path_plot))
                return false;
    return true;
}

void CDustMixture::preCalcEffProperties(parameters & param)
{
    if(mixed_component != 0)
        for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            mixed_component[i_mixture].preCalcEffProperties(param);
}

void CDustMixture::preCalcRelWeight()
{
    if(mixed_component != 0)
        for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            mixed_component[i_mixture].preCalcRelWeight();
}

string CDustMixture::getPhaseFunctionStr(uint i_mixture)
{
    string str_res = ERROR_LINE;
    str_res += "Phase function is undefined!";
    if(mixed_component != 0)
    {
        switch(mixed_component[i_mixture].getPhaseFunctionID())
        {
            case PH_ISO:
                str_res = "Isotropic scattering";
                break;

            case PH_HG:
                str_res = "Henyey and Greenstein";
                break;
            
            case PH_DHG:
                str_res = "Draine Henyey and Greenstein";
                break;
            
            case PH_TTHG:
                str_res = "Three parameter Henyey and Greenstein";
                break;

            case PH_MIE:
                str_res = "Mie scattering";
                break;
        }
    }
    return str_res;
}

void CDustMixture::killSingleComponents()
{
    // Show information
    cout << "                \r";
    cout << "-> Dust cleanup started            \r";

    // If single component was used -> delete it
    if(single_component != 0)
    {
        // Delete the single component
        delete[] single_component;
        single_component = 0;
    }
}

uint CDustMixture::getNrOfDustSpecies(uint i_mixture)
{
    return mixed_component[i_mixture].getNrOfDustSpecies();
}

uint CDustMixture::getNrOfStochasticSizes(uint i_mixture)
{
    return mixed_component[i_mixture].getNrOfStochasticSizes();
}

void CDustMixture::calcTemperature(CGridBasic * grid, cell_basic * cell, bool use_energy_density)
{
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, *cell);
            mixed_component[i_mixture].calcTemperature(grid, cell, 0, use_energy_density);
        }
        else
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                mixed_component[i_mixture].calcTemperature(grid, cell, i_mixture, use_energy_density);
    }
}

void CDustMixture::calcStochasticHeatingPropabilities(CGridBasic * grid,
                                        cell_basic * cell,
                                        dlist & wl_list) const
{
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, *cell);
            mixed_component[i_mixture].calcStochasticHeatingPropabilities(
                grid, cell, 0, wl_list);
        }
        else
        {
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                mixed_component[i_mixture].calcStochasticHeatingPropabilities(
                    grid, cell, i_mixture, wl_list);
        }
    }
}

double CDustMixture::getForegroundExtinction(double wavelength)
{
    if(extinction_magnitude == 0)
        return 1;

    // Get extinction optical depth at extinction_magnitude_wavelength
    double extinction_tau = (extinction_magnitude / 1.086);

    if(extinction_dust_mixture == MAX_UINT)
    {
        // Init extinction curve
        spline extinction_curve(100);

        // Extinction curve for ISM grains (silicate and graphite, 5nm to 250nm)
        extinction_curve.setValue(0, 5.000000074505806e-08, 4.216777503871183e-016);
        extinction_curve.setValue(1, 5.564875528216362e-08, 4.621585931369242e-016);
        extinction_curve.setValue(2, 6.193568184971808e-08, 5.721308394426405e-016);
        extinction_curve.setValue(3, 6.893286854028702e-08, 7.124210419085484e-016);
        extinction_curve.setValue(4, 7.672057300806045e-08, 7.674365729865441e-016);
        extinction_curve.setValue(5, 8.538808673620224e-08, 7.177044099170471e-016);
        extinction_curve.setValue(6, 9.503481537103653e-08, 6.998261913784690e-016);
        extinction_curve.setValue(7, 1.05771400034428e-07, 6.237277455455178e-016);
        extinction_curve.setValue(8, 1.1772090196609499e-07, 5.435211899183833e-016);
        extinction_curve.setValue(9, 1.3102050125598898e-07, 4.750205237431847e-016);
        extinction_curve.setValue(10, 1.45822495222092e-07, 4.006340109059733e-016);
        extinction_curve.setValue(11, 1.62296801805496e-07, 3.590026327160000e-016);
        extinction_curve.setValue(12, 1.8063229322433499e-07, 3.539562486442753e-016);
        extinction_curve.setValue(13, 2.01039299368858e-07, 3.874130706789749e-016);
        extinction_curve.setValue(14, 2.2375169396400497e-07, 3.996385089649330e-016);
        extinction_curve.setValue(15, 2.49030098319054e-07, 3.285992056374746e-016);
        extinction_curve.setValue(16, 2.77164310216904e-07, 2.829757611670918e-016);
        extinction_curve.setValue(17, 3.0847701430320697e-07, 2.536181550806294e-016);
        extinction_curve.setValue(18, 3.43327194452286e-07, 2.287595115539490e-016);
        extinction_curve.setValue(19, 3.82114589214325e-07, 2.065593074020058e-016);
        extinction_curve.setValue(20, 4.25284087657928e-07, 1.845888079974461e-016);
        extinction_curve.setValue(21, 4.7333058714866596e-07, 1.634203032101485e-016);
        extinction_curve.setValue(22, 5.26805222034454e-07, 1.431729732227695e-016);
        extinction_curve.setValue(23, 5.86320996284485e-07, 1.250360258234746e-016);
        extinction_curve.setValue(24, 6.525607109069819e-07, 1.078952024050406e-016);
        extinction_curve.setValue(25, 7.2628378868103e-07, 9.266067295651665e-017);
        extinction_curve.setValue(26, 8.08335781097412e-07, 7.881166441657414e-017);
        extinction_curve.setValue(27, 8.99657726287842e-07, 6.579841269777561e-017);
        extinction_curve.setValue(28, 1.00129699707031e-06, 5.458088737076844e-017);
        extinction_curve.setValue(29, 1.1144180297851599e-06, 4.498859327641566e-017);
        extinction_curve.setValue(30, 1.2403199672698998e-06, 3.660451227198234e-017);
        extinction_curve.setValue(31, 1.3804450035095198e-06, 2.970563074866552e-017);
        extinction_curve.setValue(32, 1.5364010334014899e-06, 2.401290428990413e-017);
        extinction_curve.setValue(33, 1.70997595787048e-06, 1.921721053770342e-017);
        extinction_curve.setValue(34, 1.90316104888916e-06, 1.523635586331645e-017);
        extinction_curve.setValue(35, 2.11817002296448e-06, 1.204889560936278e-017);
        extinction_curve.setValue(36, 2.35747098922729e-06, 9.562793524986687e-018);
        extinction_curve.setValue(37, 2.6238059997558596e-06, 7.653670597381650e-018);
        extinction_curve.setValue(38, 2.92023110389709e-06, 6.185056501137134e-018);
        extinction_curve.setValue(39, 3.25014495849609e-06, 5.054911063069209e-018);
        extinction_curve.setValue(40, 3.61733102798462e-06, 4.168021658769960e-018);
        extinction_curve.setValue(41, 4.02599906921387e-06, 3.468548840603535e-018);
        extinction_curve.setValue(42, 4.4808359146118195e-06, 2.910007332891138e-018);
        extinction_curve.setValue(43, 4.98706007003784e-06, 2.474287670482146e-018);
        extinction_curve.setValue(44, 5.5504732131958e-06, 2.148568553881987e-018);
        extinction_curve.setValue(45, 6.17753791809082e-06, 1.926424030502043e-018);
        extinction_curve.setValue(46, 6.875445842742919e-06, 1.737387824729840e-018);
        extinction_curve.setValue(47, 7.652201175689699e-06, 2.102465435646191e-018);
        extinction_curve.setValue(48, 8.51670932769775e-06, 6.233472209570475e-018);
        extinction_curve.setValue(49, 9.47888565063477e-06, 1.116189825799359e-017);
        extinction_curve.setValue(50, 1.05497598648071e-05, 8.459390008941990e-018);
        extinction_curve.setValue(51, 1.1741620063781699e-05, 5.243048743901026e-018);
        extinction_curve.setValue(52, 1.3068140029907199e-05, 3.222513686771834e-018);
        extinction_curve.setValue(53, 1.45445098876953e-05, 2.593459581332661e-018);
        extinction_curve.setValue(54, 1.6187679290771496e-05, 3.364742430016652e-018);
        extinction_curve.setValue(55, 1.80164794921875e-05, 4.028119360997947e-018);
        extinction_curve.setValue(56, 2.00519008636475e-05, 3.589028317207343e-018);
        extinction_curve.setValue(57, 2.2317260742187498e-05, 2.919286851797911e-018);
        extinction_curve.setValue(58, 2.48385601043701e-05, 2.424045340461135e-018);
        extinction_curve.setValue(59, 2.76446990966797e-05, 2.034588235351009e-018);
        extinction_curve.setValue(60, 3.0767860412597695e-05, 1.725524296987028e-018);
        extinction_curve.setValue(61, 3.42438583374023e-05, 1.464357218017165e-018);
        extinction_curve.setValue(62, 3.8112560272216795e-05, 1.239701787887817e-018);
        extinction_curve.setValue(63, 4.2418331146240196e-05, 1.028703213069295e-018);
        extinction_curve.setValue(64, 4.72105484008789e-05, 8.412356002324342e-019);
        extinction_curve.setValue(65, 5.25441703796387e-05, 6.810797851620609e-019);
        extinction_curve.setValue(66, 5.84803504943848e-05, 5.471715483013374e-019);
        extinction_curve.setValue(67, 6.50871810913086e-05, 4.374078863404243e-019);
        extinction_curve.setValue(68, 7.24404067993164e-05, 3.489064367360922e-019);
        extinction_curve.setValue(69, 8.06243667602539e-05, 2.781946435780552e-019);
        extinction_curve.setValue(70, 8.973293304443359e-05, 2.221216725155223e-019);
        extinction_curve.setValue(71, 9.987050628662108e-05, 1.777522731448490e-019);
        extinction_curve.setValue(72, 0.00011115339660644499, 1.426858510061328e-019);
        extinction_curve.setValue(73, 0.000123710998535156, 1.155888746296647e-019);
        extinction_curve.setValue(74, 0.000137687194824219, 9.449756554245537e-020);
        extinction_curve.setValue(75, 0.000153242401123047, 7.584280096005196e-020);
        extinction_curve.setValue(76, 0.000170554992675781, 6.018818168528717e-020);
        extinction_curve.setValue(77, 0.000189823501586914, 4.803205495665578e-020);
        extinction_curve.setValue(78, 0.00021126879882812498, 3.850335435269971e-020);
        extinction_curve.setValue(79, 0.000235136993408203, 3.094874612999559e-020);
        extinction_curve.setValue(80, 0.00026170159912109397, 2.491039874419238e-020);
        extinction_curve.setValue(81, 0.00029126739501953103, 2.006959239669001e-020);
        extinction_curve.setValue(82, 0.00032417330932617196, 1.618343487101766e-020);
        extinction_curve.setValue(83, 0.00036079681396484396, 1.305601603235169e-020);
        extinction_curve.setValue(84, 0.000401557891845703, 1.053083960958686e-020);
        extinction_curve.setValue(85, 0.00044692401123046894, 8.500255800817664e-021);
        extinction_curve.setValue(86, 0.000497415313720703, 6.863925335919417e-021);
        extinction_curve.setValue(87, 0.000553610778808594, 5.541157311593072e-021);
        extinction_curve.setValue(88, 0.000616155029296875, 4.473916606814938e-021);
        extinction_curve.setValue(89, 0.000685765197753906, 3.612131329544281e-021);
        extinction_curve.setValue(90, 0.0007632396240234379, 2.917011636593121e-021);
        extinction_curve.setValue(91, 0.0008494666748046879, 2.355426278709940e-021);
        extinction_curve.setValue(92, 0.000945435302734375, 1.902069506823009e-021);
        extinction_curve.setValue(93, 0.0010522459716796899, 1.538983867990721e-021);
        extinction_curve.setValue(94, 0.0011711240234375, 1.232992068728240e-021);
        extinction_curve.setValue(95, 0.0013034310302734399, 9.975786070486010e-022);
        extinction_curve.setValue(96, 0.0014506870117187499, 8.016047070202548e-022);
        extinction_curve.setValue(97, 0.0016145780029296899, 6.504697113317689e-022);
        extinction_curve.setValue(98, 0.0017969849853515598, 5.270413281701706e-022);
        extinction_curve.setValue(99, 0.002, 4.197725271491968e-022);

        // Init spline
        extinction_curve.createSpline();

        // Scaling factor
        double scaling_factor =
            extinction_tau / extinction_curve.getValue(extinction_magnitude_wavelength);

        return exp(-scaling_factor * extinction_curve.getValue(wavelength));
    }
    else
    {
        // Scaling factor
        double scaling_factor = extinction_tau / mixed_component[extinction_dust_mixture].getCextMean(
                                                        extinction_magnitude_wavelength);

        return exp(-scaling_factor * mixed_component[extinction_dust_mixture].getCextMean(wavelength));
    }
}

double CDustMixture::getCextMean(CGridBasic * grid, const photon_package & pp) const
{
    double sum = 0;
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, pp);
            sum = mixed_component[i_mixture].getCextMean(grid, pp);
        }
        else
        {
            if(grid->getCextMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID()) != MAX_DOUBLE)
                sum = grid->getCextMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID());
            else
            {
                double dens = 0;
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                {
                    double i_dens = getNumberDensity(grid, pp, i_mixture);
                    dens += i_dens;
                    sum += mixed_component[i_mixture].getCextMean(grid, pp) * i_dens;
                }
                if(dens != 0)
                    sum /= dens;
            }
        }
    }
    return sum;
}

double CDustMixture::getCabsMean(CGridBasic * grid, const photon_package & pp) const
{
    double sum = 0;
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, pp);
            sum = mixed_component[i_mixture].getCabsMean(grid, pp);
        }
        else
        {
            if(grid->getCabsMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID()) != MAX_DOUBLE)
                sum = grid->getCabsMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID());
            else
            {
                double dens = 0;
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                {
                    double i_dens = getNumberDensity(grid, pp, i_mixture);
                    dens += i_dens;
                    sum += mixed_component[i_mixture].getCabsMean(grid, pp) * i_dens;
                }
                if(dens != 0)
                    sum /= dens;
            }
        }
    }
    return sum;
}

double CDustMixture::getCscaMean(CGridBasic * grid, const photon_package & pp) const
{
    double sum = 0;
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, pp);
            sum = mixed_component[i_mixture].getCscaMean(grid, pp);
        }
        else
        {
            if(grid->getCscaMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID()) != MAX_DOUBLE)
                sum = grid->getCscaMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID());
            else
            {
                double dens = 0;
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                {
                    double i_dens = getNumberDensity(grid, pp, i_mixture);
                    dens += i_dens;
                    sum += mixed_component[i_mixture].getCscaMean(grid, pp) * i_dens;
                }
                if(dens != 0)
                    sum /= dens;
            }
        }
    }
    return sum;
}

double CDustMixture::getMaxSubTemperature() const
{
    if(mixed_component == 0)
        return 0;
            
    double T_max = 0;     
   
    for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
    {
        double T = mixed_component[i_mixture].getSublimationTemperature();

        T_max = max(T_max, T);
    }
    
    return T_max;
}

bool CDustMixture::adjustTempAndWavelengthBW(CGridBasic * grid, photon_package * pp, bool use_energy_density, CRandomGenerator * rand_gen)
{
    if(mixed_component != 0)
    {
        uint i_mixture = getEmittingMixture(grid, pp, rand_gen);
        return mixed_component[i_mixture].adjustTempAndWavelengthBW(
            grid, pp, i_mixture, use_energy_density, rand_gen);
    }
    return false;
}

void CDustMixture::calcAlignedRadii(CGridBasic * grid, cell_basic * cell)
{
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, *cell);
            mixed_component[i_mixture].calcAlignedRadii(grid, cell, 0);
        }
        else
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                mixed_component[i_mixture].calcAlignedRadii(grid, cell, i_mixture);
    }
}

void CDustMixture::addToWavelengthGrid(double wavelength)
{
    // Add wavelength to list
    wavelength_list.push_back(wavelength);
}

void CDustMixture::addToWavelengthGrid(double lam_min, double lam_max, double nr_of_wavelength, bool add_offset)
{
    dlist tmp_wavelength_list(nr_of_wavelength);

    // Add the amount of added wavelength to the offset to skip them for
    // writeComponent
    if(add_offset && wavelength_list.size() == 0)
        wavelength_offset += nr_of_wavelength;

    // Add new wavelengths first
    CMathFunctions::LogList(lam_min, lam_max, tmp_wavelength_list, 10);

    // Update wavelengths
    for(uint w = 0; w < nr_of_wavelength; w++)
        wavelength_list.push_back(tmp_wavelength_list[w]);
}

void CDustMixture::finalizeWavelengthList()
{
    sort(wavelength_list.begin(), wavelength_list.end());
    wavelength_list.erase(unique(wavelength_list.begin(), wavelength_list.end()), wavelength_list.end());

    nr_of_wavelength = wavelength_list.size();
}

void CDustMixture::setScatteringToRay(bool val)
{
    scattering_to_raytracing = val;
}

bool CDustMixture::getScatteringToRay()
{
    return scattering_to_raytracing;
}

const dlist & CDustMixture::getWavelengthList() const
{
    return wavelength_list;
}

double CDustMixture::getWavelength(uint wID)
{
    return wavelength_list[wID];
}

double CDustMixture::getWavelength(photon_package * pp)
{
    return wavelength_list[pp->getDustWavelengthID()];
}

uint CDustMixture::getNrOfWavelength()
{
    return nr_of_wavelength;
}

uint CDustMixture::getWavelengthID(double wavelength)
{
    dlist::iterator it = find(wavelength_list.begin(), wavelength_list.end(), wavelength);
    if(it != wavelength_list.end())
        return distance(wavelength_list.begin(), it);

    cout << WARNING_LINE << "Wavelength not found!" << endl;
    return 0;
}

double CDustMixture::getSizeMin(uint i_mixture)
{
    return mixed_component[i_mixture].getSizeMin();
}

double CDustMixture::getSizeMin(CGridBasic * grid, const cell_basic & cell) const
{
    double min_a = 1e200;
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, cell);
            min_a = mixed_component[i_mixture].getSizeMin();
        }
        else
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            {
                double tmp_a = mixed_component[i_mixture].getSizeMin();
                if(tmp_a < min_a)
                    min_a = tmp_a;
            }
    }
    return min_a;
}

double CDustMixture::getEffectiveRadius(uint i_mixture, uint a)
{
    return mixed_component[i_mixture].getEffectiveRadius(a);
}

double * CDustMixture::getEffectiveRadii(uint i_mixture)
{
    return mixed_component[i_mixture].getEffectiveRadii();
}

uint CDustMixture::getNrOfCalorimetryTemperatures(uint i_mixture)
{
    return mixed_component[i_mixture].getNrOfCalorimetryTemperatures();
}

// double CDustMixture::getMinDustTemp()
// {
//     double min_temp = 1e200;
//     if(mixed_component != 0)
//         for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
//         {
//             double temp = mixed_component[i_mixture].getMinDustTemp();
//             if(temp < min_temp)
//                 min_temp = temp;
//         }
//     return min_temp;
// }

double CDustMixture::getMaxDustTemp()
{
    double max_temp = 0;
    if(mixed_component != 0)
        for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
        {
            double temp = mixed_component[i_mixture].getMaxDustTemp();
            if(temp > max_temp)
                max_temp = temp;
        }
    return max_temp;
}

double CDustMixture::getMinAlignedRadius()
{
    double min_a_alig = 1e200;
    if(mixed_component != 0)
        for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
        {
            double a_alig = mixed_component[i_mixture].getMinAlignedRadius();
            if(a_alig < min_a_alig)
                min_a_alig = a_alig;
        }
    return min_a_alig;
}

double CDustMixture::getMaxAlignedRadius()
{
    double max_a_alig = 0;
    if(mixed_component != 0)
        for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
        {
            double a_alig = mixed_component[i_mixture].getMaxAlignedRadius();
            if(a_alig > max_a_alig)
                max_a_alig = a_alig;
        }
    return max_a_alig;
}

// Dust number density functions
double CDustMixture::getNumberDensity(CGridBasic * grid, const cell_basic & cell) const
{
    double sum = 0;
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, cell);
            sum = mixed_component[i_mixture].getNumberDensity(grid, cell);
        }
        else
        {
            if(grid->getNumberDensityTab(cell.getUniqueID()) != MAX_DOUBLE)
                sum = grid->getNumberDensityTab(cell.getUniqueID());
            else
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    sum += mixed_component[i_mixture].getNumberDensity(grid, cell, i_mixture);
        }
    }
    return sum;
}

double CDustMixture::getNumberDensity(CGridBasic * grid, const photon_package & pp) const
{
    return getNumberDensity(grid, *pp.getPositionCell());
}

double CDustMixture::getNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_mixture) const
{
    double sum = 0;
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            if(i_mixture == getMixtureID(grid, cell))
                sum = mixed_component[i_mixture].getNumberDensity(grid, cell);
        }
        else
            sum = mixed_component[i_mixture].getNumberDensity(grid, cell, i_mixture);
    }
    return sum;
}

double CDustMixture::getNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_mixture) const
{
    return getNumberDensity(grid, *pp.getPositionCell(), i_mixture);
}

// Dust mass density functions
double CDustMixture::getMassDensity(CGridBasic * grid, const cell_basic & cell) const
{
    double sum = 0;
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, cell);
            sum = mixed_component[i_mixture].getMassDensity(grid, cell);
        }
        else
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                sum += mixed_component[i_mixture].getMassDensity(grid, cell, i_mixture);
    }
    return sum;
}

double CDustMixture::getMassDensity(CGridBasic * grid, const photon_package & pp) const
{
    return getMassDensity(grid, *pp.getPositionCell());
}

double CDustMixture::getMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_mixture) const
{
    double sum = 0;
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            if(i_mixture == getMixtureID(grid, cell))
                sum = mixed_component[i_mixture].getMassDensity(grid, cell);
        }
        else
            sum = mixed_component[i_mixture].getMassDensity(grid, cell, i_mixture);
    }
    return sum;
}

double CDustMixture::getMassDensity(CGridBasic * grid, const photon_package & pp, uint i_mixture) const
{
    return getMassDensity(grid, *pp.getPositionCell(), i_mixture);
}

double CDustMixture::getRelativeDustNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const
{
    double dens = getNumberDensity(grid, cell);
    if(dens != 0)
        return getNumberDensity(grid, cell, i_density) / dens;
    else
        return 0;
}

double CDustMixture::getRelativeDustNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_density) const
{
    return getRelativeDustNumberDensity(grid, *pp.getPositionCell(), i_density);
}

double CDustMixture::getRelativeDustMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const
{
    double dens = getMassDensity(grid, cell);
    if(dens != 0)
        return getMassDensity(grid, cell, i_density) / dens;
    else
        return 0;
}

double CDustMixture::getRelativeDustMassDensity(CGridBasic * grid, const photon_package & pp, uint i_density) const
{
    return getRelativeDustMassDensity(grid, *pp.getPositionCell(), i_density);
}

void CDustMixture::markCells(CGridBasic * grid, parameters & param)
{
    ulong max_cells = grid->getMaxDataCells();
    int sub_status = param.getSubStatus();

    if(sub_status==0)
        return;
        
    if(mixed_component==0)
        return;
    
    cout << CLR_LINE;
    cout << " -> Initiating sublimation markers ... \r" << flush;
    for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
    {
        mixed_component[i_mixture].initMarker(max_cells);
    }
    
    dlist sources_list = param.getPointSources();
    uint nr_stars = param.getNrOfPointSources();

    for(uint i_star = 0; i_star < sources_list.size(); i_star += NR_OF_POINT_SOURCES)
    {
        cout << CLR_LINE;
        cout << " -> Marking sublimation cells for star " << i_star +1 <<" of " << nr_stars << " : 0%  \r" << flush;
        
        uint index = i_star / NR_OF_POINT_SOURCES;
        cell_basic * cell_star = 0;
        uint per_counter=0;
        uint last_percentage=0;

        double r_sub = sources_list[index + 5];
        Vector3D pos_star = Vector3D(sources_list[index + 0], sources_list[index + 1] ,sources_list[index + 2]);

        photon_package pp;
        pp.setPosition(pos_star);
            
        if(grid->findStartingPoint(&pp))
        {
            cell_star = pp.getPositionCell();
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
        {
            bool mark = false;
            const cell_basic * current_cell = grid->getCellFromIndex(i_cell);
            
            
            
            per_counter++;

            // Calculate percentage of total progress per source
            float percentage = 100 * float(per_counter) / float(max_cells);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {   
                #pragma omp critical
                {

                    cout << " -> Marking sublimation cells for star " << i_star +1 
                            << " of " << nr_stars << " : " << 100.0*float(per_counter)/float(max_cells) << "       \r" << flush;
                }
                
                last_percentage = percentage;
            }

            if(cell_star==current_cell)
            {
                mark=true;
            }

            if(sub_status>SUB_CENTER && !mark)
            {
                if(r_sub>0)
                {
                    Vector3D cell_center = grid->getCenter(*current_cell);
                    Vector3D diff=cell_center-pos_star;
                    double distance = diff.length();

                    double vol=grid->getVolume(*current_cell);
                    double rel_dist = cbrt((3.0 * vol) / (PIx4));

                    if(distance< (r_sub+rel_dist))
                        mark=true;
                }
            }

            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            {
                double sub_temp = mixed_component[i_mixture].getSublimationTemperature();
                const double T_dust = grid->getDustTemperature(*current_cell, i_mixture);   

                if(T_dust>=sub_temp || mark)
                {
                    mixed_component[i_mixture].setMarker(i_cell,1);
                }
            }
        }
    }
    
    cout << CLR_LINE;
    cout << "- Marking sublimation cells: done \n" << flush;
}

void CDustMixture::calcEmissivityHz(CGridBasic * grid, const photon_package & pp, StokesVector * dust_emissivity)
{
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, pp);
            mixed_component[i_mixture].calcEmissivityHz(grid, pp, 0, dust_emissivity);
        }
        else
        {
            *dust_emissivity = 0;
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            {
                mixed_component[i_mixture].calcEmissivityHz(grid, pp, i_mixture, dust_emissivity);
            }
        }
    }
}

double CDustMixture::calcEmissivity(CGridBasic * grid, const photon_package & pp)
{
    // Init variables
    double pl_abs = 0;

    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, pp);
            pl_abs = mixed_component[i_mixture].calcEmissivity(grid, pp, 0);
        }
        else
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                pl_abs += mixed_component[i_mixture].calcEmissivity(grid, pp, i_mixture);
    }

    return pl_abs;
}

void CDustMixture::calcEmissivityExt(CGridBasic * grid, const photon_package & pp, Matrix2D * dust_ext_matrix) const
{
    // Init variables
    double Cext = 0, Cpol = 0, Ccirc = 0;

    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, pp);
            double dens_dust = mixed_component[i_mixture].getNumberDensity(grid, pp);
            mixed_component[i_mixture].calcExtCrossSections(grid, pp, 0, &Cext, &Cpol, &Ccirc);
            Cext *= -dens_dust;
            Cpol *= -dens_dust;
            Ccirc *= -dens_dust;
        }
        else
        {
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            {
                double tmpCext, tmpCpol, tmpCcirc;
                double dens_dust = mixed_component[i_mixture].getNumberDensity(grid, pp, i_mixture);
                mixed_component[i_mixture].calcExtCrossSections(
                    grid, pp, i_mixture, &tmpCext, &tmpCpol, &tmpCcirc);
                Cext -= tmpCext * dens_dust;
                Cpol -= tmpCpol * dens_dust;
                Ccirc -= tmpCcirc * dens_dust;
            }
        }
    }

    // Create extinction matrix for the dust grains
    double phi = grid->getPhiMag(pp);
    double sin_2ph = sin(2.0 * phi);
    double cos_2ph = cos(2.0 * phi);

    // Reset matrix
    dust_ext_matrix->resize(4, 4);

    dust_ext_matrix->setValue(0, 0, Cext);
    dust_ext_matrix->setValue(1, 1, Cext);
    dust_ext_matrix->setValue(2, 2, Cext);
    dust_ext_matrix->setValue(3, 3, Cext);

    dust_ext_matrix->setValue(0, 1, Cpol * cos_2ph);
    dust_ext_matrix->setValue(0, 2, Cpol * sin_2ph);

    dust_ext_matrix->setValue(1, 0, Cpol * cos_2ph);
    dust_ext_matrix->setValue(2, 0, Cpol * sin_2ph);

    dust_ext_matrix->setValue(1, 3, Ccirc * sin_2ph);
    dust_ext_matrix->setValue(2, 3, -Ccirc * cos_2ph);

    dust_ext_matrix->setValue(3, 1, -Ccirc * sin_2ph);
    dust_ext_matrix->setValue(3, 2, Ccirc * cos_2ph);
}

void CDustMixture::calcEmissivityEmi(CGridBasic * grid,
                        const photon_package & pp,
                        uint i_offset,
                        uint emission_component,
                        StokesVector * dust_emissivity) const
{
    // Init variables
    double energy = 0;
    double phi = grid->getPhiMag(pp);
    Vector3D en_dir;

    // Check if radiation field is available and scattering should be included
    if(scattering_to_raytracing &&
        (emission_component == DUST_EMI_FULL || emission_component == DUST_EMI_SCAT))
    {
        // Get wavelength of photon package
        uint w = pp.getDustWavelengthID();

        // Get radiation field and calculate angle to the photon package direction
        if(grid->isRadiationFieldAvailable())
            grid->getRadiationFieldInterp(pp, wavelength_list[w], &energy, &en_dir);
        else if(i_offset != MAX_UINT)
        {
            // Use the rad field as stokes vector if all necessary things were already calculated
            *dust_emissivity += grid->getStokesFromRadiationField(pp, i_offset);
        }
        else
            grid->getRadiationField(pp, w, &energy, &en_dir);
    }

    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, pp);
            *dust_emissivity += mixed_component[i_mixture].calcEmissivityEmi(
                grid, pp, 0, emission_component, phi, energy, en_dir);
        }
        else
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                *dust_emissivity += mixed_component[i_mixture].calcEmissivityEmi(
                    grid, pp, i_mixture, emission_component, phi, energy, en_dir);
    }
}

StokesVector CDustMixture::getRadFieldScatteredFraction(CGridBasic * grid,
                                            const photon_package & pp,
                                            const Vector3D & en_dir,
                                            double energy) const
{
    StokesVector tmp_stokes;

    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, pp);
            tmp_stokes =
                mixed_component[i_mixture].getRadFieldScatteredFraction(grid, pp, 0, en_dir, energy);
        }
        else
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                tmp_stokes += mixed_component[i_mixture].getRadFieldScatteredFraction(
                    grid, pp, i_mixture, en_dir, energy);
    }
    return tmp_stokes;
}

uint CDustMixture::getNrOfMixtures() const
{
    return dust_choices.size();
}

void CDustMixture::convertTempInQB(CGridBasic * grid, cell_basic * cell, double min_gas_density, bool use_gas_temp)
{
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, *cell);
            mixed_component[i_mixture].convertTempInQB(grid, cell, 0, min_gas_density, use_gas_temp);
        }
        else
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                mixed_component[i_mixture].convertTempInQB(
                    grid, cell, i_mixture, min_gas_density, use_gas_temp);
    }
}

uint CDustMixture::getScatteringMixture(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen) const
{
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, *pp);
            return i_mixture;
        }
        else
        {
            double rnd = rand_gen->getRND();
            double dens = 0;
            double pb[getNrOfMixtures()];

            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            {
                double i_dens = getNumberDensity(grid, *pp, i_mixture);
                pb[i_mixture] = i_dens * mixed_component[i_mixture].getCscaMean(grid, *pp);
                if(i_mixture > 0)
                    pb[i_mixture] += pb[i_mixture-1];
                dens += i_dens;
            }

            double cscamean_tmp = getCscaMean(grid, *pp);

            if(dens!=0 && cscamean_tmp != 0)
            {
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                {
                    pb[i_mixture] /= dens * cscamean_tmp;
                    if(pb[i_mixture] > (1-rnd))
                        return i_mixture;
                }
            }
        }
    }
    return 0;
}

uint CDustMixture::getEmittingMixture(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen) const
{
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, *pp);
            return i_mixture;
        }
        else
        {
            double rnd = rand_gen->getRND();
            double dens = 0;
            double pb[getNrOfMixtures()];

            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            {
                double i_dens = getNumberDensity(grid, *pp, i_mixture);
                pb[i_mixture] = i_dens * mixed_component[i_mixture].getCabsMean(grid, *pp);
                if(i_mixture > 0)
                    pb[i_mixture] += pb[i_mixture-1];
                dens += i_dens;
            }

            double cabsmean_tmp = getCabsMean(grid, *pp);

            if(dens!=0 && cabsmean_tmp != 0)
            {
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                {
                    pb[i_mixture] /= dens * cabsmean_tmp;
                    if(pb[i_mixture] > (1-rnd))
                        return i_mixture;
                }
            }
        }
    }
    return 0;
}

void CDustMixture::scatter(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen, bool adjust_stokes)
{
    if(mixed_component != 0)
    {
        uint i_mixture = getScatteringMixture(grid, pp, rand_gen);
        mixed_component[i_mixture].scatter(grid, pp, rand_gen);

        // Reduce the Stokes vector by the mean albedo of the particles
        if(adjust_stokes)
        {
            if(getCextMean(grid, *pp) > 0.0)
                *pp->getStokesVector() *= getCscaMean(grid, *pp) / getCextMean(grid, *pp);
            else
            {
                cout << INFO_LINE << "Mean cross section for extinction is zero or negative!" << endl;
                pp->getStokesVector()->clear();
            }
        }
    }
}

void CDustMixture::setGridRequirements(CGridBasic * grid, parameters & param)
{
    uint nr_of_mixtures = getNrOfMixtures();
    uint * nr_dust_temp_sizes = new uint[nr_of_mixtures];
    uint * nr_stochastic_temps = new uint[nr_of_mixtures];
    uint * nr_stochastic_sizes = new uint[nr_of_mixtures];
    if(mixed_component != 0)
    {
        for(uint i_mixture = 0; i_mixture < nr_of_mixtures; i_mixture++)
        {
            // Maximum amount of dust grain sizes for which the grid has to contain a
            // temperature
            nr_dust_temp_sizes[i_mixture] = getNrOfDustSpecies(i_mixture);
            // Temperatures for which the griad has to contain propabilities
            nr_stochastic_temps[i_mixture] = getNrOfCalorimetryTemperatures(i_mixture);
            // Maximum amount of dust grain sizes affected by stochastic heating
            nr_stochastic_sizes[i_mixture] = getNrOfStochasticSizes(i_mixture);
        }
    }
    grid->setDustInformation(
        nr_of_mixtures, nr_dust_temp_sizes, nr_stochastic_sizes, nr_stochastic_temps);
}

void CDustMixture::setIDs(CDustComponent & component,
            uint i_comp,
            uint nr_of_components,
            uint i_mixture,
            uint nr_of_mixtures)
{
    component.setIDs(i_comp, nr_of_components, i_mixture, nr_of_mixtures);
}

double *** CDustMixture::getSizeFractions()
{
    // Get number of dust size bins
    uint nr_of_dust_species = single_component[0].getNrOfDustSpecies();

    // Init pointer array
    double *** size_fraction = new double **[nr_of_components];
    for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
    {
        size_fraction[i_comp] = new double *[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            size_fraction[i_comp][a] = new double[2];
    }

    // Put the mass weights into relation to each other
    double ref_mass_weight = single_component[0].getMassWeight() / single_component[0].getFraction();
    for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
    {
        double mass_weight_ratio = ref_mass_weight * single_component[i_comp].getFraction() /
                                    single_component[i_comp].getMassWeight();

        // First column is the updated size distribution
        for(uint a = 0; a < nr_of_dust_species; a++)
            if(single_component[i_comp].sizeIndexUsed(a))
            {
                size_fraction[i_comp][a][0] =
                    mass_weight_ratio * single_component[i_comp].getGrainSizeDistribution(a);
            }
    }

    for(uint a = 0; a < nr_of_dust_species; a++)
    {
        // Calulcate the sum for each size bin
        double sum = 0;
        for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
            if(single_component[i_comp].sizeIndexUsed(a))
                sum += size_fraction[i_comp][a][0];

        // Second column is the mixing ration between the dust components
        for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
            if(single_component[i_comp].sizeIndexUsed(a))
                size_fraction[i_comp][a][1] = size_fraction[i_comp][a][0] / sum;
    }
    return size_fraction;
}

void CDustMixture::getEscapePhoton(CGridBasic * grid,
                        photon_package * pp,
                        Vector3D obs_ex,
                        Vector3D dir_obs,
                        photon_package * pp_escape,
                        CRandomGenerator * rand_gen) const
{
    if(mixed_component != 0)
    {
        uint i_mixture = getScatteringMixture(grid, pp, rand_gen);
        uint a = mixed_component[i_mixture].getInteractingDust(grid, pp, rand_gen, CROSS_SCA);
        mixed_component[i_mixture].getEscapePhoton(grid, pp, a, obs_ex, dir_obs, pp_escape);

        // Reduce the Stokes vector by the mean albedo of the particles
        if(getCextMean(grid, *pp) > 0.0)
            *pp_escape->getStokesVector() *= getCscaMean(grid, *pp) / getCextMean(grid, *pp);
        else
        {
            cout << INFO_LINE << "Mean cross section for extinction is zero or negative!" << endl;
            pp_escape->getStokesVector()->clear();
        }
    }
}

double CDustMixture::getCellEmission(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen) const
{
    if(mixed_component != 0)
    {
        uint i_mixture = getEmittingMixture(grid, pp, rand_gen);
        if(grid->useDustChoice())
            return mixed_component[i_mixture].getCellEmission(grid, *pp, 0);
        else
            return mixed_component[i_mixture].getCellEmission(grid, *pp, i_mixture);
    }
    return 0;
}

double CDustMixture::getTotalCellEmission(CGridBasic * grid, const photon_package & pp) const
{
    double sum = 0;
    if(mixed_component != 0)
    {
        if(grid->useDustChoice())
        {
            uint i_mixture = getMixtureID(grid, pp);
            sum = mixed_component[i_mixture].getCellEmission(grid, pp, 0);
        }
        else
        {
            if(grid->getTotalCellEmissionTab(pp.getPositionCell()->getUniqueID()) != MAX_DOUBLE)
                sum = grid->getTotalCellEmissionTab(pp.getPositionCell()->getUniqueID());
            else
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    sum += mixed_component[i_mixture].getCellEmission(grid, pp, i_mixture);
        }
    }
    return sum;
}

double CDustMixture::getPlanck(uint w, double temp)
{
    return mixed_component[0].getPlanck(w, temp);
}
