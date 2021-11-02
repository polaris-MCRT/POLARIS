#include "RadiativeTransfer.h"
#include "CommandParser.h"
#include "Detector.h"
#include "Stokes.h"
#include <complex>
#include <limits>

#define XAxis 10
#define YAxis 20
#define ZAxis 30

#define FULL 1
#define PARA 2
#define ORTH 3

bool CRadiativeTransfer::initiateDustRaytrace(parameters & param)
{
    // Check if grid was loaded
    if(grid == 0)
    {
        cout << "\nERROR: No Grid loaded!" << endl;
        return false;
    }

    // Check if dust was loaded
    if(dust == 0)
    {
        cout << "\nERROR: No dust model!" << endl;
        return false;
    }

    // Get list of dust raytracing detector parameters
    dlist dust_ray_detectors = param.getDustRayDetectors();

    // No raytracing without detectors
    if(dust_ray_detectors.size() == 0)
    {
        cout << "\nERROR: No sequence defined!" << endl;
        return false;
    }

    // Get number of dust raytracing detectors
    nr_ray_detectors = uint(dust_ray_detectors.size()) / NR_OF_RAY_DET;

    // Get start and stop id of detectors
    start = param.getStart();
    stop = param.getStop();

    // Get maximum length of the simulation model
    double max_length = grid->getMaxLength();

    // Init array of tracer base class pointer
    tracer = new CRaytracingBasic *[nr_ray_detectors];

    // Init 2D list to get the offset index to obtain data from grid
    detector_wl_index.push_back(0);

    // Check other parameters for raytracing
    for(uint i_det = 0; i_det < nr_ray_detectors; i_det++)
    {
        // Calculate the starting position of each sequence
        uint pos = i_det * NR_OF_RAY_DET;

        // Push back the number of wavelengths of the current detector (plus the one before)
        detector_wl_index.push_back(detector_wl_index.back() + uint(dust_ray_detectors[pos + 2]));

        uint nr_source = uint(dust_ray_detectors[pos + 3]);

        // Get the ID of the chosen source
        uint detector_id = uint(dust_ray_detectors[pos + NR_OF_RAY_DET - 3]);

        // Check for wrong choice of emission source
        if(nr_source > sources_ray.size())
        {
            cout << "\nERROR: ID of source (" << nr_source << ") larger than max. amount ("
                 << sources_ray.size() << ") of defined sources!" << endl;
            return false;
        }

        // Polar raytracing background grid can only be used with cylindrical or spherical grids
        if(detector_id == DET_POLAR && grid->getDataID() != GRID_ID_SPH && grid->getDataID() != GRID_ID_CYL)
        {
            cout << "\nERROR: Polar RT grid can only be used with spherical and "
                    "cylindrical grids!"
                 << endl;
            return false;
        }

        // Enable Stokes rad field for all but HealPix detector
        stokes_dust_rad_field = true;

        // Create detector for current simulation
        switch(detector_id)
        {
            case DET_PLANE:
                tracer[i_det] = new CRaytracingCartesian(grid);
                if(!tracer[i_det]->setDustDetector(pos,
                                                   dust_ray_detectors,
                                                   max_length,
                                                   pathOutput,
                                                   param.getMaxSubpixelLvl(),
                                                   param.getAxis1(),
                                                   param.getAxis2(),
                                                   param.getAlignmentMechanism()))
                    return false;
                break;

            case DET_SPHER:
                tracer[i_det] = new CRaytracingHealPix(grid);
                if(!tracer[i_det]->setDustDetector(pos,
                                                   dust_ray_detectors,
                                                   max_length,
                                                   pathOutput,
                                                   param.getAlignmentMechanism(),
                                                   param.getHealpixOrientation()))
                    return false;
                // Disable Stokes rad field for HealPix detector
                stokes_dust_rad_field = false;
                break;

            case DET_POLAR:
                tracer[i_det] = new CRaytracingPolar(grid);
                if(!tracer[i_det]->setDustDetector(pos,
                                                   dust_ray_detectors,
                                                   max_length,
                                                   pathOutput,
                                                   param.getMaxSubpixelLvl(),
                                                   param.getAxis1(),
                                                   param.getAxis2(),
                                                   param.getAlignmentMechanism()))
                    return false;
                break;

            case DET_SLICE:
                tracer[i_det] = new CRaytracingSlice(grid);
                if(!tracer[i_det]->setDustDetector(pos,
                                                   dust_ray_detectors,
                                                   max_length,
                                                   pathOutput,
                                                   param.getAxis1(),
                                                   param.getAxis2(),
                                                   param.getAlignmentMechanism()))
                    return false;
                break;
        }
    }

    // Initiate RKF coefficients for numerical solving the radiative transfer equation
    initiateRungeKuttaFehlberg();

    return true;
}

bool CRadiativeTransfer::initiateSyncRaytrace(parameters & param)
{
    if(grid == 0)
    {
        cout << "\nERROR: No Grid loaded!" << endl;
        return false;
    }

    if(dust == 0)
    {
        cout << "\nERROR: No dust model!" << endl;
        return false;
    }

    dlist sync_ray_detectors = param.getSyncRayDetectors();

    if(sync_ray_detectors.size() == 0)
    {
        cout << "\nERROR: No sequence defined!" << endl;
        return false;
    }

    nr_ray_detectors = uint(sync_ray_detectors.size()) / NR_OF_RAY_DET;

    // Get start and stop id of detectors
    start = param.getStart();
    stop = param.getStop();

    // Get maximum length of the simulation model
    double max_length = grid->getMaxLength();

    // Init array of tracer base class pointer
    tracer = new CRaytracingBasic *[nr_ray_detectors];

    for(uint i_det = 0; i_det < nr_ray_detectors; i_det++)
    {
        // Calculate the starting position of each sequence
        uint pos = i_det * NR_OF_RAY_DET;

        uint nr_source = uint(sync_ray_detectors[pos + 3]);

        // Get the ID of the chosen source
        uint detector_id = uint(sync_ray_detectors[pos + NR_OF_RAY_DET - 3]);

        if(nr_source > sources_ray.size())
        {
            cout << "\nERROR: ID of source (" << nr_source << ") larger than max. amount ("
                 << sources_ray.size() << ") of defined sources!" << endl;
            return false;
        }

        if(detector_id == DET_POLAR && grid->getDataID() != GRID_ID_SPH && grid->getDataID() != GRID_ID_CYL)
        {
            cout << "\nERROR: Polar RT grid can only be used with spherical and "
                    "cylindrical grids!"
                 << endl;
            return false;
        }

        // Create detector for current simulation
        switch(detector_id)
        {
            case DET_PLANE:
                tracer[i_det] = new CRaytracingCartesian(grid);
                if(!tracer[i_det]->setSyncDetector(pos,
                                                   sync_ray_detectors,
                                                   max_length,
                                                   pathOutput,
                                                   param.getMaxSubpixelLvl(),
                                                   param.getAxis1(),
                                                   param.getAxis2()))
                    return false;
                break;

            case DET_SPHER:
                tracer[i_det] = new CRaytracingHealPix(grid);
                if(!tracer[i_det]->setSyncDetector(
                       pos, sync_ray_detectors, max_length, pathOutput, param.getHealpixOrientation()))
                    return false;
                break;

            case DET_POLAR:
                tracer[i_det] = new CRaytracingPolar(grid);
                if(!tracer[i_det]->setSyncDetector(pos,
                                                   sync_ray_detectors,
                                                   max_length,
                                                   pathOutput,
                                                   param.getMaxSubpixelLvl(),
                                                   param.getAxis1(),
                                                   param.getAxis2()))
                    return false;
                break;

            case DET_SLICE:
                tracer[i_det] = new CRaytracingSlice(grid);
                if(!tracer[i_det]->setSyncDetector(
                       pos, sync_ray_detectors, max_length, pathOutput, param.getAxis1(), param.getAxis2()))
                    return false;
                break;
        }
    }

    synchrotron = new CSynchrotron();

    initiateRungeKuttaFehlberg();

    return true;
}

bool CRadiativeTransfer::initiateLineRaytrace(parameters & param)
{
    if(grid == 0)
    {
        cout << "\nERROR: No Grid loaded!" << endl;
        return false;
    }

    if(dust == 0)
    {
        cout << "\nERROR: No dust model defined!" << endl;
        return false;
    }

    if(gas == 0)
    {
        cout << "\nERROR: No gas model loaded!" << endl;
        return false;
    }

    uint nrOfGasSpecies = param.getNrOfGasSpecies();

    // Get start and stop id of detectors
    start = param.getStart();
    stop = param.getStop();

    // Get maximum length of the simulation model
    double max_length = grid->getMaxLength();

    // Init array of tracer base class pointer
    nr_ray_detectors = param.getTotalNrOfLineTransitions();
    tracer = new CRaytracingBasic *[nr_ray_detectors];

    if(nrOfGasSpecies == 0)
    {
        cout << "\nERROR: No gas species transition defined!" << endl;
        return false;
    }

    uint i_det = 0;
    for(uint i_species = 0; i_species < nrOfGasSpecies; i_species++)
    {
        dlist line_ray_detectors = param.getLineRayDetector(i_species);

        for(uint i_line = 0; i_line < param.getNrOfGasSpeciesTransitions(i_species); i_line++)
        {
            // Calculate the starting position of each detector for a gas species
            uint pos = NR_OF_LINE_DET * i_line;

            uint nr_source = uint(line_ray_detectors[pos + 1]);

            // Get the ID of the chosen detector
            uint detector_id = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 4]);

            if(nr_source > sources_ray.size())
            {
                cout << "\nERROR: ID of source (" << nr_source << ") larger than max. amount ("
                     << sources_ray.size() << ") of defined sources!" << endl;
                return false;
            }

            if(detector_id == DET_POLAR && grid->getDataID() != GRID_ID_SPH &&
               grid->getDataID() != GRID_ID_CYL)
            {
                cout << "\nERROR: Polar RT grid can only be used with spherical and "
                        "cylindrical grids!"
                     << endl;
                return false;
            }

            // Create detector for current simulation
            switch(detector_id)
            {
                case DET_PLANE:
                    tracer[i_det] = new CRaytracingCartesian(grid);
                    if(!tracer[i_det]->setLineDetector(pos,
                                                       line_ray_detectors,
                                                       pathOutput,
                                                       max_length,
                                                       param.getVelMaps(),
                                                       param.getMaxSubpixelLvl(),
                                                       param.getAxis1(),
                                                       param.getAxis2()))
                        return false;
                    break;

                case DET_SPHER:
                    tracer[i_det] = new CRaytracingHealPix(grid);
                    if(!tracer[i_det]->setLineDetector(pos,
                                                       line_ray_detectors,
                                                       pathOutput,
                                                       max_length,
                                                       param.getVelMaps(),
                                                       param.getHealpixOrientation()))
                        return false;
                    break;

                case DET_POLAR:
                    tracer[i_det] = new CRaytracingPolar(grid);
                    if(!tracer[i_det]->setLineDetector(pos,
                                                       line_ray_detectors,
                                                       pathOutput,
                                                       max_length,
                                                       param.getVelMaps(),
                                                       param.getMaxSubpixelLvl(),
                                                       param.getAxis1(),
                                                       param.getAxis2()))
                        return false;
                    break;

                case DET_SLICE:
                    tracer[i_det] = new CRaytracingSlice(grid);
                    if(!tracer[i_det]->setLineDetector(pos,
                                                       line_ray_detectors,
                                                       pathOutput,
                                                       max_length,
                                                       param.getVelMaps(),
                                                       param.getAxis1(),
                                                       param.getAxis2()))
                        return false;
                    break;
            }

            // Increment index for line RT since its distributed over species and transitions
            i_det++;
        }
    }

    initiateRungeKuttaFehlberg();

    return true;
}

bool CRadiativeTransfer::initiateOPIATE(parameters & param)
{
    if(grid == 0)
    {
        cout << "\nERROR: No Grid loaded!" << endl;
        return false;
    }

    if(dust == 0)
    {
        cout << "\nHINT: No dust model defined!" << endl;
        return false;
    }

    // Get start and stop id of detectors
    start = param.getStart();
    stop = param.getStop();

    initiateRungeKuttaFehlberg();

    return true;
}

void CRadiativeTransfer::initiateDustMC(parameters & param)
{
    nr_mc_detectors = param.getNrOfDustMCDetectors();
    b_forced = param.getEnfScattering();
    peel_off = param.getPeelOff();
    mrw_step = param.getMRW();
}

void CRadiativeTransfer::initiateRadFieldMC(parameters & param)
{
    b_forced = param.getEnfScattering();
    mrw_step = param.getMRW();
    adjTgas = param.getAdjTgas();
}

bool CRadiativeTransfer::doMRWStepBWWithoutHeating(photon_package * pp)
{
    return false;
}

bool CRadiativeTransfer::doMRWStepBW(photon_package * pp)
{
    if(!mrw_step)
        return false;

    return true;
}

bool CRadiativeTransfer::calcMonteCarloRadiationField(uint command,
                                                      bool use_energy_density,
                                                      bool disable_reemission)
{
    // Init variables
    ullong nr_of_photons = 0;
    uint nr_used_wavelengths = 1;
    ullong kill_counter = 0;
    uint mrw_counter = 0;
    uint max_source = uint(sources_mc.size());

    // A loop for each source
    for(uint s = 0; s < max_source; s++)
    {
        // Init variables
        ullong per_counter = 0;
        float last_percentage = 0;

        // Init source from sources list
        CSourceBasic * tm_source = sources_mc[s];

        // Init luminosity and probability for a given wavelength to be chosen
        if(!tm_source->initSource(s, max_source, use_energy_density))
            continue;

        // Number of photons per wavelength
        nr_of_photons = tm_source->getNrOfPhotons();

        // Number of wavelength or no loop over wavelength (pick from planck distribution)
        if(use_energy_density)
            nr_used_wavelengths = dust->getNrOfWavelength();

        // Init progress visualization
        per_counter = 0;
        cout << CLR_LINE;
        switch(command)
        {
            case CMD_TEMP_RAT:
                cout << "-> MC temp. and RAT distribution: 0 [%], max. temp. " << dust->getMaxDustTemp()
                     << " [K]      \r" << flush;
                break;

            case CMD_TEMP:
                cout << "-> MC temp. distribution: 0 [%], max. temp. " << dust->getMaxDustTemp()
                     << " [K]      \r" << flush;
                break;

            default:
                cout << "-> MC radiation field: 0 [%]      \r" << flush;
                break;
        }
            // A loop for each wavelength
#pragma omp parallel for schedule(dynamic) collapse(2)
        for(int wID = 0; wID < int(nr_used_wavelengths); wID++)
        {
            // A loop for each photon
            for(llong r = 0; r < llong(nr_of_photons); r++)
            {
                // Init variables
                double end_tau, Cext, Csca;
                Vector3D old_pos;

                // Increase counter used to show progress
                per_counter++;

                // Calculate percentage of total progress per source
                float percentage = 100 * float(per_counter) / float(nr_of_photons * nr_used_wavelengths);

                // Show only new percentage number if it changed
                if((percentage - last_percentage) > PERCENTAGE_STEP)
                {
#pragma omp critical
                    {
                        switch(command)
                        {
                            case CMD_TEMP_RAT:
                                cout << "-> MC temp. and RAT distribution: " << percentage
                                     << " [%], max. temp. " << dust->getMaxDustTemp() << " [K]      \r"
                                     << flush;
                                break;

                            case CMD_TEMP:
                                cout << "-> MC temp. distribution: " << percentage << " [%], max. temp. "
                                     << dust->getMaxDustTemp() << " [K]      \r" << flush;
                                break;

                            default:
                                cout << "-> MC radiation field: " << percentage << " [%]      \r" << flush;
                                break;
                        }
                        last_percentage = percentage;
                    }
                }

                // Init photon package
                photon_package * pp = new photon_package();

                // Launch a new photon package from the source
                if(use_energy_density)
                {
                    pp->setWavelengthID(wID);
                    tm_source->createNextRay(pp, ullong(nr_of_photons * wID + r));
                }
                else
                    tm_source->createNextRay(pp, ullong(time(NULL)*r));

                if(pp->getStokesVector().I() < 1e-200)
                {
                    kill_counter++;
                    delete pp;
                    continue;
                }

                if(!grid->positionPhotonInGrid(pp))
                    if(!grid->findStartingPoint(pp))
                    {
                        kill_counter++;
                        delete pp;
                        continue;
                    }

                // Get tau for first interaction
                end_tau = -log(1.0 - pp->getRND());

                // Save the old position to use it again
                old_pos = pp->getPosition();

                // Init current number of interactions
                ullong interactions = 0;

                // Init variables
                double tmp_tau, len, dens;

                // Transfer photon through grid
                while(grid->next(pp))
                {
                    // If max interactions is reached, end photon transfer
                    if(interactions >= MAX_INTERACTION)
                    {
                        kill_counter++;
                        break;
                    }

                    // Get necessary quantities
                    dens = dust->getNumberDensity(grid, pp);

                    // skip cells without density
                    if(dens == 0)
                    {
                        old_pos = pp->getPosition();
                        continue;
                    }

                    // Get distance to next cell
                    len = pp->getTmpPathLength();

                    // Calculate the dust cross sections (for random alignment)
                    Cext = dust->getCextMean(grid, pp);

                    // Calculate optical depth that is obtained on the path_length
                    tmp_tau = Cext * len * dens;

                    // optical depth for interaction is reached or not
                    if(tmp_tau > end_tau)
                    {
                        // Increase interaction counter
                        interactions++;

                        // Reduce the photon position to match the exact
                        // interaction position
                        pp->adjustPosition(old_pos, len * end_tau / tmp_tau);

                        // Update data in grid like spectral length or radiation field
                        updateRadiationField(pp);

                        if(!doMRWStepBW(pp))
                        {
                            // Calculate the dust scattering cross section (for random
                            // alignment)
                            Csca = dust->getCscaMean(grid, pp);

                            // Calculate albedo and check if absorption or
                            // scattering occurs
                            double albedo = Csca / Cext;

                            if(pp->getRND() < albedo)
                            {
                                // Perform simple photon scattering without
                                // changing the Stokes vectors
                                dust->scatter(grid, pp);
                            }
                            else
                            {
                                // Calculate the temperature of the absorbing cell
                                // and change the wavelength of the photon
                                if(!disable_reemission &&
                                   dust->adjustTempAndWavelengthBW(grid, pp, use_energy_density))
                                {
                                    // Send this photon into a new random direction
                                    pp->calcRandomDirection();
                                }
                                else
                                    break;

                                // ---> Only an idea (not really sufficient) <---
                                // Check if new photon can escape even the smallest cell
                                // if(dust->getCextMean(grid, pp) * grid->getMinLength() * dens > 100)
                                // {
                                //     kill_counter++;
                                //     break;
                                // }
                            }
                        }
                        // Calculate new optical depth for next interaction
                        end_tau = -log(1.0 - pp->getRND());
                    }
                    else
                    {
                        // Update data in grid like spectral length or radiation field
                        updateRadiationField(pp);

                        // Remove the traveled distance from optical depth
                        end_tau -= tmp_tau;
                    }
                    // Save photon position to adjust it if necessary
                    old_pos = pp->getPosition();
                }
                // Delete the pp pointer
                delete pp;
            }
        }
    }

    // Format prints
    cout << CLR_LINE;
    cout << SEP_LINE;

    if(mrw_step)
        cout << "MRW steps: " << mrw_counter << "   " << endl;

    // Show amount of killed photons
    if(kill_counter > 0)
        cout << "- Photons killed                    : " << kill_counter << endl;

    switch(command)
    {
        case CMD_TEMP_RAT:
            cout << "- MC calc. of temperatures and RATs : done                          "
                    "      "
                 << endl;
            break;

        case CMD_TEMP:
            cout << "- MC calculation of temperatures    : done                          "
                    "      "
                 << endl;
            break;

        default:
            cout << "- MC calculation of radiation field : done                          "
                    "      "
                 << endl;
            break;
    }
    return true;
}


// Temperature Calculation with polychromatic photon packages
bool CRadiativeTransfer::calcMonteCarloRadiationFieldPoly(uint command,
                                                      bool use_energy_density,
                                                      bool disable_reemission)
{
    // Init variables
    ullong nr_of_photons = 0;
    ullong kill_counter = 0;
    uint mrw_counter = 0;
    uint max_source = uint(sources_mc.size());

    // A loop for each source
    for(uint s = 0; s < max_source; s++)
    {
        // Init variables
        ullong per_counter = 0;
        float last_percentage = 0;

        // Init source from sources list
        CSourceBasic * tm_source = sources_mc[s];

        // Init luminosity and probability for a given wavelength to be chosen
        if(!tm_source->initSource(s, max_source, use_energy_density))
            continue;

        // Number of photons per wavelength
        nr_of_photons = tm_source->getNrOfPhotons();

        // Number of wavelength or no loop over wavelength (pick from planck distribution)
        uint nr_of_wavelengths = dust->getNrOfWavelength();

        // Init progress visualization
        per_counter = 0;
        cout << CLR_LINE;
        switch(command)
        {
            case CMD_TEMP_RAT:
                cout << "-> MC temp. and RAT distribution: 0 [%], max. temp. " << dust->getMaxDustTemp()
                     << " [K]      \r" << flush;
                break;

            case CMD_TEMP:
                cout << "-> MC temp. distribution: 0 [%], max. temp. " << dust->getMaxDustTemp()
                     << " [K]      \r" << flush;
                break;
                
            case CMD_TEMP_POLY:
                cout << "-> MC polychromatic temp. distribution: 0 [%], max. temp. " << dust->getMaxDustTemp()
                     << " [K]      \r" << flush;
                break;

            default:
                cout << "-> MC radiation field: 0 [%]      \r" << flush;
                break;
        }

#pragma omp parallel for schedule(dynamic)
        // A loop for each photon
        for(llong r = 0; r < llong(nr_of_photons); r++)
        {
            // Init variables
            double end_tau, Cext, Csca;
            Vector3D old_pos;
            

            // Increase counter used to show progress
            per_counter++;

            // Calculate percentage of total progress per source
            float percentage = 100 * float(per_counter) / float(nr_of_photons);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
#pragma omp critical
                {
                    switch(command)
                    {
                        case CMD_TEMP_RAT:
                            cout << "-> MC temp. and RAT distribution: " << percentage
                                    << " [%], max. temp. " << dust->getMaxDustTemp() << " [K]      \r"
                                    << flush;
                            break;

                        case CMD_TEMP:
                            cout << "-> MC temp. distribution: " << percentage << " [%], max. temp. "
                                    << dust->getMaxDustTemp() << " [K]      \r" << flush;
                            break;
                            
                        case CMD_TEMP_POLY:
                            cout << "-> MC polychromatic temp. distribution: " << percentage
                                    << " [%], max. temp. " << dust->getMaxDustTemp() << " [K]      \r"
                                    << flush;
                            break;

                        default:
                            cout << "-> MC radiation field: " << percentage << " [%]      \r" << flush;
                            break;
                    }
                    last_percentage = percentage;
                }
            }

            // Init photon package
            photon_polychrom * pp = new photon_polychrom();

            // Launch a new photon package from the source
            tm_source->createNextRay(pp, ullong(time(NULL)*r));
            tm_source->createSpectrum(pp);
            
            // Set spectrum of photon
            //double T = tm_source->getTemperature();
            //double R = tm_source->getRadius();
            //double L = tm_source->getLuminosity();
            
            //for(uint w = 0; w < nr_of_wavelengths-1; w++)
            //{
            //    double pl = CMathFunctions::planck(dust->getWavelength(w), T);
            //    double sp_energy = 4.0 * PI * PI * (R * R_sun) * (R * R_sun) * pl;
            //    
            //    pp->addToSpectrum(sp_energy/L*(dust->getWavelength(w+1)-dust->getWavelength(w)));
            //}

            if(pp->getStokesVector().I() < 1e-200)
            {
                kill_counter++;
                delete pp;
                continue;
            }

            if(!grid->positionPhotonInGrid(pp))
                if(!grid->findStartingPoint(pp))
                {
                    kill_counter++;
                    delete pp;
                    continue;
                }

            // Get tau for first interaction
            end_tau = -log(1.0 - pp->getRND());

            // Save the old position to use it again
            old_pos = pp->getPosition();

            // Init current number of interactions
            ullong interactions = 0;

            // Init variables
            double tmp_tau, len, dens;

            // Transfer photon through grid
            while(grid->next(pp))
            {
                // If max interactions is reached, end photon transfer
                if(interactions >= MAX_INTERACTION)
                {
                    kill_counter++;
                    break;
                }

                // Get necessary quantities
                dens = dust->getNumberDensity(grid, pp);

                // skip cells without density
                if(dens == 0)
                {
                    old_pos = pp->getPosition();
                    continue;
                }

                // Get distance to next cell
                len = pp->getTmpPathLength();

                // Calculate the dust cross sections (for random alignment)
                Cext = dust->getCextMean(grid, pp);

                // Calculate optical depth that is obtained on the path_length
                tmp_tau = Cext * len * dens;
                
                // optical depth for interaction is reached or not
                if(tmp_tau > end_tau)
                {
                    // Increase interaction counter
                    interactions++;

                    // Reduce the photon position to match the exact
                    // interaction position
                    pp->adjustPosition(old_pos, len * end_tau / tmp_tau);

                    // Update data in grid like spectral length or radiation field for every wavelength
                    double tmp_w = pp->getWavelengthID();
                        
                    for(uint w = 0; w < nr_of_wavelengths; w++)
                    {
                        pp->setWavelengthID(w);
                        updateRadiationField(pp);
                    }
                    
                    pp->setWavelengthID(tmp_w);

                    // Calculate the dust scattering cross section (for random
                    // alignment)
                    Csca = dust->getCscaMean(grid, pp);

                    // Calculate albedo and check if absorption or
                    // scattering occurs
                    double albedo = Csca / Cext;
                    
                    // Weight spectrum of polychromatic photon packages
//                     double tau_weight = 0;
//                     double tau_ref = Cext * len * dens;
//                     tmp_w = pp->getWavelengthID();
//                     
//                     for(uint i = 0; i < nr_of_wavelengths; i++)
//                     {
//                         pp->setWavelengthID(i);
//                         tau_weight = exp(tau_ref*(1-dust->getCextMean(grid, pp)/Cext))*(dust->getCextMean(grid, pp)/Cext);
//                         pp->weightSingleIntensity(i,tau_weight);
//                     }
//                     
//                     pp->setWavelengthID(tmp_w);

                    if(pp->getRND() < albedo)
                    {
                        // Perform simple photon scattering without
                        // changing the Stokes vectors
                        dust->scatter(grid, pp);
                    }
                    else
                    {
                        // Calculate the temperature of the absorbing cell
                        // and change the wavelength of the photon
                        
                        if(!disable_reemission &&
                            dust->adjustTempAndWavelengthBWPoly(grid, pp, use_energy_density))
                            // Send this photon into a new random direction
                            pp->calcRandomDirection();
                        else
                            break;
                    }
                    // Calculate new optical depth for next interaction
                    end_tau = -log(1.0 - pp->getRND());
                }
                else
                {
                    // Update data in grid like spectral length or radiation field for all wavelength
                    double tmp_w = pp->getWavelengthID();
                        
                    for(uint w = 0; w < nr_of_wavelengths; w++)
                    {
                        pp->setWavelengthID(w);
                        updateRadiationField(pp);
                    }
                    
                    pp->setWavelengthID(tmp_w);

                    // Remove the traveled distance from optical depth
                    end_tau -= tmp_tau;
                }
                // Save photon position to adjust it if necessary
                old_pos = pp->getPosition();
            }
            // Delete the pp pointer
            delete pp;
        }
    }

    // Format prints
    cout << CLR_LINE;
    cout << SEP_LINE;

    if(mrw_step)
        cout << "MRW steps: " << mrw_counter << "   " << endl;

    // Show amount of killed photons
    if(kill_counter > 0)
        cout << "- Photons killed                    : " << kill_counter << endl;

    switch(command)
    {
        case CMD_TEMP_RAT:
            cout << "- MC calc. of temperatures and RATs : done                          "
                    "      "
                 << endl;
            break;

        case CMD_TEMP:
            cout << "- MC calculation of temperatures    : done                          "
                    "      "
                 << endl;
            break;
            
        case  CMD_TEMP_POLY:
            cout << "- MC calculation of polychromatic temperatures    : done                          "
                    "      "
                 << endl;
            break;

        default:
            cout << "- MC calculation of radiation field : done                          "
                    "      "
                 << endl;
            break;
    }
    return true;
}

bool CRadiativeTransfer::setTemperatureDistribution()
{
    ulong per_counter = 0;
    ulong max_cells = grid->getMaxDataCells();

    const double tol = 0.00001;
    const double GNewton = 0.0000000667;

    const double muinv = 4.70758517 * pow(10., 23.);
    const double pi = 3.141592653589793238;
    const double k_B = 1.38065812 * pow(10., -16.);

    cout << CLR_LINE;

#pragma omp parallel for schedule(dynamic)
    for(long c_i = 0; c_i < long(max_cells); c_i++)
    {
        cell_basic * cell = grid->getCellFromIndex(c_i);

        int max_newton = 2000;
        int l = 0;
        double error = 0.0001;

        double rho = dust->getNumberDensity(grid, cell) * 1e-16; // dens_data[c];
        double tdust = 1000; // grid->getGasTemperature(cell); //tdust_orig;
        double tdust_orig = grid->getGasTemperature(cell);
        double chi, tnew, tau, lambda_dust, det, eqfunc, kappa;
        double lambda_J = sqrt(pi * k_B * tdust_orig * muinv / (GNewton * rho));

        while((l < max_newton) && (error > tol))
        {
            chi = 3.3 * pow(10., -26.) * (0.054833588 * tdust) * (0.054833588 * tdust) * rho * muinv;

            if(tdust > 200.)
            {
                if(tdust < 1500.)
                {
                    chi = 3.9688735 * pow(10., -24) * rho * muinv;
                }
                else
                {
                    chi = 3.9688735 * pow(10., -24) * rho * muinv * (tdust * 0.00066666667) * pow(10., -12.);
                }
            }

            kappa = min(chi, 1.0 / lambda_J) * 6.8533297 * pow(10., 28.) / (rho * muinv);
            lambda_dust = kappa * pow(tdust, 4.);
            tau = chi * lambda_J;

            eqfunc = 3.9 * pow(10., 5.) + 0.63245552 * rho * muinv * (tdust_orig - tdust) * sqrt(tdust_orig) -
                     lambda_dust;
            if(tau > 1.0)
            {
                det = -4.0 * lambda_dust / tdust - 0.63245552 * rho * muinv * sqrt(tdust_orig);
            }
            else
            {
                det = -6.0 * lambda_dust / tdust - 0.63245552 * rho * muinv * sqrt(tdust_orig);
                if(tdust > 200.)
                {
                    det = -4.0 * lambda_dust / tdust - 0.63245552 * rho * muinv * sqrt(tdust_orig);
                }
            }

            tnew = tdust - eqfunc / det;

            if(tnew > 2000.)
            {
                tnew = 2000.;
            }

            error = abs(tnew - tdust) / tdust;
            tdust = tnew;
            l++;
        }

        tdust *= 2;

        if(tdust > 2000.)
        {
            tdust = 2000.;
        }

        grid->setDustTemperature(cell, tdust);

        if(adjTgas > 0)
            grid->setGasTemperature(cell, adjTgas * tdust);

#pragma omp critical
        {
            per_counter++;
            if(per_counter % 5000 == 0)
                cout << "-> Estimating of temperatures: " << 100.0 * float(per_counter) / float(max_cells)
                     << " [%] " << dust->getMinDustTemp() << " " << dust->getMaxDustTemp() << "        \r";
        }
    }

    cout << CLR_LINE;
    cout << " " << dust->getMinDustTemp() << " " << dust->getMaxDustTemp() << endl;
    cout << "- Estimation of final temperatures: done" << endl;
    return true;
}

bool CRadiativeTransfer::calcPolMapsViaMC()
{
    // Init variables
    ullong nr_of_photons;
    ullong per_counter, ph_max, nr_of_wavelength;
    float last_percentage;
    uint mrw_counter = 0;
    ullong kill_counter = 0;
    uint max_source = uint(sources_mc.size());

    // Perform Monte-Carlo radiative transfer for each chosen source
    for(uint s = 0; s < max_source; s++)
    {
        // Init source object
        CSourceBasic * tm_source = sources_mc[s];
        cout << CLR_LINE;

        // Get number of photons that have to be emitted from the current source source
        nr_of_photons = tm_source->getNrOfPhotons();

        // Get number of wavelength
        nr_of_wavelength = dust->getNrOfWavelength();

        // Perform Monte-Carlo radiative transfer for
        // each chosen wavelength of the chosen source
        for(uint wID = 0; wID < nr_of_wavelength; wID++)
        {
            // Init source parameters for scattering maps
            if(!tm_source->initSource(wID))
                continue;

            // Init progress visualization
            cout << "-> MC pol. maps (source ID: " << s + 1 << ", wavelength: " << dust->getWavelength(wID)
                 << " [m], photons: " << nr_of_photons << ") 0 [%]   \r" << flush;

            // Init counter and percentage to show progress
            per_counter = 0;
            last_percentage = 0;

            // Perform radiative transfer through the model for each photon
#pragma omp parallel for schedule(dynamic)
            for(llong r = 0; r < llong(nr_of_photons); r++)
            {
                // Init cross sections
                double Cext;

                // Init photon package
                photon_package * rays;

                // If forced scattering is enabled, send out two photons and let
                // at least one of them interact
                if(b_forced)
                {
                    rays = new photon_package[2];
                    ph_max = 1;
                }
                else
                {
                    rays = new photon_package[1];
                    ph_max = 0;
                }

                // Increase counter used to show progress
                per_counter++;

                // Calculate percentage of total progress per source
                float percentage = 100.0 * float(per_counter) / float(nr_of_photons);

                // Show only new percentage number if it changed
                if((percentage - last_percentage) > PERCENTAGE_STEP)
                {
#pragma omp critical
                    {
                        cout << "-> MC pol. maps (source ID: " << s + 1
                             << ", wavelength: " << dust->getWavelength(wID)
                             << " [m], photons: " << nr_of_photons << ") " << percentage << " [%]   \r"
                             << flush;
                        last_percentage = percentage;
                    }
                }

                // Execute the following once (standard mode) or
                // twice (only if enforced scattering is enabled)
                for(uint ph_i = 0; ph_i <= ph_max; ph_i++)
                {
                    // Init the photon_package pp with the specific ray
                    photon_package * pp = &rays[ph_i];

                    // Init variables
                    double end_tau;
                    Vector3D old_pos, start_pos;

                    // If this photon is the first one (if enforced scattering
                    // is enabled, the first one has to scatter once)
                    if(ph_i == 0)
                    {
                        // Set current wavelength
                        pp->setWavelengthID(wID);

                        // Launch a new photon package from the source
                        tm_source->createNextRay(pp, ullong(r));

                        // Position the photon inside the grid
                        if(!grid->positionPhotonInGrid(pp))
                            if(!grid->findStartingPoint(pp))
                            {
                                kill_counter++;
                                break;
                            }

                        // Set starting position for enforced scattering
                        start_pos = pp->getPosition();

                        // Save position of last interaction to know to which pixel
                        // the photon belongs, if it is not scattered on its further
                        // path through the model
                        pp->updatePositionLastInteraction();

                        if(b_forced)
                        {
                            // Get tau for first interaction, if the interaction is forced
                            end_tau = getEscapeTauForced(rays);

                            // If optical depth is exactly zero, send only one photon
                            // package as without enfsca
                            if(end_tau == 0)
                            {
                                ph_max = 0;
                                end_tau = -log(1.0 - pp->getRND());
                            }
                        }
                        else
                        {
                            // Get tau for first interaction
                            end_tau = -log(1.0 - pp->getRND());
                        }
                    }
                    else
                    {
                        // Get tau for first interaction
                        end_tau = -log(1.0 - pp->getRND());
                    }

                    // Init variables
                    ullong interactions = 0;
                    double tmp_tau;
                    double len, dens;

                    // Save photon position to adjust it if necessary
                    old_pos = pp->getPosition();

                    // Transfer photon through grid
                    while(grid->next(pp))
                    {
                        // If max interactions is reached or the photon intensity
                        // is too low, end photon transfer
                        if(interactions >= MAX_INTERACTION || pp->getStokesVector().I() < 1e-200)
                        {
                            if(ph_i == 0)
                                kill_counter++;
                            break;
                        }

                        // Get path length through current cell
                        len = pp->getTmpPathLength();
                        
                        // Update total pathlength of photon if time-dependent
                        if(SCA_DT > 0)
                            pp->updateTotalPathLength(len);

                        // Get dust number density of the current cell
                        dens = dust->getNumberDensity(grid, pp);

                        // If the dust density is too low, skip this cell
                        if(dens < 1e-200)
                        {
                            old_pos = pp->getPosition();
                            continue;
                        }

                        // Calculate the dust absorption cross section (for random
                        // alignment)
                        Cext = dust->getCextMean(grid, pp);

                        // Calculate optical depth that is obtained on the path_length
                        tmp_tau = Cext * len * dens;

                        // optical depth for interaction is reached or not
                        if(tmp_tau > end_tau)
                        {
                            // Increase interaction counter
                            interactions++;

                            // Reduce the photon position to match the exact
                            // interaction position
                            
                            // Correct total pathlength if necessary
                            if(SCA_DT > 0)
                                pp->updateTotalPathLength(-len*(1 - end_tau / tmp_tau));

                            len = len * end_tau / tmp_tau;
                            pp->adjustPosition(old_pos, len);

                            // Modify second photon if enforced scattering is used
                            if(b_forced && interactions == 1 && ph_i == 0)
                            {
                                rays[1].initRandomGenerator(int(r * pp->getRND()));
                                rays[1].setWavelengthID(pp->getWavelengthID());
                                rays[1].setPosition(pp->getPosition());
                                rays[1].setPositionCell(pp->getPositionCell());
                                rays[1].setPositionLastInteraction(start_pos);
                                rays[1].setDirection(pp->getDirection());
                                rays[1].initCoordSystem();
                                pp->setStokesVector(pp->getStokesVector() - rays[1].getStokesVector());
                            }

                            // Save position of last interaction to know to which pixel
                            // the photon belongs, if it is not scattered on its further
                            // path through the model
                            pp->updatePositionLastInteraction();

                            if(!doMRWStepBWWithoutHeating(pp))
                            {
                                // Calculate the dust scattering cross section (for random
                                // alignment)
                                // Csca = dust->getCscaMean(grid, pp);

                                // If peel-off is used, add flux to the detector
                                if(peel_off)
                                {
                                    // Transport a separate photon to each detector
                                    for(uint d = 0; d < nr_mc_detectors; d++)
                                    {
                                        // Get index of wavelength in current detector
                                        uint wID_det =
                                            detector[d].getDetectorWavelengthID(dust->getWavelength(wID));

                                        // Only calculate for detectors with the
                                        // corresponding wavelengths
                                        if(wID_det != MAX_UINT)
                                        {
                                            // Create an escaping photon into the
                                            // direction of the detector
                                            photon_package tmp_pp = dust->getEscapePhoton(
                                                grid, pp, detector[d].getEX(), detector[d].getDirection());

                                            // Convert the flux into Jy and consider
                                            // the distance to the observer
                                            CMathFunctions::lum2Jy(tmp_pp.getStokesVector(),
                                                                   dust->getWavelength(wID),
                                                                   detector[d].getDistance());

                                            // Consider foreground extinction
                                            tmp_pp.getStokesVector() *=
                                                dust->getForegroundExtinction(dust->getWavelength(wID));
                                                
                                            if(SCA_DT > 0)
                                            {
                                                // Add path length from grid to detector to photon
                                                tmp_pp.updateTotalPathLength(((detector[d].getDistance()*detector[d].getDirection()-pp->getPosition())*detector[d].getDirection())+pp->getTotalPathLength());
                                                // Add the photon package to the detector if in time
                                                if(tmp_pp.getTotalPathLength()/con_c+s*SCA_DT < (SCA_DT*(d+1)+detector[d].getDistance()/con_c))
                                                {
                                                    detector[d].addToMonteCarloDetector(
                                                        &tmp_pp, wID_det, SCATTERED_DUST);
                                                    break;
                                                }
                                            }
                                            else
                                                // Add the photon package to the detector
                                                detector[d].addToMonteCarloDetector(
                                                    &tmp_pp, wID_det, SCATTERED_DUST);
                                        }
                                    }
                                }
                            }
                            else
                            {
                                mrw_counter++;
                                if(mrw_counter % 500000 == 0)
                                {
#pragma omp critical
                                    {
                                        cout << " -> mrw  " << mrw_counter << "\r";
                                    }
                                }
                            }

                            // Perform scattering interaction into a new direction
                            dust->scatter(grid, pp, true);

                            // Calculate new optical depth for next interaction
                            end_tau = -log(1.0 - pp->getRND());
                        }
                        else
                        {
                            // Remove the traveled distance from optical depth
                            end_tau -= tmp_tau;
                        }
                        // Save photon position to adjust it if necessary
                        old_pos = pp->getPosition();
                    }

                    // If peel-off is not used, use classic Monte-Carlo method
                    // Now, the photon has left the model space
                    if(!peel_off && pp->getStokesVector().I() > 1e-200 && interactions <= MAX_INTERACTION)
                    {
                        // Move photon back to the point of last interaction
                        pp->resetPositionToLastInteraction();

                        // Transport photon to observer for each detector
                        for(uint d = 0; d < nr_mc_detectors; d++)
                        {
                            // Get index of wavelength in current detector
                            uint wID_det = detector[d].getDetectorWavelengthID(dust->getWavelength(wID));

                            // Only calculate for detectors with the corresponding
                            // wavelengths
                            if(wID_det != MAX_UINT)
                            {
                                // Get direction to the current detector
                                Vector3D dir_obs = detector[d].getDirection();

                                // Calculate the angle between the photon travel direction
                                // and the detector direction
                                double cos_angle = pp->getDirection() * dir_obs;

                                // Get acceptance angle from detector (minimum 1°)
                                double cos_acceptance_angle = detector[d].getAcceptanceAngle();

                                // If the angle angle between the photon direction and the
                                // detector direction is smaller than the acceptance angle
                                if(cos_angle >= cos_acceptance_angle)
                                {
                                    // Get the angle to rotate the photon space into the
                                    // detector space
                                    double rot_angle_phot_obs = CMathFunctions::getRotationAngleObserver(
                                        detector[d].getEX(), pp->getEX(), pp->getEY());

                                    // Rotate the Stokes vector
                                    pp->getStokesVector().rot(rot_angle_phot_obs);

                                    // Consider the greater solid angle due
                                    // to the acceptance angle
                                    pp->getStokesVector() *= 1.0 / ((1.0 - cos_acceptance_angle) * PIx2);

                                    // Convert the flux into Jy and consider
                                    // the distance to the observer
                                    CMathFunctions::lum2Jy(pp->getStokesVector(),
                                                           dust->getWavelength(wID),
                                                           detector[d].getDistance());

                                    // Consider foreground extinction
                                    pp->getStokesVector() *=
                                        dust->getForegroundExtinction(dust->getWavelength(wID));

                                    if(interactions == 0)
                                    {
                                        if(SCA_DT > 0)
                                        {
                                            // Add path length from grid to detector to photon
                                            pp->updateTotalPathLength((detector[d].getDistance()*detector[d].getDirection()-pp->getPosition())*detector[d].getDirection());
                                            // Add the photon package to the detector if in time
                                            if(pp->getTotalPathLength()/con_c+s*SCA_DT < (SCA_DT*(d+1)+detector[d].getDistance()/con_c))
                                                detector[d].addToMonteCarloDetector(pp, wID_det, DIRECT_STAR);
                                        }
                                        else
                                            // Add the photon package to the detector
                                            detector[d].addToMonteCarloDetector(pp, wID_det, DIRECT_STAR);
                                    }
                                    else
                                    {
                                        if(SCA_DT > 0)
                                        {
                                            // Add path length from grid to detector to photon
                                            pp->updateTotalPathLength((detector[d].getDistance()*detector[d].getDirection()-pp->getPosition())*detector[d].getDirection());
                                            // Add the photon package to the detector if in time
                                            if(pp->getTotalPathLength()/con_c+s*SCA_DT < (SCA_DT*(d+1)+detector[d].getDistance()/con_c))
                                                detector[d].addToMonteCarloDetector(pp, wID_det, SCATTERED_DUST);
                                        }
                                        else
                                            // Add the photon package to the detector
                                            detector[d].addToMonteCarloDetector(pp, wID_det, SCATTERED_DUST);
                                    }
                                }
                            }
                        }
                    }
                    // Do not launch the secondary photon if the first one has not
                    // scattered
                    if(b_forced && interactions == 0 && ph_i == 0)
                        break;
                }
                // Delete the Rays pointer
                delete[] rays;
            }

            // If peel-off is used, add direct source emission to each detector
            if(peel_off && tm_source->getStringID() != "dust source")
            {
                // Transport photon to observer for each detector
                for(uint d = 0; d < nr_mc_detectors; d++)
                {
                    // Get index of wavelength in current detector
                    uint wID_det = detector[d].getDetectorWavelengthID(dust->getWavelength(wID));

                    // Only calculate for detectors with the corresponding wavelengths
                    if(wID_det != MAX_UINT)
                    {
                        // Create temporary photon package
                        photon_package tmp_pp;

                        // Set current wavelength to temporary photon package
                        tmp_pp.setWavelengthID(wID);

                        // Get direction to the current detector
                        Vector3D dir_obs = detector[d].getDirection();

                        // Launch a new photon package from the source
                        tm_source->createDirectRay(&tmp_pp, dir_obs);

                        // Position the photon inside the grid
                        grid->positionPhotonInGrid(&tmp_pp);

                        // Init a variable to save the optical depth
                        double tau_obs = 0;

                        // Transport photon package through model to obtain the optical
                        // depth
                        while(grid->next(&tmp_pp))
                        {
                            // Get necessary information about the current cell
                            double len = tmp_pp.getTmpPathLength();
                            double dens = dust->getNumberDensity(grid, &tmp_pp);
                            double Cext = dust->getCextMean(grid, &tmp_pp);

                            // Increase the optical depth by the current cell
                            tau_obs += Cext * len * dens;
                        }

                        // Rotate photon package into the coordinate space of the detector
                        double rot_angle_phot_obs = CMathFunctions::getRotationAngleObserver(
                            detector[d].getEX(), tmp_pp.getEX(), tmp_pp.getEY());
                        tmp_pp.getStokesVector().rot(rot_angle_phot_obs);

                        // Calculate the source emission and reduce it by the optical
                        // depth
                        tmp_pp.getStokesVector() *= exp(-tau_obs);

                        // Convert the flux into Jy and consider the distance to the
                        // observer
                        CMathFunctions::lum2Jy(
                            tmp_pp.getStokesVector(), dust->getWavelength(wID), detector[d].getDistance());

                        // Consider foreground extinction
                        tmp_pp.getStokesVector() *= dust->getForegroundExtinction(dust->getWavelength(wID));

                        if(SCA_DT > 0)
                        {
                            // Add total pathlength to detector if time-dependent
                            Vector3D dir_obs = detector[d].getDirection();
                            tmp_pp.updateTotalPathLength((detector[d].getDistance()*detector[d].getDirection()-tm_source->getPosition())*dir_obs);
                            // Add the photon package to the detector if in time
                            if(tmp_pp.getTotalPathLength()/con_c+s*SCA_DT < (SCA_DT*(d+1)+detector[d].getDistance()/con_c))
                            {
                                detector[d].addToMonteCarloDetector(&tmp_pp, wID_det, DIRECT_STAR);
                                break;
                            }
                        }
                        else
                            // Add the photon package to the detector
                            detector[d].addToMonteCarloDetector(&tmp_pp, wID_det, DIRECT_STAR);
                    }
                }
            }
        }
        if(peel_off && tm_source->getStringID() == "dust source")
        {
            cout << CLR_LINE;
            cout << "\nHINT: MC simulations with dust source and peel-off include only "
                    "the scattered radiation.\n"
                 << "Add results from Raytracing simulations for full dust emission!" << endl;
        }
    }

    // Write results either as text or fits file
    for(uint d = 0; d < nr_mc_detectors; d++)
    {
        if(!detector[d].writeMap(d, RESULTS_MC))
            return false;
        if(!detector[d].writeSed(d, RESULTS_MC))
            return false;
    }

    // Format prints
    cout << CLR_LINE;
    cout << SEP_LINE;

    if(mrw_step)
        cout << "- MRW steps                           : " << mrw_counter << endl;

    // Show amount of killed photons
    if(kill_counter > 0)
        cout << "- Photons killed                   : " << kill_counter << endl;
    cout << "- Calculation of MC polarization maps (photons: " << nr_of_photons << "): done" << endl;

    return true;
}

void CRadiativeTransfer::convertTempInQB(double min_gas_density, bool use_gas_temp)
{
    ulong max_cells = grid->getMaxDataCells();
    ulong pos_counter = 0;

    cout << "-> Converting emissivities: 0 [%]        \r" << flush;

#pragma omp parallel for
    for(long c_i = 0; c_i < long(max_cells); c_i++)
    {
        cell_basic * cell = grid->getCellFromIndex(c_i);

        // Get gas density of current cell
        double gas_dens = grid->getGasNumberDensity(cell);
        if(gas_dens < min_gas_density)
            continue;

        dust->convertTempInQB(grid, cell, min_gas_density, use_gas_temp);

        pos_counter++;

        if(pos_counter % 10000 == 0)
        {
#pragma omp critical
            {
                cout << "-> Converting emissivities: " << 100.0 * float(pos_counter) / float(max_cells)
                     << " [%]       \r";
            }
        }
    }

    cout << CLR_LINE;
    // cout << "- Converting emissivities             : done" << endl;
}

void CRadiativeTransfer::calcAlignedRadii()
{
    ulong max_cells = grid->getMaxDataCells();

    ullong per_counter = 0;
    float last_percentage = 0;

    cout << CLR_LINE;
    cout << " -> Calc. RAT dust alig. radius: 0.0 [%]  (min: 0 [m]; max: 0 [m])        \r" << flush;

#pragma omp parallel for schedule(dynamic)
    for(long c_i = 0; c_i < long(max_cells); c_i++)
    {
        cell_basic * cell = grid->getCellFromIndex(c_i);
        dust->calcAlignedRadii(grid, cell);

        per_counter++;
        float percentage = 100.0 * float(per_counter) / float(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
#pragma omp critical
            {
                cout << " -> Calc. RAT dust alig. radius: " << 100.0 * float(per_counter) / float(max_cells)
                     << " [%] (min: " << dust->getMinAlignedRadius()
                     << " [m]; max: " << dust->getMaxAlignedRadius() << " [m])"
                     << "          \r";
                last_percentage = percentage;
            }
        }
    }

    cout << CLR_LINE;
    cout << "- Calculation of aligned radii      : done" << endl;
    cout << "    aligned radii [" << float(dust->getMinAlignedRadius()) << ", "
         << float(dust->getMaxAlignedRadius()) << "] [m]" << endl;
}

void CRadiativeTransfer::calcFinalTemperature(bool use_energy_density)
{
    ullong per_counter = 0;
    float last_percentage = 0;

    ulong max_cells = grid->getMaxDataCells();

    cout << CLR_LINE;
    cout << "-> Calculation of final temperatures : 0.0[%]                        \r";

#pragma omp parallel for schedule(dynamic)
    for(long c_i = 0; c_i < long(max_cells); c_i++)
    {
        cell_basic * cell = grid->getCellFromIndex(c_i);
        dust->calcTemperature(grid, cell, use_energy_density);

        if(adjTgas > 0)
            grid->setGasTemperature(cell, adjTgas * dust->getDustTemperature(grid, cell));

        per_counter++;
        float percentage = 100.0 * float(per_counter) / float(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
#pragma omp critical
            {
                cout << "-> Calculation of final temperatures : "
                     << 100.0 * float(per_counter) / float(max_cells) << " [%]              \r";

                last_percentage = percentage;
            }
        }
    }

    cout << CLR_LINE;
    cout << "- Calculation of final temperatures : done" << endl;
}

void CRadiativeTransfer::calcStochasticHeating()
{
    ullong per_counter = 0;
    float last_percentage = 0;

    // Resize wavelength list for stochastic heating propabilities
    dlist wl_list(WL_STEPS);
    CMathFunctions::LogList(WL_MIN, WL_MAX, wl_list, 10);

    ulong max_cells = grid->getMaxDataCells();

    cout << CLR_LINE;
    cout << "-> Calculation of stochastic heating: 0.0[%]    \r" << flush;

#pragma omp parallel for schedule(dynamic)
    for(long c_i = 0; c_i < long(max_cells); c_i++)
    {
        cell_basic * cell = grid->getCellFromIndex(c_i);

        dust->calcStochasticHeatingPropabilities(grid, cell, wl_list);

        per_counter++;
        float percentage = 100.0 * float(per_counter) / float(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
#pragma omp critical
            {
                cout << "-> Calculation of stochastic heating: "
                     << 100.0 * float(per_counter) / float(max_cells) << " [%]          \r" << flush;
                last_percentage = percentage;
            }
        }
    }

    cout << CLR_LINE;
    cout << "- Calculation of stochastic heating    : done" << endl;
}

double CRadiativeTransfer::getEscapeTauForced(photon_package * rays)
{
    double len, dens, Cext, enf_tau = 0;
    double rnd = rays[0].getRND();
    StokesVector stokes;
    double factor;

    // Save the old position to use it again
    Vector3D old_pos = rays[0].getPosition();
    cell_basic * tmp_cell = rays[0].getPositionCell();

    while(grid->next(&rays[0]))
    {
        len = rays[0].getTmpPathLength();
        dens = dust->getNumberDensity(grid, &rays[0]);
        Cext = dust->getCextMean(grid, &rays[0]);
        enf_tau += Cext * len * dens;
    }

    // Reset the photon position
    rays[0].setPosition(old_pos);
    rays[0].setPositionCell(tmp_cell);

    stokes = rays[0].getStokesVector();
    factor = exp(-enf_tau);

    rays[0].setStokesVector(stokes);
    rays[1].setStokesVector(stokes * factor);

    enf_tau *= 0.999;
    return -log(1.0 - rnd * (1.0 - exp(-enf_tau)));
}

// -------------------------------------------------
// ------ Calculation of Synchotron transfer -------
// -------------------------------------------------

bool CRadiativeTransfer::calcSyncMapsViaRaytracing(parameters & param)
{
    // Init math environment
    CMathFunctions mf;

    // Get list of detectors/sequences that will be simulated with the raytracer
    dlist sync_ray_detectors = param.getSyncRayDetectors();

    // Get maximum length of the simulation model
    double max_length = grid->getMaxLength();

    if(!sync_ray_detectors.empty())
    {
        for(uint i_det = start; i_det <= stop; i_det++)
        {
            // Init source object
            CSourceBasic * tmp_source;
            uint sID = tracer[i_det]->getSourceIndex();
            tmp_source = sources_ray[sID];
            tmp_source->setSideLength(max_length);

            // Calculate total number of pixel
            uint per_max = tracer[i_det]->getNpix();

            // Init counter and percentage to show progress
            ullong per_counter = 0;
            float last_percentage = 0;

            // Show information about the current detector
            cout << CLR_LINE;
            cout << "-> Raytracing synchrotron maps (Seq. " << i_det + 1 << ", source: " << sID + 1
                 << ") 0.0 [%]   \r" << flush;

            // Calculate pixel intensity for each pixel
#pragma omp parallel for schedule(dynamic)
            for(int i_pix = 0; i_pix < int(per_max); i_pix++)
            {
                double cx = 0, cy = 0;
                if(!tracer[i_det]->getRelPosition(i_pix, cx, cy))
                    continue;

                getSyncPixelIntensity(tmp_source, cx, cy, i_det, 0, i_pix);

                // Increase counter used to show progress
                per_counter++;

                // Calculate percentage of total progress per source
                float percentage = 100.0 * float(per_counter) / float(per_max);

                // Show only new percentage number if it changed
                if((percentage - last_percentage) > PERCENTAGE_STEP)
                {
#pragma omp critical
                    {
                        cout << "-> Raytracing synchrotron maps (Seq. " << i_det + 1
                             << ", source: " << sID + 1 << ")  "
                             << float(100.0 * float(per_counter) / float(per_max)) << " [%]         \r"
                             << flush;
                        last_percentage = percentage;
                    }
                }
            }

            // Show final progress
            cout << "-> Raytracing synchrotron maps (Seq. " << i_det + 1 << ", source: " << sID + 1
                 << ") 100 [%]       \r" << flush;

            // post-process raytracing simulation
            if(!tracer[i_det]->postProcessing())
                return false;

            // Write results either as text or fits file
            if(!tracer[i_det]->writeSyncResults())
                return false;
        }
    }

    // Show that raytracing is finished
    cout << CLR_LINE;
    cout << "- Raytracing synchrotron map    : done" << endl;

    return true;
}

void CRadiativeTransfer::getSyncPixelIntensity(CSourceBasic * tmp_source,
                                               double cx,
                                               double cy,
                                               uint i_det,
                                               uint subpixel_lvl,
                                               int i_pix)
{
    bool subpixel = false;

    subpixel = tracer[i_det]->getUseSubpixel(cx, cy, subpixel_lvl);

    // If any subpixel traveled through other cells, perform subpixelling
    if(subpixel == false)
    {
        // Init variables
        uint nr_used_wavelengths = tracer[i_det]->getNrSpectralBins();
        // cout << i_pix << endl;

        // Create new photon package
        photon_package * pp1 = new photon_package();
        photon_package * pp2 = new photon_package();

        // Init photon package
        pp1->initMultiStokesVector(nr_used_wavelengths);
        pp2->initMultiStokesVector(nr_used_wavelengths);

        tracer[i_det]->preparePhoton(pp1, cx, cy);
        tracer[i_det]->preparePhoton(pp2, cx, cy);

        // Calculate continuum emission along one path
        getSyncIntensity(pp1, pp2, tmp_source, cx, cy, i_det, subpixel_lvl);

        tracer[i_det]->addToDetector(pp1, pp2, i_pix);

        // Delete photon package after usage
        delete pp1;
        delete pp2;
    }
    else
    {
        // Repeat this function for each subpixel
        for(int i_sub_x = -1; i_sub_x <= 1; i_sub_x += 2)
        {
            for(int i_sub_y = -1; i_sub_y <= 1; i_sub_y += 2)
            {
                // Calculate positions of each subpixel
                double tmp_cx, tmp_cy;
                tracer[i_det]->getSubPixelCoordinates(subpixel_lvl, cx, cy, i_sub_x, i_sub_y, tmp_cx, tmp_cy);
                // Calculate radiative transfer of the current pixel
                // and add it to the detector at the corresponding position
                getSyncPixelIntensity(tmp_source, tmp_cx, tmp_cy, i_det, (subpixel_lvl + 1), i_pix);
            }
        }
    }
}

void CRadiativeTransfer::getSyncIntensity(photon_package * pp1,
                                          photon_package * pp2,
                                          CSourceBasic * tmp_source,
                                          double cx,
                                          double cy,
                                          uint i_det,
                                          uint subpixel_lvl)
{
    // Set amount of radiation coming from this pixel
    double subpixel_fraction = pow(4.0, -double(subpixel_lvl));

    // Get chosen wavelength parameter
    uint nr_used_wavelengths = tracer[i_det]->getNrSpectralBins();

    // Init multiple Stokes vectors
    MultiStokesVector WMap_cr(nr_used_wavelengths);
    MultiStokesVector WMap_ca(nr_used_wavelengths);

    // Update Stokes vectors with emission from background source
    for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
    {
        // Get current wavelength index
        uint wID = dust->getWavelengthID(tracer[i_det]->getWavelength(i_wave));

        // Set wavelength index in photon package
        pp1->setWavelengthID(wID);

        // Set related index of the multi wavelength map to the wavelength ID

        WMap_cr.setS(tmp_source->getStokesVector(pp1), i_wave);

        WMap_ca.setS(tmp_source->getStokesVector(pp1), i_wave);
    };

    // Find starting point inside the model and travel through it
    if(grid->findStartingPoint(pp1))
    {
        while(grid->next(pp1) && tracer[i_det]->isNotAtCenter(pp1, cx, cy))
        {
            double n_th = grid->getThermalElectronDensity(pp1);
            double T_e = 0.0; // grid->getElectronTemperature(pp1); //reserved for later use

            double n_cr = grid->getCRElectronDensity(pp1);
            double g_min = grid->getGammaMin(pp1);
            double g_max = grid->getGammaMax(pp1);
            double pow_p = grid->getPowerLawIndex(pp1);

            double B = grid->getMagField(pp1).length();

            // If the all the electron densities are far too low, skip the current cell
            if(n_cr + n_th >= 1e-200)
            {
                // Get path length through current cell
                double len = pp1->getTmpPathLength();

                // Calculate orientation of the Stokes vector in relation to the detector
                double phi = grid->getPhiMag(pp1);
                double theta = grid->getThetaMag(pp1);

                double sin_2ph = sin(2.0 * phi);
                double cos_2ph = cos(2.0 * phi);

                for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
                {
                    // Get wavelength/frequency of the photon package
                    double wavelength = tracer[i_det]->getWavelength(i_wave);

                    syn_param syn_th = synchrotron->get_Thermal_Parameter(n_th, T_e, wavelength, B, theta);
                    syn_param syn_cr =
                        synchrotron->get_Power_Law_Parameter(n_cr, wavelength, B, theta, g_min, g_max, pow_p);
                    syn_param syn_ca = syn_cr + syn_th;

                    // Get matrixes with the absorption's and conversions
                    Matrix2D alpha_cr = -1 * syn_cr.getSyncMatrix();
                    Matrix2D alpha_ca = -1 * syn_ca.getSyncMatrix();

                    // Init a variable to sum up path lengths until cell is crossed
                    double cell_sum = 0.0;

                    // First path length is path through cell
                    double cell_d_l = len;

                    // Save direction of the photon package
                    Vector3D dir_map_xyz = pp1->getDirection();

                    // Save the cell entry position of the photon package
                    // (Hint: next(pp) puts the photon onto the border to the next cell)
                    Vector3D pos_xyz_cell = pp1->getPosition() - (len * dir_map_xyz);

                    // Set stokes vector with emission
                    StokesVector S_em_cr(syn_cr.j_I, syn_cr.j_Q * cos_2ph, syn_cr.j_Q * sin_2ph, syn_cr.j_V);

                    StokesVector S_em_ca(syn_ca.j_I, syn_ca.j_Q * cos_2ph, syn_ca.j_Q * sin_2ph, syn_ca.j_V);

                    // Initiating the kill counter
                    ullong kill_counter = 0;

                    // Initiating the fail-save
                    bool fail = false;

                    // Make sub steps until cell is completely crossed
                    // If the error of a sub step is too high, make the step smaller
                    while(cell_sum < len)
                    {
                        // Increase the kill counter
                        kill_counter++;

                        // If too many sub steps are needed, kill the photon
                        if(kill_counter > 2.0 * MAX_SOLVER_STEPS)
                        {
#pragma omp critical
                            {
                                cout << CLR_LINE;
                                cout << "\nWARNING: Solver steps > " << 2.0 * MAX_SOLVER_STEPS
                                     << ". Too many steps!" << endl
                                     << flush;
                                cout << "         Skipping entire cell!" << endl << flush;
                            }
                            break;
                        }

                        if(kill_counter > 1.0 * MAX_SOLVER_STEPS)
                        {
                            if(fail == false)
                            {
                                fail = true;
#pragma omp critical
                                {
                                    cout << CLR_LINE;
                                    cout << "\nWARNING: Solver steps > " << 1.0 * MAX_SOLVER_STEPS
                                         << ". Too many steps!" << endl
                                         << flush;
                                    cout << "         Switching to approximate solver!" << endl << flush;
                                }
                            }
                        }

                        // Init Runge-Kutta parameters and set it to zero
                        // (see https://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method)
                        StokesVector * RK_k_cr = new StokesVector[6];
                        StokesVector * RK_k_ca = new StokesVector[6];

                        // Calculate result of the radiative transfer equation at each
                        // Runge-Kutta sub position
                        for(uint k = 0; k < 6; k++)
                        {
                            // Init scalar product
                            StokesVector scalar_product_cr, scalar_product_ca;

                            // Calculate multiplication between Runge-Kutta parameter
                            for(uint i = 0; i < 6; i++)
                            {
                                scalar_product_cr += (RK_k_cr[i] * RK_a(i, k));
                                scalar_product_ca += (RK_k_ca[i] * RK_a(i, k));
                            }

                            // Calculate new Runge-Kutta parameters as the result of the
                            // radiative transfer equation at the Runge-Kutta sub
                            // positions
                            RK_k_cr[k] =
                                alpha_cr * (scalar_product_cr * cell_d_l + WMap_cr.S(i_wave)) + S_em_cr;
                            RK_k_ca[k] =
                                alpha_ca * (scalar_product_ca * cell_d_l + WMap_ca.S(i_wave)) + S_em_ca;
                        }

                        // Init two temporary Stokes vectors
                        StokesVector stokes_new_cr = WMap_cr.S(i_wave);
                        StokesVector stokes_new2_cr = WMap_cr.S(i_wave);

                        StokesVector stokes_new_ca = WMap_ca.S(i_wave);
                        StokesVector stokes_new2_ca = WMap_ca.S(i_wave);

                        // Add the result at each Runge-Kutta sub position
                        // to the total Stokes vector (Using Runge-Kutta Fehlberg)
                        for(uint i = 0; i < 6; i++)
                        {
                            stokes_new_cr += (RK_k_cr[i] * cell_d_l) * RK_b1[i];
                            stokes_new2_cr += (RK_k_cr[i] * cell_d_l) * RK_b2[i];

                            if(!fail)
                            {
                                stokes_new_ca += (RK_k_ca[i] * cell_d_l) * RK_b1[i];
                                stokes_new2_ca += (RK_k_ca[i] * cell_d_l) * RK_b2[i];
                            }
                        }

                        // Delete the Runge-Kutta pointer
                        delete[] RK_k_cr;
                        delete[] RK_k_ca;

                        // Calculate the difference between the results with two
                        // different precisions to see if smaller steps are needed
                        // (see Reissl)

                        double epsi, dz_new;
                        if(stokes_new_cr.I() < 0 || stokes_new2_cr.I() < 0 || stokes_new_ca.I() < 0 ||
                           stokes_new2_ca.I() < 0)
                        {
                            epsi = 2.0;
                            dz_new = 0.0;
                        }
                        else
                        {
                            double epsi_cr = abs(stokes_new2_cr.I() - stokes_new_cr.I()) /
                                             (rel_err * abs(stokes_new_cr.I()) + abs_err);
                            double epsi_ca = abs(stokes_new2_ca.I() - stokes_new_ca.I()) /
                                             (rel_err * abs(stokes_new_ca.I()) + abs_err);

                            double epsi_Q = abs(abs(stokes_new2_ca.Q()) - abs(stokes_new_ca.Q())) /
                                            (rel_err * abs(stokes_new_ca.Q()) + abs_err);

                            double epsi_U = abs(abs(stokes_new2_ca.U()) - abs(stokes_new_ca.U())) /
                                            (rel_err * abs(stokes_new_ca.U()) + abs_err);

                            double epsi_V = abs(abs(stokes_new2_ca.V()) - abs(stokes_new_ca.V())) /
                                            (rel_err * abs(stokes_new_ca.V()) + abs_err);

                            double dz_new_cr = 0.9 * cell_d_l * pow(epsi_cr, -0.2);
                            double dz_new_ca = 0.9 * cell_d_l * pow(epsi_ca, -0.2);

                            double dz_new_Q = 0.9 * cell_d_l * pow(epsi_Q, -0.2);
                            double dz_new_U = 0.9 * cell_d_l * pow(epsi_U, -0.2);
                            double dz_new_V = 0.9 * cell_d_l * pow(epsi_V, -0.2);

                            // Do approximate solution
                            if(fail)
                            {
                                epsi = epsi_cr;
                                dz_new = dz_new_cr;

                                double j_Q = syn_ca.j_Q * cos_2ph;
                                double j_U = syn_ca.j_Q * sin_2ph;
                                double j_V = syn_ca.j_V;

                                double k_Q = syn_ca.kappa_Q;
                                double k_V = syn_ca.kappa_V;

                                double k = sqrt(k_Q * k_Q + k_V * k_V);
                                double k_2 = k * k;
                                double k_3 = k * k * k;

                                double rot_ang = k * cell_d_l;

                                double Q = k_Q / k_2 * (j_Q * k_Q + j_V * k_V) * cell_d_l;
                                Q += k_V / k_3 * (j_V * k_Q + j_Q * k_V) * sin(rot_ang);
                                Q -= j_U * k_V / k_2 * (1.0 - cos(rot_ang));

                                double U = (j_Q * k_V - j_V * k_Q) / k_2 * (1.0 - cos(rot_ang));
                                U += j_U / k * sin(rot_ang);

                                double V = k_V / k_2 * (j_Q * k_Q + j_V * k_V) * cell_d_l;
                                V += k_V / k_3 * (j_Q * k_V + j_V * k_Q) * sin(rot_ang);
                                V -= j_U * k_Q / k_2 * (1.0 - cos(rot_ang));

                                stokes_new_ca = stokes_new_cr + StokesVector(0, Q, U, V);
                            }
                            else
                            {
                                epsi = max(epsi_cr, max(epsi_ca, max(epsi_Q, max(epsi_U, epsi_V))));
                                dz_new =
                                    min(dz_new_cr, min(dz_new_ca, min(dz_new_Q, min(dz_new_U, dz_new_V))));
                            }

                            if(epsi == 0)
                            {
                                dz_new = 0.01 * len;
                            }
                        }

                        // Is a smaller step width needed
                        if(epsi <= 1.0)
                        {
                            // Add additional data to stokes vector
                            stokes_new_cr.setT(WMap_cr.T(i_wave));
                            stokes_new_cr.setSp(WMap_cr.Sp(i_wave));

                            // Add the temporary Stokes vector to the total one
                            WMap_cr.setS(stokes_new_cr, i_wave);

                            // tau is set to Farraday rotation*lambda^2 for CR electrons
                            // (currently its set to zero)
                            WMap_cr.addT(syn_cr.kappa_V * cell_d_l, i_wave);

                            // Sp is set to the thermal electron column
                            WMap_cr.addSp(n_th * cell_d_l, i_wave); //

                            // Add additional data to stokes vector
                            stokes_new_ca.setT(WMap_ca.T(i_wave));
                            stokes_new_ca.setSp(WMap_ca.Sp(i_wave));

                            WMap_ca.setS(stokes_new_ca, i_wave);

                            // tau is set to Farraday rotation*lambda^2 for all
                            WMap_ca.addT(syn_ca.kappa_V * cell_d_l, i_wave);

                            // Sp is set to the CR electron column
                            WMap_ca.addSp(n_cr * cell_d_l, i_wave); //

                            // Update the position of the photon package
                            pos_xyz_cell += cell_d_l * dir_map_xyz;

                            // Increase the sum of the cell path lengths
                            cell_sum += cell_d_l;

                            // Find a new path length
                            cell_d_l = min(dz_new, 4 * cell_d_l);

                            // If the new step would exceed the cell, make it smaller
                            if(cell_sum + cell_d_l > len)
                                cell_d_l = len - cell_sum;
                        }
                        else
                        {
                            // Find a smaller path length
                            cell_d_l = max(dz_new, 0.25 * cell_d_l);
                        }
                    }
                }
            }
        }
    }

    // Update the multi Stokes vectors for each wavelength
    for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
    {
        // Get wavelength/frequency of the photon package
        double wavelength = tracer[i_det]->getWavelength(i_wave);
        double frequency = con_c / wavelength;
        double mult = 1.0e+26 * subpixel_fraction * tracer[i_det]->getDistanceFactor() * con_c /
                      (frequency * frequency);

        // Include foreground extinction if necessary
        mult *= dust->getForegroundExtinction(tracer[i_det]->getWavelength(wavelength));

        // Update the photon package with the multi Stokes vectors

        if(WMap_cr.S(i_wave).I() < 0)
            WMap_cr.S(i_wave).setI(0);

        WMap_cr.S(i_wave) *= mult;
        WMap_cr.setT(WMap_cr.T(i_wave) * subpixel_fraction, i_wave);
        WMap_cr.setSp(WMap_cr.Sp(i_wave) * subpixel_fraction, i_wave);

        pp1->setMultiStokesVector(WMap_cr.S(i_wave), i_wave);

        if(WMap_ca.S(i_wave).I() < 0)
            WMap_ca.S(i_wave).setI(0);

        WMap_ca.S(i_wave) *= mult;
        WMap_ca.setT(WMap_ca.T(i_wave) * subpixel_fraction, i_wave);
        WMap_ca.setSp(WMap_ca.Sp(i_wave) * subpixel_fraction, i_wave);

        pp2->setMultiStokesVector(WMap_ca.S(i_wave), i_wave);
    }
}

// -------------------------------------------
// ------ Calculation of Dust transfer -------
// -------------------------------------------

bool CRadiativeTransfer::calcPolMapsViaRaytracing(parameters & param)
{
    // Get list of detectors/sequences that will be simulated with the raytracer
    dlist dust_ray_detectors = param.getDustRayDetectors();

    // Get maximum length of the simulation model
    double max_length = grid->getMaxLength();

    if(!dust_ray_detectors.empty())
    {
        for(uint i_det = start; i_det <= stop; i_det++)
        {
            // Init source object
            CSourceBasic * tmp_source;
            uint sID = tracer[i_det]->getSourceIndex();
            tmp_source = sources_ray[sID];

            // Calculate total number of pixel
            uint per_max = tracer[i_det]->getNpix();

            // Init counter and percentage to show progress
            ullong per_counter = 0;
            float last_percentage = 0;

            // Show information about the current detector
            cout << CLR_LINE;
            cout << "-> Raytracing dust maps (Seq. " << i_det + 1 << ", source: " << sID + 1 << ") 0 [%]   \r"
                 << flush;

            // Calculate pixel intensity for each pixel
#pragma omp parallel for schedule(dynamic)
            for(int i_pix = 0; i_pix < int(per_max); i_pix++)
            {
                double cx = 0, cy = 0;
                if(!tracer[i_det]->getRelPosition(i_pix, cx, cy))
                    continue;

                getDustPixelIntensity(tmp_source, cx, cy, i_det, 0, i_pix);

                // Increase counter used to show progress
                per_counter++;

                // Calculate percentage of total progress per source
                float percentage = 100.0 * float(per_counter) / float(per_max);

                // Show only new percentage number if it changed
                if((percentage - last_percentage) > PERCENTAGE_STEP)
                {
#pragma omp critical
                    {
                        cout << "-> Raytracing dust maps (Seq. " << i_det + 1 << ", source: " << sID + 1
                             << ") " << percentage << " [%]       \r" << flush;
                        last_percentage = percentage;
                    }
                }
            }

            // Include stellar emission, if chosen
            if(tracer[i_det]->considerPointSources())
                calcStellarEmission(i_det);

            // Show final progress
            cout << "-> Raytracing dust maps (Seq. " << i_det + 1 << ", source: " << sID + 1
                 << ") 100 [%]       \r" << flush;

            // post-process raytracing simulation
            if(!tracer[i_det]->postProcessing())
                return false;

            // Does the raytracing simulation considered the scattered light?
            uint ray_result_type = RESULTS_RAY;
            if(dust->getScatteringToRay())
                ray_result_type = RESULTS_FULL;

            // Write results either as text or fits file
            if(!tracer[i_det]->writeDustResults(ray_result_type))
                return false;
            
            // Break if time-dependent
            if(RAY_DT > 0)
            {
                for(uint i_det = start+1; i_det <= stop; i_det++)
                {
                    // Do final steps for the rest of the detectors
                    if(!tracer[i_det]->postProcessing())
                        return false;
                    uint ray_result_type = RESULTS_RAY;
                    if(dust->getScatteringToRay())
                        ray_result_type = RESULTS_FULL;
                    if(!tracer[i_det]->writeDustResults(ray_result_type))
                        return false;
                }
                break;
            }
        }
    }

    // Show that raytracing is finished
    cout << CLR_LINE;
    cout << "- Raytracing dust map           : done" << endl;

    return true;
}

void CRadiativeTransfer::getDustPixelIntensity(CSourceBasic * tmp_source,
                                               double cx,
                                               double cy,
                                               uint i_det,
                                               uint subpixel_lvl,
                                               int i_pix)
{
    bool subpixel = false;
    photon_package * pp;

    subpixel = tracer[i_det]->getUseSubpixel(cx, cy, subpixel_lvl);

    // If any subpixel traveled through other cells, perform subpixelling
    if(subpixel == false)
    {
        // Init variables
        uint nr_used_wavelengths = tracer[i_det]->getNrOfSpectralBins();

        // Create new photon package
        pp = new photon_package();

        // Init photon package
        pp->initMultiStokesVector(nr_used_wavelengths);
        tracer[i_det]->preparePhoton(pp, cx, cy);

        // Calculate continuum emission along one path
        getDustIntensity(pp, tmp_source, cx, cy, i_det, subpixel_lvl, i_pix);

        if(RAY_DT == 0)
            // Add the photon package to the detector
            tracer[i_det]->addToDetector(pp, i_pix);

        // Delete photon package after usage
        delete pp;
    }
    else
    {
        // Repeat this function for each subpixel
        for(int i_sub_x = -1; i_sub_x <= 1; i_sub_x += 2)
        {
            for(int i_sub_y = -1; i_sub_y <= 1; i_sub_y += 2)
            {
                // Calculate positions of each subpixel
                double tmp_cx, tmp_cy;
                tracer[i_det]->getSubPixelCoordinates(subpixel_lvl, cx, cy, i_sub_x, i_sub_y, tmp_cx, tmp_cy);
                // Calculate radiative transfer of the current pixel
                // and add it to the detector at the corresponding position
                getDustPixelIntensity(tmp_source, tmp_cx, tmp_cy, i_det, (subpixel_lvl + 1), i_pix);
            }
        }
    }
}

void CRadiativeTransfer::getDustIntensity(photon_package * pp,
                                          CSourceBasic * tmp_source,
                                          double cx,
                                          double cy,
                                          uint i_det,
                                          uint subpixel_lvl,
                                          int i_pix)
{
    // Set amount of radiation coming from this pixel
    double subpixel_fraction = pow(4.0, -double(subpixel_lvl));

    // Get chosen wavelength parameter
    uint nr_used_wavelengths = tracer[i_det]->getNrSpectralBins();

    // Init multiple Stokes vectors
    MultiStokesVector WMap(nr_used_wavelengths);
    
    // Reset next results counter for time-dependent raytracing
    double t_nextres = RAY_DT;
    double tmp_len = 0.0;
    double len_diff = 0.0;
    uint i_next = 0;
    Vector3D new_pos_xyz;
    Vector3D old_pos_xyz;
    cell_basic * old_cell;
    
    // Update Stokes vectors with emission from background source
    for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
    {
        // Get current wavelength index
        double wID = dust->getWavelengthID(tracer[i_det]->getWavelength(i_wave));

        // Set wavelength index in photon package
        pp->setWavelengthID(wID);

        // Get emission from background source
        WMap.setS(tmp_source->getStokesVector(pp), i_wave);
    };
    

    // Find starting point inside the model and travel through it
    if(grid->findStartingPoint(pp))
    {
        if(RAY_DT > 0)
            pp->updateTotalPathLength(pp->getTmpPathLength());
        
        while((grid->next(pp) && tracer[i_det]->isNotAtCenter(pp, cx, cy)) || len_diff > 0)
        {
            // Get path length through current cell
            double len = pp->getTmpPathLength();
            
            // Set diff length if necessary
            if(len_diff > 0)
            {
                pp->setPosition(old_pos_xyz);
                pp->setPositionCell(old_cell);
                len = len_diff;
            }
            
            // Get necessary quantities from current cell
            double dens_gas = grid->getGasNumberDensity(pp);
            double dens_dust = dust->getNumberDensity(grid, pp);
                
            // Update total pathlength of photon if time-dependent
            if(RAY_DT > 0)
            {
                tmp_len = len/con_c > RAY_DT ? RAY_DT*con_c : len;
                pp->updateTotalPathLength(tmp_len);
                len_diff = len - tmp_len;
                old_pos_xyz = pp->getPosition();
                old_cell = pp->getPositionCell();
                new_pos_xyz = old_pos_xyz - ((len - tmp_len)*pp->getDirection());
            }

            // If the dust density is far too low, skip the current cell
            if(dens_dust >= 1e-200)
            {
                // Init dust matrix
                Matrix2D alpha_dust;
                StokesVector S_dust;
                    
#ifdef CAMPS_BENCHMARK
                // Part to perform Camps et. al (2015) benchmark.
                ofstream myfile;
                string filename = pathOutput + "/stochastic_heating_test.dat";
                myfile.open(filename.c_str());

                double corr_factor = 0;
                if(pathOutput.find("_sil") != std::string::npos)
                    corr_factor = 1.471288e-7;
                else if(pathOutput.find("_gra") != std::string::npos)
                    corr_factor = 1.905816e-7;
                else if(pathOutput.find("_pah") != std::string::npos)
                    corr_factor = 2.227433e-7;

                myfile << "# nr_of wavelength = " << nr_used_wavelengths << ", corr_factor = " << corr_factor
                       << endl;
                for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
                {
                    // Get current wavelength
                    double wavelength = tracer[i_det]->getWavelength(i_wave);

                    // Get current wavelength index
                    double wID = dust->getWavelengthID(wavelength);

                    // Set wavelength index in photon package
                    pp->setWavelengthID(wID);

                    StokesVector tmp_stokes = dust->calcEmissivitiesEmi(grid, pp);
                    myfile << wavelength << TAB << corr_factor * wavelength * tmp_stokes.I() << TAB
                           << corr_factor * wavelength * tmp_stokes.Q() << TAB
                           << corr_factor * wavelength * (tmp_stokes.I() + tmp_stokes.Q()) << endl;
                }
                myfile.close();
                break;
#endif

                for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
                {
                    // Get current wavelength index
                    double wID = dust->getWavelengthID(tracer[i_det]->getWavelength(i_wave));

                    // Set wavelength index in photon package
                    pp->setWavelengthID(wID);

                    // Set stokes vector with thermal emission of the dust grains and
                    // radiation scattered at dust grains
                    S_dust = dust->calcEmissivitiesEmi(
                        grid, pp, stokes_dust_rad_field ? detector_wl_index[i_det] + wID : MAX_UINT);

                    // Get the extinction matrix of the dust grains in the current cell
                    alpha_dust = -1 * dust->calcEmissivitiesExt(grid, pp);

                    // Init a variable to sum up path lengths until cell is crossed
                    double cell_sum = 0.0;

                    // First path length is path through cell
                    double cell_d_l = len;

                    // Save direction of the photon package
                    Vector3D dir_map_xyz = pp->getDirection();

                    // Save the cell entry position of the photon package
                    // (Hint: next(pp) puts the photon onto the border to the next cell)
                    Vector3D pos_xyz_cell = pp->getPosition() - (len * dir_map_xyz);

                    // Set new length if time-dependent
                    if(RAY_DT > 0)
                    {
                        len = tmp_len;
                        cell_d_l = tmp_len;
                    }
                    
                    // Init number of steps
                    ullong kill_counter = 0;

                    // Make sub steps until cell is completely crossed
                    // If the error of a sub step is too high, make the step smaller
                    while(cell_sum < len)
                    {
                        // Increase the kill counter
                        kill_counter++;

                        // If too many sub steps are needed, kill the photon
                        if(kill_counter >= MAX_SOLVER_STEPS)
                        {
                            cout << "\nWARNING: Solver steps > " << MAX_SOLVER_STEPS << ". Too many steps!"
                                 << endl;
                            break;
                        }

                        // Init Runge-Kutta parameters and set it to zero
                        // (see
                        // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method)
                        StokesVector * RK_k = new StokesVector[6];

                        // Calculate result of the radiative transfer equation at each
                        // Runge-Kutta sub position
                        for(uint k = 0; k < 6; k++)
                        {
                            // Init scalar product
                            StokesVector scalar_product;

                            // Calculate multiplication between Runge-Kutta parmeter
                            for(uint i = 0; i < 6; i++)
                                scalar_product += (RK_k[i] * RK_a(i, k));

                            // Calculate new Runge-Kutta parameters as the result of the
                            // radiative transfer equation at the Runge-Kutta sub
                            // positions
                            RK_k[k] = alpha_dust * (scalar_product * cell_d_l + WMap.S(i_wave)) + S_dust;
                        }

                        // Init two temporary Stokes vectors
                        StokesVector stokes_new = WMap.S(i_wave);
                        StokesVector stokes_new2 = WMap.S(i_wave);

                        // Add the result at each Runge-Kutta sub position
                        // to the total Stokes vector (Using Runge-Kutta Fehlberg)
                        for(uint i = 0; i < 6; i++)
                        {
                            stokes_new += (RK_k[i] * cell_d_l) * RK_b1[i];
                            stokes_new2 += (RK_k[i] * cell_d_l) * RK_b2[i];
                        }

                        // Delete the Runge-Kutta pointer
                        delete[] RK_k;

                        // Ignore very small values
                        if(abs(stokes_new.I()) < 1e-200)
                            stokes_new.clearIntensity();
                        if(abs(stokes_new2.I()) < 1e-200)
                            stokes_new2.clearIntensity();

                        // Calculate the difference between the results with two
                        // different precisions to see if smaller steps are needed
                        double epsi, dz_new;
                        calcStepWidth(stokes_new, stokes_new2, cell_d_l, epsi, dz_new);

                        // Is a smaller step width needed
                        if(epsi <= 1.0)
                        {
                            // Add the temporary Stokes vector to the total one
                            WMap.setS(stokes_new, i_wave);

                            // Add optical depth
                            WMap.addT(-alpha_dust(0, 0) * cell_d_l, i_wave);

                            // Add to column density
                            WMap.addSp(dens_gas * cell_d_l, i_wave);

                            // Update the position of the photon package
                            pos_xyz_cell += cell_d_l * dir_map_xyz;

                            // Increase the sum of the cell path lengths
                            cell_sum += cell_d_l;

                            // Find a new path length
                            cell_d_l = min(dz_new, 4 * cell_d_l);

                            // If the new step would exceed the cell, make it smaller
                            if(cell_sum + cell_d_l > len)
                                cell_d_l = len - cell_sum;
                        }
                        else
                        {
                            // Find a smaller path length
                            cell_d_l = max(dz_new, 0.25 * cell_d_l);
                        }
                    }
                }
                
            }
            
            if(RAY_DT > 0)
                // Set photon to adjusted position
                pp->setPosition(new_pos_xyz);

            // Add photon_package to detector if in time
            if(RAY_DT > 0 && (pp->getTotalPathLength()/con_c)/t_nextres >= 1)
            {
                // Calculate detector to add to
                double dt = RAY_DT;
                
                // Check if enough detectors are available
                if(floor(pp->getTotalPathLength()/con_c/dt) > stop)
                {
                    cout << "CRITICAL ERROR: More detectors needed for time dependent emission!" << endl;
                    cout << "Last index of detector " << floor(pp->getTotalPathLength()/con_c/dt) << endl;
                    exit(0);
                }
                else
                {
                    // Calculate next detector
                    i_next = stop - floor(pp->getTotalPathLength()/con_c/dt);
                }
                
                // Init temporary multiple Stokes vectors
                MultiStokesVector WTmp(nr_used_wavelengths);
                
                // Update the multi Stokes vectors for each wavelength
                for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
                {
                    // Copy MultiStokesVector
                    WTmp.setS(WMap.S(i_wave), i_wave);
                    WTmp.setT(WMap.T(i_wave), i_wave);
                    WTmp.setSp(WMap.Sp(i_wave), i_wave);
                    
                    // Get frequency at background grid position
                    double frequency = con_c / tracer[i_det]->getWavelength(i_wave);
                    double mult = 1.0e+26 * subpixel_fraction * tracer[i_det]->getDistanceFactor() * con_c /
                                (frequency * frequency);

                    // Include foreground extinction if necessary
                    mult *= dust->getForegroundExtinction(tracer[i_det]->getWavelength(i_wave));
                    
                    if(WTmp.S(i_wave).I() < 0)
                        WTmp.S(i_wave).setI(0);

                    WTmp.S(i_wave) *= mult;
                    WTmp.setT(WTmp.T(i_wave) * subpixel_fraction, i_wave);
                    WTmp.setSp(WTmp.Sp(i_wave) * subpixel_fraction, i_wave);
                    
                    // Update the photon package with the multi Stokes vectors
                    pp->setMultiStokesVector(WTmp.S(i_wave), i_wave);
                }
                
                // Save old position and StokesVector
                Vector3D old_pos = pp->getPosition();
                
                // Add to detector
                tracer[i_next]->addToDetector(pp, i_pix);
                
                // Set back to old position
                pp->setPosition(old_pos);
                
                // Increase t_nextres
                t_nextres += dt;
            }
        }
        
        // Add to rest of detectors if time-dependent
        if(RAY_DT > 0 && i_next != start)
        {
            // Init temporary multiple Stokes vectors
            MultiStokesVector WTmp(nr_used_wavelengths);
            
            // Update the multi Stokes vectors for each wavelength
            for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
            {
                // Copy MultiStokesVector
                WTmp.setS(WMap.S(i_wave), i_wave);
                WTmp.setT(WMap.T(i_wave), i_wave);
                WTmp.setSp(WMap.Sp(i_wave), i_wave);
                
                // Get frequency at background grid position
                double frequency = con_c / tracer[i_det]->getWavelength(i_wave);
                double mult = 1.0e+26 * subpixel_fraction * tracer[i_det]->getDistanceFactor() * con_c /
                            (frequency * frequency);

                // Include foreground extinction if necessary
                mult *= dust->getForegroundExtinction(tracer[i_det]->getWavelength(i_wave));
                
                if(WTmp.S(i_wave).I() < 0)
                    WTmp.S(i_wave).setI(0);

                WTmp.S(i_wave) *= mult;
                WTmp.setT(WTmp.T(i_wave) * subpixel_fraction, i_wave);
                WTmp.setSp(WTmp.Sp(i_wave) * subpixel_fraction, i_wave);
                
                // Update the photon package with the multi Stokes vectors
                pp->setMultiStokesVector(WTmp.S(i_wave), i_wave);
            }
            
            // Save old position
            Vector3D old_pos = pp->getPosition();
            
            // Loop rest of detectors
            for(uint det = start; det < i_next; det++)
            {
                tracer[det]->addToDetector(pp, i_pix);
                pp->setPosition(old_pos);
                for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
                    pp->setMultiStokesVector(WTmp.S(i_wave), i_wave);
            }
        }
        
    }
    
    // End function here if time-dependent
    if(RAY_DT > 0)
        return;
        
    // Update the multi Stokes vectors for each wavelength
    for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
    {
        // Get frequency at background grid position
        double frequency = con_c / tracer[i_det]->getWavelength(i_wave);
        double mult = 1.0e+26 * subpixel_fraction * tracer[i_det]->getDistanceFactor() * con_c /
                      (frequency * frequency);

        // Include foreground extinction if necessary
        mult *= dust->getForegroundExtinction(tracer[i_det]->getWavelength(i_wave));

        if(WMap.S(i_wave).I() < 0)
            WMap.S(i_wave).setI(0);

        WMap.S(i_wave) *= mult;
        WMap.setT(WMap.T(i_wave) * subpixel_fraction, i_wave);
        WMap.setSp(WMap.Sp(i_wave) * subpixel_fraction, i_wave);

        // Update the photon package with the multi Stokes vectors
        pp->setMultiStokesVector(WMap.S(i_wave), i_wave);
    }
}

void CRadiativeTransfer::calcStellarEmission(uint i_det)
{
    // Init variables for photon positioning on detector
    int i_pix;

    // Transport photon to observer for each detector
    for(uint s = 0; s < sources_mc.size(); s++)
    {
        // Ignore dust as Monte-Carle radiation source
        if(sources_mc[s]->getStringID() == "dust source")
            continue;

        // Get chosen wavelength parameter
        uint nr_used_wavelengths = tracer[i_det]->getNrOfSpectralBins();

        // Init multiple Stokes vectors
        MultiStokesVector WMap(nr_used_wavelengths);

        // Create new photon package
        photon_package * pp = new photon_package();

        // Init photon package
        pp->initMultiStokesVector(nr_used_wavelengths);

        // Get position of source
        Vector3D source_pos = sources_mc[s]->getPosition();

        // Update Stokes vectors with emission from background source
        for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
        {
            // Get frequency at background grid position
            double wavelength = tracer[i_det]->getWavelength(i_wave);
            double frequency = con_c / wavelength;
            double mult = 1.0e+26 * con_c / (frequency * frequency);

            // Get current wavelength index
            double wID = dust->getWavelengthID(wavelength);

            // Set wavelength index in photon package
            pp->setWavelengthID(wID);

            // Set direction of the photon package to the observer
            tracer[i_det]->preparePhotonWithPosition(pp, source_pos, i_pix);

            // Launch photon package
            sources_mc[s]->createDirectRay(pp);

            // Consider the distance to the source
            WMap.setS(pp->getStokesVector(), i_wave);

            // Position the photon inside the grid
            grid->positionPhotonInGrid(pp);

            // Init a variable to save the optical depth
            double tau_obs = 0;

            // Transport photon package through model to obtain the optical depth
            while(grid->next(pp))
            {
                // Get necessary information about the current cell
                double len = pp->getTmpPathLength();
                double dens = dust->getNumberDensity(grid, pp);
                double Cext = dust->getCextMean(grid, pp);

                // Increase the optical depth by the current cell
                tau_obs += Cext * len * dens;
            }
            mult *= exp(-tau_obs) * tracer[i_det]->getDistanceFactor(source_pos);

            // Update the photon package with the multi Stokes vectors
            pp->setMultiStokesVector(WMap.S(i_wave) * mult, i_wave);
        }

        tracer[i_det]->addToDetector(pp, i_pix, true);

        // Delete photon package after usage
        delete pp;
    }
}

// -------------------------------------------
// ------ Calculation of Line transfer -------
// -------------------------------------------

bool CRadiativeTransfer::calcChMapsViaRaytracing(parameters & param)
{
    // Get maximum length of the simulation model
    double max_length = grid->getMaxLength();

    // Create a list for all gas species
    maplist line_ray_detector_list = param.getLineRayDetectors();
    maplist::iterator it;

    // Index to the current detector/tracer
    uint i_det = 0;

    // Perform radiative transfer for each chosen gas species
    for(it = line_ray_detector_list.begin(); it != line_ray_detector_list.end(); ++it)
    {
        // Get ID of the current gas species
        uint i_species = it->first;

        // Get number of spectral line transitions that have to be simulated
        uint nr_of_transitions = param.getNrOfGasSpeciesTransitions(i_species);

        // If no spectral line transitions are chosen for the current gas species, skip
        if(nr_of_transitions == 0)
            continue;

        // Get the detectors of the current gas species
        dlist line_ray_detectors = it->second;

        // Start and Stop define the index of the first and last gas species
        if(i_species < start || i_species > stop)
            continue;

        // Calculate the level populations for each cell
        if(!gas->calcLevelPopulation(grid, i_species))
        {
            cout << "\nERROR: Level population cannot be calculated!";
            return false;
        }

        // Calculate the line broadening for each cell
        gas->calcLineBroadening(grid, i_species);

        // Perform radiative transfer for each chosen spectral line transition
        for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
        {
            uint nr_velocity_channels = tracer[i_det]->getNrSpectralBins();
            if(gas->getZeemanSplitting(i_species) && nr_velocity_channels < 6)
            {
                cout << "\nERROR: The magnetic field information requires at least 6 "
                        "channels\n"
                     << "    for simulations with Zeeman splitting" << endl;
                return false;
            }

            // Select source as specified in the detector
            uint sID = tracer[i_det]->getSourceIndex();
            CSourceBasic * tmp_source;
            tmp_source = sources_ray[sID];

            // Show progress of the current sequence and gas species
            cout << CLR_LINE;
            cout << "-> Channel maps: gas species " << i_species + 1 << " of " << stop + 1 << ", line "
                 << i_line + 1 << " of " << nr_of_transitions << ": 0.0[%]  \r" << flush;

            uint per_max = tracer[i_det]->getNpix();

            // Init counter and percentage to show progress
            ullong per_counter = 0;
            float last_percentage = 0;

            // Calculate pixel intensity for each pixel
#pragma omp parallel for schedule(dynamic)
            for(int i_pix = 0; i_pix < int(per_max); i_pix++)
            {
                double cx = 0, cy = 0;
                if(!tracer[i_det]->getRelPosition(i_pix, cx, cy))
                    continue;

                // Increase counter used to show progress
                per_counter++;

                // Calculate percentage of total progress per source
                float percentage = 100.0 * float(per_counter) / float(per_max);

                // Show only new percentage number if it changed
                if((percentage - last_percentage) > PERCENTAGE_STEP)
                {
#pragma omp critical
                    {
                        cout << "-> Channel maps: gas species " << i_species + 1 << " of " << stop + 1
                             << ", line " << i_line + 1 << " of " << nr_of_transitions << ": " << percentage
                             << " [%]      \r" << flush;
                        last_percentage = percentage;
                    }
                }
                getLinePixelIntensity(tmp_source, cx, cy, i_det, i_species, i_line, uint(0), i_pix);
            }

            // post-process raytracing simulation
            if(!tracer[i_det]->postProcessing())
                return false;

            if(!tracer[i_det]->writeLineResults(gas, i_species, i_line))
                return false;

            // Increment index for line RT
            i_det++;
        }
    }

    cout << CLR_LINE;
    cout << "- Raytracing channel maps       : done" << endl;

    return true;
}

void CRadiativeTransfer::getLinePixelIntensity(CSourceBasic * tmp_source,
                                               double cx,
                                               double cy,
                                               uint i_species,
                                               uint i_line,
                                               uint i_det,
                                               uint subpixel_lvl,
                                               int i_pix)
{
    bool subpixel = false;
    photon_package * pp;

    subpixel = tracer[i_det]->getUseSubpixel(cx, cy, subpixel_lvl);

    // If any subpixel traveled through other cells, perform subpixelling
    if(subpixel == false)
    {
        // Create new photon package
        pp = new photon_package();

        // Init photon package
        uint nr_velocity_channels = tracer[i_det]->getNrSpectralBins();
        pp->initMultiStokesVector(nr_velocity_channels);
        tracer[i_det]->preparePhoton(pp, cx, cy);

        // Calculate line emission along one path
        getLineIntensity(pp, tmp_source, cx, cy, i_det, subpixel_lvl, i_species, i_line);

        tracer[i_det]->addToDetector(pp, i_pix);

        // Delete photon package after usage
        delete pp;
    }
    else
    {
        // Repeat this function for each subpixel
        for(int i_sub_x = -1; i_sub_x <= 1; i_sub_x += 2)
        {
            for(int i_sub_y = -1; i_sub_y <= 1; i_sub_y += 2)
            {
                // Calculate positions of each subpixel
                double tmp_cx, tmp_cy;
                tracer[i_det]->getSubPixelCoordinates(subpixel_lvl, cx, cy, i_sub_x, i_sub_y, tmp_cx, tmp_cy);

                // Calculate radiative transfer of the current pixel
                // and add it to the detector at the corresponding position
                getLinePixelIntensity(
                    tmp_source, tmp_cx, tmp_cy, i_species, i_line, i_det, (subpixel_lvl + 1), i_pix);
            }
        }
    }
}

// "getLineIntensity" is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

void CRadiativeTransfer::getLineIntensity(photon_package * pp,
                                          CSourceBasic * tmp_source,
                                          double cx,
                                          double cy,
                                          uint i_det,
                                          uint subpixel_lvl,
                                          uint i_species,
                                          uint i_line)
{
    // Set amount of radiation coming from this pixel
    double subpixel_fraction = pow(4.0, -double(subpixel_lvl));
    Vector3D obs_vel = tracer[i_det]->getObserverVelocity();

    // Get transition frequency from current species and line
    double frequency = gas->getTransitionFrequencyFromIndex(i_species, i_line);

    // Calculate wavelength and index
    double wavelength = con_c / frequency;
    uint wID = dust->getWavelengthID(wavelength);

    // velocity interpolation
    spline vel_field;
    bool zero_vel_field = true;
    bool grid_has_vel_field = grid->hasVelocityField();

    uint nr_velocity_channels = tracer[i_det]->getNrSpectralBins();
    MultiStokesVector CHMap(nr_velocity_channels);

    // Set photon package to current wavelength
    pp->setWavelengthID(wID);

    for(uint vch = 0; vch < nr_velocity_channels; vch++)
        CHMap.setS(tmp_source->getStokesVector(pp) * con_c / (frequency * frequency), vch);

    if(grid_has_vel_field && gas->getKeplerStarMass() == 0)
    {
        photon_package * tmp_pp = new photon_package();
        tracer[i_det]->preparePhoton(tmp_pp, cx, cy);
        if(grid->findStartingPoint(tmp_pp))
        {
            Vector3D pos_in_grid_0 = tmp_pp->getPosition();
            double spline_x_old = 0;
            while(grid->next(tmp_pp))
            {
                Vector3D dir_map_xyz = tmp_pp->getDirection();
                Vector3D pos_xyz_cell = tmp_pp->getPosition() - (tmp_pp->getTmpPathLength() * dir_map_xyz);
                Vector3D rel_pos = pos_xyz_cell - pos_in_grid_0;

                cell_basic * tmp_cell_pos = tmp_pp->getPositionCell();
                Vector3D cell_center = grid->getCenter(tmp_cell_pos);

                double length_on_line =
                    CMathFunctions::getClosestLinePoint(pos_xyz_cell, tmp_pp->getPosition(), cell_center);
                double spline_x = length_on_line + rel_pos.length();

                if((spline_x_old - spline_x) != 0.0)
                {
                    double spline_y = grid->getVelocityField(tmp_pp) * dir_map_xyz;

                    spline_x_old = spline_x;

                    if(spline_y != 0)
                        zero_vel_field = false;

                    vel_field.setDynValue(spline_x, spline_y);
                }
            }
            vel_field.createDynSpline();
        }
        delete tmp_pp;
    }

    tracer[i_det]->preparePhoton(pp, cx, cy);
    if(grid->findStartingPoint(pp))
    {
        // Save the starting position, to obtain relative positions
        Vector3D pos_in_grid_0 = pp->getPosition();

        // Set photon package to current wavelength
        pp->setWavelengthID(wID);

        // Init matrices
        Matrix2D alpha_ges;

        // Transport the photon package through the model
        while(grid->next(pp) && tracer[i_det]->isNotAtCenter(pp, cx, cy))
        {
            // Get gas temperature from grid
            double temp_gas = grid->getGasTemperature(pp);

            // Perform radiative transfer only if the temperature of the current species
            // are not negligible
            if(temp_gas >= 1e-200)
            {
                double cos_theta = 0, sin_theta = 0, cos_2_phi = 0, sin_2_phi = 0;
                Vector3D mag_field;
                if(gas->getZeemanSplitting(i_species))
                    grid->getMagFieldOrientation(pp, mag_field, cos_theta, sin_theta, cos_2_phi, sin_2_phi);

                // Get the path length through the current cell
                double len = pp->getTmpPathLength();

                // Get necessary quantities from the current cell
                double dens_gas = grid->getGasNumberDensity(pp);
                double dens_species = gas->getNumberDensity(grid, pp, i_species);

                // Calculate the emission of the dust grains
                StokesVector S_dust = dust->calcEmissivitiesHz(grid, pp);

                // Calculate the emission of the gas particles
                StokesVector S_gas;
                if(dens_species > 0)
                {
                    S_gas = gas->calcEmissivities(grid, pp, i_species, i_line) * dens_species;
                    S_gas.setT(S_gas.T() * dens_species);
                }

                // Perform radiative transfer for each velocity channel separately
                for(uint vch = 0; vch < nr_velocity_channels; vch++)
                {
                    // Init variables
                    double velo_dir, cell_sum = 0, cell_d_l = len;
                    ullong kill_counter = 0;

                    // Get direction and entry position of the current cell
                    Vector3D dir_map_xyz = pp->getDirection();
                    Vector3D pos_xyz_cell = pp->getPosition() - (len * dir_map_xyz);
                    Vector3D tmp_pos;

                    // Make sub steps until cell is completely crossed
                    // If the error of a sub step is too high, make the step smaller
                    while(cell_sum < len)
                    {
                        // Increase the kill counter
                        kill_counter++;

                        // If too many sub steps are needed, kill the photon
                        if(kill_counter >= MAX_SOLVER_STEPS)
                        {
                            cout << "\nWARNING: Solver steps > " << MAX_SOLVER_STEPS << ". Too many steps."
                                 << endl;
                            break;
                        }

                        // Init Runge-Kutta parameters and set it to zero
                        // (see
                        // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method)
                        StokesVector * RK_k = new StokesVector[6];

                        // Calculate result of the radiative transfer equation at each
                        // Runge-Kutta sub position
                        for(uint k = 0; k < 6; k++)
                        {
                            // Calculate Runge-Kutta position
                            tmp_pos = pos_xyz_cell + cell_d_l * dir_map_xyz * RK_c[k];

                            // Get the velocity in the photon
                            // direction of the current position
                            if(gas->getKeplerStarMass() > 0)
                            {
                                // Get velocity from Kepler rotation
                                velo_dir =
                                    CMathFunctions::calcKeplerianVelocity(tmp_pos, gas->getKeplerStarMass()) *
                                    dir_map_xyz;
                            }
                            else if(grid_has_vel_field)
                            {
                                // If the velocity field is zero, set velocity to zero
                                // If not, use veloctiy field of the grid
                                if(zero_vel_field)
                                    velo_dir = 0.0;
                                else
                                {
                                    Vector3D rel_pos = tmp_pos - pos_in_grid_0;
                                    velo_dir = vel_field.getValue(rel_pos.length());

                                    // Old version without spline interpolation
                                    // velo_dir = grid->getVelocityField(pp) *
                                    // dir_map_xyz;
                                }
                            }
                            else
                                velo_dir = 0.0;

                            if(obs_vel.length() > 0)
                                velo_dir += obs_vel * dir_map_xyz;

                            const Matrix2D & line_matrix =
                                gas->getLineMatrix(grid,
                                                   pp,
                                                   i_species,
                                                   i_line,
                                                   tracer[i_det]->getVelocityChannel(vch) - velo_dir,
                                                   mag_field,
                                                   cos_theta,
                                                   sin_theta,
                                                   cos_2_phi,
                                                   sin_2_phi);

                            // Combine the Stokes vectors from gas and dust
                            // for emission and extinction
                            StokesVector S_ges = line_matrix * S_gas + S_dust;
                            alpha_ges = S_gas.T() * line_matrix;

                            if(S_dust.T() != 0)
                                for(uint i = 0; i < 4; i++)
                                    alpha_ges(i, i) += S_dust.T();
                            alpha_ges *= -1;

                            // Init scalar product
                            StokesVector scalar_product;

                            // Calculate multiplication between Runge-Kutta parameter
                            for(uint i = 0; i <= k; i++)
                                scalar_product += (RK_k[i] * RK_a(i, k));

                            // Calculate new Runge-Kutta parameters as the result of the
                            // radiative transfer equation at the Runge-Kutta sub
                            // positions
                            RK_k[k] = alpha_ges * (scalar_product * cell_d_l + CHMap.S(vch)) + S_ges;
                        }
                        // Init two temporary Stokes vectors
                        StokesVector stokes_new = CHMap.S(vch);
                        StokesVector stokes_new2 = CHMap.S(vch);

                        for(uint i = 0; i < 6; i++)
                        {
                            stokes_new += (RK_k[i] * cell_d_l) * RK_b1[i];
                            stokes_new2 += (RK_k[i] * cell_d_l) * RK_b2[i];
                        }

                        // Delete the Runge-Kutta pointer
                        delete[] RK_k;

                        // Ignore very small values
                        if(abs(stokes_new.I()) < 1e-200)
                            stokes_new.clearIntensity();
                        if(abs(stokes_new2.I()) < 1e-200)
                            stokes_new2.clearIntensity();

                        // Calculate the difference between the results with two
                        // different precisions to see if smaller steps are needed
                        double epsi, dz_new;
                        calcStepWidth(stokes_new, stokes_new2, cell_d_l, epsi, dz_new);

                        // Is a smaller step width needed
                        if(epsi <= 1.0)
                        {
                            // Backup old intensity
                            double old_stokes = CHMap.I(vch);

                            // Stokes_new is the current flux of this line-of-sight
                            CHMap.setS(stokes_new, vch);

                            // Columns density
                            double column_density = dens_gas * cell_d_l;

                            if(gas->getZeemanSplitting(i_species))
                            {
                                // Total increase of the intensity along the line-of-sight
                                double column_flux = (stokes_new.I() - old_stokes);

                                // Total magnetic field strength of the current cell
                                double mag_strength = mag_field.length();

                                // LOS magnetic field strength of the current cell
                                double los_mag_strength = (dir_map_xyz * mag_field);

                                // Magnetic field strength in the line-of-sight direction
                                // weighted with the intensity increase of the current
                                // cell
                                double column_int_mag_field_los = los_mag_strength * column_flux;

                                // Total magnetic field strength weighted with the
                                // intensity increase of the current cell
                                double column_int_mag_field = mag_strength * column_flux;

                                // Intensity weighted LOS magnetic field
                                CHMap.addSp(column_int_mag_field_los, 0);

                                // Intensity weighted total magnetic field
                                CHMap.addSp(column_int_mag_field, 1);

                                // Flux component for weighting
                                CHMap.addSp(column_flux, 2);

                                if(vch == 0)
                                {
                                    // Magnetic field strength in the line-of-sight
                                    // direction weighted with the gas density of the
                                    // current cell
                                    double column_dens_mag_field_los = los_mag_strength * column_density;

                                    // Total magnetic field strength weighted with the
                                    // gas density of the current cell
                                    double column_dens_mag_field = mag_strength * column_density;

                                    // Density weighted LOS magnetic field
                                    CHMap.addSp(column_dens_mag_field_los, 3);

                                    // Density weighted magnetic field
                                    CHMap.addSp(column_dens_mag_field, 4);

                                    // Column density of the total gas
                                    CHMap.addSp(column_density, 5);
                                }
                            }
                            else if(vch == 0)
                            {
                                // Column density of the total gas
                                CHMap.addSp(column_density, 0);
                            }

                            // Save the optical depth of each velocity channel, if
                            // magnetic field analysis is not chosen
                            CHMap.addT(-alpha_ges(0, 0) * cell_d_l, vch);

                            // Update the position of the photon package
                            pos_xyz_cell += cell_d_l * dir_map_xyz;

                            // Increase the sum of the cell path lengths
                            cell_sum += cell_d_l;

                            // Find a new path length
                            cell_d_l = min(dz_new, 4 * cell_d_l);

                            // If the new step would exceed the cell, make it smaller
                            if(cell_sum + cell_d_l > len)
                                cell_d_l = len - cell_sum;
                        }
                        else
                        {
                            // Find a smaller path length
                            cell_d_l = max(dz_new, 0.25 * cell_d_l);
                        }
                    }
                }
            }
        }
    }
    // Update the multi Stokes vectors for each velocity channel
    for(uint vch = 0; vch < nr_velocity_channels; vch++)
    {
        double mult = 1.0e+26 * subpixel_fraction * tracer[i_det]->getDistanceFactor();

        // Include foreground extinction if necessary
        mult *= dust->getForegroundExtinction(tracer[i_det]->getWavelength(wavelength));

        if(CHMap.S(vch).I() < 0)
            CHMap.S(vch).setI(0);

        CHMap.S(vch) *= mult;
        CHMap.setT(CHMap.T(vch) * subpixel_fraction, vch);
        CHMap.setSp(CHMap.Sp(vch) * subpixel_fraction, vch);

        // Update the photon package with the multi Stokes vectors
        pp->setMultiStokesVector(CHMap.S(vch), vch);
    }
}

// ---------------------------------------------------------------
// ------ Calculation of time-dependent radiation transfer -------
// ---------------------------------------------------------------


// Full time-dependent transfer with temperature calculations

bool CRadiativeTransfer::calcMonteCarloTimeTransfer(uint command,
                                                    parameters & param,
                                                    bool use_energy_density,
                                                    bool disable_reemission)
{
    
    // Init variables
    ullong nr_of_photons;
    ullong ph_max, nr_of_wavelength;
    
    ullong per_counter = 0;
    float last_percentage = 0;
    
    ullong kill_counter = 0;
    uint max_source = uint(sources_mc.size());
    double dt, tend;
    dt = param.getTimeStep();
    tend = param.getTotalTime();
    
    // Init arrays for emission, absorption and inner energy
    ulong max_cells = grid->getMaxDataCells();
    dlist dust_abs(max_cells, 0.0);
    dlist dust_em(max_cells);
    dlist dust_u(max_cells);
    
    // Init photon pointer stack
    vector<photon_basic*> pp_stack;

    // Init enthalpies
    if (!dust->initMeanEnthalpy())
        return false;
    
    // Set or load initial temperature dist  
    if(dust->getDustTemperature(grid,grid->getCellFromIndex(1)) > 0)
    {
        // Start initial dust photons in every cell according to given temp dist
        if(!dust->initLamCdf())
            return false;
        
        // Set enthalpy from temperature
        for(long c_i = 0; c_i < long(max_cells); c_i++)
        {
            cell_basic * cell = grid->getCellFromIndex(c_i);
            dust_u[c_i] = dust->getEnthalpyMean(grid,cell,dust->getDustTemperature(grid,cell))*grid->getVolume(cell);
        }
        
        //if (!startInitialDustPhotons(dt, dust_em, pp_stack))
        //    return false;
    }
    else
    {
        // Set minimum enthalpy
        for(long c_i = 0; c_i < long(max_cells); c_i++)
        {
            cell_basic * cell = grid->getCellFromIndex(c_i);
            dust_u[c_i] = dust->getEnthalpyMean(grid,cell,uint(128))*grid->getVolume(cell);
        }
    
        // Calc initial temperature dist
        if (!setTemperatureFromU(dust_u))
            return false;
    }
    
    // Number of photons per timestep (nr of photons first source in cmd-file)
    llong nr_of_photons_step = sources_mc[0]->getNrOfPhotons();
    
    // Set points in time to print out results
    double t_results = param.getTimeOut();
    
    // Progress output
    cout << CLR_LINE;
    cout << "Time step in seconds: " << dt << endl;
    cout << "Number of photons per time step: " << nr_of_photons_step << endl;
    cout << "Length of simulation in seconds: " << tend << endl;
    cout << "Output time step in seconds: " << t_results << endl;
    cout << "-> Calculation of time-dependent transfer : 0.0[%]                        \r";
    
    pp_stack.reserve(nr_of_photons_step*(tend/dt));
    
    // Init photon deletion marker stack
    uilist pp_del;
    pp_del.reserve(nr_of_photons_step);
    
    // Set double for next output
    double t_nextres = t_results;
    
    // Loop over time series
    for(double t = 0; t < tend; t+=dt)
    {
        // Number of photons emittet from dust or source
        ullong N_d = 0;
        ullong N_s = 0;
        
        // Calc emission rates
#pragma omp parallel for schedule(dynamic)
        for(long c_i = 0; c_i < long(max_cells); c_i++)
            dust_em[c_i] = dust->getCellEmissionRate(grid,grid->getCellFromIndex(c_i));
        
        // Calc dust luminosity
        double L_d = 0;
        for(long c_i = 0; c_i < long(max_cells); c_i++)
        {
            cell_basic * cell = grid->getCellFromIndex(c_i);
            L_d += dust_em[c_i] * dust->getNumberDensity(grid, cell);
        }
        
        // Init (star)source and dust (tbd: multiple sources)
        CSourceBasic * source = sources_mc[0];
        CSourceBasic * dust_source = sources_mc[1];
        
        // Init star source and dust wave cdf
        if (t == 0)
        {
            if(!source->initSource(0, 2, false))
                return false;
            if(!dust->initLamCdf())
                return false;
        }
        
        // Calc emission probability of source and dust
        double L_s = source->getLuminosity();
        //if (t >= 50000 && t <= 60000)
        //    L_s *= 4;
        double p_d = L_d/(L_s + L_d);
        
        // Calc cumulative probability dist for cell emission
        dlist p_i(max_cells);
        double temp = 0.0;

        for(long c_i = 0; c_i < long(max_cells); c_i++)
        {
            cell_basic * cell = grid->getCellFromIndex(c_i);
            temp += dust_em[c_i] * dust->getNumberDensity(grid, cell);
            p_i[c_i] = temp/L_d;
        }
        
        // Define last index of stored photon packages
        ullong last = (pp_stack.size() > 0) ? pp_stack.size() : 0;
        
        // Start photons from dust and/or source and store in photon stack
        for (llong i = 0; i < nr_of_photons_step; i++)
        {
            pp_stack.push_back(new photon_basic());
            // Generation marker could be set here
            // pp_stack[i]->setTime(t); tbd
            
            // Decide wether photons are emitted from source or dust and emitt
            pp_stack[last+i]->initRandomGenerator(last+i);
            double rnd = pp_stack[last+i]->getRND();
            
            if(rnd < p_d)
            {
                // Emission from dust source
                
                // Get random number for cell index
                rnd = pp_stack[last+i]->getRND();
                
                // Find cell index from probability dist for cell emission
                ulong i_cell = CMathFunctions::biListIndexSearch(rnd, p_i);
                
                // Set min or max
                if (rnd < p_i[0])
                    i_cell = 0;
                if (rnd > p_i[max_cells-1])
                    i_cell = max_cells-1;
                
                // Get current cell
                cell_basic * cell = grid->getCellFromIndex(i_cell);
                
                // Set photon package into cell
                pp_stack[last+i]->setPositionCell(cell);
                
                // Set random direction, position and coordinate system for scattering
                dust_source->createNextRay(pp_stack[last+i], ullong(i));
                N_d++;
            }
            else
            {
                // Emission from source (star)
                source->createNextRay(pp_stack[last+i], ullong(i));
                N_s++;
            }
        }
        
        // Correct photon energy and set wavelength for dust photons
#pragma omp parallel for schedule(dynamic)
        for (llong i = 0; i < nr_of_photons_step; i++)
        {
            // Check if dust emission
            if(pp_stack[last+i]->getWavelengthID() == MAX_UINT)
            {
                // Get temperature ID and find wavelength
                cell_basic * cell = pp_stack[last+i]->getPositionCell();
                uint wID = dust->getWavelengthIDfromCdf(grid, cell, pp_stack[last+i]);
                
                // Set wavelength
                pp_stack[last+i]->setWavelengthID(wID);
                
                // Set Stokes vector of photon package
                double energy = L_d*dt/N_d;
                pp_stack[last+i]->setStokesVector(StokesVector(energy, 0, 0, 0));
            }
            else
            {
                // Correct Stokes vector of (stellar) source emission
                double energy = L_s*dt/N_s;
                // Scale for theta bias (see Source.cpp)
                uint exponentThetaBias = 3;
                if(exponentThetaBias != 1)
                {
                    double theta = pp_stack[last+i]->getDirection().getSphericalCoord().Theta();
                    energy *= exponentThetaBias * pow(abs(cos(theta)),(double) 1.-1./exponentThetaBias);
                }
                // Assign Stokes vector
                pp_stack[last+i]->setStokesVector(StokesVector(energy, 0, 0, 0));
            }
            // Get tau for first interaction
            pp_stack[last+i]->setTmpPathLength(-log(1.0 - pp_stack[i]->getRND()));
        }
        
        // Perform radiative transfer of all photons in photon stack
#pragma omp parallel for schedule(guided)
        for (llong i = 0; i<pp_stack.size(); i++)
        {
            // Init variables
            double end_tau, Cext, Csca;
            Vector3D old_pos;
    
            // Get tau of last time step
            end_tau = pp_stack[i]->getTmpPathLength();
            
            // Save the old position to use it again
            old_pos = pp_stack[i]->getPosition();
            
            // Init current number of interactions
            ullong interactions = 0;
            
            // Init variables for tau, path length, density
            double tmp_tau, len, dens;
            
            // Set dl = 0 and t = current_time (tbd)
            // Init variable for overall light travel distance
            double dl = 0;
            
            // Transfer photon through grid while still in time
            while(dl < con_c*dt)
            {
                // If end of grid reached, end photon transfer
                if(!grid->next(pp_stack[i]))
                    break;
                
                // If max interactions is reached or no energy left, end photon transfer
                if(interactions >= MAX_INTERACTION || pp_stack[i]->getStokesVector().I() < 1e-200)
                {
#pragma omp critical
                    pp_del.push_back(uint(i));
                    break;
                }
                
                // Get path length through current cell
                len = pp_stack[i]->getTmpPathLength();
                
                // Check rest of light travel time and adjust postion
                if (dl + len > con_c*dt)
                {
                    len = con_c*dt - dl;
                    pp_stack[i]->adjustPosition(old_pos, len);
                }
                
                // Get dust number density of the current cell
                dens = dust->getNumberDensity(grid, pp_stack[i]);
                
                // Calculate the dust absorption cross section (for random
                // alignment)
                Cext = dust->getCextMean(grid, pp_stack[i]);
                
                // Calculate optical depth that is obtained on the path_length
                tmp_tau = Cext * len * dens;
                
                // optical depth for interaction is reached or not
                if(tmp_tau > end_tau)
                {
                    // Increase interaction counter
                    interactions++;
                    
                    // Reduce the photon position to match the exact
                    // interaction position
                    len = len * end_tau / tmp_tau;
                    
                    pp_stack[i]->adjustPosition(old_pos, len);
                    
                    // Update data in grid like spectral length or radiation field
                    updateRadiationField(pp_stack[i], tend);
                    
                    // Update photon Stokes vector (only I, therefore pol results wrong, tbd)
                    double energy = len * pp_stack[i]->getStokesVector().I() * 
                                    dust->getCabsMean(grid, pp_stack[i]); 
                                    
                    // Check if energy left
                    //if (energy*dens > pp_stack[i]->getStokesVector().I())
                    //    energy = pp_stack[i]->getStokesVector().I()/dens;
                    
                    // Add to absorption estimate
                    uint c_id = (pp_stack[i]->getPositionCell())->getID();
#pragma omp atomic update
                    dust_abs[c_id] += energy/dt;
                    
                    // Calculate the dust scattering cross section (for random
                    // alignment)
                    Csca = dust->getCscaMean(grid, pp_stack[i]);
                    
                    // Calculate albedo and check if absorption or
                    // scattering occurs
                    double albedo = Csca / Cext;
                    
                    if(pp_stack[i]->getRND() < albedo)
                    {
                        
                        // Perform simple photon scattering without
                        // changing the Stokes vectors
                        dust->scatter(grid, pp_stack[i], true);
                        
                        // Subtract energy from photon package
                        //pp_stack[i]->addStokesVector(StokesVector(-(energy*dens), 0, 0, 0));
                    }
                    else
                    {
                        // Remove photon
#pragma omp critical
                        pp_del.push_back(uint(i));
                        break;
                    }
                    
                    // Add up light travel distance
                    dl += len;
                    
                    // Calculate new optical depth for next interaction
                    end_tau = -log(1.0 - pp_stack[i]->getRND());
                }
                else
                {
                    // Update data in grid like spectral length or radiation field or energy array
                    updateRadiationField(pp_stack[i], tend);
                    
                    if(dens > 1e-200)
                    {
                        // Update photon Stokes vector (only I, therefore pol results wrong, tbd)
                        double energy = len * pp_stack[i]->getStokesVector().I() * 
                                        dust->getCabsMean(grid, pp_stack[i]); 
                                                                
                        // Check if energy left
                        //if (energy*dens > pp_stack[i]->getStokesVector().I())
                        //    energy = pp_stack[i]->getStokesVector().I()/dens;
                            
                        // Subtract from photon package
                        //pp_stack[i]->addStokesVector(StokesVector(-(energy*dens), 0, 0, 0));
                        
                        // Add to absorption estimate
                        uint c_id = (pp_stack[i]->getPositionCell())->getID();
#pragma omp atomic update
                        dust_abs[c_id] += energy/dt;
                    }
                    
                    // Add up light travel time
                    dl += len;

                    // Remove the traveled distance from optical depth
                    end_tau -= tmp_tau;
                }
                // Save photon position to adjust it if necessary
                old_pos = pp_stack[i]->getPosition();
            }
            
            // Save rest of tau for next timestep in TmpPathLength
            pp_stack[i]->setTmpPathLength(end_tau);
            
            // Remove photon if outside grid
            if(!grid->positionPhotonInGrid(pp_stack[i]) && pp_stack[i]->getStokesVector().I() > 1e-200 && interactions <= MAX_INTERACTION)
            {
#pragma omp critical
                pp_del.push_back(uint(i));  
            }
        
        }
        
//         
//         // Erase photons from stack that left the grid or are absorbed
//         for (llong i = pp_del.size()-1; i>0; i--)
//         {
//                 kill_counter++;
//                 delete pp_stack[pp_del[i]];
//                 pp_stack.erase(pp_stack.begin()+pp_del[i]);
//         }

        // Erase photons from stack (remove_if)
        uint del_size = pp_del.size();
        uint pp_size = pp_stack.size();
        
        if(del_size > 0)
        {
            // Sort deletion marker in descending order with reverse iterators
            sort(pp_del.rbegin(), pp_del.rend());
            
            for (llong i = 0; i < del_size-1; i++)
            {
                kill_counter++;
                delete pp_stack[pp_del[i]];
                swap(pp_stack[pp_del[i]], pp_stack[pp_size-(i+1)]);
            }
            pp_stack.erase(pp_stack.end()-del_size,pp_stack.end());
            
            // Reset photon deletion marker
            pp_del.clear();
        }
        
        // Calc new inner energy for all cells e = e + (A-E)*dt
#pragma omp parallel for schedule(dynamic)
        for(long c_i = 0; c_i < long(max_cells); c_i++)
        {
            // Check for sublimation radius
            cell_basic * cell = grid->getCellFromIndex(c_i);
            if(dust->getNumberDensity(grid,cell) > 1e-200)
                dust_u[c_i] = dust_u[c_i] + (dust_abs[c_i] - dust_em[c_i]) * dt;
            
            // Reset absorption estimator
            dust_abs[c_i] = 0.0;
            
            // Set minimum enthalpy
            double u_min = dust->getEnthalpyMean(grid,cell,uint(128))*grid->getVolume(cell);
            if(dust_u[c_i] < u_min)
                dust_u[c_i] = u_min;
        }
        
        // Calc new temperatures
        if (!setTemperatureFromU(dust_u))
            return false;
        
        // Write out temperature dist
        if(t > 0 && (t/t_nextres >= 1))
        {
            ostringstream s;
            s << t;
            grid->saveBinaryGridFile(param.getPathOutput() + "grids/grid_temp_"+s.str()+".dat");
            t_nextres += t_results;
        }
        
        // Move to next timestept and or set new luminosities
        
        per_counter++;
        float percentage = 100.0 * float(per_counter) / float(tend/dt);
        
        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
                cout << "-> Calculation of time-dependent transfer : "
                     << 100.0 * float(per_counter) / float(tend/dt) << " [%]              \r" << flush;

                last_percentage = percentage;
        }
    }
    
    cout << CLR_LINE;
    cout << "- Calculation of time-dependent transfer : done" << endl;
    
    // Show amount of killed photons
    cout << "- Photons killed                    : " << kill_counter << endl;
    
    // Deallocate photons in stack
    for (llong i = 0; i<pp_stack.size(); i++)
        delete pp_stack[i];
    
    return true;
}

bool CRadiativeTransfer::setTemperatureFromU(dlist dust_u)
{
    // Get Temperature from the the dust mix
    // Use integraded enthalpy of dust mix
    // Tbd: add more than one dust mix
    
    ulong max_cells = grid->getMaxDataCells();
    
#pragma omp parallel for schedule(dynamic)
    for(long c_i = 0; c_i < long(max_cells); c_i++)
    {
        cell_basic * cell = grid->getCellFromIndex(c_i);
        uint i_mixture = dust->getMixtureID(grid, cell);
        
        // Get mean enthalpy for dust mixture in given cell
        dlist Enthalpy(dust->getNrOfCalorimetryTemperatures(grid,cell));
        for(uint T = 0; T < dust->getNrOfCalorimetryTemperatures(grid,cell); T++)
            Enthalpy[T] = dust->getEnthalpyMean(grid,cell,T);
            
        // Convert dust_u to enthalpy per grain
        double hd2 = (dust->getNumberDensity(grid,cell) > 1e-200) ? (dust_u[c_i]/grid->getVolume(cell)) : 0;
        
        // Check if hd2 exceeds max or min
        if (hd2 < Enthalpy[0])
            hd2 = Enthalpy[0];
        if (hd2 > Enthalpy[Enthalpy.size()-1])
            hd2 = Enthalpy[Enthalpy.size()-1];
        
        // Find the corresponding temperature
        uint t_i = CMathFunctions::biListIndexSearch(hd2, Enthalpy);
        double temp = dust->getCalorimetricTemperature(grid,cell,t_i);
        
        // Set min temperature
        if(temp < 2.72503)
            temp = 2.72503;
        
        // Set termperature in cell
        grid->setDustTemperature(cell, temp);
    }
    
    return true;
}

bool CRadiativeTransfer::startInitialDustPhotons(double dt, dlist dust_em, vector<photon_basic*> &pp_stack)
{
    // Starts dust photons from every cell to get initial conditions

    ulong max_cells = grid->getMaxDataCells();
    ullong nr_dust_photons = 100;
    ullong nr_source_photons = 1000;
    double L_d = 0;
    
    dt = (100 * con_AU) / con_c;
    
    // Write N_z - 2 to take empty cells out
    uint N_z = grid->getNumberZ()-2;
    
    // Index of last outer cell
    ulong last_cell = (grid->getDataID() == GRID_ID_CYL) ? (max_cells-grid->getNumberZ()-1) : (max_cells-2);
    
    CSourceBasic * source = sources_mc[0];
    CSourceBasic * dust_source = sources_mc[1];
    
    if(!source->initSource(0, 2, false))
        return false;
    
    ullong per_counter = 0;
    float last_percentage = 0;
    
    float total_size = nr_source_photons*max_cells+nr_dust_photons*(max_cells-(max_cells-last_cell));
    pp_stack.reserve(uint(total_size));
    
    // Progress output
    cout << CLR_LINE;
    cout << "-> Start initial photons from temp dist : 0.0[%]                        \r";
    
    // Set last position in photon stack
    ullong last = 0;
    
    // Set photons in inner cells (r<Rin)
    for(long c_i = long(last_cell); c_i < long(max_cells)-1; c_i++)
    {
        cell_basic * cell = grid->getCellFromIndex(c_i);
       
        // Place source photons in center cell
        for (llong i = 0; i < nr_source_photons; i++)
        {
            pp_stack.push_back(new photon_basic());
            
            // Set random direction, position and coordinate system for scattering
            source->createNextRay(pp_stack[last+i], i);
            
            // Set photon package into cell
            pp_stack[last+i]->setPositionCell(cell);
            // Set random position in cell
            grid->setRndPositionInCell(pp_stack[last+i]);
            
            // Adjust direction
            Vector3D dir = (pp_stack[last+i]->getPosition()-source->getPosition()).normalized();
            pp_stack[last+i]->setDirection(dir);
            
            // Correct Stokes vector of (stellar) source emission
            double energy = ((source->getLuminosity())*dt/(nr_source_photons*max_cells))*grid->getSolidAngle(cell);
            pp_stack[last+i]->setStokesVector(StokesVector(energy, 0, 0, 0));
            
            // Get tau for first interaction
            pp_stack[last+i]->setTmpPathLength(-log(1.0 - pp_stack[i]->getRND()));
            
            // Increase progress counter
            per_counter++;
            float percentage = 100.0 * float(per_counter) / total_size;
        
            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
                    cout << "-> Start initial photons from temp dist : "
                        << percentage << " [%]              \r" << flush;
                    last_percentage = percentage;
            }
        }
        
        // Update last position in stack
        last += nr_source_photons;
    }
    
    // Set photons in rest of grid
    for(long c_i = 0; c_i < long(last_cell); c_i++)
    {
        cell_basic * cell = grid->getCellFromIndex(c_i);
        
        // Get radiation field of the star
        
        // Init variables for optical depth calculation
        double len, dens, Cext, tau = 0;
        
        // Start temporary photon from source in direction of cell 
        photon_basic * pp_depth = new photon_basic;
        source->createNextRay(pp_depth, c_i);
        
        // Adjust direction
        Vector3D dir = (grid->getCenter(cell)-pp_depth->getPosition()).normalized();
        pp_depth->setDirection(dir);
        
        // Transport through grid to get optical depth
        while (grid->next(pp_depth) && ((pp_depth->getPositionCell())->getID() != c_i))
        {
            // Get the traveled distance
            len = pp_depth->getTmpPathLength();

            // Get the current density
            dens = dust->getNumberDensity(grid, pp_depth);

            // Get the current mean extinction cross-section
            Cext = dust->getCextMean(grid, pp_depth);

            // Add the optical depth of the current path to the total optical depth
            tau += Cext * len * dens;
        }
        
        // Delete temporary photon package
        delete pp_depth;
        
        // Start source photons in every cell
        for (llong i = 0; i < nr_source_photons; i++)
        {
            pp_stack.push_back(new photon_basic());
            
            // Set random direction, position and coordinate system for scattering
            source->createNextRay(pp_stack[last+i], i);
            
            // Set photon package into cell
            pp_stack[last+i]->setPositionCell(cell);
            // Set random position in cell
            grid->setRndPositionInCell(pp_stack[last+i]);
            
            // Adjust direction
            Vector3D dir = (pp_stack[last+i]->getPosition()-source->getPosition()).normalized();
            pp_stack[last+i]->setDirection(dir);
            
            // Correct Stokes vector of (stellar) source emission
            double energy = ((source->getLuminosity())*dt/(nr_source_photons*max_cells))*grid->getSolidAngle(cell);
            pp_stack[last+i]->setStokesVector(StokesVector(energy, 0, 0, 0));
            
            // Reduce Stokes vector by the optical depth
            pp_stack[last+i]->getStokesVector() *= exp(-tau);
            
            // Get tau for first interaction
            pp_stack[last+i]->setTmpPathLength(-log(1.0 - pp_stack[i]->getRND()));
            
            // Increase progress counter
            per_counter++;
            float percentage = 100.0 * float(per_counter) / total_size;
        
            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
                    cout << "-> Start initial photons from temp dist : "
                        << percentage << " [%]              \r" << flush;
                    last_percentage = percentage;
            }
        }
        
        // Update last position in stack
        last += nr_source_photons;
        
        // Get emission rate
        dust_em[c_i] = dust->getCellEmissionRate(grid,cell);
        
        // Start photons from every cell
        L_d = dust_em[c_i] * dust->getNumberDensity(grid, cell);
        
        for (llong i = 0; i < nr_dust_photons; i++)
        {
            pp_stack.push_back(new photon_basic());
            
            // Set photon package into cell
            pp_stack[last+i]->setPositionCell(cell);
            
            // Set random direction, position and coordinate system for scattering
            dust_source->createNextRay(pp_stack[last+i], ullong(i));
            
            // Get temperature ID and find wavelength
            uint wID = dust->getWavelengthIDfromCdf(grid, cell, pp_stack[last+i]);
            
            // Set wavelength
            pp_stack[last+i]->setWavelengthID(wID);
            
            // Set Stokes vector of photon package
            // double energy = L_d*dt/nr_dust_photons;
            double energy = (L_d*dt)/nr_dust_photons;
            pp_stack[last+i]->setStokesVector(StokesVector(energy, 0, 0, 0));
            
            // Get tau for first interaction
            pp_stack[last+i]->setTmpPathLength(-log(1.0 - pp_stack[i]->getRND()));
            
            // Increase progress counter
            per_counter++;
            float percentage = 100.0 * float(per_counter) / total_size;
        
            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
                    cout << "-> Start initial photons from temp dist : "
                        << percentage << " [%]              \r" << flush;
                    last_percentage = percentage;
            }
        }
        
        // Update last position in stack
        last += nr_dust_photons;
    }
    
    cout << CLR_LINE;
    cout << "- Start initial photons from temp dist : done" << endl;
    
    return true;
}
