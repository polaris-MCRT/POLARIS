#include "RadiativeTransfer.h"
#include "CommandParser.h"
#include "Detector.h"
#include "Stokes.h"
#include "GasSpecies.h"
#include "MathFunctions.h"
#include "OPIATE.h"

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
                break;

            case DET_SPHER:
                tracer[i_det] = new CRaytracingHealPix(grid);
                // Disable Stokes rad field for HealPix detector
                stokes_dust_rad_field = false;
                break;

            case DET_POLAR:
                tracer[i_det] = new CRaytracingPolar(grid);
                break;

            case DET_SLICE:
                tracer[i_det] = new CRaytracingSlice(grid);
                break;
        }

        if(!tracer[i_det]->setDustDetector(pos, param, dust_ray_detectors, max_length, pathOutput))
            return false;
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
                break;

            case DET_SPHER:
                tracer[i_det] = new CRaytracingHealPix(grid);
                break;

            case DET_POLAR:
                tracer[i_det] = new CRaytracingPolar(grid);
                break;

            case DET_SLICE:
                tracer[i_det] = new CRaytracingSlice(grid);
                break;

            default:
                cout << "ERROR: Wrong detector ID" << endl;
                return false;
                break;
        }

        if(!tracer[i_det]->setSyncDetector(pos, param, sync_ray_detectors, max_length, pathOutput))
            return false;
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

    uint nr_gas_species = param.getNrOfGasSpecies();

    // Get start and stop id of detectors
    start = param.getStart();
    stop = param.getStop();

    // Get maximum length of the simulation model
    double max_length = grid->getMaxLength();

    // Init array of tracer base class pointer
    nr_ray_detectors = gas->getTotalNrOfSpectralLines();
    tracer = new CRaytracingBasic *[nr_ray_detectors];

    if(nr_gas_species == 0)
    {
        cout << "\nERROR: No gas species transition defined!" << endl;
        return false;
    }

    uint i_det = 0;
    for(uint i_species = 0; i_species < nr_gas_species; i_species++)
    {
        dlist line_ray_detectors = param.getLineRayDetector(i_species);

        for(uint i_line = 0; i_line < param.getNrOfSpectralLines(i_species); i_line++)
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
                    break;

                case DET_SPHER:
                    tracer[i_det] = new CRaytracingHealPix(grid);
                    break;

                case DET_POLAR:
                    tracer[i_det] = new CRaytracingPolar(grid);
                    break;

                case DET_SLICE:
                    tracer[i_det] = new CRaytracingSlice(grid);
                    break;
            }
            if(!tracer[i_det]->setLineDetector(pos, param, line_ray_detectors, pathOutput, max_length))
                return false;

            // Increment index for line RT since its distributed over species and transitions
            i_det++;
        }
    }

    initiateRungeKuttaFehlberg();

    return true;
}

bool CRadiativeTransfer::initiateOPIATERaytrace(parameters & param)
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

    if(op == 0)
    {
        cout << "\nERROR: No OPIATE database loaded!" << endl;
        return false;
    }

    // Get start and stop id of detectors
    start = param.getStart();
    stop = param.getStop();

    // Get maximum length of the simulation model
    double max_length = grid->getMaxLength();

    // Init array of tracer base class pointer
    nr_ray_detectors = param.getNrOfOPIATESpecies();
    tracer = new CRaytracingBasic *[nr_ray_detectors];

    dlist op_ray_detectors = param.getOPIATERayDetectors();

    for(uint i_det = 0; i_det < nr_ray_detectors; i_det++)
    {
        // Calculate the starting position of each detector for a gas species
        uint pos = NR_OF_OPIATE_DET * i_det;

        uint nr_source = uint(op_ray_detectors[pos]);

        // Get the ID of the chosen detector
        uint detector_id = uint(op_ray_detectors[pos + NR_OF_OPIATE_DET - 4]);

        if(nr_source > sources_ray.size())
        {
            cout << "\nERROR: ID of source (" << nr_source << ") larger than max. amount ("
                 << sources_ray.size() << ") of defined sources!" << endl;
            return false;
        }

        // Create detector for current simulation
        switch(detector_id)
        {
            case DET_PLANE:
                tracer[i_det] = new CRaytracingCartesian(grid);
                break;

            case DET_SPHER:
                tracer[i_det] = new CRaytracingHealPix(grid);
                break;

            case DET_SLICE:
                //tracer[i_det] = new CRaytracingSlice(grid);
                cout << "ERROR: Slice detector not yet fully implemented!" << endl;
                break;

            default:
                //tracer[i_det] = new CRaytracingSlice(grid);
                cout << "ERROR: Detector not yet fully implemented!" << endl;
                break;
        }

        if(!tracer[i_det]->setOPIATEDetector(pos, param, op_ray_detectors, pathOutput, max_length))
            return false;
    }

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

        #if (USE_PRECALC_TABLE)
            int nr_of_dust_wavelengths = dust->getNrOfWavelength();
            ulong nr_of_cells = grid->getMaxDataCells();

            grid->initPreCalcTables(nr_of_dust_wavelengths);

#pragma omp parallel for schedule(dynamic) collapse(2)
            for(int wID = 0; wID < nr_of_dust_wavelengths; wID++)
            {
                for(long i_cell = 0; i_cell < long(nr_of_cells); i_cell++)
                {
                    photon_package pp = photon_package();

                    pp.setWavelength(dust->getWavelength(wID), wID);
                    // Put photon package into current cell
                    pp.setPositionCell(grid->getCellFromIndex(i_cell));

                    double i_cext = dust->getCextMean(grid, pp);
                    double i_cabs = dust->getCabsMean(grid, pp);
                    double i_csca = dust->getCscaMean(grid, pp);
                    grid->setCextMeanTab(i_cext, i_cell, wID);
                    grid->setCabsMeanTab(i_cabs, i_cell, wID);
                    grid->setCscaMeanTab(i_csca, i_cell, wID);

                    if(wID == 0)
                    {
                        double number_density = dust->getNumberDensity(grid, *pp.getPositionCell());
                        grid->setNumberDensityTab(number_density, i_cell);
                        double cell_emission = dust->getTotalCellEmission(grid, pp);
                        grid->setTotalCellEmissionTab(cell_emission, i_cell);
                    }
                }
            }
        #endif

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

        CRandomGenerator rand_gen = CRandomGenerator();
        // just an arbitrary, random number for the RNG seed
        // this is NOT the seed for KISS
        ullong seed = -1ULL;

        // Create parallel region with individual RNGs and seeds for each thread
        #pragma omp parallel firstprivate(rand_gen, seed)
        {
        #ifdef _OPENMP
            seed *= -1 * (omp_get_thread_num() + 1);
        #endif
        rand_gen.init(seed);

        // A loop for each wavelength
#pragma omp for schedule(dynamic) collapse(2)
        for(int wID = 0; wID < int(nr_used_wavelengths); wID++)
        {
            // A loop for each photon
            for(llong i_phot = 0; i_phot < llong(nr_of_photons); i_phot++)
            {
                // Init variables
                double end_tau, Cext, Csca;
                Vector3D old_pos;

                // Increase counter used to show progress
#pragma omp atomic update
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
                photon_package pp = photon_package();
                pp.setPhotonID(i_phot);

                if(use_energy_density)
                    pp.setWavelength(dust->getWavelength(wID), wID);

                // Launch a new photon package from the source
                tm_source->createNextRay(&pp, &rand_gen);

                if(pp.getStokesVector()->I() < 1e-200)
                {
                    #pragma omp atomic update
                    kill_counter++;
                    continue;
                }

                if(!grid->positionPhotonInGrid(&pp))
                    if(!grid->findStartingPoint(&pp))
                    {
                        #pragma omp atomic update
                        kill_counter++;
                        continue;
                    }

                // Get tau for first interaction
                end_tau = -log(1.0 - rand_gen.getRND());

                // Save the old position to use it again
                old_pos = pp.getPosition();

                // Init current number of interactions
                ullong interactions = 0;

                // Init variables
                double tmp_tau, len, dens;

                // Transfer photon through grid
                while(grid->next(&pp))
                {
                    // If max interactions is reached, end photon transfer
                    if(interactions >= MAX_INTERACTION_RADFIELD)
                    {
                        #pragma omp atomic update
                        kill_counter++;
                        break;
                    }

                    // Get necessary quantities
                    dens = dust->getNumberDensity(grid, pp);

                    // skip cells without density
                    if(dens == 0)
                    {
                        old_pos = pp.getPosition();
                        continue;
                    }

                    // Get distance to next cell
                    len = pp.getTmpPathLength();

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
                        pp.adjustPosition(old_pos, len * end_tau / tmp_tau);

                        // Update data in grid like spectral length or radiation field
                        updateRadiationField(&pp);

                        if(!doMRWStepBW(&pp))
                        {
                            // Calculate the dust scattering cross section (for random
                            // alignment)
                            Csca = dust->getCscaMean(grid, pp);

                            // Calculate albedo and check if absorption or
                            // scattering occurs
                            double albedo = Csca / Cext;

                            if(rand_gen.getRND() < albedo)
                            {
                                // Perform simple photon scattering without
                                // changing the Stokes vectors
                                dust->scatter(grid, &pp, &rand_gen);
                            }
                            else
                            {
                                // Calculate the temperature of the absorbing cell
                                // and change the wavelength of the photon
                                if(!disable_reemission &&
                                   dust->adjustTempAndWavelengthBW(grid, &pp, use_energy_density, &rand_gen))
                                {
                                    // Send this photon into a new random direction
                                    pp.setRandomDirection(rand_gen.getRND(), rand_gen.getRND());
                                    // And reset the coord system
                                    pp.initCoordSystem();

                                    // reemission is unpolarised
                                    // reset Stokes Q, U and V
                                    pp.getStokesVector()->setQ(0.0);
                                    pp.getStokesVector()->setU(0.0);
                                    pp.getStokesVector()->setV(0.0);
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
                        end_tau = -log(1.0 - rand_gen.getRND());
                    }
                    else
                    {
                        // Update data in grid like spectral length or radiation field
                        updateRadiationField(&pp);

                        // Remove the traveled distance from optical depth
                        end_tau -= tmp_tau;
                    }
                    // Save photon position to adjust it if necessary
                    old_pos = pp.getPosition();
                }
            } // end of photon loop
        } // end of wavelength loop
        } //end of parallel block
    } // end of source loop

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

bool CRadiativeTransfer::calcMonteCarloLvlPopulation(uint i_species, uint global_seed)
{
    // Check the source
    if(sources_mc.size() != 1 || sources_mc[0]->getID() != SRC_GAS_LVL)
    {
        cout << "Error: Level population can only be calculated with the GAS MC source!" << endl;
        return false;
    }

    // Init variables
    uint kill_counter = 0;

    // Iteration counter
    uint global_iteration_counter = 0;

    // Are level populations for this transition and cell converged (globally)
    bool global_converged = false;

    // Number of cells in grid
    ulong nr_of_cells = grid->getMaxDataCells();

    // Init source from sources list
    CSourceBasic * tm_source = sources_mc[0];

    // Init luminosity and probability for a given wavelength to be chosen
    cout << CLR_LINE;
    if(!tm_source->initSource(0, 1))
        return false;

    // Number of photons per wavelength
    ullong nr_of_photons = tm_source->getNrOfPhotons();

    // Get number of spectral line transitions that have to be simulated
    uint nr_of_transitions = gas->getNrOfTransitions(i_species);
    // With Zeeman subtransitions
    uint nr_of_total_transitions = gas->getNrOfTotalTransitions(i_species);

    // If no spectral line transitions are chosen for the current gas species, skip
    if(nr_of_total_transitions == 0)
        return false;

    // Calculate the level populations for each cell (initial guess, LTE)
    if(!gas->calcLevelPopulation(grid, i_species))
    {
        cout << "\nERROR: Level population cannot be calculated!";
        return false;
    }

    // Calculate the line broadening for each cell
    gas->calcLineBroadening(grid, i_species);

    // Init progress visualization
    cout << CLR_LINE;
    cout << "-> Calculating MC level population (global = 1, max local = 1) : 0.0 [%]    \r" << flush;

    // Make global iterations till converged
    while(!global_converged && global_iteration_counter < MC_LVL_POP_MAX_GLOBAL_ITER)
    {
        // Init variables
        ullong per_counter = 0;
        int max_local_iterations = 1;

        // Increase counter
        global_iteration_counter++;

        CRandomGenerator rand_gen = CRandomGenerator();
        // just an arbitrary, random number for the RNG seed
        // this is NOT the seed for KISS
        ullong seed = -1ULL;

        // Create parallel region with individual RNGs and seeds for each thread
        #pragma omp parallel firstprivate(rand_gen, seed)
        {
        #ifdef _OPENMP
            seed *= -1 * (omp_get_thread_num() + 1);
        #endif
        rand_gen.init(seed);

        // A loop for each cell
#pragma omp for schedule(dynamic)
        for(long i_cell = 0; i_cell < long(nr_of_cells); i_cell++)
        {
            // Increase counter used to show progress
#pragma omp atomic update
            per_counter++;

            // Pointer to final cell
            cell_basic * final_cell = grid->getCellFromIndex(i_cell);

            // Maximum velocity to either side (2 x Gauss width)
            double max_velocity = 2.0 / grid->getGaussA(*final_cell);

            // Obtain magnetic field strength and orientation
            double cos_theta = 0, sin_theta = 0;
            //double cos_2_phi = 0, sin_2_phi = 0;
            Vector3D mag_field;

            // Each cell automatically converged except the last one for local iteration
            // See Pavlyuchenkov & Shustov (2003)
            double * J_nu_total = new double[nr_of_total_transitions];
            double * J_nu_in = new double[nr_of_total_transitions];
            double * J_nu_in_old = new double[nr_of_total_transitions];

            // Init pointer arrays with zero
            for(uint i_trans = 0; i_trans < nr_of_total_transitions; i_trans++)
            {
                J_nu_total[i_trans] = 0;
                J_nu_in[i_trans] = 0;
                J_nu_in_old[i_trans] = 0;
            }

            // Are level populations for this transition and cell converged
            bool local_converged = false;

            // Will be true when only the inner rad field will be calculated
            bool only_J_in = false;

            // Local iteration counter (first two runs required for initial guess of rad field)
            int local_iteration_counter = -2;

            while(!local_converged && local_iteration_counter < MC_LVL_POP_MAX_LOCAL_ITER)
            {
                // Increase counter
                local_iteration_counter++;

                // Show only new percentage number
#pragma omp critical
                {
                    cout << "-> Calculating MC level population (global = " << global_iteration_counter
                         << ", max local = " << max_local_iterations
                         << ") : " << 100 * float(per_counter) / float(nr_of_cells) << " [%]    \r" << flush;
                }

                // Perform radiative transfer for each chosen spectral line transition
                // -> including Zeeman sublevel!
                for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
                {
                    // Get rest frequency of current transition
                    double trans_frequency = gas->getTransitionFrequency(i_species, i_trans);

                    // A loop for each photon
                    for(llong i_phot = 0; i_phot < llong(nr_of_photons); i_phot++)
                    {
                        // Init photon package
                        photon_package pp =
                            photon_package(trans_frequency, dust->getWavelengthID(con_c / trans_frequency));
                        pp.setPhotonID(i_phot);

                        // Init photon package outside the grid or at the border of final cell
                        tm_source->createNextRayToCell(&pp, &rand_gen, i_cell, only_J_in);

                        // Position photon at the grid border
                        if(!grid->positionPhotonInGrid(&pp))
                        {
                            if(!grid->findStartingPoint(&pp))
                            {
#pragma omp atomic update
                                kill_counter++;
                                continue;
                            }
                        }

                        // Calculate random frequency
                        double velocity = (rand_gen.getRND() * 2 - 1) * max_velocity +
                                          gas->getProjCellVelocity(grid, pp, pp.getBackupPosition());

                        // Set frequency of the photon package
                        pp.setVelocity(velocity);

                        // Set external starting energy to CMB for total rad field
                        if(!only_J_in)
                        {
                            double tmp_energy = CMathFunctions::planck_hz(trans_frequency, 2.75);
                            pp.setStokesVector(StokesVector(tmp_energy, 0, 0, 0));
                        }

                        // Set/Get Necessary information about velocity field interpolation
                        VelFieldInterp vel_field_interp;
                        vel_field_interp.start_pos = pp.getPosition();

                        // Precalculate the velocity interpolation
                        preCalcVelocityInterp(grid, pp, &vel_field_interp);

                        // Transport the photon package through the model to the current final cell
                        while(!pp.reachedBackupPosition() && grid->next(&pp))
                        {
                            // Solve radiative transfer equation by raytracing through a cell
                            rayThroughCellForLvlPop(&pp, i_species, i_trans, vel_field_interp);
                        }

                        // Add energy to the corresponding radiation field entry
                        gas->applyRadiationFieldFactor(i_species,
                                                       i_trans,
                                                       sin_theta,
                                                       cos_theta,
                                                       pp.getStokesVector()->I() / nr_of_photons,
                                                       only_J_in ? J_nu_in : J_nu_total);
                    }
                }

                if(local_iteration_counter >= 1)
                {
                    // Start with converged and set to false if error too high
                    local_converged = true;

                    // Each transition need to converge
                    for(uint i_trans_tot = 0; i_trans_tot < nr_of_total_transitions; i_trans_tot++)
                    {
                        // Check if limit is reached
                        if(abs(J_nu_in[i_trans_tot] - J_nu_in_old[i_trans_tot]) >
                           MC_LVL_POP_DIFF_LIMIT * abs(J_nu_in[i_trans_tot] + J_nu_in_old[i_trans_tot]) +
                               MC_LVL_POP_LIMIT)
                        {
                            // Not converged
                            local_converged = false;
                            break;
                        }
                    }

                    // If not converged, update full radiation field
                    if(!local_converged)
                    {
                        for(uint i_trans_tot = 0; i_trans_tot < nr_of_total_transitions; i_trans_tot++)
                            J_nu_total[i_trans_tot] += J_nu_in[i_trans_tot] - J_nu_in_old[i_trans_tot];
                    }
                }

                // Not at the first time, since J_in has to use the initial guess at the first time
                if(only_J_in)
                {
                    // Update level populations
                    gas->updateLevelPopulation(grid, final_cell, i_species, J_nu_total);
                }
                else
                {
                    // external radiation field was calculated, only J_in from now on
                    only_J_in = true;
                }

                // Backup last iteration and reset J_in for next one
                for(uint i_trans_tot = 0; i_trans_tot < nr_of_total_transitions; i_trans_tot++)
                {
                    J_nu_in_old[i_trans_tot] = J_nu_in[i_trans_tot];
                    J_nu_in[i_trans_tot] = 0;
                }
            }

            // Find maximum local iterations
            if(local_iteration_counter > max_local_iterations)
                max_local_iterations = local_iteration_counter;

            // Delete the pointer arrays
            delete[] J_nu_total;
            delete[] J_nu_in;
            delete[] J_nu_in_old;
        }
        }//end of parallel block

        // If no local iterations were necessary, globally converged
        if(max_local_iterations == 1)
            global_converged = true;
    }

    if(global_iteration_counter == MC_LVL_POP_MAX_GLOBAL_ITER)
    {
        cout << CLR_LINE;
        cout << "HINT: Global iteration not reached! Continue, but results might be wrong!" << endl;
    }

    // Format prints
    cout << CLR_LINE;
    cout << SEP_LINE;

    // Show amount of killed photons
    if(kill_counter > 0)
        cout << "- Photons killed                    : " << kill_counter << endl;
    cout << "-> Calculating MC level population : done                          " << endl;
    cout << "    Number of global iterations used  : " << global_iteration_counter << endl;
    cout << "    Seed used for MC level population : " << global_seed << endl;

    return true;
}

void CRadiativeTransfer::rayThroughCellForLvlPop(photon_package * pp,
                                                 uint i_species,
                                                 uint i_trans,
                                                 const VelFieldInterp & vel_field_interp)
{

    // Get gas species density from grid
    double dens_species = gas->getNumberDensity(grid, *pp, i_species);

    // Perform radiative transfer only if the temperature of the current species
    // are not negligible
    if(dens_species > 1e-200)
    {
        // Init matrix for absorption and dust emissivity
        Matrix2D total_absorption_matrix(4, 4);
        StokesVector dust_emi_and_ext;

        // Get extra information about the magnetic field and ine broadening
        MagFieldInfo mag_field_info;
        LineBroadening line_broadening;
        uint i_zeeman = gas->getZeemanSplitIndex(i_species, i_trans);
        if(i_zeeman != MAX_UINT)
        {
            grid->getMagFieldInfo(*pp, &mag_field_info);
            grid->getLineBroadening(*pp, i_zeeman, &line_broadening);
        }
        else
        {
            // Set only gauss_a if not zeeman split
            line_broadening.gauss_a = grid->getGaussA(*pp);
        }

        // Get the path length through the current cell
        double len = pp->getTmpPathLength();

        // Get necessary quantities from the current cell
        // double dens_gas = grid->getGasNumberDensity(*pp);

        // Calculate the emission of the dust grains
        dust->calcEmissivityHz(grid, *pp, &dust_emi_and_ext);

        // Init variables
        double cell_sum = 0, cell_d_l = len;
        ullong kill_counter = 0;

        // Get entry position of the current cell
        Vector3D pos_xyz_cell = pp->getPosition() - (len * pp->getDirection());

        // Make sub steps until cell is completely crossed
        // If the error of a sub step is too high, make the step smaller
        while(cell_sum < len)
        {
            // Increase the kill counter
            kill_counter++;

            // If too many sub steps are needed, kill the photon
            if(kill_counter >= MAX_SOLVER_STEPS)
            {
                cout << "\nWARNING: Solver steps > " << MAX_SOLVER_STEPS << ". Too many steps!" << endl;
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
                // Calculate projected velocity on each Runge-Kutta position
                double rel_velocity =
                    pp->getVelocity() -
                    gas->getProjCellVelocityInterp(pos_xyz_cell + cell_d_l * pp->getDirection() * RK_c[k],
                                                   pp->getDirection(),
                                                   vel_field_interp);

                // Init emission
                StokesVector total_emission;

                // Get line emissivity (also combined Zeeman lines)
                gas->calcEmissivity(grid,
                                    *pp,
                                    i_species,
                                    i_trans,
                                    rel_velocity,
                                    line_broadening,
                                    mag_field_info,
                                    &total_emission,
                                    &total_absorption_matrix);

                // Combine the Stokes vectors from gas and dust for emission
                total_emission *= dens_species;
                total_emission += dust_emi_and_ext;
                // and extinction
                total_absorption_matrix *= dens_species;
                if(dust_emi_and_ext.T() != 0)
                    for(uint i = 0; i < 4; i++)
                        total_absorption_matrix(i, i) += dust_emi_and_ext.T();
                total_absorption_matrix *= -1;

                // Init scalar product
                StokesVector scalar_product;

                // Calculate multiplication between Runge-Kutta parameter
                for(uint i = 0; i <= k; i++)
                    scalar_product += (RK_k[i] * RK_a(i, k));

                // Calculate new Runge-Kutta parameters as the result of the
                // radiative transfer equation at the Runge-Kutta sub positions
                RK_k[k] = total_absorption_matrix * (scalar_product * cell_d_l + *pp->getStokesVector()) +
                          total_emission;
            }

            // Init two temporary Stokes vectors
            StokesVector stokes_new = *pp->getStokesVector();
            StokesVector stokes_new2 = *pp->getStokesVector();

            for(uint i = 0; i < 6; i++)
            {
                stokes_new += RK_k[i] * cell_d_l * RK_b1[i];
                stokes_new2 += RK_k[i] * cell_d_l * RK_b2[i];
            }

            // Delete the Runge-Kutta pointer
            delete[] RK_k;

            // Ignore very small values
            if(abs(stokes_new.I()) < 1e-200)
                stokes_new.resetIntensity();
            if(abs(stokes_new2.I()) < 1e-200)
                stokes_new2.resetIntensity();

            // Calculate the difference between the results with two
            // different precisions to see if smaller steps are needed
            double epsi, dz_new;
            calcStepWidth(stokes_new, stokes_new2, cell_d_l, &epsi, &dz_new);

            // Is a smaller step width needed
            if(epsi <= 1.0)
            {
                // Stokes_new is the current flux of this line-of-sight
                pp->setStokesVector(stokes_new);

                // Update the position of the photon package
                pos_xyz_cell += cell_d_l * pp->getDirection();

                // Increase the sum of the cell path lengths
                cell_sum += cell_d_l;

                // Find a new path length
                cell_d_l = min(dz_new, 4 * cell_d_l);

                // If the new step would exceed the cell, make it smaller
                if(cell_sum + cell_d_l > len)
                    cell_d_l = len - cell_sum;

                // Local iteration for final cell (start converging)
                if(pp->reachedBackupPosition(pos_xyz_cell))
                    break;
            }
            else
            {
                // Find a smaller path length
                cell_d_l = max(dz_new, 0.25 * cell_d_l);
            }
        }
    }
}

/*bool CRadiativeTransfer::setTemperatureDistribution()
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

        double rho = dust->getNumberDensity(grid, *cell) * 1e-16; // dens_data[c];
        double tdust = 1000; // grid->getGasTemperature(cell); //tdust_orig;
        double tdust_orig = grid->getGasTemperature(*cell);
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
}*/

bool CRadiativeTransfer::calcPolMapsViaMC()
{
    // Init variables
    ullong nr_of_photons;
    ullong per_counter;
    ullong nr_of_wavelength = dust->getNrOfWavelength();
    float last_percentage;
    uint mrw_counter = 0;
    ullong kill_counter = 0;
    uint max_source = uint(sources_mc.size());
    uint nr_laser_sources = 0;

    // First, show optical depth from each source to each detector at each wavelength
    for(uint s = 0; s < max_source; s++)
    {
        for(uint wID = 0; wID < nr_of_wavelength; wID++)
        {
            double wavelength = dust->getWavelength(wID);
            for(uint d = 0; d < nr_mc_detectors; d++)
            {
                // Get index of wavelength in current detector
                uint wID_det = detector[d].getDetectorWavelengthID(wavelength);
                // Only calculate for detectors with the corresponding wavelengths
                if(wID_det != MAX_UINT)
                {
                    photon_package pp_tau = photon_package();
                    pp_tau.setWavelength(wavelength, wID);
                    pp_tau.setPosition(sources_mc[s]->getPosition());
                    pp_tau.setDirection(detector[d].getDirection());

                    // Position the photon inside the grid
                    if(!grid->positionPhotonInGrid(&pp_tau))
                        if(!grid->findStartingPoint(&pp_tau))
                            continue;

                    double tau = getOpticalDepthAlongPath(&pp_tau);
                    cout << CLR_LINE;
                    cout << "- Optical depth (wavelength = " << wavelength << " m) from source " << s+1
                        << " to detector " << d+1 << ": " << tau << endl;
                }
            }
        }
    }

    // Perform Monte-Carlo radiative transfer for each chosen source
    for(uint s = 0; s < max_source; s++)
    {
        // Init source object
        CSourceBasic * tm_source = sources_mc[s];
        cout << CLR_LINE;

        #if (USE_PRECALC_TABLE)
            grid->initPreCalcTables(nr_of_wavelength);
        #endif

        // Perform Monte-Carlo radiative transfer for
        // each chosen wavelength of the chosen source
        for(uint wID = 0; wID < nr_of_wavelength; wID++)
        {
            // Init source parameters for scattering maps
            if(!tm_source->initSource(wID))
                continue;

            if(tm_source->getID() == SRC_LASER)
                nr_laser_sources += 1;

            // Get number of photons that have to be emitted from the current source
            nr_of_photons = tm_source->getNrOfPhotons();

            #if (USE_PRECALC_TABLE)
                ulong nr_of_cells = grid->getMaxDataCells();

#pragma omp parallel for schedule(dynamic)
                for(long i_cell = 0; i_cell < long(nr_of_cells); i_cell++)
                {
                    photon_package pp = photon_package();

                    pp.setWavelength(dust->getWavelength(wID), wID);
                    // Put photon package into current cell
                    pp.setPositionCell(grid->getCellFromIndex(i_cell));

                    double i_cext = dust->getCextMean(grid, pp);
                    double i_cabs = dust->getCabsMean(grid, pp);
                    double i_csca = dust->getCscaMean(grid, pp);
                    grid->setCextMeanTab(i_cext, i_cell, wID);
                    grid->setCabsMeanTab(i_cabs, i_cell, wID);
                    grid->setCscaMeanTab(i_csca, i_cell, wID);

                    if(wID == 0)
                    {
                        double number_density = dust->getNumberDensity(grid, *pp.getPositionCell());
                        grid->setNumberDensityTab(number_density, i_cell);
                        double cell_emission = dust->getTotalCellEmission(grid, pp);
                        grid->setTotalCellEmissionTab(cell_emission, i_cell);
                    }
                }
            #endif

            // Init progress visualization
            cout << "-> MC pol. map(s) (source ID: " << s + 1 << ", wavelength: " << dust->getWavelength(wID)
                 << " [m], photons: " << float(nr_of_photons) << ") 0 [%]   \r" << flush;

            // Init counter and percentage to show progress
            per_counter = 0;
            last_percentage = 0;

            CRandomGenerator rand_gen = CRandomGenerator();
            // just an arbitrary, random number for the RNG seed
            // this is NOT the seed for KISS
            ullong seed = -1ULL;

            // Create parallel region with individual RNGs and seeds for each thread
            #pragma omp parallel firstprivate(rand_gen, seed)
            {
            #ifdef _OPENMP
                seed *= -1 * (omp_get_thread_num() + 1);
            #endif
            rand_gen.init(seed);

            // Perform radiative transfer through the model for each photon
            #pragma omp for schedule(dynamic)
            for(llong i_phot = 0; i_phot < llong(nr_of_photons); i_phot++)
            {
                // Init cross sections
                double Cext;

                // Increase counter used to show progress
                #pragma omp atomic update
                per_counter++;

                // Calculate percentage of total progress per source
                float percentage = 100.0 * float(per_counter) / float(nr_of_photons);

                // Show only new percentage number if it changed
                if((percentage - last_percentage) >= PERCENTAGE_STEP)
                {
                    #pragma omp critical
                    {
                        cout << "-> MC pol. map(s) (source ID: " << s + 1
                             << ", wavelength: " << dust->getWavelength(wID)
                             << " [m], photons: " << float(nr_of_photons) << ") " << percentage << " [%]   \r"
                             << flush;
                        last_percentage = percentage;
                    }
                }

                // Init the photon_package pp
                photon_package pp = photon_package();
                pp.setPhotonID(i_phot);

                // Init variables
                double end_tau, tau_tot;
                Vector3D old_pos, start_pos;

                // Set current wavelength
                pp.setWavelength(dust->getWavelength(wID), wID);

                // Launch a new photon package from the source
                tm_source->createNextRay(&pp, &rand_gen);

                // Position the photon inside the grid
                if(!grid->positionPhotonInGrid(&pp))
                    if(!grid->findStartingPoint(&pp))
                    {
                        #pragma omp atomic update
                        kill_counter++;
                        continue;
                    }

                // Set starting position for enforced scattering
                start_pos = pp.getPosition();

                // Save position of last interaction to know to which pixel
                // the photon belongs, if it is not scattered on its further
                // path through the model
                pp.updateBackupPosition();

                // if enforced first scattering
                if(b_forced)
                {
                    // Get tau for first interaction, if the interaction is forced
                    tau_tot = getOpticalDepthAlongPath(&pp);
                    end_tau = -log(1.0 - rand_gen.getRND() * (1.0 - exp(-tau_tot)));

                    // If optical depth is exactly zero, send photon package as without enfsca
                    if(end_tau == 0)
                        end_tau = -log(1.0 - rand_gen.getRND());
                }
                else
                {
                    // Get tau for first interaction
                    end_tau = -log(1.0 - rand_gen.getRND());
                }

                // Init variables
                ullong interactions = 0;
                double tmp_tau;
                double len, dens;

                // Save photon position to adjust it if necessary
                old_pos = pp.getPosition();

                // Transfer photon through grid
                while(grid->next(&pp))
                {
                    // If max interactions is reached or the photon intensity
                    // is too low, end photon transfer
                    if(interactions >= MAX_INTERACTION_DUST_MC || pp.getStokesVector()->I() < 1e-200)
                    {
                        #pragma omp atomic update
                        kill_counter++;
                        break;
                    }

                    // Get dust number density of the current cell
                    dens = dust->getNumberDensity(grid, pp);

                    // If the dust density is too low, skip this cell
                    if(dens < 1e-200)
                    {
                        old_pos = pp.getPosition();
                        continue;
                    }

                    // Calculate the dust absorption cross section (for random
                    // alignment)
                    Cext = dust->getCextMean(grid, pp);

                    // Get path length through current cell
                    len = pp.getTmpPathLength();

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
                        pp.adjustPosition(old_pos, len);

                        // Modify second photon if enforced scattering is used
                        if(b_forced && interactions == 1)
                        {
                            if(!peel_off)
                            {
                                for(uint d = 0; d < nr_mc_detectors; d++)
                                {
                                    // Get index of wavelength in current detector
                                    uint wID_det = detector[d].getDetectorWavelengthID(pp.getWavelength());

                                    // Only calculate for detectors with the corresponding wavelengths
                                    if(wID_det != MAX_UINT)
                                    {
                                        // Check photon direction to observer for each detector
                                        if(photonInDetectorDir(&pp, &detector[d]))
                                        {
                                            photon_package pp_enf = photon_package();
                                            pp_enf.setWavelength(pp.getWavelength(), pp.getDustWavelengthID());
                                            pp_enf.setPosition(pp.getPosition());
                                            pp_enf.setPositionCell(pp.getPositionCell());
                                            pp_enf.setBackupPosition(start_pos);
                                            pp_enf.setDirection(pp.getDirection());
                                            pp_enf.initCoordSystem();
                                            pp_enf.setStokesVector(*pp.getStokesVector() * exp(-tau_tot));
                                            scaleAddToDetector(&pp_enf, &detector[d], 0, tm_source->getID());
                                        }
                                    }
                                }
                            }
                            *pp.getStokesVector() *= (1.0 - exp(-tau_tot));
                        }

                        // Save position of last interaction to know to which pixel
                        // the photon belongs, if it is not scattered on its further
                        // path through the model
                        pp.updateBackupPosition();

                        if(!doMRWStepBWWithoutHeating(&pp))
                        {
                            // Calculate the dust scattering cross section (for random
                            // alignment)
                            // Csca = dust->getCscaMean(grid, pp);

                            // If peel-off is used, add flux to the detector
                            // except it is the first interaction of the non-forced photon
                            if(peel_off)
                            {
                                // Transport a separate photon to each detector
                                for(uint d = 0; d < nr_mc_detectors; d++)
                                {
                                    // Get index of wavelength in current detector
                                    uint wID_det =
                                        detector[d].getDetectorWavelengthID(pp.getWavelength());

                                    // Only calculate for detectors with the
                                    // corresponding wavelengths
                                    if(wID_det != MAX_UINT)
                                    {
                                        // Init photon package
                                        photon_package pp_escape = photon_package();

                                        // Create an escaping photon into the
                                        // direction of the detector
                                        dust->getEscapePhoton(
                                            grid, &pp,
                                            detector[d].getEX(),
                                            detector[d].getDirection(),
                                            &pp_escape,
                                            &rand_gen);

                                        // optical depth towards observer
                                        double tau_obs = getOpticalDepthAlongPath(&pp_escape);

                                        // Reduce the Stokes vector by the optical depth
                                        *pp_escape.getStokesVector() *= exp(-tau_obs);

                                        // Convert the flux into Jy (not if laser is used) and consider
                                        // the distance to the observer
                                        double distance = detector[d].getDistance();
                                        if(tm_source->getID() == SRC_LASER)
                                            *pp_escape.getStokesVector() /= distance * distance;
                                        else
                                            CMathFunctions::lum2Jy(
                                                pp_escape.getStokesVector(),
                                                pp.getWavelength(),
                                                distance);

                                        // Consider foreground extinction
                                        *pp_escape.getStokesVector() *=
                                            dust->getForegroundExtinction(pp.getWavelength());

                                        // If the photon intensity is too low, end
                                        // photon transfer (for better R and VOV statistics)
                                        if(pp_escape.getStokesVector()->I() < 1e-200)
                                            continue;

                                        // Add the photon package to the detector
                                        detector[d].addToMonteCarloDetector(
                                            pp_escape, wID_det, SCATTERED_DUST, tm_source->getID());
                                    }
                                }
                            }
                        }
                        else
                        {
                            #pragma omp atomic update
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
                        dust->scatter(grid, &pp, &rand_gen, true);

                        // Calculate new optical depth for next interaction
                        end_tau = -log(1.0 - rand_gen.getRND());
                    }
                    else
                    {
                        // Remove the traveled distance from optical depth
                        end_tau -= tmp_tau;
                    }
                    // Save photon position to adjust it if necessary
                    old_pos = pp.getPosition();
                }

                // If peel-off is not used, use classic Monte-Carlo method
                // Now, the photon has left the model space
                if(!peel_off && interactions <= MAX_INTERACTION_DUST_MC)
                {
                    // Move photon back to the point of last interaction
                    pp.resetPositionToBackupPos();

                    // Transport photon to observer for each detector
                    for(uint d = 0; d < nr_mc_detectors; d++)
                    {
                        // Get index of wavelength in current detector
                        uint wID_det = detector[d].getDetectorWavelengthID(pp.getWavelength());

                        // Only calculate for detectors with the corresponding wavelengths
                        if(wID_det != MAX_UINT)
                        {
                            // Check photon direction to observer for each detector
                            if(photonInDetectorDir(&pp, &detector[d]))
                                scaleAddToDetector(&pp, &detector[d], interactions, tm_source->getID());
                        }
                    }
                }
            }//end of photon loop
            }//end of parallel block

            // If peel-off is used, add direct source emission to each detector
            if(peel_off && tm_source->getID() != SRC_DUST)
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
                        photon_package pp_direct = photon_package();

                        // Set current wavelength to temporary photon package
                        pp_direct.setWavelength(dust->getWavelength(wID), wID);

                        // Get direction to the current detector
                        Vector3D dir_obs = detector[d].getDirection();

                        // Launch a new photon package from the source
                        tm_source->createDirectRay(&pp_direct, &rand_gen, dir_obs);

                        // Position the photon inside the grid
                        if(!grid->positionPhotonInGrid(&pp_direct))
                            if(!grid->findStartingPoint(&pp_direct))
                                continue;

                        // Rotate photon package into the coordinate space of the detector
                        double rot_angle_phot_obs = CMathFunctions::getRotationAngleObserver(
                            detector[d].getEX(), pp_direct.getEX(), -1*pp_direct.getEY());
                        pp_direct.getStokesVector()->rot(rot_angle_phot_obs);

                        // The scattering part is based on O. Fischer (1993)
                        // But on our detectors, U is defined the other way round
                        pp_direct.getStokesVector()->multU(-1);

                        // optical depth towards observer
                        double tau_obs = getOpticalDepthAlongPath(&pp_direct);

                        // Reduce the Stokes vector by the optical depth
                        *pp_direct.getStokesVector() *= exp(-tau_obs);

                        // Convert the flux into Jy and consider the distance to the
                        // observer (not if laser is used)
                        if(tm_source->getID() != SRC_LASER)
                            CMathFunctions::lum2Jy(
                                pp_direct.getStokesVector(),
                                pp_direct.getWavelength(),
                                detector[d].getDistance());

                        // Consider foreground extinction
                        *pp_direct.getStokesVector() *=
                            dust->getForegroundExtinction(pp_direct.getWavelength());
                            
                        // Add the photon package to the detector
                        detector[d].addToMonteCarloDetector(pp_direct, wID_det, DIRECT_STAR, tm_source->getID());
                    }
                }
            }
        }

        if(peel_off && tm_source->getID() == SRC_DUST)
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
        if(!detector[d].writeMap(d, RESULTS_MC, nr_laser_sources))
            return false;
        if(!detector[d].writeMapStats(d, RESULTS_MC, nr_laser_sources))
            return false;
        if(!detector[d].writeSed(d, RESULTS_MC, nr_laser_sources))
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
    cout << "- Calculation of MC polarization maps (photons: " << float(nr_of_photons) << "): done" << endl;

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
        double gas_dens = grid->getGasNumberDensity(*cell);
        if(gas_dens < min_gas_density)
            continue;

        dust->convertTempInQB(grid, cell, min_gas_density, use_gas_temp);

#pragma omp atomic update
        pos_counter++;

//         if(pos_counter % 10000 == 0)
//         {
// #pragma omp critical
//             {
//                 cout << "-> Converting emissivities: " << 100.0 * float(pos_counter) / float(max_cells)
//                      << " [%]       \r";
//             }
//         }
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

#pragma omp atomic update
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
            grid->setGasTemperature(cell, adjTgas * grid->getDustTemperature(*cell));

#pragma omp atomic update
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

#pragma omp atomic update
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

double CRadiativeTransfer::getOpticalDepthAlongPath(photon_package * pp)
{
    photon_package pp_new = photon_package();
    pp_new.setWavelength(pp->getWavelength(), pp->getDustWavelengthID());
    pp_new.setPosition(pp->getPosition());
    pp_new.setPositionCell(pp->getPositionCell());
    pp_new.setDirection(pp->getDirection());
    pp_new.setDirectionID(pp->getDirectionID());

    double len, dens, Cext, tau = 0.0;
    while(grid->next(&pp_new))
    {
        // Get the traveled distance
        len = pp_new.getTmpPathLength();
        // Get the current density
        dens = dust->getNumberDensity(grid, pp_new);
        // Get the current mean extinction cross-section
        Cext = dust->getCextMean(grid, pp_new);
        // Add the optical depth of the current path to the total optical depth
        tau += Cext * len * dens;
    }
    return tau;
}

bool CRadiativeTransfer::photonInDetectorDir(photon_package * pp, CDetector * detector)
{
    // Get direction to the current detector
    Vector3D dir_obs = detector->getDirection();

    // Calculate the angle between the photon travel direction
    // and the detector direction
    double cos_angle = pp->getDirection() * dir_obs;

    // Get acceptance angle from detector (minimum 1°)
    double cos_acceptance_angle = detector->getAcceptanceAngle();

    // If the angle angle between the photon direction and the
    // detector direction is smaller than the acceptance angle
    if(cos_angle >= cos_acceptance_angle)
        return true;

    return false;
}

void CRadiativeTransfer::scaleAddToDetector(photon_package * pp, CDetector * detector, ullong interactions, uint sourceID)
{
    // Get the angle to rotate the photon space into the
    // detector space
    double rot_angle_phot_obs = CMathFunctions::getRotationAngleObserver(
        detector->getEX(), pp->getEX(), -1*pp->getEY());

    // Rotate the Stokes vector
    pp->getStokesVector()->rot(rot_angle_phot_obs);

    // The scattering part is based on O. Fischer (1993)
    // But on our detectors, U is defined the other way round
    pp->getStokesVector()->multU(-1);

    // Consider the greater solid angle due
    // to the acceptance angle
    if(sourceID != SRC_LASER || interactions > 0)
    {
        double cos_acceptance_angle = detector->getAcceptanceAngle();
        *pp->getStokesVector() *= 1.0 / ((1.0 - cos_acceptance_angle) * PIx2);
    }

    // Convert the flux into Jy and consider
    // the distance to the observer (not if laser is used)
    double distance = detector->getDistance();
    if(sourceID == SRC_LASER)
    {
        if(interactions > 0)
            *pp->getStokesVector() /= distance * distance;
    }
    else
        CMathFunctions::lum2Jy(pp->getStokesVector(), pp->getWavelength(), distance);

    // Consider foreground extinction
    *pp->getStokesVector() *=
        dust->getForegroundExtinction(pp->getWavelength());

    uint wID_det = detector->getDetectorWavelengthID(pp->getWavelength());
    
    // If the photon intensity is too low, end
    // photon transfer (for better R and VOV statistics)
    if(pp->getStokesVector()->I() < 1e-200)
        return;

    // Add the photon package to the detector
    if(interactions == 0)
        detector->addToMonteCarloDetector(*pp, wID_det, DIRECT_STAR, sourceID);
    else
        detector->addToMonteCarloDetector(*pp, wID_det, SCATTERED_DUST, sourceID);
}

// -------------------------------------------------
// ------ Calculation of Synchotron transfer -------
// -------------------------------------------------

bool CRadiativeTransfer::calcSyncMapsViaRaytracing(parameters & param)
{
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
            cout << "-> Ray tracing synchrotron map(s) (Seq. " << i_det + 1 << ", source: " << sID + 1
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
#pragma omp atomic update
                per_counter++;

                // Calculate percentage of total progress per source
                float percentage = 100.0 * float(per_counter) / float(per_max);

                // Show only new percentage number if it changed
                if((percentage - last_percentage) > PERCENTAGE_STEP)
                {
#pragma omp critical
                    {
                        cout << "-> Ray tracing synchrotron map(s) (Seq. " << i_det + 1
                             << ", source: " << sID + 1 << ")  "
                             << float(100.0 * float(per_counter) / float(per_max)) << " [%]         \r"
                             << flush;
                        last_percentage = percentage;
                    }
                }
            }

            // Show final progress
            cout << "-> Ray tracing synchrotron map(s) (Seq. " << i_det + 1 << ", source: " << sID + 1
                 << ") 100 [%]       \r" << flush;

            // post-process raytracing simulation
            if(!tracer[i_det]->postProcessing())
                return false;

            // Write results either as text or fits file
            if(!tracer[i_det]->writeSyncResults())
                return false;

            if(tracer[i_det]->getSubpixelWarning())
            {
                cout << "WARNING: level of subpixeling (" << param.getMaxSubpixelLvl() << ") might be too low" << endl;
                cout << "if required, increase the maximum level of subpixeling with <max_subpixel_lvl> in the command file" << endl;
            }
            if(tracer[i_det]->getDetectorShape() == DET_PLANE && (grid->getDataID() == GRID_ID_SPH || grid->getDataID() == GRID_ID_CYL))
                cout << "HINT: a 'polar' detector should be used for a spherical or cylindrical grid" << endl;
        }
    }

    // Show that raytracing is finished
    cout << CLR_LINE;
    cout << "- Ray tracing synchrotron map    : done" << endl;

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

        // Create new photon package
        photon_package pp(tracer[i_det]->getNrExtra() * nr_used_wavelengths);

        tracer[i_det]->preparePhoton(&pp, cx, cy);

        // Calculate continuum emission along one path
        getSyncIntensity(&pp, tmp_source, cx, cy, i_det, subpixel_lvl);

        tracer[i_det]->addToDetector(&pp, i_pix);
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

void CRadiativeTransfer::getSyncIntensity(photon_package * pp,
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

    // Update Stokes vectors with emission from background source
    for(uint i_extra = 0; i_extra < tracer[i_det]->getNrExtra(); i_extra++)
    {
        for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
        {
            // Set current index in photon package
            pp->setSpectralID(i_wave + i_extra * nr_used_wavelengths);

            // Set wavelength index in photon package
            double wavelength = tracer[i_det]->getWavelength(i_wave);
            uint wID = dust->getWavelengthID(wavelength);
            pp->setWavelength(wavelength, wID);

            // Set related index of the multi wavelength map to the wavelength ID
            pp->getStokesVector()->set(tmp_source->getStokesVector(pp));
        }
    }

    // Find starting point inside the model and travel through it
    if(grid->findStartingPoint(pp))
    {
        while(grid->next(pp) && tracer[i_det]->isNotAtCenter(pp, cx, cy))
        {
            rayThroughCellSync(pp, i_det, nr_used_wavelengths);
        }
    }

    for(uint i_extra = 0; i_extra < tracer[i_det]->getNrExtra(); i_extra++)
    {
        for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
        {
            // Set current index in photon package
            pp->setSpectralID(i_wave + i_extra * nr_used_wavelengths);

            // Convert W/m2/Hz/sr to Jy
            double mult = 1e+26 * subpixel_fraction * tracer[i_det]->getDistanceFactor();

            // Include foreground extinction if necessary
            mult *= dust->getForegroundExtinction(tracer[i_det]->getWavelength(pp->getWavelength()));

            // Update the photon package with the Stokes vectors
            if(pp->getStokesVector()->I() < 0)
                pp->getStokesVector()->setI(0);

            pp->getStokesVector()->multS(mult);
            pp->getStokesVector()->multT(subpixel_fraction);
            pp->getStokesVector()->multSp(subpixel_fraction);
        }
    }
}

void CRadiativeTransfer::rayThroughCellSync(photon_package * pp, uint i_det, uint nr_used_wavelengths)
{
    double n_th = grid->getThermalElectronDensity(*pp);
    double T_e = 0.0; // grid->getElectronTemperature(pp); //reserved for later use

    double n_cr = grid->getCRElectronDensity(*pp);
    double g_min = grid->getGammaMin(*pp);
    double g_max = grid->getGammaMax(*pp);
    double pow_p = grid->getPowerLawIndex(*pp);

    double B = grid->getMagField(*pp).length();

    // If the all the electron densities are far too low, skip the current cell
    if(n_cr + n_th >= 1e-200)
    {
        // Get path length through current cell
        double len = pp->getTmpPathLength();

        // Calculate orientation of the Stokes vector in relation to the detector
        double phi = grid->getPhiMag(*pp);
        double theta = grid->getThetaMag(*pp);

        double sin_2ph = sin(2.0 * phi);
        double cos_2ph = cos(2.0 * phi);

        for(uint i_extra = 0; i_extra < tracer[i_det]->getNrExtra(); i_extra++)
        {
            for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
            {
                // Set current index in photon package
                pp->setSpectralID(i_wave + i_extra * nr_used_wavelengths);

                // Init variables
                syn_param syn_th, syn_cr, syn_ca;
                Matrix2D alpha_cr, alpha_ca;
                StokesVector S_em_cr, S_em_ca;

                // Get syn parameter
                syn_th = synchrotron->get_Thermal_Parameter(n_th, T_e, pp->getWavelength(), B, theta);
                syn_cr = synchrotron->get_Power_Law_Parameter(
                    n_cr, pp->getWavelength(), B, theta, g_min, g_max, pow_p);

                // Get matrixes with the absorption's and conversions
                alpha_cr = -1 * syn_cr.getSyncMatrix();

                // Set stokes vector with emission
                S_em_cr = StokesVector(syn_cr.j_I, syn_cr.j_Q * cos_2ph, syn_cr.j_Q * sin_2ph, syn_cr.j_V);

                if(i_extra == 1)
                {
                    syn_ca = syn_cr + syn_th;
                    alpha_ca = -1 * syn_ca.getSyncMatrix();
                    S_em_ca =
                        StokesVector(syn_ca.j_I, syn_ca.j_Q * cos_2ph, syn_ca.j_Q * sin_2ph, syn_ca.j_V);
                }

                // Init a variable to sum up path lengths until cell is crossed
                double cell_sum = 0.0;

                // First path length is path through cell
                double cell_d_l = len;

                // Save the cell entry position of the photon package
                // (Hint: next(pp) puts the photon onto the border to the next cell)
                Vector3D pos_xyz_cell = pp->getPosition() - (len * pp->getDirection());

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
                        if(!fail)
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
                    StokesVector * RK_k = new StokesVector[6];

                    // Calculate result of the radiative transfer equation at each
                    // Runge-Kutta sub position
                    for(uint k = 0; k < 6; k++)
                    {
                        // Init scalar product
                        StokesVector scalar_product;

                        // Calculate multiplication between Runge-Kutta parameter
                        for(uint i = 0; i < 6; i++)
                        {
                            scalar_product += (RK_k[i] * RK_a(i, k));
                        }

                        // Calculate new Runge-Kutta parameters as the result of the
                        // radiative transfer equation at the Runge-Kutta sub
                        // positions
                        if(i_extra == 0)
                            RK_k[k] =
                                alpha_cr * (scalar_product * cell_d_l + *pp->getStokesVector()) + S_em_cr;
                        else
                            RK_k[k] =
                                alpha_ca * (scalar_product * cell_d_l + *pp->getStokesVector()) + S_em_ca;
                    }

                    // Init two temporary Stokes vectors
                    StokesVector stokes_new = *pp->getStokesVector();
                    StokesVector stokes_new2 = *pp->getStokesVector();

                    // Add the result at each Runge-Kutta sub position
                    // to the total Stokes vector (Using Runge-Kutta Fehlberg)
                    for(uint i = 0; i < 6; i++)
                    {
                        stokes_new += RK_k[i] * cell_d_l * RK_b1[i];
                        stokes_new2 += RK_k[i] * cell_d_l * RK_b2[i];
                    }

                    // Delete the Runge-Kutta pointer
                    delete[] RK_k;

                    // Calculate the difference between the results with two
                    // different precisions to see if smaller steps are needed
                    // (see Reissl)

                    double epsi=0, dz_new;
                    // Do approximate solution
                    if(!fail)
                        calcStepWidth(stokes_new, stokes_new2, cell_d_l, &epsi, &dz_new);
                    else if(i_extra == 1)
                    {
                        calcStepWidth(stokes_new, stokes_new2, cell_d_l, &epsi, &dz_new);

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

                        stokes_new = *pp->getStokesVector(i_wave) + StokesVector(0, Q, U, V);
                    }

                    if(epsi == 0)
                    {
                        dz_new = 0.01 * len;
                    }

                    // Is a smaller step width needed
                    if(epsi <= 1.0)
                    {
                        // Add the temporary Stokes vector to the total one
                        pp->setStokesVector(stokes_new);

                        if(i_extra == 0)
                        {
                            // tau is set to Farraday rotation*lambda^2 for CR electrons
                            // (currently its set to zero)
                            pp->getStokesVector()->addT(syn_cr.kappa_V * cell_d_l);

                            // Sp is set to the thermal electron column
                            pp->getStokesVector()->addSp(n_th * cell_d_l);
                        }
                        else
                        {
                            // tau is set to Farraday rotation*lambda^2 for all
                            pp->getStokesVector()->addT(syn_ca.kappa_V * cell_d_l);

                            // Sp is set to the CR electron column
                            pp->getStokesVector()->addSp(n_cr * cell_d_l);
                        }

                        // Update the position of the photon package
                        pos_xyz_cell += cell_d_l * pp->getDirection();

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

// -------------------------------------------
// ------ Calculation of Dust transfer -------
// -------------------------------------------

bool CRadiativeTransfer::calcPolMapsViaRaytracing(parameters & param)
{
    // Get list of detectors/sequences that will be simulated with the raytracer
    dlist dust_ray_detectors = param.getDustRayDetectors();

    CRandomGenerator rand_gen = CRandomGenerator();
    // just an arbitrary, random number for the RNG seed
    // this is NOT the seed for KISS
    ullong seed = -1ULL;
    rand_gen.init(seed);

    if(!dust_ray_detectors.empty())
    {
        for(uint i_det = start; i_det <= stop; i_det++)
        {
            // Init source object
            CSourceBasic * tmp_source;
            uint sID = tracer[i_det]->getSourceIndex();
            tmp_source = sources_ray[sID];

            #if (USE_PRECALC_TABLE)
                uint nr_used_wavelengths = tracer[i_det]->getNrOfSpectralBins();
                grid->initPreCalcTables(nr_used_wavelengths);

                ulong nr_of_cells = grid->getMaxDataCells();

#pragma omp parallel for schedule(dynamic) collapse (2)
                for(uint wID = 0; wID < nr_used_wavelengths; wID++)
                {
                    for(long i_cell = 0; i_cell < long(nr_of_cells); i_cell++)
                    {
                        photon_package pp = photon_package();

                        pp.setWavelength(dust->getWavelength(wID), wID);
                        // Put photon package into current cell
                        pp.setPositionCell(grid->getCellFromIndex(i_cell));

                        double i_cext = dust->getCextMean(grid, pp);
                        double i_cabs = dust->getCabsMean(grid, pp);
                        double i_csca = dust->getCscaMean(grid, pp);
                        grid->setCextMeanTab(i_cext, i_cell, wID);
                        grid->setCabsMeanTab(i_cabs, i_cell, wID);
                        grid->setCscaMeanTab(i_csca, i_cell, wID);

                        if(wID == 0)
                        {
                            double number_density = dust->getNumberDensity(grid, *pp.getPositionCell());
                            grid->setNumberDensityTab(number_density, i_cell);
                            double cell_emission = dust->getTotalCellEmission(grid, pp);
                            grid->setTotalCellEmissionTab(cell_emission, i_cell);
                        }
                    }
                }
            #endif

            // Calculate total number of pixel
            uint per_max = tracer[i_det]->getNpix();

            // Init counter and percentage to show progress
            ullong per_counter = 0;
            float last_percentage = 0;

            // Show information about the current detector
            cout << CLR_LINE;
            cout << "-> Ray tracing dust map(s) (Seq. " << i_det + 1 << ", source: " << sID + 1 << ") 0 [%]   \r"
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
#pragma omp atomic update
                per_counter++;

                // Calculate percentage of total progress per source
                float percentage = 100.0 * float(per_counter) / float(per_max);

                // Show only new percentage number if it changed
                if((percentage - last_percentage) > PERCENTAGE_STEP)
                {
#pragma omp critical
                    {
                        cout << "-> Ray tracing dust map(s) (Seq. " << i_det + 1 << ", source: " << sID + 1
                             << ") " << percentage << " [%]       \r" << flush;
                        last_percentage = percentage;
                    }
                }
            }

            // Include stellar emission, if chosen
            if(tracer[i_det]->considerPointSources())
                calcStellarEmission(i_det, &rand_gen);

            // Show final progress
            cout << "-> Ray tracing dust map(s) (Seq. " << i_det + 1 << ", source: " << sID + 1
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

            if(tracer[i_det]->getSubpixelWarning())
            {
                cout << "WARNING: level of subpixeling (" << param.getMaxSubpixelLvl() << ") might be too low" << endl;
                cout << "if required, increase the maximum level of subpixeling with <max_subpixel_lvl> in the command file" << endl;
            }
            if(tracer[i_det]->getDetectorShape() == DET_PLANE && (grid->getDataID() == GRID_ID_SPH || grid->getDataID() == GRID_ID_CYL))
                cout << "HINT: a 'polar' detector should be used for a spherical or cylindrical grid" << endl;
        }
    }

    // Show that raytracing is finished
    cout << CLR_LINE;
    cout << "- Ray tracing dust map          : done" << endl;

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

    subpixel = tracer[i_det]->getUseSubpixel(cx, cy, subpixel_lvl);

    // If any subpixel traveled through other cells, perform subpixelling
    if(subpixel == false)
    {
        // Init variables
        uint nr_used_wavelengths = tracer[i_det]->getNrOfSpectralBins();

        // Create new photon package
        photon_package pp =
            photon_package(nr_used_wavelengths * (tracer[i_det]->splitDustEmission() ? 4 : 1));

        // Init photon package
        tracer[i_det]->preparePhoton(&pp, cx, cy);

        // Calculate continuum emission along one path
        getDustIntensity(&pp, tmp_source, cx, cy, i_det, subpixel_lvl);

        tracer[i_det]->addToDetector(&pp, i_pix);
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
                                          uint subpixel_lvl)
{
    // Set amount of radiation coming from this pixel
    double subpixel_fraction = pow(4.0, -double(subpixel_lvl));

    // Get chosen wavelength parameter
    uint nr_used_wavelengths = tracer[i_det]->getNrSpectralBins();

    // Update Stokes vectors with emission from background source
    for(uint i_extra = 0; i_extra < tracer[i_det]->getNrExtra(); i_extra++)
    {
        for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
        {
            // Set current index in photon package
            pp->setSpectralID(i_wave + i_extra * nr_used_wavelengths);

            // Set wavelength index in photon package
            double wavelength = tracer[i_det]->getWavelength(i_wave);
            pp->setWavelength(wavelength, dust->getWavelengthID(wavelength));

            // Get emission from background source
            pp->getStokesVector()->set(tmp_source->getStokesVector(pp));
        }
    }

    // Find starting point inside the model and travel through it
    if(grid->findStartingPoint(pp))
        while(grid->next(pp) && tracer[i_det]->isNotAtCenter(pp, cx, cy))
            rayThroughCellDust(pp, i_det, nr_used_wavelengths);

    // Update the multi Stokes vectors for each wavelength
    for(uint i_extra = 0; i_extra < tracer[i_det]->getNrExtra(); i_extra++)
    {
        for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
        {
            // Set current index in photon package
            pp->setSpectralID(i_wave + i_extra * nr_used_wavelengths);

            // Get frequency at background grid position
            double mult = 1e+26 * subpixel_fraction * tracer[i_det]->getDistanceFactor() *
                          pp->getWavelength() / pp->getFrequency();

            // Include foreground extinction if necessary
            mult *= dust->getForegroundExtinction(pp->getWavelength());

            if(pp->getStokesVector()->I() < 0)
                pp->getStokesVector()->setI(0);

            pp->getStokesVector()->multS(mult);
            pp->getStokesVector()->multT(subpixel_fraction);
            pp->getStokesVector()->multSp(subpixel_fraction);
        }
    }
}

void CRadiativeTransfer::rayThroughCellDust(photon_package * pp, uint i_det, uint nr_used_wavelengths)
{
    // Get necessary quantities from current cell
    double dens_gas = grid->getGasNumberDensity(*pp);
    double dens_dust = dust->getNumberDensity(grid, *pp);

    // If the dust density is far too low, skip the current cell
    if(dens_dust >= 1e-200)
    {
        // Get path length through current cell
        double len = pp->getTmpPathLength();

        // Init dust extinction matrix
        Matrix2D dust_extinction_matrix(4, 4);

#if BENCHMARK == CAMPS
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

        myfile << "# nr_of wavelength = " << nr_used_wavelengths << ", corr_factor = " << corr_factor << endl;
        for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
        {
            // Set wavelength index in photon package
            double wavelength = tracer[i_det]->getWavelength(i_wave);
            pp->setWavelength(wavelength, dust->getWavelengthID(wavelength), i_wave);

            StokesVector dust_emissivity;
            dust->calcEmissivityEmi(grid, pp, false, dust_emissivity);
            myfile << wavelength << TAB << corr_factor * wavelength * dust_emissivity.I() << TAB
                   << corr_factor * wavelength * dust_emissivity.Q() << TAB
                   << corr_factor * wavelength * (dust_emissivity.I() + dust_emissivity.Q()) << endl;
        }
        myfile.close();
        break;
#endif

        for(uint i_extra = 0; i_extra < tracer[i_det]->getNrExtra(); i_extra++)
        {
            for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
            {
                // Set current index in photon package
                pp->setSpectralID(i_wave + i_extra * nr_used_wavelengths);

                // Init dust emission matrix
                StokesVector dust_emissivity;

                // Set stokes vector with thermal emission of the dust grains and
                // radiation scattered at dust grains
                dust->calcEmissivityEmi(
                    grid,
                    *pp,
                    stokes_dust_rad_field ? detector_wl_index[i_det] + pp->getDustWavelengthID() : MAX_UINT,
                    i_extra,
                    &dust_emissivity);

                // Get the extinction matrix of the dust grains in the current cell
                dust->calcEmissivityExt(grid, *pp, &dust_extinction_matrix);

                // Init a variable to sum up path lengths until cell is crossed
                double cell_sum = 0.0;

                // First path length is path through cell
                double cell_d_l = len;

                // Save the cell entry position of the photon package
                // (Hint: next(pp) puts the photon onto the border to the next cell)
                Vector3D pos_xyz_cell = pp->getPosition() - (len * pp->getDirection());

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
                        RK_k[k] =
                            dust_extinction_matrix * (scalar_product * cell_d_l + *pp->getStokesVector()) +
                            dust_emissivity;
                    }

                    // Init two temporary Stokes vectors
                    StokesVector stokes_new = *pp->getStokesVector();
                    StokesVector stokes_new2 = *pp->getStokesVector();

                    // Add the result at each Runge-Kutta sub position
                    // to the total Stokes vector (Using Runge-Kutta Fehlberg)
                    for(uint i = 0; i < 6; i++)
                    {
                        stokes_new += RK_k[i] * cell_d_l * RK_b1[i];
                        stokes_new2 += RK_k[i] * cell_d_l * RK_b2[i];
                    }

                    // Delete the Runge-Kutta pointer
                    delete[] RK_k;

                    // Ignore very small (and negative) values
                    if(stokes_new.I() < 1e-200)
                        stokes_new.resetIntensity();
                    if(stokes_new2.I() < 1e-200)
                        stokes_new2.resetIntensity();

                    // Calculate the difference between the results with two
                    // different precisions to see if smaller steps are needed
                    double epsi, dz_new;
                    calcStepWidth(stokes_new, stokes_new2, cell_d_l, &epsi, &dz_new);

                    // If epsi > 1, a smaller step width is needed
                    if(epsi <= 1.0)
                    {
                        // Add the temporary Stokes vector to the total one
                        pp->setStokesVector(stokes_new);

                        // Add optical depth
                        pp->getStokesVector()->addT(-dust_extinction_matrix(0, 0) * cell_d_l);

                        // Add to column density
                        pp->getStokesVector()->addSp(dens_gas * cell_d_l);

                        // Update the position of the photon package
                        pos_xyz_cell += cell_d_l * pp->getDirection();

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

void CRadiativeTransfer::calcStellarEmission(uint i_det, CRandomGenerator * rand_gen)
{
    // Init variables for photon positioning on detector
    int i_pix;

    cout << CLR_LINE;
    cout << "Adding stellar contributions ...              \r";

    // Transport photon to observer for each detector
    for(uint s = 0; s < sources_mc.size(); s++)
    {
        // Ignore dust as Monte-Carle radiation source
        if(sources_mc[s]->getID() == SRC_DUST)
            continue;

        // Get chosen wavelength parameter
        uint nr_used_wavelengths = tracer[i_det]->getNrOfSpectralBins();

        // Create new photon package
        photon_package * pp = new photon_package(nr_used_wavelengths);

        // Get position of source
        Vector3D source_pos = sources_mc[s]->getPosition();

        cout << "Processing source: " << s << " of " << sources_mc.size()
             << "                              \r";

        // Update Stokes vectors with emission from background source
        for(uint i_wave = 0; i_wave < nr_used_wavelengths; i_wave++)
        {
            // Set current index in photon package
            pp->setSpectralID(i_wave);

            // Get frequency at background grid position
            double wavelength = tracer[i_det]->getWavelength(i_wave);
            pp->setWavelength(wavelength, dust->getWavelengthID(wavelength));

            double mult = 1e+26 * con_c / (pp->getFrequency() * pp->getFrequency());

            // Set direction of the photon package to the observer
            tracer[i_det]->preparePhotonWithPosition(pp, source_pos, i_pix);

            // Launch photon package
            sources_mc[s]->createDirectRay(pp, rand_gen);

            // Position the photon inside the grid
            grid->positionPhotonInGrid(pp);

            // Init a variable to save the optical depth
            double tau_obs = 0;

            // Transport photon package through model to obtain the optical depth
            while(grid->next(pp))
            {
                // Get necessary information about the current cell
                double len = pp->getTmpPathLength();
                double dens = dust->getNumberDensity(grid, *pp);
                double Cext = dust->getCextMean(grid, *pp);

                // Increase the optical depth by the current cell
                tau_obs += Cext * len * dens;
            }
            mult *= exp(-tau_obs) * tracer[i_det]->getDistanceFactor(source_pos);

            // Update the photon package with the multi Stokes vectors
            pp->getStokesVector()->multS(mult);
        }

        tracer[i_det]->addToDetector(pp, i_pix, true);
    }

    cout << CLR_LINE;
}


// -------------------------------------------
// ------ Calculation of OPIATE transfer -------
// -------------------------------------------

bool CRadiativeTransfer::calcOPIATEMapsViaRaytracing(parameters& param)
{
    // Create a list for all gas species
    dlist op_ray_detector_list = param.getOPIATERayDetectors();

    // Perform radiative transfer for each chosen gas species
    for(uint i_det = start; i_det <= stop; i_det++)
    {
        // Get total number of pixel
        uint per_max = tracer[i_det]->getNpix();
        string spec_name=param.getOpiateSpec(i_det);

        if(!op->findIndexByName(spec_name))
            return false;

        // Get BG source
        uint sID = tracer[i_det]->getSourceIndex()-1;
        CSourceBasic * tmp_source;
        tmp_source = sources_ray[sID];

        // Calculate the line broadening for each cell
        op->calcLineBroadening(grid);

        // Show progress of the current sequence and gas species
        cout << CLR_LINE;
        cout << "-> Ray tracing OPIATE map(s): species " << i_det + 1 << " of " << stop + 1
             << " : 0.0 [%]       \r" << flush;

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
# pragma omp atomic update
            per_counter++;

            // Calculate percentage of total progress per source
            float percentage = 100.0 * float(per_counter) / float(per_max);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
#pragma omp critical
                {
                    cout << "-> Ray tracing OPIATE map(s): species " << i_det + 1 << " of " << stop + 1
                         << " : " << percentage << " [%]       \r" << flush;
                    last_percentage = percentage;
                }
            }

            // Get radiative transfer results for one pixel/ray
            getOPIATEPixelIntensity(tmp_source, cx, cy, 0, 0, i_det, uint(0), i_pix);
            //getLinePixelIntensity(tmp_source, cx, cy, i_species, i_trans, i_det, uint(0), i_pix);
        }

        // post-process raytracing simulation
        //if(!tracer[i_det]->postProcessing())
        //    return false;

        if(!tracer[i_det]->writeOpiateResults(op))
            return false;

        if(tracer[i_det]->getSubpixelWarning())
        {
            cout << "WARNING: level of subpixeling (" << param.getMaxSubpixelLvl() << ") might be too low" << endl;
            cout << "if required, increase the maximum level of subpixeling with <max_subpixel_lvl> in the command file" << endl;
        }
        if(tracer[i_det]->getDetectorShape() == DET_PLANE && (grid->getDataID() == GRID_ID_SPH || grid->getDataID() == GRID_ID_CYL))
            cout << "HINT: a 'polar' detector should be used for a spherical or cylindrical grid" << endl;
    }

    cout << CLR_LINE;
    cout << "- Ray tracing channel map(s)    : done" << endl;

    return true;
}

// -------------------------------------------
// ------ Calculation of Line transfer -------
// -------------------------------------------

bool CRadiativeTransfer::calcChMapsViaRaytracing(parameters & param)
{
    // Create a list for all gas species
    maplist line_ray_detector_list = param.getLineRayDetectors();
    maplist::iterator it;

    // Index to the current detector/tracer
    uint i_det = 0;

    // Select source as specified in the detector
    uint sID = tracer[i_det]->getSourceIndex();
    CSourceBasic * tmp_source;
    tmp_source = sources_ray[sID];

    // Get total number of pixel
    uint per_max = tracer[i_det]->getNpix();

    // Perform radiative transfer for each chosen gas species
    for(it = line_ray_detector_list.begin(); it != line_ray_detector_list.end(); ++it)
    {
        // Get ID of the current gas species
        uint i_species = it->first;

        // Get number of spectral line transitions that have to be simulated
        uint nr_of_spectral_lines = gas->getNrOfSpectralLines(i_species);

        // If no spectral line transitions are chosen for the current gas species, skip
        if(nr_of_spectral_lines == 0)
            continue;

        // Get the detectors of the current gas species
        dlist line_ray_detectors = it->second;

        // Start and Stop define the index of the first and last gas species
        if(i_species < start || i_species > stop)
            continue;

        // Calculate the level populations for each cell
        bool lvl_pop_error = gas->getLevelPopType(i_species) == POP_MC
                                 ? !calcMonteCarloLvlPopulation(i_species, param.getMCLvlPopSeed())
                                 : !gas->calcLevelPopulation(grid, i_species);
        if(lvl_pop_error)
        {
            cout << "\nERROR: Level population cannot be calculated!";
            return false;
        }

        // Calculate the line broadening for each cell
        gas->calcLineBroadening(grid, i_species);

        // Perform radiative transfer for each chosen spectral line transition
        for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
        {
            // Get transition index to spectral line
            uint i_trans = gas->getTransitionFromSpectralLine(i_species, i_line);

            // Get number of velocity channels.
            uint nr_velocity_channels = tracer[i_det]->getNrSpectralBins();
            if(gas->isLineZeemanSplit(i_species, i_line) && nr_velocity_channels < 6)
            {
                cout << "\nERROR: The magnetic field information requires at least 6 "
                        "channels\n"
                     << "    for simulations with Zeeman splitting" << endl;
                return false;
            }

            // Show progress of the current sequence and gas species
            cout << CLR_LINE;
            cout << "-> Channel map(s): gas species " << i_species + 1 << " of " << stop + 1 << ", line "
                 << i_line + 1 << " of " << nr_of_spectral_lines << ": 0.0[%]  \r" << flush;

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
#pragma omp atomic update
                per_counter++;

                // Calculate percentage of total progress per source
                float percentage = 100.0 * float(per_counter) / float(per_max);

                // Show only new percentage number if it changed
                if((percentage - last_percentage) > PERCENTAGE_STEP)
                {
#pragma omp critical
                    {
                        cout << "-> Channel map(s): gas species " << i_species + 1 << " of " << stop + 1
                             << ", line " << i_line + 1 << " of " << nr_of_spectral_lines << ": "
                             << percentage << " [%]      \r" << flush;
                        last_percentage = percentage;
                    }
                }

                // Get radiative transfer results for one pixel/ray
                getLinePixelIntensity(tmp_source, cx, cy, i_species, i_trans, i_det, uint(0), i_pix);
            }

            // post-process raytracing simulation
            if(!tracer[i_det]->postProcessing())
                return false;

            if(!tracer[i_det]->writeLineResults(gas, i_species, i_line))
                return false;

            if(tracer[i_det]->getSubpixelWarning())
            {
                cout << "WARNING: level of subpixeling (" << param.getMaxSubpixelLvl() << ") might be too low" << endl;
                cout << "if required, increase the maximum level of subpixeling with <max_subpixel_lvl> in the command file" << endl;
            }
            if(tracer[i_det]->getDetectorShape() == DET_PLANE && (grid->getDataID() == GRID_ID_SPH || grid->getDataID() == GRID_ID_CYL))
                cout << "HINT: a 'polar' detector should be used for a spherical or cylindrical grid" << endl;

            // Increment index for line RT
            i_det++;
        }
    }

    cout << CLR_LINE;
    cout << "- Ray tracing channel map(s)    : done" << endl;

    return true;
}


void CRadiativeTransfer::getOPIATEPixelIntensity(CSourceBasic * tmp_source,
                                               double cx,
                                               double cy,
                                               uint i_species,
                                               uint i_trans,
                                               uint i_det,
                                               uint subpixel_lvl,
                                               int i_pix)
{
    bool subpixel = false;

    subpixel = tracer[i_det]->getUseSubpixel(cx, cy, subpixel_lvl);

    // If any subpixel traveled through other cells, perform subpixelling
    if(subpixel == false)
    {
        // Get rest frequency of current transition
        double trans_frequency = op->getCurrentFrequency();

        // Create new photon package
        photon_package pp = photon_package(trans_frequency,
                                           dust->getWavelengthID(con_c / trans_frequency),
                                           tracer[i_det]->getNrSpectralBins());

        // Init photon package
        tracer[i_det]->preparePhoton(&pp, cx, cy);

        // Calculate line emission along one path
        getOPIATEIntensity(&pp, tmp_source, cx, cy, i_species, i_trans, i_det, subpixel_lvl);

        tracer[i_det]->addToDetector(&pp, i_pix);
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
                getOPIATEPixelIntensity(tmp_source, tmp_cx, tmp_cy, i_species, i_trans, i_det, (subpixel_lvl + 1), i_pix);
            }
        }
    }
}

void CRadiativeTransfer::getOPIATEIntensity(photon_package * pp,
                                          CSourceBasic * tmp_source,
                                          double cx,
                                          double cy,
                                          uint i_species,
                                          uint i_trans,
                                          uint i_det,
                                          uint subpixel_lvl)
{
    // Set amount of radiation coming from this pixel
    double subpixel_fraction = pow(4.0, -double(subpixel_lvl));

    // Number of velocity channels from tracer
    uint nr_velocity_channels = tracer[i_det]->getNrSpectralBins();

    // Set background radiation field
    for(uint vch = 0; vch < nr_velocity_channels; vch++)
    {
        // Set current index in photon package
        pp->setSpectralID(vch);

        // Set velocity of the different velocity channels
        pp->setVelocity(tracer[i_det]->getVelocityChannel(vch));

        // Set background emission
        StokesVector st=tmp_source->getStokesVector(pp);
        double lam=pp->getWavelength();
        double freq=pp->getFrequency();

        pp->setStokesVector(st * lam / freq);
    }

    tracer[i_det]->preparePhoton(pp, cx, cy);
    if(grid->findStartingPoint(pp))
    {
        // Set/Get Necessary information about velocity field interpolation
        VelFieldInterp vel_field_interp;
        vel_field_interp.start_pos = pp->getPosition();

        // Precalculate the velocity interpolation
        preCalcVelocityInterp(grid, *pp, &vel_field_interp);

        // Transport the photon package through the model
        while(grid->next(pp) && tracer[i_det]->isNotAtCenter(pp, cx, cy))
        {
            rayThroughCellOPIATE(pp, i_species, i_trans, i_det, nr_velocity_channels, vel_field_interp);
        }
    }
    // Update the multi Stokes vectors for each velocity channel
    for(uint vch = 0; vch < nr_velocity_channels; vch++)
    {
        // Set current index in photon package
        pp->setSpectralID(vch);

        // Convert W/m2/Hz/sr to Jy
        double mult = 1e+26 * subpixel_fraction * tracer[i_det]->getDistanceFactor();

        // Include foreground extinction if necessary
        mult *= dust->getForegroundExtinction(pp->getWavelength());

        if(pp->getStokesVector()->I() < 0)
            pp->getStokesVector()->setI(0);

        pp->getStokesVector()->multS(mult);
        pp->getStokesVector()->multT(subpixel_fraction);
        pp->getStokesVector()->multSp(subpixel_fraction);
    }
}


void CRadiativeTransfer::getLinePixelIntensity(CSourceBasic * tmp_source,
                                               double cx,
                                               double cy,
                                               uint i_species,
                                               uint i_trans,
                                               uint i_det,
                                               uint subpixel_lvl,
                                               int i_pix)
{
    bool subpixel = false;

    subpixel = tracer[i_det]->getUseSubpixel(cx, cy, subpixel_lvl);

    // If any subpixel traveled through other cells, perform subpixelling
    if(subpixel == false)
    {
        // Get rest frequency of current transition
        double trans_frequency = gas->getTransitionFrequency(i_species, i_trans);

        // Create new photon package
        photon_package pp = photon_package(trans_frequency,
                                           dust->getWavelengthID(con_c / trans_frequency),
                                           tracer[i_det]->getNrSpectralBins());

        // Init photon package
        tracer[i_det]->preparePhoton(&pp, cx, cy);

        // Calculate line emission along one path
        getLineIntensity(&pp, tmp_source, cx, cy, i_species, i_trans, i_det, subpixel_lvl);

        tracer[i_det]->addToDetector(&pp, i_pix);
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
                    tmp_source, tmp_cx, tmp_cy, i_species, i_trans, i_det, (subpixel_lvl + 1), i_pix);
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
                                          uint i_species,
                                          uint i_trans,
                                          uint i_det,
                                          uint subpixel_lvl)
{
    // Set amount of radiation coming from this pixel
    double subpixel_fraction = pow(4.0, -double(subpixel_lvl));

    // Number of velocity channels from tracer
    uint nr_velocity_channels = tracer[i_det]->getNrSpectralBins();

    // Set background radiation field
    for(uint vch = 0; vch < nr_velocity_channels; vch++)
    {
        // Set current index in photon package
        pp->setSpectralID(vch);

        // Set velocity of the different velocity channels
        pp->setVelocity(tracer[i_det]->getVelocityChannel(vch));

        // Set background emission
        pp->setStokesVector(tmp_source->getStokesVector(pp) * pp->getWavelength() / pp->getFrequency());
    }

    tracer[i_det]->preparePhoton(pp, cx, cy);
    if(grid->findStartingPoint(pp))
    {
        // Set/Get Necessary information about velocity field interpolation
        VelFieldInterp vel_field_interp;
        vel_field_interp.start_pos = pp->getPosition();

        // Precalculate the velocity interpolation
        preCalcVelocityInterp(grid, *pp, &vel_field_interp);

        // Transport the photon package through the model
        while(grid->next(pp) && tracer[i_det]->isNotAtCenter(pp, cx, cy))
        {
            rayThroughCellLine(pp, i_species, i_trans, i_det, nr_velocity_channels, vel_field_interp);
        }
    }
    // Update the multi Stokes vectors for each velocity channel
    for(uint vch = 0; vch < nr_velocity_channels; vch++)
    {
        // Set current index in photon package
        pp->setSpectralID(vch);

        // Convert W/m2/Hz/sr to Jy
        double mult = 1e+26 * subpixel_fraction * tracer[i_det]->getDistanceFactor();

        // Include foreground extinction if necessary
        mult *= dust->getForegroundExtinction(pp->getWavelength());

        if(pp->getStokesVector()->I() < 0)
            pp->getStokesVector()->setI(0);

        pp->getStokesVector()->multS(mult);
        pp->getStokesVector()->multT(subpixel_fraction);
        pp->getStokesVector()->multSp(subpixel_fraction);
    }
}


void CRadiativeTransfer::rayThroughCellOPIATE(photon_package * pp,
                                            uint i_species,
                                            uint i_trans,
                                            uint i_det,
                                            uint nr_velocity_channels,
                                            const VelFieldInterp & vel_field_interp)
{
    // Get gas species density from grid
    double dens_gas = grid->getGasDensity(*pp);

    // Perform radiative transfer only if the density of the current species are not negligible
    if(dens_gas > 1e-200)
    {
        // Init matrix for absorption and dust emissivity
        Matrix2D total_absorption_matrix(4, 4);
        StokesVector dust_emi_and_ext;

        // Get extra information about the magnetic field and ine broadening
        MagFieldInfo mag_field_info;
        LineBroadening line_broadening;

        // Set only gauss_a if not zeeman split
        line_broadening.gauss_a = grid->getGaussA(*pp);

        // Get the path length through the current cell
        double len = pp->getTmpPathLength();

        // Get necessary quantities from the current cell
        // double dens_gas = grid->getGasNumberDensity(*pp);

        // Calculate the emission of the dust grains
        dust->calcEmissivityHz(grid, *pp, &dust_emi_and_ext);

        // Perform radiative transfer for each velocity channel separately
        for(uint vch = 0; vch < nr_velocity_channels; vch++)
        {
            // Set current index in photon package
            pp->setSpectralID(vch);

            // Init variables
            double cell_sum = 0, cell_d_l = len;
            ullong kill_counter = 0;

            // Get direction and entry position of the current cell
            Vector3D pos_xyz_cell = pp->getPosition() - (len * pp->getDirection());

            // Make sub steps until cell is completely crossed
            // If the error of a sub step is too high, make the step smaller
            while(cell_sum < len)
            {
                // Increase the kill counter
                kill_counter++;

                // If too many sub steps are needed, kill the photon
                if(kill_counter >= MAX_SOLVER_STEPS)
                {
                    cout << "\nWARNING: Solver steps > " << MAX_SOLVER_STEPS << ". Too many steps!" << endl;
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
                    double rel_velocity =
                        pp->getVelocity() -
                        op->getProjCellVelocityInterp(pos_xyz_cell + cell_d_l * pp->getDirection() * RK_c[k],
                                                       pp->getDirection(),
                                                       vel_field_interp);

                    Vector3D obs_vel = tracer[i_det]->getObserverVelocity();
                    if(obs_vel.length() > 0)
                        rel_velocity += obs_vel * pp->getDirection();

                    // Init emission
                    StokesVector total_emission;

                    // Get line emissivity (also combined Zeeman lines)
                    op->getMatrices(grid, pp, i_species, i_trans, rel_velocity, line_broadening, mag_field_info, &total_emission, &total_absorption_matrix);

                    // Combine the Stokes vectors from gas and dust for emission
                    //total_emission *= dens_species;

                    total_emission += dust_emi_and_ext;

                    // and extinction
                    //total_absorption_matrix *= dens_species;

                    if(dust_emi_and_ext.T() != 0)
                        for(uint i = 0; i < 4; i++)
                            total_absorption_matrix(i, i) += dust_emi_and_ext.T();
                    total_absorption_matrix *= -1;

                    // Init scalar product
                    StokesVector scalar_product;

                    // Calculate multiplication between Runge-Kutta parameter
                    for(uint i = 0; i <= k; i++)
                        scalar_product += (RK_k[i] * RK_a(i, k));

                    // Calculate new Runge-Kutta parameters as the result of the
                    // radiative transfer equation at the Runge-Kutta sub
                    // positions
                    RK_k[k] = total_absorption_matrix * (scalar_product * cell_d_l + *pp->getStokesVector()) +
                              total_emission;
                }
                // Init two temporary Stokes vectors
                StokesVector stokes_new = *pp->getStokesVector();
                StokesVector stokes_new2 = *pp->getStokesVector();

                for(uint i = 0; i < 6; i++)
                {
                    stokes_new += RK_k[i] * cell_d_l * RK_b1[i];
                    stokes_new2 += RK_k[i] * cell_d_l * RK_b2[i];
                }

                // Delete the Runge-Kutta pointer
                delete[] RK_k;

                // Ignore very small values
                if(abs(stokes_new.I()) < 1e-200)
                    stokes_new.resetIntensity();
                if(abs(stokes_new2.I()) < 1e-200)
                    stokes_new2.resetIntensity();

                // Calculate the difference between the results with two
                // different precisions to see if smaller steps are needed
                double epsi, dz_new;
                calcStepWidth(stokes_new, stokes_new2, cell_d_l, &epsi, &dz_new);

                // Is a smaller step width needed
                if(epsi <= 1.0)
                {
                    // Backup old intensity
                    double old_stokes = pp->getStokesVector()->I();

                    // Stokes_new is the current flux of this line-of-sight
                    pp->setStokesVector(stokes_new);

                    // Columns density
                    //double column_density = dens_gas * cell_d_l;
                    double column_density = dens_gas * cell_d_l;

                    if(vch == 0)
                    {
                        // Column density of the total gas
                        pp->getStokesVector()->addSp(column_density);
                    }

                    // Save the optical depth of each velocity channel, if
                    // magnetic field analysis is not chosen
                    pp->getStokesVector()->addT(-total_absorption_matrix(0, 0) * cell_d_l);

                    // Update the position of the photon package
                    pos_xyz_cell += cell_d_l * pp->getDirection();

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

void CRadiativeTransfer::rayThroughCellLine(photon_package * pp,
                                            uint i_species,
                                            uint i_trans,
                                            uint i_det,
                                            uint nr_velocity_channels,
                                            const VelFieldInterp & vel_field_interp)
{
    // Get gas species density from grid
    double dens_species = gas->getNumberDensity(grid, *pp, i_species);

    // Perform radiative transfer only if the temperature of the current species
    // are not negligible
    if(dens_species > 1e-200)
    {
        // Init matrix for absorption and dust emissivity
        Matrix2D total_absorption_matrix(4, 4);
        StokesVector dust_emi_and_ext;

        // Get extra information about the magnetic field and ine broadening
        MagFieldInfo mag_field_info;
        LineBroadening line_broadening;
        uint i_zeeman = gas->getZeemanSplitIndex(i_species, i_trans);
        if(i_zeeman != MAX_UINT)
        {
            grid->getMagFieldInfo(*pp, &mag_field_info);
            grid->getLineBroadening(*pp, i_zeeman, &line_broadening);
        }
        else
        {
            // Set only gauss_a if not zeeman split
            line_broadening.gauss_a = grid->getGaussA(*pp);
        }

        // Get the path length through the current cell
        double len = pp->getTmpPathLength();

        // Get necessary quantities from the current cell
        double dens_gas = grid->getGasNumberDensity(*pp);

        // Calculate the emission of the dust grains
        dust->calcEmissivityHz(grid, *pp, &dust_emi_and_ext);

        // Perform radiative transfer for each velocity channel separately
        for(uint vch = 0; vch < nr_velocity_channels; vch++)
        {
            // Set current index in photon package
            pp->setSpectralID(vch);

            // Init variables
            double cell_sum = 0, cell_d_l = len;
            ullong kill_counter = 0;

            // Get direction and entry position of the current cell
            Vector3D pos_xyz_cell = pp->getPosition() - (len * pp->getDirection());

            // Make sub steps until cell is completely crossed
            // If the error of a sub step is too high, make the step smaller
            while(cell_sum < len)
            {
                // Increase the kill counter
                kill_counter++;

                // If too many sub steps are needed, kill the photon
                if(kill_counter >= MAX_SOLVER_STEPS)
                {
                    cout << "\nWARNING: Solver steps > " << MAX_SOLVER_STEPS << ". Too many steps!" << endl;
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
                    double rel_velocity =
                        pp->getVelocity() -
                        gas->getProjCellVelocityInterp(pos_xyz_cell + cell_d_l * pp->getDirection() * RK_c[k],
                                                       pp->getDirection(),
                                                       vel_field_interp);

                    Vector3D obs_vel = tracer[i_det]->getObserverVelocity();
                    if(obs_vel.length() > 0)
                        rel_velocity += obs_vel * pp->getDirection();

                    // Init emission
                    StokesVector total_emission;

                    // Get line emissivity (also combined Zeeman lines)
                    gas->calcEmissivity(grid,
                                        *pp,
                                        i_species,
                                        i_trans,
                                        rel_velocity,
                                        line_broadening,
                                        mag_field_info,
                                        &total_emission,
                                        &total_absorption_matrix);

                    // Combine the Stokes vectors from gas and dust for emission
                    total_emission *= dens_species;
                    total_emission += dust_emi_and_ext;
                    // and extinction
                    total_absorption_matrix *= dens_species;
                    if(dust_emi_and_ext.T() != 0)
                        for(uint i = 0; i < 4; i++)
                            total_absorption_matrix(i, i) += dust_emi_and_ext.T();
                    total_absorption_matrix *= -1;

                    // Init scalar product
                    StokesVector scalar_product;

                    // Calculate multiplication between Runge-Kutta parameter
                    for(uint i = 0; i <= k; i++)
                        scalar_product += (RK_k[i] * RK_a(i, k));

                    // Calculate new Runge-Kutta parameters as the result of the
                    // radiative transfer equation at the Runge-Kutta sub
                    // positions
                    RK_k[k] = total_absorption_matrix * (scalar_product * cell_d_l + *pp->getStokesVector()) +
                              total_emission;
                }
                // Init two temporary Stokes vectors
                StokesVector stokes_new = *pp->getStokesVector();
                StokesVector stokes_new2 = *pp->getStokesVector();

                for(uint i = 0; i < 6; i++)
                {
                    stokes_new += RK_k[i] * cell_d_l * RK_b1[i];
                    stokes_new2 += RK_k[i] * cell_d_l * RK_b2[i];
                }

                // Delete the Runge-Kutta pointer
                delete[] RK_k;

                // Ignore very small values
                if(abs(stokes_new.I()) < 1e-200)
                    stokes_new.resetIntensity();
                if(abs(stokes_new2.I()) < 1e-200)
                    stokes_new2.resetIntensity();

                // Calculate the difference between the results with two
                // different precisions to see if smaller steps are needed
                double epsi, dz_new;
                calcStepWidth(stokes_new, stokes_new2, cell_d_l, &epsi, &dz_new);

                // Is a smaller step width needed
                if(epsi <= 1.0)
                {
                    // Backup old intensity
                    double old_stokes = pp->getStokesVector()->I();

                    // Stokes_new is the current flux of this line-of-sight
                    pp->setStokesVector(stokes_new);

                    // Columns density
                    double species_column_density = dens_species * cell_d_l;
                    double gas_column_density = dens_gas * cell_d_l;

                    if(gas->isTransZeemanSplit(i_species, i_trans))
                    {
                        // Total increase of the intensity along the line-of-sight
                        double column_flux = (stokes_new.I() - old_stokes);

                        // Total magnetic field strength of the current cell
                        double mag_strength = mag_field_info.mag_field.length();

                        // LOS magnetic field strength of the current cell
                        double los_mag_strength = (pp->getDirection() * mag_field_info.mag_field);

                        // Magnetic field strength in the line-of-sight direction
                        // weighted with the intensity increase of the current
                        // cell
                        double column_int_mag_field_los = los_mag_strength * column_flux;

                        // Total magnetic field strength weighted with the
                        // intensity increase of the current cell
                        double column_int_mag_field = mag_strength * column_flux;

                        // Intensity weighted LOS magnetic field
                        pp->getStokesVector(0)->addSp(column_int_mag_field_los);

                        // Intensity weighted total magnetic field
                        pp->getStokesVector(1)->addSp(column_int_mag_field);

                        // Flux component for weighting
                        pp->getStokesVector(2)->addSp(column_flux);

                        if(vch == 0)
                        {
                            // Magnetic field strength in the line-of-sight
                            // direction weighted with the gas density of the
                            // current cell
                            double column_dens_mag_field_los = los_mag_strength * gas_column_density;

                            // Total magnetic field strength weighted with the
                            // gas density of the current cell
                            double column_dens_mag_field = mag_strength * gas_column_density;

                            // Density weighted LOS magnetic field
                            pp->getStokesVector(3)->addSp(column_dens_mag_field_los);

                            // Density weighted magnetic field
                            pp->getStokesVector(4)->addSp(column_dens_mag_field);

                            // Column density of the total gas
                            pp->getStokesVector(5)->addSp(gas_column_density);

                            // Column density of the species
                            pp->getStokesVector(6)->addSp(species_column_density);
                        }
                    }
                    else if(vch == 0)
                    {
                        // Column density of the total gas

                        pp->getStokesVector(0)->addSp(gas_column_density);

                        // Column density of the species
                        pp->getStokesVector(1)->addSp(species_column_density);
                    }

                    // Save the optical depth of each velocity channel, if
                    // magnetic field analysis is not chosen
                    pp->getStokesVector()->addT(-total_absorption_matrix(0, 0) * cell_d_l);

                    // Update the position of the photon package
                    pos_xyz_cell += cell_d_l * pp->getDirection();

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

void CRadiativeTransfer::preCalcVelocityInterp(CGridBasic * grid,
                                               const photon_package & pp,
                                               VelFieldInterp * vel_field_interp)
{
    vel_field_interp->zero_vel_field = true;
    vel_field_interp->start_pos = pp.getPosition();

    if(grid->hasVelocityField() == true)
    {
        photon_package pp_interp = photon_package();
        pp_interp.setPosition(pp.getPosition());
        pp_interp.setDirection(pp.getDirection());
        if(grid->findStartingPoint(&pp_interp))
        {
            Vector3D pos_in_grid_0 = pp_interp.getPosition();
            double spline_x_old = 0;
            while(grid->next(&pp_interp))
            {
                Vector3D pos_xyz_cell =
                    pp_interp.getPosition() - (pp_interp.getTmpPathLength() * pp_interp.getDirection());
                Vector3D rel_pos = pos_xyz_cell - pos_in_grid_0;

                const cell_basic & tmp_cell_pos = *pp_interp.getPositionCell();
                Vector3D cell_center = grid->getCenter(tmp_cell_pos);

                double length_on_line =
                    CMathFunctions::getClosestLinePoint(pos_xyz_cell, pp_interp.getPosition(), cell_center);
                double spline_x = length_on_line + rel_pos.length();

                if((spline_x_old - spline_x) != 0.0)
                {
                    double spline_y = grid->getVelocityField(pp_interp) * pp_interp.getDirection();

                    spline_x_old = spline_x;

                    if(spline_y != 0)
                        vel_field_interp->zero_vel_field = false;

                    vel_field_interp->vel_field.setDynValue(spline_x, spline_y);
                }
            }
        }

        vel_field_interp->vel_field.createDynSpline();
    }
}