#include "GasSpecies.h"
#include "CommandParser.h"
#include "Grid.h"
#include "MathFunctions.h"
#include "typedefs.h"

#define TRANS_SIGMA_P +1
#define TRANS_PI 0
#define TRANS_SIGMA_M -1

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

bool CGasSpecies::calcLTE(CGridBasic * grid, bool full)
{
    uint nr_of_total_energy_levels = getNrOfTotalEnergyLevels();
    uint nr_of_spectral_lines = getNrOfSpectralLines();
    float last_percentage = 0;
    long cell_count = 0;
    long max_cells = grid->getMaxDataCells();

    cout << CLR_LINE;

#pragma omp parallel for schedule(dynamic)
    for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
    {
        cell_basic * cell = grid->getCellFromIndex(i_cell);
        // Calculate percentage of total progress per source
        float percentage = 100 * float(cell_count) / float(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
#pragma omp critical
            {
                cout << "-> Calculating LTE level population  : " << percentage << " [%]               \r";
                last_percentage = percentage;
            }
        }

        cell_count++;
        double ** tmp_lvl_pop = new double *[nr_of_energy_level];
        for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
        {
            tmp_lvl_pop[i_lvl] = new double[getNrOfSublevel(i_lvl)];
        }

        double temp_gas = grid->getGasTemperature(cell);
        if(temp_gas == 0)
        {
            for(uint i_lvl = 0; i_lvl < nr_of_total_energy_levels; i_lvl++)
            {
                for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
                    tmp_lvl_pop[i_lvl][i_sublvl] = 0;
            }
            tmp_lvl_pop[0][0] = 1;
        }
        else
        {
            double sum = 0;
            uint tmp_offset = 0;
            for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
            {
                // Calculate LTE level pop also for Zeeman sublevel
                for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
                {
                    tmp_lvl_pop[i_lvl][i_sublvl] =
                        getGLevel(i_lvl) *
                        exp(-con_h * getEnergylevel(i_lvl) * con_c * 100.0 / (con_kB * temp_gas));
                    sum += tmp_lvl_pop[i_lvl][i_sublvl];
                }

                tmp_offset += getNrOfSublevel(i_lvl);
            }

            for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
            {
                for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
                {
                    tmp_lvl_pop[i_lvl][i_sublvl] /= sum;
                }
            }
        }

        if(full)
        {
            for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
            {
                for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
                {
                    grid->setLvlPop(cell, i_lvl, i_sublvl, tmp_lvl_pop[i_lvl][i_sublvl]);
                }
            }
        }
        else
        {
            for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
            {
                uint i_trans = getTransitionFromSpectralLine(i_line);

                uint i_lvl_l = getLowerEnergyLevel(i_trans);
                for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl_l); i_sublvl++)
                {
                    grid->setLvlPopLower(cell, i_line, i_sublvl, tmp_lvl_pop[i_lvl_l][i_sublvl]);
                }

                uint i_lvl_u = getUpperEnergyLevel(i_trans);
                for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl_u); i_sublvl++)
                {
                    grid->setLvlPopUpper(cell, i_line, i_sublvl, tmp_lvl_pop[i_lvl_u][i_sublvl]);
                }
            }
        }

        delete[] tmp_lvl_pop;
    }
    return true;
}

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

bool CGasSpecies::calcFEP(CGridBasic * grid, bool full)
{
    uint nr_of_total_energy_levels = getNrOfTotalEnergyLevels();
    uint nr_of_total_transitions = getNrOfTotalTransitions();
    uint nr_of_spectral_lines = getNrOfSpectralLines();
    float last_percentage = 0;
    long max_cells = grid->getMaxDataCells();
    double * J_mid = new double[nr_of_total_transitions];
    bool no_error = true;

    cout << CLR_LINE;

    for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
    {
        uint i_lvl_l = getLowerEnergyLevel(i_trans);
        uint i_lvl_u = getUpperEnergyLevel(i_trans);

        for(uint i_sublvl_u = 0; i_sublvl_u < getNrOfSublevel(i_lvl_u); i_sublvl_u++)
        {
            for(uint i_sublvl_l = 0; i_sublvl_l < getNrOfSublevel(i_lvl_l); i_sublvl_l++)
            {
                uint i_tmp_trans = getUniqueTransitionIndex(i_trans, i_sublvl_u, i_sublvl_l);
                J_mid[i_tmp_trans] = CMathFunctions::planck_hz(getTransitionFrequency(i_trans), 2.75);
            }
        }
    }

#pragma omp parallel for schedule(dynamic)
    for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
    {
        cell_basic * cell = grid->getCellFromIndex(i_cell);
        double * tmp_lvl_pop = new double[nr_of_total_energy_levels];
        double gas_number_density = grid->getGasNumberDensity(cell);

        // Calculate percentage of total progress per source
        float percentage = 100 * float(i_cell) / float(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
#pragma omp critical
            {
                cout << "-> Calculating FEP level population  : " << last_percentage
                     << " [%]                       \r";
                last_percentage = percentage;
            }
        }

        if(gas_number_density > 1e-200)
        {
            Matrix2D A(nr_of_total_energy_levels, nr_of_total_energy_levels);
            double * b = new double[nr_of_total_energy_levels];

            double *** final_col_para = calcCollisionParameter(grid, cell);
            createMatrix(J_mid, A, b, final_col_para);
            CMathFunctions::gauss(A, b, tmp_lvl_pop, nr_of_total_energy_levels);

            delete[] b;
            for(uint i_col_partner = 0; i_col_partner < nr_of_col_partner; i_col_partner++)
            {
                for(uint i_col_transition = 0; i_col_transition < nr_of_col_transition[i_col_partner];
                    i_col_transition++)
                {
                    delete[] final_col_para[i_col_partner][i_col_transition];
                }
                delete[] final_col_para[i_col_partner];
            }
            delete[] final_col_para;

            double sum_p = 0;
            for(uint i_lvl_tot = 0; i_lvl_tot < nr_of_total_energy_levels; i_lvl_tot++)
            {
                if(tmp_lvl_pop[i_lvl_tot] >= 0)
                    sum_p += tmp_lvl_pop[i_lvl_tot];
                else
                {
                    cout << "WARNING: Level population element not greater than zero! Level = " << i_lvl_tot
                         << ", Level pop = " << tmp_lvl_pop[i_lvl_tot] << endl;
                    no_error = false;
                }
            }
            for(uint i_lvl_tot = 0; i_lvl_tot < nr_of_total_energy_levels; i_lvl_tot++)
                tmp_lvl_pop[i_lvl_tot] /= sum_p;
        }
        else
        {
            for(uint i_lvl_tot = 0; i_lvl_tot < nr_of_total_energy_levels; i_lvl_tot++)
                tmp_lvl_pop[i_lvl_tot] = 0;
            tmp_lvl_pop[0] = 1;
        }

        if(full)
        {
            for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
            {
                for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
                {
                    grid->setLvlPop(cell, i_lvl, i_sublvl, tmp_lvl_pop[getUniqueLevelIndex(i_lvl, i_sublvl)]);
                }
            }
        }
        else
        {
            for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
            {
                for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
                {
                    // Get indices
                    uint i_trans = getTransitionFromSpectralLine(i_line);
                    uint i_lvl_l = getLowerEnergyLevel(i_trans);
                    uint i_lvl_u = getUpperEnergyLevel(i_trans);

                    if(i_lvl == i_lvl_l)
                    {
                        for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
                            grid->setLvlPopLower(
                                cell, i_line, i_sublvl, tmp_lvl_pop[getUniqueLevelIndex(i_lvl, i_sublvl)]);
                    }
                    else if(i_lvl == i_lvl_u)
                    {
                        for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
                            grid->setLvlPopUpper(
                                cell, i_line, i_sublvl, tmp_lvl_pop[getUniqueLevelIndex(i_lvl, i_sublvl)]);
                    }
                }
            }
        }

        delete[] tmp_lvl_pop;
    }

    delete[] J_mid;
    return no_error;
}

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de
// and based on Pavlyuchenkov et. al (2007)

bool CGasSpecies::calcLVG(CGridBasic * grid, double kepler_star_mass, bool full)
{
    uint nr_of_total_energy_levels = getNrOfTotalEnergyLevels();
    uint nr_of_total_transitions = getNrOfTotalTransitions();
    uint nr_of_spectral_lines = getNrOfSpectralLines();
    float last_percentage = 0;
    long max_cells = grid->getMaxDataCells();
    bool no_error = true;

    cout << CLR_LINE;

    double * J_ext = new double[nr_of_total_transitions];
    for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
    {
        uint i_lvl_l = getLowerEnergyLevel(i_trans);
        uint i_lvl_u = getUpperEnergyLevel(i_trans);

        for(uint i_sublvl_u = 0; i_sublvl_u < getNrOfSublevel(i_lvl_u); i_sublvl_u++)
        {
            for(uint i_sublvl_l = 0; i_sublvl_l < getNrOfSublevel(i_lvl_l); i_sublvl_l++)
            {
                uint i_tmp_trans = getUniqueTransitionIndex(i_trans, i_sublvl_u, i_sublvl_l);
                J_ext[i_tmp_trans] = CMathFunctions::planck_hz(getTransitionFrequency(i_trans), 2.75);
            }
        }
    }

#pragma omp parallel for schedule(dynamic)
    for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
    {
        cell_basic * cell = grid->getCellFromIndex(i_cell);
        Matrix2D A(nr_of_total_energy_levels, nr_of_total_energy_levels);

        double * J_mid = new double[nr_of_total_transitions];
        double * b = new double[nr_of_total_energy_levels];
        double * tmp_lvl_pop = new double[nr_of_total_energy_levels];
        double * old_pop = new double[nr_of_total_energy_levels];

        for(uint i = 1; i < nr_of_total_energy_levels; i++)
        {
            tmp_lvl_pop[i] = 0.0;
            old_pop[i] = 0.0;
        }
        old_pop[0] = 1.0;
        tmp_lvl_pop[0] = 1.0;

        double turbulent_velocity = grid->getTurbulentVelocity(cell);

        // Calculate percentage of total progress per source
        float percentage = 100 * float(i_cell) / float(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
#pragma omp critical
            {
                cout << "-> Calculating LVG level population : " << last_percentage
                     << " [%]                       \r";
                last_percentage = percentage;
            }
        }

        // Check temperature and density
        double temp_gas = grid->getGasTemperature(cell);
        double dens_species = getNumberDensity(grid, cell);
        if(temp_gas < 1e-200 || dens_species < 1e-200)
            continue;

        Vector3D pos_xyz_cell = grid->getCenter(cell);
        double abs_vel;
        if(kepler_star_mass > 0)
        {
            Vector3D velo = CMathFunctions::calcKeplerianVelocity(pos_xyz_cell, kepler_star_mass);
            abs_vel = sqrt(pow(velo.X(), 2) + pow(velo.Y(), 2) + pow(velo.Z(), 2));
        }
        else if(grid->hasVelocityField())
        {
            Vector3D velo = grid->getVelocityField(cell);
            abs_vel = sqrt(pow(velo.X(), 2) + pow(velo.Y(), 2) + pow(velo.Z(), 2));
        }
        else
            abs_vel = 0;

        if(getGaussA(temp_gas, turbulent_velocity) * abs_vel < 1e-16)
            continue;

        double R_mid = sqrt(pow(pos_xyz_cell.X(), 2) + pow(pos_xyz_cell.Y(), 2));
        double L = R_mid * sqrt(2.0 / 3.0 / (getGaussA(temp_gas, turbulent_velocity) * abs_vel));
        uint i_iter = 0;

        double *** final_col_para = calcCollisionParameter(grid, cell);

        for(i_iter = 0; i_iter < MAX_LVG_ITERATIONS; i_iter++)
        {
            for(uint i_lvl = 0; i_lvl < nr_of_total_energy_levels; i_lvl++)
                old_pop[i_lvl] = tmp_lvl_pop[i_lvl];

            for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
            {
                uint i_lvl_l = getLowerEnergyLevel(i_trans);
                uint i_lvl_u = getUpperEnergyLevel(i_trans);

                for(uint i_sublvl_u = 0; i_sublvl_u < getNrOfSublevel(i_lvl_u); i_sublvl_u++)
                {
                    for(uint i_sublvl_l = 0; i_sublvl_l < getNrOfSublevel(i_lvl_l); i_sublvl_l++)
                    {
                        uint i_tmp_trans = getUniqueTransitionIndex(i_trans, i_sublvl_u, i_sublvl_l);
                        J_mid[i_tmp_trans] = elem_LVG(grid,
                                                      dens_species,
                                                      tmp_lvl_pop,
                                                      getGaussA(temp_gas, turbulent_velocity),
                                                      L,
                                                      J_ext[i_trans],
                                                      i_trans,
                                                      i_sublvl_u,
                                                      i_sublvl_l);
                    }
                }
            }

            createMatrix(J_mid, A, b, final_col_para);
            CMathFunctions::gauss(A, b, tmp_lvl_pop, nr_of_total_energy_levels);

            double sum_p = 0;
            for(uint j_lvl = 0; j_lvl < nr_of_total_energy_levels; j_lvl++)
            {
                if(tmp_lvl_pop[j_lvl] >= 0)
                    sum_p += tmp_lvl_pop[j_lvl];
                else
                {
                    cout << "WARNING: Level population element not greater than zero! Level = " << j_lvl
                         << ", Level pop = " << tmp_lvl_pop[j_lvl] << endl;
                    no_error = false;
                }
            }
            for(uint j_lvl = 0; j_lvl < nr_of_total_energy_levels; j_lvl++)
                tmp_lvl_pop[j_lvl] /= sum_p;

            uint j_lvl = 0;
            for(uint i_lvl = 1; i_lvl < nr_of_total_energy_levels; i_lvl++)
                if(tmp_lvl_pop[i_lvl] > tmp_lvl_pop[j_lvl])
                    j_lvl = i_lvl;

            if(i_iter > 1)
            {
                if(abs(tmp_lvl_pop[j_lvl] - old_pop[j_lvl]) /
                       (old_pop[j_lvl] + numeric_limits<double>::epsilon()) <
                   1e-2)
                {
                    break;
                }
            }
        }
        if(i_iter == MAX_LVG_ITERATIONS)
            cout << "WARNING: Maximum iteration needed in cell: " << i_cell << endl;

        if(full)
        {
            for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
            {
                for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
                {
                    grid->setLvlPop(cell, i_lvl, i_sublvl, tmp_lvl_pop[getUniqueLevelIndex(i_lvl, i_sublvl)]);
                }
            }
        }
        else
        {
            for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
            {
                for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
                {
                    // Get indices
                    uint i_trans = getTransitionFromSpectralLine(i_line);
                    uint i_lvl_l = getLowerEnergyLevel(i_trans);
                    uint i_lvl_u = getUpperEnergyLevel(i_trans);

                    if(i_lvl == i_lvl_l)
                    {
                        for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
                            grid->setLvlPopLower(
                                cell, i_line, i_sublvl, tmp_lvl_pop[getUniqueLevelIndex(i_lvl, i_sublvl)]);
                    }
                    else if(i_lvl == i_lvl_u)
                    {
                        for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
                            grid->setLvlPopUpper(
                                cell, i_line, i_sublvl, tmp_lvl_pop[getUniqueLevelIndex(i_lvl, i_sublvl)]);
                    }
                }
            }
        }

        delete[] J_mid;
        delete[] b;
        delete[] tmp_lvl_pop;
        delete[] old_pop;
        for(uint i_col_partner = 0; i_col_partner < nr_of_col_partner; i_col_partner++)
        {
            for(uint i_col_transition = 0; i_col_transition < nr_of_col_transition[i_col_partner];
                i_col_transition++)
            {
                delete[] final_col_para[i_col_partner][i_col_transition];
            }
            delete[] final_col_para[i_col_partner];
        }
        delete[] final_col_para;
    }
    delete[] J_ext;
    return no_error;
}

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

double CGasSpecies::elem_LVG(CGridBasic * grid,
                             double dens_species,
                             double * tmp_lvl_pop,
                             double gauss_a,
                             double L,
                             double J_ext,
                             uint i_trans,
                             uint i_sublvl_u,
                             uint i_sublvl_l)
{
    uint i_lvl_l = getLowerEnergyLevel(i_trans);
    uint i_lvl_u = getUpperEnergyLevel(i_trans);

    double lvl_pop_l = tmp_lvl_pop[getUniqueLevelIndex(i_lvl_l, i_sublvl_l)];
    double lvl_pop_u = tmp_lvl_pop[getUniqueLevelIndex(i_lvl_u, i_sublvl_u)];

    double Einst_A = getEinsteinA(i_trans);
    double Einst_B_u = getEinsteinBu(i_trans);
    double Einst_B_l = getEinsteinBl(i_trans);
    double j, alpha, tau, beta, J_mid, S;

    j = dens_species * lvl_pop_u * Einst_A * gauss_a * con_eps / sqrt(PI);
    alpha = dens_species * (lvl_pop_l * Einst_B_l - lvl_pop_u * Einst_B_u) * gauss_a * con_eps / sqrt(PI);

    if(alpha < 1e-20)
    {
        S = 0.0;
        alpha = 0.0;
    }
    else
        S = j / alpha;

    tau = alpha * L;

    if(tau < 1e-6)
        beta = 1.0 - 0.5 * tau;
    else
        beta = (1.0 - exp(-tau)) / tau;

    J_mid = (1.0 - beta) * S + beta * J_ext;
    return J_mid;
}

bool CGasSpecies::updateLevelPopulation(CGridBasic * grid, cell_basic * cell, double * J_total)
{
    // uint nr_of_energy_level = getNrOfEnergyLevels();
    uint nr_of_total_energy_levels = getNrOfTotalEnergyLevels();
    uint nr_of_spectral_lines = getNrOfSpectralLines();

    double * tmp_lvl_pop = new double[nr_of_total_energy_levels];

    double gas_number_density = grid->getGasNumberDensity(cell);
    if(gas_number_density > 1e-200)
    {
        Matrix2D A(nr_of_total_energy_levels, nr_of_total_energy_levels);
        double * b = new double[nr_of_total_energy_levels];

        double *** final_col_para = calcCollisionParameter(grid, cell);
        createMatrix(J_total, A, b, final_col_para);
        CMathFunctions::gauss(A, b, tmp_lvl_pop, nr_of_total_energy_levels);

        delete[] b;
        for(uint i_col_partner = 0; i_col_partner < nr_of_col_partner; i_col_partner++)
        {
            for(uint i_col_transition = 0; i_col_transition < nr_of_col_transition[i_col_partner];
                i_col_transition++)
            {
                delete[] final_col_para[i_col_partner][i_col_transition];
            }
            delete[] final_col_para[i_col_partner];
        }
        delete[] final_col_para;

        double sum_p = 0;
        for(uint i_lvl_tot = 0; i_lvl_tot < nr_of_total_energy_levels; i_lvl_tot++)
        {
            if(tmp_lvl_pop[i_lvl_tot] >= 0)
                sum_p += tmp_lvl_pop[i_lvl_tot];
            else
            {
                cout << "WARNING: Level population element not greater than zero!" << endl;
                return false;
            }
        }
        for(uint i_lvl_tot = 0; i_lvl_tot < nr_of_total_energy_levels; i_lvl_tot++)
            tmp_lvl_pop[i_lvl_tot] /= sum_p;
    }
    else
    {
        for(uint i_lvl_tot = 0; i_lvl_tot < nr_of_total_energy_levels; i_lvl_tot++)
            tmp_lvl_pop[i_lvl_tot] = 0;
        tmp_lvl_pop[0] = 1;
    }

    for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
    {
        for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
        {
            grid->setLvlPop(cell, i_lvl, i_sublvl, tmp_lvl_pop[getUniqueLevelIndex(i_lvl, i_sublvl)]);
        }
    }

    delete[] tmp_lvl_pop;
    return true;
}

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

double *** CGasSpecies::calcCollisionParameter(CGridBasic * grid, cell_basic * cell)
{
    // Init variables
    double temp_gas = grid->getGasTemperature(cell);
    uint nr_of_col_partner = getNrOfCollisionPartner();
    double *** final_col_para = new double **[nr_of_col_partner];

    for(uint i_col_partner = 0; i_col_partner < nr_of_col_partner; i_col_partner++)
    {
        uint nr_of_col_transition = getNrCollisionTransitions(i_col_partner);

        // Resize Matrix for lu and ul and all transitions
        final_col_para[i_col_partner] = new double *[nr_of_col_transition];

        // Get density of collision partner
        double dens = getColPartnerDensity(grid, cell, i_col_partner);

        bool skip = false;
        for(uint i_col_transition = 0; i_col_transition < nr_of_col_transition; i_col_transition++)
            if(getUpperCollision(i_col_partner, i_col_transition) == 0)
                skip = true;

        if(!skip)
        {
            uint hi_i = 0;
            if(temp_gas < getCollisionTemp(i_col_partner, 0))
                hi_i = 0;
            else if(temp_gas > getCollisionTemp(i_col_partner, getNrCollisionTemps(i_col_partner) - 1))
                hi_i = getNrCollisionTemps(i_col_partner) - 2;
            else
                hi_i = CMathFunctions::biListIndexSearch(
                    temp_gas, collision_temp[i_col_partner], getNrCollisionTemps(i_col_partner));

            for(uint i_col_transition = 0; i_col_transition < nr_of_col_transition; i_col_transition++)
            {
                // Calculate collision rates
                dlist rates = calcCollisionRate(i_col_partner, i_col_transition, hi_i, temp_gas);

                final_col_para[i_col_partner][i_col_transition] = new double[2];
                final_col_para[i_col_partner][i_col_transition][0] = rates[0] * dens;
                final_col_para[i_col_partner][i_col_transition][1] = rates[1] * dens;
            }
        }
    }
    return final_col_para;
}

double CGasSpecies::getColPartnerDensity(CGridBasic * grid, cell_basic * cell, uint i_col_partner)
{
    double dens = 0;
    switch(getOrientation_H2(i_col_partner))
    {
        case(COL_H2_FULL):
            dens = grid->getGasNumberDensity(cell);
            break;
        case(COL_H2_PARA):
            dens = grid->getGasNumberDensity(cell) * 0.25;
            break;
        case(COL_H2_ORTH):
            dens = grid->getGasNumberDensity(cell) * 0.75;
            break;
        case(COL_HE_FULL):
            dens = grid->getGasNumberDensity(cell) * max(0.0, grid->getMu() - 1.0);
            break;
        default:
            dens = 0;
            break;
    }
    return dens;
}

dlist CGasSpecies::calcCollisionRate(uint i_col_partner, uint i_col_transition, uint hi_i, double temp_gas)
{
    dlist col_mtr_tmp(2);

    double interp = CMathFunctions::interpolate(getCollisionTemp(i_col_partner, hi_i),
                                                getCollisionTemp(i_col_partner, hi_i + 1),
                                                getCollisionMatrix(i_col_partner, i_col_transition, hi_i),
                                                getCollisionMatrix(i_col_partner, i_col_transition, hi_i + 1),
                                                temp_gas);

    if(interp > 0)
    {
        col_mtr_tmp[0] = interp;

        col_mtr_tmp[1] =
            col_mtr_tmp[0] * getGLevel(getUpperCollision(i_col_partner, i_col_transition)) /
            getGLevel(getLowerCollision(i_col_partner, i_col_transition)) *
            exp(-con_h *
                (getEnergylevel(getUpperCollision(i_col_partner, i_col_transition)) * con_c * 100.0 -
                 getEnergylevel(getLowerCollision(i_col_partner, i_col_transition)) * con_c * 100.0) /
                (con_kB * temp_gas));
    }

    return col_mtr_tmp;
}

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

void CGasSpecies::createMatrix(double * J_mid, Matrix2D & A, double * b, double *** final_col_para)
{
    // Get number of colision partner and energy levels
    uint nr_of_col_partner = getNrOfCollisionPartner();
    uint nr_of_energy_level = getNrOfEnergyLevels();

    for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
    {
        for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
        {
            for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
            {
                uint i_tmp_lvl = getUniqueLevelIndex(i_lvl, i_sublvl);

                if(i_lvl == getUpperEnergyLevel(i_trans))
                {
                    uint k_lvl = getLowerEnergyLevel(i_trans);
                    for(uint k_sublvl = 0; k_sublvl < getNrOfSublevel(k_lvl); k_sublvl++)
                    {
                        uint k_tmp_lvl = getUniqueLevelIndex(k_lvl, k_sublvl);
                        uint i_tmp_trans = getUniqueTransitionIndex(i_trans, i_sublvl, k_sublvl);
                        A(i_tmp_lvl, k_tmp_lvl) += getEinsteinBl(i_trans) * J_mid[i_tmp_trans];
                        A(i_tmp_lvl, i_tmp_lvl) -=
                            getEinsteinA(i_trans) + getEinsteinBu(i_trans) * J_mid[i_tmp_trans];
                    }
                }
                else if(i_lvl == getLowerEnergyLevel(i_trans))
                {
                    uint k_lvl = getUpperEnergyLevel(i_trans);
                    for(uint k_sublvl = 0; k_sublvl < getNrOfSublevel(k_lvl); k_sublvl++)
                    {
                        uint k_tmp_lvl = getUniqueLevelIndex(k_lvl, k_sublvl);
                        uint i_tmp_trans = getUniqueTransitionIndex(i_trans, i_sublvl, k_sublvl);
                        A(i_tmp_lvl, k_tmp_lvl) +=
                            getEinsteinA(i_trans) + getEinsteinBu(i_trans) * J_mid[i_tmp_trans];
                        A(i_tmp_lvl, i_tmp_lvl) -= getEinsteinBl(i_trans) * J_mid[i_tmp_trans];
                    }
                }
            }
        }
    }

    for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
    {
        for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
        {
            uint i_tmp_lvl = getUniqueLevelIndex(i_lvl, i_sublvl);
            for(uint i_col_partner = 0; i_col_partner < nr_of_col_partner; i_col_partner++)
            {
                // Get number of transitions covered by the collisions
                uint nr_of_col_transition = getNrCollisionTransitions(i_col_partner);

                // Add collision rates to the system of equations for each transition
                for(uint i_col_transition = 0; i_col_transition < nr_of_col_transition; i_col_transition++)
                {
                    if(i_lvl == getLowerCollision(i_col_partner, i_col_transition))
                    {
                        uint k_lvl = getUpperCollision(i_col_partner, i_col_transition);
                        for(uint k_sublvl = 0; k_sublvl < getNrOfSublevel(k_lvl); k_sublvl++)
                        {
                            uint k_tmp_lvl = getUniqueLevelIndex(k_lvl, k_sublvl);
                            A(i_tmp_lvl, k_tmp_lvl) += final_col_para[i_col_partner][i_col_transition][0];
                            A(i_tmp_lvl, i_tmp_lvl) -= final_col_para[i_col_partner][i_col_transition][1];
                        }
                    }
                    else if(i_lvl == getUpperCollision(i_col_partner, i_col_transition))
                    {
                        uint k_lvl = getLowerCollision(i_col_partner, i_col_transition);
                        for(uint k_sublvl = 0; k_sublvl < getNrOfSublevel(k_lvl); k_sublvl++)
                        {
                            uint k_tmp_lvl = getUniqueLevelIndex(k_lvl, k_sublvl);
                            A(i_tmp_lvl, k_tmp_lvl) += final_col_para[i_col_partner][i_col_transition][1];
                            A(i_tmp_lvl, i_tmp_lvl) -= final_col_para[i_col_partner][i_col_transition][0];
                        }
                    }
                }
            }
        }
    }

    for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
    {
        for(uint i_sublvl = 0; i_sublvl < getNrOfSublevel(i_lvl); i_sublvl++)
        {
            uint i_tmp_lvl = getUniqueLevelIndex(i_lvl, i_sublvl);
            A(0, i_tmp_lvl) = 1.0;
            b[i_tmp_lvl] = 0;
        }
    }
    b[0] = 1.0;
}

StokesVector CGasSpecies::calcEmissivity(CGridBasic * grid, cell_basic * cell, uint i_line, uint i_sublvl)
{
    // Get width of spectral line
    double gauss_a = grid->getGaussA(cell);

    // Get index of transition
    uint i_trans = getTransitionFromSpectralLine(i_line);

    // Calculate the optical depth of the gas particles in the current cell
    double alpha_gas = (grid->getLvlPopLower(cell, i_line, i_sublvl) * getEinsteinBl(i_trans) -
                        grid->getLvlPopUpper(cell, i_line, i_sublvl) * getEinsteinBu(i_trans)) *
                       con_eps * gauss_a;
    // Calculate the Emissivity of the gas particles in the current cell
    double S_gas = grid->getLvlPopUpper(cell, i_line, i_sublvl) * getEinsteinA(i_trans) * con_eps * gauss_a;

    return StokesVector(S_gas, 0, 0, 0, alpha_gas);
}

StokesVector CGasSpecies::calcEmissivityForTransition(CGridBasic * grid,
                                                      cell_basic * cell,
                                                      uint i_trans,
                                                      uint i_sublvl_u,
                                                      uint i_sublvl_l)
{
    // Get width of spectral line
    double gauss_a = grid->getGaussA(cell);

    // Calculate the optical depth of the gas particles in the current cell
    double alpha_gas =
        (grid->getLvlPop(cell, getLowerEnergyLevel(i_trans), i_sublvl_l) * getEinsteinBl(i_trans) -
         grid->getLvlPop(cell, getUpperEnergyLevel(i_trans), i_sublvl_u) * getEinsteinBu(i_trans)) *
        con_eps * gauss_a;

    // Calculate the Emissivity of the gas particles in the current cell
    double S_gas = grid->getLvlPop(cell, getUpperEnergyLevel(i_trans), i_sublvl_u) * getEinsteinA(i_trans) *
                   con_eps * gauss_a;

    return StokesVector(S_gas, 0, 0, 0, alpha_gas);
}

StokesVector CGasSpecies::calcEmissivityForTransition(CGridBasic * grid, cell_basic * cell, uint i_trans)
{

    // Init resulting stokes vector
    StokesVector tmp_stokes;

    // Get all indices and number of sublevels
    uint i_lvl_u = getUpperEnergyLevel(i_trans);
    uint i_lvl_l = getLowerEnergyLevel(i_trans);
    uint nr_of_sublevel_u = getNrOfSublevel(i_lvl_u);
    uint nr_of_sublevel_l = getNrOfSublevel(i_lvl_l);

    for(uint i_sublvl_u = 0; i_sublvl_u < nr_of_sublevel_u; i_sublvl_u++)
    {
        for(uint i_sublvl_l = 0; i_sublvl_l < nr_of_sublevel_l; i_sublvl_l++)
        {
            tmp_stokes += calcEmissivityForTransition(grid, cell, i_trans, i_sublvl_u, i_sublvl_l);
        }
    }
    return tmp_stokes;
}

Matrix2D CGasSpecies::getGaussLineMatrix(CGridBasic * grid, photon_package * pp, double velocity)
{
    // Init line matrix
    Matrix2D line_matrix(4, 4);

    // Calculate gaussian shape
    double line_amplitude = getGaussLineShape(grid, pp, velocity);

    // Only diagonal without polarization rotation matrix elements
    for(uint i = 0; i < 4; i++)
        line_matrix(i, i) = line_amplitude;

    return line_matrix;
}

double CGasSpecies::getGaussLineShape(CGridBasic * grid, photon_package * pp, double velocity)
{
    double gauss_a = grid->getGaussA(pp);
    return exp(-(pow(velocity, 2) * pow(gauss_a, 2))) / PIsq;
}

Matrix2D CGasSpecies::getZeemanSplittingMatrix(photon_package * pp,
                                               uint i_trans,
                                               double velocity,
                                               double Gamma,
                                               double doppler_width,
                                               double voigt_a,
                                               Vector3D mag_field,
                                               double cos_theta,
                                               double sin_theta,
                                               double cos_2_phi,
                                               double sin_2_phi)
{
    // Init variables
    Matrix2D line_matrix(4, 4);
    double line_strength, freq_shift;
    double f_doppler, mult_A, mult_B;
    double frequency = getTransitionFrequency(i_trans);
    uint i_pi = 0, i_sigma_p = 0, i_sigma_m = 0;

    // Init the current value of the line function as a complex value
    complex<double> line_function;

    // Calculate the contribution of each allowed transition between Zeeman sublevels
    for(uint i_sublvl_u = 0; i_sublvl_u < getNrOfSublevelUpper(i_trans); i_sublvl_u++)
    {
        float sublvl_u = -getMaxMUpper(i_trans) + i_sublvl_u;
        for(uint i_sublvl_l = 0; i_sublvl_l < getNrOfSublevelLower(i_trans); i_sublvl_l++)
        {
            float sublvl_l = max(sublvl_u - 1, -getMaxMLower(i_trans)) + i_sublvl_l;

            // Calculate the frequency shift in relation to the unshifted line peak
            // Delta nu = (B * mu_Bohr) / h * (m' * g' - m'' * g'')
            freq_shift = mag_field.length() * con_mb / con_h *
                         (sublvl_u * getLandeUpper(i_trans) - sublvl_l * getLandeLower(i_trans));

            // Calculate the frequency value of the current velocity channel in
            // relation to the peak of the line function
            f_doppler = CMathFunctions::Velo2Freq(velocity + CMathFunctions::Freq2Velo(freq_shift, frequency),
                                                  frequency) /
                        doppler_width;

            // Calculate the line function value at the frequency
            // of the current velocity channel
            line_function = getLineShape_AB(f_doppler, voigt_a);

            // Multiply the line function value by PIsq for normalization
            mult_A = real(line_function) / PIsq;

            // Multiply the line function value by PIsq for normalization
            // and divide it by 2 to take the difference between the faddeeva
            // and Faraday-Voigt function into account
            mult_B = imag(line_function) / (PIsq * 2.0);

            // Use the correct propagation matrix and relative line strength that
            // depends on the current type of Zeeman transition (pi, sigma_-, sigma_+)
            switch(int(sublvl_l - sublvl_u))
            {
                case TRANS_SIGMA_P:
                    // Get the relative line strength from Zeeman file
                    line_strength = getLineStrengthSigmaP(i_trans, i_sigma_p);

                    // Get the propagation matrix for extinction/emission
                    line_matrix += CMathFunctions::getPropMatrixASigmaP(
                        cos_theta, sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_A);

                    // Get the propagation matrix for Faraday rotation
                    line_matrix += CMathFunctions::getPropMatrixBSigmaP(
                        cos_theta, sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_B);

                    // Increase the sigma_+ counter to circle through the line strengths
                    i_sigma_p++;
                    break;
                case TRANS_PI:
                    // Get the relative line strength from Zeeman file
                    line_strength = getLineStrengthPi(i_trans, i_pi);

                    // Get the propagation matrix for extinction/emission
                    line_matrix += CMathFunctions::getPropMatrixAPi(
                        cos_theta, sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_A);

                    // Get the propagation matrix for Faraday rotation
                    line_matrix += CMathFunctions::getPropMatrixBPi(
                        cos_theta, sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_B);

                    // Increase the pi counter to circle through the line strengths
                    i_pi++;
                    break;
                case TRANS_SIGMA_M:
                    // Get the relative line strength from Zeeman file
                    line_strength = getLineStrengthSigmaM(i_trans, i_sigma_m);

                    // Get the propagation matrix for extinction/emission
                    line_matrix += CMathFunctions::getPropMatrixASigmaM(
                        cos_theta, sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_A);

                    // Get the propagation matrix for Faraday rotation
                    line_matrix += CMathFunctions::getPropMatrixBSigmaM(
                        cos_theta, sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_B);

                    // Increase the sigma_- counter to circle through the line strengths
                    i_sigma_m++;
                    break;
            }
        }
    }
    return line_matrix;
}

bool CGasSpecies::readGasParamaterFile(string _filename, uint id, uint max)
{
    uint line_counter, cmd_counter;
    uint pos_counter = 0;
    uint row_offset, i_col_partner, i_col_transition;
    CCommandParser ps;
    fstream reader(_filename.c_str());
    unsigned char ru[4] = { '|', '/', '-', '\\' };
    string line;
    dlist values;

    cout << CLR_LINE;

    if(reader.fail())
    {
        cout << "\nERROR: Cannot open gas_species catalog:" << endl;
        cout << _filename << endl;
        return false;
    }

    line_counter = 0;
    cmd_counter = 0;

    row_offset = 0;
    i_col_partner = 0;
    i_col_transition = 0;

    uint char_counter = 0;

    while(getline(reader, line))
    {
        line_counter++;

        if(line_counter % 20 == 0)
        {
            char_counter++;
            cout << "-> Reading gas species file nr. " << id + 1 << " of " << max << " : "
                 << ru[(uint)char_counter % 4] << "           \r";
        }

        ps.formatLine(line);

        if(line.size() == 0)
            continue;

        if(cmd_counter != 0)
        {
            values = ps.parseValues(line);
            if(values.size() == 0)
                continue;
        }

        cmd_counter++;

        if(cmd_counter == 1)
            stringID = line;
        else if(cmd_counter == 2)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!"
                     << endl;
                return false;
            }
            molecular_weight = values[0];
        }
        else if(cmd_counter == 3)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!"
                     << endl;
                return false;
            }
            nr_of_energy_level = uint(values[0]);

            energy_level = new double[nr_of_energy_level];

            // Set number of sublevel to one and increase it in case of Zeeman
            nr_of_sublevel = new int[nr_of_energy_level];
            for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
            {
                nr_of_sublevel[i_lvl] = 1;
            }

            g_level = new double[nr_of_energy_level];
            j_level = new double[nr_of_energy_level];
        }
        else if(cmd_counter < 4 + nr_of_energy_level && cmd_counter > 3)
        {
            if(values.size() < 4)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!"
                     << endl;
                return false;
            }

            // ENERGIES(cm^-1)
            energy_level[pos_counter] = values[1];
            // WEIGHT
            g_level[pos_counter] = values[2];
            // J
            j_level[pos_counter] = values[3];

            pos_counter++;
        }
        else if(cmd_counter == 4 + nr_of_energy_level)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!"
                     << endl;
                return false;
            }
            nr_of_transitions = uint(values[0]);
            pos_counter = 0;

            upper_level = new int[nr_of_transitions];
            lower_level = new int[nr_of_transitions];
            trans_einstA = new double[nr_of_transitions];
            trans_freq = new double[nr_of_transitions];
            trans_inner_energy = new double[nr_of_transitions];

            trans_einstB_u = new double[nr_of_transitions];
            trans_einstB_l = new double[nr_of_transitions];
        }
        else if(cmd_counter < 5 + nr_of_energy_level + nr_of_transitions &&
                cmd_counter > 4 + nr_of_energy_level)
        {
            if(values.size() != 6)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!"
                     << endl;
                return false;
            }

            for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
            {
                if(spectral_lines[i_line] == int(values[0] - 1))
                {
                    unique_spectral_lines.push_back(i_line);
                    break;
                }
            }

            // UP (as index starting with 0)
            upper_level[pos_counter] = int(values[1] - 1);
            // LOW (as index starting with 0)
            lower_level[pos_counter] = int(values[2] - 1);
            // EINSTEINA(s^-1)
            trans_einstA[pos_counter] = values[3];
            // FREQ(GHz -> Hz)
            trans_freq[pos_counter] = values[4] * 1e9;
            // E_u(K)
            trans_inner_energy[pos_counter] = values[5];

            trans_einstB_u[pos_counter] = trans_einstA[pos_counter] *
                                          pow(con_c / trans_freq[pos_counter], 2.0) /
                                          (2.0 * con_h * trans_freq[pos_counter]);

            trans_einstB_l[pos_counter] = g_level[upper_level[pos_counter]] /
                                          g_level[lower_level[pos_counter]] * trans_einstB_u[pos_counter];

            pos_counter++;
        }
        else if(cmd_counter == 5 + nr_of_energy_level + nr_of_transitions)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!"
                     << endl;
                return false;
            }
            nr_of_col_partner = uint(values[0]);

            nr_of_col_transition = new int[nr_of_col_partner];
            nr_of_col_temp = new int[nr_of_col_partner];
            orientation_H2 = new int[nr_of_col_partner];

            collision_temp = new double *[nr_of_col_partner];
            col_upper = new int *[nr_of_col_partner];
            col_lower = new int *[nr_of_col_partner];

            col_matrix = new double **[nr_of_col_partner];

            for(uint i = 0; i < nr_of_col_partner; i++)
            {
                collision_temp[i] = 0;
                col_upper[i] = 0;
                col_lower[i] = 0;

                nr_of_col_transition[i] = 0;
                orientation_H2[i] = 0;
                nr_of_col_temp[i] = 0;
                col_matrix[i] = 0;
            }
        }
        else if(cmd_counter == 6 + nr_of_energy_level + nr_of_transitions)
        {
            if(values[0] < 1)
            {
                cout << "\nERROR: Line " << line_counter
                     << " wrong orientation of H2 collision partner (gas species file)!" << endl;
                return false;
            }
            orientation_H2[i_col_partner] = int(values[0]);
        }
        else if(cmd_counter == 7 + nr_of_energy_level + nr_of_transitions + row_offset)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!"
                     << endl;
                return false;
            }

            nr_of_col_transition[i_col_partner] = int(values[0]);
            col_upper[i_col_partner] = new int[int(values[0])];
            col_lower[i_col_partner] = new int[int(values[0])];
            col_matrix[i_col_partner] = new double *[int(values[0])];

            for(uint i = 0; i < uint(values[0]); i++)
            {
                col_matrix[i_col_partner][i] = 0;
                col_upper[i_col_partner][i] = 0;
                col_lower[i_col_partner][i] = 0;
            }
        }
        else if(cmd_counter == 8 + nr_of_energy_level + nr_of_transitions + row_offset)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!"
                     << endl;
                return false;
            }

            nr_of_col_temp[i_col_partner] = int(values[0]);

            collision_temp[i_col_partner] = new double[int(values[0])];

            for(uint i = 0; i < uint(values[0]); i++)
                collision_temp[i_col_partner][i] = 0;
        }
        else if(cmd_counter == 9 + nr_of_energy_level + nr_of_transitions + row_offset)
        {
            if(values.size() != uint(nr_of_col_temp[i_col_partner]))
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!"
                     << endl;
                return false;
            }

            for(uint i = 0; i < uint(nr_of_col_temp[i_col_partner]); i++)
                collision_temp[i_col_partner][i] = values[i];
        }
        else if(cmd_counter < 10 + nr_of_energy_level + nr_of_transitions +
                                  nr_of_col_transition[i_col_partner] + row_offset &&
                cmd_counter > 9 + nr_of_energy_level + row_offset)
        {
            if(values.size() != uint(nr_of_col_temp[i_col_partner] + 3))
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!"
                     << endl;
                return false;
            }

            i_col_transition = cmd_counter - 10 - nr_of_energy_level - nr_of_transitions - row_offset;
            col_upper[i_col_partner][i_col_transition] = int(values[1] - 1);
            col_lower[i_col_partner][i_col_transition] = int(values[2] - 1);

            col_matrix[i_col_partner][i_col_transition] = new double[nr_of_col_temp[i_col_partner]];

            for(uint i = 0; i < uint(nr_of_col_temp[i_col_partner]); i++)
                col_matrix[i_col_partner][i_col_transition][i] = values[3 + i] * 1e-6;
        }
        else if(i_col_partner < nr_of_col_partner)
        {
            if(values[0] < 1)
            {
                cout << "\nERROR: Line " << line_counter
                     << " wrong orientation of H2 collision partner (gas species file)!" << endl;
                return false;
            }

            i_col_partner++;
            orientation_H2[i_col_partner] = int(values[0]);

            row_offset = cmd_counter - (6 + nr_of_energy_level + nr_of_transitions);
        }
    }
    reader.close();

    return true;
}

bool CGasSpecies::readZeemanParamaterFile(string _filename)
{
    // Init basic variables
    fstream reader(_filename.c_str());
    CCommandParser ps;
    string line;
    dlist values;

    // Init variables
    uint sublevels_upper_nr = 0, sublevels_lower_nr = 0;
    uint offset_pi = 0, offset_sigma = 0;
    uint i_sigma_p_transition = 0, i_sigma_m_transition = 0;
    uint i_pi_transition = 0, i_trans_zeeman = 0;

    line_strength_pi = new double *[nr_of_transitions];
    line_strength_sigma_p = new double *[nr_of_transitions];
    line_strength_sigma_m = new double *[nr_of_transitions];

    lande_factor = new double[nr_of_energy_level];
    for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
    {
        lande_factor[i_lvl] = 0;
    }

    if(reader.fail())
    {
        cout << "\nERROR: Cannot open Zeeman splitting catalog:" << endl;
        cout << _filename << endl;
        return false;
    }

    uint line_counter = 0;
    uint cmd_counter = 0;

    while(getline(reader, line))
    {
        line_counter++;

        ps.formatLine(line);

        if(line.size() == 0)
            continue;

        if(cmd_counter != 0)
        {
            values = ps.parseValues(line);
            if(values.size() == 0)
                continue;
        }

        cmd_counter++;

        if(cmd_counter == 1)
        {
            if(line != stringID)
            {
                cout << "wrong Zeeman splitting catalog file chosen!" << endl;
                return false;
            }
        }
        else if(cmd_counter == 2)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            gas_species_radius = values[0];
        }
        else if(cmd_counter == 3)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            nr_zeeman_spectral_lines = uint(values[0]);
        }
        else if(cmd_counter == 4)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
            {
                if(i_trans == int(values[0] - 1))
                {
                    // Set current zeeman transition index
                    i_trans_zeeman = int(values[0] - 1);
                    break;
                }
            }
        }
        else if(cmd_counter == 5)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            // Save the lande factor of the upper energy level
            lande_factor[getUpperEnergyLevel(i_trans_zeeman)] = values[0];
        }
        else if(cmd_counter == 6)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            // Save the lande factor of the lower energy level
            lande_factor[getLowerEnergyLevel(i_trans_zeeman)] = values[0];
        }
        else if(cmd_counter == 7)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            // Save the number of sublevel per energy level
            nr_of_sublevel[getUpperEnergyLevel(i_trans_zeeman)] = int(values[0]);
        }
        else if(cmd_counter == 8)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            // Save the number of sublevel per energy level
            nr_of_sublevel[getLowerEnergyLevel(i_trans_zeeman)] = int(values[0]);

            // Set local number of sublevel for the involved energy levels
            uint nr_of_sublevel_upper = nr_of_sublevel[getUpperEnergyLevel(i_trans_zeeman)];
            uint nr_of_sublevel_lower = nr_of_sublevel[getLowerEnergyLevel(i_trans_zeeman)];

            // Calculate the numbers of transitions possible for sigma and pi transitions
            uint nr_pi_spectral_lines = min(nr_of_sublevel_upper, nr_of_sublevel_lower);
            uint nr_sigma_spectral_lines = nr_of_sublevel_upper - 1;
            if(nr_of_sublevel_upper != nr_of_sublevel_lower)
            {
                nr_sigma_spectral_lines = min(nr_of_sublevel_upper, nr_of_sublevel_lower);
            }

            // Save offset to jump of values in the text file
            offset_pi = nr_pi_spectral_lines;
            offset_sigma = nr_sigma_spectral_lines;

            // Init pointer array for pi line strengths
            line_strength_pi[i_trans_zeeman] = new double[nr_pi_spectral_lines];
            for(int i_pi = 0; i_pi < nr_pi_spectral_lines; i_pi++)
                line_strength_pi[i_trans_zeeman][i_pi] = 0;

            // Init pointer array for sigma +/- line strengths
            line_strength_sigma_p[i_trans_zeeman] = new double[nr_sigma_spectral_lines];
            line_strength_sigma_m[i_trans_zeeman] = new double[nr_sigma_spectral_lines];
            for(int i_sigma = 0; i_sigma < nr_sigma_spectral_lines; i_sigma++)
            {
                line_strength_sigma_p[i_trans_zeeman][i_sigma] = 0;
                line_strength_sigma_m[i_trans_zeeman][i_sigma] = 0;
            }

            // Init indices to address the correct trantitions
            i_pi_transition = 0;
            i_sigma_p_transition = 0;
            i_sigma_m_transition = 0;
        }
        else if(cmd_counter <= 8 + offset_pi && cmd_counter > 8)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            line_strength_pi[i_trans_zeeman][i_pi_transition] = values[0];
            i_pi_transition++;
        }
        else if(cmd_counter <= 8 + offset_pi + offset_sigma && cmd_counter > 8 + offset_pi)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            line_strength_sigma_p[i_trans_zeeman][i_sigma_p_transition] = values[0];
            i_sigma_p_transition++;
        }
        else if(cmd_counter <= 8 + offset_pi + 2 * offset_sigma && cmd_counter > 8 + offset_pi + offset_sigma)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            line_strength_sigma_m[i_trans_zeeman][i_sigma_m_transition] = values[0];
            i_sigma_m_transition++;
        }
        else if(cmd_counter == 9 + offset_pi + 2 * offset_sigma)
        {
            cmd_counter = 3;
        }
    }
    reader.close();

    for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
        if(lande_factor[getTransitionFromSpectralLine(i_line)] == 0)
        {
            cout << SEP_LINE;
            cout << "\nERROR: For transition number " << uint(getTransitionFromSpectralLine(i_line) + 1)
                 << " exists no Zeeman splitting data" << endl;
            return false;
        }
    return true;
}

void CGasSpecies::calcLineBroadening(CGridBasic * grid)
{
    long max_cells = grid->getMaxDataCells();
#pragma omp parallel for schedule(dynamic)
    for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
    {
        cell_basic * cell = grid->getCellFromIndex(i_cell);

        // Get necessary quantities from the current cell
        double temp_gas = grid->getGasTemperature(cell);
        double dens_gas = grid->getGasNumberDensity(cell);
        double dens_species = getNumberDensity(grid, cell);
        double turbulent_velocity = grid->getTurbulentVelocity(cell);

        for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
        {
            uint i_trans = getTransitionFromSpectralLine(i_line);
            double frequency = getTransitionFrequency(i_trans);

            double gauss_a = getGaussA(temp_gas, turbulent_velocity);
            double doppler_width = frequency / (con_c * gauss_a);
            double Gamma = 0, voigt_a = 0;
            if(isTransZeemanSplit(i_trans))
            {
                Gamma = getGamma(i_trans, dens_gas, dens_species, temp_gas, turbulent_velocity);
                voigt_a = Gamma / (4 * PI * doppler_width);
            }
            grid->setLineBroadening(cell, i_line, gauss_a, doppler_width, Gamma, voigt_a);
        }
    }
}

void CGasSpecies::applyRadiationFieldFactor(uint i_trans,
                                            double sin_theta,
                                            double cos_theta,
                                            double energy,
                                            double * J_nu)
{
    if(!isTransZeemanSplit(i_trans))
    {
        // Update radiation field
        uint i_tmp_trans = getUniqueTransitionIndex(i_trans);
        J_nu[i_tmp_trans] += energy;
        return;
    }

    // Get indices to involved energy level
    uint i_lvl_u = getUpperEnergyLevel(i_trans);
    uint i_lvl_l = getLowerEnergyLevel(i_trans);

    // Init indices
    uint i_pi = 0, i_sigma_p = 0, i_sigma_m = 0;

    // Two loops for each possible transition between sublevels
    for(uint i_sublvl_u = 0; i_sublvl_u < getNrOfSublevel(i_lvl_u); i_sublvl_u++)
    {
        // Get indices to involved sublevel
        float sublvl_u = -getMaxM(i_lvl_u) + i_sublvl_u;

        for(uint i_sublvl_l = 0; i_sublvl_l < getNrOfSublevel(i_lvl_l); i_sublvl_l++)
        {
            // Get indices to involved sublevel
            float sublvl_l = max(sublvl_u - 1, -getMaxM(i_lvl_l)) + i_sublvl_l;

            // Init resulting radiation field multiplier
            double radiation_field_factor = 0;

            // Find correct transition
            switch(int(sublvl_l - sublvl_u))
            {
                case TRANS_SIGMA_P:
                    // Line strength times projection
                    radiation_field_factor =
                        getLineStrengthSigmaP(i_trans, i_sigma_p) * (1 + cos_theta * cos_theta);

                    // Increase the sigma_+ counter to circle through the line strengths
                    i_sigma_p++;
                    break;
                case TRANS_PI:
                    // Line strength times projection
                    radiation_field_factor = getLineStrengthPi(i_trans, i_pi) * sin_theta * sin_theta;

                    // Increase the pi counter to circle through the line strengths
                    i_pi++;
                    break;
                case TRANS_SIGMA_M:
                    // Line strength times projection
                    radiation_field_factor =
                        getLineStrengthSigmaM(i_trans, i_sigma_m) * (1 + cos_theta * cos_theta);

                    // Increase the sigma_- counter to circle through the line strengths
                    i_sigma_m++;
                    break;
            }

            // Update radiation field
            uint i_tmp_trans = getUniqueTransitionIndex(i_trans, i_sublvl_u, i_sublvl_l);
            J_nu[i_tmp_trans] += radiation_field_factor * energy;
        }
    }
}

bool CGasMixture::createGasSpecies(parameters & param)
{
    nr_of_species = param.getNrOfGasSpecies();
    single_species = new CGasSpecies[nr_of_species];

    for(uint i_species = 0; i_species < nr_of_species; i_species++)
    {
        single_species[i_species].setAbundance(param.getGasSpeciesAbundance(i_species));
        single_species[i_species].setLevelPopType(param.getGasSpeciesLevelPopType(i_species));
        single_species[i_species].setNrOfSpectralLines(param.getNrOfSpectralLines(i_species));
        single_species[i_species].setSpectralLines(param.getSpectralLines(i_species));

        string path = param.getGasSpeciesCatalogPath(i_species);
        if(!single_species[i_species].readGasParamaterFile(path, i_species, nr_of_species))
            return false;

        if(param.getZeemanCatalog(i_species) != "")
            if(!single_species[i_species].readZeemanParamaterFile(param.getZeemanCatalog(i_species)))
                return false;

        single_species[i_species].initReferenceLists();

        if((single_species[i_species].getLevelPopType() == POP_FEP ||
            single_species[i_species].getLevelPopType() == POP_LVG) &&
           single_species[i_species].getNrOfCollisionPartner() == 0)
        {
            cout << "\nERROR: FEP and LVG level population approximations require a gas "
                    "parameters file \n"
                    "       with collisional data (e.g. from Leiden Atomic and Molecular "
                    "Database)"
                 << endl;
            return false;
        }
    }

    setKeplerStarMass(param.getKeplerStarMass());
    return true;
}

bool CGasMixture::calcLevelPopulation(CGridBasic * grid, uint i_species)
{
    // Set way of level population calculation, if not set by function call
    uint lvl_pop_type = getLevelPopType(i_species);

    // Let the grid know where to put the level populations
    grid->setGasInformation(level_to_pos[i_species], line_to_pos[i_species]);

    switch(lvl_pop_type)
    {
        case POP_MC:
            if(!single_species[i_species].calcFEP(grid, true))
                return false;
            break;
        case POP_LTE:
            if(!single_species[i_species].calcLTE(grid))
                return false;
            break;
        case POP_FEP:
            if(!single_species[i_species].calcFEP(grid))
                return false;
            break;
        case POP_LVG:
            if(!single_species[i_species].calcLVG(grid, getKeplerStarMass()))
                return false;
            break;
        default:
            return false;
            break;
    }

    return true;
}

bool CGasMixture::updateLevelPopulation(CGridBasic * grid,
                                        cell_basic * cell,
                                        uint i_species,
                                        double * J_total)
{
    uint lvl_pop_type = getLevelPopType(i_species);

    // Only used for MC level population calculations
    switch(lvl_pop_type)
    {
        case POP_MC:
            return single_species[i_species].updateLevelPopulation(grid, cell, J_total);
            break;
        default:
            return false;
            break;
    }
    return true;
}

void CGasMixture::printParameter(parameters & param, CGridBasic * grid)
{
    cout << CLR_LINE;
    cout << "Gas parameters                             " << endl;
    cout << SEP_LINE;
    cout << "- Velocity field                : ";
    if(getKeplerStarMass() > 0)
        cout << "kepler rotation, M_star: " << getKeplerStarMass() << " [M_sun]\n"
             << "    HINT: only available with one central star" << endl;
    else if(grid->getVelocityFieldAvailable())
        cout << "velocity field of the grid is used" << endl;
    else
        cout << "velocity field is zero" << endl;
    cout << "- Turbulent Velocity            : ";
    if(param.getTurbulentVelocity() > 0)
        cout << param.getTurbulentVelocity() << " [m/s]" << endl;
    else
        cout << "turbulent velocity of the grid is used" << endl;

    for(uint i_species = 0; i_species < nr_of_species; i_species++)
    {
        cout << SEP_LINE;
        cout << "Gas species " << (i_species + 1) << " (" << getGasSpeciesName(i_species) << ")" << endl;

        if(single_species[i_species].getNrOfSpectralLines() == 0)
        {
            cout << "\nWARNING: No spectral lines selected!" << endl;
            return;
        }

        stringstream transition_str, vel_channels_str, max_vel_str;
        dlist line_ray_detectors = param.getLineRayDetector(i_species);
        for(uint i = 0; i < line_ray_detectors.size(); i += NR_OF_LINE_DET)
        {
            uint pos = i / NR_OF_LINE_DET;

            transition_str << uint(line_ray_detectors[i] + 1);
            max_vel_str << line_ray_detectors[i + 2];
            vel_channels_str << uint(line_ray_detectors[i + NR_OF_LINE_DET - 1]);
            if(i < line_ray_detectors.size() - NR_OF_LINE_DET)
            {
                transition_str << ", ";
                max_vel_str << ", ";
                vel_channels_str << ", ";
            }
        }
        cout << "- Line transition(s)            : " << transition_str.str() << endl;
        cout << "- Number of velocity channels   : " << vel_channels_str.str() << endl;
        cout << "- Velocity limit(s)             : " << max_vel_str.str() << " [m/s]" << endl;
        cout << "- Level population              : ";
        uint lvl_pop_type = getLevelPopType(i_species);
        switch(lvl_pop_type)
        {
            case POP_MC:
                cout << "Monte-Carlo" << endl;
                break;
            case POP_LTE:
                cout << "LTE" << endl;
                break;
            case POP_FEP:
                cout << "FEP" << endl;
                break;
            case POP_LVG:
                cout << "LVG" << endl;
                break;
            default:
                cout << "\nERROR: UNKNOWN!" << endl;
        }

        if(isZeemanSplit(i_species))
            cout << "- Particle radius (collisions)  : " << getCollisionRadius(i_species) << " [m]" << endl;

        cout << "- Molecular weight              : " << getMolecularWeight(i_species) << endl;
        double ab = getAbundance(i_species);
        if(ab > 0)
            cout << "- Abundance                     : " << ab << endl;
        else
        {
            double dens_species, min_dens_species = 1e200, max_dens_species = 0;
            for(long i_cell = 0; i_cell < grid->getMaxDataCells(); i_cell++)
            {
                cell_basic * cell = grid->getCellFromIndex(i_cell);

                // Get abundance of a certain gas species
                double dens_species =
                    grid->getCellAbundance(cell, uint(-ab - 1)) * grid->getGasNumberDensity(cell);

                if(dens_species < min_dens_species)
                    min_dens_species = dens_species;
                if(dens_species > max_dens_species)
                    max_dens_species = dens_species;
            }
            cout << "- Abundance from grid ID nr.    : " << int(-ab) << endl;
            cout << "                      (min,max) : [" << min_dens_species << ", " << max_dens_species
                 << "] [m^-3]" << endl;
        }

        double total_species_mass = 0;
        for(long i_cell = 0; i_cell < grid->getMaxDataCells(); i_cell++)
        {
            cell_basic * cell = grid->getCellFromIndex(i_cell);
            total_species_mass += getMassDensity(grid, cell, i_species) * grid->getVolume(cell);
        }
        cout << "- Total mass                    : " << total_species_mass / M_sun << " [M_sun], "
             << total_species_mass << " [kg]" << endl;

        for(uint i = 0; i < getUniqueTransitions(i_species).size(); i++)
        {
            uint i_line = getUniqueTransitions(i_species, i);
            uint i_trans = getTransitionFromSpectralLine(i_species, i_line);

            cout << SEP_LINE;
            cout << "Line transition " << (i_line + 1) << " (gas species " << (i_species + 1) << ")" << endl;
            cout << "- Transition number             : "
                 << uint(getTransitionFromSpectralLine(i_species, i_line) + 1) << endl;
            cout << "- Involved energy levels        : "
                 << getUpperEnergyLevel(i_species, getTransitionFromSpectralLine(i_species, i_line)) + 1
                 << " -> "
                 << getLowerEnergyLevel(i_species, getTransitionFromSpectralLine(i_species, i_line)) + 1
                 << endl;
            cout << "- Transition frequency          : " << getSpectralLineFrequency(i_species, i_line)
                 << " [Hz]" << endl;
            cout << "- Transition wavelength         : "
                 << (con_c / getSpectralLineFrequency(i_species, i_line)) << " [m]" << endl;
            if(isTransZeemanSplit(i_species, i_trans))
            {
                cout << CLR_LINE;
                cout << "Zeeman splitting parameters                " << endl;
                cout << "- Lande factor of upper level   : " << getLandeUpper(i_species, i_trans) << endl;
                cout << "- Lande factor of lower level   : " << getLandeLower(i_species, i_trans) << endl;
                cout << "- Sublevels in upper level      : " << getNrOfSublevelUpper(i_species, i_trans)
                     << endl;
                cout << "- Sublevels in lower level      : " << getNrOfSublevelLower(i_species, i_trans)
                     << endl;

                uint i_pi = 0, i_sigma_p = 0, i_sigma_m = 0;

                cout << "- Sigma+ line strength\t\tm(upper)\t\tm(lower)" << endl;
                for(uint i_sublvl_u = 0; i_sublvl_u < getNrOfSublevelUpper(i_species, i_trans); i_sublvl_u++)
                {
                    float sublvl_u = -getMaxMUpper(i_species, i_trans) + i_sublvl_u;
                    float sublvl_l = sublvl_u + 1;
                    if(abs(sublvl_l) <= getMaxMLower(i_species, i_trans))
                    {
                        char LineStrengthTmp[16];
#ifdef WINDOWS
                        _snprintf_s(LineStrengthTmp,
                                    sizeof(LineStrengthTmp),
                                    "%.3f",
                                    getLineStrengthSigmaP(i_species, i_trans, i_sigma_p));
#else
                        snprintf(LineStrengthTmp,
                                 sizeof(LineStrengthTmp),
                                 "%.3f",
                                 getLineStrengthSigmaP(i_species, i_trans, i_sigma_p));
#endif

                        cout << "\t" << LineStrengthTmp << "\t\t\t   " << float(sublvl_u) << "\t\t\t   "
                             << float(sublvl_l) << endl;
                        i_sigma_p++;
                    }
                }
                cout << "- Pi line strength\t\tm(upper)\t\tm(lower)" << endl;
                for(uint i_sublvl_u = 0; i_sublvl_u < getNrOfSublevelUpper(i_species, i_trans); i_sublvl_u++)
                {
                    float sublvl_u = -getMaxMUpper(i_species, i_trans) + i_sublvl_u;
                    float sublvl_l = sublvl_u;
                    if(abs(sublvl_l) <= getMaxMLower(i_species, i_trans))
                    {
                        char LineStrengthTmp[16];
#ifdef WINDOWS
                        _snprintf_s(LineStrengthTmp,
                                    sizeof(LineStrengthTmp),
                                    "%.3f",
                                    getLineStrengthPi(i_species, i_trans, i_pi));
#else
                        snprintf(LineStrengthTmp,
                                 sizeof(LineStrengthTmp),
                                 "%.3f",
                                 getLineStrengthPi(i_species, i_trans, i_pi));
#endif
                        cout << "\t" << LineStrengthTmp << "\t\t\t   " << float(sublvl_u) << "\t\t\t   "
                             << float(sublvl_l) << endl;
                        i_pi++;
                    }
                }
                cout << "- Sigma- line strength\t\tm(upper)\t\tm(lower)" << endl;
                for(uint i_sublvl_u = 0; i_sublvl_u < getNrOfSublevelUpper(i_species, i_trans); i_sublvl_u++)
                {
                    float sublvl_u = -getMaxMUpper(i_species, i_trans) + i_sublvl_u;
                    float sublvl_l = sublvl_u - 1;
                    if(abs(sublvl_l) <= getMaxMLower(i_species, i_trans))
                    {
                        char LineStrengthTmp[16];
#ifdef WINDOWS
                        _snprintf_s(LineStrengthTmp,
                                    sizeof(LineStrengthTmp),
                                    "%.3f",
                                    getLineStrengthSigmaM(i_species, i_trans, i_sigma_m));
#else
                        snprintf(LineStrengthTmp,
                                 sizeof(LineStrengthTmp),
                                 "%.3f",
                                 getLineStrengthSigmaM(i_species, i_trans, i_sigma_m));
#endif

                        cout << "\t" << LineStrengthTmp << "\t\t\t   " << float(sublvl_u) << "\t\t\t   "
                             << float(sublvl_l) << endl;
                        i_sigma_m++;
                    }
                }
            }
        }
    }
    cout << SEP_LINE;
}
