/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "GasMixture.hpp"
#include "CommandParser.hpp"

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

string CGasMixture::getGasSpeciesName(uint i_species)
{
    return single_species[i_species].getGasSpeciesName();
}

bool CGasMixture::updateLevelPopulation(CGridBasic * grid, photon_package * pp, uint i_species, double * J_total)
{
    return updateLevelPopulation(grid, pp->getPositionCell(), i_species, J_total);
}

bool CGasMixture::updateZeemanLevelPopulation(CGridBasic * grid,
                                              cell_basic * cell,
                                              uint i_species,
                                              uint i_line,
                                              double * sublvl_fraction)
{
    return true;
}

bool CGasMixture::updateZeemanLevelPopulation(CGridBasic * grid,
                                              photon_package * pp,
                                              uint i_species,
                                              uint i_line,
                                              double * sublvl_fraction)
{
    return updateZeemanLevelPopulation(grid, pp->getPositionCell(), i_species, i_line, sublvl_fraction);
}

bool CGasMixture::isLineZeemanSplit(uint i_species, uint i_line)
{
    return single_species[i_species].isLineZeemanSplit(i_line);
}

uint CGasMixture::getZeemanSplitIndex(uint i_species, uint i_trans)
{
    return single_species[i_species].getZeemanSplitIndex(i_trans);
}

bool CGasMixture::isTransZeemanSplit(uint i_species, uint i_trans)
{
    return single_species[i_species].isTransZeemanSplit(i_trans);
}

bool CGasMixture::isZeemanSplit(uint i_species)
{
    return single_species[i_species].isZeemanSplit();
}

bool CGasMixture::hasZeemanLines()
{
    for(uint i_species = 0; i_species < nr_of_species; i_species++)
    {
        if(single_species[i_species].isZeemanSplit())
            return true;
    }
    
    return false;
}

uilist CGasMixture::getUniqueTransitions(uint i_species)
{
    return single_species[i_species].getUniqueTransitions();
}

uint CGasMixture::getNrOfSpecies()
{
    return nr_of_species;
}

uint CGasMixture::getUniqueTransitions(uint i_species, uint i)
{
    return single_species[i_species].getUniqueTransitions(i);
}

uint CGasMixture::getLevelPopType(uint i_species)
{
    return single_species[i_species].getLevelPopType();
}

uint CGasMixture::getNrOfSpectralLines(uint i_species)
{
    return single_species[i_species].getNrOfSpectralLines();
}

uint CGasMixture::getNrOfEnergyLevels(uint i_species)
{
    return single_species[i_species].getNrOfEnergyLevels();
}

uint CGasMixture::getUpperEnergyLevel(uint i_species, uint i_trans)
{
    return single_species[i_species].getUpperEnergyLevel(i_trans);
}

uint CGasMixture::getLowerEnergyLevel(uint i_species, uint i_trans)
{
    return single_species[i_species].getLowerEnergyLevel(i_trans);
}

int CGasMixture::getNrOfSublevelUpper(uint i_species, uint i_trans)
{
    return single_species[i_species].getNrOfSublevelUpper(i_trans);
}

int CGasMixture::getNrOfSublevelLower(uint i_species, uint i_trans)
{
    return single_species[i_species].getNrOfSublevelLower(i_trans);
}

float CGasMixture::getMaxMUpper(uint i_species, uint i_trans)
{
    return single_species[i_species].getMaxMUpper(i_trans);
}

float CGasMixture::getMaxMLower(uint i_species, uint i_trans)
{
    return single_species[i_species].getMaxMLower(i_trans);
}

float CGasMixture::getMaxM(uint i_species, uint i_lvl)
{
    return single_species[i_species].getMaxM(i_lvl);
}

double CGasMixture::getAbundance(uint i_species)
{
    return single_species[i_species].getAbundance();
}

double CGasMixture::getLandeUpper(uint i_species, uint i_trans)
{
    return single_species[i_species].getLandeUpper(i_trans);
}

double CGasMixture::getLandeLower(uint i_species, uint i_trans)
{
    return single_species[i_species].getLandeLower(i_trans);
}

double CGasMixture::getLande(uint i_species, uint i_lvl)
{
    return single_species[i_species].getLande(i_lvl);
}

double CGasMixture::getCollisionRadius(uint i_species)
{
    return single_species[i_species].getCollisionRadius();
}

double CGasMixture::getLineStrength(uint i_species, uint i_trans, uint i_sublvl_u, uint i_sublvl_l)
{
    return single_species[i_species].getLineStrength(i_trans, i_sublvl_u, i_sublvl_l);
}

double CGasMixture::getMolecularWeight(uint i_species)
{
    return single_species[i_species].getMolecularWeight();
}

double CGasMixture::getTransitionFrequency(uint i_species, uint i_trans)
{
    // tr is transition from 1 to ...
    return single_species[i_species].getTransitionFrequency(i_trans);
}

double CGasMixture::getSpectralLineFrequency(uint i_species, uint i_line)
{
    uint i_trans = single_species[i_species].getTransitionFromSpectralLine(i_line);
    return single_species[i_species].getTransitionFrequency(i_trans);
}

uint CGasMixture::getNrOfTransitions(uint i_species)
{
    return single_species[i_species].getNrOfTransitions();
}

uint CGasMixture::getNrOfTotalTransitions(uint i_species)
{
    return single_species[i_species].getNrOfTotalTransitions();
}

uint CGasMixture::getTransitionFromSpectralLine(uint i_species, uint i_line)
{
    return single_species[i_species].getTransitionFromSpectralLine(i_line);
}

double CGasMixture::getGaussA(uint i_species, double temp_gas, double v_turb)
{
    return single_species[i_species].getGaussA(temp_gas, v_turb);
}

double CGasMixture::getKeplerStarMass() const
{
    return single_species[0].getKeplerStarMass();
}

double CGasMixture::getNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_species) const
{
    return single_species[i_species].getNumberDensity(grid, pp);
}

double CGasMixture::getNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_species) const
{
    return single_species[i_species].getNumberDensity(grid, cell);
}

double CGasMixture::getMassDensity(CGridBasic * grid, const photon_package & pp, uint i_species) const
{
    return single_species[i_species].getMassDensity(grid, pp);
}

double CGasMixture::getMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_species) const
{
    return single_species[i_species].getMassDensity(grid, cell);
}

void CGasMixture::setKeplerStarMass(double val)
{
    for(uint i_species = 0; i_species < nr_of_species; i_species++)
        single_species[i_species].setKeplerStarMass(val);
}

void CGasMixture::calcLineBroadening(CGridBasic * grid, uint i_species)
{
    single_species[i_species].calcLineBroadening(grid);
}

void CGasMixture::applyRadiationFieldFactor(uint i_species,
                                uint i_trans,
                                double sin_theta,
                                double cos_theta,
                                double energy,
                                double * J_nu) const
{
    single_species[i_species].applyRadiationFieldFactor(i_trans, sin_theta, cos_theta, energy, J_nu);
}

void CGasMixture::calcEmissivity(CGridBasic * grid,
                    const photon_package & pp,
                    uint i_species,
                    uint i_trans,
                    double velocity,
                    const LineBroadening & line_broadening,
                    const MagFieldInfo & mag_field_info,
                    StokesVector * line_emissivity,
                    Matrix2D * line_absorption_matrix) const
{
    single_species[i_species].calcEmissivity(grid,
                                                pp,
                                                i_trans,
                                                velocity,
                                                line_broadening,
                                                mag_field_info,
                                                line_emissivity,
                                                line_absorption_matrix);
}

uint CGasMixture::getNrOfSublevel(uint i_species, uint i_lvl)
{
    return single_species[i_species].getNrOfSublevel(i_lvl);
}

uint CGasMixture::getTotalNrOfSpectralLines()
{
    uint total_spectral_lines = 0;
    for(uint i_species = 0; i_species < nr_of_species; i_species++)
        total_spectral_lines += getNrOfSpectralLines(i_species);
    return total_spectral_lines;
}

uint CGasMixture::getMaxNrOfSublevel(uint i_species)
{
    uint max_zeeman_sublevels = 0;
    for(uint i_lvl = 0; i_lvl < getNrOfEnergyLevels(i_species); i_lvl++)
    {
        uint nr_zeeman_sublevels = getNrOfSublevel(i_species, i_lvl);

        if(nr_zeeman_sublevels > max_zeeman_sublevels)
            max_zeeman_sublevels = nr_zeeman_sublevels;
    }
    return max_zeeman_sublevels;
}

uint CGasMixture::getUniqueLevelIndex(uint i_species, uint i_lvl, uint i_sublvl)
{
    return single_species[i_species].getUniqueLevelIndex(i_lvl, i_sublvl);
}

uint CGasMixture::getNrOffsetEntries(CGridBasic * grid, parameters & param)
{
    // Init variables and pointer arrays
    uint offset_entries = 0;
    uint zeeman_sublvl_offset = 0;

    // 1x Gauss_a + voigt_a for each spectral line to simulate
    uint line_broadening_offset = 1;

    // Arrays to link energy levels, simulated spectral lines and position in the grid cells
    level_to_pos = new uint **[nr_of_species];
    line_to_pos = new uint ***[nr_of_species];

    for(uint i_species = 0; i_species < nr_of_species; i_species++)
    {
        // Increase the line broadening offset for each Zeeman split spectral line
        for(uint i_trans = 0; i_trans < getNrOfTransitions(i_species); i_trans++)
        {
            if(isTransZeemanSplit(i_species, i_trans))
                line_broadening_offset++;
        }

        uint nr_of_energy_level = getNrOfEnergyLevels(i_species);
        uint nr_of_spectral_lines = getNrOfSpectralLines(i_species);

        level_to_pos[i_species] = new uint *[nr_of_energy_level];
        line_to_pos[i_species] = new uint **[nr_of_spectral_lines];

        // Two entries for upper and lower level of a certain spectral line
        for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
        {
            line_to_pos[i_species][i_line] = new uint *[2];

            uint i_trans = getTransitionFromSpectralLine(i_species, i_line);

            uint i_lvl_l = getLowerEnergyLevel(i_species, i_trans);
            line_to_pos[i_species][i_line][0] = new uint[getNrOfSublevel(i_species, i_lvl_l)];

            uint i_lvl_u = getUpperEnergyLevel(i_species, i_trans);
            line_to_pos[i_species][i_line][1] = new uint[getNrOfSublevel(i_species, i_lvl_u)];
        }

        uint used_level_populations = 0;
        for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
        {
            level_to_pos[i_species][i_lvl] = new uint[getMaxNrOfSublevel(i_species)];

            uint nr_of_sublevel = getNrOfSublevel(i_species, i_lvl);
            for(uint i_sublvl = 0; i_sublvl < nr_of_sublevel; i_sublvl++)
            {
                // Initial set to MAX_UINT
                level_to_pos[i_species][i_lvl][i_sublvl] = MAX_UINT;

                // Found a spectral line using the energy level
                bool found = false;

                // Add all energy levels to grid for MC level pop calculation
                if(param.isGasSpeciesLevelPopMC())
                {
                    level_to_pos[i_species][i_lvl][i_sublvl] =
                        line_broadening_offset + used_level_populations;
                    found = true;
                }

                for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
                {
                    uint i_trans = getTransitionFromSpectralLine(i_species, i_line);

                    if(i_lvl == getLowerEnergyLevel(i_species, i_trans))
                    {
                        level_to_pos[i_species][i_lvl][i_sublvl] =
                            line_broadening_offset + used_level_populations;
                        line_to_pos[i_species][i_line][0][i_sublvl] =
                            line_broadening_offset + used_level_populations;
                        found = true;
                    }
                    else if(i_lvl == getUpperEnergyLevel(i_species, i_trans))
                    {
                        level_to_pos[i_species][i_lvl][i_sublvl] =
                            line_broadening_offset + used_level_populations;
                        line_to_pos[i_species][i_line][1][i_sublvl] =
                            line_broadening_offset + used_level_populations;
                        found = true;
                    }

                    if(found)
                        used_level_populations++;
                }
            }
        }
        if(used_level_populations > offset_entries)
            offset_entries = used_level_populations;
    }

    return line_broadening_offset + offset_entries + zeeman_sublvl_offset;
}

double CGasMixture::getProjCellVelocity(CGridBasic * grid, const photon_package & pp, const Vector3D & tmp_pos) const
{
    return getCellVelocity(grid, *pp.getPositionCell(), tmp_pos) * pp.getDirection();
}

Vector3D CGasMixture::getCellVelocity(CGridBasic * grid, const photon_package & pp, const Vector3D & tmp_pos) const
{
    return getCellVelocity(grid, *pp.getPositionCell(), tmp_pos);
}

Vector3D CGasMixture::getCellVelocity(CGridBasic * grid, const cell_basic & cell, const Vector3D & tmp_pos) const
{
    return single_species[0].getCellVelocity(grid, cell, tmp_pos);
}

double CGasMixture::getProjCellVelocityInterp(const Vector3D & tmp_pos,
                                    const Vector3D & dir_map_xyz,
                                    const VelFieldInterp & vel_field_interp)
{
    return single_species[0].getProjCellVelocityInterp(tmp_pos, dir_map_xyz, vel_field_interp);
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
        {
            if(!single_species[i_species].readZeemanParamaterFile(param.getZeemanCatalog(i_species)))
                return false;
        }

        if((single_species[i_species].getLevelPopType() == POP_FEP ||
            single_species[i_species].getLevelPopType() == POP_LVG) &&
           single_species[i_species].getNrOfCollisionPartner() == 0)
        {
            cout << ERROR_LINE << "FEP and LVG level population approximations require a gas "
                    "parameters file \n"
                    "       with collisional data (e.g. from Leiden Atomic and Molecular "
                    "Database)"
                 << endl;
            return false;
        }

        // Create reference list
        single_species[i_species].initReferenceLists();
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
            if(!single_species[i_species].calcLVG(grid))
                return false;
            break;
        case POP_DEGUCHI_LVG:
            if(!single_species[i_species].calcDeguchiWatsonLVG(grid))
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

void CGasMixture::printParameters(parameters & param, CGridBasic * grid)
{
    cout << CLR_LINE;
    cout << "Gas parameters                             " << endl;
    cout << SEP_LINE;
    cout << "- Velocity field                : ";
    if(getKeplerStarMass() > 0)
    {
        cout << "kepler rotation, M_star: " << getKeplerStarMass() << " [M_sun]" << endl;
        cout << INFO_LINE << "Only available with one central star" << endl;
    }
    else if(grid->isVelocityFieldAvailable())
        cout << "taken from grid" << endl;
    else
        cout << "none" << endl;

    cout << "- Turbulent Velocity            : ";
    if(param.getTurbulentVelocity() > 0)
        cout << param.getTurbulentVelocity() << " [m/s]" << endl;
    else if(grid->isTurbulentVelocityAvailable())
        cout << "taken from grid" << endl;
    else
        cout << "none" << endl;

    for(uint i_species = 0; i_species < nr_of_species; i_species++)
    {
        cout << SEP_LINE;
        cout << "Gas species " << (i_species + 1) << " (" << getGasSpeciesName(i_species) << ")" << endl;

        if(single_species[i_species].getNrOfSpectralLines() == 0)
        {
            cout << WARNING_LINE << "No spectral lines selected!" << endl;
            return;
        }

        stringstream transition_str, vel_channels_str, min_vel_str, max_vel_str;
        dlist line_ray_detectors = param.getLineRayDetector(i_species);
        for(uint i = 0; i < line_ray_detectors.size(); i += NR_OF_LINE_DET)
        {
            transition_str << uint(line_ray_detectors[i] + 1);
            min_vel_str << line_ray_detectors[i + 2];
            max_vel_str << line_ray_detectors[i + 3];
            vel_channels_str << uint(line_ray_detectors[i + NR_OF_LINE_DET - 1]);
            if(i < line_ray_detectors.size() - NR_OF_LINE_DET)
            {
                transition_str << ", ";
                min_vel_str << ", ";
                max_vel_str << ", ";
                vel_channels_str << ", ";
            }
        }
        cout << "- Line transition(s)            : " << transition_str.str() << endl;
        cout << "- Number of velocity channels   : " << vel_channels_str.str() << endl;
        cout << "- Min. Velocity limit(s)        : " << min_vel_str.str() << " [m/s]" << endl;
        cout << "- Max. Velocity limit(s)        : " << max_vel_str.str() << " [m/s]" << endl;
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
            case POP_DEGUCHI_LVG:
                cout << "LVG (Deguchi & Watson 1984)" << endl;
                break;
            default:
                cout << ERROR_LINE << "UNKNOWN!" << endl;
        }

        if(isZeemanSplit(i_species))
            cout << "- Particle radius (collisions)  : " << getCollisionRadius(i_species) << " [m]" << endl;

        cout << "- Molecular weight              : " << getMolecularWeight(i_species) << endl;
        double ab = getAbundance(i_species);
        if(ab > 0)
            cout << "- Abundance                     : " << ab << endl;
        else
        {
            double min_dens_species = 1e200, max_dens_species = 0;
            for(ulong i_cell = 0; i_cell < grid->getMaxDataCells(); i_cell++)
            {
                const cell_basic & cell = *grid->getCellFromIndex(i_cell);

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
        for(ulong i_cell = 0; i_cell < grid->getMaxDataCells(); i_cell++)
        {
            cell_basic * cell = grid->getCellFromIndex(i_cell);
            total_species_mass += getMassDensity(grid, *cell, i_species) * grid->getVolume(*cell);
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

                cout << "- Transitions between sublevels : " << endl;
                cout << "    transition type\tline strength\t\tm(upper)\t\tm(lower)" << endl;
                for(uint i_sublvl_u = 0; i_sublvl_u < getNrOfSublevelUpper(i_species, i_trans); i_sublvl_u++)
                {
                    // Calculate the quantum number of the upper energy level
                    float sublvl_u = -getMaxMUpper(i_species, i_trans) + i_sublvl_u;

                    for(uint i_sublvl_l = 0; i_sublvl_l < getNrOfSublevelLower(i_species, i_trans);
                        i_sublvl_l++)
                    {
                        // Calculate the quantum number of the lower energy level
                        float sublvl_l = -getMaxMLower(i_species, i_trans) + i_sublvl_l;

                        char LineStrengthTmp[16];
#ifdef WINDOWS
                        _snprintf_s(LineStrengthTmp,
                                    sizeof(LineStrengthTmp),
                                    "%.3f",
                                    getLineStrength(i_species, i_trans, i_sublvl_u, i_sublvl_l));
#else
                        snprintf(LineStrengthTmp,
                                 sizeof(LineStrengthTmp),
                                 "%.3f",
                                 getLineStrength(i_species, i_trans, i_sublvl_u, i_sublvl_l));
#endif

                        // Use the correct propagation matrix and relative line strength that
                        // depends on the current type of Zeeman transition (pi, sigma_-, sigma_+)
                        switch(int(sublvl_l - sublvl_u))
                        {
                            case TRANS_SIGMA_P:
                                cout << "\tSigma+\t\t    " << LineStrengthTmp << "\t\t   " << float(sublvl_u)
                                     << "\t\t\t   " << float(sublvl_l) << endl;
                                break;
                            case TRANS_PI:
                                cout << "\tPi    \t\t    " << LineStrengthTmp << "\t\t   " << float(sublvl_u)
                                     << "\t\t\t   " << float(sublvl_l) << endl;
                                break;
                            case TRANS_SIGMA_M:
                                cout << "\tSigma-\t\t    " << LineStrengthTmp << "\t\t   " << float(sublvl_u)
                                     << "\t\t\t   " << float(sublvl_l) << endl;
                                break;
                        }
                    }
                }
            }
        }
    }
    cout << SEP_LINE;
}
