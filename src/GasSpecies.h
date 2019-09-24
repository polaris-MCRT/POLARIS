#pragma once
#include "Faddeeva.hh"
#include "Grid.h"
#include "Matrix2D.h"
#include "Typedefs.h"
#include <complex>

#ifndef CMOLECULE
#define CMOLECULE

// Additional Structure
struct VelFieldInterp
{
    spline vel_field;
    bool zero_vel_field;
    Vector3D start_pos;
};

class CGasSpecies
{
  public:
    CGasSpecies()
    {
        molecular_weight = 0;
        lvl_pop_type = POP_LTE;

        nr_of_energy_level = 0;
        nr_of_transitions = 0;
        nr_of_col_partner = 0;
        nr_of_spectral_lines = 0;

        abundance = 0;

        gas_species_radius = 0;

        kepler_star_mass = 0;

        nr_zeeman_spectral_lines = 0;

        quantum_numbers = 0;
        g_level = 0;
        energy_level = 0;
        nr_of_col_transition = 0;

        spectral_lines = 0;

        level_to_index = 0;
        trans_to_index = 0;

        upper_level = 0;
        lower_level = 0;

        trans_einstA = 0;
        trans_einstB_lu = 0;
        trans_einstB_ul = 0;

        trans_freq = 0;
        trans_inner_energy = 0;

        trans_is_zeeman_split = 0;

        orientation_H2 = 0;
        nr_of_col_temp = 0;

        lande_factor = 0;
        nr_of_sublevel = 0;

        collision_temp = 0;

        col_upper = 0;
        col_lower = 0;

        col_matrix = 0;

        stringID = "";

        i_species = 0;
    }

    ~CGasSpecies()
    {
        if(col_matrix != 0)
        {
            for(uint i = 0; i < nr_of_col_partner; i++)
            {
                for(int j = 0; j < nr_of_col_transition[i]; j++)
                    delete[] col_matrix[i][j];
                delete[] col_matrix[i];
            }
            delete[] col_matrix;
        }

        if(quantum_numbers != 0)
            delete[] quantum_numbers;

        if(g_level != 0)
            delete[] g_level;

        if(spectral_lines != 0)
            delete[] spectral_lines;

        if(level_to_index != 0)
        {
            for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
                delete[] level_to_index[i_lvl];
            delete[] level_to_index;
        }

        if(trans_to_index != 0)
        {
            for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
            {
                uint nr_of_sublevel_u = nr_of_sublevel[upper_level[i_trans]];
                for(uint i_sublvl_u = 0; i_sublvl_u < nr_of_sublevel_u; i_sublvl_u++)
                    delete[] trans_to_index[i_trans][i_sublvl_u];
                delete[] trans_to_index[i_trans];
            }
            delete[] trans_to_index;
        }

        if(energy_level != 0)
            delete[] energy_level;

        if(nr_of_col_transition != 0)
            delete[] nr_of_col_transition;

        if(upper_level != 0)
            delete[] upper_level;

        if(lower_level != 0)
            delete[] lower_level;

        if(trans_einstA != 0)
        {
            for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
                delete[] trans_einstA[i_trans];
            delete[] trans_einstA;
        }

        if(trans_einstB_lu != 0)
        {
            for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
                delete[] trans_einstB_lu[i_trans];
            delete[] trans_einstB_lu;
        }

        if(trans_einstB_ul != 0)
        {
            for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
                delete[] trans_einstB_ul[i_trans];
            delete[] trans_einstB_ul;
        }

        if(orientation_H2 != 0)
            delete[] orientation_H2;

        if(nr_of_sublevel != 0)
            delete[] nr_of_sublevel;

        if(nr_of_col_temp != 0)
            delete[] nr_of_col_temp;

        if(collision_temp != 0)
        {
            for(uint i = 0; i < nr_of_col_partner; i++)
                delete[] collision_temp[i];
            delete[] collision_temp;
        }

        if(col_upper != 0)
        {
            for(uint i = 0; i < nr_of_col_partner; i++)
                delete[] col_upper[i];
            delete[] col_upper;
        }

        if(col_lower != 0)
        {
            for(uint i = 0; i < nr_of_col_partner; i++)
                delete[] col_lower[i];
            delete[] col_lower;
        }
    }

    string getGasSpeciesName() const
    {
        return stringID;
    }

    double getLineStrength(uint i_trans, uint i_sublvl_u, uint i_sublvl_l) const
    {
        uint i_sublvl = getSublevelIndex(i_trans, i_sublvl_u, i_sublvl_l);
        if(getEinsteinA(i_trans) > 0)
            return getEinsteinA(i_trans, i_sublvl_u, i_sublvl_l) / getEinsteinA(i_trans);
        return 0;
    }

    double getAbundance() const
    {
        return abundance;
    }

    double getKeplerStarMass() const
    {
        return kepler_star_mass;
    }

    int getTransitionFromSpectralLine(uint i_line) const
    {
        return spectral_lines[i_line];
    }

    double getTransitionFrequency(uint i_trans) const
    {
        return trans_freq[i_trans];
    }

    double getSpectralLineFrequency(uint i_line) const
    {
        uint i_trans = getTransitionFromSpectralLine(i_line);
        return getTransitionFrequency(i_trans);
    }

    uint getNrOfSpectralLines() const
    {
        return nr_of_spectral_lines;
    }

    double getEinsteinA(uint i_trans) const
    {
        return trans_einstA[i_trans][0];
    }

    double getEinsteinBul(uint i_trans) const
    {
        return trans_einstB_ul[i_trans][0];
    }

    double getEinsteinBlu(uint i_trans) const
    {
        return trans_einstB_lu[i_trans][0];
    }

    uint getSublevelIndex(uint i_trans, uint i_sublvl_u, uint i_sublvl_l) const
    {
        return i_sublvl_u * getNrOfSublevelLower(i_trans) + i_sublvl_l;
    }

    double getEinsteinA(uint i_trans, uint i_sublvl_u, uint i_sublvl_l) const
    {
        // If not Zeeman split, use total value
        if(!isTransZeemanSplit(i_trans))
            return getEinsteinA(i_trans) / getNrOfSublevelLower(i_trans);

        uint i_sublvl = getSublevelIndex(i_trans, i_sublvl_u, i_sublvl_l);
        return trans_einstA[i_trans][i_sublvl + 1];
    }

    double getEinsteinBul(uint i_trans, uint i_sublvl_u, uint i_sublvl_l) const
    {
        // If not Zeeman split, use total value
        if(!isTransZeemanSplit(i_trans))
            return getEinsteinBul(i_trans) / getNrOfSublevelLower(i_trans);

        uint i_sublvl = getSublevelIndex(i_trans, i_sublvl_u, i_sublvl_l);
        return trans_einstB_ul[i_trans][i_sublvl + 1];
    }

    double getEinsteinBlu(uint i_trans, uint i_sublvl_u, uint i_sublvl_l) const
    {
        // If not Zeeman split, use total value
        if(!isTransZeemanSplit(i_trans))
            return getEinsteinBlu(i_trans) / getNrOfSublevelUpper(i_trans);

        uint i_sublvl = getSublevelIndex(i_trans, i_sublvl_u, i_sublvl_l);
        return trans_einstB_lu[i_trans][i_sublvl + 1];
    }

    double getJLevel(uint i_lvl) const
    {
        return quantum_numbers[i_lvl];
    }

    double getEnergyOfLevel(uint i_lvl) const
    {
        return energy_level[i_lvl];
    }

    int getUpperCollisionLevel(uint m, uint n) const
    {
        return col_upper[m][n];
    }

    int getLowerCollisionLevel(uint m, uint n) const
    {
        return col_lower[m][n];
    }

    double getCollisionTemp(uint m, uint n) const
    {
        return collision_temp[m][n];
    }

    double getCollisionMatrix(uint m, uint n, uint k) const
    {
        return col_matrix[m][n][k];
    }

    double getLandeUpper(uint i_trans) const
    {
        return lande_factor[upper_level[i_trans]];
    }

    double getLandeLower(uint i_trans) const
    {
        return lande_factor[lower_level[i_trans]];
    }

    double getLande(uint i_lvl) const
    {
        return lande_factor[i_lvl];
    }

    double getCollisionRadius() const
    {
        return gas_species_radius;
    }

    double getMolecularWeight() const
    {
        return molecular_weight;
    }

    double getNumberDensity(CGridBasic * grid, const photon_package & pp) const
    {
        return getNumberDensity(grid, *pp.getPositionCell());
    }

    double getNumberDensity(CGridBasic * grid, const cell_basic & cell) const
    {
        double dens_species = 0;
        // If the abundance is negative, get its value from the grid
        if(abundance <= 0)
        {
            uint fr_id = uint(-abundance - 1);
            dens_species = grid->getCellAbundance(cell, fr_id);
        }
        else
            dens_species = abundance;
        dens_species *= grid->getGasNumberDensity(cell);
        return dens_species;
    }

    double getMassDensity(CGridBasic * grid, const photon_package & pp) const
    {
        return getMassDensity(grid, *pp.getPositionCell());
    }

    double getMassDensity(CGridBasic * grid, const cell_basic & cell) const
    {
        double dens_species = 0;
        // If the abundance is negative, get its value from the grid
        if(abundance <= 0)
        {
            uint fr_id = uint(-abundance - 1);
            dens_species = grid->getCellAbundance(cell, fr_id);
        }
        else
            dens_species = abundance;
        dens_species *= grid->getGasNumberDensity(cell);
        dens_species *= molecular_weight * m_H;
        return dens_species;
    }

    Vector3D getCellVelocity(CGridBasic * grid, const cell_basic & cell, const Vector3D & tmp_pos) const
    {
        Vector3D cell_velocity;

        // Get the velocity in the photon
        // direction of the current position
        if(kepler_star_mass > 0)
        {
            // Get velocity from Kepler rotation
            cell_velocity = CMathFunctions::calcKeplerianVelocity(tmp_pos, kepler_star_mass);
        }
        else if(grid->hasVelocityField())
        {
            // Get velocity from grid cell
            cell_velocity = grid->getVelocityField(cell);
        }
        return cell_velocity;
    }

    double getProjCellVelocityInterp(const Vector3D & tmp_pos,
                                     const Vector3D & dir_map_xyz,
                                     const VelFieldInterp & vel_field_interp)
    {
        double cell_velocity = 0;

        // Get the velocity in the photon direction of the current position
        if(kepler_star_mass > 0)
        {
            // Get velocity from Kepler rotation
            cell_velocity = CMathFunctions::calcKeplerianVelocity(tmp_pos, kepler_star_mass) * dir_map_xyz;
        }
        else if(vel_field_interp.vel_field.size() > 0 && !vel_field_interp.zero_vel_field)
        {
            // Get velocity from grid cell with interpolation
            Vector3D rel_pos = tmp_pos - vel_field_interp.start_pos;
            cell_velocity = vel_field_interp.vel_field.getValue(rel_pos.length());
        }
        return cell_velocity;
    }

    void initReferenceLists()
    {
        // Init first dimension of 2D array
        level_to_index = new uint *[nr_of_energy_level];

        uint i_lvl_unique = 0;
        for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
        {
            // Init second dimension of 2D array
            level_to_index[i_lvl] = new uint[nr_of_sublevel[i_lvl]];

            for(uint i_sublvl = 0; i_sublvl < nr_of_sublevel[i_lvl]; i_sublvl++)
            {
                level_to_index[i_lvl][i_sublvl] = i_lvl_unique;
                i_lvl_unique++;
            }
        }

        // Init first dimension of 2D array
        trans_to_index = new uint **[nr_of_transitions];

        uint i_trans_unique = 0;
        for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
        {
            // Get all indices and number of sublevels
            uint i_lvl_u = getUpperEnergyLevel(i_trans);
            uint i_lvl_l = getLowerEnergyLevel(i_trans);
            uint nr_of_sublevel_u = nr_of_sublevel[i_lvl_u];
            uint nr_of_sublevel_l = nr_of_sublevel[i_lvl_l];

            // Init second dimension of 2D array
            trans_to_index[i_trans] = new uint *[nr_of_sublevel_u];

            for(uint i_sublvl_u = 0; i_sublvl_u < nr_of_sublevel_u; i_sublvl_u++)
            {
                // Init third dimension of 2D array
                trans_to_index[i_trans][i_sublvl_u] = new uint[nr_of_sublevel_l];

                for(uint i_sublvl_l = 0; i_sublvl_l < nr_of_sublevel_l; i_sublvl_l++)
                {
                    trans_to_index[i_trans][i_sublvl_u][i_sublvl_l] = i_trans_unique;
                    i_trans_unique++;
                }
            }
        }
    }

    uint getZeemanSplitIndex(uint i_trans)
    {
        if(trans_is_zeeman_split != 0)
        {
            uint i_zeeman = 0;
            for(uint i = 0; i < i_trans; i++)
            {
                if(trans_is_zeeman_split[i])
                    i_zeeman++;
            }
            if(trans_is_zeeman_split[i_trans])
                return i_zeeman;
            else
                return MAX_UINT;
        }
        return MAX_UINT;
    }

    bool isTransZeemanSplit(uint i_trans) const
    {
        if(trans_is_zeeman_split != 0)
            return trans_is_zeeman_split[i_trans];
        return false;
    }

    bool isLineZeemanSplit(uint i_line) const
    {
        uint i_trans = getTransitionFromSpectralLine(i_line);
        return isTransZeemanSplit(i_trans);
    }

    bool isZeemanSplit() const
    {
        for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
        {
            if(isTransZeemanSplit(i_trans))
                return true;
        }
        return false;
    }

    uint getUniqueLevelIndex(uint i_lvl, uint i_sublvl) const
    {
        return level_to_index[i_lvl][i_sublvl];
    }

    uint getUniqueTransitionIndex(uint i_trans, uint i_sublvl_u = 0, uint i_sublvl_l = 0) const
    {
        return trans_to_index[i_trans][i_sublvl_u][i_sublvl_l];
    }

    uilist getUniqueTransitions() const
    {
        return unique_spectral_lines;
    }

    uint getUniqueTransitions(uint i) const
    {
        return unique_spectral_lines[i];
    }

    uint getUpperEnergyLevel(uint i_trans) const
    {
        return upper_level[i_trans];
    }

    uint getLowerEnergyLevel(uint i_trans) const
    {
        return lower_level[i_trans];
    }

    uint getNrOfEnergyLevels() const
    {
        return nr_of_energy_level;
    }

    uint getNrOfTotalEnergyLevels() const
    /*
    Including Zeeman sublevels
    */
    {
        uint res = 0;
        for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
        {
            // Includes the 1 for no Zeeman split energy levels
            res += nr_of_sublevel[i_lvl];
        }
        return res;
    }

    uint getNrOfTotalTransitions() const
    /*
    Including Zeeman sublevels
    */
    {
        uint res = 0;
        for(uint i_trans = 0; i_trans < nr_of_transitions; i_trans++)
        {
            uint i_lvl_u = getUpperEnergyLevel(i_trans);
            uint i_lvl_l = getLowerEnergyLevel(i_trans);

            // Includes the 1 for no Zeeman split energy levels
            res += nr_of_sublevel[i_lvl_u] * nr_of_sublevel[i_lvl_l];
        }
        return res;
    }

    uint getNrOfTransitions() const
    {
        return nr_of_transitions;
    }

    uint getNrOfCollisionPartner() const
    {
        return nr_of_col_partner;
    }

    uint getNrCollisionTransitions(uint i_col_partner) const
    {
        return nr_of_col_transition[i_col_partner];
    }

    uint getNrCollisionTemps(uint i_col_partner) const
    {
        return nr_of_col_temp[i_col_partner];
    }

    uint getLevelPopType() const
    {
        return lvl_pop_type;
    }

    uint getOrientation_H2(uint i_col_partner) const
    {
        return orientation_H2[i_col_partner];
    }

    int getNrOfSublevelUpper(uint i_trans) const
    {
        return nr_of_sublevel[upper_level[i_trans]];
    }

    int getNrOfSublevelLower(uint i_trans) const
    {
        return nr_of_sublevel[lower_level[i_trans]];
    }

    int getNrOfSublevel(uint i_lvl) const
    {
        return nr_of_sublevel[i_lvl];
    }

    int getNrOfTransBetweenSublevels(uint i_trans) const
    {
        return getNrOfSublevelUpper(i_trans) * getNrOfSublevelLower(i_trans);
    }

    float getMaxMUpper(uint i_trans) const
    {
        return float((getNrOfSublevelUpper(i_trans) - 1) / 2.0);
    }

    float getMaxMLower(uint i_trans) const
    {
        return float((getNrOfSublevelLower(i_trans) - 1) / 2.0);
    }

    float getMaxM(uint i_lvl) const
    {
        return float((nr_of_sublevel[i_lvl] - 1) / 2.0);
    }

    void setKeplerStarMass(double val)
    {
        kepler_star_mass = val;
    }

    void setLevelPopType(uint type)
    {
        lvl_pop_type = type;
    }

    void setMolecularWeight(double val)
    {
        molecular_weight = val;
    }

    void setAbundance(double val)
    {
        abundance = val;
    }

    void setNrOfSpectralLines(uint val)
    {
        nr_of_spectral_lines = val;
    }

    void setSpectralLines(int * lines)
    {
        spectral_lines = lines;
    }

    double getGamma(uint i_trans, double dens_gas, double dens_species, double temp_gas, double v_turb)
    {
        double gamma = getEinsteinA(i_trans);

        // "http://chem.libretexts.org/Core/Physical_and_Theoretical_Chemistry/
        //  ->
        //  Kinetics/Modeling_Reaction_Kinetics/Collision_Theory/Collisional_Cross_Section"
        // "http://www.phy.ohiou.edu/~mboett/astro401_fall12/broadening.pdf
        double v_th = sqrt(2.0 * con_kB * temp_gas / (molecular_weight * 1e-3 / con_Na));
        double col_param =
            dens_gas * PI * pow(con_r_bohr + gas_species_radius, 2) * sqrt(pow(v_th, 2) + pow(v_turb, 2));

        return gamma + 2 * col_param;
    }

    double getGaussA(double temp_gas, double v_turb)
    {
        double v_th = sqrt(2.0 * con_kB * temp_gas / (molecular_weight * 1e-3 / con_Na));
        double gauss_a = 1.0 / sqrt(pow(v_th, 2) + pow(v_turb, 2));
        return gauss_a;
    }

    void calcLineBroadening(CGridBasic * grid);

    static complex<double> getLineShape_AB(double f_doppler, double a)
    {
        // Faddeeva::w routine:
        // Copyright © 2012 Massachusetts Institute of Technology
        // Permission is hereby granted, free of charge, to any person obtaining a copy
        // of this software and associated documentation files (the "Software"), to deal
        // in the Software without restriction, including without limitation the rights
        // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        // copies of the Software, and to permit persons to whom the Software is
        // furnished to do so, subject to the following conditions: The above copyright
        // notice and this permission notice shall be included in all copies or
        // substantial portions of the Software.

        complex<double> z(f_doppler, a);
        complex<double> w = Faddeeva::w(z, 0);

        return w;
    }

    void getGaussLineMatrix(CGridBasic * grid,
                            const cell_basic & cell,
                            double velocity,
                            Matrix2D * line_absorption_matrix) const
    {

        // Calculate gaussian shape
        double line_amplitude = getGaussLineShape(grid, cell, velocity);

        // Only diagonal without polarization rotation matrix elements
        for(uint i = 0; i < 4; i++)
            line_absorption_matrix->setValue(i, i, line_amplitude);
    }

    void getGaussLineMatrix(CGridBasic * grid,
                            const photon_package & pp,
                            double velocity,
                            Matrix2D * line_absorption_matrix) const
    {
        getGaussLineMatrix(grid, *pp.getPositionCell(), velocity, line_absorption_matrix);
    }

    double getGaussLineShape(CGridBasic * grid, const cell_basic & cell, double velocity) const
    {
        double gauss_a = grid->getGaussA(cell);
        return exp(-(pow(velocity, 2) * pow(gauss_a, 2))) / PIsq;
    }

    double getGaussLineShape(CGridBasic * grid, const photon_package & pp, double velocity) const
    {
        return getGaussLineShape(grid, *pp.getPositionCell(), velocity);
    }

    void calcEmissivityZeeman(CGridBasic * grid,
                              const photon_package & pp,
                              uint i_line,
                              double velocity,
                              const LineBroadening & line_broadening,
                              const MagFieldInfo & mfo,
                              StokesVector * line_emissivity,
                              Matrix2D * line_absorption_matrix) const;

    void calcEmissivity(CGridBasic * grid,
                        const photon_package & pp,
                        uint i_line,
                        double velocity,
                        const LineBroadening & line_broadening,
                        const MagFieldInfo & mag_field_info,
                        StokesVector * line_emissivity,
                        Matrix2D * line_absorption_matrix) const;

    bool calcLTE(CGridBasic * grid, bool full = false);
    bool calcFEP(CGridBasic * grid, bool full = false);
    bool calcLVG(CGridBasic * grid, bool full = false);
    bool calcDeguchiWatsonLVG(CGridBasic * grid, bool full = false);
    bool updateLevelPopulation(CGridBasic * grid, cell_basic * cell, double * J_total);

    double elem_LVG(CGridBasic * grid,
                    double dens_species,
                    double * tmp_lvl_pop,
                    double doppler,
                    double L,
                    double J_ext,
                    uint i_trans,
                    uint i_sublvl_u,
                    uint i_sublvl_l);

    void createMatrix(double * J_mid, Matrix2D * A, double * b, double *** final_col_para);

    double *** calcCollisionParameter(CGridBasic * grid, cell_basic * cell);
    double getColPartnerDensity(CGridBasic * grid, cell_basic * cell, uint i_col_partner);
    dlist calcCollisionRate(uint i_col_partner, uint i_col_transition, uint hi_i, double temp_gas);

    bool readGasParamaterFile(string filename, uint id, uint max);
    bool readZeemanParamaterFile(string filename);

    void applyRadiationFieldFactor(uint i_trans,
                                   double sin_theta,
                                   double cos_theta,
                                   double energy,
                                   double * J_nu) const;

  private:
    double ** collision_temp;
    int **col_upper, **col_lower;
    double **trans_einstA, **trans_einstB_lu, **trans_einstB_ul;

    double *** col_matrix;

    prob_list frequency_prob;

    double molecular_weight;
    double abundance;
    double max_velocity;
    double gas_species_radius;
    double kepler_star_mass;

    uint i_species;

    uint nr_of_transitions;
    uint nr_of_col_partner;
    uint nr_of_energy_level;
    uint nr_of_spectral_lines;
    uint nr_zeeman_spectral_lines;
    uint lvl_pop_type;

    uint ** level_to_index;
    uint *** trans_to_index;

    uilist unique_spectral_lines;

    bool * trans_is_zeeman_split;

    int * nr_of_col_transition;
    int * nr_of_col_temp;
    int * nr_of_sublevel;
    int * spectral_lines;
    int *upper_level, *lower_level;
    int * orientation_H2;

    double * energy_level;
    double * g_level;
    double * quantum_numbers;
    double * lande_factor;
    double *trans_freq, *trans_inner_energy;

    string stringID;
    string filename;
    string catalog_path;

    ostringstream velocity_information;
};

class CGasMixture
{
  public:
    CGasMixture()
    {
        nr_of_species = 0;

        level_to_pos = 0;
        line_to_pos = 0;

        single_species = 0;
    }

    ~CGasMixture(void)
    {
        if(single_species != 0)
            delete[] single_species;

        if(level_to_pos != 0)
        {
            for(uint i_species = 0; i_species < nr_of_species; i_species++)
            {
                for(uint i_lvl = 0; i_lvl < getNrOfEnergyLevels(i_species); i_lvl++)
                    delete[] level_to_pos[i_species][i_lvl];
                delete[] level_to_pos[i_species];
            }
            delete[] level_to_pos;
        }

        if(line_to_pos != 0)
        {
            for(uint i_species = 0; i_species < nr_of_species; i_species++)
            {
                for(uint i_line = 0; i_line < getNrOfSpectralLines(i_species); i_line++)
                {
                    for(uint i = 0; i < 2; i++)
                        delete[] line_to_pos[i_species][i_line][i];
                    delete[] line_to_pos[i_species][i_line];
                }
                delete[] line_to_pos[i_species];
            }
            delete[] line_to_pos;
        }
    }

    string getGasSpeciesName(uint i_species)
    {
        return single_species[i_species].getGasSpeciesName();
    }

    bool createGasSpecies(parameters & param);

    bool calcLevelPopulation(CGridBasic * grid, uint i_species);

    bool updateLevelPopulation(CGridBasic * grid, cell_basic * cell, uint i_species, double * J_total);
    bool updateLevelPopulation(CGridBasic * grid, photon_package * pp, uint i_species, double * J_total)
    {
        return updateLevelPopulation(grid, pp->getPositionCell(), i_species, J_total);
    }

    bool updateZeemanLevelPopulation(CGridBasic * grid,
                                     cell_basic * cell,
                                     uint i_species,
                                     uint i_line,
                                     double * sublvl_fraction);
    bool updateZeemanLevelPopulation(CGridBasic * grid,
                                     photon_package * pp,
                                     uint i_species,
                                     uint i_line,
                                     double * sublvl_fraction)
    {
        return updateZeemanLevelPopulation(grid, pp->getPositionCell(), i_species, i_line, sublvl_fraction);
    }

    bool isLineZeemanSplit(uint i_species, uint i_line)
    {
        return single_species[i_species].isLineZeemanSplit(i_line);
    }

    uint getZeemanSplitIndex(uint i_species, uint i_trans)
    {
        return single_species[i_species].getZeemanSplitIndex(i_trans);
    }

    bool isTransZeemanSplit(uint i_species, uint i_trans)
    {
        return single_species[i_species].isTransZeemanSplit(i_trans);
    }

    bool isZeemanSplit(uint i_species)
    {
        return single_species[i_species].isZeemanSplit();
    }

    uilist getUniqueTransitions(uint i_species)
    {
        return single_species[i_species].getUniqueTransitions();
    }

    uint getNrOfSpecies()
    {
        return nr_of_species;
    }

    uint getUniqueTransitions(uint i_species, uint i)
    {
        return single_species[i_species].getUniqueTransitions(i);
    }

    uint getLevelPopType(uint i_species)
    {
        return single_species[i_species].getLevelPopType();
    }

    uint getNrOfSpectralLines(uint i_species)
    {
        return single_species[i_species].getNrOfSpectralLines();
    }

    uint getNrOfEnergyLevels(uint i_species)
    {
        return single_species[i_species].getNrOfEnergyLevels();
    }

    uint getUpperEnergyLevel(uint i_species, uint i_trans)
    {
        return single_species[i_species].getUpperEnergyLevel(i_trans);
    }

    uint getLowerEnergyLevel(uint i_species, uint i_trans)
    {
        return single_species[i_species].getLowerEnergyLevel(i_trans);
    }

    int getNrOfSublevelUpper(uint i_species, uint i_trans)
    {
        return single_species[i_species].getNrOfSublevelUpper(i_trans);
    }

    int getNrOfSublevelLower(uint i_species, uint i_trans)
    {
        return single_species[i_species].getNrOfSublevelLower(i_trans);
    }

    float getMaxMUpper(uint i_species, uint i_trans)
    {
        return single_species[i_species].getMaxMUpper(i_trans);
    }

    float getMaxMLower(uint i_species, uint i_trans)
    {
        return single_species[i_species].getMaxMLower(i_trans);
    }

    float getMaxM(uint i_species, uint i_lvl)
    {
        return single_species[i_species].getMaxM(i_lvl);
    }

    double getAbundance(uint i_species)
    {
        return single_species[i_species].getAbundance();
    }

    double getLandeUpper(uint i_species, uint i_trans)
    {
        return single_species[i_species].getLandeUpper(i_trans);
    }

    double getLandeLower(uint i_species, uint i_trans)
    {
        return single_species[i_species].getLandeLower(i_trans);
    }

    double getLande(uint i_species, uint i_lvl)
    {
        return single_species[i_species].getLande(i_lvl);
    }

    double getCollisionRadius(uint i_species)
    {
        return single_species[i_species].getCollisionRadius();
    }

    double getLineStrength(uint i_species, uint i_trans, uint i_sublvl_u, uint i_sublvl_l)
    {
        return single_species[i_species].getLineStrength(i_trans, i_sublvl_u, i_sublvl_l);
    }

    double getMolecularWeight(uint i_species)
    {
        return single_species[i_species].getMolecularWeight();
    }

    double getTransitionFrequency(uint i_species, uint i_trans)
    {
        // tr is transition from 1 to ...
        return single_species[i_species].getTransitionFrequency(i_trans);
    }

    double getSpectralLineFrequency(uint i_species, uint i_line)
    {
        uint i_trans = single_species[i_species].getTransitionFromSpectralLine(i_line);
        return single_species[i_species].getTransitionFrequency(i_trans);
    }

    uint getNrOfTransitions(uint i_species)
    {
        return single_species[i_species].getNrOfTransitions();
    }

    uint getNrOfTotalTransitions(uint i_species)
    {
        return single_species[i_species].getNrOfTotalTransitions();
    }

    uint getTransitionFromSpectralLine(uint i_species, uint i_line)
    {
        return single_species[i_species].getTransitionFromSpectralLine(i_line);
    }

    double getGaussA(uint i_species, double temp_gas, double v_turb)
    {
        return single_species[i_species].getGaussA(temp_gas, v_turb);
    }

    double getKeplerStarMass() const
    {
        return single_species[0].getKeplerStarMass();
    }

    double getNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_species) const
    {
        return single_species[i_species].getNumberDensity(grid, pp);
    }

    double getNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_species) const
    {
        return single_species[i_species].getNumberDensity(grid, cell);
    }

    double getMassDensity(CGridBasic * grid, const photon_package & pp, uint i_species) const
    {
        return single_species[i_species].getMassDensity(grid, pp);
    }

    double getMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_species) const
    {
        return single_species[i_species].getMassDensity(grid, cell);
    }

    void setKeplerStarMass(double val)
    {
        for(uint i_species = 0; i_species < nr_of_species; i_species++)
            single_species[i_species].setKeplerStarMass(val);
    }

    void calcLineBroadening(CGridBasic * grid, uint i_species)
    {
        single_species[i_species].calcLineBroadening(grid);
    }

    void applyRadiationFieldFactor(uint i_species,
                                   uint i_trans,
                                   double sin_theta,
                                   double cos_theta,
                                   double energy,
                                   double * J_nu) const
    {
        single_species[i_species].applyRadiationFieldFactor(i_trans, sin_theta, cos_theta, energy, J_nu);
    }

    void calcEmissivity(CGridBasic * grid,
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

    uint getNrOfSublevel(uint i_species, uint i_lvl)
    {
        return single_species[i_species].getNrOfSublevel(i_lvl);
    }

    uint getTotalNrOfSpectralLines()
    {
        uint total_spectral_lines = 0;
        for(uint i_species = 0; i_species < nr_of_species; i_species++)
            total_spectral_lines += getNrOfSpectralLines(i_species);
        return total_spectral_lines;
    }

    uint getMaxNrOfSublevel(uint i_species)
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

    uint getUniqueLevelIndex(uint i_species, uint i_lvl, uint i_sublvl)
    {
        return single_species[i_species].getUniqueLevelIndex(i_lvl, i_sublvl);
    }

    uint getNrOffsetEntries(CGridBasic * grid, parameters & param)
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

    double getProjCellVelocity(CGridBasic * grid, const photon_package & pp, const Vector3D & tmp_pos) const
    {
        return getCellVelocity(grid, *pp.getPositionCell(), tmp_pos) * pp.getDirection();
    }

    Vector3D getCellVelocity(CGridBasic * grid, const photon_package & pp, const Vector3D & tmp_pos) const
    {
        return getCellVelocity(grid, *pp.getPositionCell(), tmp_pos);
    }

    Vector3D getCellVelocity(CGridBasic * grid, const cell_basic & cell, const Vector3D & tmp_pos) const
    {
        return single_species[0].getCellVelocity(grid, cell, tmp_pos);
    }

    double getProjCellVelocityInterp(const Vector3D & tmp_pos,
                                     const Vector3D & dir_map_xyz,
                                     const VelFieldInterp & vel_field_interp)
    {
        return single_species[0].getProjCellVelocityInterp(tmp_pos, dir_map_xyz, vel_field_interp);
    }

    void printParameter(parameters & param, CGridBasic * grid);

  private:
    CGasSpecies * single_species;

    uint *** level_to_pos;
    uint **** line_to_pos;

    uint nr_of_species;
};

#endif
