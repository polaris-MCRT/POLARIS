#pragma once
#include "Faddeeva.hh"
#include "Grid.h"
#include "Matrix2D.h"
#include "chelper.h"
#include "typedefs.h"
#include <complex>

#ifndef CMOLECULE
#define CMOLECULE

class CGasSpecies
{
  public:
    CGasSpecies()
    {
        molecular_weight = 0;
        lvl_pop_type = POP_LTE;

        nr_of_energy_levels = 0;
        nr_of_transitions = 0;
        nr_of_col_partner = 0;
        nr_of_spectral_lines = 0;

        abundance = 0;

        gas_species_radius = 0;

        nr_zeeman_spectral_lines = 0;
        zeeman_spectral_lines = 0;
        zeeman_splitting = false;

        j_level = 0;
        g_level = 0;
        energy_level = 0;
        nr_of_col_transition = 0;

        spectral_lines = 0;

        upper_level = 0;
        lower_level = 0;

        trans_einstA = 0;
        trans_einstB_l = 0;
        trans_einstB_u = 0;

        trans_freq = 0;
        trans_inner_energy = 0;

        orientation_H2 = 0;
        nr_of_col_temp = 0;

        lande_upper = 0;
        lande_lower = 0;

        nr_sublevels_upper = 0;
        nr_sublevels_lower = 0;

        nr_pi_spectral_lines = 0;
        nr_sigma_spectral_lines = 0;

        collision_temp = 0;

        col_upper = 0;
        col_lower = 0;

        line_strength_pi = 0;
        line_strength_sigma_p = 0;
        line_strength_sigma_m = 0;

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
                if(nr_of_col_transition != 0)
                {
                    for(int j = 0; j < nr_of_col_transition[i]; j++)
                    {
                        delete[] col_matrix[i][j];
                        col_matrix[i][j] = 0;
                    }
                }
            }
            col_matrix = 0;
        }

        if(j_level != 0)
            delete[] j_level;

        if(g_level != 0)
            delete[] g_level;

        if(energy_level != 0)
            delete[] energy_level;

        if(nr_of_col_transition != 0)
            delete[] nr_of_col_transition;

        if(upper_level != 0)
            delete[] upper_level;

        if(lower_level != 0)
            delete[] lower_level;

        if(trans_einstA != 0)
            delete[] trans_einstA;

        if(trans_einstB_l != 0)
            delete[] trans_einstB_l;

        if(trans_einstB_u != 0)
            delete[] trans_einstB_u;

        if(orientation_H2 != 0)
            delete[] orientation_H2;

        if(lande_upper != 0)
            delete[] lande_upper;

        if(lande_lower != 0)
            delete[] lande_lower;

        if(nr_sublevels_upper != 0)
            delete[] nr_sublevels_upper;

        if(nr_sublevels_lower != 0)
            delete[] nr_sublevels_lower;

        if(nr_pi_spectral_lines != 0)
            delete[] nr_pi_spectral_lines;

        if(nr_sigma_spectral_lines != 0)
            delete[] nr_sigma_spectral_lines;

        if(nr_of_col_temp != 0)
            delete[] nr_of_col_temp;

        if(collision_temp != 0)
        {
            for(uint i = 0; i < nr_of_col_partner; i++)
            {
                if(collision_temp[i] != 0)
                {
                    delete[] collision_temp[i];
                    collision_temp[i] = 0;
                }
            }
            collision_temp = 0;
        }

        if(col_upper != 0)
        {
            for(uint i = 0; i < nr_of_col_partner; i++)
                if(col_upper[i] != 0)
                {
                    delete[] col_upper[i];
                    col_upper[i] = 0;
                }
            col_upper = 0;
        }

        if(col_lower != 0)
        {
            for(uint i = 0; i < nr_of_col_partner; i++)
                if(col_lower[i] != 0)
                {
                    delete[] col_lower[i];
                    col_lower[i] = 0;
                }
            col_lower = 0;
        }

        if(line_strength_pi != 0)
        {
            for(uint i = 0; i < nr_of_spectral_lines; i++)
                if(line_strength_pi[i] != 0)
                {
                    delete[] line_strength_pi[i];
                    line_strength_pi[i] = 0;
                }
            line_strength_pi = 0;
        }

        if(line_strength_sigma_p != 0)
        {
            for(uint i = 0; i < nr_of_spectral_lines; i++)
                if(line_strength_sigma_p[i] != 0)
                {
                    delete[] line_strength_sigma_p[i];
                    line_strength_sigma_p[i] = 0;
                }
            line_strength_sigma_p = 0;
        }

        if(line_strength_sigma_m != 0)
        {
            for(uint i = 0; i < nr_of_spectral_lines; i++)
                if(line_strength_sigma_m[i] != 0)
                {
                    delete[] line_strength_sigma_m[i];
                    line_strength_sigma_m[i] = 0;
                }
            line_strength_sigma_m = 0;
        }
    }

    string getGasSpeciesName()
    {
        return stringID;
    }

    double getLineStrengthSigmaP(uint i_line, uint i_sublevel)
    {
        return line_strength_sigma_p[i_line][i_sublevel];
    }

    double getLineStrengthSigmaM(uint i_line, uint i_sublevel)
    {
        return line_strength_sigma_m[i_line][i_sublevel];
    }

    double getLineStrengthPi(uint i_line, uint i_sublevel)
    {
        return line_strength_pi[i_line][i_sublevel];
    }

    double getAbundance()
    {
        return abundance;
    }

    int getTransitionFromSpectralLine(uint i_line)
    {
        return spectral_lines[i_line];
    }

    double getTransitionFrequency(uint i_trans)
    {
        return trans_freq[i_trans];
    }

    double getSpectralLineFrequency(uint i_line)
    {
        uint i_trans = getTransitionFromSpectralLine(i_line);
        return getTransitionFrequency(i_trans);
    }

    uint getNrOfSpectralLines()
    {
        return nr_of_spectral_lines;
    }

    double getEinsteinA(uint i_trans)
    {
        return trans_einstA[i_trans];
    }

    double getEinsteinBu(uint i_trans)
    {
        return trans_einstB_u[i_trans];
    }

    double getEinsteinBl(uint i_trans)
    {
        return trans_einstB_l[i_trans];
    }

    double getGLevel(uint i_energy)
    {
        return g_level[i_energy];
    }

    double getJLevel(uint i_energy)
    {
        return j_level[i_energy];
    }

    double getEnergylevel(uint i_energy)
    {
        return energy_level[i_energy];
    }

    int getUpperCollision(uint m, uint n)
    {
        return col_upper[m][n];
    }

    int getLowerCollision(uint m, uint n)
    {
        return col_lower[m][n];
    }

    double getCollisionTemp(uint m, uint n)
    {
        return collision_temp[m][n];
    }

    double getCollisionMatrix(uint m, uint n, uint k)
    {
        return col_matrix[m][n][k];
    }

    double getLandeUpper(uint i_line)
    {
        return lande_upper[i_line];
    }

    double getLandeLower(uint i_line)
    {
        return lande_lower[i_line];
    }

    double getCollisionRadius()
    {
        return gas_species_radius;
    }

    double getMolecularWeight()
    {
        return molecular_weight;
    }

    double getNumberDensity(CGridBasic * grid, photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getNumberDensity(grid, cell);
    }

    double getNumberDensity(CGridBasic * grid, cell_basic * cell)
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

    double getMassDensity(CGridBasic * grid, photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getMassDensity(grid, cell);
    }

    double getMassDensity(CGridBasic * grid, cell_basic * cell)
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

    uilist getUniqueTransitions()
    {
        return unique_spectral_lines;
    }

    uint getUniqueTransitions(uint i)
    {
        return unique_spectral_lines[i];
    }

    uint getUpperEnergyLevel(uint i_trans)
    {
        return upper_level[i_trans];
    }

    uint getLowerEnergyLevel(uint i_trans)
    {
        return lower_level[i_trans];
    }

    uint getNrOfEnergyLevels()
    {
        return nr_of_energy_levels;
    }

    uint getNrOfTransitions()
    {
        return nr_of_transitions;
    }

    uint getNrOfCollisionPartner()
    {
        return nr_of_col_partner;
    }

    uint getNrCollisionTransitions(uint i_col_partner)
    {
        return nr_of_col_transition[i_col_partner];
    }

    uint getNrCollisionTemps(uint i_col_partner)
    {
        return nr_of_col_temp[i_col_partner];
    }

    uint getLevelPopType()
    {
        return lvl_pop_type;
    }

    uint getOrientation_H2(uint i_col_partner)
    {
        return orientation_H2[i_col_partner];
    }

    int getNrSublevelsUpper(uint i_line)
    {
        return nr_sublevels_upper[i_line];
    }

    int getNrSublevelsLower(uint i_line)
    {
        return nr_sublevels_lower[i_line];
    }

    float getMaxMUpper(uint i_line)
    {
        return float((nr_sublevels_upper[i_line] - 1) / 2.0);
    }

    float getMaxMLower(uint i_line)
    {
        return float((nr_sublevels_lower[i_line] - 1) / 2.0);
    }

    bool getZeemanSplitting()
    {
        return zeeman_splitting;
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
        double gamma = 0;
        for(uint i = 0; i < i_trans; i++)
            gamma += trans_einstA[i];

        // "http://chem.libretexts.org/Core/Physical_and_Theoretical_Chemistry/
        //  ->
        //  Kinetics/Modeling_Reaction_Kinetics/Collision_Theory/Collisional_Cross_Section"
        // "http://www.phy.ohiou.edu/~mboett/astro401_fall12/broadening.pdf
        double v_th = sqrt(2.0 * con_kB * temp_gas / (molecular_weight * 1.0e-3 / con_Na));
        double col_param =
            dens_gas * PI * pow(con_r_bohr + gas_species_radius, 2) * sqrt(pow(v_th, 2) + pow(v_turb, 2));

        return gamma + 2 * col_param;
    }

    double getGaussA(double temp_gas, double v_turb)
    {
        double v_th = sqrt(2.0 * con_kB * temp_gas / (molecular_weight * 1.0e-3 / con_Na));
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

    Matrix2D getLineMatrix(CGridBasic * grid,
                           photon_package * pp,
                           uint i_line,
                           double velocity,
                           Vector3D mag_field,
                           double cos_theta,
                           double sin_theta,
                           double cos_2_phi,
                           double sin_2_phi)
    {
        if(getZeemanSplitting())
            return getZeemanSplittingMatrix(
                grid, pp, i_line, velocity, mag_field, cos_theta, sin_theta, cos_2_phi, sin_2_phi);
        else
            return getGaussLineMatrix(grid, pp, i_line, velocity);
    }

    Matrix2D getGaussLineMatrix(CGridBasic * grid, photon_package * pp, uint i_line, double velocity);
    double getGaussLineShape(CGridBasic * grid, photon_package * pp, double velocity);

    Matrix2D getZeemanSplittingMatrix(CGridBasic * grid,
                                      photon_package * pp,
                                      uint i_line,
                                      double velocity,
                                      Vector3D mag_field,
                                      double cos_theta,
                                      double sin_theta,
                                      double cos_2_phi,
                                      double sin_2_phi);

    StokesVector calcEmissivity(CGridBasic * grid, cell_basic * cell, uint i_line);
    StokesVector calcEmissivityForTransition(CGridBasic * grid, cell_basic * cell, uint i_trans);

    bool calcLTE(CGridBasic * grid, bool full = false);
    bool calcFEP(CGridBasic * grid, bool full = false);
    bool calcLVG(CGridBasic * grid, double kepler_star_mass, bool full = false);
    bool updateLevelPopulation(CGridBasic * grid, cell_basic * cell, double * J_total);

    double elem_LVG(CGridBasic * grid,
                    double dens_species,
                    double * tmp_lvl_pop,
                    double doppler,
                    double L,
                    double J_ext,
                    uint i_trans);

    void createMatrix(double * J_mid, Matrix2D & A, double * b, double *** final_col_para);

    double *** calcCollisionParameter(CGridBasic * grid, cell_basic * cell);
    double getColPartnerDensity(CGridBasic * grid, cell_basic * cell, uint i_col_partner);
    dlist calcCollisionRate(uint i_col_partner, uint i_col_transition, uint hi_i, double temp_gas);

    bool readGasParamaterFile(string filename, uint id, uint max);
    bool readZeemanParamaterFile(string filename);

  private:
    double ** collision_temp;
    int **col_upper, **col_lower;
    double ** line_strength_pi;
    double ** line_strength_sigma_p;
    double ** line_strength_sigma_m;

    double *** col_matrix;

    prob_list frequency_prob;

    double molecular_weight;
    double abundance;
    double max_velocity;
    double gas_species_radius;

    uint i_species;

    uint nr_of_transitions;
    uint nr_of_col_partner;
    uint nr_of_energy_levels;
    uint nr_of_spectral_lines;
    uint nr_zeeman_spectral_lines;
    uint lvl_pop_type;

    uilist unique_spectral_lines;

    int * nr_of_col_transition;
    int * nr_of_col_temp;
    int *nr_sublevels_upper, *nr_sublevels_lower;
    int *nr_pi_spectral_lines, *nr_sigma_spectral_lines;
    int * zeeman_spectral_lines;
    int * spectral_lines;

    double * energy_level;
    double * g_level;
    double * j_level;
    int *upper_level, *lower_level;
    double *trans_einstA, *trans_einstB_l, *trans_einstB_u;
    double *trans_freq, *trans_inner_energy;
    int * orientation_H2;
    double *lande_upper, *lande_lower;

    string stringID;
    string filename;
    string catalog_path;

    ostringstream velocity_information;

    bool zeeman_splitting;
};

class CGasMixture
{
  public:
    CGasMixture()
    {
        nr_of_species = 0;

        kepler_star_mass = 0;

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
                delete[] level_to_pos[i_species];
            delete[] level_to_pos;
        }

        if(line_to_pos != 0)
        {
            for(uint i_species = 0; i_species < nr_of_species; i_species++)
            {
                for(uint i_line = 0; i_line < getNrOfSpectralLines(i_species); i_line++)
                    delete[] line_to_pos[i_species][i_line];
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

    bool updateLevelPopulation(CGridBasic * grid, photon_package * pp, uint i_species, double * J_total);
    bool updateLevelPopulation(CGridBasic * grid, cell_basic * cell, uint i_species, double * J_total);

    bool getZeemanSplitting(uint i_species)
    {
        return single_species[i_species].getZeemanSplitting();
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

    int getNrSublevelsUpper(uint i_species, uint i_line)
    {
        return single_species[i_species].getNrSublevelsUpper(i_line);
    }

    int getNrSublevelsLower(uint i_species, uint i_line)
    {
        return single_species[i_species].getNrSublevelsLower(i_line);
    }

    float getMaxMUpper(uint i_species, uint i_line)
    {
        return single_species[i_species].getMaxMUpper(i_line);
    }

    float getMaxMLower(uint i_species, uint i_line)
    {
        return single_species[i_species].getMaxMLower(i_line);
    }

    double getAbundance(uint i_species)
    {
        return single_species[i_species].getAbundance();
    }

    double getLandeUpper(uint i_species, uint i_line)
    {
        return single_species[i_species].getLandeUpper(i_line);
    }

    double getLandeLower(uint i_species, uint i_line)
    {
        return single_species[i_species].getLandeLower(i_line);
    }

    double getCollisionRadius(uint i_species)
    {
        return single_species[i_species].getCollisionRadius();
    }

    double getLineStrengthSigmaP(uint i_species, uint i_line, uint i_sublevel)
    {
        return single_species[i_species].getLineStrengthSigmaP(i_line, i_sublevel);
    }

    double getLineStrengthSigmaM(uint i_species, uint i_line, uint i_sublevel)
    {
        return single_species[i_species].getLineStrengthSigmaM(i_line, i_sublevel);
    }

    double getLineStrengthPi(uint i_species, uint i_line, uint i_sublevel)
    {
        return single_species[i_species].getLineStrengthPi(i_line, i_sublevel);
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

    uint getTransitionFromSpectralLine(uint i_species, uint i_line)
    {
        return single_species[i_species].getTransitionFromSpectralLine(i_line);
    }

    double getGaussA(uint i_species, double temp_gas, double v_turb)
    {
        return single_species[i_species].getGaussA(temp_gas, v_turb);
    }

    double getKeplerStarMass()
    {
        return kepler_star_mass;
    }

    double getNumberDensity(CGridBasic * grid, photon_package * pp, uint i_species)
    {
        return single_species[i_species].getNumberDensity(grid, pp);
    }

    double getNumberDensity(CGridBasic * grid, cell_basic * cell, uint i_species)
    {
        return single_species[i_species].getNumberDensity(grid, cell);
    }

    double getMassDensity(CGridBasic * grid, photon_package * pp, uint i_species)
    {
        return single_species[i_species].getMassDensity(grid, pp);
    }

    double getMassDensity(CGridBasic * grid, cell_basic * cell, uint i_species)
    {
        return single_species[i_species].getMassDensity(grid, cell);
    }

    void setKeplerStarMass(double val)
    {
        kepler_star_mass = val;
    }

    void calcLineBroadening(CGridBasic * grid, uint i_species)
    {
        single_species[i_species].calcLineBroadening(grid);
    }

    Matrix2D getLineMatrix(CGridBasic * grid,
                           photon_package * pp,
                           uint i_species,
                           uint i_line,
                           double velocity,
                           Vector3D mag_field,
                           double cos_theta,
                           double sin_theta,
                           double cos_2_phi,
                           double sin_2_phi)
    {
        return single_species[i_species].getLineMatrix(
            grid, pp, i_line, velocity, mag_field, cos_theta, sin_theta, cos_2_phi, sin_2_phi);
    }

    Matrix2D getGaussLineMatrix(CGridBasic * grid,
                                photon_package * pp,
                                uint i_species,
                                uint i_line,
                                double velocity)
    {
        return single_species[i_species].getGaussLineMatrix(grid, pp, i_line, velocity);
    }

    double getGaussLineShape(CGridBasic * grid, photon_package * pp, uint i_species, double velocity)
    {
        return single_species[i_species].getGaussLineShape(grid, pp, velocity);
    }

    StokesVector calcEmissivity(CGridBasic * grid, cell_basic * cell, uint i_species, uint i_line)
    {
        return single_species[i_species].calcEmissivity(grid, cell, i_line);
    }

    StokesVector calcEmissivity(CGridBasic * grid, photon_package * pp, uint i_species, uint i_line)
    {
        return single_species[i_species].calcEmissivity(grid, pp->getPositionCell(), i_line);
    }

    StokesVector calcEmissivityForTransition(CGridBasic * grid,
                                             cell_basic * cell,
                                             uint i_species,
                                             uint i_trans)
    {
        return single_species[i_species].calcEmissivityForTransition(grid, cell, i_trans);
    }

    StokesVector calcEmissivityForTransition(CGridBasic * grid,
                                             photon_package * pp,
                                             uint i_species,
                                             uint i_trans)
    {
        return single_species[i_species].calcEmissivityForTransition(grid, pp->getPositionCell(), i_trans);
    }

    uint getTotalNrOfSpectralLines()
    {
        uint total_spectral_lines = 0;
        for(uint i_species = 0; i_species < nr_of_species; i_species++)
            total_spectral_lines += getNrOfSpectralLines(i_species);
        return total_spectral_lines;
    }

    uint getMaxNrOfSpectralLines()
    {
        uint max_spectral_lines = 0;
        for(uint i_species = 0; i_species < nr_of_species; i_species++)
            if(getNrOfSpectralLines(i_species) > max_spectral_lines)
                max_spectral_lines = getNrOfSpectralLines(i_species);
        return max_spectral_lines;
    }

    uint getNrOffsetEntries(CGridBasic * grid, parameters & param)
    {
        // Init variables and pointer arrays
        uint offset_entries = 0;

        // 1x Gauss_a + doppler_width, Gamma, voigt_a for each spectral line to simulate
        uint line_broadening_offset = 1 + 3 * getMaxNrOfSpectralLines();

        // Arrays to link energy levels, simulated spectral lines and position in the grid cells
        level_to_pos = new uint *[nr_of_species];
        line_to_pos = new uint **[nr_of_species];

        for(uint i_species = 0; i_species < nr_of_species; i_species++)
        {
            uint nr_of_energy_level = getNrOfEnergyLevels(i_species);
            uint nr_of_spectral_lines = getNrOfSpectralLines(i_species);

            level_to_pos[i_species] = new uint[nr_of_energy_level];
            line_to_pos[i_species] = new uint *[nr_of_spectral_lines];

            for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
                line_to_pos[i_species][i_line] = new uint[2];

            uint used_level_populations = 0;
            for(uint i_lvl = 0; i_lvl < nr_of_energy_level; i_lvl++)
            {
                // Initial set to MAX_UINT
                level_to_pos[i_species][i_lvl] = MAX_UINT;

                // Found a spectral line using the energy level
                bool found = false;

                // Add all energy levels to grid for MC level pop calculation
                if(param.gasSpeciesLevelPopTypeIsMC())
                {
                    level_to_pos[i_species][i_lvl] = line_broadening_offset + used_level_populations;
                    found = true;
                }

                for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
                {
                    uint i_trans = getTransitionFromSpectralLine(i_species, i_line);

                    if(i_lvl == getLowerEnergyLevel(i_species, i_trans))
                    {
                        level_to_pos[i_species][i_lvl] = line_broadening_offset + used_level_populations;
                        line_to_pos[i_species][i_line][0] = line_broadening_offset + used_level_populations;
                        found = true;
                    }
                    else if(i_lvl == getUpperEnergyLevel(i_species, i_trans))
                    {
                        level_to_pos[i_species][i_lvl] = line_broadening_offset + used_level_populations;
                        line_to_pos[i_species][i_line][1] = line_broadening_offset + used_level_populations;
                        found = true;
                    }
                }

                if(found)
                    used_level_populations++;
            }
            if(used_level_populations > offset_entries)
                offset_entries = used_level_populations;
        }

        return line_broadening_offset + offset_entries;
    }

    void printParameter(parameters & param, CGridBasic * grid);

  private:
    CGasSpecies * single_species;

    uint ** level_to_pos;
    uint *** line_to_pos;

    double kepler_star_mass;
    uint nr_of_species;
};

#endif
