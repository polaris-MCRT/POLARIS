#pragma once
#include "typedefs.h"
#include "chelper.h"
#include "Matrix2D.h"
#include "Grid.h"
#include "Faddeeva.hh"
#include <complex>

#ifndef CMOLECULE
#define CMOLECULE

class CGasSpecies
{
public:

    CGasSpecies()
    {
        molecular_weight = 0;
        lvl_pop_type = 1;

        nr_of_energy_levels = 0;
        nr_of_total_transitions = 0;
        nr_of_col_partner = 0;
        nr_of_transitions = 0;

        abundance = 0;

        gas_species_radius = 0;

        nr_zeeman_transitions = 0;
        zeeman_transitions = 0;
        zeeman_splitting = false;

        j_level = 0;
        g_level = 0;
        energy_level = 0;
        nr_of_collision_transition = 0;

        transitions = 0;

        trans_upper = 0;
        trans_lower = 0;

        trans_einstA = 0;
        trans_einstB_l = 0;
        trans_einstB_u = 0;

        trans_freq = 0;
        trans_inneregy = 0;

        orientation_H2 = 0;
        nr_of_col_temp = 0;

        lande_upper = 0;
        lande_lower = 0;

        nr_sublevels_upper = 0;
        nr_sublevels_lower = 0;

        nr_pi_transitions = 0;
        nr_sigma_transitions = 0;

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
                if(nr_of_collision_transition != 0)
                {
                    for(int j = 0; j < nr_of_collision_transition[i]; j++)
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

        if(nr_of_collision_transition != 0)
            delete[] nr_of_collision_transition;

        if(trans_upper != 0)
            delete[] trans_upper;

        if(trans_lower != 0)
            delete[] trans_lower;

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

        if(nr_pi_transitions != 0)
            delete[] nr_pi_transitions;

        if(nr_sigma_transitions != 0)
            delete[] nr_sigma_transitions;

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
            for(uint i = 0; i < nr_of_transitions; i++)
                if(line_strength_pi[i] != 0)
                {
                    delete[] line_strength_pi[i];
                    line_strength_pi[i] = 0;
                }
            line_strength_pi = 0;
        }

        if(line_strength_sigma_p != 0)
        {
            for(uint i = 0; i < nr_of_transitions; i++)
                if(line_strength_sigma_p[i] != 0)
                {
                    delete[] line_strength_sigma_p[i];
                    line_strength_sigma_p[i] = 0;
                }
            line_strength_sigma_p = 0;
        }

        if(line_strength_sigma_m != 0)
        {
            for(uint i = 0; i < nr_of_transitions; i++)
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

    int getTransition(uint i_line)
    {
        return transitions[i_line];
    }

    double getTransitionFrequency(uint i_trans)
    {
        return trans_freq[i_trans];
    }

    double getTransitionFrequencyFromIndex(uint i_line)
    {
        uint i_trans = getTransition(i_line);
        return getTransitionFrequency(i_trans);
    }

    uint getNrOfTransitions()
    {
        return nr_of_transitions;
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
        return unique_transitions;
    }

    uint getUniqueTransitions(uint i)
    {
        return unique_transitions[i];
    }

    uint getUpperTransition(uint i_trans)
    {
        return trans_upper[i_trans];
    }

    uint getLowerTransition(uint i_trans)
    {
        return trans_lower[i_trans];
    }

    uint getNrEnergyLevels()
    {
        return nr_of_energy_levels;
    }

    uint getNrOfTotalTransitions()
    {
        return nr_of_total_transitions;
    }

    uint getNrCollisionPartner()
    {
        return nr_of_col_partner;
    }

    uint getNrCollisionTransitions(uint i_col_partner)
    {
        return nr_of_collision_transition[i_col_partner];
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

    void setNrTransitions(uint val)
    {
        nr_of_transitions = val;
    }

    void setTransitions(int * trans)
    {
        transitions = trans;
    }

    double getGamma(uint i_trans, double dens_gas, double dens_species, double temp_gas, double v_turb)
    {
        double gamma = 0;
        for(uint i = 0; i < i_trans; i++)
            gamma += trans_einstA[i];

        // "http://chem.libretexts.org/Core/Physical_and_Theoretical_Chemistry/
        //  -> Kinetics/Modeling_Reaction_Kinetics/Collision_Theory/Collisional_Cross_Section"
        // "http://www.phy.ohiou.edu/~mboett/astro401_fall12/broadening.pdf
        double v_th = sqrt(2.0 * con_kB * temp_gas / (molecular_weight * 1.0e-3 / con_Na));
        double col_param = dens_gas * PI * pow(con_r_bohr + gas_species_radius, 2)
                * sqrt(pow(v_th, 2) + pow(v_turb, 2));

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

    Matrix2D getLineMatrix(CGridBasic * grid, photon_package * pp, uint i_line, double velocity,
        Vector3D mag_field, double cos_theta, double sin_theta, double cos_2_phi, double sin_2_phi)
    {
        if(getZeemanSplitting())
            return getZeemanSplittingMatrix(grid, pp, i_line, velocity,
                mag_field, cos_theta, sin_theta, cos_2_phi, sin_2_phi);
        else
            return getGaussLineMatrix(grid, pp, i_line, velocity);
    }

    Matrix2D  getGaussLineMatrix(CGridBasic * grid, photon_package * pp, uint i_line, double velocity);

    Matrix2D getZeemanSplittingMatrix(CGridBasic * grid, photon_package * pp, uint i_line, double velocity,
        Vector3D mag_field, double cos_theta, double sin_theta, double cos_2_phi, double sin_2_phi);

    StokesVector calcEmissivities(CGridBasic *grid, photon_package * pp, uint i_line);

    bool calcLTE(CGridBasic * grid);
    bool calcFEP(CGridBasic * grid);
    bool calcLVG(CGridBasic * grid, double kepler_star_mass);

    double elem_LVG(CGridBasic * grid, double dens_species, double * new_pop, double doppler,
            double L, double J_ext, uint i_trans);

    void createMatrix(double * J_mid, Matrix2D & A, double * b, Matrix2D final_col_para);
    Matrix2D calc_collision_parameter(CGridBasic * grid, double temp_gas, double gas_number_density);

    bool readGasParamaterFile(string filename, uint id, uint max);
    bool readZeemanParamaterFile(string filename);

private:

    double ** collision_temp;
    int ** col_upper, ** col_lower;
    double ** line_strength_pi;
    double ** line_strength_sigma_p;
    double ** line_strength_sigma_m;

    double *** col_matrix;

    double molecular_weight;
    double abundance;
    double max_velocity;
    double gas_species_radius;

    uint i_species;

    uint nr_of_total_transitions;
    uint nr_of_col_partner;
    uint nr_of_energy_levels;
    uint nr_of_transitions;
    uint nr_zeeman_transitions;
    uint lvl_pop_type;

    uilist unique_transitions;

    int * nr_of_collision_transition;
    int * nr_of_col_temp;
    int * nr_sublevels_upper, *nr_sublevels_lower;
    int * nr_pi_transitions, *nr_sigma_transitions;
    int * zeeman_transitions;
    int * transitions;

    double * energy_level;
    double * g_level;
    double * j_level;
    int * trans_upper, *trans_lower;
    double * trans_einstA, *trans_einstB_l, *trans_einstB_u;
    double * trans_freq, *trans_inneregy;
    int * orientation_H2;
    double * lande_upper, *lande_lower;

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

        single_species = 0;
    }

    ~CGasMixture(void)
    {
        if(single_species != 0)
            delete[] single_species;
    }

    string getGasSpeciesName(uint i_species)
    {
        return single_species[i_species].getGasSpeciesName();
    }

    bool createGasSpecies(parameters & param);

    bool calcLevelPopulation(CGridBasic * grid, uint i_species);

    bool getZeemanSplitting(uint i_species)
    {
        return single_species[i_species].getZeemanSplitting();
    }

    uilist getUniqueTransitions(uint i_species)
    {
        return single_species[i_species].getUniqueTransitions();
    }

    uint getUniqueTransitions(uint i_species, uint i)
    {
        return single_species[i_species].getUniqueTransitions(i);
    }

    uint getLevelPopType(uint i_species)
    {
        return single_species[i_species].getLevelPopType();
    }

    uint getNrOfTransitions(uint i_species)
    {
        return single_species[i_species].getNrOfTransitions();
    }

    uint getUpperTransition(uint i_species, uint i_trans)
    {
        return single_species[i_species].getUpperTransition(i_trans);
    }

    uint getLowerTransition(uint i_species, uint i_trans)
    {
        return single_species[i_species].getLowerTransition(i_trans);
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

    double getTransitionFrequencyFromIndex(uint i_species, uint i_line)
    {
        uint i_trans = single_species[i_species].getTransition(i_line);
        return single_species[i_species].getTransitionFrequency(i_trans);
    }

    uint getNrOfTotalTransitions(uint i_species)
    {
        return single_species[i_species].getNrOfTotalTransitions();
    }

    uint getTransition(uint i_species, uint i_line)
    {
        return single_species[i_species].getTransition(i_line);
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

    Matrix2D getLineMatrix(CGridBasic * grid, photon_package * pp, uint i_species, uint i_line, double velocity,
        Vector3D mag_field, double cos_theta, double sin_theta, double cos_2_phi, double sin_2_phi)
    {
        return single_species[i_species].getLineMatrix(grid, pp, i_line, velocity,
            mag_field, cos_theta, sin_theta, cos_2_phi, sin_2_phi);
    }

    StokesVector calcEmissivities(CGridBasic * grid, photon_package * pp, uint i_species, uint i_line)
    {
        return single_species[i_species].calcEmissivities(grid, pp, i_line);
    }

    void printParameter(parameters & param, CGridBasic * grid);

private:
    CGasSpecies * single_species;
    double kepler_star_mass;
    uint nr_of_species;
};

#endif
