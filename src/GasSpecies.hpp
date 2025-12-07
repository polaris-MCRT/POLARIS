/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CGAS_SPECIES_H
#define CGAS_SPECIES_H

#include "Faddeeva.hh"
#include "GridBasic.hpp"
#include "Matrix2D.hpp"
#include "Typedefs.hpp"
#include "MathFunctions.hpp"
#include "Parameters.hpp"
#include "Photon.hpp"
#include "Vector3D.hpp"
//#include "GrainCharge.hpp"

#define TRANS_SIGMA_P +1
#define TRANS_PI 0
#define TRANS_SIGMA_M -1

class CGasSpecies
{
public:

    //uint tmp_counter;

    CGasSpecies()
    {
        molecular_weight = 0;
        lvl_pop_type = POP_LTE;

        nr_of_energy_level = 0;
        nr_of_transitions = 0;
        nr_of_col_partner = 0;
        nr_of_spectral_lines = 0;

        //tmp_counter=0;

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

    string getGasSpeciesName() const;

    double getLineStrength(uint i_trans, uint i_sublvl_u, uint i_sublvl_l) const;

    double getAbundance() const;

    double getKeplerStarMass() const;

    int getTransitionFromSpectralLine(uint i_line) const;

    double getTransitionFrequency(uint i_trans) const;

    double getSpectralLineFrequency(uint i_line) const;

    uint getNrOfSpectralLines() const;

    double getEinsteinA(uint i_trans) const;

    double getEinsteinBul(uint i_trans) const;

    double getEinsteinBlu(uint i_trans) const;

    uint getSublevelIndex(uint i_trans, uint i_sublvl_u, uint i_sublvl_l) const;

    double getEinsteinA(uint i_trans, uint i_sublvl_u, uint i_sublvl_l) const;

    double getEinsteinBul(uint i_trans, uint i_sublvl_u, uint i_sublvl_l) const;

    double getEinsteinBlu(uint i_trans, uint i_sublvl_u, uint i_sublvl_l) const;

    double getJLevel(uint i_lvl) const;

    double getEnergyOfLevel(uint i_lvl) const;

    uint getUpperCollisionLevel(uint m, uint n) const;

    uint getLowerCollisionLevel(uint m, uint n) const;

    double getCollisionTemp(uint m, uint n) const;

    double getCollisionMatrix(uint m, uint n, uint k) const;

    double getLandeUpper(uint i_trans) const;

    double getLandeLower(uint i_trans) const;

    double getLande(uint i_lvl) const;

    double getCollisionRadius() const;

    double getMolecularWeight() const;

    double getNumberDensity(CGridBasic * grid, const photon_package & pp) const;

    double getNumberDensity(CGridBasic * grid, const cell_basic & cell) const;

    double getMassDensity(CGridBasic * grid, const photon_package & pp) const;

    double getMassDensity(CGridBasic * grid, const cell_basic & cell) const;

    Vector3D getCellVelocity(CGridBasic * grid, const cell_basic & cell, const Vector3D & tmp_pos) const;

    double getProjCellVelocityInterp(const Vector3D & tmp_pos,
                                     const Vector3D & dir_map_xyz,
                                     const VelFieldInterp & vel_field_interp);

    void initReferenceLists();

    uint getZeemanSplitIndex(uint i_trans);

    bool isTransZeemanSplit(uint i_trans) const;

    bool isLineZeemanSplit(uint i_line) const;

    bool isZeemanSplit() const;

    uint getUniqueLevelIndex(uint i_lvl, uint i_sublvl) const;

    uint getUniqueTransitionIndex(uint i_trans, uint i_sublvl_u = 0, uint i_sublvl_l = 0) const;

    uilist getUniqueTransitions() const;

    uint getUniqueTransitions(uint i) const;

    uint getUpperEnergyLevel(uint i_trans) const;

    uint getLowerEnergyLevel(uint i_trans) const;

    uint getNrOfEnergyLevels() const;

    uint getNrOfTotalEnergyLevels() const;

    uint getNrOfTotalTransitions() const;

    uint getNrOfTransitions() const;

    uint getNrOfCollisionPartner() const;

    uint getNrCollisionTransitions(uint i_col_partner) const;

    uint getNrCollisionTemps(uint i_col_partner) const;

    uint getLevelPopType() const;

    uint getOrientation_H2(uint i_col_partner) const;

    int getNrOfSublevelUpper(uint i_trans) const;

    int getNrOfSublevelLower(uint i_trans) const;

    int getNrOfSublevel(uint i_lvl) const;

    int getNrOfTransBetweenSublevels(uint i_trans) const;

    float getMaxMUpper(uint i_trans) const;

    float getMaxMLower(uint i_trans) const;

    float getMaxM(uint i_lvl) const;

    void setKeplerStarMass(double val);

    void setLevelPopType(uint type);

    void setMolecularWeight(double val);

    void setAbundance(double val);

    void setNrOfSpectralLines(uint val);

    void setSpectralLines(int * lines);

    double getGamma(uint i_trans, double dens_gas, double dens_species, double temp_gas, double v_turb);

    double getGaussA(double temp_gas, double v_turb);

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
    };

    void getGaussLineMatrix(CGridBasic * grid,
                            const cell_basic & cell,
                            double velocity,
                            Matrix2D * line_absorption_matrix) const;

    void getGaussLineMatrix(CGridBasic * grid,
                            const photon_package & pp,
                            double velocity,
                            Matrix2D * line_absorption_matrix) const;

    double getGaussLineShape(CGridBasic * grid, const cell_basic & cell, double velocity) const;

    double getGaussLineShape(CGridBasic * grid, const photon_package & pp, double velocity) const;

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

    void calcEmissivityFromLvlPop(uint i_trans,
                                  uint i_sublvl_u,
                                  uint i_sublvl_l,
                                  double dens_species,
                                  double gauss_a,
                                  double * tmp_lvl_pop,
                                  double * j,
                                  double * alpha);

    double calcJFromInteractionLength(double j, double alpha, double J_ext, double L);

    double calcJFromOpticalDepth(double j, double alpha, double J_ext, double tau);

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
    uint **col_upper, **col_lower;
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

#endif /* CGAS_SPECIES_H */
