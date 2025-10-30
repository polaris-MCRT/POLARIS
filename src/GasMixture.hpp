/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CGAS_MIXTURE_H
#define CGAS_MIXTURE_H

#include "GasSpecies.hpp"

#define TRANS_SIGMA_P +1
#define TRANS_PI 0
#define TRANS_SIGMA_M -1

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

    string getGasSpeciesName(uint i_species);

    bool createGasSpecies(parameters & param);

    bool calcLevelPopulation(CGridBasic * grid, uint i_species);

    bool updateLevelPopulation(CGridBasic * grid, cell_basic * cell, uint i_species, double * J_total);
    bool updateLevelPopulation(CGridBasic * grid, photon_package * pp, uint i_species, double * J_total);

    bool updateZeemanLevelPopulation(CGridBasic * grid,
                                     cell_basic * cell,
                                     uint i_species,
                                     uint i_line,
                                     double * sublvl_fraction);
    bool updateZeemanLevelPopulation(CGridBasic * grid,
                                     photon_package * pp,
                                     uint i_species,
                                     uint i_line,
                                     double * sublvl_fraction);

    bool isLineZeemanSplit(uint i_species, uint i_line);

    uint getZeemanSplitIndex(uint i_species, uint i_trans);

    bool isTransZeemanSplit(uint i_species, uint i_trans);

    bool isZeemanSplit(uint i_species);
    
    bool hasZeemanLines();

    uilist getUniqueTransitions(uint i_species);

    uint getNrOfSpecies();

    uint getUniqueTransitions(uint i_species, uint i);

    uint getLevelPopType(uint i_species);

    uint getNrOfSpectralLines(uint i_species);

    uint getNrOfEnergyLevels(uint i_species);

    uint getUpperEnergyLevel(uint i_species, uint i_trans);

    uint getLowerEnergyLevel(uint i_species, uint i_trans);

    int getNrOfSublevelUpper(uint i_species, uint i_trans);

    int getNrOfSublevelLower(uint i_species, uint i_trans);

    float getMaxMUpper(uint i_species, uint i_trans);

    float getMaxMLower(uint i_species, uint i_trans);

    float getMaxM(uint i_species, uint i_lvl);

    double getAbundance(uint i_species);

    double getLandeUpper(uint i_species, uint i_trans);

    double getLandeLower(uint i_species, uint i_trans);

    double getLande(uint i_species, uint i_lvl);

    double getCollisionRadius(uint i_species);

    double getLineStrength(uint i_species, uint i_trans, uint i_sublvl_u, uint i_sublvl_l);

    double getMolecularWeight(uint i_species);

    double getTransitionFrequency(uint i_species, uint i_trans);

    double getSpectralLineFrequency(uint i_species, uint i_line);

    uint getNrOfTransitions(uint i_species);

    uint getNrOfTotalTransitions(uint i_species);

    uint getTransitionFromSpectralLine(uint i_species, uint i_line);

    double getGaussA(uint i_species, double temp_gas, double v_turb);

    double getKeplerStarMass() const;

    double getNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_species) const;

    double getNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_species) const;

    double getMassDensity(CGridBasic * grid, const photon_package & pp, uint i_species) const;

    double getMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_species) const;

    void setKeplerStarMass(double val);

    void calcLineBroadening(CGridBasic * grid, uint i_species);

    void applyRadiationFieldFactor(uint i_species,
                                   uint i_trans,
                                   double sin_theta,
                                   double cos_theta,
                                   double energy,
                                   double * J_nu) const;

    void calcEmissivity(CGridBasic * grid,
                        const photon_package & pp,
                        uint i_species,
                        uint i_trans,
                        double velocity,
                        const LineBroadening & line_broadening,
                        const MagFieldInfo & mag_field_info,
                        StokesVector * line_emissivity,
                        Matrix2D * line_absorption_matrix) const;

    uint getNrOfSublevel(uint i_species, uint i_lvl);

    uint getTotalNrOfSpectralLines();

    uint getMaxNrOfSublevel(uint i_species);

    uint getUniqueLevelIndex(uint i_species, uint i_lvl, uint i_sublvl);

    uint getNrOffsetEntries(CGridBasic * grid, parameters & param);

    double getProjCellVelocity(CGridBasic * grid, const photon_package & pp, const Vector3D & tmp_pos) const;

    Vector3D getCellVelocity(CGridBasic * grid, const photon_package & pp, const Vector3D & tmp_pos) const;

    Vector3D getCellVelocity(CGridBasic * grid, const cell_basic & cell, const Vector3D & tmp_pos) const;

    double getProjCellVelocityInterp(const Vector3D & tmp_pos,
                                     const Vector3D & dir_map_xyz,
                                     const VelFieldInterp & vel_field_interp);

    void printParameters(parameters & param, CGridBasic * grid);

private:
    CGasSpecies * single_species;

    uint *** level_to_pos;
    uint **** line_to_pos;

    uint nr_of_species;

};

#endif /* CGAS_MIXTURE_H */
