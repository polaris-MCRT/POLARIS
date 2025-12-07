/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CDUST_MIXTURE_H
#define CDUST_MIXTURE_H

#include "GridBasic.hpp"
#include "Matrix2D.hpp"
#include "Typedefs.hpp"
#include "CellBasic.hpp"
#include "MathFunctions.hpp"
#include "DustComponent.hpp"
#include "Photon.hpp"
#include "Stokes.hpp"
#include "Vector3D.hpp"

class CDustMixture
{
public:
    CDustMixture(void)
    {
        single_component = 0;
        mixed_component = 0;

        nr_of_dust_species = 0;

        scattering_to_raytracing = false;

        extinction_magnitude = 0;
        extinction_magnitude_wavelength = 0;
        extinction_dust_mixture = MAX_UINT;

        nr_of_components = 0;
        nr_of_wavelength = 0;
        wavelength_offset = 0;
        
        mark_counter = 0;
    }

    ~CDustMixture(void)
    {
        if(single_component != 0)
            delete[] single_component;
            
        if(mixed_component != 0)
            delete[] mixed_component;
    }

    uint getMixtureID(CGridBasic * grid, const cell_basic & cell) const;

    uint getMixtureID(CGridBasic * grid, const photon_package & pp) const;

    double getAvgMass(uint i_mixture);

    bool writeComponentData(string path_data);

    bool writeComponentPlot(string path_plot);

    void preCalcEffProperties(parameters & param);

    void preCalcRelWeight();

    string getPhaseFunctionStr(uint i_mixture);

    void killSingleComponents();

    uint getNrOfDustSpecies(uint i_mixture);

    uint getNrOfStochasticSizes(uint i_mixture);
    
    uint getNrOfNanoSizes(uint i_mixture);

    void calcTemperature(CGridBasic * grid, cell_basic * cell, bool use_energy_density);

    void calcStochasticHeatingPropabilities(CGridBasic * grid,
                                            cell_basic * cell,
                                            dlist & wl_list) const;

    double getForegroundExtinction(double wavelength);

    double getCextMean(CGridBasic * grid, const photon_package & pp) const;

    double getCabsMean(CGridBasic * grid, const photon_package & pp) const;

    double getCscaMean(CGridBasic * grid, const photon_package & pp) const;
    
    double getMaxSubTemperature() const;

    bool adjustTempAndWavelengthBW(CGridBasic * grid, photon_package * pp, bool use_energy_density, CRandomGenerator * rand_gen);

    void calcAlignedRadii(CGridBasic * grid, cell_basic * cell);
    
    void calcAME(CGridBasic * grid, cell_basic * cell);

    void addToWavelengthGrid(double wavelength);

    void addToWavelengthGrid(double lam_min, double lam_max, double nr_of_wavelength, bool add_offset = false);

    void finalizeWavelengthList();

    void setScatteringToRay(bool val);

    bool getScatteringToRay();

    const dlist & getWavelengthList() const;

    double getWavelength(uint wID);

    double getWavelength(photon_package * pp);
    
    void calc_dust_emi_ame(CGridBasic * grid, const cell_basic * cell, double lambda, double & j_ame, double & nd);

    uint getNrOfWavelength();

    uint getWavelengthID(double wavelength);

    double getSizeMin(uint i_mixture);

    double getSizeMin(CGridBasic * grid, const cell_basic & cell) const;

    double getEffectiveRadius(uint i_mixture, uint a);

    double * getEffectiveRadii(uint i_mixture);

    uint getNrOfCalorimetryTemperatures(uint i_mixture);

    // double getMinDustTemp();

    double getMaxDustTemp();

    double getMinAlignedRadius();

    double getMaxAlignedRadius();

    // Dust number density functions
    double getNumberDensity(CGridBasic * grid, const cell_basic & cell) const;

    double getNumberDensity(CGridBasic * grid, const photon_package & pp) const;

    double getNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_mixture) const;

    double getNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_mixture) const;

    // Dust mass density functions
    double getMassDensity(CGridBasic * grid, const cell_basic & cell) const;

    double getMassDensity(CGridBasic * grid, const photon_package & pp) const;

    double getMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_mixture) const;

    double getMassDensity(CGridBasic * grid, const photon_package & pp, uint i_mixture) const;

    double getRelativeDustNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const;

    double getRelativeDustNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_density) const;

    double getRelativeDustMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const;

    double getRelativeDustMassDensity(CGridBasic * grid, const photon_package & pp, uint i_density) const;

    void calcEmissivityHz(CGridBasic * grid, const photon_package & pp, StokesVector * dust_emissivity);

    double calcEmissivity(CGridBasic * grid, const photon_package & pp);

    void calcEmissivityExt(CGridBasic * grid, const photon_package & pp, Matrix2D * dust_ext_matrix) const;

    void calcEmissivityEmi(CGridBasic * grid,
                           const photon_package & pp,
                           uint i_offset,
                           uint emission_component,
                           StokesVector * dust_emissivity) const;

    StokesVector getRadFieldScatteredFraction(CGridBasic * grid,
                                              const photon_package & pp,
                                              const Vector3D & en_dir,
                  
                                              double energy) const;
    
    void markCells(CGridBasic * grid, parameters & param);
    
    uint getNrOfMixtures() const;

    void convertTempInQB(CGridBasic * grid, cell_basic * cell, double min_gas_density, bool use_gas_temp);

    uint getScatteringMixture(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen) const;

    uint getEmittingMixture(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen) const;

    void scatter(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen, bool adjust_stokes = false);

    void setGridRequirements(CGridBasic * grid, parameters & param);

    void setIDs(CDustComponent & component,
                uint mix_id, uint i_comp,
                uint nr_of_components,
                uint i_mixture,
                uint nr_of_mixtures);

    double *** getSizeFractions();

    void getEscapePhoton(CGridBasic * grid,
                         photon_package * pp,
                         Vector3D obs_ex,
                         Vector3D dir_obs,
                         photon_package * pp_escape,
                         CRandomGenerator * rand_gen) const;

    double getCellEmission(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen) const;

    double getTotalCellEmission(CGridBasic * grid, const photon_package & pp) const;

    double getPlanck(uint w, double temp);

    bool createDustMixtures(parameters & param, string path_data, string path_plot);
    bool readScatteringMatrices(parameters & param);
    bool readColarimetry(parameters & param);
    bool mixComponents(parameters & param, uint i_mixture);
    bool preCalcDustProperties(parameters & param, uint i_mixture);

    void printParameters(parameters & param, CGridBasic * grid);

    void getNrOfUniqueScatTheta(uint ** & nr_of_scat_theta, double *** & scat_theta);

private:
    CDustComponent * single_component;
    CDustComponent * mixed_component;
    
    uint mark_counter;

    bool scattering_to_raytracing;

    double extinction_magnitude;
    double extinction_magnitude_wavelength;
    uint extinction_dust_mixture;

    uint nr_of_components;
    uint nr_of_wavelength;
    uint wavelength_offset;
    uint nr_of_dust_species;

    uilist unique_IDs;
    //uilist dust_choices_to_index;

    dlist wavelength_list;
    
    spline diff_y; // diff_y as a function of z
};

#endif /* CDUST_MIXTURE_H */
