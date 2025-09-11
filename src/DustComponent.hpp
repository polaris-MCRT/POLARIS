/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CDUST_COMPONENT_H
#define CDUST_COMPONENT_H

#include "GridBasic.hpp"
#include "Matrix2D.hpp"
#include "CrossSections.hpp"
#include "Typedefs.hpp"
#include "CellBasic.hpp"
#include "MathFunctions.hpp"
#include "MathInterp.hpp"
#include "Photon.hpp"
#include "Stokes.hpp"
#include "Vector3D.hpp"

class CDustComponent
{
public:
    CDustComponent()
    {
        avg_scattering_frac = 0;
        phase_pdf = 0;

        HG_g_factor = 0;
        HG_g2_factor = 0;
        HG_g3_factor = 0;
        Qtrq = 0;
        tab_em = 0;
        tab_em_inv = 0;
        tab_planck = 0;
        enthalpy = 0;

        avg_planck_frac = 0;

        dust_prob = 0;
        sca_prob = 0;
        abs_prob = 0;

        sca_mat = 0;
        Qext1 = 0;
        Qext2 = 0;
        Qabs1 = 0;
        Qabs2 = 0;
        Qsca1 = 0;
        Qsca2 = 0;
        Qcirc = 0;
        HGg = 0;
        HGg2 = 0;
        HGg3 = 0;
        a_eff = 0;
        grain_distribution_x_aeff_sq = 0;
        grain_size_distribution = 0;
        a_eff_squared = 0;
        mass = 0;
        relWeightTab = 0;
        fraction = 0;
        calorimetry_temperatures = 0;

        CextMean = 0;
        CabsMean = 0;
        CscaMean = 0;

        tCext1 = 0;
        tCext2 = 0;
        tCabs1 = 0;
        tCabs2 = 0;
        tCsca1 = 0;
        tCsca2 = 0;
        tCcirc = 0;
        tHGg = 0;
        tHGg2 = 0;
        tHGg3 = 0;

        stringID = "";
        size_keyword = "";

        stochastic_heating_max_size = 0;
        delta0 = 8.28e23 * 2.5e-12 * 1e8 * 1e-6 * 1e6;
        larm_f = 4.1e-19;
        aspect_ratio = 0;
        sub_temp = 1e6;
        material_density = 0;
        gold_g_factor = 0;
        dust_mass_fraction = 0;
        R_rayleigh = 1.0;

        Q_ref = 0.4;
        alpha_Q = 3.0;

        // min_temp = 0;
        max_temp = 0;
        min_a_alig = 1e200;
        max_a_alig = 0;
        f_highJ = 0.25;
        f_cor = 0.6;
        delta_rat = 2;
        mu = 0;
        avg_mass = 0;

        dust_offset = false;
        scat_loaded = false;
        calorimetry_loaded = false;
        sublimate = false;
        is_align = false;
        is_mixture = false;
        individual_dust_fractions = false;

        // Connection between mat_elem_counter and position in scattering matrix
        elements[0] = int(1);   // S11
        elements[1] = int(2);   // S12
        elements[2] = int(0);   // S13
        elements[3] = int(0);   // S14
        elements[4] = int(2);   // S21
        elements[5] = int(1);   // S22
        elements[6] = int(0);   // S23
        elements[7] = int(0);   // S24
        elements[8] = int(0);   // S31
        elements[9] = int(0);   // S32
        elements[10] = int(3);  // S33
        elements[11] = int(4);  // S34
        elements[12] = int(0);  // S41
        elements[13] = int(0);  // S42
        elements[14] = int(-4); // S43
        elements[15] = int(3);  // S44

        i_component = 0;
        nr_of_components = 0;
        nr_of_wavelength = 0;
        wavelength_offset = 0;
        i_mixture = 0;
        nr_of_mixtures = 0;
        calorimetry_type = 0;
        alignment = ALIG_PA;
        phID = 0;
        nr_of_dust_species = 0;
        nr_of_incident_angles = 0;
        nr_of_scat_theta = 0;
        scat_theta = 0;
        nr_of_scat_phi = 0;
        nr_of_scat_mat_elements = 0;
        nr_of_calorimetry_temperatures = 0;

        a_min_global = 1e200;
        a_max_global = 0;
        
        component_id = 0;
    }

    ~CDustComponent()
    {
        if(Qext1 != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] Qext1[a];
            delete[] Qext1;
        }
        if(Qext2 != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] Qext2[a];
            delete[] Qext2;
        }
        if(Qabs1 != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] Qabs1[a];
            delete[] Qabs1;
        }
        if(Qabs2 != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] Qabs2[a];
            delete[] Qabs2;
        }
        if(Qsca1 != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] Qsca1[a];
            delete[] Qsca1;
        }
        if(Qsca2 != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] Qsca2[a];
            delete[] Qsca2;
        }
        if(Qcirc != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] Qcirc[a];
            delete[] Qcirc;
        }
        if(HGg != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] HGg[a];
            delete[] HGg;
        }
        if(HGg2 != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] HGg2[a];
            delete[] HGg2;
        }
        if(HGg3 != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] HGg3[a];
            delete[] HGg3;
        }
        if(avg_scattering_frac != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] avg_scattering_frac[a];
            delete[] avg_scattering_frac;
        }
        if(phase_pdf != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] phase_pdf[a];
            delete[] phase_pdf;
        }
        if(enthalpy != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] enthalpy[a];
            delete[] enthalpy;
        }
        if(sca_mat != 0)
            cleanScatteringData();

        if(nr_of_scat_theta != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] nr_of_scat_theta[a];
            delete[] nr_of_scat_theta;
        }
        if(scat_theta != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
            {
                for(uint w = 0; w < nr_of_wavelength; w++)
                    delete[] scat_theta[a][w];
                delete[] scat_theta[a];
            }
            delete[] scat_theta;
        }
        if(Qtrq != 0)
            delete[] Qtrq;
        if(HG_g_factor != 0)
            delete[] HG_g_factor;
        if(HG_g2_factor != 0)
            delete[] HG_g2_factor;
        if(HG_g3_factor != 0)
            delete[] HG_g3_factor;

        if(CextMean != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] CextMean[a];
            delete[] CextMean;
        }
        if(CabsMean != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] CabsMean[a];
            delete[] CabsMean;
        }
        if(CscaMean != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] CscaMean[a];
            delete[] CscaMean;
        }

        if(tCext1 != 0)
            delete[] tCext1;
        if(tCext2 != 0)
            delete[] tCext2;
        if(tCabs1 != 0)
            delete[] tCabs1;
        if(tCabs2 != 0)
            delete[] tCabs2;
        if(tCsca1 != 0)
            delete[] tCsca1;
        if(tCsca2 != 0)
            delete[] tCsca2;
        if(tCcirc != 0)
            delete[] tCcirc;
        if(tHGg != 0)
            delete[] tHGg;
        if(tHGg2 != 0)
            delete[] tHGg2;
        if(tHGg3 != 0)
            delete[] tHGg3;

        if(tab_em != 0)
            delete[] tab_em;
        if(tab_em_inv != 0)
            delete[] tab_em_inv;
        if(tab_planck != 0)
            delete[] tab_planck;
        if(calorimetry_temperatures != 0)
            delete[] calorimetry_temperatures;
        if(avg_planck_frac != 0)
            delete[] avg_planck_frac;

        if(dust_prob != 0)
            delete[] dust_prob;
        if(sca_prob != 0)
            delete[] sca_prob;
        if(abs_prob != 0)
            delete[] abs_prob;

        if(a_eff != 0)
            delete[] a_eff;
        if(grain_distribution_x_aeff_sq != 0)
            delete[] grain_distribution_x_aeff_sq;
        if(grain_size_distribution != 0)
            delete[] grain_size_distribution;
        if(a_eff_squared != 0)
            delete[] a_eff_squared;
        if(mass != 0)
            delete[] mass;
        if(relWeightTab != 0)
            delete[] relWeightTab;
    }

    // ----------------------------------------------------------------------
    // ----------- Efficiencies for grain size and wavelength ---------------
    // ----------------------------------------------------------------------
    inline double getQext1(uint a, uint w) const;

    inline double getQext2(uint a, uint w) const;

    inline double getQabs1(uint a, uint w) const;

    inline double getQabs2(uint a, uint w) const;

    inline double getQsca1(uint a, uint w) const;

    inline double getQsca2(uint a, uint w) const;

    inline double getQcirc(uint a, uint w) const;

    inline double getHGg(uint a, uint w) const;

    inline double getHGg2(uint a, uint w) const;

    inline double getHGg3(uint a, uint w) const;

    void setHGg(uint a, uint w, double val);

    void setHGg2(uint a, uint w, double val);

    void setHGg3(uint a, uint w, double val);

    // ------------------------------------------------------------------------------------
    // ----------- Add values to efficiencies for grain size and wavelength
    // ---------------
    // ------------------------------------------------------------------------------------
    void addQext1(uint a, uint w, double val);

    void addQext2(uint a, uint w, double val);

    void addQabs1(uint a, uint w, double val);

    void addQabs2(uint a, uint w, double val);

    void addQsca1(uint a, uint w, double val);

    void addQsca2(uint a, uint w, double val);

    void addQcirc(uint a, uint w, double val);

    void addHGg(uint a, uint w, double val);

    void addHGg2(uint a, uint w, double val);

    void addHGg3(uint a, uint w, double val);

    // ------------------------------------------------------------------------
    // ----------- Cross-sections for grain size and wavelength ---------------
    // ------------------------------------------------------------------------
    double getCext1(uint a, uint w) const;

    double getCext2(uint a, uint w) const;

    double getCabs1(uint a, uint w) const;

    double getCabs2(uint a, uint w) const;

    double getCsca1(uint a, uint w) const;

    double getCsca2(uint a, uint w) const;

    double getCcirc(uint a, uint w) const;

    // -------------------------------------------------------------------------------
    // ----------- Cross-sections for wavelength mixed in current cell ---------------
    // -------------------------------------------------------------------------------
    double getCext1(CGridBasic * grid, const photon_package & pp) const;

    double getCext2(CGridBasic * grid, const photon_package & pp) const;

    double getCabs1(CGridBasic * grid, const photon_package & pp) const;

    double getCabs2(CGridBasic * grid, const photon_package & pp) const;

    double getCsca1(CGridBasic * grid, const photon_package & pp) const;

    double getCsca2(CGridBasic * grid, const photon_package & pp) const;

    double getCcirc(CGridBasic * grid, const photon_package & pp) const;

    double getHGg(CGridBasic * grid, const photon_package & pp) const;

    double getHGg2(CGridBasic * grid, const photon_package & pp) const;

    double getHGg3(CGridBasic * grid, const photon_package & pp) const;

    // -----------------------------------------------------------------------------
    // ----------- Average cross-sections for grain size and wavelength ------------
    // -----------------------------------------------------------------------------

    inline double getCextMean(uint a, uint w) const;

    inline double getCabsMean(uint a, uint w) const;

    inline double getCscaMean(uint a, uint w) const;

    // -----------------------------------------------------------------------------
    // ----------- Average cross-sections for wavelength mixed in current cell -----
    // -----------------------------------------------------------------------------
    double getCextMean(CGridBasic * grid, const photon_package & pp) const;

    double getCabsMean(CGridBasic * grid, const photon_package & pp) const;

    double getCscaMean(CGridBasic * grid, const photon_package & pp) const;

    // ------------------------------------------------------------------------------------
    // ----------- Average cross-sections for wavelength mixed with global conditions -----
    // ------------------------------------------------------------------------------------
    double getCextMean(double w) const;

    double getCabsMean(double w) const;

    double getCscaMean(double w) const;

    // ------------------------------------------------------------------------------
    // ----------- Cross-sections for wavelength mixed with global limits -----------
    // ------------------------------------------------------------------------------
    double getCext1(uint w) const;

    double getCext2(uint w) const;

    double getCabs1(uint w) const;

    double getCabs2(uint w) const;

    double getCsca1(uint w) const;

    double getCsca2(uint w) const;

    double getCcirc(uint w) const;

    double getHGg(uint w) const;

    double getHGg2(uint w) const;

    double getHGg3(uint w) const;

    // ------------------------------------------------------------------------------
    // ----------- Efficiencies for wavelength mixed with global limits -----------
    // ------------------------------------------------------------------------------
    double getQext1(uint w) const;

    double getQext2(uint w) const;

    double getQabs1(uint w) const;

    double getQabs2(uint w) const;

    double getQsca1(uint w) const;

    double getQsca2(uint w) const;

    double getQcirc(uint w) const;

    // -----------------------------------------------------------------------------------
    // ----------- Mass cross-sections for wavelength mixed with global limits -----------
    // -----------------------------------------------------------------------------------
    double getKappaExt1(uint w) const;

    double getKappaExt2(uint w) const;

    double getKappaAbs1(uint w) const;

    double getKappaAbs2(uint w) const;

    double getKappaSca1(uint w) const;

    double getKappaSca2(uint w) const;

    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------

    double getDeltaRat() const;

    double * getRelWeight(double a_min, double a_max, double size_param = 0) const;

    double getRelWeightTab(uint a) const;

    void preCalcRelWeight();

    double getWeight() const;

    double getMassWeight();

    double getGoldFactor();

    double getLarmF();

    double getQref();

    double getAlphaQ();

    double getRayleighReductionFactor();

    void setAlignmentMechanism(uint al);

    void setPhaseFunctionID(uint ph);

    uint getPhaseFunctionID();

    uint getNrOfStochasticSizes();

    // double getMinDustTemp();

    double getMaxDustTemp();

    double getMinAlignedRadius();

    double getMaxAlignedRadius();

    double getScatteringMatrixElement(uint a,
                                      uint w,
                                      uint incID,
                                      uint sphID,
                                      uint sthID,
                                      uint i_mat,
                                      uint j_mat) const;

    const Matrix2D & getScatteringMatrix(uint a, uint w, uint incID, uint sphID, uint sthID) const;

    void cleanScatteringData();

    StokesVector getRadFieldScatteredFraction(CGridBasic * grid,
                                              const photon_package & pp,
                                              uint i_density,
                                              const Vector3D & en_dir,
                                              double energy) const;

    double getScatteredFraction(uint a, uint w, double theta) const;

    double getScatteredFractionMie(uint a, uint w, double theta) const;

    double getScatteredFractionMie(uint a, uint w, uint sth) const;

    void scatter(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen);

    double getStochasticHeatingMaxSize() const;

    void setStochasticHeatingMaxSize(double val);

    double getFraction() const;

    void setFraction(double val);

    double getDustMassFraction() const;

    void setDustMassFraction(double val);

    void setSizeParameter(string size_key, dlist size_parameter_list);

    string getDustSizeKeyword();

    string getDustSizeParameterString();

    double convDensityToNumber(CGridBasic * grid, const cell_basic & cell, bool from_gas = false) const;

    double convDensityToMass(CGridBasic * grid, const cell_basic & cell, bool from_gas = false) const;

    double getNumberDensity(CGridBasic * grid, const cell_basic & cell) const;

    double getNumberDensity(CGridBasic * grid, const photon_package & pp) const;

    double getNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const;

    double getNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_density) const;

    double getMassDensity(CGridBasic * grid, const cell_basic & cell) const;

    uint getMassDensity(CGridBasic * grid, const photon_package & pp) const;

    double getMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const;

    void setMu(double mu_);

    void setMaterialDensity(double dens);

    bool checkGrainSizeLimits(double a_min, double a_max);

    double combinedRFactor(double R1, double R2, double R3) const;

    void setDelta0(double val);

    void setLarmF(double val);

    double getAspectRatio();

    double getSizeParam(CGridBasic * grid, const cell_basic & cell) const;

    double getSizeParam(CGridBasic * grid, const photon_package & pp) const;

    double getSizeParam() const;

    double getSizeMin(CGridBasic * grid, const cell_basic & cell) const;

    double getSizeMin(CGridBasic * grid, const photon_package & pp) const;

    double getSizeMin() const;

    void setSizeMin(double val);

    double getSizeMax(CGridBasic * grid, const cell_basic & cell) const;

    double getSizeMax(CGridBasic * grid, const photon_package & pp) const;

    double getSizeMax() const;

    void setSizeMax(double val);

    bool sizeIndexUsed(uint a, double a_min, double a_max) const;

    bool sizeIndexUsed(uint a) const;

    double getInternalIDG(double Td, double Tg) const;

    double getInternalGOLD(double Td, double Tg, double vg) const;

    double getInternalRAT() const;

    double findTemperature(uint a, double qb) const;

    double findTemperature(CGridBasic * grid, cell_basic * cell, double qb) const;

    uint findTemperatureID(double t) const;

    double getPlanck(uint w, double temp) const;

    double getAbsRate(CGridBasic * grid, const cell_basic & cell, uint a, bool use_energy_density) const;

    double getAbsRate(CGridBasic * grid, const photon_package & pp, uint a, bool use_energy_density) const;

    void setFHighJ(double val);

    void setFcorr(double val);

    void setQref(double val);

    void setAlphaQ(double val);

    void setRayleighReductionFactor(double val);

    double getQrat(uint a, uint w, double theta) const;

    double getQratApproximation(uint a, uint w) const;

    double getQB(uint a, uint tID) const;

    uint findWavelengthID(uint a, uint tIDnew, double rnd) const;

    double findTheta(uint a, uint w, double rnd) const;

    uint findSizeID(photon_package * pp, prob_list & prob, double a_min, double a_max, CRandomGenerator * rand_gen) const;

    double getAvgMass(CGridBasic * grid, const cell_basic & cell) const;

    double getAvgMass() const;

    double getPahMass(uint a);

    double getEnthalpy(uint a, uint t) const;

    double getEnthalpyBinWidth(uint a, uint t) const;

    double getCalorimetricTemperature(uint i) const;

    double getMinCalorimetricTemperature() const;

    double getMaxCalorimetricTemperature() const;

    uint getNrOfCalorimetryTemperatures() const;

    uint getNrOfDustSpecies() const;

    void setNrOfDustSpecies(uint val);

    void setNrOfWavelength(uint val);

    uint getNrOfWavelength() const;

    uint getWavelengthID(double wavelength);

    uint getNrOfIncidentAngles() const;

    uint getNrOfScatTheta(uint a, uint w) const;

    uint ** getNrOfScatTheta();

    uint getNrOfScatPhi() const;

    uint getNrOfScatMatElements() const;

    string getStringID() const;

    void createStringID(CDustComponent * comp);

    void setIndividualDustMassFractions(bool val);

    bool getIndividualDustMassFractions();

    void setWavelengthList(dlist _wavelength_list, uint _wavelength_offset);

    bool calcWavelengthDiff();

    void updateStokesVector(photon_package * pp, uint wnew) const;

    double getEffectiveRadius(uint a) const;

    double * getEffectiveRadii();

    double getGrainDistributionXRadiusSq(uint a) const;

    double getGrainSizeDistribution(uint a) const;

    double * getGrainSizeDistribution();

    double getEffectiveRadiusSquared(uint a) const;

    double getMass(uint a) const;

    double getVolume(uint a) const;

    double getMaterialDensity(uint a) const;

    double getMaterialDensity() const;

    double getFHighJ() const;

    double getFcorr() const;

    bool isAligned() const;

    void setIsAligned(bool val);

    void setIsMixture(bool val);

    double getSublimationTemperature();

    void getQtrq(uint i, uint j, double & x, double & y);

    void getHG_g_factor(uint i, uint j, double & x, double & y);

    void getHG_g2_factor(uint i, uint j, double & x, double & y);

    void getHG_g3_factor(uint i, uint j, double & x, double & y);

    bool getCalorimetryLoaded();

    void setCalorimetryLoaded(bool val);

    void setIDs(uint i_comp, uint nr_components, uint i_mix, uint nr_mixtures);

    void printIDs();

    bool getScatLoaded();

    void setScatLoaded(bool val);

    void setSublimate(bool val);

    uint getComponentId();

    uint getNrOfComponents();

    void SetNrOfScatTheta(uint ** nr_of_scat_theta_tmp);

    void SetScatTheta(double *** scat_theta_tmp);

    double * getScatTheta(uint a, uint w);

    double getScatTheta(uint a, uint w, uint sth) const;

    uint getScatThetaID(double theta, uint a, uint w) const;

    void preCalcEffProperties(parameters & param);

    void henyeygreen(photon_package * pp, uint a, CRandomGenerator * rand_gen);
    void drainehenyeygreen(photon_package * pp, uint a, CRandomGenerator * rand_gen);
    void threeparamhenyeygreen(photon_package * pp, uint a, CRandomGenerator * rand_gen);
    void miesca(photon_package * pp, uint a, CRandomGenerator * rand_gen);

    void preCalcTemperatureLists(double _minTemp, double _maxTemp, uint _nr_of_temperatures);
    void preCalcAbsorptionRates();
    void preCalcWaveProb();
    void preCalcMieScatteringProb();
    void preCalcCalorimetry();

    void convertTempInQB(CGridBasic * grid,
                         cell_basic * cell,
                         uint i_density,
                         double min_gas_density,
                         bool use_gas_temp);
    bool adjustTempAndWavelengthBW(CGridBasic * grid,
                                   photon_package * pp,
                                   uint i_density,
                                   bool use_energy_density,
                                   CRandomGenerator * rand_gen);
    double updateDustTemperature(CGridBasic * grid,
                                 const photon_package & pp,
                                 uint i_density,
                                 uint a,
                                 bool use_energy_density);
    void calcTemperature(CGridBasic * grid, cell_basic * cell, uint i_density, bool use_energy_density);
    void calcStochasticHeatingPropabilities(CGridBasic * grid,
                                            cell_basic * cell,
                                            uint i_density,
                                            dlist & wl_list);

    void calcAlignedRadii(CGridBasic * grid, cell_basic * cell, uint i_density);

    void initDustProperties();
    void initScatteringMatrixArray();
    void initScatteringMatrixArray(uint nr_of_scat_theta);
    void initNrOfScatThetaArray();
    void initScatThetaArray();
    void initCalorimetry();

    bool readDustParameterFile(parameters & param, uint dust_component_choice);
    bool readDustRefractiveIndexFile(parameters & param,
                                     uint dust_component_choice,
                                     double a_min_mixture,
                                     double a_max_mixture);
    bool readScatteringMatrices(string path, uint nr_of_wavelength_dustcat, dlist wavelength_list_dustcat);
    bool readCalorimetryFile(parameters & param, uint dust_component_choice);

    bool writeComponentData(string path_data);
    bool writeComponentPlot(string path_plot);
    bool calcSizeDistribution(dlist values, double * mass);
    bool add(double ** size_fraction, CDustComponent * comp, uint ** nr_of_scat_theta_tmp, double *** scat_theta_tmp);

    uint getInteractingDust(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen, uint cross_section = CROSS_ABS) const;

    void calcPACrossSections(uint a, uint w, cross_sections & cs, double theta) const;
    void calcNONPACrossSections(uint a, uint w, cross_sections & cs, double theta, double Rent) const;
    void calcExtCrossSections(CGridBasic * grid,
                              const photon_package & pp,
                              uint i_density,
                              double * avg_Cext,
                              double * avg_Cpol,
                              double * avg_Ccirc) const;
    void calcCrossSections(CGridBasic * grid,
                           const photon_package & pp,
                           uint i_density,
                           uint a,
                           double mag_field_theta,
                           cross_sections & cs) const;
    double calcGoldReductionFactor(const Vector3D & v, const Vector3D & B) const;

    void calcEmissivityHz(CGridBasic * grid,
                          const photon_package & pp,
                          uint i_density,
                          StokesVector * dust_emissivity) const;
    double calcEmissivity(CGridBasic * grid, const photon_package & pp, uint i_density) const;
    StokesVector calcEmissivityEmi(CGridBasic * grid,
                                   const photon_package & pp,
                                   uint i_density,
                                   uint emission_component,
                                   double phi,
                                   double energy,
                                   Vector3D en_dir) const;

    double getCalorimetryA(uint a, uint f, uint i, const spline & abs_rate_per_wl) const;
    long double * getStochasticProbability(uint a, const spline & abs_rate_per_wl) const;

    void getEscapePhoton(CGridBasic * grid,
                         photon_package * pp,
                         uint a,
                         Vector3D obs_ex,
                         Vector3D dir_obs,
                         photon_package * pp_escape) const;
    void getEscapePhotonMie(CGridBasic * grid,
                            photon_package * pp,
                            uint a,
                            Vector3D obs_ex,
                            Vector3D dir_obs,
                            photon_package * pp_escape) const;
    double getCellEmission(CGridBasic * grid, const photon_package & pp, uint i_density) const;

private:
    interp ** avg_scattering_frac;
    interp ** phase_pdf;

    spline *HG_g_factor, *HG_g2_factor, *HG_g3_factor, *Qtrq;
    spline * tab_planck;
    spline * tab_em;
    spline * tab_em_inv;
    spline tab_temp, tab_em_eff;

    prob_list * avg_planck_frac;
    prob_list *dust_prob, *sca_prob, *abs_prob;

    Matrix2D ***** sca_mat;
    double **Qext1, **Qext2, **Qabs1, **Qabs2, **Qsca1, **Qsca2, **Qcirc, **HGg, **HGg2, **HGg3;
    double ** enthalpy;
    double *a_eff, *grain_distribution_x_aeff_sq, *grain_size_distribution, *a_eff_squared;
    double * calorimetry_temperatures;
    double * mass;
    double *tCext1, *tCext2, *tCabs1, *tCabs2, *tCsca1, *tCsca2, *tCcirc, *tHGg, *tHGg2, *tHGg3;
    double * relWeightTab;
    double **CextMean, **CabsMean, **CscaMean;

    string stringID;
    string size_keyword;

    double stochastic_heating_max_size;
    double delta0;
    double aspect_ratio;
    double sub_temp;
    double material_density;
    double dust_mass_fraction;
    // double min_temp;
    double max_temp;
    double min_a_alig, max_a_alig;

    // alignment paramaters
    bool is_align;
    double f_highJ;
    double f_cor;
    double larm_f;
    double gold_g_factor;
    double Q_ref;
    double alpha_Q;
    double R_rayleigh;

    double delta_rat;
    double mu;
    double fraction;
    double avg_mass;
    double a_min_global;
    double a_max_global;

    bool dust_offset;
    bool scat_loaded, calorimetry_loaded;
    bool sublimate;
    bool is_mixture;
    bool individual_dust_fractions;

    int elements[16];

    uint i_component, nr_of_components;
    uint i_mixture, nr_of_mixtures;
    uint nr_of_wavelength;
    uint wavelength_offset;
    uint calorimetry_type;
    uint alignment;
    uint phID;
    uint nr_of_dust_species;
    uint nr_of_incident_angles;
    uint **nr_of_scat_theta;
    double ***scat_theta;
    uint nr_of_scat_phi;
    uint nr_of_scat_mat_elements;
    uint nr_of_calorimetry_temperatures;

    dlist wavelength_list;
    dlist wavelength_diff;
    dlist size_parameter;

    uilist dust_mixtures;
    uilist dust_choices_to_index;
    
    uint component_id;
};

#endif /* CDUST_COMPONENT_H */
