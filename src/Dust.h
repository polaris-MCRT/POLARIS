#pragma once
#include "Grid.h"
#include "Matrix2D.h"
#include "Typedefs.h"
#include "Cell.h"
#include "MathFunctions.h"
#include "Photon.h"
#include "Stokes.h"
#include "Vector.h"

class parameters;

#ifndef CDUST
#define CDUST

class CDustComponent
{
  public:
    CDustComponent()
    {
        avg_scattering_frac = 0;
        phase_pdf = 0;

        HG_g_factor = 0;
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
        a_eff = 0;
        a_eff_1_5 = 0;
        a_eff_3_5 = 0;
        a_eff_2 = 0;
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
        if(a_eff_1_5 != 0)
            delete[] a_eff_1_5;
        if(a_eff_3_5 != 0)
            delete[] a_eff_3_5;
        if(a_eff_2 != 0)
            delete[] a_eff_2;
        if(mass != 0)
            delete[] mass;
        if(relWeightTab != 0)
            delete[] relWeightTab;
    }

    // ----------------------------------------------------------------------
    // ----------- Efficiencies for grain size and wavelength ---------------
    // ----------------------------------------------------------------------
    inline double getQext1(uint a, uint w) const
    {
        return Qext1[a][w];
    }

    inline double getQext2(uint a, uint w) const
    {
        return Qext2[a][w];
    }

    inline double getQabs1(uint a, uint w) const
    {
        return Qabs1[a][w];
    }

    inline double getQabs2(uint a, uint w) const
    {
        return Qabs2[a][w];
    }

    inline double getQsca1(uint a, uint w) const
    {
        return Qsca1[a][w];
    }

    inline double getQsca2(uint a, uint w) const
    {
        return Qsca2[a][w];
    }

    inline double getQcirc(uint a, uint w) const
    {
        return Qcirc[a][w];
    }

    inline double getHGg(uint a, uint w) const
    {
        return HGg[a][w];
    }

    // ------------------------------------------------------------------------------------
    // ----------- Add values to efficiencies for grain size and wavelength
    // ---------------
    // ------------------------------------------------------------------------------------
    void addQext1(uint a, uint w, double val)
    {
        Qext1[a][w] += val;
    }

    void addQext2(uint a, uint w, double val)
    {
        Qext2[a][w] += val;
    }

    void addQabs1(uint a, uint w, double val)
    {
        Qabs1[a][w] += val;
    }

    void addQabs2(uint a, uint w, double val)
    {
        Qabs2[a][w] += val;
    }

    void addQsca1(uint a, uint w, double val)
    {
        Qsca1[a][w] += val;
    }

    void addQsca2(uint a, uint w, double val)
    {
        Qsca2[a][w] += val;
    }

    void addQcirc(uint a, uint w, double val)
    {
        Qcirc[a][w] += val;
    }

    void addHGg(uint a, uint w, double val)
    {
        HGg[a][w] += val;
    }

    // ------------------------------------------------------------------------
    // ----------- Cross-sections for grain size and wavelength ---------------
    // ------------------------------------------------------------------------
    double getCext1(uint a, uint w) const
    {
        return PI * a_eff_2[a] * getQext1(a, w);
    }

    double getCext2(uint a, uint w) const
    {
        return PI * a_eff_2[a] * getQext2(a, w);
    }

    double getCabs1(uint a, uint w) const
    {
        return PI * a_eff_2[a] * getQabs1(a, w);
    }

    double getCabs2(uint a, uint w) const
    {
        return PI * a_eff_2[a] * getQabs2(a, w);
    }

    double getCsca1(uint a, uint w) const
    {
        return PI * a_eff_2[a] * getQsca1(a, w);
    }

    double getCsca2(uint a, uint w) const
    {
        return PI * a_eff_2[a] * getQsca2(a, w);
    }

    double getCcirc(uint a, uint w) const
    {
        return PI * a_eff_2[a] * getQcirc(a, w);
    }

    // -------------------------------------------------------------------------------
    // ----------- Cross-sections for wavelength mixed in current cell ---------------
    // -------------------------------------------------------------------------------
    double getCext1(CGridBasic * grid, const photon_package & pp) const
    {
        // Get wavelength of photon package
        uint w = pp.getDustWavelengthID();

        // Return precalculated value if available
        if(tCext1 != 0)
            return tCext1[w];

        // Get local min and max grain sizes
        double a_min = getSizeMin(grid, pp);
        double a_max = getSizeMax(grid, pp);

        // Get local size parameter for size distribution
        double size_param = getSizeParam(grid, pp);

        // Get integration over the dust size distribution
        double * rel_weight = getRelWeight(a_min, a_max, size_param);

        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] *= a_eff_2[a] * getQext1(a, w);
        double res =
            PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
        delete[] rel_weight;
        return res;
    }

    double getCext2(CGridBasic * grid, const photon_package & pp) const
    {
        // Get wavelength of photon package
        uint w = pp.getDustWavelengthID();

        // Return precalculated value if available
        if(tCext2 != 0)
            return tCext2[w];

        // Get local min and max grain sizes
        double a_min = getSizeMin(grid, pp);
        double a_max = getSizeMax(grid, pp);

        // Get local size parameter for size distribution
        double size_param = getSizeParam(grid, pp);

        // Get integration over the dust size distribution
        double * rel_weight = getRelWeight(a_min, a_max, size_param);

        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] *= a_eff_2[a] * getQext2(a, w);
        double res =
            PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
        delete[] rel_weight;
        return res;
    }

    double getCabs1(CGridBasic * grid, const photon_package & pp) const
    {
        // Get wavelength of photon package
        uint w = pp.getDustWavelengthID();

        // Return precalculated value if available
        if(tCabs1 != 0)
            return tCabs1[w];

        // Get local min and max grain sizes
        double a_min = getSizeMin(grid, pp);
        double a_max = getSizeMax(grid, pp);

        // Get local size parameter for size distribution
        double size_param = getSizeParam(grid, pp);

        // Get integration over the dust size distribution
        double * rel_weight = getRelWeight(a_min, a_max, size_param);

        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] *= a_eff_2[a] * getQabs1(a, w);
        double res =
            PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
        delete[] rel_weight;
        return res;
    }

    double getCabs2(CGridBasic * grid, const photon_package & pp) const
    {
        // Get wavelength of photon package
        uint w = pp.getDustWavelengthID();

        // Return precalculated value if available
        if(tCabs2 != 0)
            return tCabs2[w];

        // Get local min and max grain sizes
        double a_min = getSizeMin(grid, pp);
        double a_max = getSizeMax(grid, pp);

        // Get local size parameter for size distribution
        double size_param = getSizeParam(grid, pp);

        // Get integration over the dust size distribution
        double * rel_weight = getRelWeight(a_min, a_max, size_param);

        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] *= a_eff_2[a] * getQabs2(a, w);
        double res =
            PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
        delete[] rel_weight;
        return res;
    }

    double getCsca1(CGridBasic * grid, const photon_package & pp) const
    {
        // Get wavelength of photon package
        uint w = pp.getDustWavelengthID();

        // Return precalculated value if available
        if(tCsca1 != 0)
            return tCsca1[w];

        // Get local min and max grain sizes
        double a_min = getSizeMin(grid, pp);
        double a_max = getSizeMax(grid, pp);

        // Get local size parameter for size distribution
        double size_param = getSizeParam(grid, pp);

        // Get integration over the dust size distribution
        double * rel_weight = getRelWeight(a_min, a_max, size_param);

        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] *= a_eff_2[a] * getQsca1(a, w);
        double res =
            PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
        delete[] rel_weight;
        return res;
    }

    double getCsca2(CGridBasic * grid, const photon_package & pp) const
    {
        // Get wavelength of photon package
        uint w = pp.getDustWavelengthID();

        // Return precalculated value if available
        if(tCsca2 != 0)
            return tCsca2[w];

        // Get local min and max grain sizes
        double a_min = getSizeMin(grid, pp);
        double a_max = getSizeMax(grid, pp);

        // Get local size parameter for size distribution
        double size_param = getSizeParam(grid, pp);

        // Get integration over the dust size distribution
        double * rel_weight = getRelWeight(a_min, a_max, size_param);

        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] *= a_eff_2[a] * getQsca2(a, w);
        double res =
            PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
        delete[] rel_weight;
        return res;
    }

    double getCcirc(CGridBasic * grid, const photon_package & pp) const
    {
        // Get wavelength of photon package
        uint w = pp.getDustWavelengthID();

        // Return precalculated value if available
        if(tCcirc != 0)
            return tCcirc[w];

        // Get local min and max grain sizes
        double a_min = getSizeMin(grid, pp);
        double a_max = getSizeMax(grid, pp);

        // Get local size parameter for size distribution
        double size_param = getSizeParam(grid, pp);

        // Get integration over the dust size distribution
        double * rel_weight = getRelWeight(a_min, a_max, size_param);

        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] *= a_eff_2[a] * getQcirc(a, w);
        double res =
            PI * CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
        delete[] rel_weight;
        return res;
    }

    double getHGg(CGridBasic * grid, const photon_package & pp) const
    {
        // Get wavelength of photon package
        uint w = pp.getDustWavelengthID();

        // Return precalculated value if available
        if(tHGg != 0)
            return tHGg[w];

        // Get local min and max grain sizes
        double a_min = getSizeMin(grid, pp);
        double a_max = getSizeMax(grid, pp);

        // Get local size parameter for size distribution
        double size_param = getSizeParam(grid, pp);

        // Get integration over the dust size distribution
        double * rel_weight = getRelWeight(a_min, a_max, size_param);

        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] *= getHGg(a, w);
        double res = CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);
        delete[] rel_weight;
        return res;
    }

    // -----------------------------------------------------------------------------
    // ----------- Average cross-sections for grain size and wavelength ------------
    // -----------------------------------------------------------------------------

    inline double getCextMean(uint a, uint w) const
    {
        return CextMean[a][w];
    }

    inline double getCabsMean(uint a, uint w) const
    {
        return CabsMean[a][w];
    }

    inline double getCscaMean(uint a, uint w) const
    {
        return CscaMean[a][w];
    }

    // -----------------------------------------------------------------------------
    // ----------- Average cross-sections for wavelength mixed in current cell -----
    // -----------------------------------------------------------------------------
    double getCextMean(CGridBasic * grid, const photon_package & pp) const
    {
        return (2.0 * getCext1(grid, pp) + getCext2(grid, pp)) / 3.0;
    }

    double getCabsMean(CGridBasic * grid, const photon_package & pp) const
    {
        return (2.0 * getCabs1(grid, pp) + getCabs2(grid, pp)) / 3.0;
    }

    double getCscaMean(CGridBasic * grid, const photon_package & pp) const
    {
        return (2.0 * getCsca1(grid, pp) + getCsca2(grid, pp)) / 3.0;
    }

    // ------------------------------------------------------------------------------------
    // ----------- Average cross-sections for wavelength mixed with global conditions -----
    // ------------------------------------------------------------------------------------
    double getCextMean(double w) const
    {
        return (2.0 * getCext1(w) + getCext2(w)) / 3.0;
    }

    double getCabsMean(double w) const
    {
        return (2.0 * getCabs1(w) + getCabs2(w)) / 3.0;
    }

    double getCscaMean(double w) const
    {
        return (2.0 * getCsca1(w) + getCsca2(w)) / 3.0;
    }

    // ------------------------------------------------------------------------------
    // ----------- Cross-sections for wavelength mixed with global limits -----------
    // ------------------------------------------------------------------------------
    double getCext1(uint w) const
    {
        double * Cext1 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Cext1[a] = a_eff_1_5[a] * getQext1(a, w);
        double res =
            PI / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Cext1, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Cext1;
        return res;
    }

    double getCext2(uint w) const
    {
        double * Cext2 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Cext2[a] = a_eff_1_5[a] * getQext2(a, w);
        double res =
            PI / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Cext2, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Cext2;
        return res;
    }

    double getCabs1(uint w) const
    {
        double * Cabs1 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Cabs1[a] = a_eff_1_5[a] * getQabs1(a, w);
        double res =
            PI / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Cabs1, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Cabs1;
        return res;
    }

    double getCabs2(uint w) const
    {
        double * Cabs2 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Cabs2[a] = a_eff_1_5[a] * getQabs2(a, w);
        double res =
            PI / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Cabs2, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Cabs2;
        return res;
    }

    double getCsca1(uint w) const
    {
        double * Csca1 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Csca1[a] = a_eff_1_5[a] * getQsca1(a, w);
        double res =
            PI / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Csca1, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Csca1;
        return res;
    }

    double getCsca2(uint w) const
    {
        double * Csca2 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Csca2[a] = a_eff_1_5[a] * getQsca2(a, w);
        double res =
            PI / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Csca2, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Csca2;
        return res;
    }

    double getCcirc(uint w) const
    {
        double * Ccirc = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Ccirc[a] = a_eff_1_5[a] * getQcirc(a, w);
        double res =
            PI / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Ccirc, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Ccirc;
        return res;
    }

    double getHGg(uint w) const
    {
        double * HGg = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            HGg[a] = a_eff_3_5[a] * getHGg(a, w);
        double res =
            1.0 / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, HGg, nr_of_dust_species, a_min_global, a_max_global);
        delete[] HGg;
        return res;
    }

    // ------------------------------------------------------------------------------
    // ----------- Efficiencies for wavelength mixed with global limits -----------
    // ------------------------------------------------------------------------------
    double getQext1(uint w) const
    {
        double * Qext1 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Qext1[a] = a_eff_3_5[a] * getQext1(a, w);
        double res =
            1.0 / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Qext1, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Qext1;
        return res;
    }

    double getQext2(uint w) const
    {
        double * Qext2 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Qext2[a] = a_eff_3_5[a] * getQext2(a, w);
        double res =
            1.0 / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Qext2, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Qext2;
        return res;
    }

    double getQabs1(uint w) const
    {
        double * Qabs1 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Qabs1[a] = a_eff_3_5[a] * getQabs1(a, w);
        double res =
            1.0 / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Qabs1, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Qabs1;
        return res;
    }

    double getQabs2(uint w) const
    {
        double * Qabs2 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Qabs2[a] = a_eff_3_5[a] * getQabs2(a, w);
        double res =
            1.0 / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Qabs2, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Qabs2;
        return res;
    }

    double getQsca1(uint w) const
    {
        double * Qsca1 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Qsca1[a] = a_eff_3_5[a] * getQsca1(a, w);
        double res =
            1.0 / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Qsca1, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Qsca1;
        return res;
    }

    double getQsca2(uint w) const
    {
        double * Qsca2 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Qsca2[a] = a_eff_3_5[a] * getQsca2(a, w);
        double res =
            1.0 / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Qsca2, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Qsca2;
        return res;
    }

    double getQcirc(uint w) const
    {
        double * Qcirc = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Qcirc[a] = a_eff_3_5[a] * getQcirc(a, w);
        double res =
            1.0 / getWeight() *
            CMathFunctions::integ_dust_size(a_eff, Qcirc, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Qcirc;
        return res;
    }

    // -----------------------------------------------------------------------------------
    // ----------- Mass cross-sections for wavelength mixed with global limits -----------
    // -----------------------------------------------------------------------------------
    double getKappaExt1(uint w) const
    {
        double * Cext1 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Cext1[a] = a_eff_1_5[a] * getQext1(a, w);
        double res =
            PI / getWeight() / getAvgMass() *
            CMathFunctions::integ_dust_size(a_eff, Cext1, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Cext1;
        return res;
    }

    double getKappaExt2(uint w) const
    {
        double * Cext2 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Cext2[a] = a_eff_1_5[a] * getQext2(a, w);
        double res =
            PI / getWeight() / getAvgMass() *
            CMathFunctions::integ_dust_size(a_eff, Cext2, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Cext2;
        return res;
    }

    double getKappaAbs1(uint w) const
    {
        double * Cabs1 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Cabs1[a] = a_eff_1_5[a] * getQabs1(a, w);
        double res =
            PI / getWeight() / getAvgMass() *
            CMathFunctions::integ_dust_size(a_eff, Cabs1, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Cabs1;
        return res;
    }

    double getKappaAbs2(uint w) const
    {
        double * Cabs2 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Cabs2[a] = a_eff_1_5[a] * getQabs2(a, w);
        double res =
            PI / getWeight() / getAvgMass() *
            CMathFunctions::integ_dust_size(a_eff, Cabs2, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Cabs2;
        return res;
    }

    double getKappaSca1(uint w) const
    {
        double * Csca1 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Csca1[a] = a_eff_1_5[a] * getQsca1(a, w);
        double res =
            PI / getWeight() / getAvgMass() *
            CMathFunctions::integ_dust_size(a_eff, Csca1, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Csca1;
        return res;
    }

    double getKappaSca2(uint w) const
    {
        double * Csca2 = new double[nr_of_dust_species];
        for(uint a = 0; a < nr_of_dust_species; a++)
            Csca2[a] = a_eff_1_5[a] * getQsca2(a, w);
        double res =
            PI / getWeight() / getAvgMass() *
            CMathFunctions::integ_dust_size(a_eff, Csca2, nr_of_dust_species, a_min_global, a_max_global);
        delete[] Csca2;
        return res;
    }

    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------

    double getDeltaRat() const
    {
        return delta_rat;
    }

    double * getRelWeight(double a_min, double a_max, double size_param = 0) const
    {
        double * rel_weight = new double[nr_of_dust_species];

        if(getRelWeightTab(0) != MAX_DOUBLE)
            for(uint a = 0; a < nr_of_dust_species; a++)
                rel_weight[a] = getRelWeightTab(a);
        else
        {
            // Create normalization factor
            for(uint a = 0; a < nr_of_dust_species; a++)
                rel_weight[a] = a_eff_3_5[a] * pow(a_eff[a], size_param);

            double weight = CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);

            // Create the final relative mass distribution
            for(uint a = 0; a < nr_of_dust_species; a++)
            {
                rel_weight[a] /= weight;
            }
        }

        return rel_weight;
    }

    double getRelWeightTab(uint a) const
    {
        if(relWeightTab != 0)
            return relWeightTab[a];
        return MAX_DOUBLE;
    }

    void preCalcRelWeight()
    {
        relWeightTab = new double[nr_of_dust_species];
        fill(relWeightTab,relWeightTab+nr_of_dust_species,MAX_DOUBLE);
        relWeightTab = getRelWeight(a_min_global, a_max_global);
    }

    double getWeight() const
    {
        double weight =
            CMathFunctions::integ_dust_size(a_eff, a_eff_3_5, nr_of_dust_species, a_min_global, a_max_global);
        return weight;
    }

    double getMassWeight()
    {
        // Init pointer array of relative mass of each dust grain size bin
        double * rel_mass = new double[nr_of_dust_species];

        // Set the relative mass of each dust grain size bin
        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_mass[a] = a_eff_3_5[a] * mass[a];

        // Calculate the average mass of the dust grains in the current cell
        double mass_weight =
            CMathFunctions::integ_dust_size(a_eff, rel_mass, nr_of_dust_species, a_min_global, a_max_global);

        // Delete pointer array
        delete[] rel_mass;

        return mass_weight;
    }

    double getGoldFactor()
    {
        return gold_g_factor;
    }

    double getLarmF()
    {
        return larm_f;
    }

    double getQref()
    {
        return Q_ref;
    }

    double getAlphaQ()
    {
        return alpha_Q;
    }

    double getRayleighReductionFactor()
    {
        return R_rayleigh;
    }

    void setAlignmentMechanism(uint al)
    {
        alignment = al;
    }

    void setPhaseFunctionID(uint ph)
    {
        phID = ph;
    }

    uint getPhaseFunctionID()
    {
        return phID;
    }

    uint getNrOfStochasticSizes()
    {
        uint nr_stochastic_sizes = 0;
        if(stochastic_heating_max_size > 0)
            for(uint a = 0; a < nr_of_dust_species; a++)
                if(a_eff[a] <= stochastic_heating_max_size)
                    nr_stochastic_sizes++;
        return nr_stochastic_sizes;
    }

    // double getMinDustTemp()
    // {
    //     return min_temp;
    // }

    double getMaxDustTemp()
    {
        return max_temp;
    }

    double getMinAlignedRadius()
    {
        return min_a_alig;
    }

    double getMaxAlignedRadius()
    {
        return max_a_alig;
    }

    double getScatteringMatrixElement(uint a,
                                      uint w,
                                      uint incID,
                                      uint sphID,
                                      uint sthID,
                                      uint i_mat,
                                      uint j_mat) const
    {
        if(sca_mat == 0)
            return 0;
            
        return sca_mat[a][w][incID][sphID][sthID](i_mat, j_mat);
    }

    const Matrix2D & getScatteringMatrix(uint a, uint w, uint incID, uint sphID, uint sthID) const
    {
        return sca_mat[a][w][incID][sphID][sthID];
    }

    void cleanScatteringData()
    {
        if(sca_mat != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
            {
                if(sizeIndexUsed(a))
                {
                    for(uint w = 0; w < nr_of_wavelength; w++)
                    {
                        if(nr_of_scat_theta[a][w] != 0)
                        {
                            for(uint inc = 0; inc < nr_of_incident_angles; inc++)
                            {
                                for(uint sph = 0; sph < nr_of_scat_phi; sph++)
                                    delete[] sca_mat[a][w][inc][sph];
                                delete[] sca_mat[a][w][inc];
                            }
                            delete[] sca_mat[a][w];
                        }
                    }
                    delete[] sca_mat[a];
                }
            }
            delete[] sca_mat;
            sca_mat = 0;
        }
    }

    StokesVector getRadFieldScatteredFraction(CGridBasic * grid,
                                              const photon_package & pp,
                                              uint i_density,
                                              const Vector3D & en_dir,
                                              double energy) const;

    double getScatteredFraction(uint a, uint w, double theta) const
    {
        double res = 1;

        switch(phID)
        {
            case PH_ISO:
                res = 1 / PIx4;
                break;

            case PH_HG:
                double g;
                g = getHGg(a, w);
                res = CMathFunctions::phaseFunctionHG(g, theta);
                break;

            case PH_MIE:
                res = getScatteredFractionMie(a, w, theta);
                break;
        }
        return res;
    }

    double getScatteredFractionMie(uint a, uint w, double theta) const
    {
        return phase_pdf[a][w].getValue(theta);
    }

    double getScatteredFractionMie(uint a, uint w, uint sth) const
    {
        return phase_pdf[a][w].getValue(sth);
    }

    void scatter(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen)
    {
        switch(phID)
        {
            case PH_HG:
            {
                uint a = getInteractingDust(grid, pp, rand_gen, CROSS_SCA);
                henyeygreen(pp, a, rand_gen);
                break;
            }

            case PH_MIE:
            {
                uint a = getInteractingDust(grid, pp, rand_gen, CROSS_SCA);
                miesca(pp, a, rand_gen);
                break;
            }

            default:
                pp->setRandomDirection(rand_gen->getRND(), rand_gen->getRND());
                pp->updateCoordSystem();
                break;
        }
    }

    double getStochasticHeatingMaxSize() const
    {
        return stochastic_heating_max_size;
    }

    void setStochasticHeatingMaxSize(double val)
    {
        stochastic_heating_max_size = val;
    }

    double getFraction() const
    {
        return fraction;
    }

    void setFraction(double val)
    {
        fraction = val;
    }

    double getDustMassFraction() const
    {
        return dust_mass_fraction;
    }

    void setDustMassFraction(double val)
    {
        dust_mass_fraction = val;
    }

    void setSizeParameter(string size_key, dlist size_parameter_list)
    {
        size_keyword = size_key;
        size_parameter = size_parameter_list;
    }

    string getDustSizeKeyword()
    {
        return size_keyword;
    }

    string getDustSizeParameterString()
    {
        stringstream str_stream;
        for(uint i = 0; i < NR_OF_SIZE_DIST_PARAM; i++)
            if(size_parameter[i] != 0)
            {
                if(i > 0)
                    str_stream << ", ";
                str_stream << size_parameter[i];
            }
            else
            {
                if(i == 0)
                    return "flat";
            }
        return str_stream.str();
    }

    double convDensityToNumber(CGridBasic * grid, const cell_basic & cell, bool from_gas = false) const
    {
        double avg_mass = getAvgMass(grid, cell);

        if(avg_mass == 0)
            return 0;

        double conversion_factor = 1.0;
        if(from_gas)
        {
            if(!grid->getGasIsMassDensity())
                conversion_factor *= grid->getMu() * m_H;
            conversion_factor *= getDustMassFraction();
            conversion_factor /= avg_mass;
        }
        else if(grid->getDustIsMassDensity())
            conversion_factor = 1 / avg_mass;

        return conversion_factor;
    }

    double convDensityToMass(CGridBasic * grid, const cell_basic & cell, bool from_gas = false) const
    {
        double conversion_factor = 1.0;
        if(from_gas)
        {
            if(!grid->getGasIsMassDensity())
                conversion_factor *= grid->getMu() * m_H;
            conversion_factor *= getDustMassFraction();
        }
        else if(!grid->getDustIsMassDensity())
            conversion_factor *= getAvgMass(grid, cell);

        return conversion_factor;
    }

    double getNumberDensity(CGridBasic * grid, const cell_basic & cell) const
    {
        if(grid->useDustDensities())
            return convDensityToNumber(grid, cell) * grid->getDustDensity(cell);
        else
            return convDensityToNumber(grid, cell, true) * grid->getGasDensity(cell);
    }

    double getNumberDensity(CGridBasic * grid, const photon_package & pp) const
    {
        return getNumberDensity(grid, *pp.getPositionCell());
    }

    double getNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const
    {
        if(grid->useDustDensities())
            return convDensityToNumber(grid, cell) * grid->getDustDensity(cell, i_density);
        else
            return convDensityToNumber(grid, cell, true) * grid->getGasDensity(cell, i_density);
    }

    double getNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_density) const
    {
        return getNumberDensity(grid, *pp.getPositionCell(), i_density);
    }

    double getMassDensity(CGridBasic * grid, const cell_basic & cell) const
    {
        if(grid->useDustDensities())
            return convDensityToMass(grid, cell) * grid->getDustDensity(cell);
        else
            return getDustMassFraction() * grid->getGasMassDensity(cell);
    }

    uint getMassDensity(CGridBasic * grid, const photon_package & pp) const
    {
        return getMassDensity(grid, *pp.getPositionCell());
    }

    double getMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const
    {
        if(grid->useDustDensities())
            return convDensityToMass(grid, cell) * grid->getDustDensity(cell, i_density);
        else
            return getDustMassFraction() * grid->getGasMassDensity(cell, i_density);
    }

    void setMu(double mu_)
    {
        mu = mu_;
    }

    void setMaterialDensity(double dens)
    {
        material_density = dens;
    }

    bool checkGrainSizeLimits(double a_min, double a_max)
    {
        if(a_min == 0)
        {
            // If no minimum grain size is set, start with the first one
            a_min_global = a_eff[0];
        }
        else if(a_min <= a_eff[nr_of_dust_species - 1])
        {
            // Check if inside boundaries
            a_min_global = a_min;
        }
        else
        {
            cout << "\nERROR: Minimum grain size must be smaller than " << a_eff[nr_of_dust_species - 1]
                 << "!        " << endl;
            return false;
        }

        if(a_max == 0)
        {
            // If no maximum grain size is set, end with the last one
            a_max_global = a_eff[nr_of_dust_species - 1];
        }
        else if(a_max >= a_eff[0])
        {
            // Check if inside boundaries
            a_max_global = a_max;
        }
        else
        {
            cout << "\nERROR: Maximum grain size must be larger than " << a_eff[0] << "!        " << endl;
            return false;
        }

        // If the minimum grain size is larger than the maximum, no calculation possible
        if(a_min_global > a_max_global)
        {
            cout << "\nERROR: Minimum grain size (" << a_min_global
                 << ") must be smaller than maximum grain size (" << a_max_global << ")!" << endl;
            return false;
        }

        return true;
    }

    double combinedRFactor(double R1, double R2, double R3) const
    {
        double res = (R1 + R2 + R3 + R1 * R2 + R1 * R3 + R2 * R3 + 3 * R1 * R2 * R3) /
                     (1 + 2 * R1 * R2 + 2 * R1 * R3 + 2 * R2 * R3 + 2 * R1 * R2 * R3);

        if(res > 1)
            res = 1;
        if(res < -0.5)
            res = -0.5;

        return res;
    }

    void setDelta0(double val)
    {
        delta0 = val;
    }

    void setLarmF(double val)
    {
        larm_f = val;
    }

    double getAspectRatio()
    {
        return aspect_ratio;
    }

    double getSizeParam(CGridBasic * grid, const cell_basic & cell) const
    {
        double size_param = grid->getGrainSizeParam(cell);
        if(size_param != 0)
            return size_param;
        return 0;
    }

    double getSizeParam(CGridBasic * grid, const photon_package & pp) const
    {
        return getSizeParam(grid, *pp.getPositionCell());
    }

    double getSizeParam() const
    {
        return 0;
    }

    double getSizeMin(CGridBasic * grid, const cell_basic & cell) const
    {
        double a_min = grid->getMinGrainRadius(cell);
        if(a_min > 0)
            return a_min;
        return a_min_global;
    }

    double getSizeMin(CGridBasic * grid, const photon_package & pp) const
    {
        return getSizeMin(grid, *pp.getPositionCell());
    }

    double getSizeMin() const
    {
        return a_min_global;
    }

    void setSizeMin(double val)
    {
        a_min_global = val;
    }

    double getSizeMax(CGridBasic * grid, const cell_basic & cell) const
    {
        double a_max = grid->getMaxGrainRadius(cell);
        if(a_max != 0 && a_max < nr_of_dust_species)
            return a_max;
        return a_max_global;
    }

    double getSizeMax(CGridBasic * grid, const photon_package & pp) const
    {
        return getSizeMax(grid, *pp.getPositionCell());
    }

    double getSizeMax() const
    {
        return a_max_global;
    }

    void setSizeMax(double val)
    {
        a_max_global = val;
    }

    bool sizeIndexUsed(uint a, double a_min, double a_max) const
    {
        if((a_eff[a] - a_min) > -a_min * 1e-5 && (a_max - a_eff[a]) > -a_max * 1e-5)
            return true;
        else
        {
            if(a < nr_of_dust_species - 1)
                if(a_eff[a] < a_min && a_eff[a + 1] >= a_min)
                    return true;
            if(a > 0)
                if(a_eff[a - 1] < a_max && a_eff[a] >= a_max)
                    return true;
        }
        return false;
    }

    bool sizeIndexUsed(uint a) const
    {
        if((a_eff[a] - a_min_global) > -a_min_global * 1e-5 &&
           (a_max_global - a_eff[a]) > -a_max_global * 1e-5)
            return true;
        else
        {
            if(a < nr_of_dust_species - 1)
                if(a_eff[a] < a_min_global && a_eff[a + 1] >= a_min_global)
                    return true;
            if(a > 0)
                if(a_eff[a - 1] < a_max_global && a_eff[a] >= a_max_global)
                    return true;
        }
        return false;
    }

    double getInternalIDG(double Td, double Tg) const
    {
        double alpha = sqrt((1 + 0.5 / (aspect_ratio * aspect_ratio)) * (1 + Td / Tg));
        double delta = 2 / (1 + aspect_ratio * aspect_ratio) - 1;

        double x = alpha * delta;
        if(x <= 1e-4)
            x = 1e-4;
        double sq_x = sqrt(x);
        double erfi = CMathFunctions::getErfi(sq_x);
        if(erfi <= 1e-4)
            erfi = 1e-4;

        double res = exp(x) / (PIsq * sq_x * erfi) - 1. / (2. * x);
        if(res != res)
            res = 1;
        if(res > 1)
            res = 1;
        if(res < 1e-5)
            res = 1e-5;

        return res;
    }

    double getInternalGOLD(double Td, double Tg, double vg) const
    {
        double h = 2 / (1 + aspect_ratio * aspect_ratio);
        double alpha = 0.5 / Td * (2 / h + 1) * (0.5 * (Td + Tg) + mu * m_H * vg * vg / (6 * con_kB));
        double delta = h - 1;

        double x = alpha * delta;
        if(x <= 1e-4)
            x = 1e-4;
        double sq_x = sqrt(x);
        double erfi = CMathFunctions::getErfi(sq_x);
        if(erfi <= 1e-4)
            erfi = 1e-4;

        double res = exp(x) / (PIsq * sq_x * erfi) - 1. / (2. * x);

        if(res != res)
            res = 1;
        if(res > 1)
            res = 1;
        if(res < 1e-5)
            res = 1e-5;

        return res;
    }

    double getInternalRAT() const
    {
        // x = delta, alpha = 1;
        double x = 2. / (1. + aspect_ratio * aspect_ratio) - 1.;
        if(x <= 1e-4)
            x = 1e-4;
        double sq_x = sqrt(x);
        double erfi = CMathFunctions::getErfi(sq_x);
        if(erfi <= 1e-4)
            erfi = 1e-4;

        double res = exp(x) / (PIsq * sq_x * erfi) - 1. / (2. * x);

        if(res != res)
            res = 1;
        if(res > 1)
            res = 1;
        if(res < 1e-5)
            res = 1e-5;

        return res;
    }

    double findTemperature(uint a, double qb) const
    {
        // Return temperature for energy (qb) value (grain size a)
        double temp = tab_em[a].getValue(qb);

        if(temp < TEMP_MIN)
            temp = double(TEMP_MIN);

        return temp;
    }

    double findTemperature(CGridBasic * grid, cell_basic * cell, double qb) const
    {
        if(tCabs1 != 0 && tCabs2 != 0)
        {
            // Get temperature from current absorbed energy
            return max(double(TEMP_MIN), tab_em_eff.getValue(qb));
        }

        // Get number of temperatures from tab_temp spline
        uint nr_of_temperatures = tab_temp.size();

        // Init spline for absorption/emission energy interpolation
        spline tmp_tab_em(nr_of_temperatures);

        // Init a temporary array for QB values
        double * tmpQB = new double[nr_of_wavelength];

        // Init a temporary array for planck values
        double * tmpCabs = new double[nr_of_wavelength];

        // init photon package for position and wavelength
        photon_package pp = photon_package();

        // Set position of the photon into current cell
        pp.setPositionCell(cell);

        for(uint w = 0; w < nr_of_wavelength; w++)
        {
            // Set wavelength of photon package
            pp.setWavelength(wavelength_list[w], w);

            // Pre calculate absorption cross-sections
            tmpCabs[w] = getCabsMean(grid, pp);
        }

        for(uint t = 0; t < nr_of_temperatures; t++)
        {
            // Get temperature from tab_temp spline
            double tmp_temp = tab_temp.getValue(t);

            // Calculate absorption cross-section times Planck function for each
            // wavelength
            for(uint w = 0; w < nr_of_wavelength; w++)
                tmpQB[w] = tmpCabs[w] * CMathFunctions::planck(wavelength_list[w], tmp_temp);

            // Calculate QB integrated over all wavelengths
            tmp_tab_em.setValue(
                t, CMathFunctions::integ(wavelength_list, tmpQB, 0, nr_of_wavelength - 1), tmp_temp);
        }

        // Delete multiple pointer
        delete[] tmpCabs;
        delete[] tmpQB;

        // Create spline for interpolation
        tmp_tab_em.createSpline();

        // Return temperature from current absorbed energy
        return max(double(TEMP_MIN), tmp_tab_em.getValue(qb));
    }

    uint findTemperatureID(double t) const
    {
        return tab_temp.getYIndex(t);
    }

    double getPlanck(uint w, double temp) const
    {
        double pl = CMathFunctions::planck(wavelength_list[w], temp);
        return max(1e-200, pl);
    }

    double getAbsRate(CGridBasic * grid, const cell_basic & cell, uint a, bool use_energy_density) const
    {
        double abs_rate = 0;
        if(use_energy_density)
        {
            double * sp_energy = new double[nr_of_wavelength];
            for(uint w = 0; w < nr_of_wavelength; w++)
                sp_energy[w] = (grid->getSpecLength(cell, w) * getCabsMean(a, w));
            abs_rate = CMathFunctions::integ(wavelength_list, sp_energy, 0, nr_of_wavelength - 1);

            delete[] sp_energy;
        }
        else
            for(uint w = 0; w < nr_of_wavelength; w++)
                abs_rate += grid->getSpecLength(cell, w) * getCabsMean(a, w);

        return abs(abs_rate) / (PIx4 * grid->getVolume(cell));
    }

    double getAbsRate(CGridBasic * grid, const photon_package & pp, uint a, bool use_energy_density) const
    {
        return getAbsRate(grid, *pp.getPositionCell(), a, use_energy_density);
    }

    void setFHighJ(double val)
    {
        f_highJ = val;
    }

    void setFcorr(double val)
    {
        f_cor = val;
    }

    void setQref(double val)
    {
        Q_ref = val;
    }

    void setAlphaQ(double val)
    {
        alpha_Q = val;
    }

    void setRayleighReductionFactor(double val)
    {
        R_rayleigh = val;
    }

    double getQrat(uint a, uint w, double theta) const
    {
        if(theta < 0)
            theta = 0;
        if(theta > PI2)
            theta = PI - theta;
        if(theta > PI2)
            theta = PI2;

        theta += PI2;

        double res = Qtrq[w * nr_of_dust_species + a].getValue(theta);

        if(res <= 0)
            res = getQratApproximation(a, w);
        if(res > 3)
            res = 3;
        return res;
    }

    double getQratApproximation(uint a, uint w) const
    {
        double ww = wavelength_list[w];
        double aa = getEffectiveRadius(a);

        if(ww / aa < 1)
            return 1;

        double res = pow(ww / aa, -3.0);
        return res;
    }

    double getQB(uint a, uint tID) const
    {
        return tab_em[a].getX(tID);
    }

    uint findWavelengthID(uint a, uint tIDnew, double rnd) const
    {
        return avg_planck_frac[tIDnew * nr_of_dust_species + a].getIndex(rnd);
    }

    double findTheta(uint a, uint w, double rnd) const
    {
        return avg_scattering_frac[a][w].getValue(rnd);
    }

    uint findSizeID(photon_package * pp, prob_list & prob, double a_min, double a_max, CRandomGenerator * rand_gen) const
    {
        uint a;

        // Ensure that the grain size index is one of the used ones
        while(true)
        {
            // Find the related grain size index
            a = prob.getIndex(rand_gen->getRND());

            // If in case the size index is outside of the used grain sizes, pick another
            // one
            if(sizeIndexUsed(a, a_min, a_max))
                break;
        }

        return a;
    }

    double getAvgMass(CGridBasic * grid, const cell_basic & cell) const
    {
        if(avg_mass != 0)
            return avg_mass;

        // Get local min and max grain sizes
        double a_min = getSizeMin(grid, cell);
        double a_max = getSizeMax(grid, cell);

        // Get local size parameter for size distribution
        double size_param = getSizeParam(grid, cell);

        // Get integration over the dust size distribution
        double * rel_weight = getRelWeight(a_min, a_max, size_param);

        // Create the final mass distribution
        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            rel_weight[a] *= mass[a];
        }

        // Calculate the average mass of the dust grains in the current cell
        double avg_mass =
            CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);

        // Delete pointer array
        delete[] rel_weight;

        return avg_mass;
    }

    double getAvgMass() const
    {
        // Get integration over the dust size distribution
        double weight = getWeight();

        // Init pointer array of relative mass of each dust grain size bin
        double * rel_mass = new double[nr_of_dust_species];

        // Set the relative mass of each dust grain size bin
        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_mass[a] = a_eff_3_5[a] * mass[a];

        // Calculate the average mass of the dust grains in the current cell
        double avg_mass =
            1 / weight *
            CMathFunctions::integ_dust_size(a_eff, rel_mass, nr_of_dust_species, a_min_global, a_max_global);

        // Delete pointer array
        delete[] rel_mass;

        return avg_mass;
    }

    double getPahMass(uint a)
    {
        // See Draine and Lee (2000)
        // Infrared Emission from Interstellar Dust. I. Stochastic Heating of Small Grains

        // Calculate the amount of carbon atoms of the PAH molecule
        double nr_carbon_atoms = round(pow(a_eff[a] * 1e9, 3) * 468.0);

        // Calculate the amount of hydrogen atoms of the PAH molecule
        double nr_hydrogen_atoms;
        if(nr_carbon_atoms <= 25)
            nr_hydrogen_atoms = floor(0.5 * nr_carbon_atoms + 0.5);
        else if(nr_carbon_atoms <= 100)
            nr_hydrogen_atoms = floor(2.5 * sqrt(nr_carbon_atoms) + 0.5);
        else
            nr_hydrogen_atoms = floor(0.25 * nr_carbon_atoms + 0.5);

        return (12 * nr_carbon_atoms + nr_hydrogen_atoms) * m_H;
    }

    double getEnthalpy(uint a, uint t) const
    {
        return enthalpy[a][t];
    }

    double getEnthalpyBinWidth(uint a, uint t) const
    {
        double diff = 0;
        if(t == nr_of_calorimetry_temperatures - 1)
            diff = enthalpy[a][t] - enthalpy[a][t - 1];
        else if(t == 0)
            diff = enthalpy[a][t + 1] - enthalpy[a][t];
        else
            diff = (enthalpy[a][t + 1] - enthalpy[a][t - 1]) / 2.0;
        return diff;
    }

    double getCalorimetricTemperature(uint i) const
    {
        return calorimetry_temperatures[i];
    }

    double getMinCalorimetricTemperature() const
    {
        return calorimetry_temperatures[0];
    }

    double getMaxCalorimetricTemperature() const
    {
        return calorimetry_temperatures[nr_of_calorimetry_temperatures - 1];
    }

    uint getNrOfCalorimetryTemperatures() const
    {
        return nr_of_calorimetry_temperatures;
    }

    uint getNrOfDustSpecies() const
    {
        return nr_of_dust_species;
    }

    uint getNrOfIncidentAngles() const
    {
        return nr_of_incident_angles;
    }

    uint getNrOfScatTheta(uint a, uint w) const
    {
        return nr_of_scat_theta[a][w];
    }

    uint ** getNrOfScatTheta()
    {
        return nr_of_scat_theta;
    }

    uint getNrOfScatPhi() const
    {
        return nr_of_scat_phi;
    }

    uint getNrOfScatMatElements() const
    {
        return nr_of_scat_mat_elements;
    }

    string getStringID() const
    {
        return stringID;
    }

    void createStringID(CDustComponent * comp)
    {
        // Init string stream
        stringstream str_stream;
        str_stream.str("");

        // Start with printing the dust components (if mixture or only single component)
        if(stringID.length() == 0 || comp->getNrOfComponents() == 1)
            str_stream << "Dust components:" << endl;

        // Use the fraction as dust to gas mass ratio (if user chosen it)
        double fraction = comp->getFraction();
        string fraction_string = "mass ratio";
        if(individual_dust_fractions)
        {
            fraction *= dust_mass_fraction;
            fraction_string = "dust-to-gas mass ratio";
        }

        // Fill the string with various parameters and format it
        char tmp_str[1024];
#ifdef WINDOWS
        sprintf_s(tmp_str,
                  "- %s\n    %s: %g, size distr. : \"%s\" (%s), size: %g [m] - %g [m]\n",
                  comp->getStringID().c_str(),
                  fraction_string.c_str(),
                  fraction,
                  comp->getDustSizeKeyword().c_str(),
                  comp->getDustSizeParameterString().c_str(),
                  comp->getSizeMin(),
                  comp->getSizeMax());
#else
        sprintf(tmp_str,
                "- %s\n    %s: %g, size distr. : \"%s\" (%s), size: %g [m] - %g [m]\n",
                comp->getStringID().c_str(),
                fraction_string.c_str(),
                fraction,
                comp->getDustSizeKeyword().c_str(),
                comp->getDustSizeParameterString().c_str(),
                comp->getSizeMin(),
                comp->getSizeMax());
#endif

        // Add formatted string to stream
        str_stream << tmp_str;

        // Add stream to stringID or replace it if only one component
        if(comp->getNrOfComponents() == 1)
            stringID = str_stream.str();
        else
            stringID += str_stream.str();
    }

    void setIndividualDustMassFractions(bool val)
    {
        individual_dust_fractions = val;
    }

    bool getIndividualDustMassFractions()
    {
        return individual_dust_fractions;
    }

    void setWavelengthList(dlist _wavelength_list, uint _wavelength_offset)
    {
        wavelength_list = _wavelength_list;
        nr_of_wavelength = wavelength_list.size();
        wavelength_offset = _wavelength_offset;
    }

    bool calcWavelengthDiff()
    {
        // Set width of each wavelength bin
        wavelength_diff.resize(nr_of_wavelength);

        // Differences need at least two wavelengths
        if(nr_of_wavelength == 1)
            return false;

        for(uint w = 0; w < nr_of_wavelength; w++)
        {
            if(w == nr_of_wavelength - 1)
                wavelength_diff[w] = wavelength_list[w] - wavelength_list[w - 1];
            else if(w == 0)
                wavelength_diff[w] = wavelength_list[w + 1] - wavelength_list[w];
            else
                wavelength_diff[w] = (wavelength_list[w + 1] - wavelength_list[w - 1]) / 2.0;
        }
        return true;
    }

    void updateStokesVector(photon_package * pp, uint wnew) const
    {
        // Get wavelength of photon package
        uint w = pp->getDustWavelengthID();
        *pp->getStokesVector() *= wavelength_diff[w] / wavelength_diff[wnew];
    }

    double getEffectiveRadius(uint a) const
    {
        return a_eff[a];
    }

    double * getEffectiveRadii()
    {
        return a_eff;
    }

    double getEffectiveRadius1_5(uint a) const
    {
        return a_eff_1_5[a];
    }

    double getEffectiveRadius3_5(uint a) const
    {
        return a_eff_3_5[a];
    }

    double * getEffectiveRadii3_5()
    {
        return a_eff_3_5;
    }

    double getEffectiveRadius_2(uint a) const
    {
        return a_eff_2[a];
    }

    double getMass(uint a) const
    {
        return mass[a];
    }

    double getVolume(uint a) const
    {
        return 4.0 / 3.0 * PI * a_eff[a] * a_eff[a] * a_eff[a];
    }

    double getMaterialDensity(uint a) const
    {
        return mass[a] / getVolume(a);
    }

    double getMaterialDensity() const
    {
        // Get local min and max grain sizes
        double a_min = getSizeMin();
        double a_max = getSizeMax();

        // Get local size parameter for size distribution
        double size_param = getSizeParam();

        // Get integration over the dust size distribution
        double * rel_weight = getRelWeight(a_min, a_max, size_param);

        for(uint a = 0; a < nr_of_dust_species; a++)
            rel_weight[a] *= mass[a] / getVolume(a);

        double res = CMathFunctions::integ_dust_size(a_eff, rel_weight, nr_of_dust_species, a_min, a_max);

        delete[] rel_weight;

        return res;
    }

    double getFHighJ() const
    {
        return f_highJ;
    }

    double getFcorr() const
    {
        return f_cor;
    }

    bool isAligned() const
    {
        return is_align;
    }

    void setIsAligned(bool val)
    {
        is_align = val;
    }

    void setIsMixture(bool val)
    {
        is_mixture = val;
    }

    double getSublimationTemperature()
    {
        return sub_temp;
    }

    void getQtrq(uint i, uint j, double & x, double & y)
    {
        x = Qtrq[i].getX(j);
        y = Qtrq[i].getY(j);
    }

    void getHG_g_factor(uint i, uint j, double & x, double & y)
    {
        x = HG_g_factor[i].getX(j);
        y = HG_g_factor[i].getY(j);
    }

    bool getCalorimetryLoaded()
    {
        return calorimetry_loaded;
    }

    void setCalorimetryLoaded(bool val)
    {
        calorimetry_loaded = val;
    }

    void setIDs(uint i_comp, uint nr_components, uint i_mix, uint nr_mixtures)
    {
        i_component = i_comp;
        nr_of_components = nr_components;

        i_mixture = i_mix;
        nr_of_mixtures = nr_mixtures;
    }

    void printIDs()
    {
        cout << "-> Dust mixture " << i_mixture + 1 << "/" << nr_of_mixtures << ", ";
        if(!is_mixture)
            cout << "component " << i_component + 1 << "/" << nr_of_components << " ";
    }

    bool getScatLoaded()
    {
        return scat_loaded;
    }

    void setScatLoaded(bool val)
    {
        scat_loaded = val;
    }

    void setSublimate(bool val)
    {
        sublimate = val;
    }

    uint getComponentId()
    {
        return i_component;
    }

    uint getNrOfComponents()
    {
        return nr_of_components;
    }

    void SetNrOfScatTheta(uint ** nr_of_scat_theta_tmp)
    {
        if(nr_of_scat_theta != 0)
        {
            for(uint a = 0; a < nr_of_dust_species; a++)
                delete[] nr_of_scat_theta[a];
            delete[] nr_of_scat_theta;
        }
        nr_of_scat_theta = nr_of_scat_theta_tmp;
    }

    void SetScatTheta(double *** scat_theta_tmp)
    {
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
        scat_theta = scat_theta_tmp;
    }

    double *getScatTheta(uint a, uint w)
    {
        return scat_theta[a][w];
    }

    double getScatTheta(uint a, uint w, uint sth) const
    {
        return scat_theta[a][w][sth];
    }

    uint getScatThetaID(double theta, uint a, uint w) const
    {
        // Returns the index of the value in scat_theta[a][w] closest to theta
        uint sth;

        // If no refinement was done, don't use binary
        // search to save (lots of) time
        if( (2*NANG-1) == nr_of_scat_theta[a][w])
            sth = uint(theta/PI * 2*(NANG-1) + 0.5);
        else
        {
            if(theta == scat_theta[a][w][0])
                return 0;

            sth = CMathFunctions::biListIndexSearch(theta,scat_theta[a][w],nr_of_scat_theta[a][w]);
            if(sth != MAX_UINT && (scat_theta[a][w][sth+1] + scat_theta[a][w][sth]) / 2. <= theta )
                sth++;
        }

        return sth;
    }

    void preCalcEffProperties(parameters & param);

    void henyeygreen(photon_package * pp, uint a, CRandomGenerator * rand_gen);
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
                                            dlist & wavelength_list_full) const;

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

    bool writeComponent(string path_data, string path_plot);
    bool calcSizeDistribution(dlist values, double * mass);
    bool add(double ** size_fraction, CDustComponent * comp, uint ** nr_of_scat_theta_tmp, double *** scat_theta_tmp);

    uint getInteractingDust(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen, uint cross_section = CROSS_ABS) const;

    void calcPACrossSections(uint a, uint w, cross_sections & cs, double theta) const;
    void calcNONPACrossSections(uint a, uint w, cross_sections & cs, double theta) const;
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

    spline *HG_g_factor, *Qtrq;
    spline * tab_planck;
    spline * tab_em;
    spline * tab_em_inv;
    spline tab_temp, tab_em_eff;

    prob_list * avg_planck_frac;
    prob_list *dust_prob, *sca_prob, *abs_prob;

    Matrix2D ***** sca_mat;
    double **Qext1, **Qext2, **Qabs1, **Qabs2, **Qsca1, **Qsca2, **Qcirc, **HGg;
    double ** enthalpy;
    double *a_eff, *a_eff_1_5, *a_eff_3_5, *a_eff_2;
    double * calorimetry_temperatures;
    double * mass;
    double *tCext1, *tCext2, *tCabs1, *tCabs2, *tCsca1, *tCsca2, *tCcirc, *tHGg;
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
};

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
    }

    ~CDustMixture(void)
    {
        if(single_component != 0)
            delete[] single_component;
        if(mixed_component != 0)
            delete[] mixed_component;
    }

    uint getMixtureID(CGridBasic * grid, const cell_basic & cell) const
    {
        uint dust_choice = grid->getDustChoiceID(cell);
        return dust_choices_to_index[dust_choice];
    }

    uint getMixtureID(CGridBasic * grid, const photon_package & pp) const
    {
        uint dust_choice = grid->getDustChoiceID(pp);
        return dust_choices_to_index[dust_choice];
    }

    double getAvgMass(uint i_mixture)
    {
        return mixed_component[i_mixture].getAvgMass();
    }

    bool writeComponent(string path_data, string path_plot)
    {
        if(mixed_component != 0)
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                if(!mixed_component[i_mixture].writeComponent(path_data, path_plot))
                    return false;
        return true;
    }

    void preCalcEffProperties(parameters & param)
    {
        if(mixed_component != 0)
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                mixed_component[i_mixture].preCalcEffProperties(param);
    }

    void preCalcRelWeight()
    {
        if(mixed_component != 0)
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                mixed_component[i_mixture].preCalcRelWeight();
    }

    string getPhaseFunctionStr(uint i_mixture)
    {
        string str_res = "\nERROR: Phase function is undefined!\n";
        if(mixed_component != 0)
        {
            switch(mixed_component[i_mixture].getPhaseFunctionID())
            {
                case PH_ISO:
                    str_res = "Isotropic scattering";
                    break;

                case PH_HG:
                    str_res = "Henyey and Greenstein";
                    break;

                case PH_MIE:
                    str_res = "Mie scattering";
                    break;
            }
        }
        return str_res;
    }

    void killSingleComponents()
    {
        // Show information
        cout << "                \r";
        cout << "-> Dust cleanup started            \r";

        // If single component was used -> delete it
        if(single_component != 0)
        {
            // Delete the single component
            delete[] single_component;
            single_component = 0;
        }
    }

    uint getNrOfDustSpecies(uint i_mixture)
    {
        return mixed_component[i_mixture].getNrOfDustSpecies();
    }

    uint getNrOfStochasticSizes(uint i_mixture)
    {
        return mixed_component[i_mixture].getNrOfStochasticSizes();
    }

    void calcTemperature(CGridBasic * grid, cell_basic * cell, bool use_energy_density)
    {
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, *cell);
                mixed_component[i_mixture].calcTemperature(grid, cell, 0, use_energy_density);
            }
            else
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    mixed_component[i_mixture].calcTemperature(grid, cell, i_mixture, use_energy_density);
        }
    }

    void calcStochasticHeatingPropabilities(CGridBasic * grid,
                                            cell_basic * cell,
                                            dlist & wavelength_list_full) const
    {
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, *cell);
                mixed_component[i_mixture].calcStochasticHeatingPropabilities(
                    grid, cell, 0, wavelength_list_full);
            }
            else
            {
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    mixed_component[i_mixture].calcStochasticHeatingPropabilities(
                        grid, cell, i_mixture, wavelength_list_full);
            }
        }
    }

    double getForegroundExtinction(double wavelength)
    {
        if(extinction_magnitude == 0)
            return 1;

        // Get extinction optical depth at extinction_magnitude_wavelength
        double extinction_tau = (extinction_magnitude / 1.086);

        if(extinction_dust_mixture == MAX_UINT)
        {
            // Init extinction curve
            spline extinction_curve(100);

            // Extinction curve for ISM grains (silicate and graphite, 5nm to 250nm)
            extinction_curve.setValue(0, 5.000000074505806e-08, 4.216777503871183e-016);
            extinction_curve.setValue(1, 5.564875528216362e-08, 4.621585931369242e-016);
            extinction_curve.setValue(2, 6.193568184971808e-08, 5.721308394426405e-016);
            extinction_curve.setValue(3, 6.893286854028702e-08, 7.124210419085484e-016);
            extinction_curve.setValue(4, 7.672057300806045e-08, 7.674365729865441e-016);
            extinction_curve.setValue(5, 8.538808673620224e-08, 7.177044099170471e-016);
            extinction_curve.setValue(6, 9.503481537103653e-08, 6.998261913784690e-016);
            extinction_curve.setValue(7, 1.05771400034428e-07, 6.237277455455178e-016);
            extinction_curve.setValue(8, 1.1772090196609499e-07, 5.435211899183833e-016);
            extinction_curve.setValue(9, 1.3102050125598898e-07, 4.750205237431847e-016);
            extinction_curve.setValue(10, 1.45822495222092e-07, 4.006340109059733e-016);
            extinction_curve.setValue(11, 1.62296801805496e-07, 3.590026327160000e-016);
            extinction_curve.setValue(12, 1.8063229322433499e-07, 3.539562486442753e-016);
            extinction_curve.setValue(13, 2.01039299368858e-07, 3.874130706789749e-016);
            extinction_curve.setValue(14, 2.2375169396400497e-07, 3.996385089649330e-016);
            extinction_curve.setValue(15, 2.49030098319054e-07, 3.285992056374746e-016);
            extinction_curve.setValue(16, 2.77164310216904e-07, 2.829757611670918e-016);
            extinction_curve.setValue(17, 3.0847701430320697e-07, 2.536181550806294e-016);
            extinction_curve.setValue(18, 3.43327194452286e-07, 2.287595115539490e-016);
            extinction_curve.setValue(19, 3.82114589214325e-07, 2.065593074020058e-016);
            extinction_curve.setValue(20, 4.25284087657928e-07, 1.845888079974461e-016);
            extinction_curve.setValue(21, 4.7333058714866596e-07, 1.634203032101485e-016);
            extinction_curve.setValue(22, 5.26805222034454e-07, 1.431729732227695e-016);
            extinction_curve.setValue(23, 5.86320996284485e-07, 1.250360258234746e-016);
            extinction_curve.setValue(24, 6.525607109069819e-07, 1.078952024050406e-016);
            extinction_curve.setValue(25, 7.2628378868103e-07, 9.266067295651665e-017);
            extinction_curve.setValue(26, 8.08335781097412e-07, 7.881166441657414e-017);
            extinction_curve.setValue(27, 8.99657726287842e-07, 6.579841269777561e-017);
            extinction_curve.setValue(28, 1.00129699707031e-06, 5.458088737076844e-017);
            extinction_curve.setValue(29, 1.1144180297851599e-06, 4.498859327641566e-017);
            extinction_curve.setValue(30, 1.2403199672698998e-06, 3.660451227198234e-017);
            extinction_curve.setValue(31, 1.3804450035095198e-06, 2.970563074866552e-017);
            extinction_curve.setValue(32, 1.5364010334014899e-06, 2.401290428990413e-017);
            extinction_curve.setValue(33, 1.70997595787048e-06, 1.921721053770342e-017);
            extinction_curve.setValue(34, 1.90316104888916e-06, 1.523635586331645e-017);
            extinction_curve.setValue(35, 2.11817002296448e-06, 1.204889560936278e-017);
            extinction_curve.setValue(36, 2.35747098922729e-06, 9.562793524986687e-018);
            extinction_curve.setValue(37, 2.6238059997558596e-06, 7.653670597381650e-018);
            extinction_curve.setValue(38, 2.92023110389709e-06, 6.185056501137134e-018);
            extinction_curve.setValue(39, 3.25014495849609e-06, 5.054911063069209e-018);
            extinction_curve.setValue(40, 3.61733102798462e-06, 4.168021658769960e-018);
            extinction_curve.setValue(41, 4.02599906921387e-06, 3.468548840603535e-018);
            extinction_curve.setValue(42, 4.4808359146118195e-06, 2.910007332891138e-018);
            extinction_curve.setValue(43, 4.98706007003784e-06, 2.474287670482146e-018);
            extinction_curve.setValue(44, 5.5504732131958e-06, 2.148568553881987e-018);
            extinction_curve.setValue(45, 6.17753791809082e-06, 1.926424030502043e-018);
            extinction_curve.setValue(46, 6.875445842742919e-06, 1.737387824729840e-018);
            extinction_curve.setValue(47, 7.652201175689699e-06, 2.102465435646191e-018);
            extinction_curve.setValue(48, 8.51670932769775e-06, 6.233472209570475e-018);
            extinction_curve.setValue(49, 9.47888565063477e-06, 1.116189825799359e-017);
            extinction_curve.setValue(50, 1.05497598648071e-05, 8.459390008941990e-018);
            extinction_curve.setValue(51, 1.1741620063781699e-05, 5.243048743901026e-018);
            extinction_curve.setValue(52, 1.3068140029907199e-05, 3.222513686771834e-018);
            extinction_curve.setValue(53, 1.45445098876953e-05, 2.593459581332661e-018);
            extinction_curve.setValue(54, 1.6187679290771496e-05, 3.364742430016652e-018);
            extinction_curve.setValue(55, 1.80164794921875e-05, 4.028119360997947e-018);
            extinction_curve.setValue(56, 2.00519008636475e-05, 3.589028317207343e-018);
            extinction_curve.setValue(57, 2.2317260742187498e-05, 2.919286851797911e-018);
            extinction_curve.setValue(58, 2.48385601043701e-05, 2.424045340461135e-018);
            extinction_curve.setValue(59, 2.76446990966797e-05, 2.034588235351009e-018);
            extinction_curve.setValue(60, 3.0767860412597695e-05, 1.725524296987028e-018);
            extinction_curve.setValue(61, 3.42438583374023e-05, 1.464357218017165e-018);
            extinction_curve.setValue(62, 3.8112560272216795e-05, 1.239701787887817e-018);
            extinction_curve.setValue(63, 4.2418331146240196e-05, 1.028703213069295e-018);
            extinction_curve.setValue(64, 4.72105484008789e-05, 8.412356002324342e-019);
            extinction_curve.setValue(65, 5.25441703796387e-05, 6.810797851620609e-019);
            extinction_curve.setValue(66, 5.84803504943848e-05, 5.471715483013374e-019);
            extinction_curve.setValue(67, 6.50871810913086e-05, 4.374078863404243e-019);
            extinction_curve.setValue(68, 7.24404067993164e-05, 3.489064367360922e-019);
            extinction_curve.setValue(69, 8.06243667602539e-05, 2.781946435780552e-019);
            extinction_curve.setValue(70, 8.973293304443359e-05, 2.221216725155223e-019);
            extinction_curve.setValue(71, 9.987050628662108e-05, 1.777522731448490e-019);
            extinction_curve.setValue(72, 0.00011115339660644499, 1.426858510061328e-019);
            extinction_curve.setValue(73, 0.000123710998535156, 1.155888746296647e-019);
            extinction_curve.setValue(74, 0.000137687194824219, 9.449756554245537e-020);
            extinction_curve.setValue(75, 0.000153242401123047, 7.584280096005196e-020);
            extinction_curve.setValue(76, 0.000170554992675781, 6.018818168528717e-020);
            extinction_curve.setValue(77, 0.000189823501586914, 4.803205495665578e-020);
            extinction_curve.setValue(78, 0.00021126879882812498, 3.850335435269971e-020);
            extinction_curve.setValue(79, 0.000235136993408203, 3.094874612999559e-020);
            extinction_curve.setValue(80, 0.00026170159912109397, 2.491039874419238e-020);
            extinction_curve.setValue(81, 0.00029126739501953103, 2.006959239669001e-020);
            extinction_curve.setValue(82, 0.00032417330932617196, 1.618343487101766e-020);
            extinction_curve.setValue(83, 0.00036079681396484396, 1.305601603235169e-020);
            extinction_curve.setValue(84, 0.000401557891845703, 1.053083960958686e-020);
            extinction_curve.setValue(85, 0.00044692401123046894, 8.500255800817664e-021);
            extinction_curve.setValue(86, 0.000497415313720703, 6.863925335919417e-021);
            extinction_curve.setValue(87, 0.000553610778808594, 5.541157311593072e-021);
            extinction_curve.setValue(88, 0.000616155029296875, 4.473916606814938e-021);
            extinction_curve.setValue(89, 0.000685765197753906, 3.612131329544281e-021);
            extinction_curve.setValue(90, 0.0007632396240234379, 2.917011636593121e-021);
            extinction_curve.setValue(91, 0.0008494666748046879, 2.355426278709940e-021);
            extinction_curve.setValue(92, 0.000945435302734375, 1.902069506823009e-021);
            extinction_curve.setValue(93, 0.0010522459716796899, 1.538983867990721e-021);
            extinction_curve.setValue(94, 0.0011711240234375, 1.232992068728240e-021);
            extinction_curve.setValue(95, 0.0013034310302734399, 9.975786070486010e-022);
            extinction_curve.setValue(96, 0.0014506870117187499, 8.016047070202548e-022);
            extinction_curve.setValue(97, 0.0016145780029296899, 6.504697113317689e-022);
            extinction_curve.setValue(98, 0.0017969849853515598, 5.270413281701706e-022);
            extinction_curve.setValue(99, 0.002, 4.197725271491968e-022);

            // Init spline
            extinction_curve.createSpline();

            // Scaling factor
            double scaling_factor =
                extinction_tau / extinction_curve.getValue(extinction_magnitude_wavelength);

            return exp(-scaling_factor * extinction_curve.getValue(wavelength));
        }
        else
        {
            // Scaling factor
            double scaling_factor = extinction_tau / mixed_component[extinction_dust_mixture].getCextMean(
                                                         extinction_magnitude_wavelength);

            return exp(-scaling_factor * mixed_component[extinction_dust_mixture].getCextMean(wavelength));
        }
    }

    double getCextMean(CGridBasic * grid, const photon_package & pp) const
    {
        double sum = 0;
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, pp);
                sum = mixed_component[i_mixture].getCextMean(grid, pp);
            }
            else
            {
                if(grid->getCextMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID()) != MAX_DOUBLE)
                    sum = grid->getCextMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID());
                else
                {
                    double dens = 0;
                    for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    {
                        double i_dens = getNumberDensity(grid, pp, i_mixture);
                        dens += i_dens;
                        sum += mixed_component[i_mixture].getCextMean(grid, pp) * i_dens;
                    }
                    if(dens != 0)
                        sum /= dens;
                }
            }
        }
        return sum;
    }

    double getCabsMean(CGridBasic * grid, const photon_package & pp) const
    {
        double sum = 0;
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, pp);
                sum = mixed_component[i_mixture].getCabsMean(grid, pp);
            }
            else
            {
                if(grid->getCabsMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID()) != MAX_DOUBLE)
                    sum = grid->getCabsMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID());
                else
                {
                    double dens = 0;
                    for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    {
                        double i_dens = getNumberDensity(grid, pp, i_mixture);
                        dens += i_dens;
                        sum += mixed_component[i_mixture].getCabsMean(grid, pp) * i_dens;
                    }
                    if(dens != 0)
                        sum /= dens;
                }
            }
        }
        return sum;
    }

    double getCscaMean(CGridBasic * grid, const photon_package & pp) const
    {
        double sum = 0;
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, pp);
                sum = mixed_component[i_mixture].getCscaMean(grid, pp);
            }
            else
            {
                if(grid->getCscaMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID()) != MAX_DOUBLE)
                    sum = grid->getCscaMeanTab(pp.getPositionCell()->getUniqueID(),pp.getDustWavelengthID());
                else
                {
                    double dens = 0;
                    for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    {
                        double i_dens = getNumberDensity(grid, pp, i_mixture);
                        dens += i_dens;
                        sum += mixed_component[i_mixture].getCscaMean(grid, pp) * i_dens;
                    }
                    if(dens != 0)
                        sum /= dens;
                }
            }
        }
        return sum;
    }

    bool adjustTempAndWavelengthBW(CGridBasic * grid, photon_package * pp, bool use_energy_density, CRandomGenerator * rand_gen)
    {
        if(mixed_component != 0)
        {
            uint i_mixture = getEmittingMixture(grid, pp, rand_gen);
            return mixed_component[i_mixture].adjustTempAndWavelengthBW(
                grid, pp, i_mixture, use_energy_density, rand_gen);
        }
        return false;
    }

    void calcAlignedRadii(CGridBasic * grid, cell_basic * cell)
    {
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, *cell);
                mixed_component[i_mixture].calcAlignedRadii(grid, cell, 0);
            }
            else
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    mixed_component[i_mixture].calcAlignedRadii(grid, cell, i_mixture);
        }
    }

    void addToWavelengthGrid(double wavelength)
    {
        // Add wavelength to list
        wavelength_list.push_back(wavelength);
    }

    void addToWavelengthGrid(double lam_min, double lam_max, double nr_of_wavelength, bool add_offset = false)
    {
        dlist tmp_wavelength_list(nr_of_wavelength);

        // Add the amount of added wavelength to the offset to skip them for
        // writeComponent
        if(add_offset && wavelength_list.size() == 0)
            wavelength_offset += nr_of_wavelength;

        // Add new wavelengths first
        CMathFunctions::LogList(lam_min, lam_max, tmp_wavelength_list, 10);

        // Update wavelengths
        for(uint w = 0; w < nr_of_wavelength; w++)
            wavelength_list.push_back(tmp_wavelength_list[w]);
    }

    void finalizeWavelengthList()
    {
        sort(wavelength_list.begin(), wavelength_list.end());
        wavelength_list.erase(unique(wavelength_list.begin(), wavelength_list.end()), wavelength_list.end());

        nr_of_wavelength = wavelength_list.size();
    }

    void setScatteringToRay(bool val)
    {
        scattering_to_raytracing = val;
    }

    bool getScatteringToRay()
    {
        return scattering_to_raytracing;
    }

    const dlist & getWavelengthList() const
    {
        return wavelength_list;
    }

    double getWavelength(uint wID)
    {
        return wavelength_list[wID];
    }

    double getWavelength(photon_package * pp)
    {
        return wavelength_list[pp->getDustWavelengthID()];
    }

    uint getNrOfWavelength()
    {
        return nr_of_wavelength;
    }

    uint getWavelengthID(double wavelength)
    {
        dlist::iterator it = find(wavelength_list.begin(), wavelength_list.end(), wavelength);
        if(it != wavelength_list.end())
            return distance(wavelength_list.begin(), it);

        cout << "\nHINT: Wavelength not found!" << endl;
        return 0;
    }

    double getSizeMin(uint i_mixture)
    {
        return mixed_component[i_mixture].getSizeMin();
    }

    double getSizeMin(CGridBasic * grid, const cell_basic & cell) const
    {
        double min_a = 1e200;
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, cell);
                min_a = mixed_component[i_mixture].getSizeMin();
            }
            else
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                {
                    double tmp_a = mixed_component[i_mixture].getSizeMin();
                    if(tmp_a < min_a)
                        min_a = tmp_a;
                }
        }
        return min_a;
    }

    double getEffectiveRadius(uint i_mixture, uint a)
    {
        return mixed_component[i_mixture].getEffectiveRadius(a);
    }

    double * getEffectiveRadii(uint i_mixture)
    {
        return mixed_component[i_mixture].getEffectiveRadii();
    }

    uint getNrOfCalorimetryTemperatures(uint i_mixture)
    {
        return mixed_component[i_mixture].getNrOfCalorimetryTemperatures();
    }

    // double getMinDustTemp()
    // {
    //     double min_temp = 1e200;
    //     if(mixed_component != 0)
    //         for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
    //         {
    //             double temp = mixed_component[i_mixture].getMinDustTemp();
    //             if(temp < min_temp)
    //                 min_temp = temp;
    //         }
    //     return min_temp;
    // }

    double getMaxDustTemp()
    {
        double max_temp = 0;
        if(mixed_component != 0)
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            {
                double temp = mixed_component[i_mixture].getMaxDustTemp();
                if(temp > max_temp)
                    max_temp = temp;
            }
        return max_temp;
    }

    double getMinAlignedRadius()
    {
        double min_a_alig = 1e200;
        if(mixed_component != 0)
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            {
                double a_alig = mixed_component[i_mixture].getMinAlignedRadius();
                if(a_alig < min_a_alig)
                    min_a_alig = a_alig;
            }
        return min_a_alig;
    }

    double getMaxAlignedRadius()
    {
        double max_a_alig = 0;
        if(mixed_component != 0)
            for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
            {
                double a_alig = mixed_component[i_mixture].getMaxAlignedRadius();
                if(a_alig > max_a_alig)
                    max_a_alig = a_alig;
            }
        return max_a_alig;
    }

    // Dust number density functions
    double getNumberDensity(CGridBasic * grid, const cell_basic & cell) const
    {
        double sum = 0;
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, cell);
                sum = mixed_component[i_mixture].getNumberDensity(grid, cell);
            }
            else
            {
                if(grid->getNumberDensityTab(cell.getUniqueID()) != MAX_DOUBLE)
                    sum = grid->getNumberDensityTab(cell.getUniqueID());
                else
                    for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                        sum += mixed_component[i_mixture].getNumberDensity(grid, cell, i_mixture);
            }
        }
        return sum;
    }

    double getNumberDensity(CGridBasic * grid, const photon_package & pp) const
    {
        return getNumberDensity(grid, *pp.getPositionCell());
    }

    double getNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_mixture) const
    {
        double sum = 0;
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                if(i_mixture == getMixtureID(grid, cell))
                    sum = mixed_component[i_mixture].getNumberDensity(grid, cell);
            }
            else
                sum = mixed_component[i_mixture].getNumberDensity(grid, cell, i_mixture);
        }
        return sum;
    }

    double getNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_mixture) const
    {
        return getNumberDensity(grid, *pp.getPositionCell(), i_mixture);
    }

    // Dust mass density functions
    double getMassDensity(CGridBasic * grid, const cell_basic & cell) const
    {
        double sum = 0;
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, cell);
                sum = mixed_component[i_mixture].getMassDensity(grid, cell);
            }
            else
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    sum += mixed_component[i_mixture].getMassDensity(grid, cell, i_mixture);
        }
        return sum;
    }

    double getMassDensity(CGridBasic * grid, const photon_package & pp) const
    {
        return getMassDensity(grid, *pp.getPositionCell());
    }

    double getMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_mixture) const
    {
        double sum = 0;
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                if(i_mixture == getMixtureID(grid, cell))
                    sum = mixed_component[i_mixture].getMassDensity(grid, cell);
            }
            else
                sum = mixed_component[i_mixture].getMassDensity(grid, cell, i_mixture);
        }
        return sum;
    }

    double getMassDensity(CGridBasic * grid, const photon_package & pp, uint i_mixture) const
    {
        return getMassDensity(grid, *pp.getPositionCell(), i_mixture);
    }

    double getRelativeDustNumberDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const
    {
        double dens = getNumberDensity(grid, cell);
        if(dens != 0)
            return getNumberDensity(grid, cell, i_density) / dens;
        else
            return 0;
    }

    double getRelativeDustNumberDensity(CGridBasic * grid, const photon_package & pp, uint i_density) const
    {
        return getRelativeDustNumberDensity(grid, *pp.getPositionCell(), i_density);
    }

    double getRelativeDustMassDensity(CGridBasic * grid, const cell_basic & cell, uint i_density) const
    {
        double dens = getMassDensity(grid, cell);
        if(dens != 0)
            return getMassDensity(grid, cell, i_density) / dens;
        else
            return 0;
    }

    double getRelativeDustMassDensity(CGridBasic * grid, const photon_package & pp, uint i_density) const
    {
        return getRelativeDustMassDensity(grid, *pp.getPositionCell(), i_density);
    }

    void calcEmissivityHz(CGridBasic * grid, const photon_package & pp, StokesVector * dust_emissivity)
    {
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, pp);
                mixed_component[i_mixture].calcEmissivityHz(grid, pp, 0, dust_emissivity);
            }
            else
            {
                *dust_emissivity = 0;
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                {
                    mixed_component[i_mixture].calcEmissivityHz(grid, pp, i_mixture, dust_emissivity);
                }
            }
        }
    }

    double calcEmissivity(CGridBasic * grid, const photon_package & pp)
    {
        // Init variables
        double pl_abs = 0;

        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, pp);
                pl_abs = mixed_component[i_mixture].calcEmissivity(grid, pp, 0);
            }
            else
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    pl_abs += mixed_component[i_mixture].calcEmissivity(grid, pp, i_mixture);
        }

        return pl_abs;
    }

    void calcEmissivityExt(CGridBasic * grid, const photon_package & pp, Matrix2D * dust_ext_matrix) const
    {
        // Init variables
        double Cext = 0, Cpol = 0, Ccirc = 0;

        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, pp);
                double dens_dust = mixed_component[i_mixture].getNumberDensity(grid, pp);
                mixed_component[i_mixture].calcExtCrossSections(grid, pp, 0, &Cext, &Cpol, &Ccirc);
                Cext *= -dens_dust;
                Cpol *= -dens_dust;
                Ccirc *= -dens_dust;
            }
            else
            {
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                {
                    double tmpCext, tmpCpol, tmpCcirc;
                    double dens_dust = mixed_component[i_mixture].getNumberDensity(grid, pp, i_mixture);
                    mixed_component[i_mixture].calcExtCrossSections(
                        grid, pp, i_mixture, &tmpCext, &tmpCpol, &tmpCcirc);
                    Cext -= tmpCext * dens_dust;
                    Cpol -= tmpCpol * dens_dust;
                    Ccirc -= tmpCcirc * dens_dust;
                }
            }
        }

        // Create extinction matrix for the dust grains
        double phi = grid->getPhiMag(pp);
        double sin_2ph = sin(2.0 * phi);
        double cos_2ph = cos(2.0 * phi);

        // Reset matrix
        dust_ext_matrix->resize(4, 4);

        dust_ext_matrix->setValue(0, 0, Cext);
        dust_ext_matrix->setValue(1, 1, Cext);
        dust_ext_matrix->setValue(2, 2, Cext);
        dust_ext_matrix->setValue(3, 3, Cext);

        dust_ext_matrix->setValue(0, 1, Cpol * cos_2ph);
        dust_ext_matrix->setValue(0, 2, Cpol * sin_2ph);

        dust_ext_matrix->setValue(1, 0, Cpol * cos_2ph);
        dust_ext_matrix->setValue(2, 0, Cpol * sin_2ph);

        dust_ext_matrix->setValue(1, 3, Ccirc * sin_2ph);
        dust_ext_matrix->setValue(2, 3, -Ccirc * cos_2ph);

        dust_ext_matrix->setValue(3, 1, -Ccirc * sin_2ph);
        dust_ext_matrix->setValue(3, 2, Ccirc * cos_2ph);
    }

    void calcEmissivityEmi(CGridBasic * grid,
                           const photon_package & pp,
                           uint i_offset,
                           uint emission_component,
                           StokesVector * dust_emissivity) const
    {
        // Init variables
        double energy = 0;
        double phi = grid->getPhiMag(pp);
        Vector3D en_dir;

        // Check if radiation field is available and scattering should be included
        if(scattering_to_raytracing &&
           (emission_component == DUST_EMI_FULL || emission_component == DUST_EMI_SCAT))
        {
            // Get wavelength of photon package
            uint w = pp.getDustWavelengthID();

            // Get radiation field and calculate angle to the photon package direction
            if(grid->isRadiationFieldAvailable())
                grid->getRadiationFieldInterp(pp, wavelength_list[w], &energy, &en_dir);
            else if(i_offset != MAX_UINT)
            {
                // Use the rad field as stokes vector if all necessary things were already calculated
                *dust_emissivity += grid->getStokesFromRadiationField(pp, i_offset);
            }
            else
                grid->getRadiationField(pp, w, &energy, &en_dir);
        }

        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, pp);
                *dust_emissivity += mixed_component[i_mixture].calcEmissivityEmi(
                    grid, pp, 0, emission_component, phi, energy, en_dir);
            }
            else
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    *dust_emissivity += mixed_component[i_mixture].calcEmissivityEmi(
                        grid, pp, i_mixture, emission_component, phi, energy, en_dir);
        }
    }

    StokesVector getRadFieldScatteredFraction(CGridBasic * grid,
                                              const photon_package & pp,
                                              const Vector3D & en_dir,
                                              double energy) const
    {
        StokesVector tmp_stokes;

        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, pp);
                tmp_stokes =
                    mixed_component[i_mixture].getRadFieldScatteredFraction(grid, pp, 0, en_dir, energy);
            }
            else
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    tmp_stokes += mixed_component[i_mixture].getRadFieldScatteredFraction(
                        grid, pp, i_mixture, en_dir, energy);
        }
        return tmp_stokes;
    }

    uint getNrOfMixtures() const
    {
        return dust_choices.size();
    }

    void convertTempInQB(CGridBasic * grid, cell_basic * cell, double min_gas_density, bool use_gas_temp)
    {
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, *cell);
                mixed_component[i_mixture].convertTempInQB(grid, cell, 0, min_gas_density, use_gas_temp);
            }
            else
                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    mixed_component[i_mixture].convertTempInQB(
                        grid, cell, i_mixture, min_gas_density, use_gas_temp);
        }
    }

    uint getScatteringMixture(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen) const
    {
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, *pp);
                return i_mixture;
            }
            else
            {
                double rnd = rand_gen->getRND();
                double dens = 0;
                double pb[getNrOfMixtures()];

                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                {
                    double i_dens = getNumberDensity(grid, *pp, i_mixture);
                    pb[i_mixture] = i_dens * mixed_component[i_mixture].getCscaMean(grid, *pp);
                    if(i_mixture > 0)
                        pb[i_mixture] += pb[i_mixture-1];
                    dens += i_dens;
                }

                double cscamean_tmp = getCscaMean(grid, *pp);

                if(dens!=0 && cscamean_tmp != 0)
                {
                    for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    {
                        pb[i_mixture] /= dens * cscamean_tmp;
                        if(pb[i_mixture] > (1-rnd))
                            return i_mixture;
                    }
                }
            }
        }
        return 0;
    }

    uint getEmittingMixture(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen) const
    {
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, *pp);
                return i_mixture;
            }
            else
            {
                double rnd = rand_gen->getRND();
                double dens = 0;
                double pb[getNrOfMixtures()];

                for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                {
                    double i_dens = getNumberDensity(grid, *pp, i_mixture);
                    pb[i_mixture] = i_dens * mixed_component[i_mixture].getCabsMean(grid, *pp);
                    if(i_mixture > 0)
                        pb[i_mixture] += pb[i_mixture-1];
                    dens += i_dens;
                }

                double cabsmean_tmp = getCabsMean(grid, *pp);

                if(dens!=0 && cabsmean_tmp != 0)
                {
                    for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                    {
                        pb[i_mixture] /= dens * cabsmean_tmp;
                        if(pb[i_mixture] > (1-rnd))
                            return i_mixture;
                    }
                }
            }
        }
        return 0;
    }

    void scatter(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen, bool adjust_stokes = false)
    {
        if(mixed_component != 0)
        {
            uint i_mixture = getScatteringMixture(grid, pp, rand_gen);
            mixed_component[i_mixture].scatter(grid, pp, rand_gen);

            // Reduce the Stokes vector by the mean albedo of the particles
            if(adjust_stokes)
            {
                if(getCextMean(grid, *pp) > 0.0)
                    *pp->getStokesVector() *= getCscaMean(grid, *pp) / getCextMean(grid, *pp);
                else
                {
                    cout << "\nHINT: Mean cross section for extinction is zero or negative!" << endl;
                    pp->getStokesVector()->clear();
                }
            }
        }
    }

    void setGridRequirements(CGridBasic * grid, parameters & param)
    {
        uint nr_of_mixtures = getNrOfMixtures();
        uint * nr_dust_temp_sizes = new uint[nr_of_mixtures];
        uint * nr_stochastic_temps = new uint[nr_of_mixtures];
        uint * nr_stochastic_sizes = new uint[nr_of_mixtures];
        if(mixed_component != 0)
        {
            for(uint i_mixture = 0; i_mixture < nr_of_mixtures; i_mixture++)
            {
                // Maximum amount of dust grain sizes for which the grid has to contain a
                // temperature
                nr_dust_temp_sizes[i_mixture] = getNrOfDustSpecies(i_mixture);
                // Temperatures for which the griad has to contain propabilities
                nr_stochastic_temps[i_mixture] = getNrOfCalorimetryTemperatures(i_mixture);
                // Maximum amount of dust grain sizes affected by stochastic heating
                nr_stochastic_sizes[i_mixture] = getNrOfStochasticSizes(i_mixture);
            }
        }
        grid->setDustInformation(
            nr_of_mixtures, nr_dust_temp_sizes, nr_stochastic_sizes, nr_stochastic_temps);
    }

    void setIDs(CDustComponent & component,
                uint i_comp,
                uint nr_of_components,
                uint i_mixture,
                uint nr_of_mixtures)
    {
        component.setIDs(i_comp, nr_of_components, i_mixture, nr_of_mixtures);
    }

    double *** getSizeFractions()
    {
        // Get number of dust size bins
        uint nr_of_dust_species = single_component[0].getNrOfDustSpecies();

        // Init pointer array
        double *** size_fraction = new double **[nr_of_components];
        for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
        {
            size_fraction[i_comp] = new double *[nr_of_dust_species];
            for(uint a = 0; a < nr_of_dust_species; a++)
                size_fraction[i_comp][a] = new double[2];
        }

        // Put the mass weights into relation to each other
        double ref_mass_weight = single_component[0].getMassWeight() / single_component[0].getFraction();
        for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
        {
            double mass_weight_ratio = ref_mass_weight * single_component[i_comp].getFraction() /
                                       single_component[i_comp].getMassWeight();

            // First column is the updated size distribution
            for(uint a = 0; a < nr_of_dust_species; a++)
                if(single_component[i_comp].sizeIndexUsed(a))
                {
                    size_fraction[i_comp][a][0] =
                        mass_weight_ratio * single_component[i_comp].getEffectiveRadius3_5(a);
                }
        }

        for(uint a = 0; a < nr_of_dust_species; a++)
        {
            // Calulcate the sum for each size bin
            double sum = 0;
            for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
                if(single_component[i_comp].sizeIndexUsed(a))
                    sum += size_fraction[i_comp][a][0];

            // Second column is the mixing ration between the dust components
            for(uint i_comp = 0; i_comp < nr_of_components; i_comp++)
                if(single_component[i_comp].sizeIndexUsed(a))
                    size_fraction[i_comp][a][1] = size_fraction[i_comp][a][0] / sum;
        }
        return size_fraction;
    }

    void getEscapePhoton(CGridBasic * grid,
                         photon_package * pp,
                         Vector3D obs_ex,
                         Vector3D dir_obs,
                         photon_package * pp_escape,
                         CRandomGenerator * rand_gen) const
    {
        if(mixed_component != 0)
        {
            uint i_mixture = getScatteringMixture(grid, pp, rand_gen);
            uint a = mixed_component[i_mixture].getInteractingDust(grid, pp, rand_gen, CROSS_SCA);
            mixed_component[i_mixture].getEscapePhoton(grid, pp, a, obs_ex, dir_obs, pp_escape);

            // Reduce the Stokes vector by the mean albedo of the particles
            if(getCextMean(grid, *pp) > 0.0)
                *pp_escape->getStokesVector() *= getCscaMean(grid, *pp) / getCextMean(grid, *pp);
            else
            {
                cout << "\nHINT: Mean cross section for extinction is zero or negative!" << endl;
                pp_escape->getStokesVector()->clear();
            }
        }
    }

    double getCellEmission(CGridBasic * grid, photon_package * pp, CRandomGenerator * rand_gen) const
    {
        if(mixed_component != 0)
        {
            uint i_mixture = getEmittingMixture(grid, pp, rand_gen);
            if(grid->useDustChoice())
                return mixed_component[i_mixture].getCellEmission(grid, *pp, 0);
            else
                return mixed_component[i_mixture].getCellEmission(grid, *pp, i_mixture);
        }
        return 0;
    }

    double getTotalCellEmission(CGridBasic * grid, const photon_package & pp) const
    {
        double sum = 0;
        if(mixed_component != 0)
        {
            if(grid->useDustChoice())
            {
                uint i_mixture = getMixtureID(grid, pp);
                sum = mixed_component[i_mixture].getCellEmission(grid, pp, 0);
            }
            else
            {
                if(grid->getTotalCellEmissionTab(pp.getPositionCell()->getUniqueID()) != MAX_DOUBLE)
                    sum = grid->getTotalCellEmissionTab(pp.getPositionCell()->getUniqueID());
                else
                    for(uint i_mixture = 0; i_mixture < getNrOfMixtures(); i_mixture++)
                        sum += mixed_component[i_mixture].getCellEmission(grid, pp, i_mixture);
            }
        }
        return sum;
    }

    double getPlanck(uint w, double temp)
    {
        return mixed_component[0].getPlanck(w, temp);
    }

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

    bool scattering_to_raytracing;

    double extinction_magnitude;
    double extinction_magnitude_wavelength;
    uint extinction_dust_mixture;

    uint nr_of_components;
    uint nr_of_wavelength;
    uint wavelength_offset;
    uint nr_of_dust_species;

    uilist dust_choices;
    uilist dust_choices_to_index;

    dlist wavelength_list;

    spline diff_y; // diff_y as a function of z
};

#endif
