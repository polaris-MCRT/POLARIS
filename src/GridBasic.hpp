/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CGRID_BASIC_H
#define CGRID_BASIC_H

#include "Parameters.hpp"
#include "Photon.hpp"
#include "Vector3D.hpp"
#include "MathFunctions.hpp"
#include "RandomGenerator.hpp"
#include "MathSpline.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"

// Additional Structure
struct VelFieldInterp
{
    spline vel_field;
    bool zero_vel_field;
    Vector3D start_pos;
};

// Additional Structures
class MagFieldInfo
{
public:
    MagFieldInfo()
    {
        cos_theta=0;
        sin_theta=0;
        cos_2_phi=0;
        sin_2_phi=0;
    }


    double cos_theta;
    double sin_theta;
    double cos_2_phi;
    double sin_2_phi;
    Vector3D mag_field;
};

class LineBroadening
{
public:
    LineBroadening()
    {
        gauss_a=0;
        voigt_a=0;
    }

    double gauss_a;
    double voigt_a;
};

class CGridBasic
{
public:
    CGridBasic()
    {
        basic_path = 0;
        buffer_size = 0;

        max_cells = 0;
        max_value = 0;
        max_data = 0;

        min_delta = 0;
        max_delta = 0;

        min_mach = 0;
        max_mach = 0;

        min_mag = 0;
        max_mag = 0;

        min_vel = 0;
        max_vel = 0;

        min_len = 0;
        max_len = 0;

        min_gas_temp = 0;
        max_gas_temp = 0;

        min_dust_temp = 0;
        max_dust_temp = 0;

        min_gas_dens = 0;
        max_gas_dens = 0;

        min_dust_dens = 0;
        max_dust_dens = 0;

        aalg_min = 0;
        aalg_max = 0;

        a_min_min = 0;
        a_min_max = 0;

        a_max_min = 0;
        a_max_max = 0;

        size_param_min = 0;
        size_param_max = 0;

        dust_id_min = 0;
        dust_id_max = 0;

        min_larm_limit = 0;
        max_larm_limit = 0;

        min_pres = 0;
        max_pres = 0;

        line_counter = 0;
        char_counter = 0;
        ru[0] = '|';
        ru[1] = '/';
        ru[2] = '-';
        ru[3] = '\\';

        conv_length_in_SI = 1;
        conv_dens_in_SI = 1;
        conv_Bfield_in_SI = 1;
        conv_Vfield_in_SI = 1;

        dust_is_mass_density = false;
        gas_is_mass_density = false;
        velocity_field_needed = false;
        spec_length_as_vector = false;

        nrOfPlotPoints = 1000;
        nrOfPlotVectors = 1000;
        maxPlotLines = 3;

        cell_list = 0;

        data_offset = 6;
        dataID = 0;
        data_len = 0;

        total_gas_mass = 0;
        mu = 0;

        nrOfDensRatios = 0;
        nrOfOpiateIDs = 0;

        nr_densities = 1;
        size_gd_list = 0;
        size_dd_list = 0;
        multi_temperature_entries = 0;
        stochastic_temperature_entries = 0;
        nr_mixtures = 0;

        nr_dust_temp_sizes = 0;
        nr_stochastic_sizes = 0;
        nr_stochastic_temps = 0;
        size_skip = 0;

        level_to_pos = 0;
        line_to_pos = 0;

        data_pos_tg = MAX_UINT;
        data_pos_mx = MAX_UINT;
        data_pos_my = MAX_UINT;
        data_pos_mz = MAX_UINT;
        data_pos_vx = MAX_UINT;
        data_pos_vy = MAX_UINT;
        data_pos_vz = MAX_UINT;
        data_pos_px = MAX_UINT;
        data_pos_py = MAX_UINT;
        data_pos_pz = MAX_UINT;
        
        data_pos_amin = MAX_UINT;
        data_pos_amax = MAX_UINT;
        data_pos_size_param = MAX_UINT;
        //data_pos_ra = MAX_UINT;
        data_pos_id = MAX_UINT; //dust id

        data_pos_vt = MAX_UINT;
        data_pos_pda = MAX_UINT;

        pos_GasSpecRatios = 0;
        pos_OpiateIDS = 0;

        rot_angle1 = 0;
        rot_angle2 = 0;

        nr_rad_field_comp = 1;

        plt_gas_dens = false;
        plt_mol_dens = false;
        plt_dust_dens = false;
        plt_gas_temp = false;
        plt_dust_temp = false;
        plt_mag = false;
        plt_vel = false;
        plt_rat = false;
        //plt_delta = false;
        plt_larm = false;
        //plt_mach = false;
        plt_dust_id = false;
        plt_amin = false;
        plt_amax = false;
        plt_size_param = false;
        plt_rad_field1 = false;
        plt_u_rad = false;
        plt_g_zero1 = false;
        plt_n_th = false;
        plt_T_e = false;
        plt_n_cr = false;
        plt_g_min = false;
        plt_g_max = false;
        plt_p = false;

        plt_avg_dir = false;
        plt_avg_th = false;

        total_volume = 0;
        cell_volume = 0;

        buffer_gas_dens = 0;
        buffer_mol_dens = 0;
        buffer_dust_dens = 0;
        buffer_gas_temp = 0;
        buffer_dust_temp = 0;
        buffer_rat = 0;
        //buffer_delta = 0;
        buffer_mag = 0;
        buffer_mag_x = 0;
        buffer_mag_y = 0;
        buffer_mag_z = 0;
        buffer_vel = 0;
        buffer_vel_x = 0;
        buffer_vel_y = 0;
        buffer_vel_z = 0;
        buffer_larm = 0;
        //buffer_mach = 0;
        buffer_dust_mixture = 0;
        buffer_dust_amin = 0;
        buffer_dust_amax = 0;
        buffer_dust_size_param = 0;
        buffer_rad_field = 0;
        buffer_g_zero1 = 0;
        buffer_u_rad = 0;
        buffer_n_th = 0;
        buffer_T_e = 0;
        buffer_n_cr = 0;
        buffer_g_min = 0;
        buffer_g_max = 0;
        buffer_p = 0;

        buffer_avg_dir = 0;
        buffer_avg_th = 0;

        turbulent_velocity = 0;

        wl_list.resize(WL_STEPS);
        CMathFunctions::LogList(WL_MIN, WL_MAX, wl_list, 10);

        CextMeanTab = 0;
        CabsMeanTab = 0;
        CscaMeanTab = 0;
        numberDensityTab = 0;
        totalCellEmissionTab = 0;
        max_wavelengths = 0;
    }

    virtual ~CGridBasic(void)
    {
        if(cell_list != 0)
        {
            delete[] cell_list;
            cell_list = 0;
        }

        if(pos_GasSpecRatios != 0)
        {
            delete[] pos_GasSpecRatios;
            pos_GasSpecRatios = 0;
        }

        if(pos_OpiateIDS != 0)
        {
            delete[] pos_OpiateIDS;
            pos_OpiateIDS = 0;
        }

        if(nr_dust_temp_sizes != 0)
            delete[] nr_dust_temp_sizes;

        if(nr_stochastic_sizes != 0)
            delete[] nr_stochastic_sizes;

        if(nr_stochastic_temps != 0)
            delete[] nr_stochastic_temps;

        if(size_skip != 0)
            delete[] size_skip;

        if(CextMeanTab != 0)
        {
            for(uint wID = 0; wID < max_wavelengths; wID++)
                delete[] CextMeanTab[wID];
            delete[] CextMeanTab;
            CextMeanTab = 0;
        }

        if(CabsMeanTab != 0)
        {
            for(uint wID = 0; wID < max_wavelengths; wID++)
                delete[] CabsMeanTab[wID];
            delete[] CabsMeanTab;
            CabsMeanTab = 0;
        }

        if(CscaMeanTab != 0)
        {
            for(uint wID = 0; wID < max_wavelengths; wID++)
                delete[] CscaMeanTab[wID];
            delete[] CscaMeanTab;
            CscaMeanTab = 0;
        }

        if(numberDensityTab != 0)
        {
            delete[] numberDensityTab;
            numberDensityTab = 0;
        }

        if(totalCellEmissionTab != 0)
        {
            delete[] totalCellEmissionTab;
            totalCellEmissionTab = 0;
        }
    }

    void printPhysicalParameters();

    void resetGridValues();

    double getCextMeanTab(uint cellID, uint wID) const;
    double getCabsMeanTab(uint cellID, uint wID) const;
    double getCscaMeanTab(uint cellID, uint wID) const;
    double getNumberDensityTab(uint cellID) const;
    double getTotalCellEmissionTab(uint cellID) const;

    void setCextMeanTab(double Cext, uint cellID, uint wID);
    void setCabsMeanTab(double Cext, uint cellID, uint wID);
    void setCscaMeanTab(double Cext, uint cellID, uint wID);
    void setNumberDensityTab(double density, uint cellID);
    void setTotalCellEmissionTab(double cell_emission, uint cellID);

    void initPreCalcTables(uint nr_used_wavelengths);

    double getTurbulentVelocity(cell_basic * cell);

    double getTurbulentVelocity(photon_package * pp);

    void updateDataRange(cell_basic * cell);

    void updateVelocity(cell_basic * cell, parameters & param);

    uint getDataOffset();

    uint getDataID();

    bool hasVelocityField();

    bool hasTurbulentVelocity();

    virtual Vector3D getCenter(const cell_basic & cell) const = 0;

    Vector3D getCenter(const photon_package & pp) const;

    uint getDataLength();

    ulong getMaxDataCells();

    uint getDataSize();

    virtual void printParameters() = 0;

    virtual bool goToNextCellBorder(photon_package * pp) = 0;

    virtual bool updateShortestDistance(photon_package * pp);

    void setOrientation(Vector3D n1, Vector3D n2, double _rot_angle1, double _rot_angle2);

    void getMagFieldInfo(const photon_package & pp, MagFieldInfo * mfo) const;

    virtual bool findStartingPoint(photon_package * pp) = 0;

    virtual void getLengths(uint bins, double & step_xy, double & off_xy) = 0;

    virtual bool positionPhotonInGrid(photon_package * pp) = 0;

    double getMinLength();

    double getMaxLength();

    virtual bool getPolarRTGridParameter(double max_len,
                                         double pixel_width,
                                         uint max_subpixel_lvl,
                                         dlist & _listR,
                                         uint & N_polar_r,
                                         uint *& N_polar_ph);

    virtual bool createCellList() = 0;

    cell_basic * getCellFromIndex(ulong i);

    void setSIConversionFactors(parameters & param);

    virtual bool writePlotFiles(string path, parameters & param) = 0;

    void setDataSize(uint sz);

    void setDustInformation(uint _nr_mixtures,
                            uint * _nr_dust_temp_sizes,
                            uint * _nr_stochastic_sizes,
                            uint * _nr_stochastic_temps);

    void setGasInformation(uint ** _level_to_pos, uint *** _line_to_pos);

    void setVelocityFieldNeeded(bool val);

    void setDataOffset(uint off);

    void setDataID(uint id);

    void setGasDensity(photon_package * pp, double dens);

    void setGasDensity(photon_package * pp, uint i_density, double dens);

    virtual bool isInside(const Vector3D & pos) const = 0;

    void setSpecLengthAsVector(bool val);

    bool specLengthIsVector();

    void updateSpecLength(photon_package * pp, double len);

    void updateSpecLength(cell_basic * cell, uint i_offset, StokesVector stokes) const;

    double getSpecLength(const cell_basic & cell, uint wID) const;

    double getSpecLength(const photon_package & pp, uint wID) const;

    void getSpecLength(const cell_basic & cell, uint wID, double * us, Vector3D * e_dir) const;

    void saveRadiationField();

    double getRadiationField(const cell_basic & cell, uint wID) const;

    double getRadiationField(const photon_package & pp, uint wID) const;

    double getRadiationFieldX(const cell_basic & cell, uint wID) const;

    double getRadiationFieldX(const photon_package & pp, uint wID) const;

    double getRadiationFieldY(const cell_basic & cell, uint wID) const;

    double getRadiationFieldY(const photon_package & pp, uint wID) const;

    double getRadiationFieldZ(const cell_basic & cell, uint wID) const;

    double getRadiationFieldZ(const photon_package & pp, uint wID) const;

    void getRadiationField(const photon_package & pp, uint w, double * us, Vector3D * e_dir) const;

    StokesVector getStokesFromRadiationField(const photon_package & pp, uint i_offset) const;

    void getRadiationFieldInterp(const photon_package & pp,
                                 double wavelength,
                                 double * us,
                                 Vector3D * e_dir) const;

    double getGZero(const cell_basic & cell) const;

    double getGZero(const photon_package & pp) const;

    double getUrad(const cell_basic & cell) const;

    double getUrad(const photon_package & pp) const;

    double getMu() const;

    void setRelOutsidePosition(photon_package * pp, double tx, double ty, double tz);

    void setRelOutsidePosition(photon_package * pp, double tx, double ty);

    void setRelDirection(photon_package * pp);

    virtual void setRndPositionInCell(photon_package * pp, CRandomGenerator * rand_gen);

    Vector3D rotateToCenter(const photon_package & pp, bool inv = false, bool phi_only = false) const;

    virtual Vector3D rotateToCenter(const photon_package & pp,
                                    Vector3D dir,
                                    bool inv = false,
                                    bool phi_only = false) const;

    void setDustTemperature(cell_basic * cell, uint i_density, uint a, double temp);

    void setDustTemperature(cell_basic * cell, uint i_density, double temp);

    void setDustTemperature(cell_basic * cell, double temp);

    void setDustTemperature(photon_package * pp, uint i_density, uint a, double temp);

    void setDustTemperature(photon_package * pp, uint i_density, double temp);

    void setDustTemperature(photon_package * pp, double temp);

    void setDustTempProbability(cell_basic * cell, uint i_density, uint a, uint t, double temp);

    void setPDAValue(cell_basic * cell, double val);

    double getMagMax();

    void setGasTemperature(cell_basic * cell, double temp);

    void setElectronTemperature(cell_basic * cell, double temp);

    void setThermalElectronDensity(cell_basic * cell, double dens);

    void setCRElectronDensity(cell_basic * cell, double dens);

    void setGammaMin(cell_basic * cell, double g_min);

    void setGammaMax(cell_basic * cell, double g_max);

    void setPowerLawIndex(cell_basic * cell, double p);

    void setAvgTheta(cell_basic * cell, double phi);

    void setAvgDir(cell_basic * cell, double dir);

    void setDustChoiceID(cell_basic * cell, uint dust_id);

    void setGasDensity(cell_basic * cell, double dens);

    void setGasDensity(cell_basic * cell, uint i_density, double dens);

    double getQBOffset(const cell_basic & cell, uint i_density) const;

    double getQBOffset(const photon_package & pp, uint i_density) const;

    double getQBOffset(const cell_basic & cell, uint i_density, uint a) const;

    double getQBOffset(const photon_package & pp, uint i_density, uint a) const;

    void setQBOffset(cell_basic * cell, uint i_density, uint a, double temp);

    void setQBOffset(cell_basic * cell, uint i_density, double temp);

    uint getNrAlignedRadii();

    double getAlignedRadius(const cell_basic & cell, uint i_density) const;

    double getAlignedRadius(const photon_package & pp, uint i_density) const;

    void setAlignedRadius(cell_basic * cell, uint i_density, double _a_alg);

    double getMinGrainRadius(const cell_basic & cell) const;

    double getMinGrainRadius(const photon_package & pp) const;

    double getMaxGrainRadius(const cell_basic & cell) const;

    double getMaxGrainRadius(const photon_package & pp) const;

    double getGrainSizeParam(const cell_basic & cell) const;

    double getGrainSizeParam(const photon_package & pp) const;

    uint getDustChoiceID(const photon_package & pp) const;

    uint getDustChoiceID(const cell_basic & cell) const;

    void getLineBroadening(const photon_package & pp, uint i_trans, LineBroadening * line_broadening) const;

    double getGaussA(const cell_basic & cell) const;

    double getGaussA(const photon_package & pp) const;

    double getVoigtA(const cell_basic & cell, uint i_line) const;

    double getVoigtA(const photon_package & pp, uint i_line) const;

    double getLvlPop(const cell_basic & cell, uint i_lvl, uint i_sublvl = 0) const;

    double getLvlPop(const photon_package & pp, uint i_lvl, uint i_sublvl = 0) const;

    void setVelocityField(cell_basic * cell, const Vector3D & vel);

    void setVelocityField(photon_package * pp, const Vector3D & vel);

    uint getCellID(cell_basic * cell);

    uint validateDataPositions(parameters & param);
    uint getDataIdsOffset(parameters & param);

    void setLvlPopLower(cell_basic * cell, uint i_line, uint i_sublvl, double lvl_lower);

    void setLvlPopLower(photon_package * pp, uint i_line, uint i_sublvl, double lvl_lower);

    void setLvlPopUpper(cell_basic * cell, uint i_line, uint i_sublvl, double lvl_upper);

    void setLvlPopUpper(photon_package * pp, uint i_line, uint i_sublvl, double lvl_upper);

    void setLvlPop(cell_basic * cell, uint i_lvl, uint i_sublvl, double lvl_pop);

    void setLvlPop(photon_package * pp, uint i_lvl, uint i_sublvl, double lvl_pop);

    void addToZeemanSublevel(cell_basic * cell, uint i_lvl, uint i_sublvl, double lvl_pop);

    void addToZeemanSublevel(photon_package * pp, uint i_lvl, uint i_sublvl, double lvl_pop);

    void setLineBroadening(cell_basic * cell, uint i_line_broad, const LineBroadening & line_broadening);

    void setGaussA(cell_basic * cell, double gauss_a);

    uint getOpiateID(cell_basic * cell);

    uint getOpiateID(const photon_package * pp);

    uint getOpiateID(photon_package pp);

    void setOpiateID(cell_basic * cell, uint id);

    double getOpiateTestData(photon_package * pp);

    double getOpiateTestData(cell_basic * cell);

    void setOpiateTestData(cell_basic * cell, double val);

    bool fillGridWithOpiateData(uint col_id);

    double getElectronTemperature(const photon_package & pp) const;

    double getElectronTemperature(const cell_basic & cell) const;

    double getThermalElectronDensity(const photon_package & pp) const;

    double getThermalElectronDensity(const cell_basic & cell) const;

    double getCRElectronDensity(const photon_package & pp) const;

    double getCRElectronDensity(const cell_basic & cell) const;

    double getGammaMin(const photon_package & pp) const;

    double getGammaMin(const cell_basic & cell) const;

    double getGammaMax(const photon_package & pp) const;

    double getGammaMax(const cell_basic & cell) const;

    double getPowerLawIndex(const photon_package & pp) const;

    double getPowerLawIndex(const cell_basic & cell) const;

    double getAvgTheta(const photon_package & pp) const;

    double getAvgTheta(const cell_basic & cell) const;

    double getAvgDir(const photon_package & pp) const;

    double getAvgDir(const cell_basic & cell) const;

    double getDustTemperature(const cell_basic & cell, uint i_density, uint a) const;

    double getDustTemperature(const photon_package & pp, uint i_density, uint a) const;

    double getDustTemperature(const cell_basic & cell, uint i_density) const;

    double getDustTemperature(const photon_package & pp, uint i_density) const;

    double getDustTemperature(const cell_basic & cell) const;

    double getDustTemperature(const photon_package & pp) const;

    double getDustTempProbability(const cell_basic & cell, uint i_density, uint a, uint t) const;

    double getDustTempProbability(const photon_package & pp, uint i_density, uint a, uint t) const;

    double getPDAValue(const cell_basic & cell) const;

    double getGasTemperature(const photon_package & pp) const;

    double getGasTemperature(const cell_basic & cell) const;

    void setPlaneParameter(uint plane_index,
                           double xy_step,
                           double off_xy,
                           double z_step,
                           double off_z,
                           double shift_z,
                           int j,
                           int k,
                           int l,
                           double & tx,
                           double & ty,
                           double & tz);

    void fillMidplaneBuffer(double tx, double ty, double tz, uint i_cell);

    bool writeSpecialLines(string path);
    bool writeAMIRAFiles(string path, parameters & param, uint bins);
    bool writeMidplaneFits(string data_path, parameters & param, uint bins, bool all = false);

    void updateMidplaneString(char * str_1, char * str_2, uint counter);

    string getDensityString(string quantity, uint counter);

    virtual bool saveBinaryGridFile(string filename) = 0;

    virtual bool saveBinaryGridFile(string filename, ushort id, ushort data_size) = 0;

    virtual bool loadGridFromBinaryFile(parameters & param) = 0;

    virtual bool loadGridFromBinaryFile(parameters & param, uint data_len) = 0;

    virtual double getVolume(const cell_basic & cell) const = 0;

    double getVolume(const photon_package & pp) const;

    double getGasDensity(const cell_basic & cell) const;

    double getGasDensity(const cell_basic & cell, uint i_density) const;

    double getGasDensity(const photon_package & pp) const;

    double getGasDensity(const photon_package & pp, uint i_density) const;

    double getGasNumberDensity(const cell_basic & cell) const;

    double getGasNumberDensity(const cell_basic & cell, uint i_density) const;

    double getGasNumberDensity(const photon_package & pp) const;

    double getGasNumberDensity(const photon_package & pp, uint i_density) const;

    double getGasMassDensity(const cell_basic & cell) const;

    double getGasMassDensity(const cell_basic & cell, uint i_density) const;

    double getGasMassDensity(const photon_package & pp) const;

    double getGasMassDensity(const photon_package & pp, uint i_density) const;

    bool useDustChoice();

    bool useConstantGrainSizes();

    bool useDustDensities();

    void setDustDensity(cell_basic * cell, double val);

    void setDustDensity(cell_basic * cell, uint i_density, double val);

    void setDustDensity(photon_package * pp, double val);

    void setDustDensity(photon_package * pp, uint i_density, double val);

    double getDustDensity(const cell_basic & cell) const;

    double getDustDensity(const cell_basic & cell, uint i_density) const;

    double getDustDensity(const photon_package & pp) const;

    double getDustDensity(const photon_package & pp, uint i_density) const;

    double getRelativeDustDensity(const cell_basic & cell, uint i_density) const;

    void adjustDustDensity(cell_basic * cell, uint i_density, double factor);

    virtual bool positionPhotonInGridTest(photon_package * pp);

    bool isVelocityFieldAvailable();

    bool isTurbulentVelocityAvailable();

    Vector3D getVelocityField(const photon_package & pp) const;

    Vector3D getVelocityField(const cell_basic & cell) const;

    double getCellAbundance(const photon_package & pp, uint id) const;

    double getCellAbundance(const cell_basic & cell, uint id) const;

    double getOpiateIDParameter(cell_basic * cell, uint id);

    Vector3D getMagField(const cell_basic & cell) const;

    Vector3D getMagField(const photon_package & pp) const;

    void setMagField(cell_basic * cell, const Vector3D & mag);

    virtual bool next(photon_package * pp) = 0;
    
    
    double getThetaSync(const photon_package & pp) const;

    double getThetaMag(const photon_package & pp) const;

    double getPhiMag(const photon_package & pp) const;

    double getTheta(const cell_basic & cell, Vector3D & dir) const;

    double getThetaPhoton(const photon_package & pp, Vector3D & dir) const;

    bool getDustIsMassDensity() const;

    bool getGasIsMassDensity() const;

    bool isRadiationFieldAvailable() const;

    double getTotalGasMass() const;

    void setDustTemperatureRange(double _min_dust_temp, double _max_dust_temp);

    void setalignedRadiusRange(double a_min, double a_max);

    bool doPDA(parameters & param, uint pda_id);

    uint getTemperatureFieldInformation() const;

    bool setDataPositionsVariable();

    bool createCompatibleTree();

    uint CheckSynchrotron(parameters & param);

    uint CheckOpiate(parameters & param);

    uint CheckTemp(parameters & param, uint & tmp_data_offset);

    uint CheckRat(parameters & param, uint & tmp_data_offset);

    uint CheckDustEmission(parameters & param);

    uint CheckDustScattering(parameters & param);

    uint CheckRadiationForce(parameters & param);

    uint CheckLineEmission(parameters & param);

    uint CheckProbing(parameters & param);

protected:
    // uint grid_type;
    ulong max_cells;
    uint max_data;
    Vector3D meanBdir;
    Vector3D meanVdir;

    char * basic_path;

    double max_gas_dens;
    double min_gas_dens;

    double max_dust_dens;
    double min_dust_dens;

    double max_gas_temp;
    double min_gas_temp;

    double max_dust_temp;
    double min_dust_temp;

    double max_larm_limit;
    double min_larm_limit;

    double min_delta;
    double max_delta;

    double max_mach;
    double min_mach;

    double aalg_min;
    double aalg_max;

    double a_min_min;
    double a_min_max;
    double a_max_min;
    double a_max_max;

    double size_param_min;
    double size_param_max;

    uint dust_id_min;
    uint dust_id_max;

    double min_pres;
    double max_pres;

    double max_vel;
    double min_vel;

    double max_mag;
    double min_mag;

    double max_value;

    double max_len;
    double min_len;

    double min_n_th;
    double min_T_e;
    double min_n_cr;
    double min_g_min;
    double min_g_max;
    double min_p;

    double max_n_th;
    double max_T_e;
    double max_n_cr;
    double max_g_min;
    double max_g_max;
    double max_p;

    double conv_length_in_SI, conv_dens_in_SI;
    double conv_Bfield_in_SI, conv_Vfield_in_SI;

    double total_gas_mass;
    double mu;

    dlist wl_list;

    int line_counter;
    uint char_counter;
    unsigned char ru[4];

    uint nrOfPlotPoints, nrOfPlotVectors;
    uint maxPlotLines;

    Vector3D ex, ey, ez;

    cell_basic ** cell_list;
    uint buffer_size;

    uint data_offset;
    uint dataID;
    uint data_len;

    uint nrOfDensRatios;
    uint nrOfOpiateIDs;

    uint nr_mixtures;
    uint nr_densities;
    uint size_gd_list;
    uint size_dd_list;
    uint multi_temperature_entries;
    uint stochastic_temperature_entries;
    uint * nr_dust_temp_sizes;
    uint * nr_stochastic_sizes;
    uint * nr_stochastic_temps;
    uint * size_skip;

    uint ** level_to_pos;
    uint *** line_to_pos;

    uilist data_pos_gd_list;
    uilist data_pos_dd_list;
    uilist data_pos_dt_list;
    uint data_pos_tg;
    uint data_pos_mx;
    uint data_pos_my;
    uint data_pos_mz;
    uint data_pos_vx;
    uint data_pos_vy;
    uint data_pos_vz;
    uint data_pos_px;
    uint data_pos_py;
    uint data_pos_pz;
    uilist data_pos_aalg_list;
    uint data_pos_amin;
    uint data_pos_amax;
    uint data_pos_size_param;
    //uint data_pos_ra;
    uint data_pos_id;

    uint data_pos_vt;
    uint data_pos_pda;
    uint data_pos_op;

    uint data_pos_n_th;
    uint data_pos_T_e;
    uint data_pos_n_cr;
    uint data_pos_g_min;
    uint data_pos_g_max;
    uint data_pos_p;

    uint nr_rad_field_comp;

    uilist data_pos_rx_list;
    uilist data_pos_ry_list;
    uilist data_pos_rz_list;
    uilist data_pos_rf_list;

    uint data_pos_avg_th;
    uint data_pos_avg_dir;

    double turbulent_velocity;

    uslist data_ids;
    uint * pos_GasSpecRatios;
    uint * pos_OpiateIDS;

    double rot_angle1, rot_angle2;

    bool plt_gas_dens;
    bool plt_mol_dens;
    bool plt_dust_dens;
    bool plt_gas_temp;
    bool plt_dust_temp;
    bool plt_mag;
    bool plt_vel;
    bool plt_rat;
    //bool plt_delta;
    bool plt_larm;
    //bool plt_mach;
    bool plt_dust_id;
    bool plt_amin;
    bool plt_amax;
    bool plt_size_param;
    bool plt_rad_field1;
    bool plt_u_rad;
    bool plt_g_zero1;
    bool plt_n_th;
    bool plt_T_e;
    bool plt_n_cr;
    bool plt_g_min;
    bool plt_g_max;
    bool plt_p;

    bool plt_avg_dir;
    bool plt_avg_th;

    bool dust_is_mass_density, gas_is_mass_density;
    bool velocity_field_needed;
    bool spec_length_as_vector;

    double delta0;
    double larm_f;

    double total_volume;
    double cell_volume;

    double ** buffer_gas_dens;
    double ** buffer_mol_dens;
    double ** buffer_dust_dens;
    double * buffer_gas_temp;
    double ** buffer_dust_temp;
    double ** buffer_rat;
    //double * buffer_delta;
    double * buffer_mag;
    double * buffer_mag_x;
    double * buffer_mag_y;
    double * buffer_mag_z;
    double * buffer_vel;
    double * buffer_vel_x;
    double * buffer_vel_y;
    double * buffer_vel_z;
    double * buffer_larm;
    //double * buffer_mach;
    double * buffer_dust_mixture;
    double * buffer_dust_amin;
    double * buffer_dust_amax;
    double * buffer_dust_size_param;
    double *** buffer_rad_field;
    double * buffer_g_zero1;
    double * buffer_u_rad;
    double * buffer_n_th;
    double * buffer_T_e;
    double * buffer_n_cr;
    double * buffer_g_min;
    double * buffer_g_max;
    double * buffer_p;

    double * buffer_avg_th;
    double * buffer_avg_dir;

    double ** CextMeanTab;
    double ** CabsMeanTab;
    double ** CscaMeanTab;
    double * numberDensityTab;
    double * totalCellEmissionTab;
    uint max_wavelengths;

    bool getPolarRTGridParameterWorker(double max_len,
                                       double pixel_width,
                                       uint max_subpixel_lvl,
                                       dlist & _listR,
                                       uint & N_polar_r,
                                       uint *& N_polar_phi,
                                       const uint &N_r,
                                       const double *listR);
};

#endif /* CGRID_BASIC_H */
