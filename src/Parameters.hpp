/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Typedefs.hpp"
#include "Vector3D.hpp"

class parameters
{
public:
    parameters()
    {
        path_grid = "";
        path_output = "";
        // phID = PH_HG;
        conv_l_in_SI = 1;
        conv_dH_in_SI = 1;
        conv_B_in_SI = 1;
        conv_V_in_SI = 1;
        //conv_mass_fraction = 1;
        align = 0;
        nr_ofThreads = 1;

        mu = 2.0;

        min_detector_pixel_x = MAX_UINT;
        max_detector_pixel_x = 0;
        min_detector_pixel_y = MAX_UINT;
        max_detector_pixel_y = 0;

        min_rot_angle_1 = 360;
        max_rot_angle_1 = 0;
        min_rot_angle_2 = 360;
        max_rot_angle_2 = 0;

        min_sidelength_x = 1e300;
        max_sidelength_x = 0;
        min_sidelength_y = 1e300;
        max_sidelength_y = 0;
        use_grid_sidelength_x = false;
        use_grid_sidelength_y = false;

        min_ray_map_shift_x = 1e300;
        max_ray_map_shift_x = 0;
        min_ray_map_shift_y = 1e300;
        max_ray_map_shift_y = 0;

        min_obs_distance = 1e300;
        max_obs_distance = -1e300;

        delta0 = 8.28e23 * 2.5e-12 * 1e8 * 1e-6 * 1e6;
        larm_f = 4.1e-19;

        nr_of_mc_lvl_pop_photons = 0;
        mc_lvl_pop_seed = 0;

        kepler_star_mass = 0;
        turbulent_velocity = 0;
        offset_min_gas_dens = 0;
        stochastic_heating_max_size = 0;

        task_id = 0;

        rt_grid_description = "";

        b_mrw = false;
        b_pda = false;
        b_enforced = false;
        peel_off = true;
        is_speed_of_sound = false;
        vel_maps_fits = false;
        
        r_sed = 0;
        dr_sed = 0;
                      
        heal_type = 0;
        
        fits_map_IDs=0;
        
        dust_offset = false;
        dust_gas_coupling = false;
        full_dust_temp = false;
        save_radiation_field = false;
        scattering_to_raytracing = false;
        split_dust_emision = false;
        sub_status = 0;
        //individual_dust_fractions = false;

        nr_ofISRFPhotons = 0;
        nr_ofDustPhotons = 0;

        cmd = -1;

        healpix_orientation = HEALPIX_YAXIS;

        start = MAX_UINT;
        stop = MAX_UINT;

        nr_ofInpAMIRAPoints = 0;
        nr_ofOutAMIRAPoints = 0;

        plot_inp_points = false;
        plot_out_points = false;
        write_radiation_field = 0;
        write_g_zero = false;
        write_dust_files = false;
        nr_ofInpMidDataPoints = 0;
        nr_ofOutMidDataPoints = 0;

        f_highJ = 0.25;
        Q_ref = 0.4;
        alpha_Q = 3.0;
        R_rayleigh = 1.0;

        f_cor = 0.6;
        adjTgas = 0;
        isrf_g_zero = 0;
        isrf_radius = 0;

        isrf_path = "";

        max_subpixel_lvl = 3;
        midplane_zoom = 1;
        max_dust_component_choice = 0;

        extinction_magnitude = 0;
        extinction_magnitude_wavelength = 0;
        extinction_i_mixture = MAX_UINT;

        reset_dust_files = false;
        acceptance_angle = 1.0;

        xymin = -1;
        xymax = 1;
        xysteps = 2;
        xy_bins = MAX_UINT;
        xylabel = "[a.u.]";
        autoscale = true;

        axis1.set(1, 0, 0);
        axis2.set(0, 1, 0);

        // opiate parmeters
        opiata_path_emi="";
        opiata_path_abs="";
    }

    ~parameters(){}

    string getOpiatePathEmission();

    string getOpiatePathAbsorption();

    string getOpiateSpec(uint pos);

    uint getNrOfOPIATESpecies();
    
    double getSEDRadius();
    
    double getSEDRange();

    const Vector3D & getAxis1() const;

    const Vector3D & getAxis2() const;

    bool plotInpMidPoints() const;

    bool plotOutMidPoints() const;

    dlist getMidplane3dParams() const;

    const uilist & getPlotList() const;

    bool isInPlotList(uint id);

    uint getInpMidDataPoints() const;

    uint getOutMidDataPoints() const;

    /*
    const dlist & getOpiateSequence() const;
    */

    uint getMidplaneZoom() const;

    int getCommand() const;

    bool isRatSimulation() const;
    
    bool isAMESimulation() const;

    bool isMonteCarloSimulation() const;

    bool isRaytracingSimulation() const;

    bool isTemperatureSimulation() const;

    double getStarMass(uint i) const;

    string getPathGrid() const;

    string getPathOutput() const;

    string getPathInput() const;

    uint getMinDetectorPixelX() const;

    uint getMaxDetectorPixelX() const;

    uint getMinDetectorPixelY() const;

    uint getMaxDetectorPixelY() const;

    double getMinDetectorAngle1() const;

    double getMaxDetectorAngle1() const;

    double getMinDetectorAngle2() const;

    double getMaxDetectorAngle2() const;

    double getMinSidelengthX() const;

    double getMaxSidelengthX() const;

    double getMinSidelengthY() const;

    double getMaxSidelengthY() const;

    bool getUseGridSidelengthX() const;

    bool getUseGridSidelengthY() const;

    double getMinMapShiftX() const;

    double getMinMapShiftY() const;

    double getMaxMapShiftX() const;

    double getMaxMapShiftY() const;

    double getSIConvLength() const;

    double getSIConvDH() const;

    double getDelta0() const;

    double getLarmF() const;

    double getSIConvBField() const;

    double getSIConvVField() const;

    bool getDustOffset() const;

    bool getDustGasCoupling() const;

    double getOffsetMinGasDensity() const;

    bool getDustTempMulti() const;

    double getSizeMin(uint i) const;

    double getSizeMax(uint i) const;

    double getMaterialDensity(uint i) const;

    bool getDustSource() const;

    bool getISRFSource() const;

    ullong getNrOfDustPhotons() const;

    dlist getDustMassFraction() const;
    
    double getDustMassFraction(uint i) const;
    
    bool getMapID_I() const;
    bool getMapID_QU() const;
    bool getMapID_V() const;
    bool getMapID_TAU() const;
    bool getMapID_N() const;
    uint getMapIDs() const;

    uint getAlign() const;

    bool getAligRANDOM() const;

    bool getAligPA() const;

    bool getAligNONPA() const;

    bool getAligIDG() const;

    bool getAligRAT() const;

    bool getAligGOLD() const;
    
    bool getAligRD() const;
    
    bool getAligKRAT() const;

    bool getAligINTERNAL() const;

    double getMu() const;

    bool getMRW() const;

    bool getPDA() const;

    bool getEnfScattering() const;

    double getStochasticHeatingMaxSize() const;

    bool getSaveRadiationField() const;

    bool getScatteringToRay() const;

    bool splitDustEmission() const;

    //bool getIndividualDustMassFractions() const;

    bool getIsSpeedOfSound() const;

    bool getPeelOff() const;

    double getForegroundExtinctionMagnitude() const;

    double getForegroundExtinctionWavelength() const;

    uint getForegroundExtinctionDustMixture() const;

    bool getVelFieldType() const;

    uint getWriteRadiationField() const;

    bool getWriteGZero() const;

    bool getWriteDustFiles() const;

    double getISRFGZero() const;

    double getISRFRadius() const;

    string getISRFPath() const;

    string getZeemanCatalog(uint i_species) const;

    uint getAlignmentMechanism() const;

    double getMinObserverDistance() const;

    double getMaxObserverDistance() const;

    double getKeplerStarMass() const;

    double getTurbulentVelocity() const;

    uint getMCLvlPopNrOfPhotons() const;

    uint getMCLvlPopSeed() const;

    uint getTaskID() const;

    uint getNrOfThreads() const;

    ullong getNrOfISRFPhotons() const;

    uint getNrOfMixtures() const;

    uilist getDustComponentChoices() const;

    /*
    uint getPhaseFunctionID() const;
    */

    uint getPhaseFunctionID(uint i) const;

    double getFHighJ() const;

    double getQref() const;

    double getAlphaQ() const;

    double getRayleighReductionFactor() const;

    double getFcorr() const;

    double getAdjTgas() const;

    uint getNrOfDiffuseSources() const;

    uint getNrOfPointSources() const;

    uint getNrOfLaserSources() const;

    uint getNrOfBackgroundSources() const;

    double getXYMin() const;

    double getXYMax() const;

    double getXYSteps() const;

    uint getXYBins() const;

    string getXYLabel() const;

    bool isAutoScale() const;

    uint getNrOfSources() const;

    uint getStart() const;

    uint getStop() const;
    
    void setSEDRadius(double r, double dr);

    void setOpiatePathEmission(string str);

    void setOpiatePathAbsorption(string str);

    void setXYMin(double val);

    void setXYMax(double val);

    void setXYSteps(double val);

    void setXYBins(uint val);

    void setXYLabel(string val);

    void setAutoScale(bool val);

    void setAxis1(double x, double y, double z);

    void setAxis2(double x, double y, double z);

    void setStart(uint val);

    void setStop(uint val);

    void setMu(double val);

    void setNrOfThreads(uint val);

    void setCommand(int val);

    void setISRF(string path, double g_zero = 0, double radius = 2);

    void setNrOfISRFPhotons(long val);

    void addOpiateSpec(string str);

    void AddDustComponentChoice(uint dust_component_choice);

    uint getMaxDustComponentChoice();

    void setTaskID(uint val);

    void addStarMass(double val);

    void setStarMass(dlist val);

    void setPathGrid(string val);

    void setPathInput(string val);

    void setPathOutput(string val);

    void setDustOffset(bool val);

    void setDustOffset(double _offset_min_gas_dens);

    void setDustGasCoupling(bool val);

    void setDustGasCoupling(double _offset_min_gas_dens);

    void setFullDustTemp(bool val);

    void setNrOfDustPhotons(long val);

    void setDelta0(double val);

    void addToPlotList(uint id);

    void setLarmF(double val);

    void setMRW(bool val);

    void setPDA(bool val);

    void setEnfScattering(bool val);

    void setIsSpeedOfSound(bool val);

    void setPeelOff(bool val);

    void setForegroundExtinction(double _extinction_magnitude,
                                 double _extinction_magnitude_wavelength = 0.55e-6,
                                 uint _extinction_i_mixture = MAX_UINT);

    void setInpAMIRAPoints(uint val);

    void setOutAMIRAPoints(uint val);

    void setInpMidPlot(bool val);

    void setOutMidPlot(bool val);

    void set3dMidplane(uint plane, uint nr_of_slices, double z_min, double z_max);

    void setInpMidDataPoints(uint val);

    void setOutMidDataPoints(uint val);

    void setMidplaneZoom(uint val);

    void setWriteRadiationField(uint val);

    void setWriteGZero(bool val);

    void setWriteDustFiles(bool val);

    void updateDetectorPixel(uint pixel_x, uint pixel_y);

    void updateDetectorAngles(double rot_angle_1, double rot_angle_2);

    void updateMapSidelength(double sidelength_x, double sidelength_y);

    void updateRayGridShift(double map_shift_x, double map_shift_y);

    void setStochasticHeatingMaxSize(double val);

    void setSaveRadiationField(bool val);

    void setScatteringToRay(bool val);

    void setSplitDustEmission(bool val);

    void setSIConvLength(double val);

    void setSIConvDH(double val);

    void setSIConvBField(double val);

    void setSIConvVField(double val);

    void updateSIConvLength(double val);

    void updateSIConvDH(double val);

    void updateSIConvBField(double val);

    void updateSIConvVField(double val);

    void setDustMassFraction(dlist val);

    //void setIndividualDustMassFractions(bool val);

    void addAlignmentMechanism(uint val);

    void updateObserverDistance(double val);

    void setKeplerStarMass(double val);

    void setMCLvlPopNrOfPhotons(uint val);

    void setMCLvlPopSeed(uint val);

    void setTurbulentVelocity(double val);

    void addZeemanCatalog(string val);

    /*
    void setPhaseFunctionID(uint val);
    */

    void setPhaseFunctionID(uint val, uint pos);

    void setFhighJ(double val);

    void setQref(double val);

    void setAlphaQ(double val);

    void setRayleighReductionFactor(double val);

    void setFcorr(double val);

    void setAdjTgas(double val);

    void setVelMapsFits(bool val);
    
    void addFitsID(uint val);
    
    void addHealType(uint type);
    
    void setProjectedFitsParam(dlist val);

    void setAcceptanceAngle(double angle);

    void setMaxSubpixelLvl(int val);

    dlist & getDustRayDetectors();

    dlist & getSyncRayDetectors();
    
    dlist & getFreeRayDetectors();
    
    dlist & getDustAMEDetectors();

    dlist & getOPIATERayDetectors();

    dlist & getPointSources();

    dlist & getLaserSources();

    strlist & getPointSourceStringList();

    string & getPointSourceString(uint i_str);

    strlist & getDiffuseSourceStringList();

    string & getDiffuseSourceString(uint i_str);

    dlist & getDiffuseSources();

    dlist & getBackgroundSources();

    strlist & getBackgroundSourceStringList();

    string & getBackgroundSourceString(uint i_str);

    void addDustRayDetector(dlist & val);

    void addOpiateRayDetector(dlist & val);

    void addSyncRayDetector(dlist & val);
    
    void addFreeRayDetector(dlist & val);
    
    void addDustAMEDetector(dlist & val);

    /*
    void addLineOpiateDetector(dlist & val);

    void setOpiateParamPath(string val);

    void setOpiateDataPath(string val);
    */

    void addLineRayDetector(dlist & val);

    void addDustMCDetector(dlist & val);

    uint getNrOfDustMCDetectors();

    dlist getDustMCDetectors();

    void addPointSource(dlist & val, string path);

    void addLaserSource(dlist & val);

    void addBackgroundSource(dlist & val);

    void addBackgroundSource(string path, dlist & val);

    void addBackgroundSource(string path);
    
    void setGauntPath(string path);
    
    string getGauntPath();
    
    void resetDustFiles();

    void addGasSpecies(string gas_species_path, string zeeman_path, dlist & val);

    uint getNrOfGasSpecies();

    uint getNrOfSpectralLines(uint i_species);

    int * getSpectralLines(uint i_species);

    uint getNrOfDustRayDetectors();

    uint getNrOfSyncRayDetectors();
    
    uint getNrOfFreeRayDetectors();
    
    uint getNrOfDustAMERayDetectors();

    void addSubStatus(int val);

    void setHealpixOrientation(uint val);

    uint getHealpixOrientation() const;

    const dlist & getLineRayDetector(uint i_species) const;

    const maplist & getLineRayDetectors() const;
    
    void addNanoGrains(dlist & pr);
    
    uint getNrOfNanoGrains();
    
    dlist findNanoGrainList(uint id);
    
    

    void addDustComponent(string path,
                          string size_key,
                          double fr,
                          double mat_dens,
                          double a_min,
                          double a_max,
                          dlist size_parameter);

    bool getVelMapsFits() const;
    uint getHealType() const;
    
    llist getProjectedFitsParam() const;
    long getProjectedFitsX() const;
    long getProjectedFitsY() const;
    bool  projectHealMaps() const;
    

    //uint getDustChoiceFromComponentId(uint i) const;

    uint getDustChoiceFromMixtureId(uint i) const;

    void printRTGridDescription();

    double getAcceptanceAngle() const;

    string getDustPath(uint i) const;

    int getSubStatus();

    double getDustFraction(uint i) const;

    string getDustSizeKeyword(uint i) const;

    dlist getDustSizeParameter(uint i_comp) const;

    string getGasSpeciesCatalogPath(uint i_species) const;

    double getGasSpeciesAbundance(uint i_species) const;

    uint getGasSpeciesLevelPopType(uint i_species) const;

    bool isGasSpeciesLevelPopMC() const;

    uint getMaxSubpixelLvl() const;

    uint getTotalNrOfDustComponents() const;

    void resetNrOfDustComponents();

    void addDiffuseSource(dlist & val, string path);

private:
    int cmd;
    string path_grid;
    string path_input;
    string path_output;

    uint max_subpixel_lvl;
    uint nr_ofThreads;
    uint task_id;

    double conv_l_in_SI;
    double conv_dH_in_SI;
    double conv_B_in_SI;
    double conv_V_in_SI;
    dlist conv_mass_fraction;
    double mu;

    double min_rot_angle_1, max_rot_angle_1;
    double min_rot_angle_2, max_rot_angle_2;

    double min_sidelength_x, max_sidelength_x;
    double min_sidelength_y, max_sidelength_y;
    bool use_grid_sidelength_x, use_grid_sidelength_y;

    double min_ray_map_shift_x, min_ray_map_shift_y;
    double max_ray_map_shift_x, max_ray_map_shift_y;

    uint align;
    uint min_detector_pixel_x, max_detector_pixel_x;
    uint min_detector_pixel_y, max_detector_pixel_y;
    uint nr_ofInpAMIRAPoints;
    uint nr_ofOutAMIRAPoints;
    uint nr_ofInpMidDataPoints;
    uint nr_ofOutMidDataPoints;
    uint midplane_zoom;
    uint max_dust_component_choice;

    bool plot_inp_points;
    bool plot_out_points;
    uint write_radiation_field;
    bool write_g_zero;
    bool write_dust_files;

    dlist midplane_3d_param;
    dlist star_mass;

    uilist plot_list;

    bool b_mrw;
    bool b_pda;
    bool b_enforced;
    bool is_speed_of_sound;
    bool peel_off;
    bool vel_maps_fits;
    
    uint heal_type;
    llist param_projection;
    
    uint fits_map_IDs;
    
    double r_sed;
    double dr_sed;
    
    bool dust_offset, dust_gas_coupling;
    bool full_dust_temp, save_radiation_field;
    bool scattering_to_raytracing;
    bool split_dust_emision;
    //bool individual_dust_fractions;

    strlist zeeman_catalog_path;

    // uint phID;
    uilist phIDs;

    double min_obs_distance, max_obs_distance;
    double kepler_star_mass, turbulent_velocity;
    double stochastic_heating_max_size;
    double delta0;
    double larm_f;
    double acceptance_angle;
    double offset_min_gas_dens;

    double extinction_magnitude;
    double extinction_magnitude_wavelength;
    uint extinction_i_mixture;

    uint nr_of_mc_lvl_pop_photons;
    uint mc_lvl_pop_seed;

    uint healpix_orientation;

    maplist line_ray_detector_list;

    string rt_grid_description;

    dlist dust_mc_detectors;
    dlist dust_ray_detectors;
    dlist sync_ray_detectors;
    dlist free_ray_detectors1;
    dlist dust_ame_detectors;

    dlist point_sources;
    dlist diffuse_sources;
    dlist laser_sources;
    dlist background_sources;
    dlist gas_species_abundance;

    uilist gas_species_level_pop_type;

    strlist gas_species_cat_path;
    strlist point_sources_str;
    strlist diffuse_sources_str;
    strlist background_sources_path;
    string isrf_path;

    double xymin, xymax, xysteps;
    uint xy_bins;
    string xylabel;
    bool autoscale;
    int sub_status;

    Vector3D axis1, axis2;

    uint start;
    uint stop;

    double f_highJ;
    double f_cor;
    double Q_ref;
    double alpha_Q;
    double R_rayleigh;

    double adjTgas;
    double isrf_g_zero;
    double isrf_radius;

    ullong nr_ofISRFPhotons;
    ullong nr_ofDustPhotons;

    dlist dust_fractions;
    dlist material_density;
    dlist a_min_global;
    dlist a_max_global;
    maplist size_parameter_map;
    uilist dust_choices;
    //uilist component_id_to_choice;

    strlist dust_paths;
    strlist size_keywords;
    
    dlist nano_grains;

    bool reset_dust_files;

    // opiate parameters
    dlist opiate_ray_detectors;
    strlist opiate_spec_ids;
    
    string gaunt_path;
    
    string opiata_path_emi;
    string opiata_path_abs;
};

#endif /* PARAMETERS_H */
