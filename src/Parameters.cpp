/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "Parameters.hpp"

string parameters::getOpiatePathEmission()
{
    return opiata_path_emi;
}

string parameters::getOpiatePathAbsorption()
{
    return opiata_path_abs;
}

string parameters::getOpiateSpec(uint pos)
{
    return opiate_spec_ids[pos];
}

uint parameters::getNrOfOPIATESpecies()
{
    return uint(opiate_spec_ids.size());
}

double parameters::getSEDRadius()
{
    return r_sed;
}

double parameters::getSEDRange()
{
    return dr_sed;
}

const Vector3D & parameters::getAxis1() const
{
    return axis1;
}

const Vector3D & parameters::getAxis2() const
{
    return axis2;
}

bool parameters::plotInpMidPoints() const
{
    return plot_inp_points;
}

bool parameters::plotOutMidPoints() const
{
    return plot_out_points;
}

dlist parameters::getMidplane3dParams() const
{
    return midplane_3d_param;
}

const uilist & parameters::getPlotList() const
{
    return plot_list;
}

bool parameters::isInPlotList(uint id)
{
    if(plot_list.empty())
        return true;

    return (find(plot_list.begin(), plot_list.end(), id) != plot_list.end());
}

uint parameters::getInpMidDataPoints() const
{
    return nr_ofInpMidDataPoints;
}

uint parameters::getOutMidDataPoints() const
{
    return nr_ofOutMidDataPoints;
}

/*const dlist & parameters::getOpiateSequence() const
{
    return line_opiate_detectors;
}*/

uint parameters::getMidplaneZoom() const
{
    return midplane_zoom;
}

int parameters::getCommand() const
{
    return cmd;
}

bool parameters::isRatSimulation() const
{
    if(getCommand() == CMD_RAT || getCommand() == CMD_TEMP_RAT)
        return true;
    return false;
}

bool parameters::isAMESimulation() const
{
    if(getCommand() == CMD_RAT || getCommand() == CMD_TEMP_RAT)
    {
        if(!nano_grains.empty())
            return true;
    }
        
    return false;
}

bool parameters::isMonteCarloSimulation() const
{
    if(getCommand() == CMD_TEMP || getCommand() == CMD_TEMP_RAT ||
        getCommand() == CMD_RAT || getCommand() == CMD_DUST_SCATTERING)
        return true;
    return false;
}

bool parameters::isRaytracingSimulation() const
{
    if(getCommand() == CMD_OPIATE || getCommand() == CMD_DUST_EMISSION ||
        getCommand() == CMD_SYNCHROTRON || getCommand() == CMD_LINE_EMISSION || 
            getCommand() == CMD_FREE_FREE || getCommand() == CMD_AME_EMISSION)
        return true;
    
    return false;
}

bool parameters::isTemperatureSimulation() const
{
    if(getCommand() == CMD_TEMP || getCommand() == CMD_TEMP_RAT)
        return true;
    return false;
}

double parameters::getStarMass(uint i) const
{
    return star_mass[i];
}

string parameters::getPathGrid() const
{
    return path_grid;
}

string parameters::getPathOutput() const
{
    return path_output;
}

string parameters::getPathInput() const
{
    return path_input;
}

uint parameters::getMinDetectorPixelX() const
{
    return min_detector_pixel_x;
}

uint parameters::getMaxDetectorPixelX() const
{
    return max_detector_pixel_x;
}

uint parameters::getMinDetectorPixelY() const
{
    return min_detector_pixel_y;
}

uint parameters::getMaxDetectorPixelY() const
{
    return max_detector_pixel_y;
}

double parameters::getMinDetectorAngle1() const
{
    return min_rot_angle_1;
}

double parameters::getMaxDetectorAngle1() const
{
    return max_rot_angle_1;
}

double parameters::getMinDetectorAngle2() const
{
    return min_rot_angle_2;
}

double parameters::getMaxDetectorAngle2() const
{
    return max_rot_angle_2;
}

double parameters::getMinSidelengthX() const
{
    return min_sidelength_x;
}

double parameters::getMaxSidelengthX() const
{
    return max_sidelength_x;
}

double parameters::getMinSidelengthY() const
{
    return min_sidelength_y;
}

double parameters::getMaxSidelengthY() const
{
    return max_sidelength_y;
}

bool parameters::getUseGridSidelengthX() const
{
    return use_grid_sidelength_x;
}

bool parameters::getUseGridSidelengthY() const
{
    return use_grid_sidelength_y;
}

double parameters::getMinMapShiftX() const
{
    return min_ray_map_shift_x;
}

double parameters::getMinMapShiftY() const
{
    return min_ray_map_shift_y;
}

double parameters::getMaxMapShiftX() const
{
    return max_ray_map_shift_x;
}

double parameters::getMaxMapShiftY() const
{
    return max_ray_map_shift_y;
}

double parameters::getSIConvLength() const
{
    return conv_l_in_SI;
}

double parameters::getSIConvDH() const
{
    return conv_dH_in_SI;
}

double parameters::getDelta0() const
{
    return delta0;
}

double parameters::getLarmF() const
{
    return larm_f;
}

double parameters::getSIConvBField() const
{
    return conv_B_in_SI;
}

double parameters::getSIConvVField() const
{
    return conv_V_in_SI;
}

bool parameters::getDustOffset() const
{
    return dust_offset;
}

bool parameters::getDustGasCoupling() const
{
    return dust_gas_coupling;
}

double parameters::getOffsetMinGasDensity() const
{
    return offset_min_gas_dens;
}

bool parameters::getDustTempMulti() const
{
    return full_dust_temp;
}

double parameters::getSizeMin(uint i) const
{
    return a_min_global[i];
}

double parameters::getSizeMax(uint i) const
{
    return a_max_global[i];
}

double parameters::getMaterialDensity(uint i) const
{
    return material_density[i];
}

bool parameters::getDustSource() const
{
    return nr_ofDustPhotons > 0;
}

bool parameters::getISRFSource() const
{
    return nr_ofISRFPhotons > 0;
}

ullong parameters::getNrOfDustPhotons() const
{
    return nr_ofDustPhotons;
}

dlist parameters::getDustMassFraction() const
{
    return conv_mass_fraction;
}

double parameters::getDustMassFraction(uint i) const
{
    return conv_mass_fraction[i];
}

uint parameters::getMapIDs() const
{
    return fits_map_IDs;
}

bool parameters::getMapID_I() const
{
    return (fits_map_IDs == 0 ) || (fits_map_IDs & MAP_ID_I) == MAP_ID_I;
}

bool parameters::getMapID_QU() const
{
    return (fits_map_IDs == 0 ) || (fits_map_IDs & MAP_ID_QU) == MAP_ID_QU;
}

bool parameters::getMapID_V() const
{
    return (fits_map_IDs == 0 ) || (fits_map_IDs & MAP_ID_V) == MAP_ID_V;
}

bool parameters::getMapID_TAU() const
{
    return (fits_map_IDs == 0 ) || (fits_map_IDs & MAP_ID_TAU) == MAP_ID_TAU;
}

bool parameters::getMapID_N() const
{
    return (fits_map_IDs == 0 ) || (fits_map_IDs & MAP_ID_N) == MAP_ID_N;
}


uint parameters::getAlign() const
{
    return align;
}

bool parameters::getAligRANDOM() const
{
    return align == 0;
}

bool parameters::getAligPA() const
{
    return (align & ALIG_PA) == ALIG_PA;
}

bool parameters::getAligNONPA() const
{
    return (align & ALIG_NONPA) == ALIG_NONPA;
}

bool parameters::getAligIDG() const
{
    return (align & ALIG_IDG) == ALIG_IDG;
}

bool parameters::getAligRAT() const
{
    return (align & ALIG_RAT) == ALIG_RAT;
}

bool parameters::getAligGOLD() const
{
    return (align & ALIG_GOLD) == ALIG_GOLD;
}

bool parameters::getAligRD() const
{
    return (align & ALIG_RD) == ALIG_RD;
}

bool parameters::getAligKRAT() const
{
    return (align & ALIG_KRAT) == ALIG_KRAT;
}

bool parameters::getAligINTERNAL() const
{
    return (align & ALIG_INTERNAL) == ALIG_INTERNAL;
}

double parameters::getMu() const
{
    return mu;
}

bool parameters::getMRW() const
{
    return b_mrw;
}

bool parameters::getPDA() const
{
    return b_pda;
}

bool parameters::getEnfScattering() const
{
    return b_enforced;
}

double parameters::getStochasticHeatingMaxSize() const
{
    return stochastic_heating_max_size;
}

bool parameters::getSaveRadiationField() const
{
    return save_radiation_field;
}

bool parameters::getScatteringToRay() const
{
    return scattering_to_raytracing;
}

bool parameters::splitDustEmission() const
{
    return split_dust_emision;
}

/*bool parameters::getIndividualDustMassFractions() const
{
    return individual_dust_fractions;
}*/

bool parameters::getIsSpeedOfSound() const
{
    return is_speed_of_sound;
}

bool parameters::getPeelOff() const
{
    return peel_off;
}

double parameters::getForegroundExtinctionMagnitude() const
{
    return extinction_magnitude;
}

double parameters::getForegroundExtinctionWavelength() const
{
    return extinction_magnitude_wavelength;
}

uint parameters::getForegroundExtinctionDustMixture() const
{
    return extinction_i_mixture;
}

bool parameters::getVelFieldType() const
{
    if(kepler_star_mass == 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

uint parameters::getWriteRadiationField() const
{
    return write_radiation_field;
}

bool parameters::getWriteGZero() const
{
    return write_g_zero;
}

bool parameters::getWriteDustFiles() const
{
    return write_dust_files;
}

double parameters::getISRFGZero() const
{
    return isrf_g_zero;
}

double parameters::getISRFRadius() const
{
    return isrf_radius;
}

string parameters::getISRFPath() const
{
    return isrf_path;
}

string parameters::getZeemanCatalog(uint i_species) const
{
    return zeeman_catalog_path[i_species];
}

uint parameters::getAlignmentMechanism() const
{
    return align;
}

double parameters::getMinObserverDistance() const
{
    return min_obs_distance;
}

double parameters::getMaxObserverDistance() const
{
    return max_obs_distance;
}

double parameters::getKeplerStarMass() const
{
    return kepler_star_mass;
}

double parameters::getTurbulentVelocity() const
{
    return turbulent_velocity;
}

uint parameters::getMCLvlPopNrOfPhotons() const
{
    return nr_of_mc_lvl_pop_photons;
}

uint parameters::getMCLvlPopSeed() const
{
    return mc_lvl_pop_seed;
}

uint parameters::getTaskID() const
{
    return task_id;
}

uint parameters::getNrOfThreads() const
{
    return nr_ofThreads;
}

ullong parameters::getNrOfISRFPhotons() const
{
    return nr_ofISRFPhotons;
}

uint parameters::getNrOfMixtures() const
{
    if(!dust_choices.empty())
        return dust_choices.size();
    else
        return 1;
}

uilist parameters::getDustComponentChoices() const
{
    return dust_choices;
}

/*
uint parameters::getPhaseFunctionID() const
{
    return phID;
}
*/

uint parameters::getPhaseFunctionID(uint i) const
{
    if(phIDs.size() == 1)
        return phIDs[0];
    if(i < phIDs.size())
        return phIDs[i];
    return PH_ISO;
}

double parameters::getFHighJ() const
{
    return f_highJ;
}

double parameters::getQref() const
{
    return Q_ref;
}

double parameters::getAlphaQ() const
{
    return alpha_Q;
}

double parameters::getRayleighReductionFactor() const
{
    return R_rayleigh;
}

double parameters::getFcorr() const
{
    return f_cor;
}

double parameters::getAdjTgas() const
{
    return adjTgas;
}

uint parameters::getNrOfDiffuseSources() const
{
    return uint(diffuse_sources.size() / NR_OF_DIFF_SOURCES);
}

uint parameters::getNrOfPointSources() const
{
    return uint(point_sources.size() / NR_OF_POINT_SOURCES);
}

uint parameters::getNrOfLaserSources() const
{
    return uint(laser_sources.size() / NR_OF_LASER_SOURCES);
}

uint parameters::getNrOfBackgroundSources() const
{
    return uint(background_sources.size() / NR_OF_BG_SOURCES);
}

double parameters::getXYMin() const
{
    return xymin;
}

double parameters::getXYMax() const
{
    return xymax;
}

double parameters::getXYSteps() const
{
    return xysteps;
}

uint parameters::getXYBins() const
{
    return xy_bins;
}

string parameters::getXYLabel() const
{
    return xylabel;
}

bool parameters::isAutoScale() const
{
    return autoscale;
}

uint parameters::getNrOfSources() const
{
    uint res = getNrOfPointSources() + getNrOfDiffuseSources() + getNrOfBackgroundSources() +
                getNrOfLaserSources();

    if(nr_ofDustPhotons > 0)
        res++;

    if(nr_ofISRFPhotons > 0)
        res++;

    // Add gas source for levl population calculation
    if(isGasSpeciesLevelPopMC())
        res++;

    return res;
}

uint parameters::getStart() const
{
    return start;
}

uint parameters::getStop() const
{
    return stop;
}

void parameters::setSEDRadius(double r, double dr)
{
    r_sed=r;
    dr_sed=dr;
}

void parameters::setOpiatePathEmission(string str)
{
    opiata_path_emi=str;
}

void parameters::setOpiatePathAbsorption(string str)
{
    opiata_path_abs=str;
}


void parameters::setXYMin(double val)
{
    xymin = val;
}

void parameters::setXYMax(double val)
{
    xymax = val;
}

void parameters::setXYSteps(double val)
{
    xysteps = val;
}

void parameters::setXYBins(uint val)
{
    xy_bins = val;
}

void parameters::setXYLabel(string val)
{
    xylabel = val;
}

void parameters::setAutoScale(bool val)
{
    autoscale = val;
}

void parameters::setAxis1(double x, double y, double z)
{
    axis1.set(x, y, z);
    axis1.normalize();
}

void parameters::setAxis2(double x, double y, double z)
{
    axis2.set(x, y, z);
    axis2.normalize();
}

void parameters::setStart(uint val)
{
    start = val;
}

void parameters::setStop(uint val)
{
    stop = val;
}

void parameters::setMu(double val)
{
    mu = val;
}

void parameters::setNrOfThreads(uint val)
{
    nr_ofThreads = val;
}

void parameters::setCommand(int val)
{
    cmd = val;
}

void parameters::setISRF(string path, double g_zero, double radius)
{
    isrf_path = path;
    isrf_g_zero = g_zero;
    isrf_radius = radius;
}

void parameters::setNrOfISRFPhotons(long val)
{
    nr_ofISRFPhotons = val;
}

void parameters::addOpiateSpec(string str)
{
    opiate_spec_ids.push_back(str);
}

void parameters::AddDustComponentChoice(uint dust_component_choice)
{
    // If the dust component choice of a dust component is already loaded,
    // add the "dust_component_choice" to the list of component_id_to_choice but not
    // to the dust_choices.
    /*if(!component_id_to_choice.empty())
        for(uint a = 0; a < component_id_to_choice.size(); a++)
            if(component_id_to_choice[a] == dust_component_choice)
            {
                component_id_to_choice.push_back(dust_component_choice);
                return;
            }*/

    // If the "dust_component_choice" is not yet used, add it to the amount of
    // available choices.
    //component_id_to_choice.push_back(dust_component_choice);
    dust_choices.push_back(dust_component_choice);

    // Sort dust choices
    //sort(dust_choices.begin(), dust_choices.end());

    // Update the highest value of the dust choice ids.
    dust_component_choice = max(dust_component_choice,max_dust_component_choice);
}

uint parameters::getMaxDustComponentChoice()
{
    return max_dust_component_choice;
}

void parameters::setTaskID(uint val)
{
    task_id = val;
}

void parameters::addStarMass(double val)
{
    star_mass.push_back(val);
}

void parameters::setStarMass(dlist val)
{
    star_mass = val;
}

void parameters::setPathGrid(string val)
{
    path_grid = val;
}

void parameters::setPathInput(string val)
{
    path_input = val;
}

void parameters::setPathOutput(string val)
{
    path_output = val;
}

void parameters::setDustOffset(bool val)
{
    dust_offset = val;
}

void parameters::setDustOffset(double _offset_min_gas_dens)
{
    dust_offset = true;
    offset_min_gas_dens = _offset_min_gas_dens;
}

void parameters::setDustGasCoupling(bool val)
{
    dust_gas_coupling = val;
}

void parameters::setDustGasCoupling(double _offset_min_gas_dens)
{
    dust_gas_coupling = true;
    offset_min_gas_dens = _offset_min_gas_dens;
}

void parameters::setFullDustTemp(bool val)
{
    full_dust_temp = val;
}

void parameters::setNrOfDustPhotons(long val)
{
    nr_ofDustPhotons = val;
}

void parameters::setDelta0(double val)
{
    delta0 = val;
}

void parameters::addToPlotList(uint id)
{
    plot_list.push_back(id);
}

void parameters::setLarmF(double val)
{
    larm_f = val;
}

void parameters::setMRW(bool val)
{
    b_mrw = val;
}

void parameters::setPDA(bool val)
{
    b_pda = val;
}

void parameters::setEnfScattering(bool val)
{
    b_enforced = val;
}

void parameters::setIsSpeedOfSound(bool val)
{
    is_speed_of_sound = val;
}

void parameters::setPeelOff(bool val)
{
    peel_off = val;
}

void parameters::setForegroundExtinction(double _extinction_magnitude,
                                double _extinction_magnitude_wavelength,
                                uint _extinction_i_mixture)
{
    extinction_magnitude = _extinction_magnitude;
    extinction_magnitude_wavelength = _extinction_magnitude_wavelength;
    extinction_i_mixture = _extinction_i_mixture;
}

void parameters::setInpAMIRAPoints(uint val)
{
    nr_ofInpAMIRAPoints = val;
}

void parameters::setOutAMIRAPoints(uint val)
{
    nr_ofOutAMIRAPoints = val;
}

void parameters::setInpMidPlot(bool val)
{
    plot_inp_points = val;
}

void parameters::setOutMidPlot(bool val)
{
    plot_out_points = val;
}

void parameters::set3dMidplane(uint plane, uint nr_of_slices, double z_min, double z_max)
{
    if(midplane_3d_param.size() == 0)
    {
        midplane_3d_param.push_back(plane);
        midplane_3d_param.push_back(nr_of_slices);
        midplane_3d_param.push_back(z_min);
        midplane_3d_param.push_back(z_max);
    }
    else
        cout << INFO_LINE << "Multiple 3D midplane plot commands found. The first command will be considered!" << endl;
}

void parameters::setInpMidDataPoints(uint val)
{
    nr_ofInpMidDataPoints = val;
}

void parameters::setOutMidDataPoints(uint val)
{
    nr_ofOutMidDataPoints = val;
}

void parameters::setMidplaneZoom(uint val)
{
    midplane_zoom = val;
}

void parameters::setWriteRadiationField(uint val)
{
    write_radiation_field = val;
}

void parameters::setWriteGZero(bool val)
{
    write_g_zero = val;
}

void parameters::setWriteDustFiles(bool val)
{
    write_dust_files = val;
}

void parameters::updateDetectorPixel(uint pixel_x, uint pixel_y)
{
    if(min_detector_pixel_x > pixel_x)
        min_detector_pixel_x = pixel_x;

    if(max_detector_pixel_x < pixel_x)
        max_detector_pixel_x = pixel_x;

    if(min_detector_pixel_y > pixel_y)
        min_detector_pixel_y = pixel_y;

    if(max_detector_pixel_y < pixel_y)
        max_detector_pixel_y = pixel_y;
}

void parameters::updateDetectorAngles(double rot_angle_1, double rot_angle_2)
{
    if(min_rot_angle_1 > rot_angle_1)
        min_rot_angle_1 = rot_angle_1;

    if(max_rot_angle_1 < rot_angle_1)
        max_rot_angle_1 = rot_angle_1;

    if(min_rot_angle_2 > rot_angle_2)
        min_rot_angle_2 = rot_angle_2;

    if(max_rot_angle_2 < rot_angle_2)
        max_rot_angle_2 = rot_angle_2;
}

void parameters::updateMapSidelength(double sidelength_x, double sidelength_y)
{
    if(sidelength_x != -1)
    {
        if(min_sidelength_x > sidelength_x)
            min_sidelength_x = sidelength_x;
        if(max_sidelength_x < sidelength_x)
            max_sidelength_x = sidelength_x;
    }
    else
        use_grid_sidelength_x = true;
    if(sidelength_y != -1)
    {
        if(min_sidelength_y > sidelength_y)
            min_sidelength_y = sidelength_y;
        if(max_sidelength_y < sidelength_y)
            max_sidelength_y = sidelength_y;
    }
    else
        use_grid_sidelength_y = true;
}

void parameters::updateRayGridShift(double map_shift_x, double map_shift_y)
{
    if(abs(min_ray_map_shift_x) > abs(map_shift_x))
        min_ray_map_shift_x = map_shift_x;

    if(abs(max_ray_map_shift_x) < abs(map_shift_x))
        max_ray_map_shift_x = map_shift_x;

    if(abs(min_ray_map_shift_y) > abs(map_shift_y))
        min_ray_map_shift_y = map_shift_y;

    if(abs(max_ray_map_shift_y) < abs(map_shift_y))
        max_ray_map_shift_y = map_shift_y;
}

void parameters::setStochasticHeatingMaxSize(double val)
{
    stochastic_heating_max_size = val;
}

void parameters::setSaveRadiationField(bool val)
{
    save_radiation_field = val;
}

void parameters::setScatteringToRay(bool val)
{
    scattering_to_raytracing = val;
}

void parameters::setSplitDustEmission(bool val)
{
    split_dust_emision = val;
}

void parameters::setSIConvLength(double val)
{
    conv_l_in_SI = val;
}

void parameters::setSIConvDH(double val)
{
    conv_dH_in_SI = val;
}

void parameters::setSIConvBField(double val)
{
    conv_B_in_SI = val;
}

void parameters::setSIConvVField(double val)
{
    conv_V_in_SI = val;
}

void parameters::updateSIConvLength(double val)
{
    if(conv_l_in_SI != 1 && val != 1)
    {
        cout << INFO_LINE << "<conv_len> may not be used multiple times!" << endl;
        cout << "  -> No problem if <path_grid_cgs> was used!" << endl;
    }
    conv_l_in_SI *= val;
}

void parameters::updateSIConvDH(double val)
{
    if(conv_dH_in_SI != 1 && val != 1)
    {
        cout << INFO_LINE << "<conv_dens> may not be used multiple times!" << endl;
        cout << "  -> No problem if <path_grid_cgs> was used!" << endl;
    }
    conv_dH_in_SI *= val;
}

void parameters::updateSIConvBField(double val)
{
    if(conv_B_in_SI != 1 && val != 1)
    {
        cout << INFO_LINE << "<conv_mag> may not be used multiple times!" << endl;
        cout << "  -> No problem if <path_grid_cgs> was used!" << endl;
    }
    conv_B_in_SI *= val;
}

void parameters::updateSIConvVField(double val)
{
    if(conv_V_in_SI != 1 && val != 1)
    {
        cout << INFO_LINE << "<conv_vel> may not be used multiple times!" << endl;
        cout << "  -> No problem if <path_grid_cgs> was used!" << endl;
    }
    conv_V_in_SI *= val;
}

void parameters::setDustMassFraction(dlist val)
{
    conv_mass_fraction = val;
}

/*void parameters::setIndividualDustMassFractions(bool val)
{
    individual_dust_fractions = val;
}*/

void parameters::addAlignmentMechanism(uint val)
{
    align |= val;
}

void parameters::updateObserverDistance(double val)
{
    if(min_obs_distance > val)
        min_obs_distance = val;

    if(max_obs_distance < val)
        max_obs_distance = val;
}

void parameters::setKeplerStarMass(double val)
{
    kepler_star_mass = val;
}

void parameters::setMCLvlPopNrOfPhotons(uint val)
{
    nr_of_mc_lvl_pop_photons = val;
}

void parameters::setMCLvlPopSeed(uint val)
{
    if(val > 0)
        mc_lvl_pop_seed = val;
}

void parameters::setTurbulentVelocity(double val)
{
    turbulent_velocity = val;
}

void parameters::addZeemanCatalog(string val)
{
    zeeman_catalog_path.push_back(val);
}

/*
void parameters::setPhaseFunctionID(uint val)
{
    phID = val;
}
*/

void parameters::setPhaseFunctionID(uint val, uint pos)
{
    if(pos >= phIDs.size())
        phIDs.resize(pos+1);
    phIDs[pos] = val;
}

void parameters::setFhighJ(double val)
{
    f_highJ = val;
}

void parameters::setQref(double val)
{
    Q_ref = val;
}

void parameters::setAlphaQ(double val)
{
    alpha_Q = val;
}

void parameters::setRayleighReductionFactor(double val)
{
    R_rayleigh = val;
}

void parameters::setFcorr(double val)
{
    f_cor = val;
}

void parameters::setAdjTgas(double val)
{
    adjTgas = val;
}

void parameters::setVelMapsFits(bool val)
{
    vel_maps_fits = val;
}

void parameters::addFitsID(uint val)
{
    fits_map_IDs |= val;
}

void parameters::addHealType(uint type)
{
    heal_type |= type;
}

void parameters::setProjectedFitsParam(dlist val)
{
    param_projection.clear();
            
    param_projection.push_back(long(val[0]));
    param_projection.push_back(long(val[1]));
}

void parameters::setAcceptanceAngle(double angle)
{
    acceptance_angle = angle;
}

void parameters::setMaxSubpixelLvl(int val)
{
    max_subpixel_lvl = val;
}

dlist & parameters::getDustRayDetectors()
{
    return dust_ray_detectors;
}

dlist & parameters::getSyncRayDetectors()
{
    return sync_ray_detectors;
}

dlist & parameters::getFreeRayDetectors()
{
    return free_ray_detectors1;
}

dlist & parameters::getDustAMEDetectors()
{
    return dust_ame_detectors;
}

dlist & parameters::getOPIATERayDetectors()
{
    return opiate_ray_detectors;
}

dlist & parameters::getPointSources()
{
    return point_sources;
}

dlist & parameters::getLaserSources()
{
    return laser_sources;
}

strlist & parameters::getPointSourceStringList()
{
    return point_sources_str;
}

string & parameters::getPointSourceString(uint i_str)
{
    return point_sources_str[i_str];
}

strlist & parameters::getDiffuseSourceStringList()
{
    return diffuse_sources_str;
}

string & parameters::getDiffuseSourceString(uint i_str)
{
    return diffuse_sources_str[i_str];
}

dlist & parameters::getDiffuseSources()
{
    return diffuse_sources;
}

dlist & parameters::getBackgroundSources()
{
    return background_sources;
}

strlist & parameters::getBackgroundSourceStringList()
{
    return background_sources_path;
}

string & parameters::getBackgroundSourceString(uint i_str)
{
    return background_sources_path[i_str];
}

void parameters::addDustRayDetector(dlist & val)
{
    // Minimum wavelength (in SI)
    dust_ray_detectors.push_back(val[0]);
    // Maximum wavelength (in SI)
    dust_ray_detectors.push_back(val[1]);
    // Number of wavelengths (all)
    dust_ray_detectors.push_back(val[2]);
    // Source index (all)
    dust_ray_detectors.push_back(val[3] - 1);
    // Rot angle #1 (cart, polar, slice) / obs. position X (healpix)
    dust_ray_detectors.push_back(val[4]);
    // Rot angle #2 (cart, polar, slice) / obs. position Y (healpix)
    dust_ray_detectors.push_back(val[5]);
    // Distance from observer to model (cart, polar, slice) / obs. position Z
    // (healpix)
    dust_ray_detectors.push_back(val[6]);
    // Side length of dust detector in x-dir (cart, polar, slice) / l_max (healpix)
    dust_ray_detectors.push_back(val[7]);
    // Side length of dust detector in y-dir (cart, polar, slice) / l_min (healpix)
    dust_ray_detectors.push_back(val[8]);
    // delta_x (cart, slice) / None (polar) / b_max (healpix)
    dust_ray_detectors.push_back(val[9]);
    // delta_y (cart, slice) / None (polar) / b_min (healpix)
    dust_ray_detectors.push_back(val[10]);
    //bubble radius (heal)
    dust_ray_detectors.push_back(val[11]);
    // dust detector type/grid (all)
    dust_ray_detectors.push_back(val[12]);
    // number of pixel in x-dir. (cart, polar, slice) / N_side (healpix)
    dust_ray_detectors.push_back(val[13]);
    // number of pixel in y-direction (cart, polar, slice)
    dust_ray_detectors.push_back(val[14]);

    switch(uint(val[NR_OF_RAY_DET - 3]))
    {
        case DET_SPHER:
            if(rt_grid_description.find("healpix") == string::npos)
            {
                rt_grid_description += "healpix";
                if(dust_ray_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_POLAR:
            if(rt_grid_description.find("polar") == string::npos)
            {
                rt_grid_description += "polar";
                if(dust_ray_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_SLICE:
            if(rt_grid_description.find("slice") == string::npos)
            {
                rt_grid_description += "slice";
                if(dust_ray_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        default:
            if(rt_grid_description.find("cartesian") == string::npos)
            {
                rt_grid_description += "cartesian";
                if(dust_ray_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;
    }
}

void parameters::addOpiateRayDetector(dlist & val)
{
    // BG surce ID
    opiate_ray_detectors.push_back(val[0]);

    // Maximum velocity (SI)
    opiate_ray_detectors.push_back(val[1]);

    // ang1 or pos_x
    opiate_ray_detectors.push_back(val[2]);

    // ang2 or pos_y
    opiate_ray_detectors.push_back(val[3]);

    // dist or pos_z
    opiate_ray_detectors.push_back(val[4]);

    // sidelength_x or l_min
    opiate_ray_detectors.push_back(val[5]);

    // sidelength_y or l_max
    opiate_ray_detectors.push_back(val[6]);

    // off_x or b_min
    opiate_ray_detectors.push_back(val[7]);

    // off_y or b_max
    opiate_ray_detectors.push_back(val[8]);

    // empty or d_vx
    opiate_ray_detectors.push_back(val[9]);

    // empty or d_vy
    opiate_ray_detectors.push_back(val[10]);

    // empty or d_vz
    opiate_ray_detectors.push_back(val[11]);

    // empty or bubble radius
    opiate_ray_detectors.push_back(val[12]);

    // det. type
    opiate_ray_detectors.push_back(val[13]);

    // N_x or n_side
    opiate_ray_detectors.push_back(val[14]);

    // N_y or n_side
    opiate_ray_detectors.push_back(val[15]);

    // N_vel
    opiate_ray_detectors.push_back(val[16]);

    switch(uint(val[13]))
    {
        case DET_SPHER:
            if(rt_grid_description.find("healpix") == string::npos)
            {
                rt_grid_description += "healpix";
                if(opiate_ray_detectors.size() > NR_OF_OPIATE_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_POLAR:
            if(rt_grid_description.find("polar") == string::npos)
            {
                rt_grid_description += "polar";
                if(opiate_ray_detectors.size() > NR_OF_OPIATE_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_SLICE:
            if(rt_grid_description.find("slice") == string::npos)
            {
                rt_grid_description += "slice";
                if(opiate_ray_detectors.size() > NR_OF_OPIATE_DET)
                    rt_grid_description += ", ";
            }
            break;

        default:
            if(rt_grid_description.find("cartesian") == string::npos)
            {
                rt_grid_description += "cartesian";
                if(opiate_ray_detectors.size() > NR_OF_OPIATE_DET)
                    rt_grid_description += ", ";
            }
            break;
    }
}

void parameters::addSyncRayDetector(dlist & val)
{
    // Minimum wavelength (in SI)
    sync_ray_detectors.push_back(val[0]);
    // Maximum wavelength (in SI)
    sync_ray_detectors.push_back(val[1]);
    // Number of wavelengths (all)
    sync_ray_detectors.push_back(val[2]);
    // Source index (all)
    sync_ray_detectors.push_back(val[3] - 1);
    // Rot angle #1 (cart, polar, slice) / obs. position X (healpix)
    sync_ray_detectors.push_back(val[4]);
    // Rot angle #2 (cart, polar, slice) / obs. position Y (healpix)
    sync_ray_detectors.push_back(val[5]);
    // Distance from observer to model (cart, polar, slice) / obs. position Z
    // (healpix)
    sync_ray_detectors.push_back(val[6]);
    // Side length of dust detector in x-dir (cart, polar, slice) / l_max (healpix)
    sync_ray_detectors.push_back(val[7]);
    // Side length of dust detector in y-dir (cart, polar, slice) / l_min (healpix)
    sync_ray_detectors.push_back(val[8]);
    // delta_x (cart, slice) / None (polar) / b_max (healpix)
    sync_ray_detectors.push_back(val[9]);
    // delta_y (cart, slice) / None (polar) / b_min (healpix)
    sync_ray_detectors.push_back(val[10]);
    // bubble size (heal)
    sync_ray_detectors.push_back(val[11]);
    // dust detector type/grid (all)
    sync_ray_detectors.push_back(val[12]);
    // number of pixel in x-dir. (cart, polar, slice) / N_side (healpix)
    sync_ray_detectors.push_back(val[13]);
    // number of pixel in y-direction (cart, polar, slice)
    sync_ray_detectors.push_back(val[14]);

    switch(uint(val[NR_OF_RAY_DET - 3]))//
    {
        case DET_SPHER:
            if(rt_grid_description.find("healpix") == string::npos)
            {
                rt_grid_description += "healpix";
                if(sync_ray_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_POLAR:
            if(rt_grid_description.find("polar") == string::npos)
            {
                rt_grid_description += "polar";
                if(sync_ray_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_SLICE:
            if(rt_grid_description.find("slice") == string::npos)
            {
                rt_grid_description += "slice";
                if(sync_ray_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        default:
            if(rt_grid_description.find("cartesian") == string::npos)
            {
                rt_grid_description += "cartesian";
                if(sync_ray_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;
    }
}

void parameters::addFreeRayDetector(dlist & val)
{
    // Minimum wavelength (in SI)
    free_ray_detectors1.push_back(val[0]);
    // Maximum wavelength (in SI)
    free_ray_detectors1.push_back(val[1]);
    // Number of wavelengths (all)
    free_ray_detectors1.push_back(val[2]);
    // Source index (all)
    free_ray_detectors1.push_back(val[3] - 1);
    // Rot angle #1 (cart, polar, slice) / obs. position X (healpix)
    free_ray_detectors1.push_back(val[4]);
    // Rot angle #2 (cart, polar, slice) / obs. position Y (healpix)
    free_ray_detectors1.push_back(val[5]);
    // Distance from observer to model (cart, polar, slice) / obs. position Z
    // (healpix)
    free_ray_detectors1.push_back(val[6]);
    // Side length of dust detector in x-dir (cart, polar, slice) / l_max (healpix)
    free_ray_detectors1.push_back(val[7]);
    // Side length of dust detector in y-dir (cart, polar, slice) / l_min (healpix)
    free_ray_detectors1.push_back(val[8]);
    // delta_x (cart, slice) / None (polar) / b_max (healpix)
    free_ray_detectors1.push_back(val[9]);
    // delta_y (cart, slice) / None (polar) / b_min (healpix)
    free_ray_detectors1.push_back(val[10]);
    // bubble size (heal)
    free_ray_detectors1.push_back(val[11]);
    // dust detector type/grid (all)
    free_ray_detectors1.push_back(val[12]);
    // number of pixel in x-dir. (cart, polar, slice) / N_side (healpix)
    free_ray_detectors1.push_back(val[13]);
    // number of pixel in y-direction (cart, polar, slice)
    free_ray_detectors1.push_back(val[14]);

    switch(uint(val[NR_OF_RAY_DET - 3]))//
    {
        case DET_SPHER:
            if(rt_grid_description.find("healpix") == string::npos)
            {
                rt_grid_description += "healpix";
                if(free_ray_detectors1.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_POLAR:
            if(rt_grid_description.find("polar") == string::npos)
            {
                rt_grid_description += "polar";
                if(free_ray_detectors1.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_SLICE:
            if(rt_grid_description.find("slice") == string::npos)
            {
                rt_grid_description += "slice";
                if(free_ray_detectors1.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        default:
            if(rt_grid_description.find("cartesian") == string::npos)
            {
                rt_grid_description += "cartesian";
                if(free_ray_detectors1.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;
    }
}

void parameters::addDustAMEDetector(dlist & val)
{
    // Minimum wavelength (in SI)
    dust_ame_detectors.push_back(val[0]);
    // Maximum wavelength (in SI)
    dust_ame_detectors.push_back(val[1]);
    // Number of wavelengths (all)
    dust_ame_detectors.push_back(val[2]);
    // Source index (all)
    dust_ame_detectors.push_back(val[3] - 1);
    // Rot angle #1 (cart, polar, slice) / obs. position X (healpix)
    dust_ame_detectors.push_back(val[4]);
    // Rot angle #2 (cart, polar, slice) / obs. position Y (healpix)
    dust_ame_detectors.push_back(val[5]);
    // Distance from observer to model (cart, polar, slice) / obs. position Z
    // (healpix)
    dust_ame_detectors.push_back(val[6]);
    // Side length of dust detector in x-dir (cart, polar, slice) / l_max (healpix)
    dust_ame_detectors.push_back(val[7]);
    // Side length of dust detector in y-dir (cart, polar, slice) / l_min (healpix)
    dust_ame_detectors.push_back(val[8]);
    // delta_x (cart, slice) / None (polar) / b_max (healpix)
    dust_ame_detectors.push_back(val[9]);
    // delta_y (cart, slice) / None (polar) / b_min (healpix)
    dust_ame_detectors.push_back(val[10]);
    // bubble size (heal)
    dust_ame_detectors.push_back(val[11]);
    // dust detector type/grid (all)
    dust_ame_detectors.push_back(val[12]);
    // number of pixel in x-dir. (cart, polar, slice) / N_side (healpix)
    dust_ame_detectors.push_back(val[13]);
    // number of pixel in y-direction (cart, polar, slice)
    dust_ame_detectors.push_back(val[14]);

    switch(uint(val[NR_OF_RAY_DET - 3]))//
    {
        case DET_SPHER:
            if(rt_grid_description.find("healpix") == string::npos)
            {
                rt_grid_description += "healpix";
                if(dust_ame_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_POLAR:
            if(rt_grid_description.find("polar") == string::npos)
            {
                rt_grid_description += "polar";
                if(dust_ame_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_SLICE:
            if(rt_grid_description.find("slice") == string::npos)
            {
                rt_grid_description += "slice";
                if(dust_ame_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;

        default:
            if(rt_grid_description.find("cartesian") == string::npos)
            {
                rt_grid_description += "cartesian";
                if(dust_ame_detectors.size() > NR_OF_RAY_DET)
                    rt_grid_description += ", ";
            }
            break;
    }
}

/*void parameters::addLineOpiateDetector(dlist & val)
{
    // ang1
    line_opiate_detectors.push_back(val[0]);
    // ang2
    line_opiate_detectors.push_back(val[1]);
    // abundance
    line_opiate_detectors.push_back(val[2]);
    // mol. weight
    line_opiate_detectors.push_back(val[3]);
    // vel. max
    line_opiate_detectors.push_back(val[4]);
    // nr_of vels
    line_opiate_detectors.push_back(val[5]);
    // OPIATE column
    line_opiate_detectors.push_back(val[6]);
    // distance
    line_opiate_detectors.push_back(val[7]);
    // bins
    line_opiate_detectors.push_back(val[8]);
}

void parameters::setOpiateParamPath(string val)
{
    opiate_param_path = val;
}

void parameters::setOpiateDataPath(string val)
{
    opiate_data_path = val;
}*/

void parameters::addLineRayDetector(dlist & val)
{
    dlist tmp_list;
    // Gas species ID (all)
    uint i_species = uint(val[0]) - 1;

    if(line_ray_detector_list.count(i_species))
        tmp_list = line_ray_detector_list[i_species];

    // transition ID (all)
    tmp_list.push_back(val[1] - 1);
    // source ID (all)
    tmp_list.push_back(val[2] - 1);
    
    // min_velocity (all)
    tmp_list.push_back(val[3]);
    // max_velocity (all)
    tmp_list.push_back(val[4]);
    
    
    // Rot angle #1 (cart, polar) / obs. position X (healpix)
    tmp_list.push_back(val[5]);
    // Rot angle #2 (cart, polar) / obs. position Y (healpix)
    tmp_list.push_back(val[6]);
    // Distance from observer to model (cart, polar) / obs. position Z (healpix)
    tmp_list.push_back(val[7]);
    // Side length of dust detector in x-dir (cart, polar) / l_max (healpix)
    tmp_list.push_back(val[8]);
    // Side length of dust detector in y-dir (cart, polar) / l_min (healpix)
    tmp_list.push_back(val[9]);
    // delta_x (cart, slice) / None (polar) / b_max (healpix)
    tmp_list.push_back(val[10]);
    // delta_y (cart, slice) / None (polar) / b_min (healpix)
    tmp_list.push_back(val[11]);
    // None (cart, polar, slice) / obs. velocity X (healpix)
    tmp_list.push_back(val[12]);
    // None (cart, polar, slice) / obs. velocity Y (healpix)
    tmp_list.push_back(val[13]);
    // None (cart, polar, slice) / obs. velocity Z (healpix)
    tmp_list.push_back(val[14]);
    // dust detector type/grid (all)
    tmp_list.push_back(val[15]);
    // number of pixel in x-dir. (cart, polar, slice) / N_side (healpix)
    tmp_list.push_back(val[16]);
    // number of pixel in y-direction (all)
    tmp_list.push_back(val[17]);
    // number of velocity channels (all)
    tmp_list.push_back(val[18]);

    line_ray_detector_list[i_species] = tmp_list;

    switch(uint(val[NR_OF_LINE_DET - 3]))
    {
        case DET_SPHER:
            if(rt_grid_description.find("healpix") == string::npos)
            {
                rt_grid_description += "healpix";
                if(tmp_list.size() > NR_OF_LINE_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_POLAR:
            if(rt_grid_description.find("polar") == string::npos)
            {
                rt_grid_description += "polar";
                if(tmp_list.size() > NR_OF_LINE_DET)
                    rt_grid_description += ", ";
            }
            break;

        case DET_SLICE:
            if(rt_grid_description.find("slice") == string::npos)
            {
                rt_grid_description += "slice";
                if(tmp_list.size() > NR_OF_LINE_DET)
                    rt_grid_description += ", ";
            }
            break;

        default:
            if(rt_grid_description.find("cartesian") == string::npos)
            {
                rt_grid_description += "cartesian";
                if(tmp_list.size() > NR_OF_LINE_DET)
                    rt_grid_description += ", ";
            }
            break;
    }
}

void parameters::addDustMCDetector(dlist & val)
{
    // Minimum wavelength (in SI)
    dust_mc_detectors.push_back(val[0]);
    // Maximum wavelength (in SI)
    dust_mc_detectors.push_back(val[1]);
    // Number of wavelengths (all)
    dust_mc_detectors.push_back(val[2]);
    // Rot angle #1
    dust_mc_detectors.push_back(val[3]);
    // Rot angle #2
    dust_mc_detectors.push_back(val[4]);
    // Distance from observer to model
    dust_mc_detectors.push_back(val[5]);
    // Side length of dust detector in x-dir
    dust_mc_detectors.push_back(val[6]);
    // Side length of dust detector in y-dir
    dust_mc_detectors.push_back(val[7]);
    // Shift length of dust detector in x-dir
    dust_mc_detectors.push_back(val[8]);
    // Shift length of dust detector in y-dir
    dust_mc_detectors.push_back(val[9]);
    // number of pixel in x-direction
    dust_mc_detectors.push_back(val[10]);
    // number of pixel in y-direction
    dust_mc_detectors.push_back(val[11]);
}

uint parameters::getNrOfDustMCDetectors()
{
    return uint(dust_mc_detectors.size() / NR_OF_MC_DET);
}

dlist parameters::getDustMCDetectors()
{
    return dust_mc_detectors;
}

void parameters::addPointSource(dlist & val, string path)
{
    // Position X
    point_sources.push_back(val[0]);
    // Position Y
    point_sources.push_back(val[1]);
    // Position Z
    point_sources.push_back(val[2]);
    // RadiuDetector(s [m]
    point_sources.push_back(val[3]);
    // Effective temperature [K]
    point_sources.push_back(val[4]);
    // sublimation radius
    point_sources.push_back(val[5]);
    // Stokes q plarization
    point_sources.push_back(val[6]);
    // Stokes u plarization
    point_sources.push_back(val[7]);
    // Number of photons
    point_sources.push_back(val[8]);
    point_sources_str.push_back(path);
}

void parameters::addLaserSource(dlist & val)
{
    // Position X
    laser_sources.push_back(val[0]);
    // Position Y
    laser_sources.push_back(val[1]);
    // Position Z
    laser_sources.push_back(val[2]);
    // Direction X
    laser_sources.push_back(val[3]);
    // Direction Y
    laser_sources.push_back(val[4]);
    // Direction Z
    laser_sources.push_back(val[5]);
    // Total power [W]
    laser_sources.push_back(val[6]);
    // Central wavelength [m]
    laser_sources.push_back(val[7]);
    // FWHM of the laser emission [m]
    laser_sources.push_back(val[8]);
    // Stokes Q plarization
    laser_sources.push_back(val[9]);
    // Stokes U plarization
    laser_sources.push_back(val[10]);
    // Number of photons
    laser_sources.push_back(val[11]);
}

void parameters::addBackgroundSource(dlist & val)
{
    // Radiating effective black body surface [m^2]
    background_sources.push_back(val[0]);
    // Effective temperature [K]
    background_sources.push_back(val[1]);
    // Stokes Q plarization
    background_sources.push_back(val[2]);
    // Stokes U plarization
    background_sources.push_back(val[3]);
    // Stokes V plarization
    background_sources.push_back(val[4]);
    // Rot angle #1
    background_sources.push_back(val[5]);
    // Rot angle #2
    background_sources.push_back(val[6]);
    // Number of background source photons
    background_sources.push_back(val[7]);
    background_sources_path.push_back("");
}

void parameters::addBackgroundSource(string path, dlist & val)
{
    // Radiating effective black body surface [m^2]
    background_sources.push_back(0);
    // Effective temperature [K]
    background_sources.push_back(0);
    // Stokes Q plarization
    background_sources.push_back(0);
    // Stokes U plarization
    background_sources.push_back(0);
    // Stokes V plarization
    background_sources.push_back(0);
    // Rot angle #1
    background_sources.push_back(val[0]);
    // Rot angle #2
    background_sources.push_back(val[1]);
    // Number of dust detector pixel per axis
    background_sources.push_back(val[2]);
    background_sources_path.push_back(path);
}

void parameters::setGauntPath(string path)
{
    gaunt_path = path;
}
    
string parameters::getGauntPath()
{
    return gaunt_path;
}

void parameters::addBackgroundSource(string path)
{
    // Radiating effective black body surface [m^2]
    background_sources.push_back(0);
    // Effective temperature [K]
    background_sources.push_back(0);
    // Stokes Q plarization
    background_sources.push_back(0);
    // Stokes U plarization
    background_sources.push_back(0);
    // Stokes V plarization
    background_sources.push_back(0);
    // Rot angle #1
    background_sources.push_back(0);
    // Rot angle #2
    background_sources.push_back(0);
    // Number of dust detector pixel per axis
    background_sources.push_back(0);
    background_sources_path.push_back(path);
}

void parameters::resetDustFiles()
{
    if(!reset_dust_files)
    {
        dust_paths.clear();
        size_keywords.clear();
        dust_fractions.clear();
        material_density.clear();
        a_min_global.clear();
        a_max_global.clear();
        for(uint i = 0; i < NR_OF_SIZE_DIST_PARAM; i++)
            size_parameter_map[i].clear();
        reset_dust_files = true;
    }
}

void parameters::addGasSpecies(string gas_species_path, string zeeman_path, dlist & val)
{
    gas_species_cat_path.push_back(gas_species_path);
    zeeman_catalog_path.push_back(zeeman_path);
    gas_species_level_pop_type.push_back(uint(val[0]));
    gas_species_abundance.push_back(val[1]);
}

uint parameters::getNrOfGasSpecies()
{
    return uint(gas_species_abundance.size());
}

uint parameters::getNrOfSpectralLines(uint i_species)
{
    return uint(line_ray_detector_list[i_species].size() / NR_OF_LINE_DET);
}

int * parameters::getSpectralLines(uint i_species)
{
    uint nr_of_spectral_lines = getNrOfSpectralLines(i_species);
    int * spectral_lines = new int[nr_of_spectral_lines];
    dlist line_ray_detector = getLineRayDetector(i_species);
    for(uint i = 0; i < line_ray_detector_list[i_species].size(); i += NR_OF_LINE_DET)
    {
        uint pos = i / NR_OF_LINE_DET;

        spectral_lines[pos] = int(line_ray_detector[i]);
    }
    return spectral_lines;
}

uint parameters::getNrOfDustRayDetectors()
{
    return uint(dust_ray_detectors.size() / NR_OF_RAY_DET);
}

uint parameters::getNrOfFreeRayDetectors()
{
    return uint(free_ray_detectors1.size() / NR_OF_RAY_DET);
}

uint parameters::getNrOfDustAMERayDetectors()
{
    return uint(dust_ame_detectors.size() / NR_OF_RAY_DET);
}

void parameters::addSubStatus(int val)
{
    sub_status |= val;
}

void parameters::setHealpixOrientation(uint val)
{
    healpix_orientation = val;
}

uint parameters::getHealpixOrientation() const
{
    return healpix_orientation;
}

const dlist & parameters::getLineRayDetector(uint i_species) const
{
    return line_ray_detector_list.at(i_species);
}

const maplist & parameters::getLineRayDetectors() const
{
    return line_ray_detector_list;
}

void parameters::addNanoGrains(dlist & pr)
{
    // minimal grain size
    nano_grains.push_back(pr[0]);
    
    // maximal grain size
    nano_grains.push_back(pr[1]);
    
    // number of grains
    nano_grains.push_back(pr[2]);
    
    // distribution ref. grain size
    nano_grains.push_back(pr[3]);
    
    // distribution sigma
    nano_grains.push_back(pr[4]);
    
    // energy of work_function
    nano_grains.push_back(pr[5]);
    
    // mass of unit cell
    nano_grains.push_back(pr[6]);
    
    // beta of el. dipole moment
    nano_grains.push_back(pr[7]);
    
    // tensile strength S_max
    nano_grains.push_back(pr[8]);
    
    // dust choice ID
    nano_grains.push_back(pr[9]);
}

uint parameters::getNrOfNanoGrains()
{
    return uint(nano_grains.size() / NR_OF_NANO_GRAIN_PRAM);
}

dlist parameters::findNanoGrainList(uint id)
{
    dlist result;
    
    for(int i=0; i<nano_grains.size(); i+=NR_OF_NANO_GRAIN_PRAM)
    {
        uint tmp_id = nano_grains[i+9];
        
        if(tmp_id==id)
        {
            int end = min(int(nano_grains.size()), i + 10);
            result = dlist(nano_grains.begin() + i, nano_grains.begin() + end);
            return result;
        }
            
    }
    
    return result;
}

void parameters::addDustComponent(string path,
                        string size_key,
                        double fr,
                        double mat_dens,
                        double a_min,
                        double a_max,
                        dlist size_parameter)
{
    dust_paths.push_back(path);
    size_keywords.push_back(size_key);
    dust_fractions.push_back(fr);
    material_density.push_back(mat_dens);
    a_min_global.push_back(a_min);
    a_max_global.push_back(a_max);
    for(uint i = 0; i < NR_OF_SIZE_DIST_PARAM; i++)
        size_parameter_map[i].push_back(size_parameter[i]);
}

bool parameters::getVelMapsFits() const
{
    return vel_maps_fits;
}

uint parameters::getHealType() const
{
    return heal_type;
}

bool  parameters::projectHealMaps() const
{
    return (param_projection.size() == 2);
}

llist parameters::getProjectedFitsParam() const
{
    return param_projection;
}

long parameters::getProjectedFitsX() const
{
    return param_projection[0];
}

long parameters::getProjectedFitsY() const
{
    return param_projection[1];
}

/*uint parameters::getDustChoiceFromComponentId(uint i) const
{
    return component_id_to_choice[i];
}*/

uint parameters::getDustChoiceFromMixtureId(uint i) const
{
    return dust_choices[i];
}

void parameters::printRTGridDescription()
{
    cout << "- background grid shape : " << rt_grid_description << endl;
}

double parameters::getAcceptanceAngle() const
{
    return acceptance_angle;
}

string parameters::getDustPath(uint i) const
{
    return dust_paths[i];
}

int parameters::getSubStatus()
{
    return sub_status;
}

double parameters::getDustFraction(uint i) const
{
    return dust_fractions[i];
}

string parameters::getDustSizeKeyword(uint i) const
{
    return size_keywords[i];
}

dlist parameters::getDustSizeParameter(uint i_comp) const
{
    dlist res;
    for(uint i = 0; i < NR_OF_SIZE_DIST_PARAM; i++)
        res.push_back(size_parameter_map.at(i).at(i_comp));
    return res;
}

string parameters::getGasSpeciesCatalogPath(uint i_species) const
{
    return gas_species_cat_path[i_species];
}

double parameters::getGasSpeciesAbundance(uint i_species) const
{
    return gas_species_abundance[i_species];
}

uint parameters::getGasSpeciesLevelPopType(uint i_species) const
{
    return gas_species_level_pop_type[i_species];
}

bool parameters::isGasSpeciesLevelPopMC() const
{
    for(uint i_species = 0; i_species < gas_species_level_pop_type.size(); i_species++)
        if(gas_species_level_pop_type[i_species] == POP_MC)
            return true;
    return false;
}

uint parameters::getMaxSubpixelLvl() const
{
    return max_subpixel_lvl;
}

uint parameters::getTotalNrOfDustComponents() const
{
    return (uint)dust_paths.size();
}

void parameters::resetNrOfDustComponents()
{
    return dust_paths.clear();
}

void parameters::addDiffuseSource(dlist & val, string path)
{
    // Position X
    diffuse_sources.push_back(val[0]);
    // Position Y
    diffuse_sources.push_back(val[1]);
    // Position Z
    diffuse_sources.push_back(val[2]);
    // Radius [m]
    diffuse_sources.push_back(val[3]);
    // Effective temperature [K]
    diffuse_sources.push_back(val[4]);
    
    // Extent of the starfield (sigma_x) [m]
    diffuse_sources.push_back(val[5]);
    // Extent of the starfield (sigma_y) [m]
    diffuse_sources.push_back(val[6]);
    // Extent of the starfield (sigma_z) [m]
    diffuse_sources.push_back(val[7]);
    
    // Extent of the ellipse (a) [m]
    diffuse_sources.push_back(val[8]);
    // Extent of the ellipse (b) [m]
    diffuse_sources.push_back(val[9]);
    // Extent of the ellipse (c) [m]
    diffuse_sources.push_back(val[10]);
    
    // rot. axis rx1 [m]
    diffuse_sources.push_back(val[11]);
    // rot. axis ry1 [m]
    diffuse_sources.push_back(val[12]);
    // rot. axis rz1 [m]
    diffuse_sources.push_back(val[13]);
    // rot. angle ang1 [deg]
    diffuse_sources.push_back(val[14]);
    
    // rot. axis rx2 [m]
    diffuse_sources.push_back(val[15]);
    // rot. axis ry2 [m]
    diffuse_sources.push_back(val[16]);
    // rot. axis rz2 [m]
    diffuse_sources.push_back(val[17]);
    // rot. angle ang2 [deg]
    diffuse_sources.push_back(val[18]);    
    
    // Stokes q plarization
    diffuse_sources.push_back(val[19]);
    // Stokes u plarization
    diffuse_sources.push_back(val[20]);
    
    // Number of photons
    diffuse_sources.push_back(val[21]);
    diffuse_sources_str.push_back(path);
}

