#include "Typedefs.h"
#include "Vector.h"

#ifndef PARAMETERS
#define PARAMETERS

class parameters
{
  public:
    parameters()
    {
        path_grid = "";
        path_output = "";
        phID = PH_HG;
        conv_l_in_SI = 1;
        conv_dH_in_SI = 1;
        conv_B_in_SI = 1;
        conv_V_in_SI = 1;
        conv_mass_fraction = 0.01;
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
        larm_f = 4.1e-21;

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
        vel_maps = false;
        dust_offset = false;
        dust_gas_coupling = false;
        full_dust_temp = false;
        save_radiation_field = false;
        scattering_to_raytracing = true;
        split_dust_emision = false;
        sublimate = false;
        individual_dust_fractions = false;

        nr_ofISRFPhotons = 0;
        nr_ofDustPhotons = 0;

        nrOfPlotPoints = 0;
        nrOfPlotVectors = 0;
        maxPlotLines = 0;
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
        nr_ofInpMidDataPoints = 0;
        nr_ofOutMidDataPoints = 0;

        f_highJ = 0.25;
        Q_ref = 0.4;
        alpha_Q = 3.0;

        f_cor = 0.6;
        adjTgas = 0;
        isrf_g_zero = 0;
        isrf_radius = 0;

        isrf_path = "";

        max_subpixel_lvl = 1;
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

    ~parameters()
    {}

    string getOpiatePathEmission()
    {
        return opiata_path_emi;
    }

    string getOpiatePathAbsorption()
    {
        return opiata_path_abs;
    }

    string getOpiateSpec(uint pos)
    {
        return opiate_spec_ids[pos];
    }

    uint getNrOfOPIATESpecies()
    {
        return uint(opiate_spec_ids.size());
    }

    const Vector3D & getAxis1() const
    {
        return axis1;
    }

    const Vector3D & getAxis2() const
    {
        return axis2;
    }

    uint getOutAMIRAPoints() const
    {
        return nr_ofOutAMIRAPoints;
    }

    uint getInpAMIRAPoints() const
    {
        return nr_ofInpAMIRAPoints;
    }

    bool plotInpMidPoints() const
    {
        return plot_inp_points;
    }

    bool plotOutMidPoints() const
    {
        return plot_out_points;
    }

    dlist getMidplane3dParams() const
    {
        return midplane_3d_param;
    }

    const uilist & getPlotList() const
    {
        return plot_list;
    }

    bool isInPlotList(uint id)
    {
        if(plot_list.empty())
            return true;

        return (find(plot_list.begin(), plot_list.end(), id) != plot_list.end());
    }

    uint getInpMidDataPoints() const
    {
        return nr_ofInpMidDataPoints;
    }

    uint getOutMidDataPoints() const
    {
        return nr_ofOutMidDataPoints;
    }

    /*const dlist & getOpiateSequence() const
    {
        return line_opiate_detectors;
    }*/

    uint getMidplaneZoom() const
    {
        return midplane_zoom;
    }

    int getCommand() const
    {
        return cmd;
    }

    bool isRatSimulation() const
    {
        if(getCommand() == CMD_RAT || getCommand() == CMD_TEMP_RAT)
            return true;
        return false;
    }

    bool isMonteCarloSimulation() const
    {
        if(getCommand() == CMD_TEMP || getCommand() == CMD_TEMP_RAT || getCommand() == CMD_RAT)
            return true;
        return false;
    }

    bool isRaytracingSimulation() const
    {
        if(getCommand() == CMD_OPIATE || getCommand() == CMD_DUST_EMISSION ||
           getCommand() == CMD_SYNCHROTRON || getCommand() == CMD_LINE_EMISSION)
            return true;
        return false;
    }

    bool isTemperatureSimulation() const
    {
        if(getCommand() == CMD_TEMP || getCommand() == CMD_TEMP_RAT)
            return true;
        return false;
    }

    double getStarMass(uint i) const
    {
        return star_mass[i];
    }

    string getPathGrid() const
    {
        return path_grid;
    }

    string getPathOutput() const
    {
        return path_output;
    }

    string getPathInput() const
    {
        return path_input;
    }

    uint getMinDetectorPixelX() const
    {
        return min_detector_pixel_x;
    }

    uint getMaxDetectorPixelX() const
    {
        return max_detector_pixel_x;
    }

    uint getMinDetectorPixelY() const
    {
        return min_detector_pixel_y;
    }

    uint getMaxDetectorPixelY() const
    {
        return max_detector_pixel_y;
    }

    double getMinDetectorAngle1() const
    {
        return min_rot_angle_1;
    }

    double getMaxDetectorAngle1() const
    {
        return max_rot_angle_1;
    }

    double getMinDetectorAngle2() const
    {
        return min_rot_angle_2;
    }

    double getMaxDetectorAngle2() const
    {
        return max_rot_angle_2;
    }

    double getMinSidelengthX() const
    {
        return min_sidelength_x;
    }

    double getMaxSidelengthX() const
    {
        return max_sidelength_x;
    }

    double getMinSidelengthY() const
    {
        return min_sidelength_y;
    }

    double getMaxSidelengthY() const
    {
        return max_sidelength_y;
    }

    bool getUseGridSidelengthX() const
    {
        return use_grid_sidelength_x;
    }

    bool getUseGridSidelengthY() const
    {
        return use_grid_sidelength_y;
    }

    double getMinMapShiftX() const
    {
        return min_ray_map_shift_x;
    }

    double getMinMapShiftY() const
    {
        return min_ray_map_shift_y;
    }

    double getMaxMapShiftX() const
    {
        return max_ray_map_shift_x;
    }

    double getMaxMapShiftY() const
    {
        return max_ray_map_shift_y;
    }

    double getSIConvLength() const
    {
        return conv_l_in_SI;
    }

    double getSIConvDH() const
    {
        return conv_dH_in_SI;
    }

    double getDelta0() const
    {
        return delta0;
    }

    double getLarmF() const
    {
        return larm_f;
    }

    double getSIConvBField() const
    {
        return conv_B_in_SI;
    }

    double getSIConvVField() const
    {
        return conv_V_in_SI;
    }

    bool getDustOffset() const
    {
        return dust_offset;
    }

    bool getDustGasCoupling() const
    {
        return dust_gas_coupling;
    }

    double getOffsetMinGasDensity() const
    {
        return offset_min_gas_dens;
    }

    bool getDustTempMulti() const
    {
        return full_dust_temp;
    }

    double getSizeMin(uint i) const
    {
        return a_min_global[i];
    }

    double getSizeMax(uint i) const
    {
        return a_max_global[i];
    }

    double getMaterialDensity(uint i) const
    {
        return material_density[i];
    }

    bool getDustSource() const
    {
        return nr_ofDustPhotons > 0;
    }

    bool getISRFSource() const
    {
        return nr_ofISRFPhotons > 0;
    }

    ullong getNrOfDustPhotons() const
    {
        return nr_ofDustPhotons;
    }

    double getDustMassFraction() const
    {
        return conv_mass_fraction;
    }

    uint getAlign() const
    {
        return align;
    }

    bool getAligRANDOM() const
    {
        return align == 0;
    }

    bool getAligPA() const
    {
        return (align & ALIG_PA) == ALIG_PA;
    }

    bool getAligIDG() const
    {
        return (align & ALIG_IDG) == ALIG_IDG;
    }

    bool getAligRAT() const
    {
        return (align & ALIG_RAT) == ALIG_RAT;
    }

    bool getAligGOLD() const
    {
        return (align & ALIG_GOLD) == ALIG_GOLD;
    }

    bool getAligINTERNAL() const
    {
        return (align & ALIG_INTERNAL) == ALIG_INTERNAL;
    }

    double getMu() const
    {
        return mu;
    }

    bool getMRW() const
    {
        return b_mrw;
    }

    bool getPDA() const
    {
        return b_pda;
    }

    bool getEnfScattering() const
    {
        return b_enforced;
    }

    double getStochasticHeatingMaxSize() const
    {
        return stochastic_heating_max_size;
    }

    bool getSaveRadiationField() const
    {
        return save_radiation_field;
    }

    bool getScatteringToRay() const
    {
        return scattering_to_raytracing;
    }

    bool splitDustEmission() const
    {
        return split_dust_emision;
    }

    bool getIndividualDustMassFractions() const
    {
        return individual_dust_fractions;
    }

    bool getIsSpeedOfSound() const
    {
        return is_speed_of_sound;
    }

    bool getPeelOff() const
    {
        return peel_off;
    }

    double getForegroundExtinctionMagnitude() const
    {
        return extinction_magnitude;
    }

    double getForegroundExtinctionWavelength() const
    {
        return extinction_magnitude_wavelength;
    }

    uint getForegroundExtinctionDustMixture() const
    {
        return extinction_i_mixture;
    }

    bool getVelFieldType() const
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

    uint getWriteRadiationField() const
    {
        return write_radiation_field;
    }

    bool getWriteGZero() const
    {
        return write_g_zero;
    }

    double getISRFGZero() const
    {
        return isrf_g_zero;
    }

    double getISRFRadius() const
    {
        return isrf_radius;
    }

    string getISRFPath() const
    {
        return isrf_path;
    }

    string getZeemanCatalog(uint i_species) const
    {
        return zeeman_catalog_path[i_species];
    }

    uint getAlignmentMechanism() const
    {
        return align;
    }

    double getMinObserverDistance() const
    {
        return min_obs_distance;
    }

    double getMaxObserverDistance() const
    {
        return max_obs_distance;
    }

    double getKeplerStarMass() const
    {
        return kepler_star_mass;
    }

    double getTurbulentVelocity() const
    {
        return turbulent_velocity;
    }

    uint getMCLvlPopNrOfPhotons() const
    {
        return nr_of_mc_lvl_pop_photons;
    }

    uint getMCLvlPopSeed() const
    {
        return mc_lvl_pop_seed;
    }

    uint getTaskID() const
    {
        return task_id;
    }

    uint getNrOfThreads() const
    {
        return nr_ofThreads;
    }

    ullong getNrOfISRFPhotons() const
    {
        return nr_ofISRFPhotons;
    }

    uint getNrOfMixtures() const
    {
        if(!dust_choices.empty())
            return dust_choices.size();
        else
            return 1;
    }

    uilist getDustComponentChoices() const
    {
        return dust_choices;
    }

    uint getPhaseFunctionID() const
    {
        return phID;
    }

    double getFHighJ() const
    {
        return f_highJ;
    }

    double getQref() const
    {
        return Q_ref;
    }

    double getAlphaQ() const
    {
        return alpha_Q;
    }

    double getFcorr() const
    {
        return f_cor;
    }

    double getAdjTgas() const
    {
        return adjTgas;
    }

    uint getNrOfDiffuseSources() const
    {
        return uint(diffuse_sources.size() / NR_OF_DIFF_SOURCES);
    }

    uint getNrOfPointSources() const
    {
        return uint(point_sources.size() / NR_OF_POINT_SOURCES);
    }

    uint getNrOfLaserSources() const
    {
        return uint(laser_sources.size() / NR_OF_LASER_SOURCES);
    }

    uint getNrOfBackgroundSources() const
    {
        return uint(background_sources.size() / NR_OF_BG_SOURCES);
    }

    double getXYMin() const
    {
        return xymin;
    }

    double getXYMax() const
    {
        return xymax;
    }

    double getXYSteps() const
    {
        return xysteps;
    }

    uint getXYBins() const
    {
        return xy_bins;
    }

    string getXYLabel() const
    {
        return xylabel;
    }

    bool isAutoScale() const
    {
        return autoscale;
    }

    uint getNrOfSources() const
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

    uint getNrOfPlotPoints() const
    {
        return nrOfPlotPoints;
    }

    uint getNrOfPlotVectors() const
    {
        return nrOfPlotVectors;
    }

    uint getMaxPlotLines() const
    {
        return maxPlotLines;
    }

    uint getStart() const
    {
        return start;
    }

    uint getStop() const
    {
        return stop;
    }


    void setOpiatePathEmission(string str)
    {
        opiata_path_emi=str;
    }

    void setOpiatePathAbsorption(string str)
    {
        opiata_path_abs=str;
    }


    void setXYMin(double val)
    {
        xymin = val;
    }

    void setXYMax(double val)
    {
        xymax = val;
    }

    void setXYSteps(double val)
    {
        xysteps = val;
    }

    void setXYBins(uint val)
    {
        xy_bins = val;
    }

    void setXYLabel(string val)
    {
        xylabel = val;
    }

    void setAutoScale(bool val)
    {
        autoscale = val;
    }

    void setAxis1(double x, double y, double z)
    {
        axis1.set(x, y, z);
        axis1.normalize();
    }

    void setAxis2(double x, double y, double z)
    {
        axis2.set(x, y, z);
        axis2.normalize();
    }

    void setStart(uint val)
    {
        start = val;
    }

    void setStop(uint val)
    {
        stop = val;
    }

    void setMu(double val)
    {
        mu = val;
    }

    void setNrOfPlotPoints(uint val)
    {
        nrOfPlotPoints = val;
    }

    void setnrOfPlotVectors(uint val)
    {
        nrOfPlotVectors = val;
    }

    void setMaxPlotLines(uint val)
    {
        maxPlotLines = val;
    }

    void setNrOfThreads(uint val)
    {
        nr_ofThreads = val;
    }

    void setCommand(int val)
    {
        cmd = val;
    }

    void setISRF(string path, double g_zero = 0, double radius = 2)
    {
        isrf_path = path;
        isrf_g_zero = g_zero;
        isrf_radius = radius;
    }

    void setNrOfISRFPhotons(long val)
    {
        nr_ofISRFPhotons = val;
    }

    void addOpiateSpec(string str)
    {
        opiate_spec_ids.push_back(str);
    }

    void AddDustComponentChoice(uint dust_component_choice)
    {
        // If the dust component choice of a dust component is already loaded,
        // add the "dust_component_choice" to the list of component_id_to_choice but not
        // to the dust_choices.
        if(!component_id_to_choice.empty())
            for(uint a = 0; a < component_id_to_choice.size(); a++)
                if(component_id_to_choice[a] == dust_component_choice)
                {
                    component_id_to_choice.push_back(dust_component_choice);
                    return;
                }

        // If the "dust_component_choice" is not yet used, add it to the amount of
        // available choices.
        component_id_to_choice.push_back(dust_component_choice);
        dust_choices.push_back(dust_component_choice);

        // Sort dust choices
        sort(dust_choices.begin(), dust_choices.end());

        // Update the highest value of the dust choice ids.
        if(dust_component_choice > max_dust_component_choice)
            max_dust_component_choice = dust_component_choice;
    }

    uint getMaxDustComponentChoice()
    {
        return max_dust_component_choice;
    }

    void setTaskID(uint val)
    {
        task_id = val;
    }

    void addStarMass(double val)
    {
        star_mass.push_back(val);
    }

    void setStarMass(dlist val)
    {
        star_mass = val;
    }

    void setPathGrid(string val)
    {
        path_grid = val;
    }

    void setPathInput(string val)
    {
        path_input = val;
    }

    void setPathOutput(string val)
    {
        path_output = val;
    }

    void setDustOffset(bool val)
    {
        dust_offset = val;
    }

    void setDustOffset(double _offset_min_gas_dens)
    {
        dust_offset = true;
        offset_min_gas_dens = _offset_min_gas_dens;
    }

    void setDustGasCoupling(bool val)
    {
        dust_gas_coupling = val;
    }

    void setDustGasCoupling(double _offset_min_gas_dens)
    {
        dust_gas_coupling = true;
        offset_min_gas_dens = _offset_min_gas_dens;
    }

    void setFullDustTemp(bool val)
    {
        full_dust_temp = val;
    }

    void setNrOfDustPhotons(long val)
    {
        nr_ofDustPhotons = val;
    }

    void setDelta0(double val)
    {
        delta0 = val;
    }

    void addToPlotList(uint id)
    {
        plot_list.push_back(id);
    }

    void setLarmF(double val)
    {
        larm_f = val;
    }

    void setMRW(bool val)
    {
        b_mrw = val;
    }

    void setPDA(bool val)
    {
        b_pda = val;
    }

    void setEnfScattering(bool val)
    {
        b_enforced = val;
    }

    void setIsSpeedOfSound(bool val)
    {
        is_speed_of_sound = val;
    }

    void setPeelOff(bool val)
    {
        peel_off = val;
    }

    void setForegroundExtinction(double _extinction_magnitude,
                                 double _extinction_magnitude_wavelength = 0.55e-6,
                                 uint _extinction_i_mixture = MAX_UINT)
    {
        extinction_magnitude = _extinction_magnitude;
        extinction_magnitude_wavelength = _extinction_magnitude_wavelength;
        extinction_i_mixture = _extinction_i_mixture;
    }

    void setInpAMIRAPoints(uint val)
    {
        nr_ofInpAMIRAPoints = val;
    }

    void setOutAMIRAPoints(uint val)
    {
        nr_ofOutAMIRAPoints = val;
    }

    void setInpMidPlot(bool val)
    {
        plot_inp_points = val;
    }

    void setOutMidPlot(bool val)
    {
        plot_out_points = val;
    }

    void set3dMidplane(uint plane, uint nr_of_slices, double z_min, double z_max)
    {
        if(midplane_3d_param.size() == 0)
        {
            midplane_3d_param.push_back(plane);
            midplane_3d_param.push_back(nr_of_slices);
            midplane_3d_param.push_back(z_min);
            midplane_3d_param.push_back(z_max);
        }
        else
            cout << "\nHINT: Multiple 3D midplane plot commands found. The first command "
                    "will be considered!";
    }

    void setInpMidDataPoints(uint val)
    {
        nr_ofInpMidDataPoints = val;
    }

    void setOutMidDataPoints(uint val)
    {
        nr_ofOutMidDataPoints = val;
    }

    void setMidplaneZoom(uint val)
    {
        midplane_zoom = val;
    }

    void setWriteRadiationField(uint val)
    {
        write_radiation_field = val;
    }

    void setWriteGZero(double val)
    {
        write_g_zero = val;
    }

    void updateDetectorPixel(uint pixel_x, uint pixel_y)
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

    void updateDetectorAngles(double rot_angle_1, double rot_angle_2)
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

    void updateMapSidelength(double sidelength_x, double sidelength_y)
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

    void updateRayGridShift(double map_shift_x, double map_shift_y)
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

    void setStochasticHeatingMaxSize(double val)
    {
        stochastic_heating_max_size = val;
    }

    void setSaveRadiationField(bool val)
    {
        save_radiation_field = val;
    }

    void setScatteringToRay(bool val)
    {
        scattering_to_raytracing = val;
    }

    void setSplitDustEmission(bool val)
    {
        split_dust_emision = val;
    }

    void setSIConvLength(double val)
    {
        conv_l_in_SI = val;
    }

    void setSIConvDH(double val)
    {
        conv_dH_in_SI = val;
    }

    void setSIConvBField(double val)
    {
        conv_B_in_SI = val;
    }

    void setSIConvVField(double val)
    {
        conv_V_in_SI = val;
    }

    void updateSIConvLength(double val)
    {
        if(conv_l_in_SI != 1 && val != 1)
            cout << "\nHINT: <conv_len> may not be used multiple times!" << endl
                 << "      -> No problem if <path_grid_cgs> was used!" << endl;
        conv_l_in_SI *= val;
    }

    void updateSIConvDH(double val)
    {
        if(conv_dH_in_SI != 1 && val != 1)
            cout << "\nHINT: <conv_dens> may not be used multiple times!" << endl
                 << "      -> No problem if <path_grid_cgs> was used!" << endl;
        conv_dH_in_SI *= val;
    }

    void updateSIConvBField(double val)
    {
        if(conv_B_in_SI != 1 && val != 1)
            cout << "\nHINT: <conv_mag> may not be used multiple times!" << endl
                 << "      -> No problem if <path_grid_cgs> was used!" << endl;
        conv_B_in_SI *= val;
    }

    void updateSIConvVField(double val)
    {
        if(conv_V_in_SI != 1 && val != 1)
            cout << "\nHINT: <conv_vel> may not be used multiple times!" << endl
                 << "      -> No problem if <path_grid_cgs> was used!" << endl;
        conv_V_in_SI *= val;
    }

    void setDustMassFraction(double val)
    {
        conv_mass_fraction = val;
    }

    void setIndividualDustMassFractions(bool val)
    {
        individual_dust_fractions = val;
    }

    void addAlignmentMechanism(uint val)
    {
        align |= val;
    }

    void updateObserverDistance(double val)
    {
        if(min_obs_distance > val)
            min_obs_distance = val;

        if(max_obs_distance < val)
            max_obs_distance = val;
    }

    void setKeplerStarMass(double val)
    {
        kepler_star_mass = val;
    }

    void setMCLvlPopNrOfPhotons(uint val)
    {
        nr_of_mc_lvl_pop_photons = val;
    }

    void setMCLvlPopSeed(uint val)
    {
        if(val > 0)
            mc_lvl_pop_seed = val;
    }

    void setTurbulentVelocity(double val)
    {
        turbulent_velocity = val;
    }

    void addZeemanCatalog(string val)
    {
        zeeman_catalog_path.push_back(val);
    }

    void setPhaseFunctionID(uint val)
    {
        phID = val;
    }

    void setFhighJ(double val)
    {
        f_highJ = val;
    }

    void setQref(double val)
    {
        Q_ref = val;
    }

    void setAlphaQ(double val)
    {
        alpha_Q = val;
    }

    void setFcorr(double val)
    {
        f_cor = val;
    }

    void setAdjTgas(double val)
    {
        adjTgas = val;
    }

    void setVelMaps(bool val)
    {
        vel_maps = val;
    }

    void setAcceptanceAngle(double angle)
    {
        acceptance_angle = angle;
    }

    void setMaxSubpixelLvl(int val)
    {
        max_subpixel_lvl = val;
    }

    dlist & getDustRayDetectors()
    {
        return dust_ray_detectors;
    }

    dlist & getSyncRayDetectors()
    {
        return sync_ray_detectors;
    }

    dlist & getOPIATERayDetectors()
    {
        return opiate_ray_detectors;
    }

    dlist & getPointSources()
    {
        return point_sources;
    }

    dlist & getLaserSources()
    {
        return laser_sources;
    }

    strlist & getPointSourceStringList()
    {
        return point_sources_str;
    }

    string & getPointSourceString(uint i_str)
    {
        return point_sources_str[i_str];
    }

    strlist & getDiffuseSourceStringList()
    {
        return diffuse_sources_str;
    }

    string & getDiffuseSourceString(uint i_str)
    {
        return diffuse_sources_str[i_str];
    }

    dlist & getDiffuseSources()
    {
        return diffuse_sources;
    }

    dlist & getBackgroundSources()
    {
        return background_sources;
    }

    strlist & getBackgroundSourceStringList()
    {
        return background_sources_path;
    }

    string & getBackgroundSourceString(uint i_str)
    {
        return background_sources_path[i_str];
    }

    void addDustRayDetector(dlist & val)
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

    void addOpiateRayDetector(dlist & val)
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

    void addSyncRayDetector(dlist & val)
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

    /*void addLineOpiateDetector(dlist & val)
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

    void setOpiateParamPath(string val)
    {
        opiate_param_path = val;
    }

    void setOpiateDataPath(string val)
    {
        opiate_data_path = val;
    }*/

    void addLineRayDetector(dlist & val)
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
        // max_velocity (all)
        tmp_list.push_back(val[3]);
        // Rot angle #1 (cart, polar) / obs. position X (healpix)
        tmp_list.push_back(val[4]);
        // Rot angle #2 (cart, polar) / obs. position Y (healpix)
        tmp_list.push_back(val[5]);
        // Distance from observer to model (cart, polar) / obs. position Z (healpix)
        tmp_list.push_back(val[6]);
        // Side length of dust detector in x-dir (cart, polar) / l_max (healpix)
        tmp_list.push_back(val[7]);
        // Side length of dust detector in y-dir (cart, polar) / l_min (healpix)
        tmp_list.push_back(val[8]);
        // delta_x (cart, slice) / None (polar) / b_max (healpix)
        tmp_list.push_back(val[9]);
        // delta_y (cart, slice) / None (polar) / b_min (healpix)
        tmp_list.push_back(val[10]);
        // None (cart, polar, slice) / obs. velocity X (healpix)
        tmp_list.push_back(val[11]);
        // None (cart, polar, slice) / obs. velocity Y (healpix)
        tmp_list.push_back(val[12]);
        // None (cart, polar, slice) / obs. velocity Z (healpix)
        tmp_list.push_back(val[13]);
        // dust detector type/grid (all)
        tmp_list.push_back(val[14]);
        // number of pixel in x-dir. (cart, polar, slice) / N_side (healpix)
        tmp_list.push_back(val[15]);
        // number of pixel in y-direction (all)
        tmp_list.push_back(val[16]);
        // number of velocity channels (all)
        tmp_list.push_back(val[17]);

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

    void addDustMCDetector(dlist & val)
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

    uint getNrOfDustMCDetectors()
    {
        return uint(dust_mc_detectors.size() / NR_OF_MC_DET);
    }

    dlist getDustMCDetectors()
    {
        return dust_mc_detectors;
    }

    void addPointSource(dlist & val, string path)
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
        // Stokes Q plarization
        point_sources.push_back(val[5]);
        // Stokes U plarization
        point_sources.push_back(val[6]);
        // Number of photons
        point_sources.push_back(val[7]);
        point_sources_str.push_back(path);
    }

    void addLaserSource(dlist & val)
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

    void addBackgroundSource(dlist & val)
    {
        // RadiatDetector(ing effective black body surface [m^2]
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

    void addBackgroundSource(string path, dlist & val)
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

    void addBackgroundSource(string path)
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

    void resetDustFiles()
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

    void addGasSpecies(string gas_species_path, string zeeman_path, dlist & val)
    {
        gas_species_cat_path.push_back(gas_species_path);
        zeeman_catalog_path.push_back(zeeman_path);
        gas_species_level_pop_type.push_back(uint(val[0]));
        gas_species_abundance.push_back(val[1]);
    }

    uint getNrOfGasSpecies()
    {
        return uint(gas_species_abundance.size());
    }

    uint getNrOfSpectralLines(uint i_species)
    {
        return uint(line_ray_detector_list[i_species].size() / NR_OF_LINE_DET);
    }

    int * getSpectralLines(uint i_species)
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

    uint getNrOfDustRayDetectors()
    {
        return uint(dust_ray_detectors.size() / NR_OF_RAY_DET);
    }

    uint getNrOfSyncRayDetectors()
    {
        return uint(sync_ray_detectors.size() / NR_OF_RAY_DET);
    }

    void setSublimate(bool val)
    {
        sublimate = val;
    }

    void setHealpixOrientation(uint val)
    {
        healpix_orientation = val;
    }

    uint getHealpixOrientation() const
    {
        return healpix_orientation;
    }

    const dlist & getLineRayDetector(uint i_species) const
    {
        return line_ray_detector_list.at(i_species);
    }

    const maplist & getLineRayDetectors() const
    {
        return line_ray_detector_list;
    }

    void addDustComponent(string path,
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

    bool getVelMaps() const
    {
        return vel_maps;
    }

    uint getDustChoiceFromComponentId(uint i) const
    {
        return component_id_to_choice[i];
    }

    uint getDustChoiceFromMixtureId(uint i) const
    {
        return dust_choices[i];
    }

    void printRTGridDescription()
    {
        cout << "- background grid shape : " << rt_grid_description << endl;
    }

    double getAcceptanceAngle() const
    {
        return acceptance_angle;
    }

    string getDustPath(uint i) const
    {
        return dust_paths[i];
    }

    bool isSublimate()
    {
        return sublimate;
    }

    double getDustFraction(uint i) const
    {
        return dust_fractions[i];
    }

    string getDustSizeKeyword(uint i) const
    {
        return size_keywords[i];
    }

    dlist getDustSizeParameter(uint i_comp) const
    {
        dlist res;
        for(uint i = 0; i < NR_OF_SIZE_DIST_PARAM; i++)
            res.push_back(size_parameter_map.at(i).at(i_comp));
        return res;
    }

    string getGasSpeciesCatalogPath(uint i_species) const
    {
        return gas_species_cat_path[i_species];
    }

    double getGasSpeciesAbundance(uint i_species) const
    {
        return gas_species_abundance[i_species];
    }

    uint getGasSpeciesLevelPopType(uint i_species) const
    {
        return gas_species_level_pop_type[i_species];
    }

    bool isGasSpeciesLevelPopMC() const
    {
        for(uint i_species = 0; i_species < gas_species_level_pop_type.size(); i_species++)
            if(gas_species_level_pop_type[i_species] == POP_MC)
                return true;
        return false;
    }

    uint getMaxSubpixelLvl() const
    {
        return max_subpixel_lvl;
    }

    uint getTotalNrOfDustComponents() const
    {
        return (uint)dust_paths.size();
    }

    void resetNrOfDustComponents()
    {
        return dust_paths.clear();
    }

    void addDiffuseSource(dlist & val, string path)
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
        // Extent of the starfield [m]
        diffuse_sources.push_back(val[5]);
        // Stokes Q plarization
        diffuse_sources.push_back(val[6]);
        // Stokes U plarization
        diffuse_sources.push_back(val[7]);
        // Number of photons
        diffuse_sources.push_back(val[8]);
        diffuse_sources_str.push_back(path);
    }

    class plot_parameter
    {
      public:
        plot_parameter()
        {
            label = "";
            abs_min_cut = double(MAX_UINT);
            abs_max_cut = double(MAX_UINT);
            rel_min_cut = double(MAX_UINT);
            rel_max_cut = double(MAX_UINT);
            int_cut = 0;
            log = false;
            plot = true;
            normalized = true;
            scale = 1;
            offset = 0;

            pixel_bins = MAX_UINT;
            vec_bins = 26;
            vec_color[0] = 255;
            vec_color[1] = 255;
            vec_color[2] = 255;
        }

        ~plot_parameter()
        {}

        void addColorBarColor(double pos, double R, double G, double B)
        {
            cbar.push_back(pos);
            cbar.push_back(R);
            cbar.push_back(G);
            cbar.push_back(B);
        }

        void addContourLine(double val, double R, double G, double B)
        {
            cline.push_back(val);
            cline.push_back(R);
            cline.push_back(G);
            cline.push_back(B);
        }

        void setVectorColor(uchar R, uchar G, uchar B)
        {
            vec_color[0] = R;
            vec_color[1] = G;
            vec_color[2] = B;
        }

        string label;
        dlist cbar, cline;
        double abs_min_cut;
        double abs_max_cut;

        double rel_min_cut;
        double rel_max_cut;

        double int_cut;

        bool log;
        bool plot;
        bool normalized;

        double scale;
        double offset;

        uint vec_bins;
        uint pixel_bins;
        uchar vec_color[3];
    };

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
    double conv_mass_fraction;
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

    dlist midplane_3d_param;
    dlist star_mass;

    uilist plot_list;

    bool b_mrw;
    bool b_pda;
    bool b_enforced;
    bool is_speed_of_sound;
    bool peel_off;
    bool vel_maps;

    bool dust_offset, dust_gas_coupling;
    bool full_dust_temp, save_radiation_field;
    bool scattering_to_raytracing;
    bool split_dust_emision;
    bool individual_dust_fractions;

    strlist zeeman_catalog_path;

    uint phID;

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

    uint nrOfPlotPoints;
    uint nrOfPlotVectors;
    uint maxPlotLines;

    uint nr_of_mc_lvl_pop_photons;
    uint mc_lvl_pop_seed;

    uint healpix_orientation;

    maplist line_ray_detector_list;

    string rt_grid_description;

    dlist dust_mc_detectors;
    dlist dust_ray_detectors;
    dlist sync_ray_detectors;

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
    bool sublimate;

    Vector3D axis1, axis2;

    uint start;
    uint stop;

    double f_highJ;
    double f_cor;
    double Q_ref;
    double alpha_Q;

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
    uilist component_id_to_choice;

    strlist dust_paths;
    strlist size_keywords;

    bool reset_dust_files;

    // opiate parameters
    dlist opiate_ray_detectors;
    strlist opiate_spec_ids;
    string opiata_path_emi;
    string opiata_path_abs;
};

#endif
