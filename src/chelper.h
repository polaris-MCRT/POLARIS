#include "MathFunctions.h"
#include "Matrix2D.h"
#include "Stokes.h"
#include "Vector.h"
#include "typedefs.h"

#ifndef CHELPER
#define CHELPER

class cell_basic
{
  public:
    cell_basic()
    {
        data = 0;
    }

    virtual ~cell_basic()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    bool isValid()
    {
        return data != 0;
    }

    double getData(uint i)
    {
        if(data == 0)
            return 0;

        if(i == MAX_UINT)
            return 0;

        if(i > 4138637280)
            return 0;

        return data[i];
    }

    void setData(uint i, double d)
    {
#pragma omp atomic write
        data[i] = d;
    }

    void updateData(uint i, double d)
    {
#pragma omp atomic update
        data[i] += d;
    }

    void convertData(uint i, double c)
    {
#pragma omp atomic update
        data[i] *= c;
    }

    void resize(uint size)
    {
        if(data != 0)
            delete[] data;

        data = new double[size];
        for(uint i = 0; i < size; i++)
            data[i] = 0;
    }

    void setID(uint _id)
    {
        id = _id;
    }

    uint getID()
    {
        return id;
    }

    virtual ulong getUniqueID()
    {
        return ulong(id);
    }

    void updateID(uint _id)
    {
        id += _id;
    }

  protected:
    double * data;
    uint id;
};

class cell_oc : public cell_basic
{
  public:
    cell_oc()
    {
        parent = 0;
        children = 0;
        level = 0;
        id = 0;
        x_min = 0;
        y_min = 0;
        z_min = 0;
        length = 0;
        data = 0;
    }

    ~cell_oc()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    void setXmin(double _x_min)
    {
        x_min = _x_min;
    }

    void setYmin(double _y_min)
    {
        y_min = _y_min;
    }

    void setZmin(double _z_min)
    {
        z_min = _z_min;
    }

    void setLength(double _length)
    {
        length = _length;
    }

    void setLevel(uchar _level)
    {
        level = _level;
    }

    void setParent(cell_oc * _parent)
    {
        parent = _parent;
    }

    void setChildren(cell_oc * _children)
    {
        children = _children;
    }

    double getXmin()
    {
        return x_min;
    }

    double getYmin()
    {
        return y_min;
    }

    double getZmin()
    {
        return z_min;
    }

    double getXmax()
    {
        return x_min + length;
    }

    double getYmax()
    {
        return y_min + length;
    }

    double getZmax()
    {
        return z_min + length;
    }

    double getLength()
    {
        return length;
    }

    uchar getLevel()
    {
        return level;
    }

    ulong getUniqueID()
    {
        uint parent_id = parent->getID();
        ulong res = 8 ^ level + 8 * parent_id + id;
        return res;
    }

    cell_oc * getParent()
    {
        return parent;
    }

    cell_oc * getChildren()
    {
        return children;
    }

    cell_oc * getChild(uint i)
    {
        return &children[i];
    }

  private:
    double x_min, y_min, z_min, length;
    cell_oc * children;
    cell_oc * parent;
    uchar level;
};

class cell_vo : public cell_basic
{
  public:
    cell_vo()
    {
        id = 0;
        data = 0;
        nr_neighbors = 0;
        neighbors = 0;
        volume = -1;
    }

    ~cell_vo()
    {
        if(data != 0)
            delete[] data;

        data = 0;

        if(neighbors != 0)
            delete[] neighbors;

        neighbors = 0;
    }

    void initNeighbors(short nr)
    {
        nr_neighbors = nr;
        neighbors = new int[nr];

        for(ushort i = 0; i < nr; i++)
            neighbors[i] = 0;
    }

    void setNeighbor(uint pos, int id)
    {
        if(pos > 500)
            return;

        neighbors[pos] = id;
    }

    void setCenter(double cx, double cy, double cz)
    {
        center.set(cx, cy, cz);
    };

    void setVolume(double v)
    {
        volume = v;
    };

    Vector3D getCenter()
    {
        return center;
    };

    double getX()
    {
        return center.X();
    }

    double getY()
    {
        return center.Y();
    }

    double getZ()
    {
        return center.Z();
    }

    double getVolume()
    {
        return volume;
    };

    int getNeighborID(uint pos)
    {
        return neighbors[pos];
    }

    ushort getNrOfNeighbors()
    {
        return nr_neighbors;
    }

  private:
    Vector3D center;
    ushort nr_neighbors;
    int * neighbors;
    double volume;
};

class cell_sp : public cell_basic
{
  public:
    cell_sp()
    {
        rID = 0;
        phID = 0;
        thID = 0;
        data = 0;
        id = 0;
    }

    cell_sp(uint _rID, uint _phID, uint _thID)
    {
        rID = _rID;
        phID = _phID;
        thID = _thID;
        data = 0;
        id = 0;
    }

    ~cell_sp()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    void setRID(uint id)
    {
        rID = id;
    }

    void setPhID(uint id)
    {
        phID = id;
    }

    void setThID(uint id)
    {
        thID = id;
    }

    uint getRID()
    {
        return rID;
    }

    uint getPhID()
    {
        return phID;
    }

    uint getThID()
    {
        return thID;
    }

  private:
    uint rID, phID, thID;
};

class cell_cyl : public cell_basic
{
  public:
    cell_cyl()
    {
        rID = 0;
        phID = 0;
        zID = 0;
        data = 0;
        id = 0;
    }

    ~cell_cyl()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    void setRID(uint id)
    {
        rID = id;
    }

    void setPhID(uint id)
    {
        phID = id;
    }

    void setZID(uint id)
    {
        zID = id;
    }

    uint getRID()
    {
        return rID;
    }

    uint getPhID()
    {
        return phID;
    }

    uint getZID()
    {
        return zID;
    }

  private:
    uint rID, phID, zID;
};

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
        sublimate = false;
        individual_dust_fractions = false;

        nr_ofISRFPhotons = 0;
        nr_ofDustPhotons = 0;

        nrOfGnuPoints = 0;
        nrOfGnuVectors = 0;
        maxGridLines = 0;
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

        // opiate parmeter
        opiate_param_path = "";
        opiate_data_path = "";
    }

    ~parameters()
    {}

    string getOpiateParamPath()
    {
        return opiate_param_path;
    }

    string getOpiateDataPath()
    {
        return opiate_data_path;
    }

    Vector3D getAxis1()
    {
        return axis1;
    }

    Vector3D getAxis2()
    {
        return axis2;
    }

    uint getOutAMIRAPoints()
    {
        return nr_ofOutAMIRAPoints;
    }

    uint getInpAMIRAPoints()
    {
        return nr_ofInpAMIRAPoints;
    }

    bool plotInpMidPoints()
    {
        return plot_inp_points;
    }

    bool plotOutMidPoints()
    {
        return plot_out_points;
    }

    dlist getMidplane3dParams()
    {
        return midplane_3d_param;
    }

    uilist & getPlotList()
    {
        return plot_list;
    }

    bool isInPlotList(uint id)
    {
        if(plot_list.empty())
            return true;

        return (find(plot_list.begin(), plot_list.end(), id) != plot_list.end());
    }

    uint getInpMidDataPoints()
    {
        return nr_ofInpMidDataPoints;
    }

    uint getOutMidDataPoints()
    {
        return nr_ofOutMidDataPoints;
    }

    dlist getOpiateSequence()
    {
        return line_opiate_detectors;
    }

    uint getMidplaneZoom()
    {
        return midplane_zoom;
    }

    int getCommand()
    {
        return cmd;
    }

    bool isRatSimulation()
    {
        if(getCommand() == CMD_RAT || getCommand() == CMD_TEMP_RAT)
            return true;
        return false;
    }

    bool isMonteCarloSimulation()
    {
        if(getCommand() == CMD_TEMP || getCommand() == CMD_TEMP_RAT || getCommand() == CMD_RAT)
            return true;
        return false;
    }

    bool isRaytracing()
    {
        if(getCommand() == CMD_OPIATE || getCommand() == CMD_DUST_EMISSION ||
           getCommand() == CMD_SYNCHROTRON || getCommand() == CMD_LINE_EMISSION)
            return true;
        return false;
    }

    bool isTemperatureSimulation()
    {
        if(getCommand() == CMD_TEMP || getCommand() == CMD_TEMP_RAT)
            return true;
        return false;
    }

    double getStarMass(uint i)
    {
        return star_mass[i];
    }

    string getPathGrid()
    {
        return path_grid;
    }

    string getPathOutput()
    {
        return path_output;
    }

    string getPathInput()
    {
        return path_input;
    }

    uint getMinDetectorPixelX()
    {
        return min_detector_pixel_x;
    }

    uint getMaxDetectorPixelX()
    {
        return max_detector_pixel_x;
    }

    uint getMinDetectorPixelY()
    {
        return min_detector_pixel_y;
    }

    uint getMaxDetectorPixelY()
    {
        return max_detector_pixel_y;
    }

    double getMinDetectorAngle1()
    {
        return min_rot_angle_1;
    }

    double getMaxDetectorAngle1()
    {
        return max_rot_angle_1;
    }

    double getMinDetectorAngle2()
    {
        return min_rot_angle_2;
    }

    double getMaxDetectorAngle2()
    {
        return max_rot_angle_2;
    }

    double getMinSidelengthX()
    {
        return min_sidelength_x;
    }

    double getMaxSidelengthX()
    {
        return max_sidelength_x;
    }

    double getMinSidelengthY()
    {
        return min_sidelength_y;
    }

    double getMaxSidelengthY()
    {
        return max_sidelength_y;
    }

    bool getUseGridSidelengthX()
    {
        return use_grid_sidelength_x;
    }

    bool getUseGridSidelengthY()
    {
        return use_grid_sidelength_y;
    }

    double getMinMapShiftX()
    {
        return min_ray_map_shift_x;
    }

    double getMinMapShiftY()
    {
        return min_ray_map_shift_y;
    }

    double getMaxMapShiftX()
    {
        return max_ray_map_shift_x;
    }

    double getMaxMapShiftY()
    {
        return max_ray_map_shift_y;
    }

    double getSIConvLength()
    {
        return conv_l_in_SI;
    }

    double getSIConvDH()
    {
        return conv_dH_in_SI;
    }

    double getDelta0()
    {
        return delta0;
    }

    double getLarmF()
    {
        return larm_f;
    }

    double getSIConvBField()
    {
        return conv_B_in_SI;
    }

    double getSIConvVField()
    {
        return conv_V_in_SI;
    }

    bool getDustOffset()
    {
        return dust_offset;
    }

    bool getDustGasCoupling()
    {
        return dust_gas_coupling;
    }

    double getOffsetMinGasDensity()
    {
        return offset_min_gas_dens;
    }

    bool getDustTempMulti()
    {
        return full_dust_temp;
    }

    double getSizeMin(uint i)
    {
        return a_min_global[i];
    }

    double getSizeMax(uint i)
    {
        return a_max_global[i];
    }

    double getMaterialDensity(uint i)
    {
        return material_density[i];
    }

    bool getDustSource()
    {
        return nr_ofDustPhotons > 0;
    }

    bool getISRFSource()
    {
        return nr_ofISRFPhotons > 0;
    }

    ullong getNrOfDustPhotons()
    {
        return nr_ofDustPhotons;
    }

    double getDustMassFraction()
    {
        return conv_mass_fraction;
    }

    uint getAlign()
    {
        return align;
    }

    bool getAligRANDOM()
    {
        return align == 0;
    }

    bool getAligPA()
    {
        return (align & ALIG_PA) == ALIG_PA;
    }

    bool getAligIDG()
    {
        return (align & ALIG_IDG) == ALIG_IDG;
    }

    bool getAligRAT()
    {
        return (align & ALIG_RAT) == ALIG_RAT;
    }

    bool getAligGOLD()
    {
        return (align & ALIG_GOLD) == ALIG_GOLD;
    }

    bool getAligINTERNAL()
    {
        return (align & ALIG_INTERNAL) == ALIG_INTERNAL;
    }

    double getMu()
    {
        return mu;
    }

    bool getMRW()
    {
        return b_mrw;
    }

    bool getPDA()
    {
        return b_pda;
    }

    bool getEnfScattering()
    {
        return b_enforced;
    }

    double getStochasticHeatingMaxSize()
    {
        return stochastic_heating_max_size;
    }

    bool getSaveRadiationField()
    {
        return save_radiation_field;
    }

    bool getScatteringToRay()
    {
        return scattering_to_raytracing;
    }

    bool getIndividualDustMassFractions()
    {
        return individual_dust_fractions;
    }

    bool getIsSpeedOfSound()
    {
        return is_speed_of_sound;
    }

    bool getPeelOff()
    {
        return peel_off;
    }

    double getForegroundExtinctionMagnitude()
    {
        return extinction_magnitude;
    }

    double getForegroundExtinctionWavelength()
    {
        return extinction_magnitude_wavelength;
    }

    uint getForegroundExtinctionDustMixture()
    {
        return extinction_i_mixture;
    }

    bool getVelFieldType()
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

    uint getWriteRadiationField()
    {
        return write_radiation_field;
    }

    bool getWriteGZero()
    {
        return write_g_zero;
    }

    double getISRFGZero()
    {
        return isrf_g_zero;
    }

    double getISRFRadius()
    {
        return isrf_radius;
    }

    string getISRFPath()
    {
        return isrf_path;
    }

    string getZeemanCatalog(uint i_species)
    {
        return zeeman_catalog_path[i_species];
    }

    uint getAlignmentMechanism()
    {
        return align;
    }

    double getMinObserverDistance()
    {
        return min_obs_distance;
    }

    double getMaxObserverDistance()
    {
        return max_obs_distance;
    }

    double getKeplerStarMass()
    {
        return kepler_star_mass;
    }

    double getTurbulentVelocity()
    {
        return turbulent_velocity;
    }

    uint getTaskID()
    {
        return task_id;
    }

    uint getNrOfThreads()
    {
        return nr_ofThreads;
    }

    ullong getNrOfISRFPhotons()
    {
        return nr_ofISRFPhotons;
    }

    uint getNrOfMixtures()
    {
        if(!dust_choices.empty())
            return dust_choices.size();
        else
            return 1;
    }

    uilist getDustComponentChoices()
    {
        return dust_choices;
    }

    uint getPhaseFunctionID()
    {
        return phID;
    }

    double getFHighJ()
    {
        return f_highJ;
    }

    double getFcorr()
    {
        return f_cor;
    }

    double getAdjTgas()
    {
        return adjTgas;
    }

    uint getNrOfDiffuseSources()
    {
        return uint(diffuse_sources.size() / NR_OF_DIFF_SOURCES);
    }

    uint getNrOfPointSources()
    {
        return uint(point_sources.size() / NR_OF_POINT_SOURCES);
    }

    uint getNrOfLaserSources()
    {
        return uint(laser_sources.size() / NR_OF_LASER_SOURCES);
    }

    uint getNrOfBackgroundSources()
    {
        return uint(background_sources.size() / NR_OF_BG_SOURCES);
    }

    double getXYMin()
    {
        return xymin;
    }

    double getXYMax()
    {
        return xymax;
    }

    double getXYSteps()
    {
        return xysteps;
    }

    uint getXYBins()
    {
        return xy_bins;
    }

    string getXYLabel()
    {
        return xylabel;
    }

    bool isAutoScale()
    {
        return autoscale;
    }

    uint getNrOfSources()
    {
        uint res = getNrOfPointSources() + getNrOfDiffuseSources() + getNrOfBackgroundSources() +
                   getNrOfLaserSources();

        if(nr_ofDustPhotons > 0)
            res++;

        if(nr_ofISRFPhotons > 0)
            res++;

        // Add gas source for levl population calculation
        if(gasSpeciesLevelPopTypeIsMC())
            res++;

        return res;
    }

    uint getNrOfGnuPoints()
    {
        return nrOfGnuPoints;
    }

    uint getNrOfGnuVectors()
    {
        return nrOfGnuVectors;
    }

    uint getmaxGridLines()
    {
        return maxGridLines;
    }

    uint getStart()
    {
        return start;
    }

    uint getStop()
    {
        return stop;
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

    void setNrOfGnuPoints(uint val)
    {
        nrOfGnuPoints = val;
    }

    void setNrOfGnuVectors(uint val)
    {
        nrOfGnuVectors = val;
    }

    void setmaxGridLines(uint val)
    {
        maxGridLines = val;
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

    void setDustOffset(double offset_min_gas_dens)
    {
        dust_offset = true;
        offset_min_gas_dens = offset_min_gas_dens;
    }

    void setDustGasCoupling(bool val)
    {
        dust_gas_coupling = val;
    }

    void setDustGasCoupling(double offset_min_gas_dens)
    {
        dust_gas_coupling = true;
        offset_min_gas_dens = offset_min_gas_dens;
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
        // dust detector type/grid (all)
        dust_ray_detectors.push_back(val[11]);
        // number of pixel in x-dir. (cart, polar, slice) / N_side (healpix)
        dust_ray_detectors.push_back(val[12]);
        // number of pixel in y-direction (cart, polar, slice)
        dust_ray_detectors.push_back(val[13]);

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
        // dust detector type/grid (all)
        sync_ray_detectors.push_back(val[11]);
        // number of pixel in x-dir. (cart, polar, slice) / N_side (healpix)
        sync_ray_detectors.push_back(val[12]);
        // number of pixel in y-direction (cart, polar, slice)
        sync_ray_detectors.push_back(val[13]);

        switch(uint(val[NR_OF_RAY_DET - 3]))
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

    void addLineOpiateDetector(dlist & val)
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
    }

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
        // number of pixel in x-direction
        dust_mc_detectors.push_back(val[8]);
        // number of pixel in y-direction
        dust_mc_detectors.push_back(val[9]);
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

    uint getHealpixOrientation()
    {
        return healpix_orientation;
    }

    dlist getLineRayDetector(uint i_species)
    {
        return line_ray_detector_list[i_species];
    }

    maplist getLineRayDetectors()
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

    bool getVelMaps()
    {
        return vel_maps;
    }

    uint getDustChoiceFromComponentId(uint i)
    {
        return component_id_to_choice[i];
    }

    uint getDustChoiceFromMixtureId(uint i)
    {
        return dust_choices[i];
    }

    void printRTGridDescription()
    {
        cout << "- background grid shape : " << rt_grid_description << endl;
    }

    double getAcceptanceAngle()
    {
        return acceptance_angle;
    }

    string getDustPath(uint i)
    {
        return dust_paths[i];
    }

    bool isSublimate()
    {
        return sublimate;
    }

    double getDustFraction(uint i)
    {
        return dust_fractions[i];
    }

    string getDustSizeKeyword(uint i)
    {
        return size_keywords[i];
    }

    dlist getDustSizeParameter(uint i_comp)
    {
        dlist res;
        for(uint i = 0; i < NR_OF_SIZE_DIST_PARAM; i++)
            res.push_back(size_parameter_map[i][i_comp]);
        return res;
    }

    string getGasSpeciesCatalogPath(uint i_species)
    {
        return gas_species_cat_path[i_species];
    }

    double getGasSpeciesAbundance(uint i_species)
    {
        return gas_species_abundance[i_species];
    }

    uint getGasSpeciesLevelPopType(uint i_species)
    {
        return gas_species_level_pop_type[i_species];
    }

    bool gasSpeciesLevelPopTypeIsMC()
    {
        for(uint i_species = 0; i_species < gas_species_level_pop_type.size(); i_species++)
            if(gas_species_level_pop_type[i_species] == POP_MC)
                return true;
        return false;
    }

    uint getMaxSubpixelLvl()
    {
        return max_subpixel_lvl;
    }

    uint getTotalNrOfDustComponents()
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

    uint nrOfGnuPoints;
    uint nrOfGnuVectors;
    uint maxGridLines;

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

    // opiate data
    string opiate_param_path;
    string opiate_data_path;
    dlist line_opiate_detectors;
};

class photon_package
{
  public:
    photon_package()
    {
        cell_pos = 0;
        tmp_path = 0;

        wavelength = 0;
        frequency = 0;
        velocity = 0;

        wID = MAX_UINT;
        dirID = MAX_UINT;

        stokes.set(1, 0, 0, 0, 0);
        multi_stokes = 0;

        hasSpare = false;
        rand1 = 0.5;
        rand2 = 0.5;

        X1 = 0;
        X2 = 0;
        call = 0;
    }

    double randn(double mu, double sigma)
    {
        double U1, U2, W, mult;
        double X1, X2;
        int call = 0;

        if(call == 1)
        {
            call = !call;
            return (mu + sigma * (double)X2);
        }

        do
        {
            U1 = -1 + getRND() * 2;
            U2 = -1 + getRND() * 2;
            W = pow(U1, 2) + pow(U2, 2);
        } while(W >= 1 || W == 0);

        mult = sqrt((-2 * log(W)) / W);
        X1 = U1 * mult;
        X2 = U2 * mult;

        call = !call;

        double res = mu + sigma * (double)X1;

        if(res < 0)
            return randn(mu, sigma);

        return res;
    }

    double generateGaussianNoise(const double & variance)
    {
        hasSpare = false;
        rand1 = getRND(), rand2 = getRND();

        if(hasSpare)
        {
            hasSpare = false;
            return sqrt(variance * rand1) * sin(rand2);
        }

        hasSpare = true;

        rand1 = getRND();
        if(rand1 < 1e-100)
            rand1 = 1e-100;
        rand1 = -2 * log(rand1);
        rand2 = (getRND() / ((double)RAND_MAX)) * PIx2;

        return sqrt(variance * rand1) * cos(rand2);
    }

    ~photon_package()
    {
        if(multi_stokes != 0)
            delete[] multi_stokes;
    }

    void initRandomGenerator(ullong seed)
    {
        rand_gen.setSeed(seed);
    }

    void setD(Matrix2D & _mD)
    {
        mD = _mD;
    }

    Matrix2D & getD()
    {
        return mD;
    }

    Vector3D & getPosition()
    {
        return pos;
    }

    Vector3D & getDirection()
    {
        return ez;
    }

    StokesVector & getStokesVector()
    {
        return stokes;
    }

    double const getRND()
    {
        return rand_gen.getValue();
    }

    double & getTmpPathLength()
    {
        return tmp_path;
    }

    cell_basic * getPositionCell()
    {
        return cell_pos;
    }

    uint getWavelengthID()
    {
        return wID;
    }

    double getWavelength()
    {
        if(wavelength > 0)
            return wavelength;
        else if(frequency > 0)
            return con_c / frequency;
        else
        {
            cout << "ERROR: Frequency not set correctly in photon package!" << endl;
            return 0;
        }
    }

    double getFrequency()
    {
        if(frequency > 0)
            return frequency;
        else if(wavelength > 0)
            return con_c / wavelength;
        else
        {
            cout << "ERROR: Wavelength not set correctly in photon package!" << endl;
            return 0;
        }
    }

    double getVelocity(double f_0 = 0)
    {
        if(velocity > 0)
            return velocity;
        else if(getFrequency() > 0 && f_0 > 0)
            return getFrequency() * con_c / f_0;
        else
        {
            cout << "ERROR: Velocity not set correctly in photon package!" << endl;
            return 0;
        }
    }

    StokesVector & getMultiStokesVector(uint vch)
    {
        return multi_stokes[vch];
    }

    void calcRandomDirection()
    {
        ez.rndDir(getRND(), getRND());
    }

    void initCoordSystem()
    {
        double phi = atan3(ez.Y(), -ez.X());
        double theta = acos(ez.Z());
        double cos_phi = cos(phi);
        double sin_phi = sin(phi);
        double cos_theta = ez.Z();
        double sin_theta = sin(theta);

        mD = CMathFunctions::getRotationMatrix(cos_phi, sin_phi, cos_theta, sin_theta);

        ex = mD * Vector3D(1, 0, 0);
        ey = mD * Vector3D(0, 1, 0);
    }

    void updateCoordSystem()
    {
        double phi = atan3(ez.Y(), -ez.X());
        double theta = acos(ez.Z());
        double cos_phi = cos(phi);
        double sin_phi = sin(phi);
        double cos_theta = ez.Z();
        double sin_theta = sin(theta);

        Matrix2D D_help = CMathFunctions::getRotationMatrix(cos_phi, sin_phi, cos_theta, sin_theta);
        mD = mD * D_help;

        ex = mD * Vector3D(1, 0, 0);
        ey = mD * Vector3D(0, 1, 0);
        ez = mD * Vector3D(0, 0, 1);
    }

    void updateCoordSystem(double phi, double theta)
    {
        double cos_phi = cos(phi);
        double sin_phi = sin(phi);
        double cos_theta = cos(theta);
        double sin_theta = sin(theta);

        Matrix2D D_help = CMathFunctions::getRotationMatrix(cos_phi, sin_phi, cos_theta, sin_theta);
        mD = mD * D_help;

        ex = mD * Vector3D(1, 0, 0);
        ey = mD * Vector3D(0, 1, 0);
        ez = mD * Vector3D(0, 0, 1);
    }

    void adjustPosition(Vector3D _pos, double _len)
    {
        tmp_path = _len;
        pos = _pos + tmp_path * ez;

        // Calculate cell by position
        if(tmp_path > 0)
            dirID = MAX_UINT;
    }

    void setPosition(Vector3D val)
    {
        pos = val;
    }

    void setBackupPosition(Vector3D val)
    {
        backup_pos = val;
    }

    void updateBackupPosition()
    {
        backup_pos = pos;
    }

    void resetPositionToBackupPos()
    {
        pos = backup_pos;
    }

    bool reachedBackupPosition()
    {
        if(ez * (backup_pos - pos) <= 0)
            return true;
        return false;
    }

    bool reachedBackupPosition(Vector3D tmp_pos)
    {
        if(ez * (backup_pos - tmp_pos) <= 0)
            return true;
        return false;
    }

    void setDetectorProjection()
    {
        double tmp_vec = ey * pos;
        pos.setX(ex * pos);
        pos.setY(tmp_vec);
    }

    void setRelativePosition(double tx, double ty, double tz)
    {
        pos = tx * ex + ty * ey + tz * ez;
    }

    void setDirection(Vector3D val)
    {
        ez = val;
    }

    void setStokesVector(StokesVector val)
    {
        stokes = val;
    }

    void initMultiStokesVector(uint nr_stokes_vector)
    {
        multi_stokes = new StokesVector[nr_stokes_vector];
    }

    void setMultiStokesVector(StokesVector st, uint i)
    {
        multi_stokes[i] = st;
    }

    void addStokesVector(StokesVector val)
    {
        stokes += val;
    }

    Vector3D getEX()
    {
        return ex;
    }

    Vector3D getEY()
    {
        return ey;
    }

    Vector3D getEZ()
    {
        return ez;
    }

    void setEX(Vector3D _e)
    {
        ex = _e;
    }

    void setEY(Vector3D _e)
    {
        ey = _e;
    }

    void setEZ(Vector3D _e)
    {
        ez = _e;
    }

    void setCoordinateSystem(Vector3D _ex, Vector3D _ey, Vector3D _ez)
    {
        ex = _ex;
        ey = _ey;
        ez = _ez;
    }

    void getCoordinateSystem(Vector3D & _ex, Vector3D & _ey, Vector3D & _ez)
    {
        _ex = ex;
        _ey = ey;
        _ez = ez;
    }

    void setTmpPathLength(double val)
    {
        tmp_path = val;
    }

    void setPositionCell(cell_basic * val)
    {
        cell_pos = val;
    }

    void setDirectionID(uint val)
    {
        dirID = val;
    }

    uint getDirectionID()
    {
        return dirID;
    }

    void setWavelength(uint _wID, double val)
    {
        wID = _wID;
        wavelength = val;
        frequency = 0;
        velocity = 0;
    }

    void setFrequency(uint _wID, double val)
    {
        wID = _wID;
        frequency = val;
        wavelength = 0;
        velocity = 0;
    }

    void setVelocity(uint _wID, double val)
    {
        wID = _wID;
        velocity = val;
        frequency = 0;
        wavelength = 0;
    }

  private:
    CRandomGenerator rand_gen;
    Vector3D pos;
    Vector3D backup_pos;
    Vector3D ex;
    Vector3D ey;
    Vector3D ez;
    Matrix2D mD;
    StokesVector stokes;
    StokesVector * multi_stokes;

    double wavelength;
    double frequency;
    double velocity;

    uint wID;
    uint dirID;

    double tmp_path;

    cell_basic * cell_pos;

    bool hasSpare;
    double rand1, rand2;

    double X1, X2;
    int call;
};
#endif
