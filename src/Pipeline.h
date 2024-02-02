#pragma once
#include "CommandParser.h"
#include "Detector.h"
#include "Dust.h"
#include "GasSpecies.h"
#include "Source.h"
#include "Typedefs.h"
#include "Parameters.h"
#include "Grid.h"
#include "Vector.h"

#ifndef CPIPELINE
#define CPIPELINE

class CPipeline
{
  public:
    CPipeline(void)
    {
        begin = 0;
        end = 0;
        len = 0;

        m = 0;
        s = 0;
        h = 0;

        ru[0] = '|';
        ru[1] = '/';
        ru[2] = '-';
        ru[3] = '\\';
    }

    ~CPipeline(void)
    {}

    bool calcMonteCarloRadiationField(parameters & param);
    bool calcPolarizationMapsViaRayTracing(parameters & param);
    bool calcPolarizationMapsViaSynchrotron(parameters & param);
    bool calcChMapsViaRayTracing(parameters & param);
    bool calcPolarizationMapsViaMC(parameters & param);
    bool calcOpiateMapsViaRayTracing(parameters & param);
    bool proberobeLineOfSight(parameters & param);

    bool assignDustMixture(parameters & param, CDustMixture * dust, CGridBasic * grid);
    bool assignGasSpecies(parameters & param, CGasMixture * gas, CGridBasic * grid);

    // bool calcRadPressure(parameter & param);

    bool solveNextVelocityField(uint itID)
    {
        return true;
    };
    // bool preparePressureData(CGridBasic * grid, CDustMixture * dust, parameters &
    // param, bool plot, uint itID);

    bool createWavelengthList(parameters & param, CDustMixture * dust, CGasMixture * gas, COpiateDataBase * op);

    void printParameters(parameters & param, uint max_id);

    bool Init(int argc, char ** argv);
    void Finish();
    void Error();
    void Run();

    bool createOutputPaths(string path)
    {
        const char * sep = SEP;
        string::size_type pos1 = 0, pos2 = 0;
        int len = -1;
        uint offset;

#ifdef WINDOWS
        offset = 0;
#else
        offset = 1;
#endif

        path_data = path + "data" + sep;
        path_plot = path + "plots" + sep;

        strlist folders;
        string folder;

        while(path.find(sep, offset) != string::npos)
        {
            pos1 = path.find(sep, offset);
            pos2 = path.find(sep, pos1 + 1);

            len = int(pos2 - pos1 + 1);

            if(len < 0)
                break;

            folder = path.substr(pos1 + 1, len - 2);
            path.erase(pos1 + 1, len - 1);
            folders.push_back(folder);
        }

        for(uint i = 0; i < folders.size(); i++)
        {
            folder = folders[i];
            path += folder + sep;

            if(!createPath(path))
            {
                cout << "\nERROR: Failed to create output folder for data!" << endl;
                cout << path << std::endl;
                return false;
            }
        }

        if(!createPath(path))
        {
            cout << "\nERROR: Failed to create output folder(s)!" << endl;
            cout << path << std::endl;
            return false;
        }

        if(!createPath(path_data))
        {
            cout << "\nERROR: Failed to create output folder for data!" << endl;
            cout << path_data << std::endl;
            return false;
        }

        if(!createPath(path_plot))
        {
            cout << "\nERROR: Failed to create output folder for plots!" << endl;
            cout << path_plot << std::endl;
            return false;
        }

        return true;
    }

    bool assignGridType(CGridBasic *& grid, parameters & param);

    bool createPath(string path)
    {
        const char * tmp_path = path.c_str();
#ifdef WINDOWS
        if(_access(tmp_path, 0) != 0)
            if(_mkdir(tmp_path))
                return false;
#else
        if(access(tmp_path, 0) != 0)
            if(mkdir(tmp_path, S_IRWXU | S_IRWXG | S_IRWXO))
                return false;
#endif

        return true;
    }

    void printConversionParameters(parameters & param)
    {
        cout << "Conversion factors" << endl;
        cout << "- Conv. gas dens. in SI        : " << param.getSIConvDH() << endl;
        cout << "- Conv. length in SI           : " << param.getSIConvLength() << endl;
        cout << "- Conv. vel. field in SI.      : " << param.getSIConvVField() << endl;
        cout << "- Conv. mag. field in SI.      : " << param.getSIConvBField() << endl;
        if(param.getIndividualDustMassFractions())
            cout << "- Mass fraction (Mdust/Mgas)   : set by the dust components" << endl;
        else
            cout << "- Mass fraction (Mdust/Mgas)   : " << param.getDustMassFraction() << endl;
        cout << "- Relative molecular mass (mu) : " << param.getMu() << endl;
    }

    void printAdditionalParameters(parameters & param)
    {
        cout << "Additional parameters" << endl;
        if(param.isRatSimulation())
        {
            cout << "- RAT efficiency ref. (Qref)   : " << param.getQref() << endl;
            cout << "- RAT efficiency exp. (alphaQ) : " << param.getAlphaQ() << endl;
        }
        else
            cout << "- None" << endl;
    }

    void printPathParameters(parameters & param)
    {
        cout << "- Path grid file   : " << param.getPathGrid() << endl;
        cout << "- Path output      : " << param.getPathOutput() << endl;
        cout << "- Number of threads: " << param.getNrOfThreads() << endl;
    }

    void printPlotParameters(parameters & param, bool input_output = false)
    {
        if(param.getNrOfPlotPoints() + param.getNrOfPlotVectors() + param.getInpMidDataPoints() +
               param.getOutMidDataPoints() + param.getInpAMIRAPoints() + param.getOutAMIRAPoints() >
           0)
            cout << "Plot parameters" << endl;

        if(param.getNrOfPlotPoints() + param.getNrOfPlotVectors() != 0)
        {
            cout << "- Raw data                     : " << param.getNrOfPlotPoints() << " points, ";
            cout << param.getNrOfPlotVectors() << " vectors, ";
            cout << param.getMaxPlotLines() << " lines" << endl;
        }

        if(param.getInpMidDataPoints() != 0 && param.getOutMidDataPoints() != 0)
        {
            cout << "- Midplane            (in,out) : " << param.getInpMidDataPoints() << ", "
                 << param.getOutMidDataPoints() << " pixel" << endl;
        }
        else
        {
            if(param.getInpMidDataPoints() != 0)
                cout << "- Midplane                (in) : " << param.getInpMidDataPoints() << " pixel"
                     << endl;

            if(param.getOutMidDataPoints() != 0)
                cout << "- Midplane               (out) : " << param.getOutMidDataPoints() << " pixel"
                     << endl;
        }

        if(param.getInpAMIRAPoints() != 0 && param.getOutAMIRAPoints() != 0)
            cout << "- Amira points (in,out)        : " << param.getInpAMIRAPoints() << ", "
                 << param.getOutAMIRAPoints() << " pixel" << endl;
        else
        {
            if(param.getInpAMIRAPoints() != 0)
                cout << "- Amira points            (in) : " << param.getInpAMIRAPoints() << " pixel" << endl;
            if(param.getOutAMIRAPoints() != 0)
                cout << "- Amira points           (out) : " << param.getOutAMIRAPoints() << " pixel" << endl;
        }
    }

    void printAlignmentParameters(parameters & param)
    {
        cout << "Dust grain alignment" << endl;

        if(param.getAligRANDOM())
            cout << "- No alignment" << endl;
        else if(param.getAligPA())
            cout << "- Perfect alignment" << endl;
        else if(param.getAligNONPA())
        {
            cout << "- Non-perfect alignment" << endl;
            cout << "    Rayleigh reduction fact. (R) : " << param.getRayleighReductionFactor() << endl;
        }
        else
        {
            if(param.getAligIDG())
                cout << "- IDG (paramagnetic alignment)" << endl;
            if(param.getAligRAT())
            {
                cout << "- RAT (radiative torques)" << endl;
                if(!param.getAligINTERNAL())
                    cout << "    Rayleigh reduction fact. (R) : " << param.getRayleighReductionFactor() << endl;
            }
            if(param.getAligGOLD())
                cout << "- GOLD (mechanical alignment)" << endl;
            if(param.getAligINTERNAL())
            {
                cout << "- Internal alignment" << endl;
                if(param.getAligIDG() || param.getAligGOLD())
                    cout << "    Correlation factor (f_c)     : " << param.getFcorr() << endl;
                if(param.getAligRAT())
                    cout << "    Fraction at high-J (f_highJ) : " << param.getFHighJ() << endl;
            }
            cout << "    Align limit prefact. (larmF) : " << param.getLarmF() << endl;
        }
    }

    void printSourceParameters(parameters & param, bool show_dust = false)
    {
        if(param.getNrOfSources() > 0 || param.isRaytracingSimulation())
        {
            cout << "Defined radiation source(s):" << endl;

            if(param.getNrOfPointSources() > 0 && (!param.isRaytracingSimulation() || param.getScatteringToRay()))
            {
                dlist sources_list = param.getPointSources();
                for(uint s = 0; s < sources_list.size(); s += NR_OF_POINT_SOURCES)
                {
                    // Index of current source
                    uint pos = s / NR_OF_POINT_SOURCES;
                    cout << "- Star " << pos + 1 << " :\n"
                         << "    Position    : "
                         << sources_list[s + 0] << ", "
                         << sources_list[s + 1] << ", "
                         << sources_list[s + 2] << " (x,y,z) [m]\n"
                         << "    Radius      : " << sources_list[s + 3] << " [R_sun]\n"
                         << "    Temperature : " << sources_list[s + 4] << " [K]\n"
                         << "    Stokes      : "
                         << sources_list[s + 5] << ", "
                         << sources_list[s + 6] << " (q,u)" << endl;
                }
            }
                // cout << "- Star(s)        : " << param.getNrOfPointSources() << endl;

            if(param.getNrOfDiffuseSources() > 0 && (!param.isRaytracingSimulation() || param.getScatteringToRay()))
            {
                dlist sources_list = param.getDiffuseSources();
                for(uint s = 0; s < sources_list.size(); s += NR_OF_DIFF_SOURCES)
                {
                    // Index of current source
                    uint pos = s / NR_OF_DIFF_SOURCES;
                    cout << "- Starfield " << pos + 1 << " :\n"
                         << "    Position    : "
                         << sources_list[s + 0] << ", "
                         << sources_list[s + 1] << ", "
                         << sources_list[s + 2] << " (x,y,z) [m]\n"
                         << "    Radius      : " << sources_list[s + 3] << " [R_sun]\n"
                         << "    Temperature : " << sources_list[s + 4] << " [K]\n"
                         << "    Variance    : " << sources_list[s + 5] << " [m]\n"
                         << "    Stokes      : "
                         << sources_list[s + 6] << ", "
                         << sources_list[s + 7] << " (q,u)" << endl;
                }
            }
                // cout << "- Star field(s)  : " << param.getNrOfDiffuseSources() << endl;

            if(param.getNrOfLaserSources() > 0 && (!param.isRaytracingSimulation() || param.getScatteringToRay()))
            {
                dlist sources_list = param.getLaserSources();
                for(uint s = 0; s < sources_list.size(); s += NR_OF_LASER_SOURCES)
                {
                    // Index of current source
                    uint pos = s / NR_OF_LASER_SOURCES;
                    cout << "- Laser " << pos + 1 << " :\n"
                         << "    Position           : "
                         << sources_list[s + 0] << ", "
                         << sources_list[s + 1] << ", "
                         << sources_list[s + 2] << " (x,y,z) [m]\n"
                         << "    Direction          : "
                         << sources_list[s + 3] << ", "
                         << sources_list[s + 4] << ", "
                         << sources_list[s + 5] << " (x,y,z) [m]\n"
                         << "    Total power        : " << sources_list[s + 6] << " [W]\n"
                         << "    Central wavelength : " << sources_list[s + 7] << " [m]\n"
                         << "    FWHM               : " << sources_list[s + 8] << " [m]\n"
                         << "    Stokes             : "
                         << sources_list[s + 9] << ", "
                         << sources_list[s + 10] << " (q,u)" << endl;
                }
            }
                // cout << "- Laser(s)  : " << param.getNrOfLaserSources() << endl;

            if(param.getNrOfBackgroundSources() > 0 || param.isRaytracingSimulation())
            {
                dlist sources_list = param.getBackgroundSources();
                for(uint s = 0; s < sources_list.size(); s += NR_OF_BG_SOURCES)
                {
                    // Index of current source
                    uint pos = s / NR_OF_BG_SOURCES;
                    cout << "- Background " << pos + 1 << " :\n"
                         << "    Scaling factor  : " << sources_list[s + 0] << "\n"
                         << "    Temperature     : " << sources_list[s + 1] << " [K]\n"
                         << "    Stokes          : "
                         << sources_list[s + 2] << ", "
                         << sources_list[s + 3] << ", "
                         << sources_list[s + 4] << " (q,u,v)\n"
                         << "    Rotation angles : "
                         << sources_list[s + 5] << ", "
                         << sources_list[s + 6] << " (1,2) [°]" << endl;
                }
            }
                // cout << "- Background(s)  : " << max(uint(1), param.getNrOfBackgroundSources()) << endl;

            if(param.getDustSource() && show_dust)
                cout << "- Dust as sources: yes" << endl;

            if(param.getISRFSource() && (!param.isRaytracingSimulation() || param.getScatteringToRay()))
                cout << "- ISRF as sources: yes" << endl;
        }
    }

    void printDetectorParameters(parameters & param, bool monte_carlo = false)
    {
        if(monte_carlo)
        {
            cout << "Monte-Carlo parameter" << endl;
            if(param.getNrOfDustMCDetectors() > 0)
                cout << "- Number of detectors   : " << param.getNrOfDustMCDetectors() << endl;
        }
        else
        {
            cout << "Raytrace parameter" << endl;
            if(param.getNrOfDustRayDetectors() > 0)
                cout << "- Number of detectors   : " << param.getNrOfDustRayDetectors() << endl;
            param.printRTGridDescription();
            cout << "- Start, Stop           : " << param.getStart() + 1 << ", " << param.getStop() + 1
                 << endl;
        }
        cout << "- Observer distance [m] : " << param.getMinObserverDistance() << " - "
             << param.getMaxObserverDistance() << endl;
        cout << "- Pixel in x-dir/nsides : " << param.getMinDetectorPixelX() << " (min) - "
             << param.getMaxDetectorPixelX() << " (max)" << endl;
        if(param.getMinDetectorPixelY() != 0 || param.getMaxDetectorPixelY() != 0)
            cout << "- Pixel in y-direction  : " << param.getMinDetectorPixelY() << " (min) - "
                 << param.getMaxDetectorPixelY() << " (max)" << endl;
        if(!monte_carlo)
            cout << "- Max subpixel level    : " << param.getMaxSubpixelLvl() << endl;

        if(param.getMaxSidelengthX() > 0)
        {
            cout << "- Sidelength X [m]      : " << param.getMinSidelengthX() << " (min) - "
                 << param.getMaxSidelengthX() << " (max)" << endl;
            if(param.getUseGridSidelengthX())
                cout << "  (some detectors are using the grid size)" << endl;
        }
        else if(param.getUseGridSidelengthX())
            cout << "- Sidelength X [m]      : only the grid size is used" << endl;
        if(param.getMaxSidelengthY() > 0)
        {
            cout << "- Sidelength Y [m]      : " << param.getMinSidelengthY() << " (min) - "
                 << param.getMaxSidelengthY() << " (max)" << endl;
            if(param.getUseGridSidelengthY())
                cout << "  (some detectors are using the grid size)" << endl;
        }
        else if(param.getUseGridSidelengthY())
            cout << "- Sidelength Y [m]      : only the grid size is used" << endl;

        if(!monte_carlo)
        {
            if(param.getMaxMapShiftX() != 0)
                cout << "- Map shift X [m]       : " << param.getMinMapShiftX() << " (min) - "
                     << param.getMaxMapShiftX() << " (max)" << endl;
            if(param.getMaxMapShiftY() != 0)
                cout << "- Map shift Y [m]       : " << param.getMinMapShiftY() << " (min) - "
                     << param.getMaxMapShiftY() << " (max)" << endl;
        }

        cout << "- Rotation axis         : n1 = ";
        param.getAxis1().printValues();
        cout << ", angle [°] = " << param.getMinDetectorAngle1() << " (min) - "
             << param.getMaxDetectorAngle1() << " (max)" << endl;
        cout << "- Rotation axis         : n2 = ";
        param.getAxis2().printValues();
        cout << ", angle [°] = " << param.getMinDetectorAngle2() << " (min) - "
             << param.getMaxDetectorAngle2() << " (max)" << endl;
    }

    void printSynchrotronParameters(parameters & param)
    {
        cout << SEP_LINE;
        cout << "Synchotron parameter" << endl;
        cout << SEP_LINE;
        cout << "Observed wavelengths" << endl;
        dlist sync_ray_detectors = param.getSyncRayDetectors();
        for(uint i = 0; i < sync_ray_detectors.size(); i += NR_OF_RAY_DET)
        {
            uint pos = i / NR_OF_RAY_DET;
            double lam_min = sync_ray_detectors[i + 0];
            double lam_max = sync_ray_detectors[i + 1];
            uint lam_skip = uint(sync_ray_detectors[i + 2]);
            cout << "    - Synchotron emission detector " << (pos + 1) << ": from wl = " << lam_min
                 << " [m] to wl = " << lam_max << " [m] with " << lam_skip << " step(s)" << endl;
        }
    }

    void deleteSourceLists()
    {
        for(uint s = 0; s < sources_mc.size(); s++)
            delete sources_mc[s];
        sources_mc.clear();

        for(uint s = 0; s < sources_ray.size(); s++)
            delete sources_ray[s];
        sources_ray.clear();
    }

    uint getNrOffsetEntriesRay(parameters & param, CDustMixture * dust, CGridBasic * grid)
    {
        uint nr_of_offset_entries = 0;
        if(param.getStochasticHeatingMaxSize() > 0)
        {
            // Add fields to store the stochastic heating propabilities for each cell
            for(uint i_mixture = 0; i_mixture < dust->getNrOfMixtures(); i_mixture++)
                nr_of_offset_entries +=
                    dust->getNrOfStochasticSizes(i_mixture) * dust->getNrOfCalorimetryTemperatures(i_mixture);
        }
        else if(dust->getScatteringToRay())
        {
            // Set that the spectral length is saved as vector
            grid->setSpecLengthAsVector(true);

            // Get list of dust raytracing detector parameters
            dlist dust_ray_detectors = param.getDustRayDetectors();

            // Get number of dust raytracing detectors
            uint nr_ray_detectors = uint(dust_ray_detectors.size()) / NR_OF_RAY_DET;

            // Check if one of the detectors is using the healpix background grid
            bool detector_uses_healpix = false;
            for(uint i = 0; i < nr_ray_detectors; i++)
            {
                uint pos = i * NR_OF_RAY_DET;
                uint detector_id = uint(dust_ray_detectors[pos + NR_OF_RAY_DET - 3]);

                if(detector_id == DET_SPHER)
                    detector_uses_healpix = true;
            }

            // Set up a list of all detector directions to use for the scattering via radiation field
            if(!detector_uses_healpix)
            {
                for(uint i_det = 0; i_det < nr_ray_detectors; i_det++)
                {
                    uint nr_spectral_bins = uint(dust_ray_detectors[i_det * NR_OF_RAY_DET + 2]);

                    // Add fields for the radiation field of each considered wavelength and detector
                    nr_of_offset_entries += 4 * nr_spectral_bins;
                }
            }
            else
            {
                nr_of_offset_entries += 4 * dust->getNrOfWavelength();
            }
        }

        return nr_of_offset_entries;
    }

    double getNrOfRayDetector(parameters & param)
    {
        dlist dust_ray_detectors = param.getDustRayDetectors();
        return uint(dust_ray_detectors.size()) / NR_OF_RAY_DET;
    }

  private:
    bool writeSources(parameters & param, CGridBasic * grid);
    void createSourceLists(parameters & param, CDustMixture * dust, CGridBasic * grid);
    CDetector * createDetectorList(parameters & param, CDustMixture * dust, CGridBasic * grid);
    string path_plot, path_data;
    slist sources_mc, sources_ray;
    long h, m, s;
    double begin, end;

    long len;
    parameter_list param_list;
    unsigned char ru[4];

    Vector3D ** det_coord_systems;
};

#endif

