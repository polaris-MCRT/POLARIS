#pragma once
#include "Detector.h"
#include "Dust.h"
#include "GasSpecies.h"
#include "Grid.h"
#include "MathFunctions.h"
#include "OcTree.h"
#include "Photon.h"
#include "Source.h"
#include "Typedefs.h"
#include "Raytracing.h"
#include "Synchrotron.h"
#include "Matrix2D.h"
#include "Parameters.h"
#include "Stokes.h"
#include "Vector.h"

class CRadiativeTransfer
{
  public:
    CRadiativeTransfer(parameters & param)
    {
        grid = 0;
        dust = 0;
        gas = 0;
        op = 0;
        tracer = 0;
        b_forced = true;
        peel_off = false;
        mrw_step = false;
        stokes_dust_rad_field = false;

        start = MAX_UINT;
        stop = MAX_UINT;

        adjTgas = 0;

        probing_points = 0;
        detector = 0;
        nr_mc_detectors = 0;
        nr_ray_detectors = 0;

        RK_c = 0;
        RK_b1 = 0;
        RK_b2 = 0;

        probing_points = 0;

        synchrotron = 0;

        pathOutput = param.getPathOutput();
    }

    ~CRadiativeTransfer(void)
    {
        if(probing_points != 0)
            delete[] probing_points;

        if(tracer != 0)
        {
            for(uint i = 0; i < nr_ray_detectors; i++)
                delete tracer[i];
            delete[] tracer;
        }

        if(RK_c != 0)
            delete[] RK_c;

        if(RK_b1 != 0)
            delete[] RK_b1;

        if(RK_b2 != 0)
            delete[] RK_b2;

        if(synchrotron != 0)
            delete synchrotron;
    }

    bool initiateDustRaytrace(parameters & param);
    bool initiateSyncRaytrace(parameters & param);
    bool initiateLineRaytrace(parameters & param);
    bool initiateOPIATERaytrace(parameters & param);
    bool initiateProbing(parameters & param);
    void initiateRadFieldMC(parameters & param);
    void initiateDustMC(parameters & param);

    void initiateRungeKuttaFehlberg()
    {
        RK_c = new double[6];
        RK_c[0] = 0.0;
        RK_c[1] = 1.0 / 4.0;
        RK_c[2] = 3.0 / 8.0;
        RK_c[3] = 12.0 / 13.0;
        RK_c[4] = 1.0;
        RK_c[5] = 0.5;

        RK_b1 = new double[6];
        RK_b1[0] = 16.0 / 135.0;
        RK_b1[1] = 0.0;
        RK_b1[2] = 6656.0 / 12825.0;
        RK_b1[3] = 28561.0 / 56430.0;
        RK_b1[4] = -9.0 / 50.0;
        RK_b1[5] = 2.0 / 55.0;

        RK_b2 = new double[6];
        RK_b2[0] = 25.0 / 216.0;
        RK_b2[1] = 0.0;
        RK_b2[2] = 1408.0 / 2565.0;
        RK_b2[3] = 2197.0 / 4104.0;
        RK_b2[4] = -1.0 / 5.0;
        RK_b2[5] = 0.0;

        RK_a.resize(6, 6);
        RK_a.addValue(0, 0, 0.0);
        RK_a.addValue(1, 0, 0.0);
        RK_a.addValue(2, 0, 0.0);
        RK_a.addValue(3, 0, 0.0);
        RK_a.addValue(4, 0, 0.0);
        RK_a.addValue(5, 0, 0.0);

        RK_a.addValue(0, 1, 1.0 / 4.0);
        RK_a.addValue(1, 1, 0.0);
        RK_a.addValue(2, 1, 0.0);
        RK_a.addValue(3, 1, 0.0);
        RK_a.addValue(4, 1, 0.0);
        RK_a.addValue(5, 1, 0.0);

        RK_a.addValue(0, 2, 3.0 / 32.0);
        RK_a.addValue(1, 2, 9.0 / 32.0);
        RK_a.addValue(2, 2, 0.0);
        RK_a.addValue(3, 2, 0.0);
        RK_a.addValue(4, 2, 0.0);
        RK_a.addValue(5, 2, 0.0);

        RK_a.addValue(0, 3, 1932.0 / 2197.0);
        RK_a.addValue(1, 3, -7200.0 / 2197.0);
        RK_a.addValue(2, 3, 7296.0 / 2197.0);
        RK_a.addValue(3, 3, 0.0);
        RK_a.addValue(4, 3, 0.0);
        RK_a.addValue(5, 3, 0.0);

        RK_a.addValue(0, 4, 439.0 / 216.0);
        RK_a.addValue(1, 4, -8.0);
        RK_a.addValue(2, 4, 3680.0 / 513.0);
        RK_a.addValue(3, 4, -845.0 / 4104.0);
        RK_a.addValue(4, 4, 0.0);
        RK_a.addValue(5, 4, 0.0);

        RK_a.addValue(0, 5, -8.0 / 27.0);
        RK_a.addValue(1, 5, 2.0);
        RK_a.addValue(2, 5, -3544.0 / 2565.0);
        RK_a.addValue(3, 5, 1859.0 / 4104.0);
        RK_a.addValue(4, 5, -11.0 / 40.0);
        RK_a.addValue(5, 5, 0.0);
    }

    // Temperature calculation and RATs
    bool calcMonteCarloRadiationField(uint command, bool use_energy_density, bool disable_reemission = false);
    bool calcMonteCarloLvlPopulation(uint i_species, uint global_seed);
    void rayThroughCellForLvlPop(photon_package * pp,
                                 uint i_species,
                                 uint i_trans,
                                 const VelFieldInterp & vel_field_interp);

    // Set temperature (old!)
    // bool setTemperatureDistribution();

    // Dust scattered light
    bool calcPolMapsViaMC();

    // Dust emission
    bool calcPolMapsViaRaytracing(parameters & param);
    void getDustPixelIntensity(CSourceBasic * tmp_source,
                               double cx,
                               double cy,
                               uint i_det,
                               uint subpixel_lvl,
                               int pos_id);
    void getDustIntensity(photon_package * pp,
                          CSourceBasic * tmp_source,
                          double cx,
                          double cy,
                          uint i_det,
                          uint subpixel_lvl);
    void rayThroughCellDust(photon_package * pp, uint i_det, uint nr_used_wavelengths);
    void calcStellarEmission(uint i_det, CRandomGenerator * rand_gen);

    // Synchrontron emission
    bool calcSyncMapsViaRaytracing(parameters & param);
    void getSyncPixelIntensity(CSourceBasic * tmp_source,
                               double cx,
                               double cy,
                               uint i_det,
                               uint subpixel_lvl,
                               int pos_id);
    void getSyncIntensity(photon_package * pp1,
                          CSourceBasic * tmp_source,
                          double cx,
                          double cy,
                          uint i_det,
                          uint subpixel_lvl);
    void rayThroughCellSync(photon_package * pp1, uint i_det, uint nr_used_wavelengths);


    //OPIATE database RT

    bool calcOPIATEMapsViaRaytracing(parameters & param);
    void getOPIATEPixelIntensity(CSourceBasic * tmp_source,
                               double cx,
                               double cy,
                               uint i_species,
                               uint i_trans,
                               uint i_det,
                               uint subpixel_lvl,
                               int pos_id);

    void getOPIATEIntensity(photon_package * pp,
                          CSourceBasic * tmp_source,
                          double cx,
                          double cy,
                          uint i_species,
                          uint i_trans,
                          uint i_det,
                          uint subpixel_lvl);

    void rayThroughCellOPIATE(photon_package * pp,
                            uint i_species,
                            uint i_trans,
                            uint i_det,
                            uint nr_velocity_channels,
                            const VelFieldInterp & vel_field_interp);


    // Line emission
    bool calcChMapsViaRaytracing(parameters & param);

    void getLinePixelIntensity(CSourceBasic * tmp_source,
                               double cx,
                               double cy,
                               uint i_species,
                               uint i_trans,
                               uint i_det,
                               uint subpixel_lvl,
                               int pos_id);

    void getLineIntensity(photon_package * pp,
                          CSourceBasic * tmp_source,
                          double cx,
                          double cy,
                          uint i_species,
                          uint i_trans,
                          uint i_det,
                          uint subpixel_lvl);
    void rayThroughCellLine(photon_package * pp,
                            uint i_species,
                            uint i_trans,
                            uint i_det,
                            uint nr_velocity_channels,
                            const VelFieldInterp & vel_field_interp);

    void preCalcVelocityInterp(CGridBasic * grid,
                               const photon_package & pp,
                               VelFieldInterp * vel_field_interp);

    // Calc radiation pressure
    // bool calcRadiativePressure(parameter & param);

    void updateRadiationField(photon_package * pp)
    {
        double energy = pp->getTmpPathLength() * pp->getStokesVector()->I();

        if(stokes_dust_rad_field)
        {
            // Rotate vector of radiation field to cell center
            Vector3D rad_field_dir = grid->rotateToCenter(*pp);

            // Backup original direction of photon
            Vector3D old_dir = pp->getDirection();

            // For each detector check if wavelength fits
            for(uint i_det = 0; i_det < nr_ray_detectors; i_det++)
            {
                // Set coordinate system of temporary photon package for the map direction
                tracer[i_det]->setCoordinateSystem(pp);

                // Go through each wavelength
                for(uint i_wave = 0; i_wave < tracer[i_det]->getNrSpectralBins(); i_wave++)
                {
                    // If the wavelengths fit, save Stokes
                    if(pp->getWavelength() == tracer[i_det]->getWavelength(i_wave))
                    {
                        // Save the scattering Stokes vector in the grid
                        grid->updateSpecLength(
                            pp->getPositionCell(),
                            detector_wl_index[i_det] + i_wave,
                            dust->getRadFieldScatteredFraction(grid, *pp, rad_field_dir, energy));
                    }
                }
            }

            // Recopy original direction of photon
            pp->setDirection(old_dir);
        }
        else
        {
            grid->updateSpecLength(pp, energy);
        }
    }

    void setGrid(CGridBasic * _grid)
    {
        grid = _grid;
    }

    void setDust(CDustMixture * _dust)
    {
        dust = _dust;
    }

    void setGas(CGasMixture * _gas)
    {
        gas = _gas;
    }

    void setOpiateDataBase(COpiateDataBase * _op)
    {
        op=_op;
    }

    void setSourcesLists(slist & _sources_mc, slist & _sources_ray)
    {
        sources_mc = _sources_mc;
        sources_ray = _sources_ray;
    }

    void setDetectors(CDetector * d)
    {
        detector = d;
    }

    double getOpticalDepthAlongPath(photon_package * pp);

    bool photonInDetectorDir(photon_package * pp, CDetector * detector);

    void scaleAddToDetector(photon_package * pp, CDetector * detector, ullong interactions, uint sourceID);

    bool doMRWStepBW(photon_package * pp);
    bool doMRWStepBWWithoutHeating(photon_package * pp);

    void calcFinalTemperature(bool use_energy_density);
    void calcStochasticHeating();
    void calcAlignedRadii();

    bool isInvalid(double val)
    {
        if(val != val)
            return true;

        if(val == numeric_limits<double>::infinity())
            return true;

        return false;
    }

    void calcStepWidth(StokesVector & stokes_new,
                       StokesVector & stokes_new2,
                       double cell_d_l,
                       double * epsi,
                       double * dz_new)
    {
        double epsi_I = abs(stokes_new2.I() - stokes_new.I()) /
                        (REL_ERROR * abs(stokes_new.I()) + ABS_ERROR);

        double epsi_Q = abs(abs(stokes_new2.Q()) - abs(stokes_new.Q())) /
                        (REL_ERROR * abs(stokes_new.Q()) + ABS_ERROR);

        double epsi_U = abs(abs(stokes_new2.U()) - abs(stokes_new.U())) /
                        (REL_ERROR * abs(stokes_new.U()) + ABS_ERROR);

        double epsi_V = abs(abs(stokes_new2.V()) - abs(stokes_new.V())) /
                        (REL_ERROR * abs(stokes_new.V()) + ABS_ERROR);

        *epsi = max(epsi_I, max(epsi_Q, max(epsi_U, epsi_V)));
        *dz_new = 0.9 * cell_d_l * pow(*epsi, -0.2);
    }

    void convertTempInQB(double min_gas_density, bool use_gas_temp);

  private:
    string pathOutput;

    int * probing_points;

    uint nr_ray_detectors;
    uint nr_mc_detectors;

    CRaytracingBasic ** tracer;
    CGridBasic * grid;
    COpiateDataBase * op;
    CDustMixture * dust;
    CGasMixture * gas;
    slist sources_mc, sources_ray;
    CDetector * detector;

    double * RK_c;
    double * RK_b1;
    double * RK_b2;
    Matrix2D RK_a;

    double adjTgas;

    uint start, stop;

    bool b_forced;
    bool peel_off;
    bool mrw_step;
    bool stokes_dust_rad_field;

    uilist detector_wl_index;

    CSynchrotron * synchrotron;
};
