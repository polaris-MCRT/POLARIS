/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CRADIATIVE_TRANSFER_H
#define CRADIATIVE_TRANSFER_H

#include "Detector.hpp"
#include "DustMixture.hpp"
#include "GasSpecies.hpp"
#include "GridBasic.hpp"
#include "MathFunctions.hpp"
#include "GridOcTree.hpp"
#include "Photon.hpp"
#include "SourceBasic.hpp"
#include "Typedefs.hpp"
#include "RaytracingBasic.hpp"
#include "RaytracingCartesian.hpp"
#include "RaytracingHealPix.hpp"
#include "RaytracingPolar.hpp"
#include "RaytracingSlice.hpp"
#include "Synchrotron.hpp"
#include "FreeFree.hpp"
#include "Matrix2D.hpp"
#include "Parameters.hpp"
#include "Stokes.hpp"
#include "Vector3D.hpp"
#include "FreeFree.hpp"

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
        freefree = 0;

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
    
    bool initiateFreeFreeRaytrace(parameters & param);

    bool initiateLineRaytrace(parameters & param); 

    bool initiateOPIATERaytrace(parameters & param);

    bool initiateProbing(parameters & param);

    void initiateRadFieldMC(parameters & param);

    void initiateDustMC(parameters & param);

    void initiateRungeKuttaFehlberg();

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

    //free free emission 
    bool calcFreeFreeMapsViaRaytracing(parameters & param);
    
    void getFreeFreePixelIntensity(CSourceBasic * tmp_source,
                               double cx,
                               double cy,
                               uint i_det,
                               uint subpixel_lvl,
                               int pos_id);
    
    void getFreeFreeIntensity(photon_package * pp1,
                          CSourceBasic * tmp_source,
                          double cx,
                          double cy,
                          uint i_det,
                          uint subpixel_lvl);
                        
    void rayThroughCellFreeFree(photon_package * pp1, uint i_det, uint nr_used_wavelengths);
    
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

    void updateRadiationField(photon_package * pp);

    void setGrid(CGridBasic * _grid);

    void setDust(CDustMixture * _dust);
    
    void setFreeFree(CFreeFree * _free);

    void setGas(CGasMixture * _gas);

    void setOpiateDataBase(COpiateDataBase * _op);

    void setSourcesLists(slist & _sources_mc, slist & _sources_ray);

    void setDetectors(CDetector * d);

    double getOpticalDepthAlongPath(photon_package * pp);

    bool photonInDetectorDir(photon_package * pp, CDetector * detector);

    void scaleAddToDetector(photon_package * pp, CDetector * detector, ullong interactions);

    bool doMRWStepBW(photon_package * pp);

    bool doMRWStepBWWithoutHeating(photon_package * pp);

    void calcFinalTemperature(bool use_energy_density);

    void calcStochasticHeating();

    void calcAlignedRadii();

    bool isInvalid(double val);

    void calcStepWidth(StokesVector & stokes_new,
                       StokesVector & stokes_new2,
                       double cell_d_l,
                       double * epsi,
                       double * dz_new);

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
    CFreeFree * freefree;
};

#endif /* CRADIATIVE_TRANSFER_H */
