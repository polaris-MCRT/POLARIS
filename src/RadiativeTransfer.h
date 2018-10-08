#pragma once
#include "typedefs.h"
#include "OcTree.h"
#include "Grid.h"
#include "Dust.h"
#include "GasSpecies.h"
#include "MathFunctions.h"
#include "chelper.h"
#include "Source.h"
#include "Detector.h"
//#include "OPIATE.h"
#include "Raytracing.h"
#include "Synchrotron.h"

class CRadiativeTransfer
{
public:

    CRadiativeTransfer(parameter & param)
    {
        grid = 0;
        dust = 0;
        gas = 0;
        //op = 0;
        tracer = 0;
        b_forced = true;
        peel_off = false;
        mrw_step = false;

        start = MAX_UINT;
        stop = MAX_UINT;

        adjTgas = 0;

        probing_points = 0;
        detector = 0;
        nr_ofMCDetectors = 0;

        axis1.set(1, 0, 0);
        axis2.set(0, 1, 0);

        RK_c = 0;
        RK_b1 = 0;
        RK_b2 = 0;

        probing_points = 0;

        synchrotron=0;

        pathOutput = param.getPathOutput();
    }

    ~CRadiativeTransfer(void)
    {
        if(probing_points != 0)
            delete[] probing_points;

        if(RK_c != 0)
            delete[] RK_c;

        if(RK_b1 != 0)
            delete[] RK_b1;

        if(RK_b2 != 0)
            delete[] RK_b2;

        if(synchrotron!=0)
            delete synchrotron;
    }

    bool initiateDustRaytrace(parameter & param);
    bool initiateSyncRaytrace(parameter & param);
    bool initiateLineRaytrace(parameter & param);
    bool initiateOPIATE(parameter & param);
    bool initiateProbing(parameter & param);
    void initiateRadFieldMC(parameter & param);
    void initiateDustMC(parameter & param);

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
    bool calcMonteCarloRadiationField(uint command, bool use_energy_density, bool disable_reemission=false);
    // Set temperature (old!)
    bool setTemperatureDistribution();

    // Dust scattered light
    bool calcPolMapsViaMC();

    // Dust emission
    bool calcPolMapsViaRaytracing(parameter & param);
    void getDustPixelIntensity(CSourceBasic * tmp_source, double cx, double cy, uint subpixel_lvl, int pos_id);
    void getDustIntensity(photon_package * pp, CSourceBasic * tmp_source,
            double cx, double cy, uint subpixel_lvl);
    void calcStellarEmission();

    // Synchrontron emission
    bool calcSyncMapsViaRaytracing(parameter & param);
    void getSyncPixelIntensity(CSourceBasic * tmp_source, double cx, double cy, uint subpixel_lvl, int pos_id);
    void getSyncIntensity(photon_package * pp1, photon_package * pp2,photon_package * pp3, CSourceBasic * tmp_source,
            double cx, double cy, uint subpixel_lvl);

    // Line meission
    bool calcChMapsViaRaytracing(parameter & param);
    void getLinePixelIntensity(CSourceBasic * tmp_source, double cx, double cy,
            const uint i_species, const uint i_line, uint subpixel_lvl, int pos_id);
    void getLineIntensity(photon_package * pp, CSourceBasic * tmp_source, double cx, double cy,
            uint subpixel_lvl, const uint i_species, const uint i_line);

    // Calc radiation pressure
    //bool calcRadiativePressure(parameter & param);

    void updateRadiationField(photon_package * pp, double len)
    {
        double energy = len * pp->getStokesVector().I();
        grid->updateSpecLength(pp, energy);
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

    void setSourcesLists(slist & _sources_mc, slist & _sources_ray)
    {
        sources_mc = _sources_mc;
        sources_ray = _sources_ray;
    }

    void setDetectors(CDetector * d)
    {
        detector = d;
    }

    double getEscapeTauForced(photon_package * pp);

    bool doMRWStepBW(photon_package * pp);
    bool doMRWStepBWWithoutHeating(photon_package * pp);

    void calcFinalTemperature(bool use_energy_density);
    void calcStochasticHeating(bool update_temperature);
    void calcAlignedRadii();

    bool isInvalid(double val)
    {
        if(val != val)
            return true;

        if(val == numeric_limits<double>::infinity())
            return true;

        return false;
    }

    void convertTempInQB(double min_gas_density, bool use_gas_temp);

private:
    string pathOutput;

    int * probing_points;

    CRaytracingBasic * tracer;
    CGridBasic * grid;
    //COpiate * op;
    CDustMixture * dust;
    CGasMixture * gas;
    slist sources_mc, sources_ray;
    CDetector * detector;
    uint nr_ofMCDetectors;

    double * RK_c;
    double * RK_b1;
    double * RK_b2;
    Matrix2D RK_a;

    double adjTgas;

    uint start, stop;

    bool b_forced;
    bool peel_off;
    bool mrw_step;

    CSynchrotron * synchrotron;

    Vector3D axis1, axis2;
};
