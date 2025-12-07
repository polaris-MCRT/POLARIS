/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CPIPELINE_H
#define CPIPELINE_H

#include "CommandParser.hpp"
#include "Detector.hpp"
#include "DustMixture.hpp"
#include "GasSpecies.hpp"
#include "SourceBasic.hpp"
#include "Typedefs.hpp"
#include "Parameters.hpp"
#include "GridBasic.hpp"
#include "Vector3D.hpp"

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

    bool calcDustMapsViaRayTracing(parameters & param);
    
    bool calcFreeFreeMapsViaRayTracing(parameters & param);
    
    bool calcAMEMapsViaRayTracing(parameters & param);

    bool calcPolarizationMapsViaSynchrotron(parameters & param);

    bool calcChMapsViaRayTracing(parameters & param);

    bool calcPolarizationMapsViaMC(parameters & param);

    bool calcOpiateMapsViaRayTracing(parameters & param);

    bool proberobeLineOfSight(parameters & param);

    bool assignDustMixture(parameters & param, CDustMixture * dust, CGridBasic * grid);

    bool assignGasSpecies(parameters & param, CGasMixture * gas, CGridBasic * grid);

    // bool calcRadPressure(parameter & param);

    bool solveNextVelocityField(uint itID);

    // bool preparePressureData(CGridBasic * grid, CDustMixture * dust, parameters & param, bool plot, uint itID);

    bool createWavelengthList(parameters & param, CDustMixture * dust, CGasMixture * gas, COpiateDataBase * op);

    void printParameters(parameters & param, uint max_id);

    bool Init(int argc, char ** argv);

    void Finish();

    void Error();

    void Run();

    bool createOutputPaths(string path);

    bool assignGridType(CGridBasic *& grid, parameters & param);

    bool createPath(string path);

    void printConversionParameters(parameters & param);

    void printAdditionalParameters(parameters & param);

    void printPathParameters(parameters & param);

    void printPlotParameters(parameters & param, bool input_output = false);

    void printAlignmentParameters(parameters & param);

    void printSourceParameters(parameters & param, bool show_dust = false);

    void printDetectorParameters(parameters & param, bool monte_carlo = false);

    void printSynchrotronParameters(parameters & param);
    
    void printFreeFree(parameters & param);

    void deleteSourceLists();

    uint getNrOffsetEntriesRay(parameters & param, CDustMixture * dust, CGridBasic * grid);

    double getNrOfRayDetector(parameters & param);

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

#endif /* CPIPELINE_H */
