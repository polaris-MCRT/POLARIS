/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CRAYTRACING_SLICE_H
#define CRAYTRACING_SLICE_H

#include "Detector.hpp"
#include "GasSpecies.hpp"
#include "GridBasic.hpp"
#include "MathFunctions.hpp"
#include "Matrix2D.hpp"
#include "Parameters.hpp"
#include "RaytracingBasic.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"
#include "Vector3D.hpp"

class CRaytracingSlice : public CRaytracingBasic
{
public:
    CRaytracingSlice(CGridBasic * _grid)
    {
        grid = _grid;
    }

    ~CRaytracingSlice(void)
    {}

    bool setDustDetector(uint pos,
                         const parameters & param,
                         dlist dust_ray_detectors,
                         double _max_length,
                         string path);

    bool getRelPositionMap(int64_t i_pix, double & cx, double & cy);

    bool setSyncDetector(uint pos,
                         const parameters & param,
                         dlist sync_ray_detectors,
                         double _max_length,
                         string path);

    bool setLineDetector(uint pos,
                         const parameters & param,
                         dlist line_ray_detectors,
                         string path,
                         double _max_length, bool hasZeeman);

    void preparePhoton(photon_package * pp, double cx, double cy);

    void resetPhotonPosition(photon_package * pp, int64_t i_pix);

    void addToDetector(photon_package * pp, int64_t i_pix, bool direct = false);

    long getNpix();
};

#endif /* CRAYTRACING_SLICE_H */
