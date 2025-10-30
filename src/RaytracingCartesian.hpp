/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CRAYTRACING_CARTESIAN_H
#define CRAYTRACING_CARTESIAN_H

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

class CRaytracingCartesian : public CRaytracingBasic
{
public:
    CRaytracingCartesian(CGridBasic * _grid)
    {
        grid = _grid;
    }

    ~CRaytracingCartesian(void)
    {}

    bool setDustDetector(uint pos,
                         const parameters & param,
                         dlist dust_ray_detectors,
                         double _max_length,
                         string path);

    bool setSyncDetector(uint pos,
                         const parameters & param,
                         dlist sync_ray_detectors,
                         double _max_length,
                         string path);
    
    bool setFreeFreeDetector(uint pos,
                         const parameters & param,
                         dlist free_ray_detectors,
                         double _max_length,
                         string path);

    bool setLineDetector(uint pos,
                         const parameters & param,
                         dlist line_ray_detectors,
                         string path,
                         double _max_length, bool hasZeeman);

    bool setOPIATEDetector(uint pos,
                         const parameters & param,
                         dlist op_ray_detectors,
                         string path,
                         double _max_length);

    void getSubPixelCoordinates(uint subpixel_lvl,
                                double cx,
                                double cy,
                                int i_sub_x,
                                int i_sub_y,
                                double & subpix_cx,
                                double & subpix_cy);

    bool getUseSubpixel(double cx, double cy, uint subpixel_lvl);

    void addToDetector(photon_package * pp, int i_pix, bool direct = false);

    long getNpix();
};

#endif /* CRAYTRACING_CARTESIAN_H */
