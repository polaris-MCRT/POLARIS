/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CRAYTRACING_POLAR_H
#define CRAYTRACING_POLAR_H

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

class CRaytracingPolar : public CRaytracingBasic
{
public:
    CRaytracingPolar(CGridBasic * _grid)
    {
        npix_ph = 0;
        npix_r = 0;
        npix_total = 0;

        tmpStokes = 0;

        grid = _grid;
    }

    ~CRaytracingPolar(void)
    {
        if(npix_ph != 0)
            delete[] npix_ph;

        if(tmpStokes != 0)
        {
            for(uint i_spectral = 0; i_spectral < nr_extra * nr_spectral_bins; i_spectral++)
            {
                for(uint i_r = 0; i_r <= npix_r; i_r++)
                    delete[] tmpStokes[i_spectral][i_r];
                delete[] tmpStokes[i_spectral];
            }
            delete[] tmpStokes;
            tmpStokes = 0;
        }
    }

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

    bool setLineDetector(uint pos,
                         const parameters & param,
                         dlist line_ray_detectors,
                         string path,
                         double _max_length, bool hasZeeman);

    bool initPolarGridParameter();

    void initTmpStokes();

    bool getRelPosition(int64_t i_pix, double & cx, double & cy);

    void addToDetector(photon_package * pp, int64_t i_pix, bool direct = false);

    /*
    during post processing, the polar detector is mapped onto a cartesian detector

    Interpolation:
      - iterate over all cartesian pixels
      - for each cartesian pixel, take Stokes vector of nearest polar pixel and its neighboring polar pixels
      - interpolate to get resulting Stokes vector
      - good if cartesian detector has more pixels than the polar detector
      - if not, polar pixels (i.e. flux) will be ignored

    Nearest:
      - iterate over all polar pixels
      - for each polar pixel, take Stokes vector and add to corresponding cartesian pixel
      - good if cartesian detector has less pixels than the polar detector
      - if not, cartesian pixels potentially have zero flux
      - however, flux is maintained

    -> if the resulting cartesian detector has more pixels than there are polar pixels,
       then use the interpolation method
    */
    bool postProcessing();

    bool postProcessingUsingNearest();

    bool postProcessingUsingInterpolation();

    long getNpix();

private:
    void getCoordinateIDs(int64_t i_pix, uint & rID, uint & phID);

    double getRingElementArea(uint rID);

    dlist listR;
    uint * npix_ph;
    uint npix_r, npix_total;
    StokesVector *** tmpStokes;
};

#endif /* CRAYTRACING_POLAR_H */
