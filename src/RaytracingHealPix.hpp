/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CRAYTRACING_HEALPIX_H
#define CRAYTRACING_HEALPIX_H

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

class CRaytracingHealPix : public CRaytracingBasic
{
public:
    CRaytracingHealPix(CGridBasic * _grid)
    {
        nside = 1;
        npix = 12;

        sx = 0;
        sy = 0;
        sz = 0;

        vx = 0;
        vy = 0;
        vz = 0;

        l_min = 0;
        l_max = PIx2;

        b_min = 0;
        b_max = PI;

        grid = _grid;
    }

    ~CRaytracingHealPix(void)
    {}

    CRaytracingHealPix(int _nside)
    {
        nside = _nside;
        npix = 12 * nside * nside;
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

    void setOrientation(uint orientation_reference);

    long getNpix();

    double getMinArea();

    Vector3D getObserverVelocity();

    bool isNotAtCenter(photon_package * pp, double cx, double cy);

    void preparePhoton(photon_package * pp, double cx, double cy);

    void preparePhotonWithPosition(photon_package * pp, Vector3D pos, int & i_pix);

    void setDirection(photon_package * pp);

    void setPosition(Vector3D pos);

    bool getRelPosition(int i_pix, double & cx, double & cy);

    double getDistance();

    double getDistance(Vector3D pos);

    void addToDetector(photon_package * pp, int i_pix, bool direct = false);

    bool writeDustResults(uint ray_result_type);

    bool writeLineResults(CGasMixture * gas, uint i_species, uint i_line);

    bool writeSyncResults();
    
    bool writeFreeFreeResults();

    void setObserverPosition(Vector3D pos);

    void pix2ang_ring(int i_pix, double * theta, double * phi);

    void ang2ring_ring(double theta, double phi, int * i_pix);

private:
    static int isqrt(int v);

    static void pix2ang_ring_z_phi(int nside_, int pix, double * z, double * phi);

    static void ang2pix_ring_z_phi(int nside_, double z, double phi, int * pix);

    Vector3D det_pos;
    Vector3D detector_angle_offset;

    double sx, sy, sz;
    double vx, vy, vz;

    double l_min;
    double l_max;

    double b_min;
    double b_max;

    int nside;
    long npix;
};

#endif /* CRAYTRACING_HEALPIX_H */
