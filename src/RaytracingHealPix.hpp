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

static const int jrll[] = { 2,2,2,2,3,3,3,3,4,4,4,4 };
static const int jpll[] = { 1,3,5,7,0,2,4,6,1,3,5,7 };

static const int NB_DX[8] = {-1,-1, 0, 0, 0, 1, 1, 1};
static const int NB_DY[8] = {-1, 0,-1, 1, 1, 1, 0,-1};

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

        l_min = PI;
        l_max = PI;

        b_min = 0;
        b_max = PIx2;

        grid = _grid;
        
        is_patch = false;
        
        detector_proj = 0;
        
        proj_x = -1;
        proj_y = -1;
    }

    ~CRaytracingHealPix(void)
    {
        if(detector_proj!=0)
        {
            delete detector_proj;
            detector_proj = 0;
        }
    }

    CRaytracingHealPix(int64_t _nside)
    {
        nside = _nside;
        npix = 12 * nside * nside;
        
        is_patch = false;
        detector_proj = 0;
    }
    
    void initIndices();
    
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
    
    bool setDustAMEDetector(uint pos,
                             const parameters & param,
                             dlist ame_ray_detectors,
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

    void preparePhotonWithPosition(photon_package * pp, Vector3D pos, int64_t & i_pix);

    void setDirection(photon_package * pp);

    void setPosition(Vector3D pos);

    bool getRelPosition(int64_t i_pix, double & cx, double & cy);
    
    bool getRefPosition(int64_t i_pix, double & cx, double & cy);

    double getDistance();

    double getDistance(Vector3D pos);

    void addToDetector(photon_package * pp, int64_t i_pix, bool direct = false);

    bool writeDustResults(uint ray_result_type);

    bool writeLineResults(CGasMixture * gas, uint i_species, uint i_line);

    bool writeSyncResults();
    
    bool writeFreeFreeResults();
    
    bool writeDustAMEResults();

    void setObserverPosition(Vector3D pos);

    //void pix2ang_ring(int i_pix, double * theta, double * phi);

    //void ang2ring_ring(double theta, double phi, int * i_pix);
    

    void pix2ang_ring64(int64_t ipix, double *theta, double *phi);

    void ang2pix_ring64(double theta, double phi, int64_t *ipix);

    long isqrt64(int64_t v);
    
    int ring_neighbors_standalone64(int64_t ipix, int64_t out8[8]);
    
    double clamp_double(double x, double a, double b);
    
    double wrap2pi(double x);
    
    double ang_distance_to_pixel_center64(int64_t ipix, double theta, double phi);
    
    int pick_three_closest64(double theta, double phi, const int64_t *cands,
                                int ncand, int64_t out3[3]);
    
    void interpolate_ring_to_regular_grid_4pt(const double* hp_map,
                                          double* out, int Nx, int Ny,
                                          double theta0, double theta1,
                                          double phi0, double phi1);
    
    int collect_neighbor_candidates64(double theta, double phi, int64_t *cands, int maxcands);
    
    void ang2vec(double theta, double phi, Vector3D & v);
    
    void pixcenter_vec_ring64(int64_t ipix, Vector3D & v);
    
    bool projectHealMaps();
    
    int64_t getHealIndex(int64_t value);

    double fmodulo (double v1, double v2);

    int64_t imodulo64 (int64_t v1, int64_t v2);

    int64_t ang2pix_ring_z_phi64 (int64_t nside_, double z, double s,  double phi);

    void pix2ang_ring_z_phi64(int64_t pix, double *z, double *s, double *phi);    
    
    int64_t xyf2ring64 (int64_t nside_, int ix, int iy, int face_num);
    
    void ring2xyf64 (int64_t nside_, int64_t pix, int *ix, int *iy, int *face_num);
    
    int64_t special_div64 (int64_t a, int64_t b);

private:
    static int isqrt(int v);

    //static void pix2ang_ring_z_phi(int nside_, int pix, double * z, double * phi);

    //static void ang2pix_ring_z_phi(int nside_, double z, double phi, int * pix);

    Vector3D det_pos;
    Vector3D detector_angle_offset;

    double sx, sy, sz;
    double vx, vy, vz;

    double l_min;
    double l_max;

    double b_min;
    double b_max;

    int64_t nside;
    int64_t npix;
    
    bool is_patch;
    ilist64 heal_indices;
    dlist arr_theta;
    dlist arr_phi;
    
    long proj_x, proj_y;
    CDetector * detector_proj;
    
};

#endif /* CRAYTRACING_HEALPIX_H */
