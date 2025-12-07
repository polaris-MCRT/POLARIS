/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CRAYTRACING_BASIC_H
#define CRAYTRACING_BASIC_H

#include "Detector.hpp"
#include "GasSpecies.hpp"
#include "GridBasic.hpp"
#include "MathFunctions.hpp"
#include "Matrix2D.hpp"
#include "Parameters.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"
#include "Vector3D.hpp"

class CRaytracingBasic
{
public:
    CRaytracingBasic()
    {
        rt_detector_shape = 0;

        dID = 0;
        sID = 0;

        map_pixel_x = 0;
        map_pixel_y = 0;

        nr_spectral_bins = 0;
        max_subpixel_lvl = 0;
        nr_extra = 0;

        subpixel_warning = false;

        sidelength_x = 0;
        sidelength_y = 0;

        step_x = 0;
        step_y = 0;

        off_x = 0;
        off_y = 0;

        map_shift_x = 0;
        map_shift_y = 0;

        rot_angle1 = 0;
        rot_angle2 = 0;

        rad_bubble = 0;

        distance = 0;
        max_length = 0;

        off_len_x = 0;
        off_len_y = 0;

        vel_maps_fits = false;
        heal_type = 0;
        project_healpix = false;

        detector = 0;
        grid = 0;
    }

    virtual ~CRaytracingBasic(void)
    {
        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }
    }

    virtual bool setDustDetector(uint pos,
                                 const parameters & param,
                                 dlist dust_ray_detectors,
                                 double _max_length,
                                 string path) = 0;

    virtual bool setSyncDetector(uint pos,
                                 const parameters & param,
                                 dlist sync_ray_detectors,
                                 double _max_length,
                                 string path);
    
    virtual bool setFreeFreeDetector(uint pos,
                                 const parameters & param,
                                 dlist free_ray_detectors,
                                 double _max_length,
                                 string path);
    
    virtual bool setDustAMEDetector(uint pos,
                                 const parameters & param,
                                 dlist ame_ray_detectors,
                                 double _max_length,
                                 string path);

    virtual bool setLineDetector(uint pos,
                                 const parameters & param,
                                 dlist line_ray_detectors,
                                 string path,
                                 double _max_length, bool hasZeeman);

    virtual bool setOPIATEDetector(uint pos,
                                 const parameters & param,
                                 dlist op_ray_detectors,
                                 string path,
                                 double _max_length);

    double getWavelength(uint i_wave);

    virtual long getNpix();

    virtual double getMinArea();

    uint getSourceIndex();

    uint considerPointSources();

    uint getNrSpectralBins();

    virtual void preparePhoton(photon_package * pp, double cx, double cy);

    virtual void preparePhotonWithPosition(photon_package * pp, Vector3D pos, int64_t & i_pix);

    virtual void setCoordinateSystem(photon_package * pp);

    virtual bool getRelPosition(int64_t i_pix, double & cx, double & cy);

    void calcMapParameter();

    virtual bool getRelPositionMap(int64_t i_pix, double & cx, double & cy);

    virtual void setDetCoordSystem(const Vector3D & n1, const Vector3D & n2);

    bool splitDustEmission();

    uint getNrExtra();

    uint getDetectorShape();

    bool getSubpixelWarning();

    virtual double getDistance();

    virtual double getDistance(Vector3D pos);

    virtual void addToDetector(photon_package * pp, int64_t i_pix, bool direct = false);

    virtual void addToDetector(photon_package * pp1, photon_package * pp2, int64_t i_pix, bool direct = false);

    virtual void setObserverPosition(Vector3D pos);

    virtual bool postProcessing();

    virtual bool writeDustResults(uint ray_result_type);

    virtual bool writeLineResults(CGasMixture * gas, uint i_species, uint i_line);

    virtual bool writeOpiateResults(COpiateDataBase * op);

    virtual bool writeSyncResults();
    
    virtual bool writeFreeFreeResults();
    
    virtual bool writeDustAMEResults();

    virtual bool getUseSubpixel(double cx, double cy, uint subpixel_lvl);

    double getDistanceFactor();

    double getDistanceFactor(Vector3D pos);

    virtual void getDetectorData(dlist & C, dlist & T);

    virtual void setPosition(Vector3D pos);

    virtual bool isNotAtCenter(photon_package * pp, double cx, double cy);

    virtual Vector3D getObserverVelocity();

    virtual void getSubPixelCoordinates(uint subpixel_lvl,
                                        double cx,
                                        double cy,
                                        int i_sub_x,
                                        int i_sub_y,
                                        double & subpix_cx,
                                        double & subpix_cy);

    double getLamMin();

    double getLamMax();

    double getChannelWidth();

    uint getNrOfSpectralBins();

    double getVelocityChannel(uint vch);

protected:
    uint rt_detector_shape;
    uint dID, sID;
    uint map_pixel_x, map_pixel_y;
    uint nr_spectral_bins;
    uint max_subpixel_lvl;
    uint nr_extra;

    double sidelength_x, sidelength_y;
    double step_x, step_y;
    double off_x, off_y;
    double map_shift_x, map_shift_y;
    double distance, max_length;
    double rot_angle1, rot_angle2;
    double rad_bubble;

    bool subpixel_warning;

    int off_len_x, off_len_y;

    bool vel_maps_fits;
    uint heal_type;
    bool project_healpix;
    bool split_emission;
    CDetector * detector;
    CGridBasic * grid;
    Vector3D ex, ey, ez;
};

#endif /* CRAYTRACING_BASIC_H */
