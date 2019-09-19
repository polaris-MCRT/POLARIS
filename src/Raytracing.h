#pragma once
#include "Detector.h"
#include "GasSpecies.h"
#include "Grid.h"
#include "MathFunctions.h"
#include "Matrix2D.h"
#include "Parameters.h"
#include "Stokes.h"
#include "Typedefs.h"
#include "Vector.h"

#ifndef RAYTRACING
#define RAYTRACING

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

        distance = 0;
        max_length = 0;

        off_len_x = 0;
        off_len_y = 0;

        vel_maps = false;

        detector = 0;
        grid = 0;
    }

    virtual ~CRaytracingBasic(void)
    {}

    virtual bool setDustDetector(uint pos,
                                 const parameters & param,
                                 dlist dust_ray_detectors,
                                 double _max_length,
                                 string path) = 0;

    // synchrotron detectors

    virtual bool setSyncDetector(uint pos,
                                 const parameters & param,
                                 dlist sync_ray_detectors,
                                 double _max_length,
                                 string path)
    {
        return false;
    }

    // end synchrotron

    virtual bool setLineDetector(uint pos,
                                 const parameters & param,
                                 dlist line_ray_detectors,
                                 string path,
                                 double _max_length)
    {
        return false;
    }

    double getWavelength(uint i_wave)
    {
        return detector->getWavelength(i_wave);
    }

    virtual long getNpix()
    {
        return 0;
    }

    virtual double getMinArea()
    {
        return (sidelength_x * sidelength_y) / double(map_pixel_x * map_pixel_y);
    }

    uint getSourceIndex()
    {
        return sID;
    }

    uint considerPointSources()
    {
        switch(rt_detector_shape)
        {
            case DET_PLANE:
            case DET_SPHER:
            case DET_POLAR:
                return true;
                break;

            default:
                return false;
                break;
        }
    }

    uint getNrSpectralBins()
    {
        return nr_spectral_bins;
    }

    virtual void preparePhoton(photon_package * pp, double cx, double cy)
    {
        Vector3D pos = cx * ex + cy * ey - max_length * ez;
        pp->setPosition(pos);
        pp->setEX(ex);
        pp->setEY(ey);
        pp->setEZ(ez);
    }

    virtual void preparePhotonWithPosition(photon_package * pp, Vector3D pos, int & i_pix)
    {
        pp->setPosition(pos);
        pp->setEX(ex);
        pp->setEY(ey);
        pp->setEZ(ez);

        i_pix = 0;
    }

    virtual void setCoordinateSystem(photon_package * pp)
    {
        pp->setEX(ex);
        pp->setEY(ey);
        pp->setEZ(ez);
    }

    virtual bool getRelPosition(int i_pix, double & cx, double & cy)
    {
        return getRelPositionMap(i_pix, cx, cy);
    }

    void calcMapParameter()
    {
        if(sidelength_x <= 0)
            sidelength_x = max_length;
        if(sidelength_y <= 0)
            sidelength_y = max_length;

        step_x = sidelength_x / double(map_pixel_x);
        off_x = step_x / 2.0;

        step_y = sidelength_y / double(map_pixel_y);
        off_y = step_y / 2.0;

        if(map_pixel_x % 2 == 0)
            off_len_x = int(0.5 * map_pixel_x);
        else
            off_len_x = int(0.5 * (map_pixel_x - 1));

        if(map_pixel_y % 2 == 0)
            off_len_y = int(0.5 * map_pixel_y);
        else
            off_len_y = int(0.5 * (map_pixel_y - 1));
    }

    bool getRelPositionMap(int i_pix, double & cx, double & cy)
    {
        int y = (i_pix % map_pixel_y);
        int x = i_pix / map_pixel_y - off_len_x;
        y -= off_len_y;

        if(map_pixel_x % 2 == 1)
            cx = x * step_x;
        else
        {
            if(x > -1)
                x++;

            cx = x * step_x - CMathFunctions::sgn(x) * off_x;
            cx += map_shift_x;
        }

        if(map_pixel_y % 2 == 1)
            cy = y * step_y;
        else
        {
            if(y > -1)
                y++;

            cy = y * step_y - CMathFunctions::sgn(y) * off_y;
            cy += map_shift_y;
        }
        return true;
    }

    virtual void setDetCoordSystem(const Vector3D & n1, const Vector3D & n2)
    {
        ex.set(1, 0, 0);
        ey.set(0, 1, 0);
        ez.set(0, 0, 1);

        double cos_a = cos(rot_angle1);
        double sin_a = sin(rot_angle1);

        ex.rot(n1, cos_a, sin_a);
        ey.rot(n1, cos_a, sin_a);
        ez.rot(n1, cos_a, sin_a);

        cos_a = cos(rot_angle2);
        sin_a = sin(rot_angle2);

        ex.rot(n2, cos_a, sin_a);
        ey.rot(n2, cos_a, sin_a);
        ez.rot(n2, cos_a, sin_a);

        ex.normalize();
        ey.normalize();
        ez.normalize();
    }

    bool splitDustEmission()
    {
        return split_emission;
    }

    uint getNrExtra()
    {
        return nr_extra;
    }

    virtual double getDistance()
    {
        return distance;
    }

    virtual double getDistance(Vector3D pos)
    {
        double proj_length = -ez * pos;
        return proj_length + getDistance();
    }

    virtual void addToDetector(photon_package * pp, int i_pix, bool direct = false)
    {}

    virtual void setObserverPosition(Vector3D pos)
    {}

    virtual bool postProcessing()
    {
        return true;
    }

    virtual bool writeDustResults(uint ray_result_type)
    {
        if(!detector->writeMap(dID, ray_result_type))
            return false;

        if(!detector->writeSed(dID, ray_result_type))
            return false;

        return true;
    }

    virtual bool writeLineResults(CGasMixture * gas, uint i_species, uint i_line)
    {
        if(vel_maps)
            if(!detector->writeVelChannelMaps(gas, i_species, i_line))
                return false;

        if(!detector->writeIntChannelMaps(gas, i_species, i_line))
            return false;

        if(!detector->writeLineSpectrum(gas, i_species, i_line))
            return false;

        return true;
    }

    virtual bool writeSyncResults()
    {
        if(!detector->writeSyncMap(dID))
            return false;

        return true;
    }

    virtual bool getUseSubpixel(double cx, double cy, uint subpixel_lvl)
    {
        return false;
    }

    double getDistanceFactor()
    {
        return 1 / (getDistance() * getDistance());
    }

    double getDistanceFactor(Vector3D pos)
    {
        return 1 / (getDistance(pos) * (getDistance(pos)));
    }

    virtual void getDetectorData(dlist & C, dlist & T)
    {}

    virtual void setPosition(Vector3D pos){};

    virtual bool isNotAtCenter(photon_package * pp, double cx, double cy)
    {
        return true;
    }

    virtual Vector3D getObserverVelocity()
    {
        return Vector3D(0, 0, 0);
    }

    virtual void getSubPixelCoordinates(uint subpixel_lvl,
                                        double cx,
                                        double cy,
                                        int i_sub_x,
                                        int i_sub_y,
                                        double & subpix_cx,
                                        double & subpix_cy)
    {
        subpix_cx = cx;
        subpix_cy = cy;
    }

    double getLamMin()
    {
        return detector->getLamMin();
    }

    double getLamMax()
    {
        return detector->getLamMax();
    }

    double getChannelWidth()
    {
        return detector->getChannelWidth();
    }

    uint getNrOfSpectralBins()
    {
        return detector->getNrOfSpectralBins();
    }

    double getVelocityChannel(uint vch)
    {
        return detector->getVelocityChannel(vch);
    }

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

    int off_len_x, off_len_y;

    bool vel_maps, split_emission;
    CDetector * detector;
    CGridBasic * grid;
    Vector3D ex, ey, ez;
};

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
                         string path)
    {
        rt_detector_shape = DET_PLANE;

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }

        dID = pos / NR_OF_RAY_DET;

        split_emission = param.splitDustEmission();

        double lam_min = dust_ray_detectors[pos + 0];
        double lam_max = dust_ray_detectors[pos + 1];
        nr_spectral_bins = uint(dust_ray_detectors[pos + 2]);
        nr_extra = (split_emission ? 4 : 1);

        sID = uint(dust_ray_detectors[pos + 3]);

        rot_angle1 = PI / 180.0 * dust_ray_detectors[pos + 4];
        rot_angle2 = PI / 180.0 * dust_ray_detectors[pos + 5];

        distance = dust_ray_detectors[pos + 6];

        sidelength_x = dust_ray_detectors[pos + 7];
        sidelength_y = dust_ray_detectors[pos + 8];

        max_length = _max_length;

        if(dust_ray_detectors[pos + 9] != -1)
            map_shift_x = dust_ray_detectors[pos + 9];
        if(dust_ray_detectors[pos + 10] != -1)
            map_shift_y = dust_ray_detectors[pos + 10];

        map_pixel_x = uint(dust_ray_detectors[pos + NR_OF_RAY_DET - 2]);
        map_pixel_y = uint(dust_ray_detectors[pos + NR_OF_RAY_DET - 1]);

        max_subpixel_lvl = param.getMaxSubpixelLvl();

        calcMapParameter();

        Vector3D n1 = param.getAxis1();
        Vector3D n2 = param.getAxis2();

        setDetCoordSystem(n1, n2);

        detector = new CDetector(rt_detector_shape,
                                 path,
                                 map_pixel_x,
                                 map_pixel_y,
                                 dID,
                                 sidelength_x,
                                 sidelength_y,
                                 map_shift_x,
                                 map_shift_y,
                                 distance,
                                 lam_min,
                                 lam_max,
                                 nr_spectral_bins,
                                 nr_extra,
                                 param.getAlignmentMechanism());
        detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

        return true;
    }

    bool setSyncDetector(uint pos,
                         const parameters & param,
                         dlist sync_ray_detectors,
                         double _max_length,
                         string path)
    {
        rt_detector_shape = DET_PLANE;

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }

        dID = pos / NR_OF_RAY_DET;

        double lam_min = sync_ray_detectors[pos + 0];
        double lam_max = sync_ray_detectors[pos + 1];
        nr_spectral_bins = uint(sync_ray_detectors[pos + 2]);
        nr_extra = 2;

        sID = uint(sync_ray_detectors[pos + 3]);

        rot_angle1 = PI / 180.0 * sync_ray_detectors[pos + 4];
        rot_angle2 = PI / 180.0 * sync_ray_detectors[pos + 5];

        distance = sync_ray_detectors[pos + 6];

        sidelength_x = sync_ray_detectors[pos + 7];
        sidelength_y = sync_ray_detectors[pos + 8];

        max_length = _max_length;

        if(sync_ray_detectors[pos + 9] != -1)
            map_shift_x = sync_ray_detectors[pos + 9];
        if(sync_ray_detectors[pos + 10] != -1)
            map_shift_y = sync_ray_detectors[pos + 10];

        map_pixel_x = uint(sync_ray_detectors[pos + NR_OF_RAY_DET - 2]);
        map_pixel_y = uint(sync_ray_detectors[pos + NR_OF_RAY_DET - 1]);

        max_subpixel_lvl = param.getMaxSubpixelLvl();

        calcMapParameter();

        Vector3D n1 = param.getAxis1();
        Vector3D n2 = param.getAxis2();

        setDetCoordSystem(n1, n2);

        detector = new CDetector(rt_detector_shape,
                                 path,
                                 map_pixel_x,
                                 map_pixel_y,
                                 dID,
                                 sidelength_x,
                                 sidelength_y,
                                 map_shift_x,
                                 map_shift_y,
                                 distance,
                                 lam_min,
                                 lam_max,
                                 nr_spectral_bins,
                                 nr_extra);
        detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

        return true;
    }

    bool setLineDetector(uint pos,
                         const parameters & param,
                         dlist line_ray_detectors,
                         string path,
                         double _max_length)
    {
        rt_detector_shape = DET_PLANE;
        vel_maps = param.getVelMaps();

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        };

        dID = pos / NR_OF_LINE_DET;

        uint i_trans = uint(line_ray_detectors[pos + 0]);
        sID = uint(line_ray_detectors[pos + 1]);
        double max_velocity = line_ray_detectors[pos + 2];

        rot_angle1 = PI / 180.0 * line_ray_detectors[pos + 3];
        rot_angle2 = PI / 180.0 * line_ray_detectors[pos + 4];

        distance = line_ray_detectors[pos + 5];

        sidelength_x = line_ray_detectors[pos + 6];
        sidelength_y = line_ray_detectors[pos + 7];

        max_length = _max_length;

        if(line_ray_detectors[pos + 8] != -1)
            map_shift_x = line_ray_detectors[pos + 8];
        if(line_ray_detectors[pos + 9] != -1)
            map_shift_y = line_ray_detectors[pos + 9];

        map_pixel_x = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 3]);
        map_pixel_y = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 2]);
        nr_spectral_bins = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 1]);
        nr_extra = 1;

        max_subpixel_lvl = param.getMaxSubpixelLvl();

        calcMapParameter();

        Vector3D n1 = param.getAxis1();
        Vector3D n2 = param.getAxis2();

        setDetCoordSystem(n1, n2);

        detector = new CDetector(rt_detector_shape,
                                 path,
                                 map_pixel_x,
                                 map_pixel_y,
                                 dID,
                                 sidelength_x,
                                 sidelength_y,
                                 map_shift_x,
                                 map_shift_y,
                                 distance,
                                 i_trans,
                                 nr_spectral_bins,
                                 max_velocity);
        detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

        return true;
    }

    void getSubPixelCoordinates(uint subpixel_lvl,
                                double cx,
                                double cy,
                                int i_sub_x,
                                int i_sub_y,
                                double & subpix_cx,
                                double & subpix_cy)
    {
        double subpixel_fraction_1D = pow(2.0, -double(subpixel_lvl + 1));
        subpix_cx = cx + i_sub_x * off_x * subpixel_fraction_1D;
        subpix_cy = cy + i_sub_y * off_y * subpixel_fraction_1D;
    }

    bool getUseSubpixel(double cx, double cy, uint subpixel_lvl)
    {
        // Init variables
        bool subpixel = false;

        if(subpixel_lvl < max_subpixel_lvl)
        {
            // Init sum of indices
            ulong ID_sum_1 = 0;

            // Create new photon package with wavelength index wID and position it in
            // model
            photon_package pp = photon_package();

            preparePhoton(&pp, cx, cy);

            if(grid->findStartingPoint(&pp))
            {
                while(grid->next(&pp))
                {
                    // Add the index of each cell along the path onto ID_sum_1
                    ID_sum_1 += pp.getPositionCell()->getUniqueID();
                }
            }

            // Delete photon package after usage
            for(int i_sub_x = -1; i_sub_x <= 1 && subpixel == false; i_sub_x += 2)
            {
                for(int i_sub_y = -1; i_sub_y <= 1 && subpixel == false; i_sub_y += 2)
                {
                    // Init sum of indices
                    ulong ID_sum_2 = 0;

                    // Calculate positions of each subpixel
                    double tmp_cx, tmp_cy;
                    getSubPixelCoordinates(subpixel_lvl, cx, cy, i_sub_x, i_sub_y, tmp_cx, tmp_cy);

                    // Create new photon package with wavelength index wID and position it
                    // in model
                    preparePhoton(&pp, tmp_cx, tmp_cy);

                    // Find starting point of the photon package
                    // and transport it through the model
                    if(grid->findStartingPoint(&pp))
                    {
                        while(grid->next(&pp))
                        {
                            // Add the index of each cell along the path onto ID_sum_2
                            ID_sum_2 += pp.getPositionCell()->getUniqueID();
                        }
                    }

                    // If any subpixel travel through other cells, perform subpixelling
                    if(ID_sum_1 != ID_sum_2)
                    {
                        subpixel = true;
                        break;
                    }
                }
            }
        }
        return subpixel;
    }

    void addToDetector(photon_package * pp, int i_pix, bool direct = false)
    {
        pp->setDetectorProjection();
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins * nr_extra; i_spectral++)
        {
            // Set wavelength of photon package
            pp->setSpectralID(i_spectral);

            // Multiply by min area if such a multiplication did not happen before
            if(!direct)
                pp->getStokesVector(i_spectral)->multS(getMinArea());

            // Add photon Stokes vector to detector
            detector->addToRaytracingDetector(*pp);
            detector->addToRaytracingSedDetector(*pp);
        }
    }

    long getNpix()
    {
        return map_pixel_x * map_pixel_y;
    }
};

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
                         string path)
    {
        rt_detector_shape = DET_SPHER;

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }

        dID = pos / NR_OF_RAY_DET;

        split_emission = param.splitDustEmission();

        double lam_min = dust_ray_detectors[pos + 0];
        double lam_max = dust_ray_detectors[pos + 1];
        nr_spectral_bins = uint(dust_ray_detectors[pos + 2]);
        nr_extra = (split_emission ? 4 : 1);

        sID = uint(dust_ray_detectors[pos + 3]);

        sx = dust_ray_detectors[pos + 4];
        sy = dust_ray_detectors[pos + 5];
        sz = dust_ray_detectors[pos + 6];

        l_min = PI * (-dust_ray_detectors[pos + 8] + 180.0) / 180.0;
        l_max = PI * (-dust_ray_detectors[pos + 7] + 180.0) / 180.0;
        b_min = PI * (-dust_ray_detectors[pos + 10] + 90.0) / 180;
        b_max = PI * (-dust_ray_detectors[pos + 9] + 90.0) / 180;

        nside = uint(dust_ray_detectors[pos + NR_OF_RAY_DET - 1]);

        npix = 12 * nside * nside;

        det_pos.setX(sx);
        det_pos.setY(sy);
        det_pos.setZ(sz);

        max_length = _max_length * 10;

        setOrientation(param.getHealpixOrientation());

        detector = new CDetector(path,
                                 npix,
                                 1,
                                 det_pos,
                                 max_length,
                                 lam_min,
                                 lam_max,
                                 nr_spectral_bins,
                                 nr_extra,
                                 param.getAlignmentMechanism());
        detector->setObsPosition(Vector3D(sx, sy, sz), Vector3D(0, 0, 0), l_min, l_max, b_min, b_max);

        return true;
    }

    bool setSyncDetector(uint pos,
                         const parameters & param,
                         dlist sync_ray_detectors,
                         double _max_length,
                         string path)
    {
        rt_detector_shape = DET_SPHER;

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }

        dID = pos / NR_OF_RAY_DET;

        double lam_min = sync_ray_detectors[pos + 0];
        double lam_max = sync_ray_detectors[pos + 1];
        nr_spectral_bins = uint(sync_ray_detectors[pos + 2]);
        nr_extra = 2;

        sID = uint(sync_ray_detectors[pos + 3]);

        sx = sync_ray_detectors[pos + 4];
        sy = sync_ray_detectors[pos + 5];
        sz = sync_ray_detectors[pos + 6];

        l_min = PI * (-sync_ray_detectors[pos + 8] + 180.0) / 180.0;
        l_max = PI * (-sync_ray_detectors[pos + 7] + 180.0) / 180.0;
        b_min = PI * (-sync_ray_detectors[pos + 10] + 90.0) / 180;
        b_max = PI * (-sync_ray_detectors[pos + 9] + 90.0) / 180;

        nside = uint(sync_ray_detectors[pos + NR_OF_RAY_DET - 1]);

        npix = 12 * nside * nside;

        det_pos.setX(sx);
        det_pos.setY(sy);
        det_pos.setZ(sz);

        max_length = _max_length * 10;

        setOrientation(param.getHealpixOrientation());

        detector =
            new CDetector(path, npix, 1, det_pos, max_length, lam_min, lam_max, nr_spectral_bins, nr_extra);
        detector->setObsPosition(Vector3D(sx, sy, sz), Vector3D(0, 0, 0), l_min, l_max, b_min, b_max);

        return true;
    }

    bool setLineDetector(uint pos,
                         const parameters & param,
                         dlist line_ray_detectors,
                         string path,
                         double _max_length)
    {
        rt_detector_shape = DET_SPHER;
        vel_maps = param.getVelMaps();

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }

        dID = pos / NR_OF_LINE_DET;

        uint i_trans = uint(line_ray_detectors[pos + 0]);
        sID = uint(line_ray_detectors[pos + 1]);
        double max_velocity = line_ray_detectors[pos + 2];

        sx = line_ray_detectors[pos + 3];
        sy = line_ray_detectors[pos + 4];
        sz = line_ray_detectors[pos + 5];

        l_min = PI * (-line_ray_detectors[pos + 7] + 180.0) / 180.0;
        l_max = PI * (-line_ray_detectors[pos + 6] + 180.0) / 180.0;
        b_min = PI * (-line_ray_detectors[pos + 9] + 90.0) / 180;
        b_max = PI * (-line_ray_detectors[pos + 8] + 90.0) / 180;

        setOrientation(param.getHealpixOrientation());

        vx = line_ray_detectors[pos + 10];
        vy = line_ray_detectors[pos + 11];
        vz = line_ray_detectors[pos + 12];

        nside = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 2]);
        nr_spectral_bins = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 1]);
        nr_extra = 1;

        npix = 12 * nside * nside;

        max_length = _max_length * 10;

        det_pos.setX(sx);
        det_pos.setY(sy);
        det_pos.setZ(sz);
        detector = new CDetector(path, npix, 1, det_pos, max_length, i_trans, nr_spectral_bins, max_velocity);
        detector->setObsPosition(Vector3D(sx, sy, sz), Vector3D(vx, vy, vz), l_min, l_max, b_min, b_max);

        return true;
    }

    void setOrientation(uint orientation_reference)
    {
        if(orientation_reference == HEALPIX_CENTER)
            detector_angle_offset = det_pos.getSphericalCoord();
        else if(orientation_reference == HEALPIX_YAXIS)
        {
            detector_angle_offset = det_pos.getSphericalCoord();
            detector_angle_offset.setTheta(PI2);
        }
        else
        {
            detector_angle_offset.setPhi(-PI);
            detector_angle_offset.setTheta(PI2);
        }
    }

    long getNpix()
    {
        return npix;
    }

    double getMinArea()
    {
        return 4 * PI / double(npix);
    }

    Vector3D getObserverVelocity()
    {
        return Vector3D(vx, vy, vz);
    }

    bool isNotAtCenter(photon_package * pp, double cx, double cy)
    {
        double theta = cx;
        double phi = cy;

        /*if(theta <(45-40)*PI/180)
            return false;

        if(theta > (45+40)*PI/180)
            return false;

        if(phi <-40*PI/180)
            return false;

        if(phi > 40*PI/180)
            return false;*/

        Vector3D ph_dir = pp->getDirection();
        Vector3D ph_pos = pp->getPosition();
        Vector3D new_dir = det_pos - ph_pos;

        double lam = ph_dir * new_dir;

        if(lam <= 0)
            return false;

        return true;
    }

    void preparePhoton(photon_package * pp, double cx, double cy)
    {
        double theta = cx + (PI2 - detector_angle_offset.Theta());
        double phi = cy + PI + detector_angle_offset.Phi();

        Vector3D start_pos, tmp_ex, tmp_ey, tmp_ez;

        tmp_ez.setX(sin(theta) * cos(phi));
        tmp_ez.setY(sin(theta) * sin(phi));
        tmp_ez.setZ(cos(theta));

        tmp_ey.setX(cos(theta) * cos(phi));
        tmp_ey.setY(cos(theta) * sin(phi));
        tmp_ey.setZ(-sin(theta));

        tmp_ex.setX(-sin(phi));
        tmp_ex.setY(cos(phi));
        tmp_ex.setZ(0);

        start_pos += max_length * tmp_ez + det_pos;

        pp->setPosition(start_pos);
        pp->setEX(tmp_ex);
        pp->setEY(-tmp_ey);
        pp->setEZ(-tmp_ez);
    }

    void preparePhotonWithPosition(photon_package * pp, Vector3D pos, int & i_pix)
    {
        pp->setPosition(pos);

        Vector3D tmp_ex, tmp_ey;
        Vector3D tmp_ez = pos - det_pos;
        tmp_ez.normalize();

        double theta = acos(tmp_ez.Z()) - (PI2 - detector_angle_offset.Theta());
        double phi = atan3(tmp_ez.X(), tmp_ez.Y()) - (PI + detector_angle_offset.Phi());

        tmp_ey.setX(cos(theta) * cos(phi));
        tmp_ey.setY(cos(theta) * sin(phi));
        tmp_ey.setZ(-sin(theta));

        tmp_ex.setX(-sin(phi));
        tmp_ex.setY(cos(phi));
        tmp_ex.setZ(0);

        ang2ring_ring(theta, phi, &i_pix);

        pp->setEX(tmp_ex);
        pp->setEY(-tmp_ey);
        pp->setEZ(-tmp_ez);
    }

    void setDirection(photon_package * pp)
    {
        Vector3D dir;
        Vector3D pos = pp->getPosition();

        dir = (pos - det_pos) / max_length;
        dir.normalize();

        pp->setEZ(dir);
    }

    void setPosition(Vector3D pos)
    {
        sx = pos.X();
        sy = pos.Y();
        sz = pos.Z();

        det_pos.setX(sx);
        det_pos.setY(sy);
        det_pos.setZ(sz);
    }

    bool getRelPosition(int i_pix, double & cx, double & cy)
    {
        pix2ang_ring(i_pix, &cx, &cy);

        if(cx < b_min)
            return false;

        if(cx > b_max)
            return false;

        if(cy < l_min)
            return false;

        if(cy > l_max)
            return false;

        return true;
    }

    double getDistance()
    {
        return 1.0;
    }

    double getDistance(Vector3D pos)
    {
        Vector3D pos_obs = Vector3D(sx, sy, sz);
        Vector3D diff = (pos - pos_obs);
        return diff.length();
    }

    void addToDetector(photon_package * pp, int i_pix, bool direct = false)
    {
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins * nr_extra; i_spectral++)
        {
            // Set wavelength of photon package
            pp->setSpectralID(i_spectral);

            // Multiply by min area if such a multiplication did not happen before
            if(!direct)
                pp->getStokesVector(i_spectral)->multS(getMinArea());

            // Add photon Stokes vector to detector
            detector->addToRaytracingDetector(*pp, i_pix);
            detector->addToRaytracingSedDetector(*pp);
        }
    }

    bool writeDustResults(uint ray_result_type)
    {
        if(!detector->writeHealMaps(dID, ray_result_type))
            return false;

        if(!detector->writeSed(dID, ray_result_type))
            return false;

        return true;
    }

    bool writeLineResults(CGasMixture * gas, uint i_species, uint i_line)
    {
        if(vel_maps)
            if(!detector->writeVelChannelHealMaps(gas, i_species, i_line))
                return false;

        if(!detector->writeIntVelChannelHealMaps(gas, i_species, i_line))
            return false;

        if(!detector->writeLineSpectrum(gas, i_species, i_line))
            return false;

        return true;
    }

    bool writeSyncResults()
    {
        if(!detector->writeSyncHealMap(dID))
            return false;

        return true;
    }

    void setObserverPosition(Vector3D pos)
    {
        sx = pos.X();
        sy = pos.Y();
        sz = pos.Z();
    }

    void pix2ang_ring(int i_pix, double * theta, double * phi)
    {
        double z;
        pix2ang_ring_z_phi(nside, i_pix, &z, phi);
        *theta = acos(z);
    }

    void ang2ring_ring(double theta, double phi, int * i_pix)
    {
        double z = cos(theta);
        ang2pix_ring_z_phi(nside, z, phi, i_pix);
    }

  private:
    static int isqrt(int v)
    {
        return (int)(sqrt(v + 0.5));
    }

    static void pix2ang_ring_z_phi(int nside_, int pix, double * z, double * phi)
    {
        long ncap_ = nside_ * (nside_ - 1) * 2;
        long npix_ = 12 * nside_ * nside_;
        double fact2_ = 4. / npix_;
        if(pix < ncap_) /* North Polar cap */
        {
            int iring = (1 + isqrt(1 + 2 * pix)) >> 1; /* counted from North pole */
            int iphi = (pix + 1) - 2 * iring * (iring - 1);

            *z = 1.0 - (iring * iring) * fact2_;
            *phi = (iphi - 0.5) * PI2 / iring;
        }
        else if(pix < (npix_ - ncap_)) /* Equatorial region */
        {
            double fact1_ = (nside_ << 1) * fact2_;
            int ip = pix - ncap_;
            int iring = ip / (4 * nside_) + nside_; /* counted from North pole */
            int iphi = ip % (4 * nside_) + 1;
            /* 1 if iring+nside is odd, 1/2 otherwise */
            double fodd = ((iring + nside_) & 1) ? 1 : 0.5;

            int nl2 = 2 * nside_;
            *z = (nl2 - iring) * fact1_;
            *phi = (iphi - fodd) * PI / nl2;
        }
        else /* South Polar cap */
        {
            int ip = npix_ - pix;
            int iring = (1 + isqrt(2 * ip - 1)) >> 1; /* counted from South pole */
            int iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));

            *z = -1.0 + (iring * iring) * fact2_;
            *phi = (iphi - 0.5) * PI2 / iring;
        }
    }

    static void ang2pix_ring_z_phi(int nside_, double z, double phi, int * pix)
    {
        long ncap_ = nside_ * (nside_ - 1) * 2;
        long npix_ = 12 * nside_ * nside_;
        double za = abs(z);
        double tt = CMathFunctions::fmodulo(phi, PIx2) * invPI2; // in [0,4)
        if(za <= TWOTHIRD)                                       // Equatorial region
        {
            double temp1 = nside_ * (0.5 + tt);
            double temp2 = nside_ * z * 0.75;
            int jp = (int)(temp1 - temp2); /* index of  ascending edge line */
            int jm = (int)(temp1 + temp2); /* index of descending edge line */

            /* ring number counted from z=2/3 */
            int ir = nside_ + 1 + jp - jm; /* in {1,2n+1} */
            int kshift = 1 - (ir & 1);     /* kshift=1 if ir even, 0 otherwise */

            int ip = (jp + jm - nside_ + kshift + 1) / 2; /* in {0,4n-1} */
            ip = CMathFunctions::imodulo(ip, 4 * nside_);

            *pix = nside_ * (nside_ - 1) * 2 + (ir - 1) * 4 * nside_ + ip;
        }
        else // North & South polar caps
        {
            double tp = tt - (int)(tt);
            double tmp = nside_ * sqrt(3 * (1 - za));

            int jp = (int)(tp * tmp);         /* increasing edge line index */
            int jm = (int)((1.0 - tp) * tmp); /* decreasing edge line index */

            int ir = jp + jm + 1;    /* ring number counted from the closest pole */
            int ip = (int)(tt * ir); /* in {0,4*ir-1} */
            ip = CMathFunctions::imodulo(ip, 4 * ir);

            if(z > 0)
                *pix = 2 * ir * (ir - 1) + ip;
            else
                *pix = 12 * nside_ * nside_ - 2 * ir * (ir + 1) + ip;
        }
    }

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
                         string path)
    {
        rt_detector_shape = DET_POLAR;

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }

        dID = pos / NR_OF_RAY_DET;

        split_emission = param.splitDustEmission();

        double lam_min = dust_ray_detectors[pos + 0];
        double lam_max = dust_ray_detectors[pos + 1];
        nr_spectral_bins = uint(dust_ray_detectors[pos + 2]);
        nr_extra = (split_emission ? 4 : 1);

        sID = uint(dust_ray_detectors[pos + 3]);

        rot_angle1 = PI / 180.0 * dust_ray_detectors[pos + 4];
        rot_angle2 = PI / 180.0 * dust_ray_detectors[pos + 5];

        distance = dust_ray_detectors[pos + 6];

        sidelength_x = dust_ray_detectors[pos + 7];
        sidelength_y = dust_ray_detectors[pos + 8];

        max_length = _max_length;

        map_pixel_x = uint(dust_ray_detectors[pos + NR_OF_RAY_DET - 2]);
        map_pixel_y = uint(dust_ray_detectors[pos + NR_OF_RAY_DET - 1]);

        calcMapParameter();

        Vector3D n1 = param.getAxis1();
        Vector3D n2 = param.getAxis2();

        setDetCoordSystem(n1, n2);

        max_subpixel_lvl = param.getMaxSubpixelLvl();

        if(!initPolarGridParameter())
            return false;

        initTmpStokes();

        detector = new CDetector(rt_detector_shape,
                                 path,
                                 map_pixel_x,
                                 map_pixel_y,
                                 dID,
                                 sidelength_x,
                                 sidelength_y,
                                 0.0,
                                 0.0,
                                 distance,
                                 lam_min,
                                 lam_max,
                                 nr_spectral_bins,
                                 nr_extra,
                                 param.getAlignmentMechanism());
        detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

        return true;
    }

    bool setSyncDetector(uint pos,
                         const parameters & param,
                         dlist sync_ray_detectors,
                         double _max_length,
                         string path)
    {
        rt_detector_shape = DET_POLAR;

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }

        dID = pos / NR_OF_RAY_DET;

        double lam_min = sync_ray_detectors[pos + 0];
        double lam_max = sync_ray_detectors[pos + 1];
        nr_spectral_bins = uint(sync_ray_detectors[pos + 2]);
        nr_extra = 2;

        sID = uint(sync_ray_detectors[pos + 3]);

        rot_angle1 = PI / 180.0 * sync_ray_detectors[pos + 4];
        rot_angle2 = PI / 180.0 * sync_ray_detectors[pos + 5];

        distance = sync_ray_detectors[pos + 6];

        sidelength_x = sync_ray_detectors[pos + 7];
        sidelength_y = sync_ray_detectors[pos + 8];

        max_length = _max_length;

        map_pixel_x = uint(sync_ray_detectors[pos + NR_OF_RAY_DET - 2]);
        map_pixel_y = uint(sync_ray_detectors[pos + NR_OF_RAY_DET - 1]);

        calcMapParameter();

        Vector3D n1 = param.getAxis1();
        Vector3D n2 = param.getAxis2();

        setDetCoordSystem(n1, n2);

        max_subpixel_lvl = param.getMaxSubpixelLvl();

        if(!initPolarGridParameter())
            return false;

        initTmpStokes();

        detector = new CDetector(rt_detector_shape,
                                 path,
                                 map_pixel_x,
                                 map_pixel_y,
                                 dID,
                                 sidelength_x,
                                 sidelength_y,
                                 0.0,
                                 0.0,
                                 distance,
                                 lam_min,
                                 lam_max,
                                 nr_spectral_bins,
                                 nr_extra);
        detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

        return true;
    }

    bool setLineDetector(uint pos,
                         const parameters & param,
                         dlist line_ray_detectors,
                         string path,
                         double _max_length)
    {
        rt_detector_shape = DET_POLAR;
        vel_maps = param.getVelMaps();

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }

        dID = pos / NR_OF_LINE_DET;

        uint i_trans = uint(line_ray_detectors[pos + 0]);
        sID = uint(line_ray_detectors[pos + 1]);
        double max_velocity = line_ray_detectors[pos + 2];

        rot_angle1 = PI / 180.0 * line_ray_detectors[pos + 3];
        rot_angle2 = PI / 180.0 * line_ray_detectors[pos + 4];

        distance = line_ray_detectors[pos + 5];

        sidelength_x = line_ray_detectors[pos + 6];
        sidelength_y = line_ray_detectors[pos + 7];

        max_length = _max_length;

        map_pixel_x = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 3]);
        map_pixel_y = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 2]);
        nr_spectral_bins = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 1]);
        nr_extra = 1;

        calcMapParameter();

        Vector3D n1 = param.getAxis1();
        Vector3D n2 = param.getAxis2();

        setDetCoordSystem(n1, n2);

        max_subpixel_lvl = param.getMaxSubpixelLvl();

        if(!initPolarGridParameter())
            return false;

        initTmpStokes();

        detector = new CDetector(rt_detector_shape,
                                 path,
                                 map_pixel_x,
                                 map_pixel_y,
                                 dID,
                                 sidelength_x,
                                 sidelength_y,
                                 0,
                                 0,
                                 distance,
                                 i_trans,
                                 nr_spectral_bins,
                                 max_velocity);
        detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

        return true;
    }

    bool initPolarGridParameter()
    {
        double pixel_width = max(sidelength_x / map_pixel_x, sidelength_y / map_pixel_y);
        double max_len = 0.5 * sqrt(sidelength_x * sidelength_x + sidelength_y * sidelength_y);
        if(!grid->getPolarRTGridParameter(max_len, pixel_width, max_subpixel_lvl, listR, npix_r, npix_ph))
        {
            cout << "\nERROR: Polar detector can only be used with spherical or "
                    "cylindrical grids!";
            return false;
        }

        for(uint i_r = 0; i_r < npix_r; i_r++)
        {
            // Make sure that the number of phi cells is even
            if(npix_ph[i_r] % 2 == 1)
                (npix_ph[i_r])++;

            // Add the amount of phi cells in the current ring to the total amount
            npix_total += npix_ph[i_r];
        }

        npix_total++;

        if(npix_total > MAX_RT_RAYS)
        {
            cout << "\nHINT: Very high amount of rays required for DUST EMISSION "
                    "simulation with polar raytracing grid!"
                 << endl
                 << "      Problem: The simulation may take a long time to process all "
                    "rays."
                 << endl
                 << "      Solutions: Decrease max subpixel level or use the cartesian "
                    "raytracing grid."
                 << endl;
        }

        return true;
    }

    void initTmpStokes()
    {
        tmpStokes = new StokesVector **[nr_extra * nr_spectral_bins];
        for(uint i_spectral = 0; i_spectral < nr_extra * nr_spectral_bins; i_spectral++)
        {
            tmpStokes[i_spectral] = new StokesVector *[npix_r + 1];
            for(uint i_r = 0; i_r < npix_r; i_r++)
                tmpStokes[i_spectral][i_r] = new StokesVector[npix_ph[i_r]];

            // Add center pixel
            tmpStokes[i_spectral][npix_r] = new StokesVector[1];
        }
    }

    bool getRelPosition(int i_pix, double & cx, double & cy)
    {
        // Return (0,0) coordinate, if central pixel
        if(i_pix == npix_total - 1)
        {
            cx = 0;
            cy = 0;
            return true;
        }

        int rID, phID;
        getCoordinateIDs(i_pix, rID, phID);

        double phi = (phID + 0.5) * PIx2 / double(npix_ph[rID]);

        double r1 = listR[rID];
        double r2 = listR[rID + 1];

        double radius = (r1 + r2) / 2.0;

        cx = radius * cos(phi);
        cy = radius * sin(phi);

        return true;
    }

    void addToDetector(photon_package * pp, int i_pix, bool direct = false)
    {
        if(direct)
        {
            pp->setDetectorProjection();
            for(uint i_spectral = 0; i_spectral < nr_spectral_bins * nr_extra; i_spectral++)
            {
                // Set wavelength of photon package
                pp->setSpectralID(i_spectral);

                // Add photon Stokes vector to detector
                detector->addToRaytracingDetector(*pp);
                detector->addToRaytracingSedDetector(*pp);
            }
        }
        else
        {
            int rID, phID;
            getCoordinateIDs(i_pix, rID, phID);
            for(uint i_spectral = 0; i_spectral < nr_spectral_bins * nr_extra; i_spectral++)
            {
                // Set wavelength of photon package
                pp->setSpectralID(i_spectral);

                // Set Stokes vector of polar ring element
                tmpStokes[pp->getSpectralID()][rID][phID] = *pp->getStokesVector();

                // Update Stokes vector of photon package
                pp->getStokesVector()->multS(getRingElementArea(rID));

                // Add photon Stokes vector SED detector
                detector->addToRaytracingSedDetector(*pp);
            }
        }
    }

    bool postProcessing()
    {
        // Init counter and percentage to show progress
        ullong per_counter = 0;
        float last_percentage = 0;

        // Create a double list with the center points of the rt grid in radial direction
        dlist r_center_pos;
        for(uint i_r = 0; i_r < npix_r; i_r++)
            r_center_pos.push_back(0.5 * (listR[i_r] + listR[i_r + 1]));

#pragma omp parallel for schedule(dynamic)
        for(int i_pix = 0; i_pix < int(map_pixel_x * map_pixel_y); i_pix++)
        {
            photon_package pp = photon_package(nr_spectral_bins * nr_extra);

            // Init variables
            Vector3D pos;
            StokesVector Q1, Q2, Q3, Q4;
            int rID, rID1, rID2;
            int phID, phID1, phID2, phID3, phID4;
            double r1, r2, p1, p2, p3, p4, dph1, dph2;
            double cx = 0, cy = 0;

            // Increase counter used to show progress
            per_counter++;

            // Calculate percentage of total progress per source
            float percentage =
                100.0 * float(per_counter) / float(map_pixel_x * map_pixel_y * nr_spectral_bins);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
#pragma omp critical
                {
                    cout << "-> Interpolating from polar grid to detector map: " << percentage
                         << " [%]            \r" << flush;
                    last_percentage = percentage;
                }
            }

            // Get cartesian coordinates of the current pixel
            if(!getRelPositionMap(i_pix, cx, cy))
                continue;

            // Set photon package position and get R and Phi
            pos = Vector3D(cx, cy, 0);
            pp.setPosition(pos);
            pos.cart2cyl();

            // Check if central pixel
            if(cx == 0 && cy == 0)
            {
                for(uint i_spectral = 0; i_spectral < nr_spectral_bins * nr_extra; i_spectral++)
                {
                    // Set wavelength of photon package
                    pp.setSpectralID(i_spectral);

                    // Get Stokes Vector
                    pp.setStokesVector(tmpStokes[i_spectral][npix_r][0] * getMinArea(), i_spectral);

                    // Transport pixel value to detector
                    detector->addToRaytracingDetector(pp);
                }
                continue;
            }

            // Get radius index from center position list (subtract first outer border)
            rID = CMathFunctions::biListIndexSearch(pos.R(), r_center_pos);

            // Calculate the radial indices and positions
            if(rID != MAX_UINT)
            {
                rID1 = rID;
                rID2 = rID + 1;
                r1 = r_center_pos[rID1];
                r2 = r_center_pos[rID2];
            }
            else if(pos.R() > r_center_pos[npix_r - 1] && pos.R() <= listR[npix_r])
            {
                rID1 = npix_r - 1;
                rID2 = npix_r - 1;
                r1 = r_center_pos[rID1];
                r2 = listR[npix_r];
            }
            else
                continue;

            // Calculate the difference between two phi cells
            dph1 = PIx2 / double(npix_ph[rID1]);
            dph2 = PIx2 / double(npix_ph[rID2]);

            // Calculate the phi indices
            phID1 = floor(pos.Phi() / dph1 - 0.5);
            phID2 = floor(pos.Phi() / dph1 + 0.5);
            phID3 = floor(pos.Phi() / dph2 - 0.5);
            phID4 = floor(pos.Phi() / dph2 + 0.5);

            // Calculate the phi positions
            p1 = phID1 * dph1;
            p2 = (phID2 + 1.0) * dph1;
            p3 = phID3 * dph2;
            p4 = (phID4 + 1.0) * dph2;

            // Modify phi index if phi is close to 0
            if(phID1 < 0)
                phID1 += npix_ph[rID1];
            if(phID3 < 0)
                phID3 += npix_ph[rID2];

            // Modify phi index if phi is close to PIx2
            if(phID2 >= npix_ph[rID1])
                phID2 -= npix_ph[rID1];
            if(phID4 >= npix_ph[rID2])
                phID4 -= npix_ph[rID2];

            for(uint i_spectral = 0; i_spectral < nr_spectral_bins * nr_extra; i_spectral++)
            {
                // Set wavelength of photon package
                pp.setSpectralID(i_spectral);

                if(rID1 == npix_r - 1)
                {
                    // If outside of last ring center interpolate only from one side
                    // (leave stokes empty)
                    Q1 = tmpStokes[i_spectral][rID1][phID1];
                    Q2 = tmpStokes[i_spectral][rID1][phID2];
                    Q3.clear();
                    Q4.clear();
                }
                else
                {
                    // If between two ring center positions consider all points
                    Q1 = tmpStokes[i_spectral][rID1][phID1];
                    Q2 = tmpStokes[i_spectral][rID1][phID2];
                    Q3 = tmpStokes[i_spectral][rID2][phID3];
                    Q4 = tmpStokes[i_spectral][rID2][phID4];
                }

                // Calculate relative position factors
                double t = (pos.R() - r1) / (r2 - r1);
                double s = (pos.Phi() - p1 - (p3 - p1) * t) / (p2 + (p4 - p2) * t - p1 - (p3 - p1) * t);

                // Do 2D linear interpolation
                double stokes_I = CMathFunctions::getPolarInterp(t, s, Q1.I(), Q2.I(), Q3.I(), Q4.I());
                double stokes_Q = CMathFunctions::getPolarInterp(t, s, Q1.Q(), Q2.Q(), Q3.Q(), Q4.Q());
                double stokes_U = CMathFunctions::getPolarInterp(t, s, Q1.U(), Q2.U(), Q3.U(), Q4.U());
                double stokes_V = CMathFunctions::getPolarInterp(t, s, Q1.V(), Q2.V(), Q3.V(), Q4.V());
                double stokes_T = CMathFunctions::getPolarInterp(t, s, Q1.T(), Q2.T(), Q3.T(), Q4.T());
                double stokes_Sp = CMathFunctions::getPolarInterp(t, s, Q1.Sp(), Q2.Sp(), Q3.Sp(), Q4.Sp());

                // Add pixel value to photon_package
                StokesVector res_stokes =
                    StokesVector(stokes_I, stokes_Q, stokes_U, stokes_V, stokes_T, stokes_Sp);
                pp.setStokesVector(res_stokes * getMinArea());
                // Transport pixel value to detector
                detector->addToRaytracingDetector(pp);
            }
        }

        cout << "-> Interpolating from polar grid to detector map: 100 [%]         \r" << flush;

        return true;
    }

    long getNpix()
    {
        return npix_total;
    }

  private:
    void getCoordinateIDs(int i_pix, int & rID, int & phID)
    {
        rID = 0;
        phID = 0;
        if(i_pix == npix_total - 1)
        {
            rID = npix_r;
            return;
        }

        while(true)
        {
            if(i_pix < npix_ph[rID])
            {
                phID = i_pix;
                return;
            }
            else
            {
                i_pix -= npix_ph[rID];
                rID++;
            }
        }
    }

    double getRingElementArea(int rID)
    {
        if(rID == npix_r)
            return PI * listR[0] * listR[0];

        double r1 = listR[rID];
        double r2 = listR[rID + 1];

        return PI * (r2 * r2 - r1 * r1) / npix_ph[rID];
    }

    dlist listR;
    uint * npix_ph;
    uint npix_r, npix_total;
    StokesVector *** tmpStokes;
};

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
                         string path)
    {
        rt_detector_shape = DET_SLICE;

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }

        dID = pos / NR_OF_RAY_DET;

        split_emission = param.splitDustEmission();

        double lam_min = dust_ray_detectors[pos + 0];
        double lam_max = dust_ray_detectors[pos + 1];
        nr_spectral_bins = uint(dust_ray_detectors[pos + 2]);
        nr_extra = (split_emission ? 4 : 1);

        sID = uint(dust_ray_detectors[pos + 3]);

        rot_angle1 = PI / 180.0 * dust_ray_detectors[pos + 4];
        rot_angle2 = PI / 180.0 * dust_ray_detectors[pos + 5];

        distance = dust_ray_detectors[pos + 6];

        sidelength_x = dust_ray_detectors[pos + 7];
        sidelength_y = dust_ray_detectors[pos + 8];

        max_length = _max_length;

        if(dust_ray_detectors[pos + 9] != -1)
            map_shift_x = dust_ray_detectors[pos + 9];
        if(dust_ray_detectors[pos + 10] != -1)
            map_shift_y = dust_ray_detectors[pos + 10];

        map_pixel_x = uint(dust_ray_detectors[pos + NR_OF_RAY_DET - 2]);
        map_pixel_y = uint(dust_ray_detectors[pos + NR_OF_RAY_DET - 1]);

        calcMapParameter();

        Vector3D n1 = param.getAxis1();
        Vector3D n2 = param.getAxis2();

        setDetCoordSystem(n1, n2);

        detector = new CDetector(rt_detector_shape,
                                 path,
                                 map_pixel_x,
                                 map_pixel_y,
                                 dID,
                                 sidelength_x,
                                 sidelength_y,
                                 map_shift_x,
                                 map_shift_y,
                                 distance,
                                 lam_min,
                                 lam_max,
                                 nr_spectral_bins,
                                 nr_extra,
                                 param.getAlignmentMechanism());
        detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

        return true;
    }

    bool setSyncDetector(uint pos,
                         const parameters & param,
                         dlist sync_ray_detectors,
                         double _max_length,
                         string path)
    {
        rt_detector_shape = DET_SLICE;

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }

        dID = pos / NR_OF_RAY_DET;

        double lam_min = sync_ray_detectors[pos + 0];
        double lam_max = sync_ray_detectors[pos + 1];
        nr_spectral_bins = uint(sync_ray_detectors[pos + 2]);
        nr_extra = 2;

        sID = uint(sync_ray_detectors[pos + 3]);

        rot_angle1 = PI / 180.0 * sync_ray_detectors[pos + 4];
        rot_angle2 = PI / 180.0 * sync_ray_detectors[pos + 5];

        distance = sync_ray_detectors[pos + 6];

        sidelength_x = sync_ray_detectors[pos + 7];
        sidelength_y = sync_ray_detectors[pos + 8];

        max_length = _max_length;

        if(sync_ray_detectors[pos + 9] != -1)
            map_shift_x = sync_ray_detectors[pos + 9];
        if(sync_ray_detectors[pos + 10] != -1)
            map_shift_y = sync_ray_detectors[pos + 10];

        map_pixel_x = uint(sync_ray_detectors[pos + NR_OF_RAY_DET - 2]);
        map_pixel_y = uint(sync_ray_detectors[pos + NR_OF_RAY_DET - 1]);

        calcMapParameter();

        Vector3D n1 = param.getAxis1();
        Vector3D n2 = param.getAxis2();

        setDetCoordSystem(n1, n2);

        detector = new CDetector(rt_detector_shape,
                                 path,
                                 map_pixel_x,
                                 map_pixel_y,
                                 dID,
                                 sidelength_x,
                                 sidelength_y,
                                 map_shift_x,
                                 map_shift_y,
                                 distance,
                                 lam_min,
                                 lam_max,
                                 nr_spectral_bins,
                                 nr_extra);
        detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

        return true;
    }

    bool setLineDetector(uint pos,
                         const parameters & param,
                         dlist line_ray_detectors,
                         string path,
                         double _max_length)
    {
        rt_detector_shape = DET_SLICE;
        vel_maps = param.getVelMaps();

        if(detector != 0)
        {
            delete detector;
            detector = 0;
        }

        dID = pos / NR_OF_LINE_DET;

        uint i_trans = uint(line_ray_detectors[pos + 0]);
        sID = uint(line_ray_detectors[pos + 1]);
        double max_velocity = line_ray_detectors[pos + 2];

        rot_angle1 = PI / 180.0 * line_ray_detectors[pos + 3];
        rot_angle2 = PI / 180.0 * line_ray_detectors[pos + 4];

        distance = line_ray_detectors[pos + 5];
        sidelength_x = line_ray_detectors[pos + 6];
        sidelength_y = line_ray_detectors[pos + 7];

        max_length = _max_length;

        if(line_ray_detectors[pos + 8] != -1)
            map_shift_x = line_ray_detectors[pos + 8];
        if(line_ray_detectors[pos + 9] != -1)
            map_shift_y = line_ray_detectors[pos + 9];

        map_pixel_x = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 3]);
        map_pixel_y = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 2]);
        nr_spectral_bins = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 1]);
        nr_extra = 1;

        calcMapParameter();

        Vector3D n1 = param.getAxis1();
        Vector3D n2 = param.getAxis2();

        setDetCoordSystem(n1, n2);

        detector = new CDetector(rt_detector_shape,
                                 path,
                                 map_pixel_x,
                                 map_pixel_y,
                                 dID,
                                 sidelength_x,
                                 sidelength_y,
                                 map_shift_x,
                                 map_shift_y,
                                 distance,
                                 i_trans,
                                 nr_spectral_bins,
                                 max_velocity);
        detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

        return true;
    }

    void preparePhoton(photon_package * pp, double cx, double cy)
    {
        Vector3D pos = cx * ex + cy * ez;
        pp->setPosition(pos);
        pp->setEX(ex);
        pp->setEY(ey);
        pp->setEZ(ez);
    }

    void resetPhotonPosition(photon_package * pp, int i_pix)
    {
        double cx, cy;
        getRelPositionMap(i_pix, cx, cy);
        pp->setPosition(Vector3D(cx, cy, 0));
    }

    void addToDetector(photon_package * pp, int i_pix, bool direct = false)
    {
        resetPhotonPosition(pp, i_pix);
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins * nr_extra; i_spectral++)
        {
            // Set wavelength of photon package
            pp->setSpectralID(i_spectral);

            // Multiply by min area if such a multiplication did not happen before
            if(!direct)
                pp->getStokesVector(i_spectral)->multS(getMinArea());

            // Add photon Stokes vector to detector
            detector->addToRaytracingDetector(*pp);
            detector->addToRaytracingSedDetector(*pp);
        }
    }

    long getNpix()
    {
        return map_pixel_x * map_pixel_y;
    }
};

#endif
