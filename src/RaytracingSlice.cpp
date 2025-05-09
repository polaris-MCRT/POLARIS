/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "RaytracingSlice.hpp"

bool CRaytracingSlice::setDustDetector(uint pos,
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
                                0,
                                0,
                                distance,
                                lam_min,
                                lam_max,
                                nr_spectral_bins,
                                nr_extra,
                                param.getAlignmentMechanism());
    detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

    return true;
}

bool CRaytracingSlice::getRelPositionMap(int i_pix, double & cx, double & cy)
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
        //cx += map_shift_x;
    }

    if(map_pixel_y % 2 == 1)
        cy = y * step_y;
    else
    {
        if(y > -1)
            y++;

        cy = y * step_y - CMathFunctions::sgn(y) * off_y;
        //cy += map_shift_y;
    }
    return true;
}

bool CRaytracingSlice::setSyncDetector(uint pos,
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
                                0,
                                0,
                                distance,
                                lam_min,
                                lam_max,
                                nr_spectral_bins,
                                nr_extra);
    detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

    return true;
}

bool CRaytracingSlice::setLineDetector(uint pos,
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
    double min_velocity = line_ray_detectors[pos + 2];
    double max_velocity = line_ray_detectors[pos + 3];

    rot_angle1 = PI / 180.0 * line_ray_detectors[pos + 4];
    rot_angle2 = PI / 180.0 * line_ray_detectors[pos + 5];

    distance = line_ray_detectors[pos + 6];
    sidelength_x = line_ray_detectors[pos + 7];
    sidelength_y = line_ray_detectors[pos + 8];

    max_length = _max_length;

    if(line_ray_detectors[pos + 9] != -1)
        map_shift_x = line_ray_detectors[pos + 9];
    if(line_ray_detectors[pos + 10] != -1)
        map_shift_y = line_ray_detectors[pos + 10];

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
                                0,
                                0,
                                distance,
                                i_trans,
                                nr_spectral_bins,
                                min_velocity,
                                max_velocity);
    detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

    return true;
}

void CRaytracingSlice::preparePhoton(photon_package * pp, double cx, double cy)
{
    Vector3D pos = cx * ex + cy * ez + map_shift_x *ey;
    pp->setPosition(pos);
    pp->setEX(ex);
    pp->setEY(ey);
    pp->setEZ(ez);
}

void CRaytracingSlice::resetPhotonPosition(photon_package * pp, int i_pix)
{
    double cx, cy;
    getRelPositionMap(i_pix, cx, cy);
    pp->setPosition(Vector3D(cx, cy, 0));
}

void CRaytracingSlice::addToDetector(photon_package * pp, int i_pix, bool direct)
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

long CRaytracingSlice::getNpix()
{
    return map_pixel_x * map_pixel_y;
}
