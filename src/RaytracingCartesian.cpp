/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "RaytracingCartesian.hpp"

bool CRaytracingCartesian::setDustDetector(uint pos,
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
    
    /*(uint _detector_id,
              string _path,
              uint _bins_x,
              uint _bins_y,
              uint _id,
              double _sidelength_x,
              double _sidelength_y,
              double _map_shift_x,
              double _map_shift_y,
              double _distance,
              double _l_min,
              double _l_max,
              uint _nr_spectral_bins,
              uint _nr_extra,
              uint _special_param,
              uint _alignment,
              uint _fits_map_IDs)*/
    
    uint special_param=0;
    uint alignment= param.getAlignmentMechanism();
    uint fits_map_IDs=param.getMapIDs();


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
                                special_param,
                                alignment,
                                fits_map_IDs);
    
    detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

    return true;
}

bool CRaytracingCartesian::setSyncDetector(uint pos,
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

    //todo: rewrite for sync
    /*detector = new CDetector(rt_detector_shape,
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
    detector->setOrientation(n1, n2, rot_angle1, rot_angle2);*/

    return true;
}

bool CRaytracingCartesian::setFreeFreeDetector(uint pos,
                                           const parameters & param,
                                           dlist free_ray_detectors,
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

    double lam_min = free_ray_detectors[pos + 0];
    double lam_max = free_ray_detectors[pos + 1];
    nr_spectral_bins = uint(free_ray_detectors[pos + 2]);
    nr_extra = 1;

    sID = uint(free_ray_detectors[pos + 3]);

    rot_angle1 = PI / 180.0 * free_ray_detectors[pos + 4];
    rot_angle2 = PI / 180.0 * free_ray_detectors[pos + 5];

    distance = free_ray_detectors[pos + 6];

    sidelength_x = free_ray_detectors[pos + 7];
    sidelength_y = free_ray_detectors[pos + 8];

    max_length = _max_length;

    if(free_ray_detectors[pos + 9] != -1)
        map_shift_x = free_ray_detectors[pos + 9];
    if(free_ray_detectors[pos + 10] != -1)
        map_shift_y = free_ray_detectors[pos + 10];

    map_pixel_x = uint(free_ray_detectors[pos + NR_OF_RAY_DET - 2]);
    map_pixel_y = uint(free_ray_detectors[pos + NR_OF_RAY_DET - 1]);

    max_subpixel_lvl = param.getMaxSubpixelLvl();

    calcMapParameter();

    Vector3D n1 = param.getAxis1();
    Vector3D n2 = param.getAxis2();

    setDetCoordSystem(n1, n2);
    
    /*(uint _detector_id,
              string _path,
              uint _bins_x,
              uint _bins_y,
              uint _id,
              double _sidelength_x,
              double _sidelength_y,
              double _map_shift_x,
              double _map_shift_y,
              double _distance,
              double _l_min,
              double _l_max,
              uint _nr_spectral_bins,
              uint _nr_extra,
              uint _special_param,
              uint _alignment,
              uint _fits_map_IDs)*/
    
    uint special_param=3;
    uint fits_map_IDs=param.getMapIDs();
    uint alignment = 0;

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
                                nr_extra, special_param,
                                alignment, fits_map_IDs);
    
    detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

    return true;
}

bool CRaytracingCartesian::setDustAMEDetector(uint pos,
                             const parameters & param,
                             dlist ame_ray_detectors,
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

    double lam_min = ame_ray_detectors[pos + 0];
    double lam_max = ame_ray_detectors[pos + 1];
    nr_spectral_bins = uint(ame_ray_detectors[pos + 2]);
    nr_extra = 1;

    sID = uint(ame_ray_detectors[pos + 3]);

    rot_angle1 = PI / 180.0 * ame_ray_detectors[pos + 4];
    rot_angle2 = PI / 180.0 * ame_ray_detectors[pos + 5];

    distance = ame_ray_detectors[pos + 6];

    sidelength_x = ame_ray_detectors[pos + 7];
    sidelength_y = ame_ray_detectors[pos + 8];

    max_length = _max_length;

    if(ame_ray_detectors[pos + 9] != -1)
        map_shift_x = ame_ray_detectors[pos + 9];
    if(ame_ray_detectors[pos + 10] != -1)
        map_shift_y = ame_ray_detectors[pos + 10];

    map_pixel_x = uint(ame_ray_detectors[pos + NR_OF_RAY_DET - 2]);
    map_pixel_y = uint(ame_ray_detectors[pos + NR_OF_RAY_DET - 1]);

    max_subpixel_lvl = param.getMaxSubpixelLvl();

    calcMapParameter();

    Vector3D n1 = param.getAxis1();
    Vector3D n2 = param.getAxis2();

    setDetCoordSystem(n1, n2);
    
    /*    CDetector(uint _detector_id,
              string _path,
              uint _bins_x,
              uint _bins_y,
              uint _id,
              double _sidelength_x,
              double _sidelength_y,
              double _map_shift_x,
              double _map_shift_y,
              double _distance,
              double _l_min,
              double _l_max,
              uint _nr_spectral_bins,
              uint _nr_extra,
              uint _special_param,
              uint _alignment,
              uint _fits_map_IDs)*/
    
    uint special_param = 2;
    uint alignment = 0;
    uint fits_map_IDs = 0;

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
                                special_param,
                                alignment,
                                fits_map_IDs);
    
    detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

    return true;
}

bool CRaytracingCartesian::setLineDetector(uint pos,
                                           const parameters & param,
                                           dlist line_ray_detectors,
                                           string path,
                                           double _max_length, bool hasZeeman)
{
    rt_detector_shape = DET_PLANE;
    vel_maps_fits = param.getVelMapsFits();

    if(detector != 0)
    {
        delete detector;
        detector = 0;
    };

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
                                min_velocity,
                                max_velocity,
                                hasZeeman);
    detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

    return true;
}

bool CRaytracingCartesian::setOPIATEDetector(uint pos,
                                             const parameters & param,
                                             dlist op_ray_detectors,
                                             string path,
                                             double _max_length)
{
    rt_detector_shape = DET_PLANE;
    vel_maps_fits = param.getVelMapsFits();

    if(detector != 0)
    {
        delete detector;
        detector = 0;
    };

    dID = pos / NR_OF_OPIATE_DET;

    uint i_trans = -1;
    sID = uint(op_ray_detectors[pos]);
    double min_velocity = 0;
    double max_velocity = op_ray_detectors[pos + 1];

    double tmp_angle1 = op_ray_detectors[pos + 2];
    double tmp_angle2 = op_ray_detectors[pos + 3];

    rot_angle1 = PI / 180.0 * tmp_angle1;
    rot_angle2 = PI / 180.0 * tmp_angle2;

    distance = op_ray_detectors[pos + 4];

    sidelength_x = op_ray_detectors[pos + 5];
    sidelength_y = op_ray_detectors[pos + 6];

    max_length = _max_length;

    if(op_ray_detectors[pos + 7] != -1)
        map_shift_x = op_ray_detectors[pos + 7];
    if(op_ray_detectors[pos + 8] != -1)
        map_shift_y = op_ray_detectors[pos + 8];

    map_pixel_x = uint(op_ray_detectors[pos + NR_OF_OPIATE_DET - 3]);
    map_pixel_y = uint(op_ray_detectors[pos + NR_OF_OPIATE_DET - 2]);
    nr_spectral_bins = uint(op_ray_detectors[pos + NR_OF_OPIATE_DET - 1]);
    nr_extra = 1;

    max_subpixel_lvl = param.getMaxSubpixelLvl();

    calcMapParameter();

    Vector3D n1 = param.getAxis1();
    Vector3D n2 = param.getAxis2();

    setDetCoordSystem(n1, n2);

    /*detector = new CDetector(rt_detector_shape,
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
                                min_velocity,
                                max_velocity);
    detector->setOrientation(n1, n2, tmp_angle1, tmp_angle2);*/

    return true;
}

void CRaytracingCartesian::getSubPixelCoordinates(uint subpixel_lvl,
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

bool CRaytracingCartesian::getUseSubpixel(double cx, double cy, uint subpixel_lvl)
{
    // Init variables
    bool subpixel = false;

    if(subpixel_lvl < max_subpixel_lvl)
    {
        // Init vector of indices
        uilist ID_1;

        // Create new photon package with wavelength index wID and position it in
        // model
        photon_package pp = photon_package();

        preparePhoton(&pp, cx, cy);

        if(grid->findStartingPoint(&pp))
        {
            while(grid->next(&pp))
            {
                // Push the index of each cell along the path into ID_1
                ID_1.push_back(pp.getPositionCell()->getUniqueID());
            }
        }

        for(int i_sub_x = -1; i_sub_x <= 1 && subpixel == false; i_sub_x += 2)
        {
            for(int i_sub_y = -1; i_sub_y <= 1 && subpixel == false; i_sub_y += 2)
            {
                // Init vector of indices
                uilist ID_2;

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
                        // Push the index of each cell along the path into ID_2
                        ID_2.push_back(pp.getPositionCell()->getUniqueID());
                    }
                }

                // If any subpixel travel through different cells, perform subpixeling
                if(ID_1 != ID_2)
                {
                    subpixel = true;
                    break;
                }
            }
        }
    }
    else
        subpixel_warning = true;

    return subpixel;
}

void CRaytracingCartesian::addToDetector(photon_package * pp, int64_t i_pix, bool direct)
{
    pp->setDetectorProjection();

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins * nr_extra; i_spectral++)
    {
        // Set wavelength of photon package
        pp->setSpectralID(i_spectral);

        // Multiply by min area if such a multiplication did not happen before
        if(!direct)
            pp->getStokesVector(i_spectral)->multStokesParam(getMinArea());

        // Add photon Stokes vector to detector
        detector->addToRaytracingDetector(*pp);
        detector->addToRaytracingSedDetector(*pp);
    }
}

long CRaytracingCartesian::getNpix()
{
    return map_pixel_x * map_pixel_y;
}
