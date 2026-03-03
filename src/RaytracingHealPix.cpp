/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "RaytracingHealPix.hpp"

bool CRaytracingHealPix::setDustDetector(uint pos,
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
    
    heal_type = param.getHealType();
    project_healpix = param.projectHealMaps();
    
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

    double tmp_l_min = dust_ray_detectors[pos + 7];
    double tmp_l_max = dust_ray_detectors[pos + 8];

    double tmp_b_min = dust_ray_detectors[pos + 9];
    double tmp_b_max = dust_ray_detectors[pos + 10];
    
    if( dust_ray_detectors[pos + 7] != -180 ||
        dust_ray_detectors[pos + 8] !=  180 ||
            dust_ray_detectors[pos + 9]  != -90 ||
            dust_ray_detectors[pos + 10] !=  90)
    {
        is_patch=true;
        
        //l_min = PI * (-dust_ray_detectors[pos + 8] + 180.0) / 180.0;
        //l_max = PI * (-dust_ray_detectors[pos + 7] + 180.0) / 180.0;

        l_min = PI * (-dust_ray_detectors[pos + 7]) / 180.0;
        l_max = PI * (-dust_ray_detectors[pos + 8] + 360.0) / 180.0;

        b_min = PI * (-dust_ray_detectors[pos + 10] + 90.0) / 180;
        b_max = PI * (-dust_ray_detectors[pos + 9] + 90.0) / 180;        
    }
    else
    {
        is_patch=false;
        
        l_min=PI;
        l_max=PI;
        
        b_min=0;
        b_max=PI;
    }
    
    if(dust_ray_detectors[pos + 11]>0)
        rad_bubble = dust_ray_detectors[pos + 11];

    nside = int64_t(dust_ray_detectors[pos + NR_OF_RAY_DET - 1]);

    npix = 12 * nside * nside;

    det_pos.setX(sx);
    det_pos.setY(sy);
    det_pos.setZ(sz);

    max_length = _max_length * 10;
    uint alignment= param.getAlignmentMechanism();
    uint fits_map_IDs=param.getMapIDs();

    setOrientation(param.getHealpixOrientation());
    


    initIndices();
    
    int64_t tmp_bins=heal_indices.size();
            
    //detector = new CDetector(
    //    path, npix, 1, det_pos, max_length, lam_min, lam_max, rad_bubble, nr_spectral_bins, nr_extra, 1,alignment,compact_fits);
    
    detector = new CDetector(
        path, tmp_bins, 1, det_pos, max_length, lam_min, lam_max, rad_bubble, nr_spectral_bins, nr_extra, 1,alignment,fits_map_IDs,heal_type);
    
    detector->setObsPosition(Vector3D(sx, sy, sz), Vector3D(0, 0, 0), tmp_l_min, tmp_l_max, tmp_b_min, tmp_b_max);
    
    if(project_healpix)
    {
        if(detector_proj != 0)
        {
            delete detector_proj;
            detector_proj = 0;
        }

        uint special_param=0;
        distance = dust_ray_detectors[pos + 6];

        sidelength_x = abs(dust_ray_detectors[pos + 7]-dust_ray_detectors[pos + 8]);
        sidelength_y = abs(dust_ray_detectors[pos + 9]-dust_ray_detectors[pos + 10]);

        max_length = _max_length;

        map_shift_x = 0;
        map_shift_y = 0;

        proj_x = param.getProjectedFitsX();
        proj_y = param.getProjectedFitsY();

        max_subpixel_lvl = 1;
        
        detector_proj = new CDetector(DET_PLANE,
                                    path,
                                    proj_x,
                                    proj_y,
                                    sID,
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
    }

    return true;
}

void CRaytracingHealPix::initIndices()
{
    // Init counter and percentage to show progress
    ullong per_counter = 0;
    float last_percentage = 0;
    bool do_print = false;
    
    // just in case...
    heal_indices.clear();
    arr_theta.clear();
    arr_phi.clear();
    
    if(nside>1024)
    {
        cout << CLR_LINE;
        cout << "-> Selecting heal. indices: 0 [%]      \r" << flush;
        do_print = true;
    }
    
    for(int64_t i_pix = 0; i_pix < npix; i_pix++)
    {   
        if(do_print)
        {
            per_counter++;

            // Calculate percentage of total progress per source
            float percentage = 100.0 * float(per_counter) / float(npix);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
                cout << "-> Selecting heal. indices: " << percentage << " [%]       \r" << flush;

                last_percentage = percentage;
            }        
        }
        
        double theta=0;
        double phi=0;
        
        if(getRefPosition(i_pix,theta,phi))
        {
            heal_indices.push_back(i_pix);
            arr_theta.push_back(theta);
            arr_phi.push_back(phi);
        }        
    }
    
    bool sorted= is_sorted(heal_indices.begin(), heal_indices.end());
    
    if(!sorted)
    {
        cout << CLR_LINE;
        cout << WARNING_LINE << "Healpix indices are not sorted!\n";    
    }
    
    cout << CLR_LINE;
    return;
}

bool CRaytracingHealPix::setSyncDetector(uint pos,
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

    double tmp_l_min = sync_ray_detectors[pos + 7];
    double tmp_l_max = sync_ray_detectors[pos + 8];

    double tmp_b_min = sync_ray_detectors[pos + 9];
    double tmp_b_max = sync_ray_detectors[pos + 10];

    //l_min = PI * (-sync_ray_detectors[pos + 8] + 180.0) / 180.0;
    //l_max = PI * (-sync_ray_detectors[pos + 7] + 180.0) / 180.0;
    
    l_min = PI * (-sync_ray_detectors[pos + 7]) / 180.0;
    l_max = PI * (-sync_ray_detectors[pos + 8] + 360.0) / 180.0;
    
    b_min = PI * (-sync_ray_detectors[pos + 10] + 90.0) / 180;
    b_max = PI * (-sync_ray_detectors[pos + 9] + 90.0) / 180;

    if(sync_ray_detectors[pos + 11]>0)
        rad_bubble = sync_ray_detectors[pos + 11];

    nside = uint(sync_ray_detectors[pos + NR_OF_RAY_DET - 1]);

    npix = 12 * nside * nside;

    det_pos.setX(sx);
    det_pos.setY(sy);
    det_pos.setZ(sz);

    max_length = _max_length * 10;

    setOrientation(param.getHealpixOrientation());
    


    detector = new CDetector(path, npix, 1, det_pos, max_length, lam_min, lam_max,rad_bubble, nr_spectral_bins, nr_extra, 1, ALIG_PA,0, false);
    detector->setObsPosition(Vector3D(sx, sy, sz), Vector3D(0, 0, 0), tmp_l_min, tmp_l_max, tmp_b_min, tmp_b_max);

    return true;
}

bool CRaytracingHealPix::setDustAMEDetector(uint pos,
                             const parameters & param,
                             dlist ame_ray_detectors,
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

    double lam_min = ame_ray_detectors[pos + 0];
    double lam_max = ame_ray_detectors[pos + 1];
    nr_spectral_bins = uint(ame_ray_detectors[pos + 2]);
    nr_extra = 1;

    sID = uint(ame_ray_detectors[pos + 3]);

    sx = ame_ray_detectors[pos + 4];
    sy = ame_ray_detectors[pos + 5];
    sz = ame_ray_detectors[pos + 6];

    double tmp_l_min = ame_ray_detectors[pos + 7];
    double tmp_l_max = ame_ray_detectors[pos + 8];

    double tmp_b_min = ame_ray_detectors[pos + 9];
    double tmp_b_max = ame_ray_detectors[pos + 10];

    //l_min = PI * (-free_ray_detectors[pos + 8] + 180.0) / 180.0;
    //l_max = PI * (-free_ray_detectors[pos + 7] + 180.0) / 180.0;
    
    l_min = PI * (-ame_ray_detectors[pos + 7]) / 180.0;
    l_max = PI * (-ame_ray_detectors[pos + 8] + 360.0) / 180.0;
    
    
    b_min = PI * (-ame_ray_detectors[pos + 10] + 90.0) / 180;
    b_max = PI * (-ame_ray_detectors[pos + 9] + 90.0) / 180;

    if(ame_ray_detectors[pos + 11]>0)
        rad_bubble = ame_ray_detectors[pos + 11];

    nside = uint(ame_ray_detectors[pos + NR_OF_RAY_DET - 1]);

    npix = 12 * nside * nside;

    det_pos.setX(sx);
    det_pos.setY(sy);
    det_pos.setZ(sz);

    max_length = _max_length * 10;

    setOrientation(param.getHealpixOrientation());

    detector = new CDetector(path, npix, 1, det_pos, max_length, lam_min, lam_max,rad_bubble, nr_spectral_bins, nr_extra, 3,ALIG_PA,0,false);
    detector->setObsPosition(Vector3D(sx, sy, sz), Vector3D(0, 0, 0), tmp_l_min, tmp_l_max, tmp_b_min, tmp_b_max);

    return true;
}

bool CRaytracingHealPix::setFreeFreeDetector(uint pos,
                                         const parameters & param,
                                         dlist free_ray_detectors,
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

    double lam_min = free_ray_detectors[pos + 0];
    double lam_max = free_ray_detectors[pos + 1];
    nr_spectral_bins = uint(free_ray_detectors[pos + 2]);
    nr_extra = 1;

    sID = uint(free_ray_detectors[pos + 3]);

    sx = free_ray_detectors[pos + 4];
    sy = free_ray_detectors[pos + 5];
    sz = free_ray_detectors[pos + 6];

    double tmp_l_min = free_ray_detectors[pos + 7];
    double tmp_l_max = free_ray_detectors[pos + 8];

    double tmp_b_min = free_ray_detectors[pos + 9];
    double tmp_b_max = free_ray_detectors[pos + 10];

    //l_min = PI * (-free_ray_detectors[pos + 8] + 180.0) / 180.0;
    //l_max = PI * (-free_ray_detectors[pos + 7] + 180.0) / 180.0;
    
    l_min = PI * (-free_ray_detectors[pos + 7]) / 180.0;
    l_max = PI * (-free_ray_detectors[pos + 8] + 360.0) / 180.0;
    
    
    b_min = PI * (-free_ray_detectors[pos + 10] + 90.0) / 180;
    b_max = PI * (-free_ray_detectors[pos + 9] + 90.0) / 180;

    if(free_ray_detectors[pos + 11]>0)
        rad_bubble = free_ray_detectors[pos + 11];

    nside = uint(free_ray_detectors[pos + NR_OF_RAY_DET - 1]);

    npix = 12 * nside * nside;

    det_pos.setX(sx);
    det_pos.setY(sy);
    det_pos.setZ(sz);

    max_length = _max_length * 10;

    setOrientation(param.getHealpixOrientation());
    


    detector = new CDetector(path, npix, 1, det_pos, max_length, lam_min, lam_max,rad_bubble, nr_spectral_bins, nr_extra, 3,ALIG_PA,0,false);
    detector->setObsPosition(Vector3D(sx, sy, sz), Vector3D(0, 0, 0), tmp_l_min, tmp_l_max, tmp_b_min, tmp_b_max);

    return true;
}

bool CRaytracingHealPix::setLineDetector(uint pos,
                                         const parameters & param,
                                         dlist line_ray_detectors,
                                         string path,
                                         double _max_length, bool hasZeeman)
{
    rt_detector_shape = DET_SPHER;
    vel_maps_fits = param.getVelMapsFits();
    heal_type = param.getHealType();

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

    sx = line_ray_detectors[pos + 4];
    sy = line_ray_detectors[pos + 5];
    sz = line_ray_detectors[pos + 6];

    //l_min = PI * (-line_ray_detectors[pos + 8] + 180.0) / 180.0;
    //l_max = PI * (-line_ray_detectors[pos + 7] + 180.0) / 180.0;
    
    l_min = PI * (-line_ray_detectors[pos + 7]) / 180.0;
    l_max = PI * (-line_ray_detectors[pos + 8] + 360.0) / 180.0;
    
    b_min = PI * (-line_ray_detectors[pos + 10] + 90.0) / 180;
    b_max = PI * (-line_ray_detectors[pos + 9] + 90.0) / 180;

    setOrientation(param.getHealpixOrientation());

    vx = line_ray_detectors[pos + 11];
    vy = line_ray_detectors[pos + 12];
    vz = line_ray_detectors[pos + 13];

    nside = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 2]);
    nr_spectral_bins = uint(line_ray_detectors[pos + NR_OF_LINE_DET - 1]);
    nr_extra = 1;

    npix = 12 * nside * nside;

    max_length = _max_length * 10;

    det_pos.setX(sx);
    det_pos.setY(sy);
    det_pos.setZ(sz);
    detector = new CDetector(path, npix, 1, det_pos, max_length, i_trans, nr_spectral_bins, min_velocity, max_velocity, hasZeeman);
    detector->setObsPosition(Vector3D(sx, sy, sz), Vector3D(vx, vy, vz), l_min, l_max, b_min, b_max);

    return true;
}

void CRaytracingHealPix::setOrientation(uint orientation_reference)
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

long CRaytracingHealPix::getNpix()
{
    if(!heal_indices.empty())
        return heal_indices.size();
        
    return npix;
}

double CRaytracingHealPix::getMinArea()
{
    return 4 * PI / double(npix);
}

Vector3D CRaytracingHealPix::getObserverVelocity()
{
    return Vector3D(vx, vy, vz);
}

bool CRaytracingHealPix::isNotAtCenter(photon_package * pp, double cx, double cy)
{
    /*double theta = cx;
    double phi = cy;

    if(theta <(45-40)*PI/180)
        return false;

    if(theta > (45+40)*PI/180)
        return false;

    if(phi <-40*PI/180)
        return false;

    if(phi > 40*PI/180)
        return false;*/

    Vector3D ph_dir = pp->getDirection();
    Vector3D ph_pos = pp->getPosition();
    Vector3D new_dir = det_pos - (rad_bubble*pp->getDirection()+ph_pos);

    double lam = ph_dir * new_dir;

    if(lam <= 0)
        return false;

    return true;
}

void CRaytracingHealPix::preparePhoton(photon_package * pp, double cx, double cy)
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

    start_pos += 1e3*max_length * tmp_ez + det_pos;

    pp->setPosition(start_pos);
    pp->setEX(tmp_ex);
    pp->setEY(-tmp_ey);
    pp->setEZ(-tmp_ez);
}

void CRaytracingHealPix::preparePhotonWithPosition(photon_package * pp, Vector3D pos, int64_t & i_pix)
{
    pp->setPosition(pos);

    Vector3D tmp_ex, tmp_ey;
    Vector3D tmp_ez = pos - det_pos;
    tmp_ez.normalize();

    double theta = acos(tmp_ez.Z()) - (PI2 - detector_angle_offset.Theta());
    double phi = Vector3D::atan3(tmp_ez.X(), tmp_ez.Y()) - (PI + detector_angle_offset.Phi());

    tmp_ey.setX(cos(theta) * cos(phi));
    tmp_ey.setY(cos(theta) * sin(phi));
    tmp_ey.setZ(-sin(theta));

    tmp_ex.setX(-sin(phi));
    tmp_ex.setY(cos(phi));
    tmp_ex.setZ(0);

    ang2pix_ring64(theta, phi, &i_pix);

    pp->setEX(tmp_ex);
    pp->setEY(-tmp_ey);
    pp->setEZ(-tmp_ez);
}

void CRaytracingHealPix::setDirection(photon_package * pp)
{
    Vector3D dir;
    Vector3D pos = pp->getPosition();

    dir = (pos - det_pos) / max_length;
    dir.normalize();

    pp->setEZ(dir);
}

void CRaytracingHealPix::setPosition(Vector3D pos)
{
    sx = pos.X();
    sy = pos.Y();
    sz = pos.Z();

    det_pos.setX(sx);
    det_pos.setY(sy);
    det_pos.setZ(sz);
}

//todo: check
bool CRaytracingHealPix::getRelPosition(int64_t i_pix, double & cx, double & cy)
{
    bool result = false;
    
    if(is_patch)
    {
        int64_t tmp_pix=heal_indices[i_pix];
        result = getRefPosition(tmp_pix, cx, cy);
    }
    else
        result = getRefPosition(i_pix, cx, cy);

    return result;
}

bool CRaytracingHealPix::getRefPosition(int64_t i_pix, double & cx, double & cy)
{
    pix2ang_ring64(i_pix, &cx, &cy);

    if(cx < b_min)
        return false;

    if(cx > b_max)
        return false;
    
     if(cy > l_min && cy < l_max)
        return false;

    /*if(cy < l_min)
        return false;

    if(cy > l_max)
        return false;*/

    return true;
}

double CRaytracingHealPix::getDistance()
{
    return 1.0;
}

double CRaytracingHealPix::getDistance(Vector3D pos)
{
    Vector3D pos_obs = Vector3D(sx, sy, sz);
    Vector3D diff = (pos - pos_obs);
    return diff.length();
}

void CRaytracingHealPix::addToDetector(photon_package * pp, int64_t i_pix, bool direct)
{
    for(uint i_spectral = 0; i_spectral < nr_spectral_bins * nr_extra; i_spectral++)
    {
        // Set wavelength of photon package
        pp->setSpectralID(i_spectral);

        // Multiply by min area if such a multiplication did not happen before
        if(!direct)
            pp->getStokesVector(i_spectral)->multStokesParam(getMinArea());

        // Add photon Stokes vector to detector
        detector->addToRaytracingDetector(*pp, i_pix);
        detector->addToRaytracingSedDetector(*pp);
    }
}

bool CRaytracingHealPix::writeDustResults(uint ray_result_type)
{
    if(heal_type==0)
    {
        if(!detector->writeDustHealMaps(dID, heal_indices, npix, ray_result_type))
            return false;
    }
    else
    {
        if((heal_type & HEAL_NONE) == HEAL_NONE)
        {
            cout << CLR_LINE;
            cout << "Skipping writing of healpix maps!\n" << flush;
        }
        else
        {
            if((heal_type & HEAL_INDEX) == HEAL_INDEX)
            {
                if(!detector->writeDustHealMapsTiny(dID, heal_indices, npix, ray_result_type))
                    return false;
            }
            
            if((heal_type & HEAL_FULL) == HEAL_FULL)
            {
                if(!detector->writeDustHealMaps(dID, heal_indices, npix, ray_result_type))
                    return false;
            }
        }
    }
    
    if(project_healpix)
    {
        projectHealMaps();
        detector_proj->writeProjMap(dID, ray_result_type);
    }
    
    
    if(!detector->writeSed(dID, ray_result_type))
        return false;

    return true;
}

bool CRaytracingHealPix::writeLineResults(CGasMixture * gas, uint i_species, uint i_line)
{
    if(vel_maps_fits)
    {
        if(!detector->writeVelChannelHealMaps(gas, i_species, i_line))
            return false;
    }
    
    //if(compact_fits1)
    {
        if(!detector->writeVelChannelHealMapsTiny(gas, i_species, i_line))
            return false;
    }

    if(!detector->writeIntVelChannelHealMaps(gas, i_species, i_line))
        return false;

    if(!detector->writeLineSpectrum(gas, i_species, i_line))
        return false;

    return true;
}

bool CRaytracingHealPix::writeFreeFreeResults()
{
    if(!detector->writeFreeFreeHealMap(dID))
        return false;

    return true;
}

bool CRaytracingHealPix::writeDustAMEResults()
{
    if(!detector->writeDustAMEHealMap(dID))
        return false;

    return true;
}


bool CRaytracingHealPix::writeSyncResults()
{
    if(!detector->writeSyncHealMap(dID))
        return false;

    return true;
}

void CRaytracingHealPix::setObserverPosition(Vector3D pos)
{
    sx = pos.X();
    sy = pos.Y();
    sz = pos.Z();
}

/*void CRaytracingHealPix::pix2ang_ring(int i_pix, double * theta, double * phi)
{
    double z;
    pix2ang_ring_z_phi(nside, i_pix, &z, phi);
    *theta = acos(z);
}

void CRaytracingHealPix::ang2ring_ring(double theta, double phi, int * i_pix)
{
    double z = cos(theta);
    ang2pix_ring_z_phi(nside, z, phi, i_pix);
}


void CRaytracingHealPix::pix2ang_ring_z_phi(int nside_, int pix, double * z, double * phi)
{
    long ncap_ = nside_ * (nside_ - 1) * 2;
    long npix_ = 12 * nside_ * nside_;
    double fact2_ = 4. / npix_;
    if(pix < ncap_) // North Polar cap 
    {
        int iring = (1 + isqrt(1 + 2 * pix)) >> 1; // counted from North pole
        int iphi = (pix + 1) - 2 * iring * (iring - 1);

        *z = 1.0 - (iring * iring) * fact2_;
        *phi = (iphi - 0.5) * PI2 / iring;
    }
    else if(pix < (npix_ - ncap_)) // Equatorial region
    {
        double fact1_ = (nside_ << 1) * fact2_;
        int ip = pix - ncap_;
        int iring = ip / (4 * nside_) + nside_; // counted from North pole 
        int iphi = ip % (4 * nside_) + 1;
        // 1 if iring+nside is odd, 1/2 otherwise 
        double fodd = ((iring + nside_) & 1) ? 1 : 0.5;

        int nl2 = 2 * nside_;
        *z = (nl2 - iring) * fact1_;
        *phi = (iphi - fodd) * PI / nl2;
    }
    else // South Polar cap 
    {
        int ip = npix_ - pix;
        int iring = (1 + isqrt(2 * ip - 1)) >> 1; // counted from South pole 
        int iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));

        *z = -1.0 + (iring * iring) * fact2_;
        *phi = (iphi - 0.5) * PI2 / iring;
    }
}

void CRaytracingHealPix::ang2pix_ring_z_phi(int nside_, double z, double phi, int * pix)
{
    double za = abs(z);
    double tt = CMathFunctions::fmodulo(phi, PIx2) * invPI2; // in [0,4)
    if(za <= TWOTHIRD)                                       // Equatorial region
    {
        double temp1 = nside_ * (0.5 + tt);
        double temp2 = nside_ * z * 0.75;
        int jp = (int)(temp1 - temp2); // index of  ascending edge line 
        int jm = (int)(temp1 + temp2); // index of descending edge line 

        /* ring number counted from z=2/3 
        int ir = nside_ + 1 + jp - jm; // in {1,2n+1} 
        int kshift = 1 - (ir & 1);     // kshift=1 if ir even, 0 otherwise 

        int ip = (jp + jm - nside_ + kshift + 1) / 2; /* in {0,4n-1} 
        ip = CMathFunctions::imodulo(ip, 4 * nside_);

        *pix = nside_ * (nside_ - 1) * 2 + (ir - 1) * 4 * nside_ + ip;
    }
    else // North & South polar caps
    {
        double tp = tt - (int)(tt);
        double tmp = nside_ * sqrt(3 * (1 - za));

        int jp = (int)(tp * tmp);         // increasing edge line index 
        int jm = (int)((1.0 - tp) * tmp); // decreasing edge line index 

        int ir = jp + jm + 1;    // ring number counted from the closest pole 
        int ip = (int)(tt * ir); // in {0,4*ir-1} 
        ip = CMathFunctions::imodulo(ip, 4 * ir);

        if(z > 0)
            *pix = 2 * ir * (ir - 1) + ip;
        else
            *pix = 12 * nside_ * nside_ - 2 * ir * (ir + 1) + ip;
    }
}*/

int CRaytracingHealPix::isqrt(int v)
{
    return (int)(sqrt(v + 0.5));
}


void CRaytracingHealPix::pix2ang_ring64(int64_t ipix, double *theta, double *phi)
{
    double z,s;
    pix2ang_ring_z_phi64 (ipix,&z,&s,phi);
    *theta= (s<-2) ? acos(z) : atan2(s,z);
}

void CRaytracingHealPix::ang2pix_ring64(double theta, double phi, int64_t *ipix)
{
    //UTIL_ASSERT((theta>=0)&&(theta<=pi),"theta out of range");
    double cth=cos(theta), sth=(fabs(cth)>0.99) ? sin(theta) : -5;
    *ipix=ang2pix_ring_z_phi64 (nside,cth,sth,phi);
}

long CRaytracingHealPix::isqrt64(int64_t v)
{
    int64_t res = sqrt(v+0.5);
    if (v<((int64_t)(1)<<50)) return (long)res;
    if (res*res>v)
        --res;
    else if ((res+1)*(res+1)<=v)
        ++res;
    
    return (long)res;
}


int CRaytracingHealPix::ring_neighbors_standalone64(int64_t ipix, int64_t out8[8])
{
    for(int i=0; i<8; i++)
        out8[i] = -1;

    int ix, iy, face;
    ring2xyf64(nside, ipix, &ix, &iy, &face);

    int nsm1 = (int)nside - 1;

    // get pixel center angles once (for fallback probes)
    double theta_c = 0.0, phi_c = 0.0;
    pix2ang_ring64(ipix, &theta_c, &phi_c);

    // small angular step roughly a pixel "radius"
    const double step = 1.05 / (std::sqrt(3.0) * (double)nside);
    const double inv_sqrt2 = 0.70710678118654752440;
    const double sth = std::sin(theta_c);
    const double dphi_unit = (std::fabs(sth) > 1e-12) ? (step / sth) : 0.0;

    // direction-aligned deltas matching NB_DX/NB_DY ordering:
    // SW,W,NW,N,NE,E,SE,S
    const double dth[8] = { +step*inv_sqrt2, 0.0, -step*inv_sqrt2, -step, -step*inv_sqrt2, 0.0, +step*inv_sqrt2, +step };
    const double dph[8] = { -dphi_unit*inv_sqrt2, -dphi_unit, -dphi_unit*inv_sqrt2, 0.0, +dphi_unit*inv_sqrt2, +dphi_unit, +dphi_unit*inv_sqrt2, 0.0 };

    for (int m=0; m<8; ++m)
    {
        int nx = ix + NB_DX[m];
        int ny = iy + NB_DY[m];

        if (nx >= 0 && nx < (int)nside && ny >= 0 && ny < (int)nside)
        {
              // same face: exact integer map, O(1)
              out8[m] = xyf2ring64(nside, nx, ny, face);
        }
        else
        {
              // crossed a face: robust tiny probe toward that neighbor
              double th = clamp_double(theta_c + dth[m], 1e-9, 3.14159265358979323846 - 1e-9);
              double ph = wrap2pi(phi_c   + dph[m]);
              int64_t nb = -1;
              ang2pix_ring64(th, ph, &nb);
              out8[m] = nb;
        }
    }
    
    return 8;
}

void CRaytracingHealPix::ang2vec(double theta, double phi, Vector3D & v)
{
    double st = sin(theta);
    v.setX(st * cos(phi));
    v.setY(st * sin(phi));
    v.setZ(cos(theta));
}

/* center-vector of a RING pixel using your (z,s,phi) routine */
void CRaytracingHealPix::pixcenter_vec_ring64(int64_t ipix, Vector3D & v)
{
    double z = 0.0;
    double s = 0.0;
    double phi = 0.0;

    pix2ang_ring_z_phi64(ipix, &z, &s, &phi);
    v.setX(s * std::cos(phi));
    v.setY(s * std::sin(phi));
    v.setZ(z);
}

/* great-circle distance between query (theta,phi) and pixel center (via vectors) */
double CRaytracingHealPix::ang_distance_to_pixel_center64(int64_t ipix, double theta, double phi)
{
    Vector3D q, p;
    ang2vec(theta, phi, q);
    
    double theta_c;
    double phi_c;
            
    pix2ang_ring64(ipix, &theta_c, &phi_c);
    ang2vec(theta_c, phi_c, p);
    
    double dotp = q.X()*p.X() + q.Y()*p.Y() + q.Z()*p.Z();
    
    
    if (dotp > 1.0) dotp = 1.0;
    if (dotp < -1.0) dotp = -1.0;
    return std::acos(dotp);
}

/* ---------------- Main interpolation (center + three nearest neighbors) ----------------
   hp_map:  input HEALPix RING map (size = 12*nside^2), double
   out:     output matrix, row-major, size Ny * Nx
   Nx, Ny:  grid resolution (x ~ phi index, y ~ theta index)
   theta0, theta1 in [0,pi], phi0, phi1 in [0,2pi)

   For each grid point (theta,phi):
     - find center pixel ip_center
     - collect up to 8 neighbor candidates around (theta,phi)
     - pick the three closest of those candidates
     - inverse-distance blend of 4 samples: center + 3 nearest neighbors*/

void CRaytracingHealPix::interpolate_ring_to_regular_grid_4pt(const double* hp_map,
                                          double* out, int Nx, int Ny,
                                          double theta0, double theta1,
                                          double phi0, double phi1)
{
    double dtheta, dphi;
    
    if(Ny > 1)
        dtheta = (theta1 - theta0) / (double)(Ny - 1);
    else
        dtheta = 0.0;

    if(Nx > 1)
        dphi   = (phi1   - phi0) / (double)(Nx - 1);
    else
        dphi   = 0.0;

    const double eps = 1e-12;

    for (int iy = 0; iy < Ny; ++iy)
    {
        double theta = theta0 + dtheta * (double)iy;
        
        if(theta < 0.0)
            theta = 0.0;
        
        if(theta > PI)
            theta = PI;

        for (int ix = 0; ix < Nx; ++ix)
        {
            double phi = phi0 + dphi * (double)ix;
            phi = wrap2pi(phi);

            /* 1) center pixel for this query direction */
            int64_t ip_center = -1;
            ang2pix_ring64(theta, phi, &ip_center);

            /* 2) collect neighbor candidates around the query direction */
            int64_t cand[8];
            int ncand = collect_neighbor_candidates64(theta, phi, cand, 8);

            /* remove any accidental center duplicates from candidates */
            int iwrite = 0;
            
            for (int k = 0; k < ncand; ++k)
            {
                if (cand[k] != ip_center)
                {
                    cand[iwrite] = cand[k];
                    iwrite += 1;
                }
            }
            ncand = iwrite;

            /* 3) choose the 3 closest neighbors */
            int64_t nn[3] = { -1, -1, -1 };
            int nnb = pick_three_closest64(theta, phi, cand, ncand, nn);

            /* 4) inverse-distance blend of center + up to 3 neighbors */
            /* If we land exactly on a pixel center, return that value directly. */
            double best_dist = ang_distance_to_pixel_center64(ip_center, theta, phi);
            
            if (best_dist < 1e-14)
            {
                out[iy * Nx + ix] = hp_map[ip_center];
                continue;
            }

            double num = 0.0;
            double den = 0.0;

            /* center weight */
            double w0 = 1.0 / (best_dist + eps);
            num += w0 * hp_map[ip_center];
            den += w0;

            /* neighbor weights */
            for (int t = 0; t < nnb; ++t)
            {
                if(nn[t] < 0)
                    continue;

                double d = ang_distance_to_pixel_center64(nn[t], theta, phi);
                double w = 1.0 / (d + eps);
                num += w * hp_map[nn[t]];
                den += w;
            }

            if (den > 0.0)
                out[iy * Nx + ix] = num / den;
            else
                out[iy * Nx + ix] = 0.0;
        }
    }
}

/* Select the three closest neighbors (by great-circle distance to the query point).
   Expects 'cands' to be unique and not including the center pixel. */
int CRaytracingHealPix::pick_three_closest64(double theta, double phi, const int64_t *cands,
                                int ncand, int64_t out3[3])
{
    double best_d[3];
    int64_t best_i[3];
    int nbest = 0;

    for (int i = 0; i < 3; ++i)
    {
        best_d[i] = 1e300;
        best_i[i] = -1;
    }

    for (int i = 0; i < ncand; ++i)
    {
        double d = ang_distance_to_pixel_center64(cands[i], theta, phi);

        /* insert in sorted triplet (smallest first) */
        if (d < best_d[0])
        {
              best_d[2] = best_d[1]; best_i[2] = best_i[1];
              best_d[1] = best_d[0]; best_i[1] = best_i[0];
              best_d[0] = d;         best_i[0] = cands[i];

              if (nbest < 3)
                  nbest += 1;
        }
        else if(d < best_d[1])
        {
              best_d[2] = best_d[1]; best_i[2] = best_i[1];
              best_d[1] = d;         best_i[1] = cands[i];
              if (nbest < 3)
                  nbest += 1;
        }
        else if(d < best_d[2])
        {
              best_d[2] = d;         best_i[2] = cands[i];
              if (nbest < 3)
                  nbest += 1;
        }
    }

    for (int k = 0; k < 3; ++k)
        out3[k] = best_i[k];
    
    return nbest;
}

/* Collect up to 8 neighbor candidates around (theta,phi) by tiny spherical offsets.
   We only need the 3 closest later, so this is perfect. */
int CRaytracingHealPix::collect_neighbor_candidates64(double theta, double phi,
                                         int64_t *cands, int maxcands)
{
    /* step ~ pixel "radius" */
    const double step = 1.05 / (std::sqrt(3.0) * (double)nside);
    const double inv_sqrt2 = 0.70710678118654752440;

    double sth = std::sin(theta);
    double dphi_unit;

    if (fabs(sth) > 1e-12)
        dphi_unit = step / sth;
    else
        dphi_unit = 0.0;

    /* Directions: W, NW, N, NE, E, SE, S, SW (8 of them) */
    double dth[8], dph[8];
    dth[0] =  0.0;                 dph[0] = -dphi_unit;            /* W  */
    dth[1] = -step*inv_sqrt2;      dph[1] = -dphi_unit*inv_sqrt2;  /* NW */
    dth[2] = -step;                dph[2] =  0.0;                  /* N  */
    dth[3] = -step*inv_sqrt2;      dph[3] =  dphi_unit*inv_sqrt2;  /* NE */
    dth[4] =  0.0;                 dph[4] =  dphi_unit;            /* E  */
    dth[5] =  step*inv_sqrt2;      dph[5] =  dphi_unit*inv_sqrt2;  /* SE */
    dth[6] =  step;                dph[6] =  0.0;                  /* S  */
    dth[7] =  step*inv_sqrt2;      dph[7] = -dphi_unit*inv_sqrt2;  /* SW */

    int nfound = 0;

    for(int k = 0; k < 8; k++)
    {
        double th = clamp_double(theta + dth[k], 1e-9, 3.14159265358979323846 - 1e-9);
        double ph = wrap2pi(phi + dph[k]);

        int64_t nb = -1;
        ang2pix_ring64(th, ph, &nb);

        int dup = 0;
        for (int j = 0; j < nfound; ++j)
        {
            if (cands[j] == nb)
            {
                dup = 1;
                break;
            }
        }

        if (dup == 0)
        {
            if(nfound < maxcands)
            {
                cands[nfound] = nb;
                nfound += 1;
            }
        }
    }
    return nfound;
}

double CRaytracingHealPix::wrap2pi(double x)
{
  double y = fmod(x, PIx2);
  
  if (y < 0.0)
      y += PIx2;
  
  return y;
}

double CRaytracingHealPix::clamp_double(double x, double a, double b)
{
  if (x < a) return a;
  if (x > b) return b;
  
  return x;
}

int64_t CRaytracingHealPix::getHealIndex(int64_t value)
{
    int64_t left  = 0;
    int64_t right = (int64_t)heal_indices.size() - 1;

    while (left <= right)
    {
        int64_t mid = left + (right - left) / 2;
        int64_t midval = heal_indices[mid];

        if (midval == value)
        {
            return mid;      // exact match
        }
        else if (midval < value)
        {
            left = mid + 1;  // search right half
        }
        else
        {
            right = mid - 1; // search left half
        }
    }
    return -1; // not found
}

bool CRaytracingHealPix::projectHealMaps()
{
    double d_phi = (l_min + (PIx2-l_max)) / (double)(proj_x - 1); 
    double d_theta = (b_max - b_min) / (double)(proj_y - 1);
    
                // Calculate total number of pixel
    long per_max = proj_y * proj_x;

    // Init counter and percentage to show progress
    ullong per_counter = 0;
    float last_percentage = 0;
    
    #pragma omp parallel for schedule(dynamic)
    for(uint i_p=0; i_p<proj_x; i_p++)
    {
        double phi = l_min - d_phi * (double)i_p;
        
        for(uint i_t=0; i_t<proj_y; i_t++)
        {
            
            #pragma omp atomic update
            per_counter++;

            // Calculate percentage of total progress per source
            float percentage = 100.0 * float(per_counter) / float(per_max);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
                #pragma omp critical
                {
                    cout << "-> Projecting map: " << percentage << " [%]       \r" << flush;
                    last_percentage = percentage;
                }
            }
            
            double theta = b_min + d_theta * (double)i_t;
            
            double theta_c;
            double phi_c;
            
            int64_t ip_center = -1;
            ang2pix_ring64(theta, phi, &ip_center);
            pix2ang_ring64(ip_center, &theta_c, &phi_c);

            //cout << ip_center << endl << flush;
            
            int64_t cand[8];
            int ncand = collect_neighbor_candidates64(theta, phi, cand, 8);
            
            /* remove any accidental center duplicates from candidates */
            int iwrite = 0;
            
            for (int k = 0; k < ncand; ++k)
            {
                if (cand[k] != ip_center)
                {
                    cand[iwrite] = cand[k];
                    iwrite += 1;
                }
            }
            ncand = iwrite;

            /* 3) choose the 3 closest neighbors */
            int64_t nn[3] = { -1, -1, -1 };
            int nnb = pick_three_closest64(theta, phi, cand, ncand, nn);

            /* 4) inverse-distance blend of center + up to 3 neighbors */
            /* If we land exactly on a pixel center, return that value directly. */
            double best_dist = ang_distance_to_pixel_center64(ip_center, theta, phi);
                        
            int64_t center_pos = getHealIndex(ip_center);
            
            StokesVector center_st;
            photon_package pp;
            
            pp.setSpectralID(0);
            pp.setPosition(Vector3D(double(i_p), double(i_t),0));
                    
            if(center_pos!=-1)
                center_st = detector->getStokesVector(center_pos,0);

            if (best_dist < 1e-14)
            {
                pp.setStokesVector(center_st, MAX_UINT);
                    
                detector_proj->addToRaytracingDetector(pp,MAX_UINT-1);
                continue;
            }

            StokesVector num;
            double den = 0;

            /* center weight */
            if(center_pos!=-1)
            {
                double w0 = 1.0 / (best_dist + EPS_DOUBLE);
                num += w0 * center_st;
                den += w0;
            }

            /* neighbor weights */
            for (int i_ne = 0; i_ne < nnb; ++i_ne)
            {
                int64_t n_ipix=nn[i_ne];
                
                if(n_ipix < 0)
                    continue;
                
                int64_t n_pos = getHealIndex(n_ipix);
                
                if(n_pos < 0)
                    continue;
                
                StokesVector n_st = detector->getStokesVector(n_pos,0);
                
                double d = ang_distance_to_pixel_center64(n_ipix, theta, phi);
                double w = 1.0 / (d + EPS_DOUBLE);
                num += w * n_st;
                den += w;
            }
            
            if (den > 0.0)
            {
                num /= den;
                num.multSp1(1/den);
                pp.setStokesVector(center_st, MAX_UINT);        
                detector_proj->addToRaytracingDetector(pp,MAX_UINT-1);
            }
        }
    }
    
    cout << CLR_LINE;
    return true;
}

double CRaytracingHealPix::fmodulo (double v1, double v2)
{
  if (v1>=0)
    return (v1<v2) ? v1 : fmod(v1,v2);

  double tmp=fmod(v1,v2)+v2;
  return (tmp==v2) ? 0. : tmp;
}

int64_t CRaytracingHealPix::imodulo64 (int64_t v1, int64_t v2)
{
    int64_t v=v1%v2;
    return (v>=0) ? v : v+v2;
}

int64_t CRaytracingHealPix::ang2pix_ring_z_phi64 (int64_t nside_, double z, double s,  double phi)
{    
    double za = fabs(z);
    double tt = fmodulo(phi,PIx2) * invPI2; /* in [0,4) */

    if (za<=TWOTHIRD) /* Equatorial region */
    {
        double temp1 = nside_*(0.5+tt);
        double temp2 = nside_*z*0.75;
        int64_t jp = (int64_t)(temp1-temp2); /* index of  ascending edge line */
        int64_t jm = (int64_t)(temp1+temp2); /* index of descending edge line */

        /* ring number counted from z=2/3 */
        int64_t ir = nside_ + 1 + jp - jm; /* in {1,2n+1} */
        int kshift = 1-(ir&1); /* kshift=1 if ir even, 0 otherwise */

        int64_t ip = (jp+jm-nside_+kshift+1)/2; /* in {0,4n-1} */
        ip = imodulo64(ip,4*nside_);

        return nside_*(nside_-1)*2 + (ir-1)*4*nside_ + ip;
    }
    else  /* North & South polar caps */
    {
        double tp = tt-(int)(tt);
        double tmp = (s>-2.) ? nside_*s/sqrt((1.+za)/3.) : nside_*sqrt(3*(1-za));

        int64_t jp = (int64_t)(tp*tmp); /* increasing edge line index */
        int64_t jm = (int64_t)((1.0-tp)*tmp); /* decreasing edge line index */

        int64_t ir = jp+jm+1; /* ring number counted from the closest pole */
        int64_t ip = (int64_t)(tt*ir); /* in {0,4*ir-1} */
        ip = imodulo64(ip,4*ir);

    if (z>0)
        return 2*ir*(ir-1) + ip;
    else
        return 12*nside_*nside_ - 2*ir*(ir+1) + ip;
    }
}

void CRaytracingHealPix::pix2ang_ring_z_phi64 (int64_t pix, double *z, double *s, double *phi)
{
    int64_t ncap_=nside*(nside-1)*2;
    int64_t npix_=12*nside*nside;
    double fact2_  = 4./npix_;
    *s=-5;
  
    if (pix<ncap_) /* North Polar cap */
    {
        int64_t iring = (1+isqrt64(1+2*pix))>>1; /* from N pole */
        int64_t iphi  = (pix+1) - 2*iring*(iring-1);
        double tmp=(iring*iring)*fact2_;

        *z = 1.0 - tmp;
        if (*z>0.99) *s=sqrt(tmp*(2.-tmp));
        *phi = (iphi-0.5) * PI2/iring;
    }
    else if (pix<(npix_-ncap_)) /* Equatorial region */
    {
        double fact1_  = (nside<<1)*fact2_;
        int64_t ip  = pix - ncap_;
        int64_t iring = ip/(4*nside) + nside; /* counted from North pole */
        int64_t iphi  = ip%(4*nside) + 1;
        /* 1 if iring+nside is odd, 1/2 otherwise */
        double fodd = ((iring+nside)&1) ? 1 : 0.5;

        int64_t nl2 = 2*nside;
        *z = (nl2-iring)*fact1_;
        *phi = (iphi-fodd) * PI/nl2;
    }
    else /* South Polar cap */
    {
        int64_t ip = npix_ - pix;
        int64_t iring = (1+isqrt64(2*ip-1))>>1; /* from S pole */
        int64_t iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

        double tmp=(iring*iring)*fact2_;
        *z = tmp - 1.0;
        if (*z<-0.99) *s=sqrt(tmp*(2.-tmp));
        *phi = (iphi-0.5) * PI2/iring;
    }
}

int64_t CRaytracingHealPix::xyf2ring64 (int64_t nside_, int ix, int iy, int face_num)
  {
  int64_t nl4 = 4*nside_;
  int64_t jr = (jrll[face_num]*nside_) - ix - iy  - 1, jp;

  int64_t nr, kshift, n_before;
  if (jr<nside_)
    {
    nr = jr;
    n_before = 2*nr*(nr-1);
    kshift = 0;
    }
  else if (jr > 3*nside_)
    {
    nr = nl4-jr;
    n_before = 12*nside_*nside_ - 2*(nr+1)*nr;
    kshift = 0;
    }
  else
    {
    int64_t ncap_=2*nside_*(nside_-1);
    nr = nside_;
    n_before = ncap_ + (jr-nside_)*nl4;
    kshift = (jr-nside_)&1;
    }

  jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
  if (jp>nl4)
    jp-=nl4;
  else
    if (jp<1) jp+=nl4;

  return n_before + jp - 1;
  }
  

int64_t CRaytracingHealPix::special_div64 (int64_t a, int64_t b)
{
	int64_t t=(a>=(b<<1));
	a-=t*(b<<1);
	return (t<<1)+(a>=b);
}

void CRaytracingHealPix::ring2xyf64 (int64_t nside_, int64_t pix, int *ix, int *iy, int *face_num)
  {
  int64_t iring, iphi, kshift, nr;
  int64_t ncap_=2*nside_*(nside_-1);
  int64_t npix_=12*nside_*nside_;
  int64_t nl2 = 2*nside_;

  if (pix<ncap_) /* North Polar cap */
    {
    iring = (1+isqrt64(1+2*pix))>>1; /* counted from North pole */
    iphi  = (pix+1) - 2*iring*(iring-1);
    kshift = 0;
    nr = iring;
    *face_num=special_div64(iphi-1,nr);
    }
  else if (pix<(npix_-ncap_)) /* Equatorial region */
    {
    int64_t ip = pix - ncap_;
    iring = (ip/(4*nside_)) + nside_; /* counted from North pole */
    iphi  = (ip%(4*nside_)) + 1;
    kshift = (iring+nside_)&1;
    nr = nside_;
    int64_t ire = iring-nside_+1;
    int64_t irm = nl2+2-ire;
    int64_t ifm = (iphi - ire/2 + nside_ -1) / nside_;
    int64_t ifp = (iphi - irm/2 + nside_ -1) / nside_;
    *face_num = (ifp==ifm) ? (ifp|4) : ((ifp<ifm) ? ifp : (ifm+8));
    }
  else /* South Polar cap */
    {
    int64_t ip = npix_ - pix;
    iring = (1+isqrt64(2*ip-1))>>1; /* counted from South pole */
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
    kshift = 0;
    nr = iring;
    iring = 2*nl2-iring;
    *face_num=8+special_div64(iphi-1,nr);
    }

  int64_t irt = iring - (jrll[*face_num]*nside_) + 1;
  int64_t ipt = 2*iphi- jpll[*face_num]*nr - kshift -1;
  if (ipt>=nl2) ipt-=8*nside_;

  *ix =  (ipt-irt) >>1;
  *iy =(-(ipt+irt))>>1;
  }



