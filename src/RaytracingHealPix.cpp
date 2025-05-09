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

    l_min = PI * (-dust_ray_detectors[pos + 8] + 180.0) / 180.0;
    l_max = PI * (-dust_ray_detectors[pos + 7] + 180.0) / 180.0;
    b_min = PI * (-dust_ray_detectors[pos + 10] + 90.0) / 180;
    b_max = PI * (-dust_ray_detectors[pos + 9] + 90.0) / 180;

    if(dust_ray_detectors[pos + 11]>0)
        rad_bubble = dust_ray_detectors[pos + 11];

    nside = uint(dust_ray_detectors[pos + NR_OF_RAY_DET - 1]);

    npix = 12 * nside * nside;

    det_pos.setX(sx);
    det_pos.setY(sy);
    det_pos.setZ(sz);

    max_length = _max_length * 10;

    setOrientation(param.getHealpixOrientation());

    detector = new CDetector(
        path, npix, 1, det_pos, max_length, lam_min, lam_max, rad_bubble, nr_spectral_bins, 1, param.getAlignmentMechanism());
    detector->setObsPosition(Vector3D(sx, sy, sz), Vector3D(0, 0, 0), tmp_l_min, tmp_l_max, tmp_b_min, tmp_b_max);


    return true;
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

    l_min = PI * (-sync_ray_detectors[pos + 8] + 180.0) / 180.0;
    l_max = PI * (-sync_ray_detectors[pos + 7] + 180.0) / 180.0;
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

    detector =
        new CDetector(path, npix, 1, det_pos, max_length, lam_min, lam_max,rad_bubble, nr_spectral_bins, nr_extra);
    detector->setObsPosition(Vector3D(sx, sy, sz), Vector3D(0, 0, 0), tmp_l_min, tmp_l_max, tmp_b_min, tmp_b_max);

    return true;
}

bool CRaytracingHealPix::setLineDetector(uint pos,
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
    
    double min_velocity = line_ray_detectors[pos + 2];
    double max_velocity = line_ray_detectors[pos + 3];

    sx = line_ray_detectors[pos + 4];
    sy = line_ray_detectors[pos + 5];
    sz = line_ray_detectors[pos + 6];

    l_min = PI * (-line_ray_detectors[pos + 8] + 180.0) / 180.0;
    l_max = PI * (-line_ray_detectors[pos + 7] + 180.0) / 180.0;
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
    detector = new CDetector(path, npix, 1, det_pos, max_length, i_trans, nr_spectral_bins, min_velocity, max_velocity);
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

    start_pos += max_length * tmp_ez + det_pos;

    pp->setPosition(start_pos);
    pp->setEX(tmp_ex);
    pp->setEY(-tmp_ey);
    pp->setEZ(-tmp_ez);
}

void CRaytracingHealPix::preparePhotonWithPosition(photon_package * pp, Vector3D pos, int & i_pix)
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

    ang2ring_ring(theta, phi, &i_pix);

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

bool CRaytracingHealPix::getRelPosition(int i_pix, double & cx, double & cy)
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

void CRaytracingHealPix::addToDetector(photon_package * pp, int i_pix, bool direct)
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

bool CRaytracingHealPix::writeDustResults(uint ray_result_type)
{
    if(!detector->writeHealMaps(dID, ray_result_type))
        return false;

    if(!detector->writeSed(dID, ray_result_type))
        return false;

    return true;
}

bool CRaytracingHealPix::writeLineResults(CGasMixture * gas, uint i_species, uint i_line)
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

void CRaytracingHealPix::pix2ang_ring(int i_pix, double * theta, double * phi)
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

int CRaytracingHealPix::isqrt(int v)
{
    return (int)(sqrt(v + 0.5));
}

void CRaytracingHealPix::pix2ang_ring_z_phi(int nside_, int pix, double * z, double * phi)
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

void CRaytracingHealPix::ang2pix_ring_z_phi(int nside_, double z, double phi, int * pix)
{
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
