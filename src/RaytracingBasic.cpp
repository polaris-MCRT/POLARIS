/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "RaytracingBasic.hpp"

bool CRaytracingBasic::setSyncDetector(uint pos,
                                       const parameters & param,
                                       dlist sync_ray_detectors,
                                       double _max_length,
                                       string path)
{
    return false;
}

bool CRaytracingBasic::setFreeFreeDetector(uint pos,
                                       const parameters & param,
                                       dlist free_ray_detectors,
                                       double _max_length,
                                       string path)
{
    return false;
}

bool CRaytracingBasic::setLineDetector(uint pos,
                                       const parameters & param,
                                       dlist line_ray_detectors,
                                       string path,
                                       double _max_length, bool hasZeeman)
{
    return false;
}

bool CRaytracingBasic::setOPIATEDetector(uint pos,
                                         const parameters & param,
                                         dlist op_ray_detectors,
                                         string path,
                                         double _max_length)
{
    return false;
}

double CRaytracingBasic::getWavelength(uint i_wave)
{
    return detector->getWavelength(i_wave);
}

long CRaytracingBasic::getNpix()
{
    return 0;
}

double CRaytracingBasic::getMinArea()
{
    return (sidelength_x * sidelength_y) / double(map_pixel_x * map_pixel_y);
}

uint CRaytracingBasic::getSourceIndex()
{
    return sID;
}

uint CRaytracingBasic::considerPointSources()
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

uint CRaytracingBasic::getNrSpectralBins()
{
    return nr_spectral_bins;
}

void CRaytracingBasic::preparePhoton(photon_package * pp, double cx, double cy)
{
    Vector3D pos = cx * ex + cy * ey - max_length * ez;
    pp->setPosition(pos);
    pp->setEX(ex);
    pp->setEY(ey);
    pp->setEZ(ez);
}

void CRaytracingBasic::preparePhotonWithPosition(photon_package * pp, Vector3D pos, int & i_pix)
{
    pp->setPosition(pos);
    pp->setEX(ex);
    pp->setEY(ey);
    pp->setEZ(ez);

    i_pix = 0;
}

void CRaytracingBasic::setCoordinateSystem(photon_package * pp)
{
    pp->setEX(ex);
    pp->setEY(ey);
    pp->setEZ(ez);
}

bool CRaytracingBasic::getRelPosition(int i_pix, double & cx, double & cy)
{
    return getRelPositionMap(i_pix, cx, cy);
}

void CRaytracingBasic::calcMapParameter()
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

bool CRaytracingBasic::getRelPositionMap(int i_pix, double & cx, double & cy)
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

void CRaytracingBasic::setDetCoordSystem(const Vector3D & n1, const Vector3D & n2)
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

bool CRaytracingBasic::splitDustEmission()
{
    return split_emission;
}

uint CRaytracingBasic::getNrExtra()
{
    return nr_extra;
}

uint CRaytracingBasic::getDetectorShape()
{
    return rt_detector_shape;
}

bool CRaytracingBasic::getSubpixelWarning()
{
    return subpixel_warning;
}

double CRaytracingBasic::getDistance()
{
    return distance;
}

double CRaytracingBasic::getDistance(Vector3D pos)
{
    double proj_length = -ez * pos;
    return proj_length + getDistance();
}

void CRaytracingBasic::addToDetector(photon_package * pp, int i_pix, bool direct)
{}

void CRaytracingBasic::addToDetector(photon_package * pp1, photon_package * pp2, int i_pix, bool direct)
{
    // pos was only traced of first photon package
    pp2->setPosition(pp1->getPosition());

    // Add first package of photons to detector
    addToDetector(pp1, i_pix, direct);

    // Add second package of photons to detector
    addToDetector(pp2, i_pix, direct);
}

void CRaytracingBasic::setObserverPosition(Vector3D pos)
{}

bool CRaytracingBasic::postProcessing()
{
    return true;
}

bool CRaytracingBasic::writeDustResults(uint ray_result_type)
{
    if(!detector->writeMap(dID, ray_result_type))
        return false;

    if(!detector->writeSed(dID, ray_result_type))
        return false;

    return true;
}

bool CRaytracingBasic::writeLineResults(CGasMixture * gas, uint i_species, uint i_line)
{
    if(vel_maps_fits)
    {
        if(!detector->writeVelChannelMaps(gas, i_species, i_line))
            return false;
    }
    
    if(compact_fits)
    {
        //todo: include tiny fits
        if(!detector->writeVelChannelMaps(gas, i_species, i_line))
            return false;
    }

    if(!detector->writeIntChannelMaps(gas, i_species, i_line))
        return false;

    if(!detector->writeLineSpectrum(gas, i_species, i_line))
        return false;

    return true;
}

bool CRaytracingBasic::writeOpiateResults(COpiateDataBase * op)
{
    /*if(vel_maps)
        if(!detector->writeOPIATEVelChannelMaps(op,dID))
            return false;

    if(!detector->writeOPIATEIntChannelMaps(op,dID))
        return false;

    if(!detector->writeOPIATESpectrum(op,dID))
        return false;*/

    return true;
}

bool CRaytracingBasic::writeSyncResults()
{
    if(!detector->writeSyncMap(dID))
        return false;

    return true;
}

bool CRaytracingBasic::writeFreeFreeResults()
{
    if(!detector->writeFreeFreeMap(dID))
        return false;

    return true;
}

bool CRaytracingBasic::getUseSubpixel(double cx, double cy, uint subpixel_lvl)
{
    return false;
}

double CRaytracingBasic::getDistanceFactor()
{
    return 1.0 / (getDistance() * getDistance());
}

double CRaytracingBasic::getDistanceFactor(Vector3D pos)
{
    return 1.0 / (getDistance(pos) * (getDistance(pos)));
}

void CRaytracingBasic::getDetectorData(dlist & C, dlist & T)
{}

void CRaytracingBasic::setPosition(Vector3D pos){};

bool CRaytracingBasic::isNotAtCenter(photon_package * pp, double cx, double cy)
{
    return true;
}

Vector3D CRaytracingBasic::getObserverVelocity()
{
    return Vector3D(0, 0, 0);
}

void CRaytracingBasic::getSubPixelCoordinates(uint subpixel_lvl,
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

double CRaytracingBasic::getLamMin()
{
    return detector->getLamMin();
}

double CRaytracingBasic::getLamMax()
{
    return detector->getLamMax();
}

double CRaytracingBasic::getChannelWidth()
{
    return detector->getChannelWidth();
}

uint CRaytracingBasic::getNrOfSpectralBins()
{
    return detector->getNrOfSpectralBins();
}

double CRaytracingBasic::getVelocityChannel(uint vch)
{
    return detector->getVelocityChannel(vch);
}
