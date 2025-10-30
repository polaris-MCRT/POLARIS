/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "RaytracingPolar.hpp"

bool CRaytracingPolar::setDustDetector(uint pos,
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

bool CRaytracingPolar::setSyncDetector(uint pos,
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

    /*detector = new CDetector(rt_detector_shape,
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
    detector->setOrientation(n1, n2, rot_angle1, rot_angle2);*/

    return true;
}

bool CRaytracingPolar::setLineDetector(uint pos,
                                       const parameters & param,
                                       dlist line_ray_detectors,
                                       string path,
                                       double _max_length, bool hasZeeman)
{
    rt_detector_shape = DET_POLAR;
    vel_maps_fits = param.getVelMapsFits();

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
                                min_velocity,
                                max_velocity,
                                hasZeeman);
    detector->setOrientation(n1, n2, rot_angle1, rot_angle2);

    return true;
}

bool CRaytracingPolar::initPolarGridParameter()
{
    double pixel_width = max(sidelength_x / map_pixel_x, sidelength_y / map_pixel_y);
    double max_len = 0.5 * sqrt(sidelength_x * sidelength_x + sidelength_y * sidelength_y);
    if(!grid->getPolarRTGridParameter(max_len, pixel_width, max_subpixel_lvl, listR, npix_r, npix_ph))
    {
        cout << ERROR_LINE << "Polar detector can only be used with spherical or "
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

    npix_total++; // last pixel is central grid cell

    if(npix_total > MAX_RT_RAYS)
    {
        cout << INFO_LINE << "Very high amount of rays required for DUST EMISSION simulation with polar raytracing grid!" << endl;
        cout << "  Problem: The simulation may take a long time to process all rays." << endl;
        cout << "  Solution: Decrease max subpixel level or use the cartesian raytracing grid." << endl;
    }

    return true;
}

void CRaytracingPolar::initTmpStokes()
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

bool CRaytracingPolar::getRelPosition(int i_pix, double & cx, double & cy)
{
    // Return (0,0) coordinate, if central pixel
    if(i_pix == int(npix_total - 1))
    {
        cx = 0;
        cy = 0;
        return true;
    }

    uint rID, phID;
    getCoordinateIDs(uint(i_pix), rID, phID);

    double phi = (phID + 0.5) * PIx2 / double(npix_ph[rID]);

    double r1 = listR[rID];
    double r2 = listR[rID + 1];

    double radius = (r1 + r2) / 2.0;

    cx = radius * cos(phi);
    cy = radius * sin(phi);

    return true;
}

void CRaytracingPolar::addToDetector(photon_package * pp, int i_pix, bool direct)
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
        uint rID, phID;
        getCoordinateIDs(uint(i_pix), rID, phID);
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins * nr_extra; i_spectral++)
        {
            // Set wavelength of photon package
            pp->setSpectralID(i_spectral);

            // Set Stokes vector of polar ring element
            tmpStokes[pp->getSpectralID()][rID][phID] = *pp->getStokesVector();

            // Update Stokes vector of photon package
            pp->getStokesVector()->multStokesParam(getRingElementArea(rID));

            // Add photon Stokes vector SED detector
            detector->addToRaytracingSedDetector(*pp);
        }
    }
}

bool CRaytracingPolar::postProcessing()
{
    if(map_pixel_x * map_pixel_y > npix_total)
    {
        detector->setProcessingMethod(INTERP);
        cout << INFO_LINE << "Using 'interpolation' method to post process from polar to cartesian detector" << endl;
        return postProcessingUsingInterpolation();
    }
    else
    {
        detector->setProcessingMethod(NEAREST);
        cout << INFO_LINE << "Using 'nearest' method to post process from polar to cartesian detector" << endl;
        return postProcessingUsingNearest();
    }

    return false;
}

bool CRaytracingPolar::postProcessingUsingNearest()
{
    // Init counter and percentage to show progress
    ullong per_counter = 0;
    float last_percentage = 0;

    #pragma omp parallel for schedule(dynamic)
    for(int i_pix = 0; i_pix < npix_total; i_pix++)
    {
        photon_package pp = photon_package(nr_spectral_bins * nr_extra);

        // Init variables
        Vector3D pos;
        uint rID, phID;
        double cx = 0, cy = 0;

        // Increase counter used to show progress
        #pragma omp atomic update
        per_counter++;

        // Calculate percentage of total progress per source
        float percentage =
            100.0 * float(per_counter) / float(npix_total * nr_spectral_bins);

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

        // Get coordinates of the current pixel
        if(!getRelPosition(i_pix, cx, cy))
            continue;

        // Set photon package position and get R and Phi
        pos = Vector3D(cx, cy, 0);
        pp.setPosition(pos);
        pos.cart2cyl();

        getCoordinateIDs(uint(i_pix), rID, phID);

        for(uint i_spectral = 0; i_spectral < nr_spectral_bins * nr_extra; i_spectral++)
        {
            // Set wavelength of photon package
            pp.setSpectralID(i_spectral);

            pp.setStokesVector(tmpStokes[i_spectral][rID][phID], i_spectral);
            pp.getStokesVector(i_spectral)->multStokesParam(getRingElementArea(rID));

            // Transport pixel value to detector
            detector->addToRaytracingDetector(pp);
        }
    }

    cout << "-> Interpolating from polar grid to detector map: 100 [%]         \r" << flush;

    return true;
}

bool CRaytracingPolar::postProcessingUsingInterpolation()
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
        uint rID, rID1, rID2;
        int phID1, phID2, phID3, phID4;
        double r1, r2, p1, p2, p3, p4, dph1, dph2;
        double cx = 0, cy = 0;

        // Increase counter used to show progress
        #pragma omp atomic update
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
                pp.setStokesVector(tmpStokes[i_spectral][npix_r][0], i_spectral);
                pp.getStokesVector(i_spectral)->multStokesParam(getMinArea());

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
        if(phID2 >= int(npix_ph[rID1]))
            phID2 -= npix_ph[rID1];
        if(phID4 >= int(npix_ph[rID2]))
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
            double stokes_Sp = CMathFunctions::getPolarInterp(t, s, Q1.Sp1(), Q2.Sp1(), Q3.Sp1(), Q4.Sp1());

            // Add pixel value to photon_package
            StokesVector res_stokes =
                StokesVector(stokes_I, stokes_Q, stokes_U, stokes_V, stokes_T, stokes_Sp);
            pp.setStokesVector(res_stokes, i_spectral);
            pp.getStokesVector(i_spectral)->multStokesParam(getMinArea());
            // Transport pixel value to detector
            detector->addToRaytracingDetector(pp);
        }
    }

    cout << "-> Interpolating from polar grid to detector map: 100 [%]         \r" << flush;

    return true;
}

long CRaytracingPolar::getNpix()
{
    return npix_total;
}

void CRaytracingPolar::getCoordinateIDs(uint i_pix, uint & rID, uint & phID)
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

double CRaytracingPolar::getRingElementArea(uint rID)
{
    if(rID == npix_r)
        return PI * listR[0] * listR[0];

    double r1 = listR[rID];
    double r2 = listR[rID + 1];

    return PI * (r2 * r2 - r1 * r1) / npix_ph[rID];
}
