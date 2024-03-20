#include "Photon.h"
#include "Vector.h"
#include "MathFunctions.h"
#include "Surface.h"
#include "Parameters.h"
#include "Stokes.h"
#include "Grid.h"

void CSurface::initSurface(
    uint _surface_refl_model, uint _surface_pol_model,
    dlist _surface_refl_param, dlist _surface_pol_param)
{
    surface_refl_model = _surface_refl_model;
    surface_pol_model = _surface_pol_model;
    surface_refl_model_loc = surface_refl_model;
    surface_pol_model_loc = surface_pol_model;
    surface_refl_param = _surface_refl_param;
    surface_pol_param = _surface_pol_param;
    if(surface_refl_model == OCEAN) {
        initOceanParam();
        cout << "-> Calculating Integral for Cox-Munk wave facet distribution     \r" << flush;
        cm_norm = CMathFunctions::integCoxMunkNorm(sigma_sq);
    }
}

void CSurface::initOceanParam()
{
    // see Zhai et al. 2010, JQSRT 111, 1025
    // see Monahan & Ó Muircheartaigh 1980, JPO 10, 2094
    sigma_sq = 0.003 + 0.00512 * surface_refl_param[1];
    white_cap_fraction = 2.95e-6 * pow(surface_refl_param[1], 3.52);
}

void CSurface::setLocalParam(photon_package * pp, CRandomGenerator * rand_gen)
{
    surface_normal = pp->getPosition();
    surface_normal.normalize();

    if(surface_refl_model == OCEAN) {
        refractive_ratio = surface_refl_param[2];
    }

    surface_refl_model_loc = surface_refl_model;
    surface_pol_model_loc = surface_pol_model;
    if(surface_refl_model == OCEAN) {
        if(rand_gen->getRND() < white_cap_fraction) {
            surface_refl_model_loc = LAMBERTIAN;
            surface_pol_model_loc = DEPOL;
        }
    }

    cos_0 = -1*pp->getDirection() * surface_normal;
}

void CSurface::getEscapePhotonSurface(
    photon_package * pp,
    Vector3D obs_ex,
    Vector3D dir_obs,
    photon_package * pp_surface)
{
    // Init temporary photon package
    // Set the photon package at the position of the current photon
    pp_surface->setPosition(pp->getPosition());
    pp_surface->setPositionCell(pp->getPositionCell());

    // Synchronize the direction and wavelength as well
    pp_surface->setD(pp->getD());
    pp_surface->setWavelength(pp->getWavelength(), pp->getDustWavelengthID());

    switch(surface_refl_model_loc) {
        case LAMBERTIAN:
        {
            getEscapePhotonSurfaceLambertian(pp, dir_obs, pp_surface);
            break;
        }
        case LOMMELSEELIGER:
        {
            getEscapePhotonSurfaceLommelSeeliger(pp, dir_obs, pp_surface);
            break;
        }
        case OCEAN:
        {
            getEscapePhotonSurfaceOcean(pp, dir_obs, pp_surface);
            break;
        }
    }

    switch(surface_pol_model_loc) {
        case DEPOL:
        {
            dePolarization(pp_surface);
            break;
        }
        case SPECULAR:
        {
            specularPolarization(pp_surface);
            break;
        }
    }

    if(!pp_surface->getStokesVector()->isConsistent()) {
        pp_surface->getStokesVector()->depolarize();
    }
    
    // Rotate photon package into the coordinate space of the detector
    // (x or r)-axis of photon is y-axis of detector
    // (y or l)-axis of photon is negative x-axis of detector
    // see Figure 12 in O. Fischer (1993)
    double rot_angle_phot_obs = CMathFunctions::getRotationAngleObserver(obs_ex, pp_surface->getEX(), -1*pp_surface->getEY());
    pp_surface->getStokesVector()->rot(rot_angle_phot_obs);

    // The scattering part is based on O. Fischer (1993)
    // But on our detectors, U is defined the other way round
    pp_surface->getStokesVector()->multU(-1);
}

void CSurface::getEscapePhotonSurfaceLambertian(
    photon_package * pp,
    Vector3D dir_obs,
    photon_package * pp_surface)
{
    // Get the Stokes vector of the current photon package
    StokesVector tmp_stokes = *pp->getStokesVector();

    Vector3D d_in = pp->getDirection();

    // Calculate the theta and phi angle to the observer
    double phi_photon_to_obs = atan3(pp->getEY() * dir_obs, -1*pp->getEX() * dir_obs);
    double theta_photon_to_obs = acos(d_in * dir_obs);

    // Update the coordinate space of the photon
    pp_surface->updateCoordSystem(phi_photon_to_obs, theta_photon_to_obs);

    Vector3D w_normal = dir_obs - d_in;
    w_normal.normalize();
    cos_r = dir_obs * surface_normal;
    cos_i = -1*d_in * w_normal;

    // Calculate the fraction that is scattered into this direction
    double scat_fraction = cos_r / PI;

    // Reduce Stokes vector by albedo and scattering propability into theta
    tmp_stokes *= scat_fraction * surface_refl_param[0];

    // Rotate Stokes vector to new photon direction
    tmp_stokes.rot(phi_photon_to_obs);

    pp_surface->setStokesVector(tmp_stokes);
}

void CSurface::getEscapePhotonSurfaceLommelSeeliger(
    photon_package * pp,
    Vector3D dir_obs,
    photon_package * pp_surface)
{
    // Get the Stokes vector of the current photon package
    StokesVector tmp_stokes = *pp->getStokesVector();

    Vector3D d_in = pp->getDirection();

    // Calculate the theta and phi angle to the observer
    double phi_photon_to_obs = atan3(pp->getEY() * dir_obs, -1*pp->getEX() * dir_obs);
    double theta_photon_to_obs = acos(d_in * dir_obs);

    // Update the coordinate space of the photon
    pp_surface->updateCoordSystem(phi_photon_to_obs, theta_photon_to_obs);

    Vector3D w_normal = dir_obs - d_in;
    w_normal.normalize();
    cos_r = dir_obs * surface_normal;
    cos_i = -1*d_in * w_normal;

    // Calculate the fraction that is scattered into this direction
    double scat_fraction = cos_r / (cos_0 + cos_r);

    // Reduce Stokes vector by albedo and scattering propability into theta
    tmp_stokes *= scat_fraction * surface_refl_param[0] / PIx4;

    // Rotate Stokes vector to new photon direction
    tmp_stokes.rot(phi_photon_to_obs);

    pp_surface->setStokesVector(tmp_stokes);
}

void CSurface::getEscapePhotonSurfaceOcean(
    photon_package * pp,
    Vector3D dir_obs,
    photon_package * pp_surface)
{
    // Get the Stokes vector of the current photon package
    StokesVector tmp_stokes = *pp->getStokesVector();

    Vector3D d_in = pp->getDirection();

    // Calculate the theta and phi angle to the observer
    double phi_photon_to_obs = atan3(pp->getEY() * dir_obs, -1*pp->getEX() * dir_obs);
    double theta_photon_to_obs = acos(d_in * dir_obs);

    // Update the coordinate space of the photon
    pp_surface->updateCoordSystem(phi_photon_to_obs, theta_photon_to_obs);

    double cos_n;
    double cos_0_tmp = -1*d_in * surface_normal;
    cos_r = dir_obs * surface_normal;
    double scat_fraction = 0.0;
    Vector3D w_normal;

    w_normal = dir_obs - d_in;
    w_normal.normalize();

    cos_n = w_normal * surface_normal;
    cos_i = -1*d_in * w_normal;

    // get probability of water surface normal
    scat_fraction = CMathFunctions::getCoxMunkProb(cos_n, sigma_sq);
    scat_fraction /= (4.0 * cos_n * cos_0_tmp);

    // get shadowing function
    scat_fraction *= CMathFunctions::getShadowing(cos_0_tmp, cos_r, sigma_sq);

    // Reduce Stokes vector by albedo and scattering propability into theta and phi
    tmp_stokes *= CMathFunctions::effFresnelReflection(refractive_ratio, cos_i) * scat_fraction;
    tmp_stokes *= cos_i / (cm_norm.getValue(cos_0_tmp) * cos_n);

    // Rotate Stokes vector to new photon direction
    tmp_stokes.rot(phi_photon_to_obs);

    pp_surface->setStokesVector(tmp_stokes);
}

void CSurface::dePolarization(photon_package * pp)
{
    // Get the Stokes vector of the current photon package
    StokesVector tmp_stokes = *pp->getStokesVector();

    // Depolarize Stokes vector
    tmp_stokes.depolarize();
    pp->setStokesVector(tmp_stokes);
}

void CSurface::specularPolarization(photon_package * pp)
{
    // Get the Stokes vector of the current photon package
    StokesVector tmp_stokes = *pp->getStokesVector();
    double stokes_1_bak = tmp_stokes.I();

    Matrix2D mat = CMathFunctions::getReflectionMatrix(refractive_ratio, cos_i);

    if(stokes_1_bak > 1e-200) {
        // Multiply Stokes vector with scattering matrix
        tmp_stokes *= mat;
        // Normalize Stokes vector to preserve total intensity
        tmp_stokes *= stokes_1_bak / tmp_stokes.I();
    } else {
        tmp_stokes.clear();
    }

    pp->setStokesVector(tmp_stokes);
}

void CSurface::readSurfaceParameter(string path)
{
    // TODO: input file for surface parameters
    // instead of defining a homogeneous surface
    // use this file to define an inhomogeneous surface
}

void CSurface::setRandomWaveNormal(photon_package * pp, CRandomGenerator * rand_gen)
{
    Vector3D dir = pp->getDirection();
    double cos_0_tmp = -1*dir * surface_normal;
    double phi, theta;
    uint run_counter = 0;

    do {
        // for wave facet probability distribution
        // see Cox & Munk 1954, JOSA 44, 838
        // see Nakajima 1983, JQSRT 29, 521
        
        // get random water surface normal (wave normal vector)
        phi = PIx2 * rand_gen->getRND();
        theta = atan( sqrt( -sigma_sq * log(rand_gen->getRND()) ) );
        // alternative method
        // theta = acos( 1.0 / sqrt( 1.0 - sigma_sq * log(rand_gen.getRND()) ) );

        // get vector in the plane of surface
        Vector3D ref_vector = cross(dir, surface_normal);

        // rotate around phi and theta to get wave normal vector
        ref_vector.rot(surface_normal, phi);
        wave_normal = surface_normal;
        wave_normal.rot(ref_vector, theta);

        run_counter++;
    // find another normal if the angle between the wave normal vector and incoming radiation is >= PI/2
    } while(-1*dir * wave_normal <= 0.0 && run_counter < 1000);

    if(run_counter == 1000) {
        cout << "\nERROR: No wave normal found!" << endl;
    }

    // for probability of interaction of the photon and facet, i.e. weighting
    // see Plass et al. 1975, ApOpt 14, 1924
    // see Zeisse 1995, JOSAA 12, 2022
    double cos_n = wave_normal * surface_normal;
    cos_i = -1*dir * wave_normal;
    *pp->getStokesVector() *= cos_i / (cm_norm.getValue(cos_0_tmp) * cos_n);
}

bool CSurface::surfaceInteraction(photon_package * pp, CRandomGenerator * rand_gen)
{   
    // return true for an additional surface interaction, else false
    int res = -1;

    switch(surface_refl_model_loc) {
        case LAMBERTIAN:
        {
            res = lambertianReflection(pp, rand_gen);
            break;
        }
        case LOMMELSEELIGER:
        {
            res = lommelSeeligerReflection(pp, rand_gen);
            break;
        }
        case OCEAN:
        {
            setRandomWaveNormal(pp, rand_gen);
            res = specularInteraction(pp, rand_gen);
            break;
        }
    }

    switch(surface_pol_model_loc) {
        case DEPOL:
        {
            dePolarization(pp);
            break;
        }
        case SPECULAR:
        {
            specularPolarization(pp);
            break;
        }
    }

    if(!pp->getStokesVector()->isConsistent()) {
        pp->getStokesVector()->depolarize();
    }

    // res = 1 reflected
    Vector3D dir = pp->getDirection();
    cos_r = dir * surface_normal;

    if(res == 1) { // reflected
        if(cos_r > 0.0) {
            // no more surface interaction needed
            if(surface_refl_model_loc == OCEAN) {
                // scale photon package with shadowing function
                *pp->getStokesVector() *= CMathFunctions::getShadowing(cos_0, cos_r, sigma_sq);
            }
            return false;
        } else {
            // do an additional interaction
            return true;
        }
    } else {
        pp->getStokesVector()->clear();
        return false;
    }
}

int CSurface::lambertianReflection(photon_package * pp, CRandomGenerator * rand_gen)
{
    Vector3D d_in = pp->getDirection();

    StokesVector stokes = *pp->getStokesVector();
    // Calculate random phi and theta for lambertian surface reflection
    double phi =  PIx2 * rand_gen->getRND();
    double cos_2 = 1.0 - 2.0 * rand_gen->getRND();
    double theta = 0.5 * acos(cos_2);
    // alternative method
    // double theta = acos(sqrt(rand_gen->getRND()));

    pp->rotateToVector(surface_normal);
    pp->updateCoordSystem(phi, theta);

    Vector3D d_out = pp->getDirection();
    Vector3D w_normal = d_out - d_in;
    w_normal.normalize();
    cos_r = d_out * surface_normal;
    cos_i = -1*d_in * w_normal;

    // Multiply stokes vector by the surface albedo and depolarize light
    stokes *= surface_refl_param[0];

    pp->setStokesVector(stokes);
    return 1;
}

int CSurface::lommelSeeligerReflection(photon_package * pp, CRandomGenerator * rand_gen)
{
    Vector3D d_in = pp->getDirection();

    StokesVector stokes = *pp->getStokesVector();
    // Calculate random phi and theta for surface reflection
    double phi =  PIx2 * rand_gen->getRND();

    cos_r = CMathFunctions::findRootBrent(0.0, 1.0, &CMathFunctions::getLommelSeeligerIntegral, {cos_0, rand_gen->getRND()});
    double theta = acos(cos_r);

    pp->rotateToVector(surface_normal);
    pp->updateCoordSystem(phi, theta);

    Vector3D d_out = pp->getDirection();
    Vector3D w_normal = d_out - d_in;
    w_normal.normalize();
    cos_i = -1*d_in * w_normal;

    // Multiply stokes vector by the surface albedo
    stokes *= 0.5 * surface_refl_param[0] * (1.0 - cos_0 * log(1.0 + 1.0 / cos_0));

    pp->setStokesVector(stokes);
    return 1;
}

int CSurface::specularInteraction(photon_package * pp, CRandomGenerator * rand_gen)
{
    Vector3D d_in = pp->getDirection();
    StokesVector stokes = *pp->getStokesVector();

    cos_i = -1*d_in * wave_normal;
    double fresnelRefl = CMathFunctions::effFresnelReflection(refractive_ratio, cos_i);

    // reflection
    d_in.reflect(wave_normal);
    pp->rotateToVector(d_in);

    // Multiply stokes vector by the albedo
    stokes *= fresnelRefl;

    Vector3D d_out = pp->getDirection();
    cos_r = d_out * surface_normal;
    pp->setStokesVector(stokes);
    return 1;
}

void CSurface::printParameters() const
{
    cout << SEP_LINE;
    switch(surface_refl_model) {
        case LAMBERTIAN:
        {
            cout << "- Reflection model: Lambertian" << endl;
            cout << "  - Albedo: " << surface_refl_param[0] << endl;
            break;
        }
        case LOMMELSEELIGER:
        {
            cout << "- Reflection model: Lommel-Seeliger" << endl;
            cout << "  - Single scattering albedo: " << surface_refl_param[0] << endl;
            break;
        }
        case OCEAN:
        {
            cout << "- Reflection model: Ocean" << endl;
            cout << "  - Wind speed: " << surface_refl_param[1] << " [m/s]," << endl;
            cout << "  - White caps: " << 100*white_cap_fraction << " [%] with albedo: " << surface_refl_param[0] << "," << endl;
            cout << "  - Ratio of refractive index: " << surface_refl_param[2] << "," << endl;
            break;
        }
    }

    switch(surface_pol_model) {
        case DEPOL:
        {
            cout << "- Polarization model: Depolarization" << endl;
            break;
        }
        case SPECULAR:
        {
            cout << "- Polarization model: Specular (Ocean)" << endl;
            break;
        }
    }
}
