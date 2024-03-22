#pragma once
#include "Photon.h"
#include "Vector.h"
#include "Parameters.h"
#include "MathFunctions.h"
#include "Typedefs.h"
#include "Grid.h"

// Surface parameter for reflection (dlist surface_refl_param)
// LAMBERTIAN:        [surface albedo]
// LOMMELSEELIGER:    [single scattering albedo]
// OCEAN:             [white caps albedo,
//                     wind speed (m/s),
//                     ratio of refractive index (n2/n1)]
// 
// Surface parameter for polarization (dlist surface_pol_param)
// DEPOLARIZATION:    []
// OCEAN -> SPECULAR: []
//
//
// LAMBERTIAN / LOMMELSEELIGER:
// see Lester et al. 1979, JRASC 73, 233
//
// OCEAN:
// see Zhai et al. 2010, JQSRT 111, 1025
//     Trees & Stam 2019, A&A 626, 129


class CSurface
{
    public:
        CSurface()
        {
            surface_refl_model = LAMBERTIAN;
            surface_refl_param.push_back(0);
            surface_pol_model = DEPOL;
        }

        void initSurface(parameters & param);

        void initOceanParam();

        void setLocalParam(photon_package * pp, CRandomGenerator * rand_gen);

        void getEscapePhotonSurface(
            photon_package * pp,
            Vector3D obs_ex,
            Vector3D dir_obs,
            photon_package * pp_surface);

        void getEscapePhotonSurfaceLambertian(
            photon_package * pp,
            Vector3D dir_obs,
            photon_package * pp_surface);
        
        void getEscapePhotonSurfaceLommelSeeliger(
            photon_package * pp,
            Vector3D dir_obs,
            photon_package * pp_surface);

        void getEscapePhotonSurfaceOcean(
            photon_package * pp,
            Vector3D dir_obs,
            photon_package * pp_surface);

        void dePolarization(photon_package * pp);

        void specularPolarization(photon_package * pp);

        bool surfaceInteraction(photon_package * pp, CRandomGenerator * rand_gen);

        int lambertianReflection(photon_package * pp, CRandomGenerator * rand_gen);

        int lommelSeeligerReflection(photon_package * pp, CRandomGenerator * rand_gen);

        int specularInteraction(photon_package * pp, CRandomGenerator * rand_gen);

        // TODO: input file for surface parameters
        // instead of defining a homogeneous surface
        // use this file to define an inhomogeneous surface
        void readSurfaceParameter(string path);

        void setRandomWaveNormal(photon_package * pp, CRandomGenerator * rand_gen);

        void printParameters() const;

        double getSurfaceReflModel() const
        {
            return surface_refl_model;
        }

        double getSurfacePolModel() const
        {
            return surface_pol_model;
        }

        double getSurfaceReflParam(uint v) const
        {
            return surface_refl_param[v];
        }

        double getSurfacePolParam(uint v) const
        {
            return surface_pol_param[v];
        }

        dlist getSurfaceReflParam() const
        {
            return surface_refl_param;
        }

        dlist getSurfacePolParam() const
        {
            return surface_pol_param;
        }

        interp getCMNorm() const
        {
            return cm_norm;
        }

        void setSurfaceReflModel(double _surface_refl_model)
        {
            surface_refl_model = _surface_refl_model;
            surface_refl_model_loc = _surface_refl_model;
        }

        void setSurfacePolModel(double _surface_pol_model)
        {
            surface_pol_model = _surface_pol_model;
            surface_pol_model_loc = _surface_pol_model;
        }

        void setSurfaceReflParam(dlist _surface_refl_param)
        {
            surface_refl_param = _surface_refl_param;
        }

        void setSurfacePolParam(dlist _surface_pol_param)
        {
            surface_pol_param = _surface_pol_param;
        }

        void setCMNorm(interp _cm_norm)
        {
            cm_norm = _cm_norm;
        }

    private:
        uint surface_refl_model, surface_pol_model;
        uint surface_refl_model_loc, surface_pol_model_loc;

        dlist surface_refl_param, surface_pol_param;

        double sigma_sq, white_cap_fraction;
        double cos_0, cos_r, cos_i;

        Vector3D surface_normal, wave_normal;

        interp cm_norm;
};
