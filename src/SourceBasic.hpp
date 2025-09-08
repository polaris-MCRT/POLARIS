/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CSOURCE_BASIC_H
#define CSOURCE_BASIC_H

#include "DustMixture.hpp"
#include "Vector3D.hpp"
#include "MathFunctions.hpp"
#include "Matrix2D.hpp"
#include "Parameters.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"

class CSourceBasic
{
public:
    CSourceBasic(void)
    {
        L = 0;
        T = 0;
        R = 0;
        
        r_sub=-1;

        q = 0;
        u = 0;

        is_ext = false;

        nr_of_photons = 0;

        dust = 0;
        grid = 0;

        source_id = SRC_BASIC;
    }

    virtual ~CSourceBasic(void)
    {}

    Vector3D getPosition();

    double getRadius();

    double getLuminosity();

    double getTemperature();

    virtual double getSublimationRadius();

    void setParameter(parameters & param, CGridBasic * _grid, CDustMixture * _dust, uint p);

    virtual void setParameter(parameters & param, uint p);

    bool setParameterFromFile(parameters & param, CGridBasic * _grid, CDustMixture * _dust, uint p);

    virtual bool setParameterFromFile(parameters & param, uint p);

    void setGrid(CGridBasic * _grid);

    void setDust(CDustMixture * _dust);

    virtual bool initSource(uint id, uint max, bool use_energy_density = false);

    virtual bool initSource(uint w);

    virtual void clean();

    uint getID();

    void setNrOfPhotons(ullong val);

    void updateNrOfPhotons(double val);

    virtual void setSideLength(double val);

    virtual void createNextRay(photon_package * pp, CRandomGenerator * rand_gen);

    virtual void createNextRayToCell(photon_package * pp,
                                     CRandomGenerator * rand_gen,
                                     ulong i_cell,
                                     bool cell_as_border = false);

    virtual void createDirectRay(photon_package * pp, CRandomGenerator * rand_gen, Vector3D dir_obs = Vector3D());

    virtual ullong getNrOfPhotons();

    uint getNrOfWavelength();

    virtual void setOrientation(Vector3D n1, Vector3D n2, double _theta, double _phi);

    virtual StokesVector getStokesVector(photon_package * pp);

    virtual uint getBins();

protected:
    CDustMixture * dust;
    CGridBasic * grid;

    ullong nr_of_photons;
    dlist wavelength_list;

    Vector3D pos;
    double R;
    double L;
    double T;
    double r_sub;

    double q;
    double u;

    spline lam_pf;
    spline sp_ext;
    spline sp_ext_q;
    spline sp_ext_u;

    bool is_ext;

    uint source_id;
};

typedef vector<CSourceBasic *> slist;

#endif /* CSOURCE_BASIC_H */
