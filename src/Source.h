#pragma once
#include "Dust.h"
#include "Grid.h"
#include "Vector.h"
#include "chelper.h"
#include "typedefs.h"

#ifndef CSOURCE
#define CSOURCE

class CSourceBasic
{
  public:
    CSourceBasic(void)
    {
        L = 0;
        T = 0;
        R = 0;

        q = 0;
        u = 0;

        is_ext = false;

        nr_of_photons = 0;

        dust = 0;
        grid = 0;

        StringID = "Basic radiation source";
    }

    virtual ~CSourceBasic(void)
    {}

    Vector3D getPosition()
    {
        return pos;
    }

    double getRadius()
    {
        return R;
    }

    double getLuminosity()
    {
        return L;
    }

    double getTemperature()
    {
        return T;
    }

    virtual double getSublimationRadius()
    {
        return 0;
    }

    void setParameter(parameters & param, CGridBasic * _grid, CDustMixture * _dust, uint p)
    {
        setGrid(_grid);
        setDust(_dust);

        // Get global wavelength grid
        wavelength_list = _dust->getWavelengthList();

        setParameter(param, p);
    }

    virtual void setParameter(parameters & param, uint p)
    {}

    bool setParameterFromFile(parameters & param, CGridBasic * _grid, CDustMixture * _dust, uint p)
    {
        setGrid(_grid);
        setDust(_dust);

        // Get global wavelength grid
        wavelength_list = _dust->getWavelengthList();

        return setParameterFromFile(param, p);
    }

    virtual bool setParameterFromFile(parameters & param, uint p)
    {
        return false;
    }

    void setGrid(CGridBasic * _grid)
    {
        grid = _grid;
    }

    void setDust(CDustMixture * _dust)
    {
        dust = _dust;
    }

    virtual bool initSource(uint id, uint max, bool use_energy_density)
    {
        return true;
    }

    virtual bool initSource(uint w)
    {
        return true;
    }

    virtual void clean()
    {}

    string getStringID()
    {
        return StringID;
    }

    void setNrOfPhotons(llong val)
    {
        nr_of_photons = val;
    }

    void updateNrOfPhotons(double val)
    {
        nr_of_photons = llong(nr_of_photons * val);
    }

    virtual void setSideLength(double val)
    {}

    virtual void createNextRay(photon_package * pp, llong i_pos, llong nr_photons = 0)
    {}

    virtual llong getNrOfPhotons()
    {
        return nr_of_photons;
    }

    uint getNrOfWavelength()
    {
        return wavelength_list.size();
    }

    virtual void setOrientation(Vector3D n1, Vector3D n2, double _theta, double _phi)
    {}

    virtual StokesVector getStokesVector(photon_package * pp)
    {
        return StokesVector(0, 0, 0, 0);
    }

    virtual uint getBins()
    {
        return 1;
    }

  protected:
    CDustMixture * dust;
    CGridBasic * grid;

    llong nr_of_photons;
    dlist wavelength_list;

    Vector3D pos;
    double R;
    double L;
    double T;

    double q;
    double u;

    spline lam_pf;
    spline sp_ext;
    spline sp_ext_q;
    spline sp_ext_u;

    bool is_ext;

    string StringID;
};

class CSourceStar : public CSourceBasic
{
  public:
    CSourceStar()
    {
        pos = 0;
        StringID = "point source";
    }

    ~CSourceStar()
    {}

    bool initSource(uint id, uint max, bool use_energy_density);

    void createNextRay(photon_package * pp, llong i_pos, llong nr_photons = 0);

    bool setParameterFromFile(parameters & param, uint p);
    void setParameter(parameters & param, uint p)
    {
        dlist values = param.getPointSources();

        pos = Vector3D(values[p], values[p + 1], values[p + 2]);
        R = values[p + 3];
        T = values[p + 4];

        q = values[p + 5];
        u = values[p + 6];

        nr_of_photons = (llong)values[p + NR_OF_POINT_SOURCES - 1];

        L = PIx4 * con_sigma * (R * R_sun) * (R * R_sun) * T * T * T * T;
    }
};

class CSourceStarField : public CSourceBasic
{
  public:
    CSourceStarField(void)
    {
        pos = 0;
        var = 0;
        StringID = "diffuse source";
    }

    ~CSourceStarField(void)
    {}

    bool initSource(uint id, uint max, bool use_energy_density);

    void createNextRay(photon_package * pp, llong i_pos, llong nr_photons = 0);

    bool setParameterFromFile(parameters & param, uint p);
    void setParameter(parameters & param, uint p)
    {
        dlist values = param.getDiffuseSources();

        pos = Vector3D(values[p + 0], values[p + 1], values[p + 2]);
        R = values[p + 3];
        T = values[p + 4];
        var = values[p + 5];

        q = values[p + 6];
        u = values[p + 7];

        nr_of_photons = llong(values[p + NR_OF_DIFF_SOURCES - 1]);

        L = R * T * T * T * T;
    }

  private:
    double var;
};

class CSourceBackground : public CSourceBasic
{
  public:
    CSourceBackground()
    {
        init = false;
        constant = false;

        step_xy = 0;
        off_xy = 0;

        grid = 0;

        bins = 0;
        sidelength = 0;

        ex.set(1, 0, 0);
        ey.set(0, 1, 0);
        ez.set(0, 0, 1);

        c_temp = 0;
        c_f = 0;
        c_q = 0;
        c_u = 0;
        c_v = 0;

        max_len = 0;

        lam_pf = 0;
        L = 0;
        rot_angle1 = 0;
        rot_angle2 = 0;

        StringID = "background source";
    }

    ~CSourceBackground()
    {
        if(L != 0)
            delete[] L;
        if(lam_pf != 0)
            delete[] lam_pf;
    }

    bool initSource(uint id, uint max, bool use_energy_density);

    StokesVector getStokesVector(photon_package * pp);

    bool setParameterFromFile(parameters & param, uint p);
    void setParameter(parameters & param, uint p)
    {
        dlist values;
        if(param.getNrOfBackgroundSources() == 0)
            values.resize(10, 0);
        else
            values = param.getBackgroundSources();

        bins = 1;
        max_len = 1;

        c_f = values[p + 0];
        c_temp = values[p + 1];
        c_q = values[p + 2];
        c_u = values[p + 3];
        c_v = values[p + 4];
        rot_angle1 = values[p + 5];
        rot_angle2 = values[p + 6];

        constant = true;
        init = true;

        nr_of_photons = llong(values[p + 7]);
    }

    llong getNrOfPhotons()
    {
        return llong(bins * bins) * nr_of_photons;
    }

    uint getBins()
    {
        return bins;
    }

    void setOrientation(Vector3D n1, Vector3D n2, double _rot_angle1, double _rot_angle2)
    {
        rot_angle1 = _rot_angle1;
        rot_angle2 = _rot_angle2;

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

  private:
    Matrix2D temp, f, q, u, v;
    Vector3D ex, ey, ez;
    uint bins;
    uint max_len;
    bool init;
    bool constant;
    double sidelength;

    double step_xy;
    double off_xy;

    double c_temp, c_f, c_q, c_u, c_v;
    double rot_angle1, rot_angle2;
    spline * lam_pf;
    double * L;
};

class CSourceISRF : public CSourceBasic
{
  public:
    CSourceISRF()
    {
        init = false;
        kill_count = 0;
        sidelength = 0;
        g_zero = 0;

        grid = 0;

        c_w = 0;
        c_f = 0;
        c_q = 0;
        c_u = 0;
        c_v = 0;

        L = 0;

        StringID = "isrf source";
    }

    ~CSourceISRF()
    {
        if(c_w != 0)
            delete[] c_w;
        if(c_f != 0)
            delete[] c_f;
    }

    bool initSource(uint id, uint max, bool use_energy_density);

    bool setParameterFromFile(parameters & param, uint p);

    void setParameter(parameters & param, uint p)
    {
        nr_of_photons = param.getNrOfISRFPhotons();
        g_zero = param.getISRFGZero();

        sp_ext.resize(getNrOfWavelength());
        for(uint w = 0; w < getNrOfWavelength(); w++)
        {
            // Get Mathis radiation field
            double rad_field = CMathFunctions::mathis_isrf(wavelength_list[w]);

            // Calculate final emission
            sp_ext.setValue(w,
                            wavelength_list[w],
                            g_zero * rad_field * dust->getForegroundExtinction(wavelength_list[w]));
        }

        sp_ext.createSpline();
    }

    void createNextRay(photon_package * pp, llong i_pos, llong nr_photons = 0);

  private:
    Vector3D e, l;

    int kill_count;
    bool init;
    double c_q, c_u, c_v;
    double sidelength, g_zero;
    double *c_w, *c_f;
};

class CSourceDust : public CSourceBasic
{
  public:
    CSourceDust(void)
    {
        total_energy = 0;
        cell_prob = 0;

        StringID = "dust source";
    }

    ~CSourceDust(void)
    {
        if(cell_prob != 0)
            delete[] cell_prob;
    }

    bool initSource(uint w);

    bool initSource(uint id, uint max, bool use_energy_density);

    void createNextRay(photon_package * pp, llong i_pos, llong nr_photons = 0);

    void setParameter(parameters & param, uint p)
    {
        nr_of_photons = param.getNrOfDustPhotons();
    }

    llong getNrOfPhotons()
    {
        return nr_of_photons;
    }

  private:
    double * total_energy;
    prob_list * cell_prob;
};

typedef vector<CSourceBasic *> slist;
#endif
