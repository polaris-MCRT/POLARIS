/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CSOURCE_STARFIELD_H
#define CSOURCE_STARFIELD_H

#include "DustMixture.hpp"
#include "Vector3D.hpp"
#include "MathFunctions.hpp"
#include "Matrix2D.hpp"
#include "Parameters.hpp"
#include "SourceBasic.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"

class CSourceStarField : public CSourceBasic
{
public:
    CSourceStarField(void)
    {
        pos = 0;
        
        sig_x = 0;
        sig_y = 0;
        sig_z = 0;
    
        a = 0;
        b = 0;
        c = 0;
    
        rot1 = Vector3D(1,0,0);
        rot1 = Vector3D(0,1,0);
    
        ang1 = 0;
        ang2 = 0;
        
        source_id = SRC_SFIELD;
        
        Npos=0;
        dist_pos=0;
    }

    ~CSourceStarField(void)
    {
        if(dist_pos!=0)
        {
            delete [] dist_pos;
            dist_pos=0;
        }
    }

    bool initSource(uint id, uint max, bool use_energy_density);

    void createNextRay(photon_package * pp, CRandomGenerator * rand_gen);
    void createDirectRay(photon_package * pp, CRandomGenerator * rand_gen, Vector3D dir_obs);

    bool setParameterFromFile(parameters & param, uint p);
    
    bool read_fits_file(string & filename);
    
    void setParameter(parameters & param, uint p);

private:
    double sig_x;
    double sig_y;
    double sig_z;
    
    double a,b,c;
    
    Vector3D rot1;
    Vector3D rot2;
    
    double ang1, ang2;
    
    long Npos;
    Vector3D * dist_pos;
};

#endif /* CSOURCE_STARFIELD_H */