/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef SYNC_PARAMETERS_H
#define SYNC_PARAMETERS_H

#include "Matrix2D.hpp"
#include "Typedefs.hpp"

// redefine some natural constants in cgs
#define syn_me 9.1093826e-28 // electron mass      [g]
#define syn_e 4.80320680e-10 // electron charge    [Statcoulomb]
#define syn_h 6.6260693e-27  // Planck constant    [erg/s]
#define syn_kB 1.380662e-16  // Boltzmann constant [erg/K]
#define syn_c 2.99792458e10  // speed of light     [cm/s]
// conversion factor 100.0 (1/cm -> 1/m) * 1e-3 (erg s^-1 cm^-2 Hz^-1 -> W m^-2 Hz^-1) 
#define syn_SI_int 0.1
// conversion factor 100.0 (1/cm -> 1/m) 
#define syn_SI_abs 100.0

// container class for the parameters of sync. RT
class syn_param
{
public:
    syn_param()
    {
        j_I = 0;
        j_Q = 0;
        j_V = 0;

        alpha_I = 0;
        alpha_Q = 0;
        alpha_V = 0;

        kappa_Q = 0;
        kappa_V = 0;
    }

    syn_param(double _j_I,
              double _j_Q,
              double _j_V,
              double _alpha_I,
              double _alpha_Q,
              double _alpha_V,
              double _kappa_Q,
              double _kappa_V)
    {
        j_I = _j_I;
        j_Q = _j_Q;
        j_V = _j_V;

        alpha_I = _alpha_I;
        alpha_Q = _alpha_Q;
        alpha_V = _alpha_V;

        kappa_Q = _kappa_Q;
        kappa_V = _kappa_V;
    }

    syn_param operator+(const syn_param & rhs);

    // (back) conversion into SI
    void scale();

    // generate matrix for synchrotron RT
    Matrix2D getSyncMatrix();

    double j_I;
    double j_Q;
    double j_V;

    double alpha_I;
    double alpha_Q;
    double alpha_V;

    double kappa_Q;
    double kappa_V;
};

#endif /* SYNC_PARAMETERS_H */
