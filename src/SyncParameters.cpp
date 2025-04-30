/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "SyncParameters.hpp"

void syn_param::scale()
{
    j_I *= syn_SI_int;
    j_Q *= syn_SI_int;
    j_V *= syn_SI_int;

    alpha_I *= syn_SI_abs;
    alpha_Q *= syn_SI_abs;
    alpha_V *= syn_SI_abs;

    kappa_Q *= syn_SI_abs;
    kappa_V *= syn_SI_abs;
}

Matrix2D syn_param::getSyncMatrix()
{
    //------------------------------------------------------
    // (  alpha_I  alpha_Q  0        alpha_V  )
    // (  alpha_Q  alpha_I  kappa_V  0        )
    // (  0       -kappa_V  alpha_I  kappa_Q  )
    // (  alpha_V  0       -kappa_Q  alpha_I  )
    //------------------------------------------------------
    Matrix2D cont_matrix(4, 4);
    cont_matrix.setValue(0, 0, alpha_I);
    cont_matrix.setValue(1, 1, alpha_I);
    cont_matrix.setValue(2, 2, alpha_I);
    cont_matrix.setValue(3, 3, alpha_I);

    cont_matrix.setValue(0, 1, alpha_Q);
    cont_matrix.setValue(1, 0, alpha_Q);

    cont_matrix.setValue(0, 3, alpha_V);
    cont_matrix.setValue(3, 0, alpha_V);

    cont_matrix.setValue(2, 3, kappa_Q);
    cont_matrix.setValue(3, 2, -kappa_Q);

    cont_matrix.setValue(1, 2, kappa_V);
    cont_matrix.setValue(2, 1, -kappa_V);

    return cont_matrix;
}

syn_param syn_param::operator+(const syn_param & rhs)
{
    return syn_param(j_I + rhs.j_I,
                     j_Q + rhs.j_Q,
                     j_V + rhs.j_V,
                     alpha_I + rhs.alpha_I,
                     alpha_Q + rhs.alpha_Q,
                     alpha_V + rhs.alpha_V,
                     kappa_Q + rhs.kappa_Q,
                     kappa_V + rhs.kappa_V);
}
