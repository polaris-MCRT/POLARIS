/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "FreeFree.hpp"

#include "CCfits/FITS.h"
#include "CCfits/FITSUtilT.h"
#include "CCfits/FitsError.h"
#include "CCfits/KeyData.h"
#include "CCfits/PHDU.h"
#include "CCfits/PHDUT.h"
using namespace CCfits;

bool CFreeFree::loadGauntFITS(const string& filename)
{
    try
    {
        unique_ptr<FITS> pInfile;
        pInfile.reset(new FITS(filename, Read, true));

        PHDU& primary = pInfile->pHDU();
        long axes = primary.axes();
        if(axes != 3)
        {
            cout << "ERROR: Expected 3D data in primary HDU of Gaunt fits!" << endl;
            return false;
        }

        N_gam2 = static_cast<int>(primary.axis(0));
        N_u = static_cast<int>(primary.axis(1));
        N_Z = static_cast<int>(primary.axis(2));

        valarray<float> data;
        primary.read(data);


        gaunt_table = new float**[N_Z];
        for (int z = 0; z < N_Z; ++z)
        {
            gaunt_table[z] = new float*[N_u];
            
            for (int y = 0; y < N_u; ++y)
            {
                gaunt_table[z][y] = new float[N_gam2];
            }
        }

        size_t index = 0;
        for (int z = 0; z < N_Z; ++z)
            for (int y = 0; y < N_u; ++y)
                for (int x = 0; x < N_gam2; ++x)
                    gaunt_table[z][y][x] = data[index++];

        readExtension<float>(pInfile.get(), "LOG_U", log_u, N_u);
        readExtension<float>(pInfile.get(), "LOG_GAM2", log_gam2, N_gam2);
        
        short* Z_short = nullptr;
        int tmp_len;
        readExtension<short>(pInfile.get(), "Z_VALUES", Z_short, tmp_len);
        Z_vals = new float[tmp_len];
        for (int i = 0; i < tmp_len; ++i)
        {
            Z_vals[i] = static_cast<float>(Z_short[i]);
        }
        delete[] Z_short;
        
        Z_min = Z_vals[0];
        Z_max = Z_vals[N_Z-1];
        
        gam2_min = log_gam2[0];
        gam2_max = log_gam2[N_gam2-1];
        
        u_min  = log_u[0];
        u_max = log_u[N_u-1];


        return true;
    }
    catch (CCfits::FitsException& e)
    {
        cout << "FITS read error: " << e.message() << endl;
        return false;
    }
}

float CFreeFree::interpolate(float Z, float log_gam2_in, float log_u_in) const
{    
    if(Z<Z_min || Z>Z_max)
    {
        cout << "ERROR: Interpolation of Gaunt factor is not possible! Z is out if range!" << endl << flush;
        return -1;
    }
    
    if(log_gam2_in<gam2_min || log_gam2_in>gam2_max)
    {
        cout << "ERROR: Interpolation of Gaunt factor is not possible! log_gam2 is out if range!" << endl << flush;
        return -1;
    }
    
    if(log_u_in<u_min || log_u_in>u_max)
    {
        cout << "ERROR: Interpolation of Gaunt factor is not possible! log_u is out if range!" << endl << flush;
        return -1;
    }

    int z_idx = findIndex(Z, Z_vals, N_Z);
    float z0 = Z_vals[z_idx];
    float z1 = Z_vals[z_idx + 1];
    float fz = (Z - z0) / (z1 - z0);

    float g00 = bilinearInterp(gaunt_table[z_idx], log_gam2_in, log_u_in);
    float g01 = bilinearInterp(gaunt_table[z_idx + 1], log_gam2_in, log_u_in);

    return lerp(g00, g01, fz);
}

void CFreeFree::printParameters()
{
    cout << "Free-Free parameters" << endl << flush;
}

inline double CFreeFree::log10_u(double lambda, double Te)
{
    // check for non-physical inputs
    if(lambda <= 0.0)
    {
        cout << "ERROR: Non positive wavelengths!" << endl << flush;
        return 0;
    }
    
    if(Te <= 0.0)
    {
        cout << "ERROR: Non positive electron temperature!" << endl << flush;
        return 0;
    }
        
    const double u = (con_h * con_c) / (lambda * con_kB * Te);
    return log10(u);
} 

inline double CFreeFree::log10_gamma2(int Z, double Te)
{
    if(Z <= 1)
    {
        cout << "ERROR: Wrong charge number!" << endl << flush;
        return 0;
    }
    
    if(Te <= 0.0)
    {
        cout << "ERROR: Non positive electron temperature!" << endl << flush;
        return 0;
    }
    
    const double g2   = static_cast<double>(Z) * static_cast<double>(Z) * (con_RyJ / (con_kB * Te));
    return log10(g2);
}

inline double CFreeFree::j_lambda(double lambda, double Te, double ne, double ni, int Z, double gff)
{
    if(lambda * Te * ne  * ni * Z  * gff <= 0.0)
    {
        cout << "ERROR: Non-physical units in j_lambda calculation!" << endl << flush;
        return 0.0;
    }        

    // Common temperature factor
    const double temp_pref = sqrt( (PIx2) / (3.0 * con_kB * con_m_e) );

    // Emissivity prefactor
    const double J_PREF = (16.0 * PI * pow(con_e, 6)) /
        (3.0 * pow(con_m_e, 2) * pow(con_c, 2)) * temp_pref;

    const double inv_sqrt_T = 1.0 / sqrt(Te);
    const double x          = (con_h * con_c) / (lambda * con_kB * Te); // hc/(λ kB Te)
    const double expo       = exp(-x);
    const double lambda2    = lambda * lambda;

    const double pre = J_PREF * double(Z) * double(Z) * ne * ni * inv_sqrt_T;
    return pre * (expo / lambda2) * gff; // [W m^-4]
}

inline double CFreeFree::alpha_lambda(double lambda, double Te, double ne, double ni, int Z, double gff)
{
    if(lambda <= 0.0 || Te <= 0.0 || ne <= 0.0 || ni <= 0.0 || Z == 0 || gff <= 0.0) return 0.0;
    {
        cout << "ERROR: Non-physical units in alpha_lambda calculation!" << endl << flush;
        return 0.0;
    }

    // Common temperature factor
    const double temp_pref = sqrt( (PIx2) / (3.0 * con_kB * con_m_e) );

    // Absorption prefactor
    const double A_PREF = (8.0 * PI * pow(con_e, 6)) /
        (3.0 * con_h * pow(con_c, 4) * pow(con_m_e, 2)) * temp_pref;

    const double inv_sqrt_T = 1.0 / sqrt(Te);
    const double x          = (con_h * con_c) / (lambda * con_kB * Te); // hc/(λ kB Te)
    const double one_minus_exp = -expm1(-x); // = 1 - exp(-x), numerically stable
    const double lambda3    = lambda * lambda * lambda;

    const double pre = A_PREF * double(Z) * double(Z) * ne * ni * inv_sqrt_T;
    return pre * (lambda3 * one_minus_exp) * gff; // [m^-1]
}

inline int CFreeFree::findIndex(float value, const float* axis, int length) const
{
    int low = 0;
    int high = length - 2;
    
    while (low <= high)
    {
        int mid = (low + high) / 2;
        
        if(axis[mid] <= value && value < axis[mid + 1])
        {
            return mid;
        }
        else if(value < axis[mid])
        {
            high = mid - 1;
        }
        else
        {
            low = mid + 1;
        }
    }
    return length - 2;
}

inline float CFreeFree::lerp(float a, float b, float t) const
{
    return a + t * (b - a);
}   

inline float CFreeFree::bilinearInterp(float** plane, float xval, float yval) const
{
    int ix = findIndex(xval, log_gam2, N_gam2);
    int iy = findIndex(yval, log_u, N_u);

    float fx = (xval - log_gam2[ix]) / (log_gam2[ix + 1] - log_gam2[ix]);
    float fy = (yval - log_u[iy]) / (log_u[iy + 1] - log_u[iy]);

    float v00 = plane[iy][ix];
    float v10 = plane[iy][ix + 1];
    float v01 = plane[iy + 1][ix];
    float v11 = plane[iy + 1][ix + 1];

    float vx0 = lerp(v00, v10, fx);
    float vx1 = lerp(v01, v11, fx);
    
    return lerp(vx0, vx1, fy);
}