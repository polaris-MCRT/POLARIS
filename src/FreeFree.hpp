/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CFREEFREE_HPP
#define CFREEFREE_HPP

#include "MathFunctions.hpp"
#include "MathSpline.hpp"
#include "Matrix2D.hpp"
#include "Typedefs.hpp"
#include <CCfits/CCfits>
#include "Stokes.hpp"

#include <valarray>
#include <memory>

using namespace CCfits;
using namespace std;

template <typename T>
void readExtension(FITS* pFITS, const string& extName, T*& buffer, int& length)
{
    ExtHDU& ext = pFITS->extension(extName);
    valarray<T> data;
    ext.read(data);
    length = static_cast<int>(data.size());
    buffer = new T[length];
    
    for (int i = 0; i < length; ++i)
    {
        buffer[i] = data[i];
    }
}

// class for the physics of sync. RT
class CFreeFree
{
public:
    CFreeFree()
    {
        Z_min =-1;
        Z_max =-1;
        
        gam2_min =-1;
        gam2_max =-1;
        
        u_min =-1;
        u_max =-1;
        
        N_gam2 =-1;
        N_u =-1;
        N_Z =-1;
        
        gaunt_table = 0;
        log_u = 0;
        log_gam2 = 0;
        Z_vals = 0;
        
        inter_gaunt = false;
    };

    ~CFreeFree()
    {
        if(gaunt_table)
        {
            for (int z = 0; z < N_Z; ++z)
            {
                for (int y = 0; y < N_u; ++y) 
                    delete[] gaunt_table[z][y];
                    
                delete[] gaunt_table[z];
            }
            delete[] gaunt_table;
        }
        
        if(log_u!=0)
            delete[] log_u;
            
        if(log_gam2!=0)
            delete[] log_gam2;
        
        if(Z_vals!=0)
            delete[] Z_vals;
    };
    
    bool loadGauntFITS(const string& filename);
    bool loadGauntFITS();
    float interpolate(float Z, float log_gam2_in, float log_u_in) const;
    
    void setGauntFile(string f);
    void printParameters();
    

    double getGauntFactor(double lambda, double Te, double Z);
    
    //inline double get_j_lambda(double lambda, double Te, double ne, double ni, double Z);
    //inline double get_alpha_lambda(double lambda, double Te, double ne, double ni, double Z);
    
    void get_coeff_lambda(Matrix2D & alpha, StokesVector & J, double lambda, double Te, double ne, double ni, double Z);
    
private:
    float Z_min, Z_max;
    float gam2_min, gam2_max;
    float u_min, u_max;
    
    float*** gaunt_table;
    float* log_u;
    float* log_gam2;
    float* Z_vals;
    int N_gam2, N_u, N_Z;
    
    string gaunt_file;
    bool inter_gaunt;
    
    inline double gaunt_ff_approx(double lambda, double Te, double Z);

    // ---------------- Emissivity and absorption (angle-integrated) ----------------
    //
    // j_λ  [W m^-4]:
    //   j_λ = (16π e^6 / (3 m_e^2 c^2)) * sqrt(2π / (3 kB m_e))
    //         * Z^2 n_e n_i T_e^(-1/2) * exp[-hc/(λ kB T_e)] / λ^2 * ḡ_ff(λ,T_e)
    //
    // α_λ  [m^-1]:
    //   α_λ = ( 8π e^6 / (3 h c^4 m_e^2)) * sqrt(2π / (3 kB m_e))
    //         * Z^2 n_e n_i T_e^(-1/2) * λ^3 * (1 - exp[-hc/(λ kB T_e)]) * ḡ_ff
    //
    // These are angle-integrated (no sr). Te is electron temperature.

    inline double j_lambda(double lambda, double Te, double ne, double ni, double Z, double gff);
    inline double alpha_lambda(double lambda, double Te, double ne, double ni, double Z, double gff);

    // log10 u, with u = (h c)/(lambda kB Te)    
    inline double log10_u(double lambda, double Te);

    // log10 gamma^2, with gamma^2 = Z^2 * Ry / (kB Te) and Ry = h c R_infty
    inline double log10_gamma2(double Z, double Te);

    // Gaunt factor interpolation functions 
    inline int findIndex(float value, const float* axis, int length) const;
    inline float lerp(float a, float b, float t) const;
    inline float bilinearInterp(float** plane, float xval, float yval) const;
};

#endif 

