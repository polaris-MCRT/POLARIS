/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CMATH_FUNCTIONS_H
#define CMATH_FUNCTIONS_H

#include "Matrix2D.hpp"
#include "MathProbList.hpp"
#include "Stokes.hpp"
#include "Vector3D.hpp"
#include "Typedefs.hpp"

class CMathFunctions
{
public:
    CMathFunctions(void)
    {
        b = 0;
        a = 0;
        r = 0;
        d = 0;
        y = 0;
        sigma = 0;

        nr_ofSeq = 0;
    }

    ~CMathFunctions(void)
    {
        if(b != 0)
            delete b;
        if(a != 0)
            delete a;
        if(d != 0)
            delete d;
        if(r != 0)
            delete r;
        if(y != 0)
            delete y;
        if(sigma != 0)
            delete sigma;
    }

    static bool isPowerOfTwo(int num);

    static Vector3D SolveEqSys(Matrix2D & inA, Vector3D b);

    static Vector3D gauss(Matrix2D inA, Vector3D inb);

    static void gauss(Matrix2D & A, double * b, int n);

    static double BE_Mass(double r, double rho0, double rc);

    static void CholDec(Matrix2D & A);

    void setHeader(uint pos,
                   uint seqID,
                   uint sourceID,
                   double wavelength,
                   double _rot_angle1,
                   double _rot_angle2);

    void updateSatistics(uint pos, StokesVector st);

    static int insertInList(double val, dlist & list);

    static uint biListIndexSearch(double val, const dlist & list);

    static uint biListIndexSearchRec(double val, const dlist & list);

    static uint biListIndexSearchRec(double val, const double * list, uint N);

    static uint biListIndexSearch(double val, const double * list, uint N);

    static double interpolate(double x_min, double x_max, double y_min, double y_max, double x_ipol);

    static uint inSphereTest(Vector3D a, Vector3D b, Vector3D c, Vector3D d, Vector3D e);

    static uint orientationTest(Vector3D a, Vector3D b, Vector3D c, Vector3D d);

    static double det5x5(double mat[5][5], int N = 5);

    static double det4x4(double mat[4][4], int N = 4);

    static Vector3D calcKeplerianVelocity(Vector3D pos, double stellar_mass);

    static double getClosestLinePoint(Vector3D line_start, Vector3D line_end, Vector3D point);

    static double maxValue(double v1, double v2);

    static uint maxValue(uint v1, uint v2);

    static double grad2rad(double grad);

    static double rad2grad(double rad);

    static bool areSame(double a, double b);

    static double generateGaussianNoise(const double & variance);

    static double diffSeries(double y, double upper_limit);

    static double getErfi(double x);

    void sleep(int milliseconds);

    static double integ(const double * x, const double * y, uint xlow, uint xup);

    static double integ(const dlist & x, double * y, uint xlow, uint xup);

    // for debugging reasons only
    static double integ1(dlist x, double * y, uint xlow, uint xup);

    static double integ(const dlist & x, const dlist & y, uint xlow, uint xup);

    static double integ(const dlist & x, const uilist & y, uint xlow, uint xup);

    static double integ(const double * x, dlist & y, uint xlow, uint xup);

    static double integ(double * x, dlist & y, uint xlow, uint xup);

    static double integ_dust_size(const double * a_eff,
                                         const double * quantity,
                                         uint nr_of_dust_species,
                                         double a_min,
                                         double a_max);

    static StokesVector integ_dust_size(const double * a_eff,
                                               const StokesVector * stokes,
                                               uint nr_of_dust_species,
                                               double a_min,
                                               double a_max);

    static double full_integ(dlist x, dlist y, double xlow, double xup);

    static void probListInteg(const double * x, const double * y, prob_list & integ_spline);

    static void probListInteg(dlist x, const double * y, prob_list & integ_spline);

    static double calc_delta(double B, double Td, double Tg, double ng);

    static double calc_mach(double vel, double Tg, double mu);

    static double calc_larm_limit(double B, double Td, double Tg, double ng, double s, double larm_f);

    static double planck(double l, double T);

    static double planck_hz(double f, double T);

    static double dplanck_dT(double l, double T);

    static double mathis_isrf(double wavelength);

    static void SinList(double start, double stop, double * list, uint N, double f);

    static void LinearList(double start, double stop, double * list, uint N);

    static void LinearList(double start, double stop, dlist & list);

    static dlist LinearList(double start, double stop, uint N);

    static void ExpList(double start, double stop, double * list, uint N, double base);

    static void ExpListSym(double start, double stop, double * list, uint N, double base);

    static void LogList(double start, double stop, double * list, uint N, double base);

    static void LogList(double start, double stop, dlist & list, double base);

    static void SymLogList(double start, double stop, double * list, uint N, double base);

    static double sgn(double x);

    static double sgn(int x);
    
    static double sgn(int64_t x);

    static Matrix2D getRotationMatrix(double cos_phi,
                                             double sin_phi,
                                             double cos_theta,
                                             double sin_theta);

    static double getRotationAngleObserver(const Vector3D & obs_ex,
                                                  const Vector3D & photon_ex,
                                                  const Vector3D & photon_ey);

    static void getPropMatrixAPi(double cos_theta,
                                        double sin_theta,
                                        double cos_2_phi,
                                        double sin_2_phi,
                                        double mult,
                                        Matrix2D * propMatrix);

    static void getPropMatrixBSigmaP(double cos_theta,
                                            double sin_theta,
                                            double cos_2_phi,
                                            double sin_2_phi,
                                            double mult,
                                            Matrix2D * propMatrix);

    static void getPropMatrixASigmaP(double cos_theta,
                                            double sin_theta,
                                            double cos_2_phi,
                                            double sin_2_phi,
                                            double mult,
                                            Matrix2D * propMatrix);

    static void getPropMatrixASigmaM(double cos_theta,
                                            double sin_theta,
                                            double cos_2_phi,
                                            double sin_2_phi,
                                            double mult,
                                            Matrix2D * propMatrix);

    static void getPropMatrixBSigmaM(double cos_theta,
                                            double sin_theta,
                                            double cos_2_phi,
                                            double sin_2_phi,
                                            double mult,
                                            Matrix2D * propMatrix);

    static void getPropMatrixBPi(double cos_theta,
                                        double sin_theta,
                                        double cos_2_phi,
                                        double sin_2_phi,
                                        double mult,
                                        Matrix2D * propMatrix);

    static Vector3D interpVector(Vector3D y_1, Vector3D y_2, double x_1, double x_2);

    static double lum2Jy(double intensity, double wavelength, double distance);

    static void lum2Jy(StokesVector * S, double wavelength, double distance);

    static uint findListIndex(double * list, uint l_limit, uint h_limit, double val);

    void initChiSq(uint _N, uint _M);

    static double getPolarInterp(double t, double s, double Q1, double Q2, double Q3, double Q4);

    /*
    Returns normalized scattering phase function of the Henyey-Greenstein function
    Henyey & Greenstein 1941, ApJ 93, 70
    */
    static double phaseFunctionHG(double g, double theta);

    /*
    Returns normalized scattering phase function of the Draine Henyey-Greenstein function
    Draine 2003, ApJ 598, 1017
    */
    static double phaseFunctionDHG(double g, double alpha, double theta);

    /*
    Returns normalized scattering phase function of the three parameter Henyey-Greenstein function
    Kattawar 1975, JQSRT 15, 839
    Witt 1977, ApJS 31, 1
    */
    static double phaseFunctionTTHG(double g1, double g2, double weight, double theta);

    /*
    Returns the integral of the phi scattering angle distribution minus random number
    see Eq. (A16) in Fischer et al 1994, A&A 284, 187
    */
    static double getPhiIntegral(double phi, dlist args);

    /*
    Returns the integral of the Draine Henyey-Greenstein phase function minus random number
    */
    static double getDHGIntegral(double cos_theta, dlist args);

    /*
    Returns the integral of the three parameter Henyey-Greenstein phase function minus random number
    see Eq. (20) in Witt 1977, ApJS 31, 1
    */
    static double getTTHGIntegral(double cos_theta, dlist args);

    /*
    Hybrid root-finding algorithm by Brent and Dekker
    Brent, R. P. (1973), "Chapter 4: An Algorithm with Guaranteed Convergence for Finding a Zero of a Function",
        Algorithms for Minimization without Derivatives, Englewood Cliffs, NJ: Prentice-Hall
    Dekker, T. J. (1969), "Finding a zero by means of successive linear interpolation", in Dejon, B.; Henrici, P. (eds.),
        Constructive Aspects of the Fundamental Theorem of Algebra, London: Wiley-Interscience

    Procedure returns root of function in the given interval [a, b], within tolerance 6 * macheps * |x| + 2t,
    where macheps is the relative machine precision and t is a positive tolerance.
    This procedure assumes that f(a) and f(b) have different signs.
    */
    static double findRootBrent(double a, double b, double (*func)(double, dlist), dlist args);

    static double Freq2Velo(double _f, double _f0);

    static double Velo2Freq(double _v, double _f0);

    static double getWavelengthStepWidth(double lam_min, double lam_max, double nr_of_spectral_bins);

    static void gauss(Matrix2D & inA, double * b, double * res, uint n);

    /*
    Basic Gaussian elimination with pivoting
    DESCRIPTION:
    - Algorithm for solving system of n linear equations
      with n unknowns (x1,x2,...,xn)
    - Gaussian elimination is algorithm for solving
      system of linear equations. This method can be
      used to calculate determinant, inverse matrix or
      find inverse matrix.
    Author: Ervin B,
    May 2013
    */
    static dlist gauss_solver(Matrix2D & A, uint n, dlist & b);

    /*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
    static double fmodulo(double v1, double v2);

    /*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
    static int imodulo(int v1, int v2);

    static dcomplex Csqrt(dcomplex z);

    /*
    Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
    to calculate scattering and absorption by a homogenous isotropic
    sphere.

    Comment:
        NANG = number of angles between 0 and 90 degrees
                (will calculate 2 * NANG - 1 directions from 0 to 180 deg.)

    Given:
        X = 2*pi*a/lambda
        REFREL = (complex refractive index of sphere) / (real index of medium)

    Returns:
        S1(1 .. 2 * NANG - 1) =  (incident E perpendicular to scattering plane,
                                  scattering E perpendicular to scattering plane)
        S2(1 .. 2 * NANG - 1) =  (incident E parallel to scattering plane,
                                  scattering E parallel to scattering plane)
        QEXT = C_ext/pi*a**2 = efficiency factor for extinction
        QSCA = C_sca/pi*a**2 = efficiency factor for scattering
        QBACK = 4*pi*(dC_sca/domega)/pi*a**2 = backscattering efficiency
        GSCA = <cos(theta)> for scattering

    Original program taken from Bohren and Huffman (1983), Appendix A
    Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
    in order to compute <cos(theta)>
    */
    static bool calcBHMie(double x,
                          dcomplex refractive_index,
                          double & qext,
                          double & qabs,
                          double & qsca,
                          double & gsca,
                          double * S11,
                          double * S12,
                          double * S33,
                          double * S34);

    /*
    Wolf & Voshchinnikov approximation of optical properties for spherical grains.
    see Wolf & Voshchinnikov (2004), Comput. Phys. Commun. 162, 113
    */
    static bool calcWVMie(double x,
                          dlist scat_angle,
                          dcomplex refractive_index,
                          double & qext,
                          double & qabs,
                          double & qsca,
                          double & gsca,
                          double * S11,
                          double * S12,
                          double * S33,
                          double * S34);

    static bool calcGeometricOptics(double x,
                                           dcomplex refractive_index,
                                           double & qext,
                                           double & qabs,
                                           double & qsca,
                                           double & gsca,
                                           double * S11,
                                           double * S12,
                                           double * S33,
                                           double * S34);

    static double calcReflectionCoefficients(dcomplex refractive_index, double theta);

    static int factorial(int n);

    int IDUM;

private:
    // statistics and SED

    void CholDec();

    void SolveEqSys();

    double ChiSq();

    Matrix2D stat_data;
    uint nr_ofSeq;

    // random number generator
    double RM;
    int IY, IFF;
    int IR[98];

    int Mr;
    int IA;
    int IC;

    // Chi squared

    double * b;
    double * a;
    double * r;
    double * d;
    double * y;
    double * sigma;
    Matrix2D A;
    Matrix2D C;
    uint N;
    uint M;

    ullong kiss_x;
    ullong kiss_y;
    ullong kiss_z;
    ullong kiss_c;
};

#endif /* CMATH_FUNCTIONS_H */
