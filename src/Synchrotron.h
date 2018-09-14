#include "typedefs.h"
#include "chelper.h"
#include <cmath>

#ifndef SYNCHROTRON_H
#define SYNCHROTRON_H

//redefine some natural constants in cgs
#define syn_me   9.1093826e-28  //electron mass      [g]            9.10938188e-31
#define syn_e    4.80320680e-10 //electron charge    [Statcoulomb]  1.602176634e-19
#define syn_h    6.6260693e-27  //Planck constant    [erg/s]        6.626176e-34
#define syn_kB   1.380662e-16   //Boltzmann constant [erg/K]
#define syn_c    2.99792458e10  //speed of light     [m/s]          299792458.0

//conversion factor 1/cm -> 1/m
#define syn_SI   100.0


//container class for the parameters of sync. RT
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


    syn_param(double _j_I, double _j_Q, double _j_V, double _alpha_I, double _alpha_Q, double _alpha_V, double _kappa_Q, double _kappa_V)
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

    syn_param operator+(const syn_param & rhs)
    {
        return syn_param(j_I + rhs.j_I, j_Q + rhs.j_Q, j_V + rhs.j_V,
                alpha_I + rhs.alpha_I, alpha_Q + rhs.alpha_Q, alpha_V + rhs.alpha_V,
                kappa_Q + rhs.kappa_Q, kappa_V+ rhs.kappa_V);
    }

    void scale()
    {
        j_I *= syn_SI;
        j_Q *= syn_SI;
        j_V *= syn_SI;

        alpha_I *= syn_SI;
        alpha_Q *= syn_SI;
        alpha_V *= syn_SI;

        kappa_Q *= syn_SI;
        kappa_V *= syn_SI;
    }

    //generate matrix for synchrotron RT
    Matrix2D getSyncMatrix()
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
        cont_matrix.setValue(3, 2,-kappa_Q);

        cont_matrix.setValue(1, 2, kappa_V);
        cont_matrix.setValue(2, 1,-kappa_V);

        //cont_matrix.printMatrix();
        return cont_matrix;
    }

    double j_I;
    double j_Q;
    double j_V;

    double alpha_I;
    double alpha_Q;
    double alpha_V;

    double kappa_Q;
    double kappa_V;
};


//class for the physics of sync. RT
class CSynchrotron
{
public:

    CSynchrotron()
    {
        initGamma();
        initBesselK();
    };

    ~CSynchrotron()
    {
        if(p != 0)
        {
            delete[] p;
            p = 0;
        }

        if(q != 0)
        {
            delete[] q;
            q = 0;
        }

        if(c != 0)
        {
            delete[] c;
            c = 0;
        }
    };


    syn_param get_Thermal_Parameter(double n_e, double T_e, double l, double B, double theta);
    syn_param get_Power_Law_Parameter(double n_e, double l, double B, double theta, double g_min, double g_max, double p);

    double BesselK(uint n, double x)
    {
        double res = 0;

        switch(n)
        {
            case 0:
                res = BesselK_0.getValue(x);
                break;

            case 1:
                res = BesselK_1.getValue(x);
                break;

            case 2:
                res = BesselK_2.getValue(x);
                break;

            default:
                //cout << "ERROR: BesselK_"<<n<<" is not defined                           \n";
                return 0;
        }

        return pow(10.0, res);
    }

private:
    // Gamma numerator coefficients for approximation over the interval (1,2)
    double * p;

    // Gamma denominator coefficients for approximation over the interval (1,2)
    double * q;

    // Asymptotic series LogGamma, Abramowitz and Stegun 6.1.41
    double * c;

    double Euler_gamma;
    double halfLogTwoPi;

    spline BesselK_0;
    spline BesselK_1;
    spline BesselK_2;

    void initGamma()
    {
        p = new double[8];
        q = new double[8];
        c = new double[8];

        p[0] = -1.71618513886549492533811E+0;
        p[1] = 2.47656508055759199108314E+1;
        p[2] = -3.79804256470945635097577E+2;
        p[3] = 6.29331155312818442661052E+2;
        p[4] = 8.66966202790413211295064E+2;
        p[5] = -3.14512729688483675254357E+4;
        p[6] = -3.61444134186911729807069E+4;
        p[7] = 6.64561438202405440627855E+4;

        q[0] = -3.08402300119738975254353E+1;
        q[1] = 3.15350626979604161529144E+2;
        q[2] = -1.01515636749021914166146E+3;
        q[3] = -3.10777167157231109440444E+3;
        q[4] = 2.25381184209801510330112E+4;
        q[5] = 4.75584627752788110767815E+3;
        q[6] = -1.34659959864969306392456E+5;
        q[7] = -1.15132259675553483497211E+5;

        c[0] = 1.0 / 12.0;
        c[1] = -1.0 / 360.0;
        c[2] = 1.0 / 1260.0;
        c[3] = -1.0 / 1680.0;
        c[4] = 1.0 / 1188.0;
        c[5] = -691.0 / 360360.0;
        c[6] = 1.0 / 156.0;
        c[7] = -3617.0 / 122400.0;

        Euler_gamma = 0.577215664901532860606512090; // Euler's gamma constant
        halfLogTwoPi = 0.91893853320467274178032973640562;
    }

    void initBesselK()
    {
        uint size = 68;
        BesselK_0.resize(size);
        BesselK_1.resize(size);
        BesselK_2.resize(size);

        BesselK_0.setValue(0, 1.0e-9, 1.318880988);
        BesselK_0.setValue(1, 1.0e-8, 1.268030366);
        BesselK_0.setValue(2, 1.0e-7, 1.210426268);
        BesselK_0.setValue(3, 1.0e-6, 1.144);
        BesselK_0.setValue(4, 2.0e-6, 1.12183);
        BesselK_0.setValue(5, 3.0e-6, 1.10832);
        BesselK_0.setValue(6, 4.0e-6, 1.09848);
        BesselK_0.setValue(7, 5.0e-6, 1.09068);
        BesselK_0.setValue(8, 6.0e-6, 1.08421);
        BesselK_0.setValue(9, 7.0e-6, 1.07866);
        BesselK_0.setValue(10, 8.0e-6, 1.07379);
        BesselK_0.setValue(11, 9.0e-6, 1.06945);
        BesselK_0.setValue(12, 0.00001, 1.06554);
        BesselK_0.setValue(13, 0.00002, 1.03885);
        BesselK_0.setValue(14, 0.00003, 1.02244);
        BesselK_0.setValue(15, 0.00004, 1.01041);
        BesselK_0.setValue(16, 0.00005, 1.00084);
        BesselK_0.setValue(17, 0.00006, 0.992867);
        BesselK_0.setValue(18, 0.00007, 0.986008);
        BesselK_0.setValue(19, 0.00008, 0.979977);
        BesselK_0.setValue(20, 0.00009, 0.974587);
        BesselK_0.setValue(21, 0.0001, 0.969708);
        BesselK_0.setValue(22, 0.0002, 0.936168);
        BesselK_0.setValue(23, 0.0003, 0.915276);
        BesselK_0.setValue(24, 0.0004, 0.899819);
        BesselK_0.setValue(25, 0.0005, 0.887439);
        BesselK_0.setValue(26, 0.0006, 0.877055);
        BesselK_0.setValue(27, 0.0007, 0.868078);
        BesselK_0.setValue(28, 0.0008, 0.860148);
        BesselK_0.setValue(29, 0.0009, 0.853032);
        BesselK_0.setValue(30, 0.001, 0.846565);
        BesselK_0.setValue(31, 0.002, 0.801441);
        BesselK_0.setValue(32, 0.003, 0.772695);
        BesselK_0.setValue(33, 0.004, 0.75108);
        BesselK_0.setValue(34, 0.005, 0.733541);
        BesselK_0.setValue(35, 0.006, 0.718666);
        BesselK_0.setValue(36, 0.007, 0.70568);
        BesselK_0.setValue(37, 0.008, 0.694108);
        BesselK_0.setValue(38, 0.009, 0.68364);
        BesselK_0.setValue(39, 0.01, 0.674057);
        BesselK_0.setValue(40, 0.02, 0.605139);
        BesselK_0.setValue(41, 0.03, 0.559132);
        BesselK_0.setValue(42, 0.04, 0.523297);
        BesselK_0.setValue(43, 0.05, 0.493351);
        BesselK_0.setValue(44, 0.06, 0.467294);
        BesselK_0.setValue(45, 0.07, 0.444016);
        BesselK_0.setValue(46, 0.08, 0.422834);
        BesselK_0.setValue(47, 0.09, 0.403295);
        BesselK_0.setValue(48, 0.1, 0.385082);
        BesselK_0.setValue(49, 0.2, 0.243709);
        BesselK_0.setValue(50, 0.3, 0.1375);
        BesselK_0.setValue(51, 0.4, 0.0470914);
        BesselK_0.setValue(52, 0.5, -0.0341311);
        BesselK_0.setValue(53, 0.6, -0.109287);
        BesselK_0.setValue(54, 0.7, -0.180114);
        BesselK_0.setValue(55, 0.8, -0.247685);
        BesselK_0.setValue(56, 0.9, -0.312712);
        BesselK_0.setValue(57, 1.0, -0.375693);
        BesselK_0.setValue(58, 2.0, -0.9435);
        BesselK_0.setValue(59, 3.0, -1.45918);
        BesselK_0.setValue(60, 4.0, -1.95235);
        BesselK_0.setValue(61, 5.0, -2.43284);
        BesselK_0.setValue(62, 6.0, -2.90518);
        BesselK_0.setValue(63, 7.0, -3.37182);
        BesselK_0.setValue(64, 8.0, -3.83425);
        BesselK_0.setValue(65, 9.0, -4.29344);
        BesselK_0.setValue(66, 10.0, -4.75007);
        BesselK_0.setValue(67, 20.0, -9.24099);

        BesselK_1.setValue(0, 1.0e-9, 9.000000000);
        BesselK_1.setValue(1, 1.0e-8, 8.000000000);
        BesselK_1.setValue(2, 1.0e-7, 7.000000000);
        BesselK_1.setValue(3, 1.0e-6, 6.);
        BesselK_1.setValue(4, 2.0e-6, 5.69897);
        BesselK_1.setValue(5, 3.0e-6, 5.52288);
        BesselK_1.setValue(6, 4.0e-6, 5.39794);
        BesselK_1.setValue(7, 5.0e-6, 5.30103);
        BesselK_1.setValue(8, 6.0e-6, 5.22185);
        BesselK_1.setValue(9, 7.0e-6, 5.1549);
        BesselK_1.setValue(10, 8.0e-6, 5.09691);
        BesselK_1.setValue(11, 9.0e-6, 5.04576);
        BesselK_1.setValue(12, 0.00001, 5.);
        BesselK_1.setValue(13, 0.00002, 4.69897);
        BesselK_1.setValue(14, 0.00003, 4.52288);
        BesselK_1.setValue(15, 0.00004, 4.39794);
        BesselK_1.setValue(16, 0.00005, 4.30103);
        BesselK_1.setValue(17, 0.00006, 4.22185);
        BesselK_1.setValue(18, 0.00007, 4.1549);
        BesselK_1.setValue(19, 0.00008, 4.09691);
        BesselK_1.setValue(20, 0.00009, 4.04576);
        BesselK_1.setValue(21, 0.0001, 4.);
        BesselK_1.setValue(22, 0.0002, 3.69897);
        BesselK_1.setValue(23, 0.0003, 3.52288);
        BesselK_1.setValue(24, 0.0004, 3.39794);
        BesselK_1.setValue(25, 0.0005, 3.30103);
        BesselK_1.setValue(26, 0.0006, 3.22185);
        BesselK_1.setValue(27, 0.0007, 3.1549);
        BesselK_1.setValue(28, 0.0008, 3.09691);
        BesselK_1.setValue(29, 0.0009, 3.04576);
        BesselK_1.setValue(30, 0.001, 3.);
        BesselK_1.setValue(31, 0.002, 2.69896);
        BesselK_1.setValue(32, 0.003, 2.52287);
        BesselK_1.setValue(33, 0.004, 2.39792);
        BesselK_1.setValue(34, 0.005, 2.301);
        BesselK_1.setValue(35, 0.006, 2.2218);
        BesselK_1.setValue(36, 0.007, 2.15484);
        BesselK_1.setValue(37, 0.008, 2.09683);
        BesselK_1.setValue(38, 0.009, 2.04566);
        BesselK_1.setValue(39, 0.01, 1.99989);
        BesselK_1.setValue(40, 0.02, 1.69858);
        BesselK_1.setValue(41, 0.03, 1.52207);
        BesselK_1.setValue(42, 0.04, 1.39661);
        BesselK_1.setValue(43, 0.05, 1.29906);
        BesselK_1.setValue(44, 0.06, 1.21916);
        BesselK_1.setValue(45, 0.07, 1.1514);
        BesselK_1.setValue(46, 0.08, 1.09252);
        BesselK_1.setValue(47, 0.09, 1.0404);
        BesselK_1.setValue(48, 0.1, 0.993606);
        BesselK_1.setValue(49, 0.2, 0.679062);
        BesselK_1.setValue(50, 0.3, 0.485152);
        BesselK_1.setValue(51, 0.4, 0.339323);
        BesselK_1.setValue(52, 0.5, 0.219176);
        BesselK_1.setValue(53, 0.6, 0.114889);
        BesselK_1.setValue(54, 0.7, 0.0213066);
        BesselK_1.setValue(55, 0.8, -0.0646028);
        BesselK_1.setValue(56, 0.9, -0.144763);
        BesselK_1.setValue(57, 1.0, -0.22047);
        BesselK_1.setValue(58, 2.0, -0.854288);
        BesselK_1.setValue(59, 3.0, -1.39624);
        BesselK_1.setValue(60, 4.0, -1.90366);
        BesselK_1.setValue(61, 5.0, -2.39312);
        BesselK_1.setValue(62, 6.0, -2.87163);
        BesselK_1.setValue(63, 7.0, -3.34277);
        BesselK_1.setValue(64, 8.0, -3.80864);
        BesselK_1.setValue(65, 9.0, -4.27054);
        BesselK_1.setValue(66, 10.0, -4.72935);
        BesselK_1.setValue(67, 20.0, -9.2304);

        BesselK_2.setValue(0, 1.0e-9, 18.30103000);
        BesselK_2.setValue(1, 1.0e-8, 16.30103000);
        BesselK_2.setValue(2, 1.0e-7, 14.30103000);
        BesselK_2.setValue(3, 1.0e-6, 12.301);
        BesselK_2.setValue(4, 2.0e-6, 11.699);
        BesselK_2.setValue(5, 3.0e-6, 11.3468);
        BesselK_2.setValue(6, 4.0e-6, 11.0969);
        BesselK_2.setValue(7, 5.0e-6, 10.9031);
        BesselK_2.setValue(8, 6.0e-6, 10.7447);
        BesselK_2.setValue(9, 7.0e-6, 10.6108);
        BesselK_2.setValue(10, 8.0e-6, 10.4949);
        BesselK_2.setValue(11, 9.0e-6, 10.3925);
        BesselK_2.setValue(12, 0.00001, 10.301);
        BesselK_2.setValue(13, 0.00002, 9.69897);
        BesselK_2.setValue(14, 0.00003, 9.34679);
        BesselK_2.setValue(15, 0.00004, 9.09691);
        BesselK_2.setValue(16, 0.00005, 8.90309);
        BesselK_2.setValue(17, 0.00006, 8.74473);
        BesselK_2.setValue(18, 0.00007, 8.61083);
        BesselK_2.setValue(19, 0.00008, 8.49485);
        BesselK_2.setValue(20, 0.00009, 8.39254);
        BesselK_2.setValue(21, 0.0001, 8.30103);
        BesselK_2.setValue(22, 0.0002, 7.69897);
        BesselK_2.setValue(23, 0.0003, 7.34679);
        BesselK_2.setValue(24, 0.0004, 7.09691);
        BesselK_2.setValue(25, 0.0005, 6.90309);
        BesselK_2.setValue(26, 0.0006, 6.74473);
        BesselK_2.setValue(27, 0.0007, 6.61083);
        BesselK_2.setValue(28, 0.0008, 6.49485);
        BesselK_2.setValue(29, 0.0009, 6.39254);
        BesselK_2.setValue(30, 0.001, 6.30103);
        BesselK_2.setValue(31, 0.002, 5.69897);
        BesselK_2.setValue(32, 0.003, 5.34679);
        BesselK_2.setValue(33, 0.004, 5.09691);
        BesselK_2.setValue(34, 0.005, 4.90309);
        BesselK_2.setValue(35, 0.006, 4.74472);
        BesselK_2.setValue(36, 0.007, 4.61083);
        BesselK_2.setValue(37, 0.008, 4.49484);
        BesselK_2.setValue(38, 0.009, 4.39254);
        BesselK_2.setValue(39, 0.01, 4.30102);
        BesselK_2.setValue(40, 0.02, 3.69893);
        BesselK_2.setValue(41, 0.03, 3.34669);
        BesselK_2.setValue(42, 0.04, 3.09674);
        BesselK_2.setValue(43, 0.05, 2.90282);
        BesselK_2.setValue(44, 0.06, 2.74434);
        BesselK_2.setValue(45, 0.07, 2.6103);
        BesselK_2.setValue(46, 0.08, 2.49416);
        BesselK_2.setValue(47, 0.09, 2.39167);
        BesselK_2.setValue(48, 0.1, 2.29995);
        BesselK_2.setValue(49, 0.2, 1.69471);
        BesselK_2.setValue(50, 0.3, 1.33737);
        BesselK_2.setValue(51, 0.4, 1.08049);
        BesselK_2.setValue(52, 0.5, 0.877958);
        BesselK_2.setValue(53, 0.6, 0.709296);
        BesselK_2.setValue(54, 0.7, 0.563639);
        BesselK_2.setValue(55, 0.8, 0.434537);
        BesselK_2.setValue(56, 0.9, 0.31786);
        BesselK_2.setValue(57, 1.0, 0.21081);
        BesselK_2.setValue(58, 2.0, -0.595577);
        BesselK_2.setValue(59, 3.0, -1.21105);
        BesselK_2.setValue(60, 4.0, -1.75942);
        BesselK_2.setValue(61, 5.0, -2.27499);
        BesselK_2.setValue(62, 6.0, -2.77161);
        BesselK_2.setValue(63, 7.0, -3.25605);
        BesselK_2.setValue(64, 8.0, -3.73209);
        BesselK_2.setValue(65, 9.0, -4.20204);
        BesselK_2.setValue(66, 10.0, -4.66736);
        BesselK_2.setValue(67, 20.0, -9.19863);

        BesselK_0.createSpline();
        BesselK_1.createSpline();
        BesselK_2.createSpline();
    }

    double LogGamma(double x);
    double Gamma(double x);

    double Gamma_I_p(double g_min, double g_max, double p)
    {
        return Gamma((3.0 * p - 1.0) / 12.0) * Gamma((3.0 * p + 19.) / 12.) / (2.0 * (p + 1.)*(pow(g_min, 1.0 - p) - pow(g_max, 1. - p)));
    }

    double Gamma_A_p(double g_min, double g_max, double p)
    {
        return Gamma((3. * p + 2.) / 12.) * Gamma((3. * p + 22.) / 12.) / (4. * (pow(g_min, 1. - p) - pow(g_max, 1. - p)));
    }

    double getI_Q_p(double p)
    {
        return(-(p + 1.0) / (p + 7.0 / 3.0));
    }

    double getI_V_p(double p, double tan_theta)
    {
        return(171.0 / 250.0)*pow(p, 49.0 / 100.0) / tan_theta;
    }

    double getA_Q_p(double p)
    {
        return(-pow((17. / 500.) * p - 43. / 1250., 43. / 500.));
    }

    double getA_V_p(double p, double sin_theta, double cos_theta)
    {
        double sign = 1;

        if(cos_theta != 0)
            sign = cos_theta / fabs(cos_theta);

        return sign * (pow((71. / 100.) * p + 22. / 625., 197. / 500.)) * pow((31. / 10.) * pow(sin_theta, -48. / 25) - 31. / 10., 64. / 125.);
    }

    double getK_Q_th(double x)
    {
        // Correction factor according to  Dexter (2016)
        double corr = (.011 * exp(-x / 47.2) - pow(2., (-1.0 / 3.0)) / pow(3.0, (23. / 6.0))
                * PI * 1.0e4 * pow(x, (-8.0 / 3.0)))
                * (0.5 + 0.5 * tanh((log(x) - log(120.0)) / 0.1));

        return 2.011 * exp(-pow(x, 1.035) / 4.7) - cos(x / 2.)
                * exp(-pow(x, 1.2) / 2.73) - .011 * exp(-x / 47.2) + corr;
    }

    double getK_V_th(double x)
    {
        return 0.43793091 * log(1.0 + 0.00185777 * pow(x, 1.50316886));
    }

    double getI_I_th(double x)
    {
        double term1 = sqrt(2.0) * PI / 27.0;
        double term2 = sqrt(x) + pow(2.0, 11.0 / 12.0) * pow(x, 1. / 6.);
        term2 *= term2;
        double term3 = exp(-pow(x, 1. / 3.));
        return term1 * term2*term3;
    }

    double getI_Q_th(double x, double Theta)
    {
        double term1 = -sqrt(2.0) * PI / 27.0;
        double term2 = (7.0 * pow(Theta, 24.0 / 25.0) + 35.0) / (10.0 * pow(Theta, 24.0 / 25.0) + 75.0);
        double term3 = sqrt(x) + term2 * pow(2., 11. / 12.) * pow(x, 1. / 6.);

        term3 *= term3;

        return term1 * term3 * exp(-pow(x, 1.0 / 3.0));
    }

    double getI_V_th(double x, double theta, double Theta)
    {
        double term1 = (37.0 - 87.0 * sin(theta - 28.0 / 25.0)) / (100.0 * (Theta + 1.0));
        double term2 = pow(1.0 + (pow(Theta, 3.0 / 5.0) / 25.0 + 7.0 / 10.0) * pow(x, 9.0 / 25.0), 5.0 / 3.0);
        return term1 * term2 * exp(-pow(x, 1. / 3.));
    }

    double planck_hz(double nu, double Theta)
    {
        double term1 = (2. * syn_h * pow(nu, 3.)) / pow(syn_c, 2.);
        double term2 = (exp(syn_h * nu / (Theta * syn_me * pow(syn_c, 2.))) - 1.);

        return term1 / term2;
    }
};
#endif
