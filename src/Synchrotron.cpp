#include "Synchrotron.h"

// Split the function domain into three intervals:
// (0, 0.001), [0.001, 12), and (12, infinity)
double CSynchrotron::Gamma(double x)
{
    //Actually Gamma(0) is infinity. However, this case shouldn't happen anyway for any 2<p<3.
    if (x <= 0.0)
    {
        cout << "ERROR: Gamma(0)                                            \n";
        return 0;
    }

    // First interval: (0, 0.001)
    // In this range, 1/Gamma(x) = x + gamma x^2 with error less than 6e-7.

    if (x < 0.001)
        return 1.0 / (x * (1.0 + Euler_gamma * x));

    // Second interval: [0.001, 12)
    // Algorithm directly with reduction identities to reduce other arguments to this interval.

    if (x < 12.0)
    {
        double y = x;
        int n = 0;
        bool arg_was_less_than_one = (y < 1.0);

        if (arg_was_less_than_one)
            y += 1.0;
        else
        {
            n = int(floor(y)) - 1;
            y -= n;
        }

        double num = 0.0;
        double den = 1.0;

        double z = y - 1;
        for (uint i = 0; i < 8; i++)
        {
            num = (num + p[i]) * z;
            den = den * z + q[i];
        }
        double result = num / den + 1.0;

        if (arg_was_less_than_one)
            result /= (y - 1.0);
        else
        {
            // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
            for (uint i = 0; i < n; i++)
                result *= y++;
        }

        return result;
    }

    // Third interval: [12, infinity)
    // Gamma(171.624) is larger than double precession. However, this shouldn't happen anyway for any realistic value of p.
    if (x > 171.624)
    {
        cout << "ERROR: Gamma(171.624)                                      \n";
        return 0;
    }

    return exp(LogGamma(x));
}

double CSynchrotron::LogGamma(double x)
{
    if (x <= 0.0)
        return 0;

    if (x < 12.0)
        return log(fabs(Gamma(x)));

    double z = 1.0 / (x * x);
    double sum = c[7];

    for (int i = 6; i >= 0; i--)
    {
        sum *= z;
        sum += c[i];
    }

    return (x - 0.5) * log(x) - x + halfLogTwoPi + sum / x;
}

syn_param CSynchrotron::get_Thermal_Parameter(double n_e, double T_e, double l, double B, double theta)
{
    syn_param res;

    if(n_e==0)
        return res;

    double Theta=syn_kB*T_e/(syn_me*syn_c*syn_c);

    n_e*=1e-6;
    B*=1e4; //B in mu Gauss

    double sin_theta=sin(theta);
    double cos_theta=cos(theta);

    double nu=con_c/l;
    double nu_c =  syn_e*B/ (PIx2 * syn_me * syn_c);
    double nu_s = (2./9.)*nu_c*sin_theta*Theta*Theta;

    double x = nu/nu_s;


    //cout << "Z: " << nu_c << "\t" << nu_s << "\t" << x << "\n";



    double prefactor = (n_e*syn_e*syn_e* nu_c)/syn_c;
    double planck = planck_hz(nu, Theta);

    //cout << "Planck 2: " << planck << "\n";

//    cout << "P1: " << prefactor << "\t" << sin_theta*getI_I(x) << "\n";
//    cout << "P2: " << prefactor << "\t" << sin_theta*getI_Q(x,Theta) << "\n";
//    cout << "P3: " << prefactor << "\t" << getI_V(x,theta,Theta) << "\n";

    res.j_I=prefactor*sin_theta*getI_I_th(x);
    res.j_Q=prefactor*sin_theta*getI_Q_th(x,Theta);
    res.j_V=prefactor*getI_V_th(x,theta,Theta);

    res.alpha_I=res.j_I/planck;
    res.alpha_Q=res.j_Q/planck;
    res.alpha_V=res.j_V/planck;

    double omega0 = syn_e*B / (syn_me*syn_c);

    double wp2 = PIx4 * n_e*syn_e*syn_e / syn_me;

  // argument for function corr. funcitons
  x = Theta * sqrt(sqrt(2.) * sin_theta*(1.e3*omega0 / (2. * PI * nu)));
  
  //cout << "Z1   omega0: " << omega0 << "\t wp2: " << wp2 << "\t x: " << x <<  "\n\n";
  
  double bessel_0=BesselK(1,1./Theta);
  double bessel_1=BesselK(1,1./Theta);
  double bessel_2=BesselK(2,1./Theta);
  
  //(bessel_1/ bessel_2 + 6. * Theta)*
  res.kappa_Q = getK_Q_th(x)*PIx2 * nu /(2. * syn_c) * wp2 * omega0*omega0 
                    / pow(PIx2 * nu, 4.);//*sin_theta*sin_theta;

//  l=0.01*l;
//  double test=n_e*syn_e*syn_e*B/(PI*syn_me*syn_me*syn_c*syn_c*syn_c*syn_c)*l*l;
  //* (bessel_0 - getK_V_th(x))/ bessel_2;
  res.kappa_V = PIx2 * nu / syn_c * wp2 * omega0 / pow((2. * PI * nu), 3.);//* cos_theta;

//  double test=n_e*pow(syn_e,4)*pow(B*/(PI*syn_me*syn_me*syn_c*syn_c*syn_c*syn_c)*l*l;  
//  res.kappa_V=res.kappa_Q;               

  if(res.kappa_Q!=res.kappa_Q)
      int tt=0;
                   
  if(res.kappa_V!=res.kappa_V)
      int tt=0;

    //converting back into SI;
    res.scale();

    return res;
}

syn_param CSynchrotron::get_Power_Law_Parameter(double n_e, double l, double B, double theta, double g_min, double g_max, double p)
{
    syn_param res;

    n_e*=1e-6;
    B*=1e4;

    double nu=con_c/l;
    double nu_c =  syn_e*B/ (PIx2 * syn_me * syn_c);
    double sin_theta=sin(theta);
    double cos_theta=cos(theta);
    double tan_theta=tan(theta);

    double gammas=Gamma_I_p(g_min, g_max, p);

    res.j_I=(n_e*syn_e*syn_e*nu_c)/syn_c*pow(3.0,0.5*p)*(p-1.0)*sin_theta*gammas*pow(nu/(nu_c*sin_theta),-(p-1.0)/2.0);
    res.j_Q=res.j_I*getI_Q_p(p);
    res.j_V=res.j_I/sqrt(nu/(3.*nu_c*sin_theta))*getI_V_p(p,tan_theta);

    gammas = Gamma_A_p(g_min,g_max,p);

    res.alpha_I= (n_e*pow(syn_e,2.))/(nu*syn_me*syn_c)*
                pow(3.0, (p+1.0)/2.0)*(p-1.0)*
                gammas*pow(nu/(nu_c*sin_theta),-(p+2.)/2.);

    res.alpha_Q=-res.alpha_I*getA_Q_p(p);

    res.alpha_V=res.alpha_I/sqrt(nu/(nu_c*sin_theta))*getA_V_p(p,sin_theta,cos_theta);

    double kappa_perp=n_e*syn_e*syn_e*(p-1)/(syn_me*syn_c*nu_c*sin_theta)/
    (pow(g_min,1-p)-pow(g_max,1-p));

    res.kappa_Q=0.0;//-kappa_perp*pow(nu_c*sin_theta/nu,3)*pow(g_min,2.0-p)*(1-pow(2*nu_c/(3*nu),0.5*p-1))/(0.5*p-1);
    res.kappa_V=0.0;//2.0*(p+2)/(p+1)*kappa_perp*pow(nu_c*sin_theta/nu,2)*pow(g_min,-(p+1)) * log(g_min)/tan_theta;

    //converting back into SI;
    res.scale();

  return res;
}
