#include "GrainCharge.hpp"
#include "DustComponent.hpp"

CGrainCharge::CGrainCharge()
{
    mix_ID = 0;
    //W_carbon_eV = 4.4;
    //W_sil_eV    = 8.0;
    W_eV = 8.0;
    
    //rho_sil = 2200;
    //rho_carb = 3300;
    
    dust = 0;
    
    rho = 2000;
    
    // 1 D  = 3.33564e−30 C m
    //beta_sil=sqrt(0.5*2.44*2.44 + 0.5*1.9*1.9)*3.33564e-30;//*3.33564e−30; //C m
    //beta_carb=0.38*3.33564e-30;//*3.33564e−30; //C m
    
    beta=sqrt(0.5*2.44*2.44 + 0.5*1.9*1.9)*3.33564e-30;//*3.33564e−30; //C m
    
    m0 = 24*con_m_p; //sil
    //m0 = 12*con_m_p; //carb
    
    // Tensile strengths
    Smax          = 5.0e9;  // Pa (material tensile strength)
    
    eV_to_J = con_e;
    J_to_eV = 1.0 / eV_to_J;
    m_to_AA = 1.0e10;
    
    //carbon = false;
    
    Tgas   = 1e3;         // K
    Tdust  = 15;            // K
    
    n_el  = 1.0e4;         // m^-3
    n_ion = 1.0e5;         // m^-3
    n_gas = 1.0e6;          // m^-3
    
    n_neu = n_gas - n_ion;
    
    stick_e    = 0.5;           // electron sticking (simplified)
    stick_ion  = 1.0;           // ion sticking
    
    Zgas = 0;
    mu_mol=0;
    
    //PhysRevA.111.012801
    alpha_gas = 3.75e-40; //H2
    //alpha_gas = 8.8e-41; //H2
            
    //alpha_gas = 7.4e-41 //neutral H
    //alpha_gas = 2.3e-41 //neutral He

    // Radiation field (demo): flat u_lambda between 613.6 eV
    E_low_eV      = 6.0;
    E_up_eV       = 13.6;

    Nla = 100;

    // Recurrence window controls
    tail_tol_ratio = 1e-6;
    Z_span_default = 60;
    
    arr_a_eff = 0;
    arr_dnda = 0;
    
    arr_a_eff_large = 0;    
    arr_dnda_large = 0;    
    
    tot_Zmin=10000;
    tot_Zmax=-10000;

    a_min=0;
    a_max=1;
    Na=0;
    a0=1;
    a_sigma=1;
    
    lambda_min = 0;
    lambda_max = 0;
}

CGrainCharge::~CGrainCharge()
{
    clear();
}

void CGrainCharge::set_density(double _n_gas, double _n_el, double _n_ion)
{
    n_el  = _n_el;  // m^-3
    n_ion = _n_ion; // m^-3
    n_gas = _n_gas; // m^-3
    
    n_neu = abs(n_gas - n_ion);
}

void CGrainCharge::set_temp(double _T_gas, double _T_el, double _T_dust)
{
    Tgas = _T_gas;
    Tdust = _T_dust;
}

void CGrainCharge::set_wavelengths(const dlist & _wavelengths_list)
{
    wavelength_list = _wavelengths_list;
    
    lambda_min = wavelength_list[0];
    lambda_max = wavelength_list[wavelength_list.size()-1];
}

void CGrainCharge::set_r(double _Zgas, double _mu)
{
    Zgas = _Zgas;
    mu_mol = _mu;
}


// ======================= Currents struct =======================
Currents CGrainCharge::total_currents(const CGridBasic * grid, const cell_basic & cell, int Z, double a_eff)
{
    Currents J;
    J.J_up   = J_ion_Hp(Z, a_eff) + J_photoelectric(grid, cell,Z, a_eff);
    J.J_down = J_electron(Z, a_eff) + J_photodetachment(grid, cell, Z, a_eff);
    return J;
}

// ======================= Hard Z-bounds =======================
inline int CGrainCharge::compute_Zmin_auto(double a_m)
{
    // Smallest Z with EA(Z,a) > 0
    int Z = -1;

    while ( Z > -10000 )
    {
        double EA = EA_eV( a_m, Z);
        if ( EA <= 0.0 )
        {
            return Z + 1;
        }
        Z = Z - 1;
    }

    return -9999; // conservative default value
}

int CGrainCharge::compute_Zmax_coulomb( double a_m, double Smax_Pa )
{
    // Stress p = Q^2 / ( 32 π^2 ε0 a^4 ); set p = Smax.
    // Z_max = Q_max / e = ( 4 √2 π a^2 / e ) √( ε0 Smax )
    double pref   = 4.0 * sqrt(2) * PI * a_m * a_m / con_e;
    double inside = con_epsilon_0 * Smax_Pa;
    if ( inside < 0.0 ) inside = 0.0;
    double Zmax_real = pref * sqrt( inside );

    int Zmax_int = int( floor( Zmax_real ) );
    return Zmax_int;
}

// ======================= Helper functions for recurrence =======================
dlist CGrainCharge::compute_f_window(const CGridBasic * grid, const cell_basic & cell, int Zmin, int Zmax, int Zhard_min, int Zhard_max, double a_eff)
{
    int N = ( Zmax - Zmin + 1 );
    dlist logf(N, -1e300);

    logf[0] = 0.0;

    for (int k = 0; k < N - 1; k++)
    {
        int Z      = Zmin + k;
        Currents Jz   = total_currents(grid, cell, Z, a_eff);
        Currents Jzp1 = total_currents(grid, cell, Z + 1, a_eff);

        double num = Jz.J_up;
        double den = Jzp1.J_down;

        if (Z >= Zhard_max) num = 0.0;
        if (Z + 1 <= Zhard_min) den = 0.0;

        double r;

        if (den > 0.0)
        {
            r = num / den;
        }
        else
        {
            r = 0.0;
        }

        if (r > 0.0)
        {
            logf[k + 1] = logf[k] + log(r);
        }
        else
        {
            logf[k + 1] = -1e300;
        }
    }

    // normalize safely
    double mval = -1e300;
    for ( int k = 0; k < N; k++)
    {
        if (logf[k] > mval)
        {
            mval = logf[k];
        }
    }

    dlist f(N, 0.0 );
    double sum = 0.0;

    for ( int k = 0; k < N; ++k )
    {
        double val = exp( logf[k] - mval );
        f[k] = val;

        sum += val;
    }

    if(sum > 0.0)
    {
        for (int k = 0; k < N; ++k )
        {
            f[k] /= sum;
        }
    }

    return f;
}

bool CGrainCharge::tails_are_small(const dlist &f)
{
    int N = f.size();

    if (N < 3)
        return false;

    double left  = f[0] + f[1];
    double right = f[N - 1] + f[N - 2];

    if ( left  < tail_tol_ratio && right < tail_tol_ratio )
        return true;

    return false;
}

int CGrainCharge::find_Zeq(const CGridBasic * grid, const cell_basic & cell, int Zmin, int Zmax, double a_eff)
{
    double best = numeric_limits<double>::infinity();
    int Zeq = 0;

    for (int Z = Zmin; Z <= Zmax; Z++)
    {
        Currents J = total_currents(grid, cell, Z, a_eff);

        double diff = fabs( J.J_up - J.J_down );

        if (diff < best)
        {
            best = diff;
            Zeq  = Z;
        }
    }

    return Zeq;
}

// Evaluate Gaussian at integer Z given parameters.
double CGrainCharge::gaussian_value(double mu, double sigma, int Z)
{
    if ( sigma <= 0.0 )
        return 0.0;

    double dz   = double( Z ) - mu;
    double arg  = -0.5 * ( dz * dz ) / ( sigma * sigma );
    double A = 1.0 / ( sqrt(PIx2sq) * sigma );

    const double val  = A * exp( arg );

    return val;
}

// Compute sum of squares error between f and model for given parameters (diagnostics).
double CGrainCharge::gaussian_sse( const ilist &Z, const dlist &f,
                            double mu, double sigma)
{
    const size_t N = Z.size();
    double sse = 0.0;

    size_t i = 0;
    while ( i < N )
    {
        const double g  = gaussian_value(mu, sigma, Z[i]);
        const double r  = g - f[ i ];
        sse += r * r;
        i = i + 1;
    }

    return sse;
}

void CGrainCharge::fit_gaussian_pdf( const ilist   &Z,
                                     const dlist &f , 
                                     double & mu, double & sigma)
{
    mu    = 0.0;
    sigma = 0.0;

    int N = Z.size();

    if ( N == 0 )
        return;

    if ( f.size() != N )
        return;

    dlist p(N, 0.0);
    double sumf = 0.0;

    int i = 0;

    while (i < N)
    {
        if ( f[ i ] > 0.0 )
        {
            p[ i ] = f[ i ];
            sumf  += f[ i ];
        }
        else
        {
            p[ i ] = 0.0;
        }
        i = i + 1;
    }

    if ( sumf <= 0.0 )
    {
        return;
    }

    i = 0;

    while ( i < N )
    {
        p[ i ] = p[ i ] / sumf;
        i = i + 1;
    }

    // Moment estimates for mu and sigma.
    double var  = 0.0;

    i = 0;
    while ( i < N )
    {
        mu += (Z[i]) * p[ i ];
        i = i + 1;
    }

    i = 0;
    while ( i < N )
    {
        double dz = double(Z[i]) - mu;
        var += dz * dz * p[i];
        i = i + 1;
    }

    double sigma_min = 1.0e-8;
    if ( var > 0.0 )
    {
        sigma = sqrt( var );
    }
    else
    {
        sigma = sigma_min;
    }

    if ( sigma < sigma_min )
    {
        sigma = sigma_min;
    }

    double num = 0.0;
    double den = 0.0;

    i = 0;
    while ( i < N )
    {
        const double dz   = double(Z[i]) - mu;
        const double arg  = -0.5 * (dz * dz) / (sigma * sigma);
        double g = 0.0;

        if ( arg > -700.0 )
        {
            g = exp(arg);
        }
        else
        {
            g = 0.0;
        }

        num += p[ i ] * g;
        den += g * g;
        i = i + 1;
    }
}

inline double CGrainCharge::Y_band( double E_eV, double a_m, int Z)
{
    double th  = Theta_eV(E_eV, a_m, Z);
    double y0v = y0_bulk(th);
    double y1v = y1_smallgrain(a_m);
    double y2v = y2_escape(E_eV, a_m, Z);

    double y01 = y0v * y1v;

    if(y01 > 1.0)
        y01 = 1.0;

    return y2v * y01;
}

// Placeholder for EUV/X-ray
inline double CGrainCharge::Y_inner( double E_eV, double a_eff, int Z)
{
    return 0.0;
}

void CGrainCharge::init(CDustComponent * _dust, double _a_min, double _a_max, uint _Na, double _a0, double _a_sigma)
{
    clear();      

    a_min = _a_min;
    a_max = _a_max;
    Na = _Na;
    a0 = _a0;
    a_sigma = _a_sigma;
    
    dust = _dust;

    Nla = max(Na,uint(20));

    arr_a_eff = new double[Na];    
    arr_dnda = new double[Na];  

    arr_a_eff_large = new double[Nla];    
    arr_dnda_large = new double[Nla];    

    double log_amin = log10(a_min);
    double log_amax = log10(a_max);
    double delta_a = (log_amax - log_amin) / double(Na-1);
    double sum = 0;

    for (int i = 0; i < Na; i++)
    {
        double tmp_a_eff= pow(10.0, log_amin + i * delta_a);
        double tmp1 = log(tmp_a_eff/a0);
        double tmp2 = 1.0 / tmp_a_eff * exp( -tmp1*tmp1 / (2*a_sigma *a_sigma) );

        sum+=tmp2;

        arr_a_eff[i] = tmp_a_eff;
        arr_dnda[i] = tmp2;
    }

    arr_a_eff[0]=a_min;
    arr_a_eff[Na-1]=a_max;

    for (int i = 0; i < Na; i++)
        arr_dnda[i] /= sum;

    delta_a = (log_amax - log_amin) / double(Nla-1);
    sum = 0;

    for (int i = 0; i < Nla; i++)
    {
        double tmp_a_eff= pow(10.0, log_amin + i * delta_a);
        double tmp1 = log(tmp_a_eff/a0);
        double tmp2 = 1.0 / tmp_a_eff * exp( -tmp1*tmp1 / (2*a_sigma *a_sigma) );

        sum+=tmp2;

        arr_a_eff_large[i] = tmp_a_eff;
        arr_dnda_large[i] = tmp2;
    }

    arr_a_eff_large[0]=a_min;
    arr_a_eff_large[Nla-1]=a_max;

    for (int i = 0; i < Nla; i++)
        arr_dnda_large[i] /= sum;       
}

double CGrainCharge::J_photodetachment(const CGridBasic * grid, const cell_basic & cell, int Z, double a_eff )
{
    if ( Z >= 0 )
        return 0.0;

    double Epdt_eV = E_pdt_eV(a_eff, Z);
    double l0      = lambda_min;
    double l1      = ( con_h * con_c ) / ( Epdt_eV * eV_to_J );

    if ( l1 > lambda_max ) l1 = lambda_max;
    if ( l1 <= l0 ) return 0.0;
        
    double sum = 0.0;

    for ( int iw = 1; iw <  wavelength_list.size(); iw++)
    {
        double lp = wavelength_list[iw];
        double E_eV_p   = (con_h * con_c / lp) * J_to_eV;

        double ulam_p   = get_u_lam(grid,cell,iw);
        double sigma_p  = sigma_pdt_m2(E_eV_p, a_eff, Z);

        double yp = sigma_p * ( lp / con_h ) * ulam_p;

        //if(lp>lambda_max)
        //    break; 

        if (lp > l1)
            break;  

        double ln = wavelength_list[iw-1];
        double E_eV_n   = (con_h * con_c / ln) * J_to_eV;
        
        double ulam_n   = get_u_lam(grid,cell,iw-1);
        double sigma_n  = sigma_pdt_m2(E_eV_n, a_eff, Z);

        double yn = sigma_n * ( ln / con_h ) * ulam_n;

        sum += (lp - ln) * yn + 0.5 * (lp - ln) * (yp - yn);
    }

    if ( sum < 0.0 )
        sum = 0.0;

    return sum;
}    

double CGrainCharge::J_photoelectric(const CGridBasic * grid, const cell_basic & cell, int Z, double a_eff)
    {
        double Epet_eV = E_pet_eV( a_eff, Z);
        double l0      = lambda_min;
        double l1      = ( con_h * con_c ) / ( Epet_eV * eV_to_J );

        if ( l1 > lambda_max ) l1 = lambda_max;
        if ( l1 <= l0 ) return 0.0;

        //double dl  = ( l1 - l0 ) / double( N_lambda );
        double sum = 0.0;

        for(int iw = 1; iw < wavelength_list.size(); ++iw )
        {
            double lp = wavelength_list[iw];
            double E_eV_p   = ( con_h * con_c / lp ) * J_to_eV;
            double Y_p      = get_Yield(E_eV_p, a_eff, Z);
            
            //double Qabs_p   = test_Qabs_lambda( lp, a_eff );
            double Qabs_p   = dust->getQabsMean(0,iw);

            double ulam_p   =  get_u_lam(grid,cell,iw);
            
            
            double yp = Y_p * Qabs_p * (lp / con_h) * ulam_p;
            
            if(lp>l1)
                break;     
            
            double ln = wavelength_list[iw-1];
            double E_eV_n   = ( con_h * con_c / ln ) * J_to_eV;
            double Y_n      = get_Yield(E_eV_n, a_eff, Z);
            //double Qabs_n   = test_Qabs_lambda( ln, a_eff );
            double Qabs_n   = dust->getQabsMean(0,iw-1);
            
            double ulam_n   = get_u_lam(grid,cell,iw-1);
            
            double yn = Y_n * Qabs_n * (ln / con_h) * ulam_n;

            sum += (lp - ln) * yn + 0.5 * (lp - ln) * (yp - yn);
        }

        double Jpe = PI * a_eff * a_eff * sum;
        
        if ( Jpe < 0.0 )
            Jpe = 0.0;
        
        return Jpe;
    }

// Placeholder for inner-shell Auger
inline double CGrainCharge::Y_auger( double E_eV, double a_eff, int Z)
{        
    return 0.0;
}

inline double CGrainCharge::vth_pref( double m )
{
    return sqrt( ( 8.0 * con_kB * Tgas ) / ( PI * m ) );
}

double CGrainCharge::J_electron( int Z, double a_eff)
{
    double q    = -con_e;
    double tau  = get_tau(a_eff, Tgas, fabs(q));
    double nu   = get_nu(Z, q);
    double Jt   = Jtilde(tau, nu);
    double rate = n_el * stick_e * PI * a_eff * a_eff * vth_pref(con_m_e) * Jt;

    if ( rate < 0.0 )
        rate = 0.0;

    return rate;
}

double CGrainCharge::J_ion_Hp( int Z, double a_eff )
{
    double q    = con_e;
    double tau  = get_tau(a_eff, Tgas, fabs(q));
    double nu   = get_nu (Z, q);
    double Jt   = Jtilde(tau, nu);
    double rate = n_ion * stick_ion * PI * a_eff * a_eff * vth_pref(con_m_p) * Jt;

    if( rate < 0.0 )
        rate = 0.0;

    return rate;
}

double CGrainCharge::sigma_pdt_m2( double E_eV, double a_m, int Z)
{
    if ( Z >= 0 ) return 0.0;

    double Epdt   = E_pdt_eV( a_m, Z);
    double DeltaE = 3.0; // eV
    double x      = ( E_eV - Epdt ) / DeltaE;

    if ( x <= 0.0 ) return 0.0;

    // 1.2e-17 cm^2 * |Z| * x / ( 1 + x^2 / 3 )^2
    double sigma_cm2 = 1.2e-17 * double(abs( Z )) * x;
    double denom     = 1.0 + ( x * x / 3.0 );
    denom = denom * denom;

    if ( denom <= 0.0 ) return 0.0;

    sigma_cm2 = sigma_cm2 / denom;

    double sigma_m2 = sigma_cm2 * 1.0e-4; // cm^2 -> m^2
    return sigma_m2;
}

// Placeholder for secondary electrons
double CGrainCharge::Y_secondary( double E_eV, double a_m, int Z)
{
    return 0.0;
}

//WD01 WDB06
double CGrainCharge::get_Yield(double E_eV, double a_m, int Z)
{
    double y  = 0.0;

    y += Y_band(E_eV, a_m, Z);
    y += Y_inner(E_eV, a_m, Z);
    y += Y_auger(E_eV, a_m, Z);
    y += Y_secondary(E_eV, a_m, Z);

    if ( y < 0.0 ) y = 0.0;
    if ( y > 1.0 ) y = 1.0;

    return y;
}

inline double CGrainCharge::Theta_eV( double E_eV, double a_m, int Z)
{
    double Epet = E_pet_eV( a_m, Z);

    if ( Z >= 0 )
    {
        double Ec = get_ECoul_eV(a_m);
        return max(0.0, E_eV - Epet + ( double( Z ) + 1.0 ) * Ec );
    }

    return max( 0.0, E_eV - Epet );
}

double CGrainCharge::y0_bulk(double Theta)
{
    if ( Theta <= 0.0 ) 
        return 0.0;

    double W;
    double TW;
    double TW5;
    double num;
    double den;

    /*if ( carbon )
    {
        W   = W_carbon_eV;
        TW5 = pow( Theta / W, 5.0 );
        num = 9.0e-3 * TW5;
        den = 1.0 + 3.7e-2 * TW5;
        return num / den;
    }*/

    //W   = W_sil_eV;

    W   = get_work_function_eV();
    TW  = Theta / W;
    num = 0.5 * TW;
    den = 1.0 + 5.0 * TW;

    if ( den <= 0.0 ) 
        return 0.0;

    return num / den;
}

// Updated WD01: fixed lengths
double CGrainCharge::y1_smallgrain( double a_m )
{

    const double le_AA = 10.0;
    const double la_AA = 100.0;

    double a_AA  = a_m * m_to_AA;
    double alpha = a_AA / la_AA + a_AA / le_AA;
    double beta  = a_AA / la_AA;

    if ( alpha <= 0.0 || beta <= 0.0 ) 
        return 0.0;

    double expa = exp( -alpha );
    double expb = exp( -beta );

    double num  = ( alpha * alpha ) - ( 2.0 * alpha ) + 2.0 - ( 2.0 * expa );
    double den  = ( beta  * beta  ) - ( 2.0 * beta  ) + 2.0 - ( 2.0 * expb );

    if ( den <= 0.0 )
        return 0.0;

    double fac  = ( beta / alpha );
    fac = fac * fac;

    double y1 = fac * ( num / den );

    return y1;
}

inline double CGrainCharge::Elow_eV( double a_m, int Z )
{
    if(Z < 0) 
        return Emin_eV( a_m, Z );

    return - ( double( Z ) + 1.0 ) * get_ECoul_eV( a_m );
}

inline double CGrainCharge::Ehigh_eV( double E_eV, double a_m, int Z)
{
    if ( Z < 0 )
    {
        return Emin_eV(a_m, Z) + E_eV - E_pet_eV(a_m, Z);
    }

    return E_eV - E_pet_eV( a_m, Z);
}

double CGrainCharge::y2_escape( double E_eV, double a_m, int Z)
{
    if ( Z < 0 )
        return 1.0;

    double Eh = Ehigh_eV( E_eV, a_m, Z);
    double El = Elow_eV ( a_m, Z );

    if ( Eh <= El )
        return 0.0;

    if ( Eh < 0.0 )
    {
        double En = fabs( Eh );
        El += En;
        Eh  = 0.0;
    }

    double diff = Eh - El;
    if(diff <= 0.0)
        return 0.0;

    double num = ( Eh * Eh ) * ( Eh - 3.0 * El );
    double den = diff * diff * diff;

    if ( den <= 0.0 )
        return 0.0;

    double y2 = num / den;
    y2 = clamp_value(y2, 0.0, 1.0 );

    return y2;
}

inline double CGrainCharge::EA_eV( double a_m, int Z )
{
    double W  = get_work_function_eV();
    double Ec = get_ECoul_eV( a_m );
    return W + (double( Z ) - 0.5 ) * Ec;
}

double CGrainCharge::IPv_eV( double a_m, int Z)
{
    double W  = get_work_function_eV();
    double Ec = get_ECoul_eV( a_m );

    if ( Z >= 0 )
    {
        return W + ( double(Z) + 0.5 ) * Ec;
    }
    else
    {
        // WD01: for Z < 0, IP_v(Z) = EA(Z+1)
        return EA_eV( a_m, Z + 1);
    }
}

double CGrainCharge::Emin_eV( double a_m, int Z )
{
    if ( Z >= -1 ) return 0.0;

    double zprime = static_cast<double>( abs( Z + 1 ) );
    double nu     = zprime;
    double theta  = nu / ( 1.0 + ( 1.0 / sqrt( nu ) ) );
    double a_AA   = a_m * m_to_AA;

    double pow1 = pow( a_AA / 10.0, -0.45 );
    double pow2 = pow( max( 1.0, zprime ), -0.26 );
    double corr = 1.0 - 0.3 * pow1 * pow2;
    corr = clamp_value(corr, 0.0, 1.0);

    return theta * get_ECoul_eV( a_m ) * corr;
}

inline double CGrainCharge::E_pet_eV( double a_m, int Z)
{
    if ( Z >= -1 )
        return IPv_eV( a_m, Z);

    return IPv_eV(a_m, Z) + Emin_eV(a_m, Z);
}

// WD01 2.3.3
inline double CGrainCharge::E_pdt_eV( double a_m, int Z)
{        
    return EA_eV( a_m, Z + 1) + Emin_eV( a_m, Z );
}

inline double CGrainCharge::get_tau( double a, double T, double qabs )
{
    return ( (PIx4 * con_epsilon_0) * a * con_kB * T ) / ( qabs * qabs );
}

inline double CGrainCharge::get_nu( int Z, double q )
{
    return ( double( Z ) * con_e ) / q;
}

inline double CGrainCharge::Jtilde_0( double tau )
{
    return 1.0 + sqrt( PI / ( 2.0 * tau ) );
}

double CGrainCharge::Jtilde_neg( double tau, double nu )
{
    // attractive (nu < 0)
    double denom = tau - 2.0 * nu;
    if ( denom <= 0.0 ) return 0.0;

    double term1 = 1.0 - ( nu / tau );
    double term2 = 1.0 + sqrt( 2.0 / denom );
    return term1 * term2;
}

double CGrainCharge::Jtilde_pos( double tau, double nu )
{
    // repulsive (nu > 0)
    double inner = 1.0 / ( 4.0 * tau + 3.0 * nu );
    double root  = 0.0;
    if ( inner > 0.0 )
    {
        root = sqrt( inner );
    }

    double theta = 0.0;
    if ( nu > 0.0 )
    {
        double sq = sqrt( nu );
        double denom = 1.0 + ( 1.0 / sq );
        theta = nu / denom;
    }

    double pref = ( 1.0 + root );
    pref = pref * pref;

    double result = pref * exp( -theta / tau );
    return result;
}

inline double CGrainCharge::Jtilde( double tau, double nu )
{
    if ( nu < 0.0 ) return Jtilde_neg( tau, nu );
    if ( nu == 0.0 ) return Jtilde_0 ( tau );
    return Jtilde_pos( tau, nu );
}

inline double CGrainCharge::get_ECoul_eV(double a_eff)
{
    double val = ( con_e * con_e ) / ( (PIx4 * con_epsilon_0) * a_eff );
    return val * J_to_eV;
}

void CGrainCharge::set_ID(uint id)
{
    mix_ID = id;
}

void CGrainCharge::set_work_function_eV(double w)
{
    W_eV=w;
}

void CGrainCharge::set_work_function_J(double w)
{
    W_eV=w*J_to_eV;
}

void CGrainCharge::set_density(double r)
{
    rho = r;
}

void CGrainCharge::set_m0_rel(double m)
{
    m0 = m*con_m_p;
}

void CGrainCharge::set_beta(double b)
{
    beta = b;
}

void CGrainCharge::set_Smax(double S)
{
    Smax = S;
}

uint CGrainCharge::getGrainSizes()
{
    return Na;
}

inline double CGrainCharge::get_work_function_eV()
{
    return W_eV;
}

inline double CGrainCharge::clamp_value( double x, double a, double b )
{
    if ( x < a ) return a;
    if ( x > b ) return b;
    return x;
}

void CGrainCharge::clear()
{
    tot_Zmin=10000;
    tot_Zmax=-10000;

    if(arr_a_eff!=0)
    {
        delete [] arr_a_eff;    
        arr_a_eff=0;
    }

    if(arr_dnda!=0)
    {
        delete [] arr_dnda; 
        arr_dnda = 0;
    }

    if(arr_a_eff_large!=0)
    {
        delete [] arr_a_eff_large;    
        arr_a_eff_large=0;
    }

    if(arr_dnda_large!=0)
    {
        delete [] arr_dnda_large; 
        arr_dnda_large = 0;
    }
}

void CGrainCharge::build_distribution(const CGridBasic * grid, const cell_basic & cell, int ia, double & Zgr, double & Zs, double & Trot)
{    
        double tmp_a_eff=arr_a_eff[ia];

        int Zhard_min = compute_Zmin_auto    (tmp_a_eff);
        int Zhard_max = compute_Zmax_coulomb (tmp_a_eff, Smax);

        int Zcoarse = find_Zeq(grid, cell, Zhard_min, Zhard_max, tmp_a_eff);

        int Zmin = Zcoarse - 6;
        if ( Zmin < Zhard_min ) Zmin = Zhard_min;

        int Zmax = Zcoarse + 6;
        if ( Zmax > Zhard_max ) Zmax = Zhard_max;

        int guard = 0;
        while ( guard < 50 )
        {
            dlist ftmp = compute_f_window(grid, cell, Zmin, Zmax, Zhard_min, Zhard_max, tmp_a_eff);
            if ( tails_are_small( ftmp ) ) break;

            Zmin -= 6;
            Zmax += 6;

            if ( Zmin < Zhard_min ) Zmin = Zhard_min;
            if ( Zmax > Zhard_max ) Zmax = Zhard_max;

            if ( ( Zmax - Zmin ) > 1200 ) break;
            if ( Zmin == Zhard_min && Zmax == Zhard_max ) break;

            guard += 1;
        }

        if ( guard >= 50 )
        {
            Zmin = Zcoarse - Z_span_default;
            Zmax = Zcoarse + Z_span_default;
            if ( Zmin < Zhard_min ) Zmin = Zhard_min;
            if ( Zmax > Zhard_max ) Zmax = Zhard_max;
        }

        dlist f = compute_f_window(grid, cell, Zmin, Zmax, Zhard_min, Zhard_max, tmp_a_eff);
        int N = ( Zmax - Zmin + 1 );

        ilist Zvals(N, 0);

        //double Trot=0;

        double Gn=0, Fn=0;
        double Gi=0, Fi=0;
        double Gp=0, Fp=0;
        double Gpe=0, Fpe=0;
        double GIR=0, FIR=0;
        double tau_H=calc_tau_H(tmp_a_eff);
        double tau_ed=0;

        for ( int k = 0; k < N; ++k )
        {
            int Z = Zmin + k;
            Zvals[k] = Z;

            double tmp_Gn=0, tmp_Fn=0;
            double tmp_Gi=0, tmp_Fi=0;
            double tmp_Gp=0, tmp_Fp=0;
            double tmp_Gpe=0, tmp_Fpe=0;

            calc_FGn(tmp_Fn, tmp_Gn, tmp_a_eff, Z);
            calc_FGi(tmp_Fi, tmp_Gi, tmp_a_eff, Z);
            calc_FGp(tmp_Fp, tmp_Gp, tmp_a_eff, Z);
            calc_FGpe(grid, cell, tmp_Fpe, tmp_Gpe, tmp_a_eff, Z);


            Fn+=f[k]*tmp_Fn;
            Gn+=f[k]*tmp_Gn;

            Fi+=f[k]*tmp_Fi;
            Gi+=f[k]*tmp_Gi;

            Fpe+=f[k]*tmp_Fpe;
            Gpe+=f[k]*tmp_Gpe;

            tau_ed = max( tau_ed, calc_tau_ed(tmp_a_eff,Z) );
        }

        calc_FGIR(grid, cell, FIR, GIR, tmp_a_eff);
        //compute_IR_coeffs(FFIR, GFIR, tmp_a_eff);

        double sec=20*tau_H/(3*tau_ed);

        double F=Fn + Fi + Fp + Fpe + FIR;
        double G=Gn + Gi + Gp + Gpe + GIR;

        double fr= 2 * G / F / (1+sqrt(1+G/(F*F) * sec)) ;

        Trot = Tgas * fr ;

        tot_Zmin = min(tot_Zmin,Zmin);
        tot_Zmax = max(tot_Zmax,Zmax);

        double mu=0;
        double sigma=0;

        fit_gaussian_pdf(Zvals,f,mu,sigma);

        Zgr = mu;
        Zs = sigma;
}

inline int CGrainCharge::findIndex(double a) const
{
    if(a<=arr_a_eff[0])
        return 0;

    if(a>=arr_a_eff[Na-1])
        return Na - 2;

    int low = 0;
    int high = Na - 2;

    while (low <= high)
    {
        int mid = (low + high) / 2;

        if(arr_a_eff[mid] <= a && a < arr_a_eff[mid + 1])
        {
            return mid;
        }
        else if(a < arr_a_eff[mid])
        {
            high = mid - 1;
        }
        else
        {
            low = mid + 1;
        }
    }

    return Na - 2;
}

void CGrainCharge::print_distribution()
{
    /*for(int iZ=tot_Zmin; iZ<=tot_Zmax; iZ++)
    {
        for(int ia=0; ia<Na; ia++)
        {
            double mu = arrZmean[ia];
            double sigma = arrZsig[ia];
            double f = gaussian_value( mu, sigma, iZ);

            arr_sum[ia]+=f;
        }
    }

    cout << "\n\n\n";

    for(int iZ=tot_Zmin; iZ<=tot_Zmax; iZ++)
    {
        cout << iZ << "\t";
        for(int ia=0; ia<Na; ia++)
        {
            double mu = arrZmean[ia];
            double sigma = arrZsig[ia];
            double f = gaussian_value( mu, sigma, iZ)/arr_sum[ia];

            cout << f << "\t";
        }

        cout << "\n";
    }

    cout << "\n\n\n";

    for(int ia=0; ia < Na; ia++)
    {
        cout << arr_a_eff[ia] << " " << arr_dnda[ia] << " " << arrZmean[ia] << " " << arrZsig[ia] << " " << arrTrot[ia] << "\n";        
    }
    cout << "\n";*/
}

inline double CGrainCharge::g1(double x)
{
    if(x<1)
        return 1-x;

    return exp(-x);    
}

inline double CGrainCharge::g2(double x)
{
    if(x<1)
        return 1-x+0.5*x*x;

    return exp(-x);    
}

double CGrainCharge::Gamma_photoelectric(const CGridBasic * grid, const cell_basic & cell, int Zgr, double a_eff)
{
    const double area = PI * a_eff * a_eff;

    double I = 0;
    double EpetJ = E_pet_eV(a_eff, Zgr) * eV_to_J;               // WD01 threshold in J

    for (int iw = 1; iw < wavelength_list.size(); iw++)
    {
        double lp   = wavelength_list[iw];
        double up  = get_u_lam(grid,cell,iw); // J m^-4 (energy density per wavelength)

        double E_eVp  = (con_h * con_c / lp) * J_to_eV;      // photon energy in eV
        double Yp     = get_Yield(E_eVp, a_eff, Zgr); // band + (placeholders for inner/Auger/secondary)
        double Qp     = dust->getQabsMean(0,iw);
        double Ekinp  = (con_h * con_c / lp) - EpetJ;       // mean kinetic energy (simple approx)

        double ln   = wavelength_list[iw-1];
        double un  = get_u_lam(grid,cell,iw-1);

        double E_eVn  = (con_h * con_c / ln) * J_to_eV;
        double Yn     = get_Yield(E_eVn, a_eff, Zgr);
        double Qn     = dust->getQabsMean(0,iw);
        double Ekinn  = (con_h * con_c / ln) - EpetJ;               

        if(Ekinp<0)
            Ekinp=0;

        if(Ekinn<0)
            Ekinn=0;

        double yp = Qp * up * lp * Yp * Ekinp;
        double yn = Qn * un * ln * Yn * Ekinn;

        I += (lp - ln) * yn + 0.5 * (lp - ln) * (yp - yn);
    }

    double Gamma_pe = (area / con_h) * I; // [W]

    if(Gamma_pe<0)
        Gamma_pe=0;

    return Gamma_pe;
}

void CGrainCharge::calc_FGpe(const CGridBasic * grid, const cell_basic & cell, double & Fpe, double & Gpe, double a_eff, int Zgr)
{
    Fpe = 0.0;
    Gpe = 0.0;

    // 1) photoelectron number rate and heating power
    const double Jpe = J_photoelectric(grid, cell, Zgr, a_eff);        // [s^-1] (already implemented in class)
    const double GpeW = Gamma_photoelectric(grid, cell, Zgr, a_eff);   // [W]    (computed above)

    /*if(Jpe>0)
        int tt=0;
    
    if(GpeW>0)
        int tt=0;*/

    const double mH = con_m_p;
    const double vF = sqrt( 2.0 * con_kB * Tgas / (PI * mH) );
    const double denomF = 2.0 * PI * a_eff * a_eff * n_gas * vF;

    const double coulombJ = ( (double(Zgr) + 1.0) * con_e * con_e ) / ( 4.0 * PI * con_epsilon_0 * a_eff ); // J = e Φ_g

    const double vG = sqrt( 8.0 * PI * mH * con_kB * Tgas );                // sqrt(8 π m_H k T)
    const double denomG = 4.0 * n_gas * a_eff * a_eff * con_kB * Tgas * vG;

    // 3) assemble
    if (denomF > 0.0)
        Fpe = (con_m_e / mH) * (Jpe / denomF);

    if (denomG > 0.0)
        Gpe = con_m_e * ( GpeW + coulombJ * Jpe ) / denomG;

    if (Fpe < 0.0) Fpe = 0.0;
    if (Gpe < 0.0) Gpe = 0.0;
}

double CGrainCharge::get_u_lam(const CGridBasic * grid, const cell_basic & cell, uint iw)
{
    double vol = grid->getVolume(cell);

    double arr_en_dens = 0;
    Vector3D en_dir;
        
    grid->getSpecLength(cell, iw, &arr_en_dens, &en_dir);
    
    // arr_en_dens = 4 * PI * vol * J -> 4 * PI / c * J
    double u_lam   = arr_en_dens / double(vol * con_c);
    
    if(u_lam>0)
        int tt=0;
    
    return u_lam;
}

// ======================================================
// IR emission recoil: F_IR, G_IR  (SI, wavelength form)
// ======================================================
void CGrainCharge::calc_FGIR(const CGridBasic * grid, const cell_basic & cell, double & FIR, double & GIR, double a_eff)
{
    FIR = 0.0;
    GIR = 0.0;

    // thermal speed of neutral H (SpDust convention)
    double vth = sqrt( 8.0 * con_kB * Tgas / (PI * con_m_p) );

    double IF = 0;//trapz_arr(arr_lambda, iF.data(), N_lambda);
    double IG = 0;//trapz_arr(arr_lambda, iG.data(), N_lambda);

    for (int iw = 1; iw < wavelength_list.size(); ++iw)
    {
        double lp = wavelength_list[iw];        // [m]

        double Qp  = dust->getQabsMean(0,iw);     // dimensionless
        double Blp = CMathFunctions::planck(lp, Tdust);            // W m^-3 sr^-1
        double wFp = (lp * lp) / (con_c * con_c); // λ^2 / c^2
        double wGp = lp / con_c;                   // λ / c

        double ln = wavelength_list[iw-1];

        double Qn  = dust->getQabsMean(0,iw-1);
        double Bln = CMathFunctions::planck(ln, Tdust);    
        double wFn = (ln * ln) / (con_c * con_c); 
        double wGn = ln / con_c;                  

        double yFp = Qp * Blp * wFp;
        double yFn = Qn * Bln * wFn;

        double yGp = Qp * Blp * wGp;
        double yGn = Qn * Bln * wGn;

        IF += (lp - ln) * yFn + 0.5 * (lp - ln) * (yFp - yFn);
        IG += (lp - ln) * yGn + 0.5 * (lp - ln) * (yGp - yGn);
    }


    //double tmp = 0.1 / ( n_gas * con_m_p * vth * a_eff * a_eff );
    double tmp = 1.0 / ( n_gas * con_m_p * vth * a_eff * a_eff );

    FIR = (3.0 / (2.0 * PI))            * tmp * IF;
    GIR = (con_hq / (8.0 * PI * con_kB * Tgas)) * tmp * IG;
}

void CGrainCharge::calc_FGn(double & Fn, double & Gn, double a_eff, double Zgr)
{

    double epsilon_n_sq = (con_e*con_e) / (PIx4 * con_epsilon_0) * Zgr * Zgr * alpha_gas /(2*pow(a_eff,4)*con_kB*Tgas);
    double epsilon_e_sq = (con_e*con_e) / (PIx4 * con_epsilon_0) * Zgr * Zgr * alpha_gas /(2*pow(a_eff,4)*con_kB*Tdust);

    double epsilon_n = sqrt(epsilon_n_sq);
    double epsilon_e = sqrt(epsilon_e_sq);

    double Tev=Tdust;

    double tmp1 = exp(-epsilon_n_sq)+2*epsilon_n_sq;
    double tmp2 = Tev / Tgas *(exp(-epsilon_n_sq) + PIsq * epsilon_n*erf(epsilon_n) ) / (exp(-epsilon_e_sq) + PIsq * epsilon_e*erf(epsilon_e) );
    double tmp3 = exp(-epsilon_e_sq)+2*epsilon_e_sq;

    //double 

    Fn =  n_neu / n_gas * sqrt(mu_mol) * (exp(-epsilon_n_sq) + PIsq * epsilon_n*erf(epsilon_n)  );

    Gn = n_neu / (2*n_gas) * sqrt(mu_mol)* (tmp1 + tmp2*tmp3) ;
}

/*void calc_FGpe(double & Fpe, double & Gpe, double a_eff, double Zgr)
{
    double Jpe = J_photoelectric(Zgr, a_eff);

    Fpe=(con_m_e / con_m_p) * Jpe / (PIx2 *a_eff*a_eff*n_gas * sqrt( 2*con_kB*Tgas / (PI*con_m_p)));
    Gpe = 0;
}*/

void CGrainCharge::calc_FGi(double & Fi, double & Gi, double a_eff, double Zgr)
{
    double Gin, Gev;

    if(Zgr*Zgas!=0)
    {
        double psi = (con_e*con_e) / (PIx4 * con_epsilon_0) * Zgr*Zgas / (a_eff*con_kB*Tgas);


        double epsilon_i_sq = (con_e*con_e) / (PIx4 * con_epsilon_0) * alpha_gas /(2*pow(a_eff,4)*con_kB*Tdust);
        double epsilon_i = sqrt(epsilon_i_sq);

        double tmp1 = exp( -Zgr*Zgr*epsilon_i_sq)+2*Zgr*Zgr*epsilon_i_sq;
        double tmp2 = exp( -Zgr*Zgr*epsilon_i_sq)+abs(Zgr)*PIsq*epsilon_i*erf(abs(Zgr)*epsilon_i);

        Fi = n_ion / n_gas * sqrt(mu_mol) * g1(psi);

        double Tev=Tdust;

        Gin= n_ion / (2*n_gas) *sqrt(mu_mol)  * g2(psi);
        Gev= Fi * ( Tev/(2*Tgas) ) * tmp1/tmp2;

        Gi = Gin + Gev;
    }
    else
    {
        double phi_sq = (con_e*con_e) / (PIx4 * con_epsilon_0) * 2*Zgas*Zgas / (a_eff*con_kB*Tgas);
        double phi = sqrt(phi_sq);

        Fi =  n_ion / n_gas * sqrt(mu_mol) * (1+PIsq/2*phi);
        Gin= n_ion / (2*n_gas) * sqrt(mu_mol) * (1+ 3*PIsq/4*phi+0.5*phi_sq);
        Gev = Tdust/(2*Tgas)*Fi;
        Gi = Gin + Gev;
    }
}

double CGrainCharge::calc_mu_dipole(double Zgr, double a_eff)
{
    double epsilon = 0.1;
    //double rho = rho_sil;
    //double m=24*con_m_p;
    //double beta = beta_sil;

    /*if(carbon)
    {
        rho = rho_carb;
        m=12*con_m_p;
        beta = beta_carb;
    }*/

    double N=PIx4 * rho * pow(a_eff,3) / m0;
    double mu_int = beta * sqrt(N);
    double mu_charge = epsilon*Zgr*con_e*a_eff;       

    double mu = sqrt(mu_int*mu_int + mu_charge*mu_charge);

    return mu;
}

void CGrainCharge::calc_FGp(double & Fp, double & Gp, double a_eff, double Zgr)
{
    /*double rho=rho_sil;

    if(carbon)
        rho=rho_carb;*/

    double xi = 1.0; //for spheres DL98 Eq. A4
    double I = 8./15. * PI * rho * pow(a_eff,5);//for spheres
    double v_th = sqrt(2*con_kB*Tgas/(mu_mol * con_m_p));
    double omega_th = sqrt(2*con_kB*Tgas/I);

    double lambda_D = sqrt( (con_epsilon_0 * con_kB) * Tgas / (n_el * con_e* con_e)  );
    double b_omega = v_th / omega_th;
    double b_q = I * v_th / con_hq; 

    double mu_dip = calc_mu_dipole(Zgr, a_eff);
    double tmp1 = n_ion / n_gas * sqrt(mu_mol) * (con_e*con_e) / (PIx4 * con_epsilon_0) * 2 * Zgas*Zgas / (3*pow(a_eff,4)*pow(con_kB*Tgas,2)) *mu_dip*mu_dip;

    // <cos²> = 1/3
    double tmp2 = log(b_omega / a_eff)+1.0/3.0*log(min(b_q,lambda_D)/b_omega);

    Gp = tmp1*tmp2;

    Fp = Gp;        
}

/*void CGrainCharge::compute_IR_coeffs(double & F_IR, double & G_IR, double a_eff)
{
    // Constants and sanity guards
    const double nH_ref = 2.0e7;   // m^-3  (20 cm^-3)
    const double T_ref  = 100.0;   // K
    const double a_unit = 1.0e-9;  // m (10^-7 cm)
    const double U_isrf = 8.64e-14;// to avoid division by zero


    // Dimensionless size a_-7 = a / (10^-9 m)
    const double a_m7 = a_eff / a_unit;

    double U_ratio = u_tot / U_isrf;

    // Common factors
    const double nH_fac = nH_ref / n_gas;
    const double T_fac_sqrt = sqrt(T_ref / Tgas);      // (100 K / T)^{1/2}
    const double T_fac_32   = pow(T_ref / Tgas, 1.5);  // (100 K / T)^{3/2}

    // ---- F_IR ----
    // Continuous-T regime
    const double F_IR_c = 60.8
        * (1.0 / a_m7)
        * pow(U_ratio, 2.0/3.0) * nH_fac * T_fac_sqrt;

    // Spike (quantized) regime
    const double F_IR_q = 4.49
        * sqrt(a_m7)
        * U_ratio
        * nH_fac
        * T_fac_sqrt;

    // ---- G_IR ----
    // Continuous-T regime
    const double G_IR_c = 7.34
        * (1.0 / a_m7)
        * pow(U_ratio, 5.0/6.0)
        * nH_fac
        * T_fac_32;

    // Spike regime
    const double G_IR_q = 2.11
        * pow(1.0 / a_m7, 0.25)
        * U_ratio
        * nH_fac
        * T_fac_32;

    // Pick regime by transition size (and min for robustness)

    F_IR = min(F_IR_q, F_IR_c);
    G_IR = min(G_IR_q, G_IR_c);
}*/

/*void calc_FGIR(double & FIR, double & GIR, double a_eff)
{


    double I = 8./15. * PI * rho * pow(a_eff,5);//for spheres
    double v_th = sqrt(8*con_kB*Tgas/(8* mu_mol * con_m_p));

    //double tau_gas

    double tmp_FIR=3* con_h*con_h /(n_gas*sqrt(con_m_p * con_kB * Tgas)) / pow(a_eff * con_kB * Tdust,2);// / v_th;
    double tmp_GIR= con_h*con_h /(16*n_gas*sqrt(con_m_p)*pow(con_kB * Tgas,3./2.)) / (a_eff * a_eff * con_kB * Tdust);// / v_th;

    double tmp=0;

    for(uint i = 1; i < N_lambda; i++)
    {
        double lp = arr_lambda[i];
        double qp=lp*lp*planck(lp,Tdust)*Qabs_lambda(lp, a_eff)/(con_c);

        double ln = arr_lambda[i-1];
        double qn=ln*ln*planck(ln,Tdust)*Qabs_lambda(ln, a_eff)/(con_c);

        tmp += (lp - ln) * qn + 0.5 * (lp - ln) * (qp - qn);
    }

    FIR = tmp * tmp_FIR;
    GIR = tmp * tmp_GIR;
}*/

inline double CGrainCharge::calc_tau_H(double a_eff)
{
    /*double rho=rho_sil;

    if(carbon)
        rho=rho_carb;*/

    double vth=sqrt(8*con_kB*Tgas / (PI*con_m_p));
    double tau_H = 3. / (4.*PIsq) * rho*a_eff /(n_gas*con_m_p)/vth;

    return tau_H;
}

inline double CGrainCharge::calc_tau_ed(double a_eff, double Zgr)
{
    /*double rho=rho_sil;

    if(carbon)
        rho=rho_carb;*/

    double I = 8./15. * PI * rho * pow(a_eff,5);//for spheres
    double mu_dip = calc_mu_dipole(Zgr, a_eff);

    // is it I or I*I ?
    double tau_ed = 3.*I*pow(con_c,3)*con_epsilon_0 / (4.*mu_dip*mu_dip * con_kB *Tgas);

    return tau_ed;
}

inline double CGrainCharge::Power(double lambda, double a_eff, double Z)
{
    double mu_dip=calc_mu_dipole(Z, a_eff);
    double res=mu_dip*mu_dip / (18*PI*con_epsilon_0*con_c*con_c*con_c)*pow(PIx2*con_c/lambda,4.0);


    return res;    
}

inline double CGrainCharge::fMW(double lambda, double Trot, double a_eff)
{
    /*double rho=rho_sil;

    if(carbon)
        rho=rho_carb;*/

    double I = 8./15. * PI * rho * pow(a_eff,5);//for spheres
    double tmp1 = I / (2 * con_kB * Trot);
    double tmp2 = (PIx2 * con_c) / lambda;
    double res=4.0 / PIsq;

    res *= pow( tmp1 , 3./2.);
    res *= tmp2*tmp2;
    res *= exp(-tmp1*tmp2*tmp2);

    return res;
}

inline double CGrainCharge::interpolate(double x, double x1, double x2, double y1, double y2)
{
    double y = y1 + (y2 - y1) * (x - x1) / (x2 - x1);
    return y;
}

double CGrainCharge::j_lambda_inter(const CGridBasic * grid, const cell_basic & cell, double lambda)
{
    double jl=0;

    for(int ia=0; ia<Nla; ia++)
    {
        double Zmin=0, Zmax=0;
        double sum=0;
        double tmp_j=0;

        double a_eff = arr_a_eff_large[ia];

        int i_ame=findIndex(a_eff);

        double a_lower = arr_a_eff[i_ame];
        double Zsig_lower = grid->getAMEZgr(cell,i_ame);
        double Zmean_lower = grid->getAMEZs(cell,i_ame);
        double Trot_lower = grid->getAMETrot(cell,i_ame);

        double a_upper = arr_a_eff[i_ame+1];
        double Zsig_upper = grid->getAMEZgr(cell,i_ame+1);
        double Zmean_upper = grid->getAMEZs(cell,i_ame+1);
        double Trot_upper = grid->getAMETrot(cell,i_ame+1);

        double Zsig = interpolate(a_eff, a_lower, a_upper, Zsig_lower, Zsig_upper);
        double Zmean = interpolate(a_eff, a_lower, a_upper, Zmean_lower, Zmean_upper);
        double Trot = interpolate(a_eff, a_lower, a_upper, Trot_lower, Trot_upper);

        /*cout << ia << " " << index << " -- " << a_eff << ", " << a_lower << ", " << a_upper << " -- "
                << Zsig << ", " << Zsig_lower << ", " << Zsig_upper << " -- "
                << Zmean << ", " << Zmean_lower << ", " << Zmean_upper << " -- "
                << Trot << ", " << Trot_lower << ", " << Trot_upper<< endl << flush;/**/

        double f_MW = fMW(lambda, Trot, a_eff);

        if(f_MW<1e-200)
            continue;

        gaussian_bounds(Zmean,Zsig,1e-6,Zmin,Zmax);

        double dnda=arr_dnda[ia];

        for(int Z=int(Zmin); Z<=int(Zmax);Z++)
        {
            double f_charge=gaussian_value(Zmean,Zsig,Z);        
            double P = Power(lambda, a_eff, double(Z));

            tmp_j += dnda * f_charge * f_MW * P * (con_c / (lambda * lambda)) / PIx4;
            sum += f_charge;
        }

        if(sum>0)
            jl += tmp_j/sum;        
    }

    return jl;
}

/*double CGrainCharge::j_lambda(double lambda)
{
    double jl=0;

    for(int ia=0; ia<Na; ia++)
    {
        double Zmin=0, Zmax=0;
        double sum=0;
        double tmp_j=0;

        double Zsig=arrZsig[ia];
        double Zmean=arrZmean[ia];
        double a_eff=arr_a_eff[ia];
        double Trot = arrTrot[ia];

        double f_MW = fMW(lambda, Trot, a_eff);

        if(f_MW<1e-200)
            continue;

        gaussian_bounds(Zmean,Zsig,1e-6,Zmin,Zmax);

        double dnda=arr_dnda[ia];

        for(int Z=int(Zmin); Z<=int(Zmax);Z++)
        {
            double f_charge=gaussian_value(Zmean,Zsig,Z);        
            double P = Power(lambda, a_eff, double(Z));

            tmp_j += dnda * f_charge * f_MW * P * (con_c / (lambda * lambda)) / PIx4;
            sum += f_charge;
        }

        if(sum>0)
            jl += tmp_j/sum;        
    }

    return jl;
}

double CGrainCharge::j_lambda_test(double lambda)
{
    double jl=0;

    for(int ia=0; ia<Na; ia++)
    {
        double Zmin=0, Zmax=0;
        double sum=0;
        double tmp_j=0;

        double Zsig=arrZsig[ia];
        double Zmean=arrZmean[ia];
        double a_eff=arr_a_eff[ia];

        double Trot = arrTrot[ia];

        double Gn=0, Fn=0;
        double Gi=0, Fi=0;
        double Gp=0, Fp=0;
        double Gpe=0, Fpe=0;
        double GIR=0, FIR=0;
        double tau_H=calc_tau_H(a_eff);
        double tau_ed=0;

        gaussian_bounds(Zmean,Zsig,1e-5,Zmin,Zmax);

        int Zhard_min = compute_Zmin_auto    (a_eff);
        int Zhard_max = compute_Zmax_coulomb (a_eff, Smax);

        if(Zmin<Zhard_min)
            Zmin=Zhard_min;

        if(Zmax>Zhard_max)
            Zmax=Zhard_max;

        double o_cr = 2. / a_eff * sqrt(Smax / rho);
        double l_cr = PIx2 * con_c / o_cr;

        if(lambda<l_cr)
            continue;

        for(int Z=int(Zmin); Z<=int(Zmax);Z++)
        {   
            double tmp_Gn=0, tmp_Fn=0;
            double tmp_Gi=0, tmp_Fi=0;
            double tmp_Gp=0, tmp_Fp=0;
            double tmp_Gpe=0, tmp_Fpe=0;

            double f_charge=gaussian_value(Zmean,Zsig,Z);  

            calc_FGn(tmp_Fn, tmp_Gn, a_eff, Z);
            calc_FGi(tmp_Fi, tmp_Gi, a_eff, Z);
            calc_FGp(tmp_Fp, tmp_Gp, a_eff, Z);
            calc_FGpe(tmp_Fpe, tmp_Gpe, a_eff, Z);

            Fn+=f_charge*tmp_Fn;
            Gn+=f_charge*tmp_Gn;

            Fi+=f_charge*tmp_Fi;
            Gi+=f_charge*tmp_Gi;

            Fp+=f_charge*tmp_Fp;
            Gp+=f_charge*tmp_Gp;

            Fpe+=f_charge*tmp_Fpe;
            Gpe+=f_charge*tmp_Gpe;

            tau_ed = max( tau_ed, calc_tau_ed(a_eff,Z) );

            sum += f_charge;
        }

        if(sum>0)
        {
            Fn /= sum;
            Gn /= sum;

            Fi /= sum;
            Gi /= sum;

            Fp /= sum;
            Gp /= sum;

            Fpe /= sum;
            Gpe /= sum;
        }

        calc_FGIR(FIR, GIR, a_eff);
        //compute_IR_coeffs(FFIR, GFIR, a_eff);

        double sec=20*tau_H/(3*tau_ed);

        double F=Fn + Fi + Fp + Fpe + FIR;
        double G=Gn + Gi + Gp + Gpe + GIR;

        double fr= 2 * G / F / (1+sqrt(1+G/(F*F) * sec)) ;

        Trot = Tgas * fr ;

        double f_MW = fMW(lambda, Trot, a_eff);

        if(f_MW<1e-200)
            continue;

        double dnda=arr_dnda[ia];

        for(int Z=int(Zmin); Z<=int(Zmax);Z++)
        {
            double f_charge=gaussian_value(Zmean,Zsig,Z);        
            double P = Power(lambda, a_eff, double(Z));

            tmp_j += dnda * f_charge * f_MW * P * (con_c / (lambda * lambda)) / PIx4;
        }

        if(sum>0)
            jl += tmp_j/sum;   

    }

    return jl;
}*/

inline void CGrainCharge::gaussian_bounds(double mu, double sigma, double rel, double & x_min, double & x_max)
{
    // |x - mu| = sigma * sqrt(-2 ln r)
    double delta = sigma * sqrt(-2.0 * log(rel));
    x_min = double(int(mu - delta-1.5));
    x_max = double(int(mu + delta+1.5));
}

