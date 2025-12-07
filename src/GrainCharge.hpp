#include "Typedefs.hpp"
#include "GridBasic.hpp"

#ifndef CGRAINCHARGE_HPP
#define CGRAINCHARGE_HPP

struct Currents
{
    double J_up;
    double J_down;
};

struct Distribution
{
    int                 Zmin;
    int                 Zmax;
    ilist               Z;
    dlist               f;
    double              mean;
    double              stddev;
};

class CDustComponent;

class CGrainCharge {
public:
    CGrainCharge();
    virtual ~CGrainCharge();
    
    //setter functions
    void set_density(double _n_gas, double _n_el, double _n_ion);
    void set_temp(double _T_gas, double _T_el, double _T_dust);
    void set_wavelengths(const dlist & _wavelengths_list);
    
    void set_r(double _Zgas, double _mu);
    
    void set_work_function_eV(double w);
    void set_work_function_J(double w);
    
    void set_density(double r);
    void set_m0_rel(double m);
    void set_beta(double b);
    void set_Smax(double S);
    
    uint getGrainSizes();
    
    void set_ID(uint id);
    
    //getter functions
    double get_work_function_eV();
    double get_ECoul_eV(double a_eff);

    // ======================= WD01 thresholds (eV) =======================
    double EA_eV( double a_m, int Z );

    double IPv_eV( double a_m, int Z);

    double Emin_eV( double a_m, int Z );

    double E_pet_eV( double a_m, int Z);

    // WD01 2.3.3
    double E_pdt_eV( double a_m, int Z);
    
    
    
    /*double J_photoelectric( int Z, double a_eff)
    {
        double Epet_eV = E_pet_eV( a_eff, Z);
        double l0      = lambda_min;
        double l1      = ( con_h * con_c ) / ( Epet_eV * eV_to_J );

        if ( l1 > lambda_max ) l1 = lambda_max;
        if ( l1 <= l0 ) return 0.0;

        double dl  = ( l1 - l0 ) / double( N_lambda );
        double sum = 0.0;

        for ( int i = 0; i < N_lambda; ++i )
        {
            double lambda = l0 + ( ( double(i) + 0.5 ) * dl);
            double E_eV   = ( con_h * con_c / lambda ) * J_to_eV;
            double Y      = Yield_WD01_WDB06( E_eV, a_eff, Z);
            double Qabs   = Qabs_lambda( lambda, a_eff );
            double ulam   = u_lambda(lambda);

            sum += Y * Qabs * (lambda / con_h) * ulam * dl;
        }

        double Jpe = PI * a_eff * a_eff * sum;
        
        if ( Jpe < 0.0 )
            Jpe = 0.0;
        
        return Jpe;
    }

    double J_photodetachment( int Z, double a_eff )
    {
        if ( Z >= 0 )
            return 0.0;

        double Epdt_eV = E_pdt_eV(a_eff, Z);
        double l0      = lambda_min;
        double l1      = ( con_h * con_c ) / ( Epdt_eV * eV_to_J );

        if ( l1 > lambda_max ) l1 = lambda_max;
        if ( l1 <= l0 ) return 0.0;

        double dl  = ( l1 - l0 ) / static_cast<double>( N_lambda );
        double sum = 0.0;

        for ( int i = 0; i < N_lambda; ++i )
        {
            double lambda = l0 + ( (double( i ) + 0.5 ) * dl );
            double E_eV   = (con_h * con_c / lambda) * J_to_eV;
            double ulam   = u_lambda(lambda);
            double sigma  = sigma_pdt_m2(E_eV, a_eff, Z);

            sum += sigma * ( lambda / con_h ) * ulam;
        }

        if ( sum < 0.0 )
            sum = 0.0;
        
        return sum * dl;
    }*/
    
    
    // ======================= Currents struct =======================
    Currents total_currents(const CGridBasic * grid, const cell_basic & cell, int Z, double a_eff);

    // ======================= Hard Z-bounds =======================
    int compute_Zmin_auto(double a_m);
    
    int compute_Zmax_coulomb( double a_m, double Smax_Pa );
    
    
    
    void init(CDustComponent * dust, double _a_min, double _a_max, uint _Na, double _a0, double _a_sigma);
    


    // ======================= Build distribution =======================
       
    void build_distribution(const CGridBasic * grid, const cell_basic & cell, int ia, double & Zgr, double & Zs, double & Trot);
    
    void print_distribution();
    
    // emission functions
    double j_lambda_inter(const CGridBasic * grid, const cell_basic & cell, double lambda);
    
    //double j_lambda(double lambda);
        
    //double j_lambda_test(double lambda);
    
    
    /*void calc_FGpe(double & Fpe, double & Gpe, double a_eff, double Zgr)
    {
        double Jpe = J_photoelectric(Zgr, a_eff);
        
        Fpe=(con_m_e / con_m_p) * Jpe / (PIx2 *a_eff*a_eff*n_gas * sqrt( 2*con_kB*Tgas / (PI*con_m_p)));
        Gpe = 0;
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
    

    
    
private:    
    double W_eV;
    
    double eV_to_J;
    double J_to_eV;
    double m_to_AA;
    
    double * arr_a_eff;
    double * arr_dnda;
    
    double * arr_a_eff_large;
    double * arr_dnda_large;
    
    double Tgas;
    double Tdust;
    double mu_mol;
    
    double rho;
    
    double beta;
    
    double Zgas;
    
    double mu_sil;
    
    double n_el;
    double n_ion;
    double n_gas;
    double n_neu;
    
    double stick_e;
    double stick_ion;
    
    double alpha_gas;

    // Radiation field (test only): flat u_lambda between 613.6 eV
    double E_low_eV;
    double E_up_eV;
    
    // Integration
    dlist wavelength_list;
    double lambda_min;
    double lambda_max;
    CDustComponent * dust;
    
    // Charge bounds
    double Smax;
    
    uint mix_ID;

    // Recurrence window controls
    double tail_tol_ratio;
    int    Z_span_default;
    
    int tot_Zmin;
    int tot_Zmax;

    double a_min;
    double a_max;
    uint Na;
    uint Nla;
    
    double a0;
    double a_sigma;
    double m0;
    
    //helper functions
    double clamp_value( double x, double a, double b);
    void clear();
    int findIndex(double a) const;
    void gaussian_bounds(double mu, double sigma, double rel, double & x_min, double & x_max);
    double g1(double x);
    double g2(double x);
    double interpolate(double x, double x1, double x2, double y1, double y2);
    
    // ======================= Helper functions for recurrence =======================
    dlist compute_f_window(const CGridBasic * grid, const cell_basic & cell, int Zmin, int Zmax, int Zhard_min, int Zhard_max, double a_eff);

    bool tails_are_small(const dlist &f);

    int find_Zeq(const CGridBasic * grid, const cell_basic & cell, int Zmin, int Zmax, double a_eff);
    
    // Evaluate Gaussian at integer Z given parameters.
    double gaussian_value(double mu, double sigma, int Z);

    // Compute sum of squares error between f and model for given parameters (diagnostics).
    double gaussian_sse( const ilist &Z, const dlist &f,
                                double mu, double sigma);

    void fit_gaussian_pdf( const ilist    &Z,
                                         const dlist &f , double & mu, double & sigma);
    
    double get_u_lam(const CGridBasic * grid, const cell_basic & cell, uint iw);
    
    // excitation and damping coeff.
    void calc_FGIR(const CGridBasic * grid, const cell_basic & cell, double & FIR, double & GIR, double a_eff);
    void calc_FGn(double & Fn, double & Gn, double a_eff, double Zgr);
    void calc_FGpe(const CGridBasic * grid, const cell_basic & cell, double & Fpe, double & Gpe, double a_eff, int Zgr);
    void calc_FGi(double & Fi, double & Gi, double a_eff, double Zgr);
    void calc_FGp(double & Fp, double & Gp, double a_eff, double Zgr);
    //void compute_IR_coeffs(double & F_IR, double & G_IR, double a_eff);
    
    // grain charge
    double Gamma_photoelectric(const CGridBasic * grid, const cell_basic & cell, int Zgr, double a_eff);
    
    // ======================= DS87 reduced rates (SI) =======================
    double get_tau( double a, double T, double qabs );

    double get_nu( int Z, double q );

    double Jtilde_0( double tau );

    double Jtilde_neg( double tau, double nu );

    double Jtilde_pos( double tau, double nu );

    double Jtilde( double tau, double nu );
    
    //auxillary parameters
        
    double calc_tau_H(double a_eff);
    
    double calc_tau_ed(double a_eff, double Zgr);
    
    // ======================= WD01 band-yield factors (updated) =======================
    double Theta_eV( double E_eV, double a_m, int Z);

    double y0_bulk(double Theta);

    // Updated WD01: fixed lengths
    double y1_smallgrain( double a_m );

    double Elow_eV( double a_m, int Z );

    double Ehigh_eV( double E_eV, double a_m, int Z);

    double y2_escape( double E_eV, double a_m, int Z);
    
    
    // ======================= WDB06 yield: four-term structure =======================
    double Y_band( double E_eV, double a_m, int Z);

    // Placeholder for EUV/X-ray
    double Y_inner( double E_eV, double a_eff, int Z);

    // Placeholder for inner-shell Auger
    double Y_auger( double E_eV, double a_eff, int Z);

    // Placeholder for secondary electrons
    double Y_secondary( double /*E_eV*/, double /*a_m*/, int /*Z*/);

    //WD01 WDB06
    double get_Yield(double E_eV, double a_m, int Z);

    // ======================= Photodetachment cross section (WD01-like) =======================
    double sigma_pdt_m2( double E_eV, double a_m, int Z);

    // ======================= Collisional currents (DS87 in SI) =======================
    double vth_pref( double m );

    double J_electron( int Z, double a_eff);

    double J_ion_Hp( int Z, double a_eff );

    // ======================= Photoelectric emission & Photodetachment (lambda-integrals, SI) =======================
    // trapezoid (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);
    double J_photoelectric(const CGridBasic * grid, const cell_basic & cell, int Z, double a_eff);

    double J_photodetachment(const CGridBasic * grid, const cell_basic & cell, int Z, double a_eff );
    
    
    // emission helper
    
    double Power(double lambda, double a_eff, double Z);
    
    double fMW(double lambda, double Trot, double a_eff);
    
    double calc_mu_dipole(double Zgr, double a_eff);
};

#endif


