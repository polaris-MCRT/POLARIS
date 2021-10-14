#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <time.h>
#include <vector>
using namespace std;

// Header and Version of POLARIS
#define PROG_ID "POLARIS    V4.06.00        "

// Flags to activate WINDOWS support, some DEBUG messages, and CAMPS BENCHMARK
//#define DEBUG
//#define WINDOWS
//#define CAMPS_BENCHMARK

// Constants taken from astropy/numpy (reference: CODATA 2014, IAU 2012 Resolution B2)
#define PI 3.1415926535897932384626433832795028841971 // PI
#define PIsq sqrt(PI)                                 // sqrt(PI)
#define PI2 (PI / 2.0)                                // PI / 2
#define PI4 (PI / 4.0)                                // PI / 4
#define invPI2 (2.0 / PI)                             // 2 / PI
#define PIx2 (PI * 2.0)                               // 2 * PI
#define PIx4 (PI * 4.0)                               // 4 * PI
#define PI4x3 (PI * 3.0 / 4.0)                        // 3 * PI / 4
#define TWOTHIRD 0.6666666666666666666666666666666    // 2 / 3
#define m_H 1.66053904e-27                            // Atomic mass unit [kg]
#define con_h 6.62607004e-34                          // Planck constant [Js]
#define con_hq 1.0545718e-34                          // Reduced Planck constant [Js]
#define con_kB 1.38064852e-23                         // Boltzmann constant [J / K]
#define con_c 299792458.0                             // Speed of light in vacuum [m / s]
#define PIx4_c (PIx4 / con_c)                         // (4 * PI / c)
#define con_e 1.6021766208e-19                        // Electron charge [C]
#define con_Ryd 10973731.568508                       // Rydberg constant [1 / m]
#define con_AU 149597870700.0                         // Astronomical Unit [m]
#define con_pc 3.0856775814671916e+16                 // Parsec [m]
#define con_ly 9460730472580800.0                     // Lightyear [m]
#define con_G 6.67408e-11                             // Gravitational constant [m^3 / (kg * s^2)]
#define con_sigma 5.670367e-08                        // Stefan-Boltzmann constant [W / (m^2 * K^4)]
#define L_sun 3.828e+26                               // Nominal solar luminosity [W]
#define M_sun 1.9884754153381438e+30                  // Solar mass [kg]
#define R_sun 695700000.0                             // Nominal solar radius [m]
#define con_Na 6.02214129e23                          // Avogadro's number [1 / mol]
#define con_mb 9.274009994e-24                        // Bohr magneton [J / T]
#define con_r_bohr 5.2917721067e-11                   // Bohr radius [m]
#define con_m_e 9.10938356e-31                        // Electron mass [kg]
#define con_m_p 1.672621898e-27                       // Proton mass [kg]
#define con_epsilon_0 8.854187817620389e-12           // Vacuum permittivity [F / m]
#define con_eps (con_h * con_c / PIx4)                // (h * c / (4 * PI))

#ifdef CAMPS_BENCHMARK
// Part to perform Camps et. al (2015) benchmark (adjust definition below, if
// necessary!!!)
#define WL_MIN 1e-9
#define WL_MAX 1e-2
#define WL_STEPS 1201
#else
// Default parameters of the global wavelength grid
#define WL_MIN 0.1e-6
#define WL_MAX 2000.0e-6
#define WL_STEPS 100
#endif

// Parameter for numerical limitations
#define MAX_LVG_ITERATIONS 200
#define MAX_INTERACTION 1e6
#define MAX_RT_RAYS 1e7
#define MIN_LEN_STEP 1e4
#define ACC_SELECT_LEVEL 1.0e-6
#define DIFF_GAMMA 7.00
#define PERCENTAGE_STEP 0.1

// Limits of the Runge-Kutta-Fehlberg raytracing method
#define rel_err 1.0e-6
#define abs_err 1.0e-30
#define MAX_SOLVER_STEPS 1500000

// Define the fits file extension
// ".fits" normal fits file
// ".fits.gz" compressed fits file
#define FITS_COMPRESS_EXT ".fits.gz"

// Number of quantities in dust component file (C_{x,abs}, C_{x,ext}, C_{x,sca}, ...)
#define NR_OF_EFF 8

// Number of dust grain size distribution parameters
#define NR_OF_SIZE_DIST_PARAM 14

// Number of entries for different detectors
#define NR_OF_MC_DET 10
#define NR_OF_RAY_DET 14
#define NR_OF_LINE_DET 17

// Number of entries for different sources
#define NR_OF_POINT_SOURCES 8
#define NR_OF_DIFF_SOURCES 9
#define NR_OF_LASER_SOURCES 12
#define NR_OF_BG_SOURCES 8

#define TEMP_MIN 2.728
#define TEMP_MAX 3000
#define TEMP_STEP 1000

// detector ids
#define DET_PLANE 0
#define DET_POLAR 1
#define DET_SPHER 2
#define DET_SLICE 3

// timestep and offset for time-dependent scattering and dust emission
#define SCA_DT 0
#define RAY_DT 0

// phase functions
#define PH_ISO 0
#define PH_HG 1
#define PH_MIE 2

// Alignment mechanisms
#define ALIG_RND 0
#define ALIG_INTERNAL 1
#define ALIG_PA 2
#define ALIG_IDG 4
#define ALIG_RAT 8
#define ALIG_GOLD 16
#define ALIG_KRAT 32

#define SUPERTHERMAL_LIMIT 3
#define MACH_LIMIT 1
#define MRW_LIMIT 7
#define PDA_LIMIT 30

#define CMD_TEMP 0
#define CMD_DUST_EMISSION 1
#define CMD_DUST_SCATTERING 2
#define CMD_PROBING 3
#define CMD_RAT 4
#define CMD_TEMP_RAT 5
#define CMD_LINE_EMISSION 6
#define CMD_FORCE 7
#define CMD_OPIATE 8
#define CMD_SYNCHROTRON 9
#define CMD_DUST_TIME 10
#define CMD_TEMP_POLY 11

// PDA IDs
#define PDA_TEMP 0
#define PDA_RAT 1
#define PDA_PRES 2

// grid data IDs
#define GRIDgas_dens 0
#define GRIDdust_dens 1
#define GRIDdust_temp 2
#define GRIDgas_temp 3
#define GRIDmx 4
#define GRIDmy 5
#define GRIDmz 6
#define GRIDvx 7
#define GRIDvy 8
#define GRIDvz 9
#define GRIDpx 10
#define GRIDpy 11
#define GRIDpz 12
#define GRIDa_alg 13
#define GRIDa_min 14
#define GRIDa_max 15
#define GRIDq 16
#define GRIDratio 17
#define GRIDv_turb 18
#define GRIDPDA 19
#define GRIDopiate 20
#define GRIDdust_id 21

#define GRIDn_th 22  // number density of thermal electrons
#define GRIDT_e 23   // Temperature of thermal electrons
#define GRIDn_cr 24  // number density of CR electrons
#define GRIDg_min 25 // gamma min for power law distribution
#define GRIDg_max 26 // gamma max for power law distribution
#define GRIDp 27     // power law exponent

#define GRIDgas_mdens 28
#define GRIDdust_mdens 29
#define GRIDradx 30
#define GRIDrady 31
#define GRIDradz 32
#define GRIDrad 33

#define GRIDavg_th 34
#define GRIDavg_dir 35

#define minGRID GRIDgas_dens
#define maxGRID GRIDavg_dir

#define MAX_UINT uint(-1)
#define MAX_DOUBLE double(uint(-1))

#define EPS_DOUBLE numeric_limits<double>::epsilon()
#define EPS_FLOAT numeric_limits<float>::epsilon()

// grid types
#define GRID_ID_OCT 20
#define GRID_ID_SPH 30
#define GRID_ID_CYL 40
#define GRID_ID_VOR 50

#define POP_LTE 1
#define POP_FEP 2
#define POP_LVG 3

#define COL_H2_FULL 1
#define COL_H2_PARA 2
#define COL_H2_ORTH 3
#define COL_HE_FULL 6

// Dust temperature cases
#define TEMP_EMPTY 0
#define TEMP_SINGLE 1
#define TEMP_EFF 2
#define TEMP_STOCH 3
#define TEMP_FULL 4

// Cross-sections IDs
#define CROSS_ABS 1
#define CROSS_SCA 2

// Extrapolation/Interpolation IDs
#define CONST 0
#define LINEAR 1
#define SPLINE 2

// Radiation types on detector
// direct starlight
#define DIRECT_STAR 0
// all scattering at dust grains
#define SCATTERED_DUST 1
// single scattering at dust grains
#define SCATTERED_DUST_1 2
// multiple scattering at dust grains
#define SCATTERED_DUST_2 3

// Type of emission
#define RESULTS_RAY 0
#define RESULTS_MC 1
#define RESULTS_FULL 2

// Type of calorimetry data
#define CALO_HEAT_CAP 0
#define CALO_ENTHALPY 1

// Type of healpix orientation
#define HEALPIX_FIXED 0
#define HEALPIX_YAXIS 1
#define HEALPIX_CENTER 2

// Mie-scattering calculation
#define MIE_SIZE_STEPS 100
// Number of angles for scattering between 0° and 90°
#define NANG 91
#define MAX_MIE_ITERATIONS 1000000
#define MIN_MIE_SIZE_PARAM 1e-6
#define MIE_ACCURACY 1e-20

// Projections for midplane files
#define PROJ_XY 1
#define PROJ_XZ 2
#define PROJ_YZ 3

#define SEP_LINE                                                                                             \
    "**********************************************************************************"                     \
    "***\n"
#define CLR_LINE                                                                                             \
    "                                                                                  "                     \
    "   \r"
#define TAB "\t"

#ifdef WINDOWS
#include <direct.h>
#include <io.h>
#include <string>
#include <windows.h>
#define SEP "\\"
#define LINE_DELAY 5
#pragma warning(disable : 4566)
#else
#include <cmath>
#include <stdlib.h>
#include <string.h>
//	#include <asm/io.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#define SEP "/"
#define LINE_DELAY 5
#endif

// data types
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned long long ullong;
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef long long llong;
typedef vector<double> dlist;
typedef map<uint, vector<double> > maplist;
typedef vector<uchar> clist;
typedef vector<uint> uilist;
typedef vector<ushort> uslist;
typedef vector<int> ilist;
typedef vector<string> strlist;
typedef complex<float> fcomplex;
typedef complex<double> dcomplex;
