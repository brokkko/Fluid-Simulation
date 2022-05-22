#ifndef FLUID_SIMULATION_CONSTANTS_H
#define FLUID_SIMULATION_CONSTANTS_H

#define ctg(x) std::cos(x)/std::sin(x)
#define clamp(x, y, z) std::max(y, std::min(x, z))
#define DT 100000
#define CFL 0.4
#define DATA "../data/bnd.nc"
#define SIZE_R 90
#define SIZE_PH 180
#define SIZE_TH 21
#define MIN_RADIUS 1.497131e10
#define MAX_RADIUS 2.28e11
#define COLOR_SCHEME 1
#define mu 1.2566e-6
#define m_div_k 8249.805773
#define Ms 1.991e30
#define SOLAR_ROTATION 2114208
#define small_P 1e-20

#define small_rho 1e-20
#define G 6.670e-11
#define gamma 5./3
#define A_SPEED 1000000
//#define USE_CONST_A
//#define PRINT_NEG
#define ARROW_LEN_MULT 0.1f

#define T_R 0
#define T_PHI 1
#define T_THETA 2

#endif //FLUID_SIMULATION_CONSTANTS_H