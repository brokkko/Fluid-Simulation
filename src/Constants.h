//
// Created by alex on 05.04.2022.
//

#ifndef FLUID_SIMULATION_CONSTANTS_H
#define FLUID_SIMULATION_CONSTANTS_H

#define ctg(x) std::cos(x)/std::sin(x)
#define clamp(x, y, z) std::max(y, std::min(x, z))
#define DT 120
#define CELL_SIZE 1
#define COLOR_SCHEME 2
#define mu 1.2566e-8
#define m_div_k 8249.805773
#define Ms 1.991e30
#define small_P 1e-10

#define small_rho 1e-10
#define G 0 //6.670e-11
#define gamma 5./3
#define A_SPEED 2000000
#define USE_CONST_A
//#define PRINT_NEG
#define ARROW_LEN_MULT 0.1f




#endif //FLUID_SIMULATION_CONSTANTS_H