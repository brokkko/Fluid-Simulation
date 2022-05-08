//
// Created by alex on 08.05.2022.
//

#ifndef FLUID_SIMULATION_TVD_H
#define FLUID_SIMULATION_TVD_H
#include "Cell.h"
double nonZeroDouble(double val);
Cell nonZeroDenom(Cell denom);
Cell SlopeLim(Cell r);
#endif //FLUID_SIMULATION_TVD_H

