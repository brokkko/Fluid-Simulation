#ifndef FLUID_SIMULATION_SIMULATION_H
#define FLUID_SIMULATION_SIMULATION_H
#include <cmath>
#include "Grid.h"


// one variable!!
// ospre flux limiter
Cell SlopeLim(Cell r);

double nonZeroDouble(double val);

Cell nonZeroDenom(Cell denom);

Cell F(Cell DX,Cell DY,Cell U);

void CalculateFlux(Grid& out, Grid& in);

Cell S(int x,int y,Cell val);

void ApplyBoundaryConditions(Grid& grid);

void RKIntegrator(Grid& grid, double dt);
#endif //FLUID_SIMULATION_SIMULATION_H
