#ifndef FLUID_SIMULATION_SIMULATION_H
#define FLUID_SIMULATION_SIMULATION_H
#include <cmath>
#include "Grid.h"
#include "SphericalGrid.h"


// one variable!!
// ospre flux limiter
Cell SlopeLim(Cell r);

double nonZeroDouble(double val);

Cell nonZeroDenom(Cell denom);

Cell F(Cell DX,Cell DY,Cell U);

void CalculateFlux(std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> out, Grid& in);

Cell S(int x,int y,Cell val);

void InitialConditions(SphericalGrid& grid);
void ApplyBoundaryConditions(SphericalGrid& grid);

void RKIntegrator(SphericalGrid& grid, double dt);
#endif //FLUID_SIMULATION_SIMULATION_H
