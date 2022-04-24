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

Cell F(Cell Dr,Cell Dtheta, Cell Dphi, Cell U, double r, double phi, double theta);

void CalculateFlux(std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> out, Grid& in);

Cell S(int x,int y,Cell val);

void InitialConditions(SphericalGrid& grid);
void ApplyBoundaryConditions(SphericalGrid& grid,double t);

void RKIntegrator(SphericalGrid& grid, double dt,double& t);
#endif //FLUID_SIMULATION_SIMULATION_H
