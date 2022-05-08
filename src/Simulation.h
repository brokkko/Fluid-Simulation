#ifndef FLUID_SIMULATION_SIMULATION_H
#define FLUID_SIMULATION_SIMULATION_H
#include <cmath>
#include "Grid.h"
#include "SphericalGrid.h"
#include "fluxes.h"
#include "Derivative.h"
#include "Conditions.h"
#include "tvd.h"

double nonZeroDouble(double val);

Cell nonZeroDenom(Cell denom);


Cell S(int x,int y,Cell val);

void RKIntegrator(SphericalGrid& grid, double dt,double& t);
#endif //FLUID_SIMULATION_SIMULATION_H
