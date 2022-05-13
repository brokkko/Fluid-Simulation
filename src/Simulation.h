#ifndef FLUID_SIMULATION_SIMULATION_H
#define FLUID_SIMULATION_SIMULATION_H
#include <cmath>
#include "Grid.h"
#include "SphericalGrid.h"
#include "fluxes.h"
#include "Derivative.h"
#include "Conditions.h"
#include "tvd.h"

class Simulation {
    SphericalGrid &grid;
    SphericalGrid k = SphericalGrid::copyGrid(grid);

    SphericalGrid fluxR = SphericalGrid::copyGrid(grid);
    SphericalGrid fluxPhi = SphericalGrid::copyGrid(grid);
    SphericalGrid fluxTheta = SphericalGrid::copyGrid(grid);

    SphericalGrid DR = SphericalGrid::copyGrid(grid);
    SphericalGrid DPhi = SphericalGrid::copyGrid(grid);
    SphericalGrid DTheta = SphericalGrid::copyGrid(grid);
public:
    Simulation(SphericalGrid &grid);
    void RKIntegrator(double &dt, double &t);
};
#endif //FLUID_SIMULATION_SIMULATION_H
