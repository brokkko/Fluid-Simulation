#ifndef FLUID_SIMULATION_SIMULATION_H
#define FLUID_SIMULATION_SIMULATION_H
#include <cmath>
#include "Grid.h"
#include "SphericalGrid.h"
#include "fluxes.h"
#include "Derivative.h"
#include "Conditions.h"
#include "tvd.h"
#include <SFML/Graphics.hpp>
class Simulation {
    SphericalGrid &grid;
    SphericalGrid k = SphericalGrid::copyGrid(grid);
    int step;
    int order[9] = {0,1,2,1,2,0,2,0,1};
   //int order[9] = {0,1,2,1,2,0,1,2};
    SphericalGrid flux = SphericalGrid::copyGrid(grid);
    SphericalGrid grad = SphericalGrid::copyGrid(grid);

    double *densities;
    double *vels;
    double *temperature;
    double *magneticField;

public:
    Simulation(SphericalGrid &grid);
    void RKIntegrator(double &dt, double &t);
};
#endif //FLUID_SIMULATION_SIMULATION_H
