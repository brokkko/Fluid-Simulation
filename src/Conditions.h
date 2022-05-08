//
// Created by alex on 08.05.2022.
//

#ifndef FLUID_SIMULATION_CONDITIONS_H
#define FLUID_SIMULATION_CONDITIONS_H

#include "SphericalGrid.h"

void InitialConditions(SphericalGrid& grid);
void ApplyBoundaryConditions(SphericalGrid& grid,double t,double* dens,double* vels);

#endif //FLUID_SIMULATION_CONDITIONS_H
