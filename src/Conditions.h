//
// Created by alex on 08.05.2022.
//

#ifndef FLUID_SIMULATION_CONDITIONS_H
#define FLUID_SIMULATION_CONDITIONS_H

#include "SphericalGrid.h"
#include "ReadData.h"

void InitialConditions(SphericalGrid& grid, double* dens,double* vels, double *temperature, double *magneticField );
void ApplyBoundaryConditions(SphericalGrid& grid,double t,double* dens,double* vels, double *temperature, double *magneticField);

double* getDensity();
double* getVelocity();
double* getTemperature();
double* getMagneticField();

#endif //FLUID_SIMULATION_CONDITIONS_H
