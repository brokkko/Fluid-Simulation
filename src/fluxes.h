//
// Created by alex on 08.05.2022.
//

#ifndef FLUID_SIMULATION_FLUXES_H
#define FLUID_SIMULATION_FLUXES_H
#include <tuple>
#include "SphericalGrid.h"
#include "Cell.h"
#include "tvd.h"


Cell FluxR(Cell U);
Cell FluxPhi(Cell U);
Cell FluxTheta(Cell U);

void CalculateFlux(std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> out, SphericalGrid& in, std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> grad);
void CalculateFluxR(SphericalGrid& out, SphericalGrid& in, SphericalGrid& grad );
void CalculateFluxPh(SphericalGrid& out, SphericalGrid& in, SphericalGrid& grad );
void CalculateFluxTh(SphericalGrid& out, SphericalGrid& in, SphericalGrid& grad );
#endif //FLUID_SIMULATION_FLUXES_H
