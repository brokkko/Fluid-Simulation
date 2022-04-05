#ifndef FLUID_SIMULATION_GRID_H
#define FLUID_SIMULATION_GRID_H
#include "Cell.h"

struct Grid
{
    Cell** mesh;
    Cell** fluxMesh;
    unsigned int sizeX;
    unsigned int sizeY;
    Grid(unsigned int sizeX,unsigned int sizeY);
    void fluxMeshInit();
    void Fill(double v);
    ~Grid();
};
#endif //FLUID_SIMULATION_GRID_H
