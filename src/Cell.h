#ifndef FLUID_SIMULATION_CELL_H
#define FLUID_SIMULATION_CELL_H
#include <iostream>

struct Cell
{
    double rho;
    double vx;
    double vy;
    double vz;
    double Bx;
    double By;
    double Bz;
    double E;
    static Cell zeros();
    Cell operator+(Cell r) const;
    Cell operator-(Cell r) const;
    Cell operator*(Cell r) const;
    Cell operator/(Cell r) const;
    Cell operator*(double r) const;
    Cell operator/(double r) const;
    friend Cell operator*(double l, Cell r);
};



#endif //FLUID_SIMULATION_CELL_H

