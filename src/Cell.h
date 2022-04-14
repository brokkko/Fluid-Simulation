#ifndef FLUID_SIMULATION_CELL_H
#define FLUID_SIMULATION_CELL_H
#include <iostream>

struct Cell
{
    double rho;
    double rhoVr;
    double rhoVphi;
    double rhoVtheta;
    double Br;
    double Bphi;
    double Btheta;
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

