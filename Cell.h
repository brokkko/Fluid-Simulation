//
// Created by alex on 18.03.2022.
//

#ifndef FLUID_SIMULATION_CELL_H
#define FLUID_SIMULATION_CELL_H

#endif //FLUID_SIMULATION_CELL_H

struct Cell
{
    double rho;
    double vx;
    double vy;
    double B;

    Cell operator+(const Cell r) const
    {
        return Cell{ rho + r.rho ,vx + r.vx,vy + r.vy,B+r.B};
    }
    Cell operator-(const Cell r) const
    {
        return Cell{ rho - r.rho ,vx - r.vx,vy - r.vy,B-r.B};
    }
    Cell operator*(const Cell r) const
    {
        return Cell{ rho * r.rho ,vx * r.vx,vy * r.vy,B*r.B};
    }
    Cell operator/(const Cell r) const
    {
        return Cell{ rho / r.rho ,vx / r.vx,vy / r.vy,B/r.B};
    }
    Cell operator*(const double r) const
    {
        return Cell{ rho * r ,vx*r,vy*r,B*r};
    }
    Cell operator/(const double r) const
    {
        return Cell{ rho / r ,vx/r,vy/r,B/r};
    }
    friend Cell operator*(const double l, const Cell r)
    {
        return { l * r.rho ,l*r.vx,l*r.vy,l*r.B};
    }
};