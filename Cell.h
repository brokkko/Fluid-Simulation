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
    double vz;
    double Bx;
    double By;
    double Bz;
    double E;
    static Cell zeros()
    {
        return Cell{ 0,0,0,0,0,0,0,0};
    }
    Cell operator+(const Cell r) const
    {
        return Cell{ rho + r.rho ,
                     vx + r.vx,
                     vy + r.vy,
                     vz + r.vy,
                     Bx+r.Bx,
                     By+r.By,
                     Bz+r.Bz,
                     E+r.E,};
    }
    Cell operator-(const Cell r) const
    {
        return Cell{ rho - r.rho ,
                     vx - r.vx,
                     vy - r.vy,
                     vz - r.vy,
                     Bx - r.Bx,
                     By - r.By,
                     Bz - r.Bz,
                     E - r.E,};
    }
    Cell operator*(const Cell r) const
    {
        return Cell{ rho * r.rho ,
                     vx * r.vx,
                     vy * r.vy,
                     vz * r.vy,
                     Bx * r.Bx,
                     By * r.By,
                     Bz * r.Bz,
                     E * r.E,};
    }
    Cell operator/(const Cell r) const
    {
        return Cell{ rho / r.rho ,
                     vx / r.vx,
                     vy / r.vy,
                     vz / r.vy,
                     Bx / r.Bx,
                     By / r.By,
                     Bz / r.Bz,
                     E / r.E,};
    }
    Cell operator*(const double r) const
    {
        return Cell{ rho * r ,
                     vx * r,
                     vy * r,
                     vz * r,
                     Bx * r,
                     By * r,
                     Bz * r,
                     E * r};
    }
    Cell operator/(const double r) const
    {
        return Cell{ rho / r ,
                     vx / r,
                     vy / r,
                     vz / r,
                     Bx / r,
                     By / r,
                     Bz / r,
                     E / r};
    }
    friend Cell operator*(const double l, const Cell r)
    {
        return Cell{ l*r.rho ,
                     l*r.vx ,
                     l*r.vy ,
                     l*r.vz ,
                     l*r.Bx ,
                     l*r.By ,
                     l*r.Bz ,
                     l*r.E};
    }
};