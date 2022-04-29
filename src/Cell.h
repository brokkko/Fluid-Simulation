#ifndef FLUID_SIMULATION_CELL_H
#define FLUID_SIMULATION_CELL_H
#include <iostream>

struct Cell
{
    double p_rho;
    double p_Vr;
    double p_Vph;
    double p_Vth;
    double p_Br;
    double p_Bph;
    double p_Bth;
    double p_P;

    double c_rho;
    double c_Mr;
    double c_Mph;
    double c_Mth;
    double c_Br;
    double c_Bph;
    double c_Bth;
    double c_E;


    double volume;
    double r;
    double phi;
    double theta;


    Cell(double volume,double r,double phi,double theta)
    {
        this->volume=volume;
        this->r=r;
        this->phi=phi;
        this->theta=theta;
        auto this_p = reinterpret_cast<double*>(this);
        for(int i=0;i<16;i++)
            this_p[i]=0.0;
    }
    Cell() { };

    void UpdatePrim();
    void UpdateCons();


    Cell& zeros();
    Cell operator+(Cell r) const;
    Cell operator-(Cell r) const;
    Cell operator*(Cell r) const;
    Cell operator/(Cell r) const;
    Cell operator*(double r) const;
    Cell operator/(double r) const;
    friend Cell operator*(double l, Cell r);
};



#endif //FLUID_SIMULATION_CELL_H

