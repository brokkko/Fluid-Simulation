#ifndef FLUID_SIMULATION_CELL_H
#define FLUID_SIMULATION_CELL_H
#include <iostream>
#include<cmath>

struct PrimitiveVector
{
    double rho;
    double Vr;
    double Vph;
    double Vth;
    double Br;
    double Bph;
    double Bth;
    double P;

    PrimitiveVector& zeros();
    PrimitiveVector operator+(PrimitiveVector r) const;
    PrimitiveVector operator-(PrimitiveVector r) const;
    PrimitiveVector operator*(PrimitiveVector r) const;
    PrimitiveVector operator/(PrimitiveVector r) const;
    PrimitiveVector operator*(double r) const;
    PrimitiveVector operator/(double r) const;
    PrimitiveVector rotate(double phi,double theta);
    friend PrimitiveVector operator*(double l, PrimitiveVector r);
};

struct ConservativeVector
{
    double m;
    double Mr;
    double Mph;
    double Mth;
    double Br;
    double Bph;
    double Bth;
    double E;

    ConservativeVector& zeros();
    ConservativeVector operator+(ConservativeVector r) const;
    ConservativeVector operator-(ConservativeVector r) const;
    ConservativeVector operator*(ConservativeVector r) const;
    ConservativeVector operator/(ConservativeVector r) const;
    ConservativeVector operator*(double r) const;
    ConservativeVector operator/(double r) const;
    ConservativeVector rotate(double phi,double theta);
    friend ConservativeVector operator*(double l, ConservativeVector r);
};



struct Cell
{
    PrimitiveVector p;
    ConservativeVector c;

    double volume;
    double r;
    double phi;
    double theta;
    double Sr;
    double Sph;
    double Sth;

    Cell(double volume,double r,double phi,double theta,double Sr,double Sph,double Sth)
    {
        this->volume=volume;
        this->r=r;
        this->phi=phi;
        this->theta=theta;
        this->Sr = Sr;
        this->Sph = Sph;
        this->Sth = Sth;

        p.zeros();
        c.zeros();
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

