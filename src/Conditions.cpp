//
// Created by alex on 08.05.2022.
//
#include "Conditions.h"
#include "Constants.h"

void InitialConditions(SphericalGrid& grid) {
    for (int x = 0; x < grid.getSizeR(); x++) {
        for (int y = 0; y < grid.getSizePhi(); y++) {
            double rho=small_rho/(grid.getRFromIndex(x)*grid.getRFromIndex(x))*(grid.getRFromIndex(0)*grid.getRFromIndex(0));
            double vx=0;
            double vy=0;
            double vz=0;
            double Bx=0.0;
            double By=0.0;
            double Bz=0.0;
            double T =10;
            // for (int y = 45; y < 55; y++){
            //        for (int x = 20; x < 40; x++){
            if(x>10 && x< 30 && y>45 && y<55)
            {
                rho=small_rho*20;
                //vx=1000000;
                T=10;
                vy=1000000;
            }



            //nk=p_rho/m m=1.6733e-27
            //k=1.38044e-23
            double E = 2 * rho * m_div_k * T/(gamma-1)
                       + rho * (vx*vx + vy*vy + vz*vz)
                       + (Bx*Bx + By*By + Bz*Bz) /(2);
            double p = (gamma - 1) * (E - 0.5 * rho * (vx*vx + vy*vy + vz*vz) - 0.5 * (Bx*Bx + By*By + Bz*Bz));
            Cell& c =grid.getCellRef(x,y,0);
            c.p.rho = rho;
            c.p.Vr = vx;
            c.p.Vph = vy;
            c.p.Vth = vz;
            c.p.Br = Bx;
            c.p.Bph = By;
            c.p.Bth = Bz;
            c.p.P = p;

        }
    }
    grid.UpdateCons();
}


void ApplyBoundaryConditions(SphericalGrid& grid,double t,double* dens,double* vels)
{
    for (int x=0;x<grid.getSizePhi();x++) {
        double vx=0;
       // double vx=vels[x];
        double vy=0;
        double vz=0;
        double Bx=0.000;
        double By=0.000;
        double Bz=0.000;
        double T = 10;
        double rho=small_rho;
        // double rho=dens[x];
        if( x>5 && x< 15)
        {
            // T=10;
             //rho=small_rho*1000;
            // vx=1000000;
        }
        double E = 2 * rho * m_div_k * T/(gamma-1)
                   + rho * (vx*vx + vy*vy + vz*vz)
                   + (Bx*Bx + By*By + Bz*Bz) /2;
        double p = (gamma - 1) * (E - 0.5 * rho * (vx*vx + vy*vy + vz*vz) - 0.5 * (Bx*Bx + By*By + Bz*Bz));
        Cell& c =grid.getCellRef(0,x,0);
        c.p.rho = rho;
        c.p.Vr = vx;
        c.p.Vph = vy;
        c.p.Vth = vz;
        c.p.Br = Bx;
        c.p.Bph = By;
        c.p.Bth = Bz;
        c.p.P = p;
    }
    grid.UpdateCons();
}