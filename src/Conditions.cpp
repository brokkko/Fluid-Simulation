//
// Created by alex on 08.05.2022.
//
#include "Conditions.h"
#include "Constants.h"

void InitialConditions(SphericalGrid& grid) {
    for (int x = 0; x < grid.getSizeR(); x++) {
        for (int y = 0; y < grid.getSizePhi(); y++) {
            double rho=small_rho;
            double vx=0;
            double vy=0;
            double vz=0;
            double Bx=0.000;
            double By=0.000;
            double Bz=0.000;
            double T =10;
            // for (int y = 45; y < 55; y++){
            //        for (int x = 20; x < 40; x++){
            if(x>30 && x< 60 && y>45 && y<55)
            {
                rho=small_rho*20;
                //vx=100000;
                T=10;
                vy=-100000;
            }


            //nk=p_rho/m m=1.6733e-27
            //k=1.38044e-23
            double E = 2 * rho * m_div_k * T/(gamma-1)
                       + rho * (vx*vx + vy*vy + vz*vz)
                       + (Bx*Bx + By*By + Bz*Bz) /(2*mu);
            double p = (gamma - 1) * (E - 0.5 * rho * (vx*vx + vy*vy + vz*vz) - 0.5 / mu * (Bx*Bx + By*By + Bz*Bz));
            Cell& c =grid.getCellRef(x,y,0);
            c.p_rho = rho;
            c.p_Vr = vx;
            c.p_Vph = vy;
            c.p_Vth = vz;
            c.p_Br = Bx;
            c.p_Bph = By;
            c.p_Bth = Bz;
            c.p_P = p;
        }
    }
    grid.UpdateCons();
}


void ApplyBoundaryConditions(SphericalGrid& grid,double t,double* dens,double* vels)
{
    for (int x=0;x<grid.getSizePhi();x++) {
        double vx=0;
        //double vx=vels[x];
        double vy=0;
        double vz=0;
        double Bx=0.000;
        double By=0.000;
        double Bz=0.000;
        double T = 500000;
        double rho=small_rho;
        // double rho=dens[x];
        if( x>5 && x< 25)
        {
            // T=500000;
            // rho=small_rho*100;
            // vx=1000000;
        }
        double E = 2 * rho * m_div_k * T/(gamma-1)
                   + rho * (vx*vx + vy*vy + vz*vz)
                   + (Bx*Bx + By*By + Bz*Bz) /(2*mu);
        double p = (gamma - 1) * (E - 0.5 * rho * (vx*vx + vy*vy + vz*vz) - 0.5 / mu * (Bx*Bx + By*By + Bz*Bz));
        Cell& c =grid.getCellRef(1,x,0);
        c.p_rho = rho;
        c.p_Vr = vx;
        c.p_Vph = vy;
        c.p_Vth = vz;
        c.p_Br = Bx;
        c.p_Bph = By;
        c.p_Bth = Bz;
        c.p_P = p;
    }
    grid.UpdateCons();
}