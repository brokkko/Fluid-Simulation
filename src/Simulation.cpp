#include <algorithm>
#include <iostream>
#include <tuple>
#include "Simulation.h"
#include "Constants.h"

double maxSpeed(Cell U){
    double p = gamma * (U.c_E - 0.5 * U.p_rho * (U.p_Vr * U.p_Vr + U.p_Vth * U.p_Vth + U.p_Vph * U.p_Vph)
                        - 0.5/mu * (U.p_Br * U.p_Br + U.p_Bph * U.p_Bph + U.p_Bth * U.p_Bth));
    double P = p + 0.5 / mu * (U.p_Br * U.p_Br + U.p_Bph * U.p_Bph + U.p_Bth * U.p_Bth);

    double B_2 = U.p_Br * U.p_Br + U.p_Bph * U.p_Bph + U.p_Bth * U.p_Bth;
    double B = std::sqrt(B_2);

    double cmax =
            std::sqrt(U.p_Vr * U.p_Vr + U.p_Vph * U.p_Vph + U.p_Vth * U.p_Vth)
            + 0.5*((gamma * p + B_2)/U.p_rho
            + std::sqrt(((gamma + B)/U.p_rho) * ((gamma + B) / U.p_rho) - 4 * (gamma * U.p_Br * U.p_Br) / (U.p_rho * U.p_rho)));

    return cmax*DT/CELL_SIZE;
}

Cell S(int x,int y,Cell val)
{
    if (x>3 && x<7 && y>32 && y< 58)
        return {0,0,0,0};
    else return {0,0,0,0};
}

void RKIntegrator(SphericalGrid& grid, double dt,double& t)
{
    t+=dt;

    double dr = grid.getRFromIndex(1)-grid.getRFromIndex(0);
    double dphi = (grid.getPhiFromIndex(1)-grid.getPhiFromIndex(0));
    double dtheta = grid.getThetaFromIndex(1)-grid.getThetaFromIndex(0);

    SphericalGrid fluxR = SphericalGrid::copyGrid(grid);
    SphericalGrid fluxPhi = SphericalGrid::copyGrid(grid);
    SphericalGrid fluxTheta = SphericalGrid::copyGrid(grid);

    SphericalGrid DR = SphericalGrid::copyGrid(grid);
    SphericalGrid DPhi = SphericalGrid::copyGrid(grid);
    SphericalGrid DTheta = SphericalGrid::copyGrid(grid);

    grid.UpdatePrim();

    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                DR.getCellRef(r,phi,theta) = (grid.getCell(r+1,phi,theta)-grid.getCell(r-1,phi,theta))/(2*1/*dr*/);
                DPhi.getCellRef(r,phi,theta) = (grid.getCell(r,phi+1,theta)-grid.getCell(r,phi-1,theta))/(2*1/*dphi * grid.getRFromIndex(r)*/);
                DTheta.getCellRef(r,phi,theta) = (grid.getCell(r,phi,theta+1)-grid.getCell(r,phi,theta-1))/(2*1/*dtheta * grid.getRFromIndex(r)*/);
            }
        }
    }

    std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> flux(fluxR,fluxPhi,fluxTheta);



    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                Cell& c=grid.getCellRef(r,phi,theta);
                if (r==30 && phi == 50)
                {
                    int a = 0;
                }
                auto Dr = DR.getCell(r, phi,theta);
                auto Dtheta = DTheta.getCell(r, phi, theta);
                auto Dphi = DPhi.getCell(r, phi, theta);
                c = c -0.5 * dt * F(Dr,Dtheta,Dphi,c);

            }
        }
    }

    grid.UpdateCons();

    SphericalGrid k = SphericalGrid::copyGrid(grid);
    CalculateFlux(flux, grid);

    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                Cell& c=grid.getCellRef(r,phi,theta);
                c = c - dt * (   (fluxR.getCell(r+1,phi,theta)     - fluxR.getCell(r,phi,theta)) /1/*dr*/     +
                        (fluxPhi.getCell(r,phi+1,theta)   - fluxPhi.getCell(r,phi,theta))/(1/*dphi * grid.getRFromIndex(r)*/)  +
                        (fluxTheta.getCell(r,phi,theta+1) - fluxTheta.getCell(r,phi,theta))/(1/*dtheta * grid.getRFromIndex(r)*/));
            }
        }
    }


}
