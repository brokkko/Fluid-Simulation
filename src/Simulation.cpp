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

void CalculateGradients(std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> grad,SphericalGrid& grid)
{
    double dr = grid.getRFromIndex(1)-grid.getRFromIndex(0);
    double dphi = (grid.getPhiFromIndex(1)-grid.getPhiFromIndex(0));
    double dtheta = grid.getThetaFromIndex(1)-grid.getThetaFromIndex(0);

    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                auto r1 =     grid.getCell(r+1,phi,theta);
                auto r_1 =    grid.getCell(r-1,phi,theta);
                auto phi1 =   grid.getCell(r,phi+1,theta);
                auto phi_1 =  grid.getCell(r,phi-1,theta);
                auto theta1 = grid.getCell(r,phi,theta+1);
                auto theta_1 =grid.getCell(r,phi,theta-1);
                auto curr =   grid.getCell(r,phi,theta);

                auto Dr=      (r1     - r_1)     / (2*dr);
                auto Dphi =   (phi1   - phi_1)   / (2*dphi * grid.getRFromIndex(r));
                auto Dtheta = (theta1 - theta_1) / (2*dtheta * grid.getRFromIndex(r));

                auto rR =     (curr - r_1)     / nonZeroDenom(r1     - curr);
                auto rPhi =   (curr - phi_1)   / nonZeroDenom(phi1   - curr);
                auto rTheta = (curr - theta_1) / nonZeroDenom(theta1 - curr);
                //TODO: Check for mass loss
                std::get<T_R>(grad).getCellRef(r,phi,theta)     = SlopeLim(rR)     * Dr;
                std::get<T_PHI>(grad).getCellRef(r,phi,theta)   = SlopeLim(rPhi)   * Dphi;
                std::get<T_THETA>(grad).getCellRef(r,phi,theta) = SlopeLim(rTheta) * Dtheta;
            }
        }
    }
}

void PredictorStep(std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> grad,SphericalGrid& grid,double dt)
{
    SphericalGrid& DR =std::get<T_R>(grad);
    SphericalGrid& DPhi =std::get<T_PHI>(grad);
    SphericalGrid& DTheta =std::get<T_THETA>(grad);
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
}

void ApplyFluxes(std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> flux,SphericalGrid& grid,double dt)
{
    double dr = grid.getRFromIndex(1)-grid.getRFromIndex(0);
    double dphi = (grid.getPhiFromIndex(1)-grid.getPhiFromIndex(0));
    double dtheta = grid.getThetaFromIndex(1)-grid.getThetaFromIndex(0);

    auto& fluxR = std::get<T_R>(flux);
    auto& fluxPhi = std::get<T_PHI>(flux);
    auto& fluxTheta = std::get<T_THETA>(flux);
    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                Cell& c=grid.getCellRef(r,phi,theta);
                c = c - dt * (
                         (fluxR.getCell(r+1,phi,theta)     - fluxR.getCell(r,phi,theta))     /  dr
                        +(fluxPhi.getCell(r,phi+1,theta)   - fluxPhi.getCell(r,phi,theta))   / (dphi * grid.getRFromIndex(r))
                        +(fluxTheta.getCell(r,phi,theta+1) - fluxTheta.getCell(r,phi,theta)) / (dtheta * grid.getRFromIndex(r))
                );
            }
        }
    }

}


void RKIntegrator(SphericalGrid& grid, double dt,double& t)
{
    t+=dt;



    SphericalGrid fluxR = SphericalGrid::copyGrid(grid);
    SphericalGrid fluxPhi = SphericalGrid::copyGrid(grid);
    SphericalGrid fluxTheta = SphericalGrid::copyGrid(grid);
    std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> flux(fluxR,fluxPhi,fluxTheta);


    SphericalGrid DR = SphericalGrid::copyGrid(grid);
    SphericalGrid DPhi = SphericalGrid::copyGrid(grid);
    SphericalGrid DTheta = SphericalGrid::copyGrid(grid);
    std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> grad(DR,DPhi,DTheta);

    grid.UpdatePrim();
    CalculateGradients(grad,grid);
    PredictorStep(grad,grid,dt);
    grid.UpdateCons();
    CalculateFlux(flux, grid, grad);
    ApplyFluxes(flux,grid,dt);

}
