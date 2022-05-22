#include <algorithm>
#include <iostream>
#include <tuple>
#include "Simulation.h"
#include "Constants.h"

double CalculateDT(SphericalGrid& grid){
    double minDT = DT;

    double dr =     grid.getRFromIndex(1)       - grid.getRFromIndex(0);
    double dphi =  (grid.getPhiFromIndex(1)     - grid.getPhiFromIndex(0));
    double dtheta = grid.getThetaFromIndex(1) - grid.getThetaFromIndex(0);

    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {

                auto curr =   grid.getCell(r,phi,theta);

                double V_r = curr.p.Vr*curr.p.Vr + curr.p.Vth*curr.p.Vth + curr.p.Vph*curr.p.Vph;
                double minDr = CFL *  (dr / (std::sqrt(gamma*curr.p.P/curr.p.rho) + std::abs(curr.p.Vr)));
                double minDph = CFL *  (dphi*curr.r / (std::sqrt(gamma*curr.p.P/curr.p.rho) + std::abs(curr.p.Vph)));
                double minDth = CFL *  (dtheta*curr.r / (std::sqrt(gamma*curr.p.P/curr.p.rho) + std::abs(curr.p.Vth)));
                minDT = std::min(minDr, std::min(minDph, std::min(minDth, minDT)));
            }
        }
    }
    return minDT;
}


/*
double CalculateDT(SphericalGrid &grid) {
    double minDT = DT;

    double dr = grid.getRFromIndex(1) - grid.getRFromIndex(0);
    double dphi = (grid.getPhiFromIndex(1) - grid.getPhiFromIndex(0));
    double dtheta = grid.getThetaFromIndex(1) - grid.getThetaFromIndex(0);

    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {

                auto curr = grid.getCell(r, phi, theta);
                double gpr=gamma * curr.p.P;
                if(std::isnan(curr.p.P)){
                    //std::cout << "HERE" << std::endl;
                }
                double Bmag2  = curr.p.Br*curr.p.Br + curr.p.Bph*curr.p.Bph+curr.p.Bth*curr.p.Bth;
                double Cf = gpr - Bmag2;
                double Cf2 = gamma * curr.p.P + Bmag2;
                double Btr =sqrt(Cf*Cf + 4.0*gpr*(Bmag2-curr.p.Br*curr.p.Br));
                double Btph =sqrt(Cf*Cf + 4.0*gpr*(Bmag2-curr.p.Bph*curr.p.Bph));
                double Btth =sqrt(Cf*Cf + 4.0*gpr*(Bmag2-curr.p.Bth*curr.p.Bth));
                double minDr = CFL * (dr / (sqrt(0.5*(Cf2+Btr)/curr.p.rho) + std::abs(curr.p.Vr)));
                double minDph =
                        CFL * (dphi * curr.r / (sqrt(0.5*(Cf2+Btph)/curr.p.rho) + std::abs(curr.p.Vph)));
                double minDth =
                        CFL * (dtheta * curr.r / (sqrt(0.5*(Cf2+Btth)/curr.p.rho) + std::abs(curr.p.Vth)));
                minDT = std::min(minDr, std::min(minDph, std::min(minDth, minDT)));
            }
        }
    }
    return minDT;
}*/

void CalculateGradients(SphericalGrid & grad, SphericalGrid &grid,int dir) {
    double dr = grid.getRFromIndex(1) - grid.getRFromIndex(0);
    double dphi = (grid.getPhiFromIndex(1) - grid.getPhiFromIndex(0));
    double dtheta = grid.getThetaFromIndex(1) - grid.getThetaFromIndex(0);

    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                auto curr = grid.getCell(r, phi, theta).p;
                if (dir==0)
                {
                    auto r1 = grid.getCell(r + 1, phi, theta).p;
                    auto r_1 = grid.getCell(r - 1, phi, theta).p;
                    auto Dr = (r1 - r_1) / (2 * dr);
                    auto rR = (curr - r_1) / nonZeroDenom(r1 - curr);
                    grad.getCellRef(r, phi, theta).p = SlopeLim(rR) * Dr;
                }
                else if (dir==1)
                {
                    auto phi1 = grid.getCell(r, phi + 1, theta).p;
                    auto phi_1 = grid.getCell(r, phi - 1, theta).p;
                    auto Dphi = (phi1 - phi_1) / (2 * dphi * grid.getRFromIndex(r));
                    auto rPhi = (curr - phi_1) / nonZeroDenom(phi1 - curr);
                    grad.getCellRef(r, phi, theta).p = SlopeLim(rPhi) * Dphi;
                }
                else if (dir ==2 )
                {
                    auto theta1 = grid.getCell(r, phi, theta + 1).p;
                    auto theta_1 = grid.getCell(r, phi, theta - 1).p;
                    auto Dtheta = (theta1 - theta_1) / (2 * dtheta * grid.getRFromIndex(r));
                    auto rTheta = (curr - theta_1) / nonZeroDenom(theta1 - curr);
                    grad.getCellRef(r, phi, theta).p = SlopeLim(rTheta) * Dtheta;
                }
            }
        }
    }
}

void PredictorStep(SphericalGrid & grad, SphericalGrid &out, SphericalGrid &grid, double dt,int dir) {

    auto Dr = grad.getCell(0, 0, 0).zeros();
    auto Dtheta = grad.getCell(0, 0, 0).zeros();
    auto Dphi = grad.getCell(0, 0, 0).zeros();
    auto& Ddir =Dr;
    if (dir == 1) Ddir =Dphi;
    else if (dir == 2) Ddir =Dtheta;

    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                Cell &c = grid.getCellRef(r, phi, theta);
                Ddir = grad.getCell(r, phi, theta);
                out.getCellRef(r, phi, theta) = c - 0.5 * dt * F(Dr, Dtheta, Dphi, c);

            }
        }
    }
}

void ApplyFluxes(SphericalGrid &flux, SphericalGrid &grid, double dt,int dir) {
    double dr = grid.getRFromIndex(1) - grid.getRFromIndex(0);
    double dphi = (grid.getPhiFromIndex(1) - grid.getPhiFromIndex(0));
    double dtheta = grid.getThetaFromIndex(1) - grid.getThetaFromIndex(0);

    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                auto &c = grid.getCellRef(r, phi, theta);
                auto f= flux.getCell(r, phi, theta);
                if (dir==0)
                {
                    auto f1=flux.getCell(r + 1, phi, theta);
                    c.c = c.c - dt * ( f1.c*(f1.Sr/c.volume) - f.c*(f.Sr/c.volume));
                }
                else if (dir==1)
                {
                    auto f1=flux.getCell(r, phi + 1, theta);
                    c.c = c.c - dt *  (f1.c.rotate(dphi/2,0)*(f1.Sph/c.volume) - f.c.rotate(-dphi/2,0)*(f.Sph/c.volume)); //
                }
                else if (dir ==2 )
                {
                    auto f1=flux.getCell(r, phi, theta+1);
                    c.c = c.c - dt *  (f1.c.rotate(0,dtheta/2)*(f1.Sth/c.volume) - f.c.rotate(0,-dtheta/2)*(f.Sth/c.volume)); //
                }
            }
        }
    }
}

Simulation::Simulation(SphericalGrid &grid):grid(grid),step(0)
{
    densities=getDensity();
    vels=getVelocity();
    temperature=getTemperature();
    magneticField=getMagneticField();
    InitialConditions(grid,densities,vels,temperature,magneticField);
    ApplyBoundaryConditions(grid,0,densities, vels, temperature, magneticField);
}


void Simulation::RKIntegrator(double &dt, double &t) {
    //grid.UpdatePrim();

    dt = CalculateDT(grid);
    t += dt;
    for (int i=0;i<3;i++) {


        grid.UpdatePrim();

        ApplyBoundaryConditions(grid,t,densities, vels, temperature, magneticField);

        int dir = order[step];
        std::cout<<"dir "<<dir<<"\n";
        step++;
        if (step == 9) step = 0;
        //dir=1;
       // if (dir==1) return;

        CalculateGradients(grad, grid, dir);
        PredictorStep(grad, k, grid, dt, dir);


        if (dir == 0)
            CalculateFluxR(flux, k, grad);
        else if (dir == 1)
            CalculateFluxPh(flux, k, grad);
        else if (dir == 2)
            CalculateFluxTh(flux, k, grad);
        ApplyFluxes(flux, grid, dt, dir);

    }
    grid.UpdatePrim();
    ApplyBoundaryConditions(grid,t,densities, vels, temperature, magneticField);
}
