#include <algorithm>
#include <iostream>
#include <tuple>
#include "Simulation.h"
#include "Constants.h"


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
                    std::cout << "HERE" << std::endl;
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
}

void CalculateGradients(std::tuple<SphericalGrid &, SphericalGrid &, SphericalGrid &> grad, SphericalGrid &grid) {
    double dr = grid.getRFromIndex(1) - grid.getRFromIndex(0);
    double dphi = (grid.getPhiFromIndex(1) - grid.getPhiFromIndex(0));
    double dtheta = grid.getThetaFromIndex(1) - grid.getThetaFromIndex(0);

    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                auto r1 = grid.getCell(r + 1, phi, theta).p;
                auto r_1 = grid.getCell(r - 1, phi, theta).p;
                auto phi1 = grid.getCell(r, phi + 1, theta).p;
                auto phi_1 = grid.getCell(r, phi - 1, theta).p;
                auto theta1 = grid.getCell(r, phi, theta + 1).p;
                auto theta_1 = grid.getCell(r, phi, theta - 1).p;
                auto curr = grid.getCell(r, phi, theta).p;

                auto Dr = (r1 - r_1) / (2 * dr);
                auto Dphi = (phi1 - phi_1) / (2 * dphi * grid.getRFromIndex(r));
                auto Dtheta = (theta1 - theta_1) / (2 * dtheta * grid.getRFromIndex(r));

                /* auto Dr=      (r1     - curr)     / (dr);
                auto Dphi =   (phi1   - curr)   / (dphi * grid.getRFromIndex(r));
                auto Dtheta = (theta1 - curr) / (dtheta * grid.getRFromIndex(r));*/


                auto rR = (curr - r_1) / nonZeroDenom(r1 - curr);
                auto rPhi = (curr - phi_1) / nonZeroDenom(phi1 - curr);
                auto rTheta = (curr - theta_1) / nonZeroDenom(theta1 - curr);
                std::get<T_R>(grad).getCellRef(r, phi, theta).p = SlopeLim(rR) * Dr;
                std::get<T_PHI>(grad).getCellRef(r, phi, theta).p = SlopeLim(rPhi) * Dphi;
                std::get<T_THETA>(grad).getCellRef(r, phi, theta).p = SlopeLim(rTheta) * Dtheta;
            }
        }
    }
}

void PredictorStep(std::tuple<SphericalGrid &, SphericalGrid &, SphericalGrid &> grad, SphericalGrid &out,
                   SphericalGrid &grid, double dt) {
    SphericalGrid &DR = std::get<T_R>(grad);
    SphericalGrid &DPhi = std::get<T_PHI>(grad);
    SphericalGrid &DTheta = std::get<T_THETA>(grad);
    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                Cell &c = grid.getCellRef(r, phi, theta);
                if (r == 30 && phi == 50) {
                    int a = 0;
                }
                auto Dr = DR.getCell(r, phi, theta);
                auto Dtheta = DTheta.getCell(r, phi, theta);
                auto Dphi = DPhi.getCell(r, phi, theta);
                out.getCellRef(r, phi, theta) = c - 0.5 * dt * F(Dr, Dtheta, Dphi, c);

            }
        }
    }
}

void ApplyFluxes(std::tuple<SphericalGrid &, SphericalGrid &, SphericalGrid &> flux, SphericalGrid &grid, double dt) {
    double dr = grid.getRFromIndex(1) - grid.getRFromIndex(0);
    double dphi = (grid.getPhiFromIndex(1) - grid.getPhiFromIndex(0));
    double dtheta = grid.getThetaFromIndex(1) - grid.getThetaFromIndex(0);

    auto &fluxR = std::get<T_R>(flux);
    auto &fluxPhi = std::get<T_PHI>(flux);
    auto &fluxTheta = std::get<T_THETA>(flux);
    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                auto &c = grid.getCellRef(r, phi, theta).c;
                c = c - dt * (
                        (fluxR.getCell(r + 1, phi, theta).c - fluxR.getCell(r, phi, theta).c) // dr
                        + (fluxPhi.getCell(r, phi + 1, theta).c.rotate(dphi/2,0) - fluxPhi.getCell(r, phi, theta).c.rotate(-dphi/2,0)) //
                          //(dphi * grid.getRFromIndex(r))
                        + (fluxTheta.getCell(r, phi, theta + 1).c.rotate(0,dtheta/2) - fluxTheta.getCell(r, phi, theta).c.rotate(0,-dtheta/2)) //
                          //(dtheta * grid.getRFromIndex(r))
                );
            }
        }
    }

}

Simulation::Simulation(SphericalGrid &grid):grid(grid)
{

}


void Simulation::RKIntegrator(double &dt, double &t) {
    dt = CalculateDT(grid);
    t += dt;
    //std::cout << CalculateDT(grid) << "                   " << dt << std::endl;

    std::tuple<SphericalGrid &, SphericalGrid &, SphericalGrid &> flux(fluxR, fluxPhi, fluxTheta);
    std::tuple<SphericalGrid &, SphericalGrid &, SphericalGrid &> grad(DR, DPhi, DTheta);

    CalculateGradients(grad, grid);
    PredictorStep(grad, k, grid, dt);
    CalculateFlux(flux, k, grad);

    ApplyFluxes(flux, grid, dt);
    grid.UpdatePrim();


}
