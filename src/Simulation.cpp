#include <algorithm>
#include <iostream>
#include <tuple>
#include "Simulation.h"
#include "Constants.h"

Cell SlopeLim(Cell r)
{
    /* return {std::max(0.0,std::max(std::min(2*r.rho,1.0),std::min(r.rho,2.0))),
             std::max(0.0,std::max(std::min(2*r.Vr,1.0),std::min(r.Vr,2.0))),
             std::max(0.0,std::max(std::min(2*r.Vphi,1.0),std::min(r.Vphi,2.0))),
                 std::max(0.0,std::max(std::min(2*r.B,1.0),std::min(r.B,2.0)))};*/
    return {std::max(0.0, std::min(1.0, r.rho)),
            std::max(0.0, std::min(1.0, r.Vr)),
            std::max(0.0, std::min(1.0, r.Vphi)),
            std::max(0.0, std::min(1.0, r.Vtheta)),
            std::max(0.0, std::min(1.0, r.Br)),
            std::max(0.0, std::min(1.0, r.Bphi)),
            std::max(0.0, std::min(1.0, r.Btheta)),
            std::max(0.0, std::min(1.0, r.E))};

    /* return Cell{std::max(0.0, 1.5 * (r.rho * r.rho + r.rho) / (r.rho * r.rho + r.rho + 1)),
                 std::max(0.0, 1.5 * (r.Vr * r.Vr + r.Vr) / (r.Vr * r.Vr + r.Vr + 1)),
                          std::max(0.0, 1.5 * (r.Vphi * r.Vphi + r.Vphi) / (r.Vphi * r.Vphi + r.Vphi + 1)),
                          std::max(0.0, 1.5 * (r.B * r.B + r.B) / (r.B * r.B + r.B + 1))};*/



}

double nonZeroDouble(double val)
{
    if (std::abs(val) < __DBL_EPSILON__)
    {
        if (val < 0)
            return -__DBL_EPSILON__;
        else
            return __DBL_EPSILON__;
    }
    return val;
}

Cell nonZeroDenom(Cell denom)
{
    return {nonZeroDouble(denom.rho),
            nonZeroDouble(denom.Vr),
            nonZeroDouble(denom.Vphi),
            nonZeroDouble(denom.Vtheta),
            nonZeroDouble(denom.Br),
            nonZeroDouble(denom.Bphi),
            nonZeroDouble(denom.Btheta),
            nonZeroDouble(denom.E)};
}

Cell F(Cell Dr,Cell Dtheta, Cell Dphi, Cell U, double r, double phi, double theta)
{
    theta=M_PI_2;
    /*double Vr = U.Vr/U.rho;
    double Vrdr = Dr.Vr - Dr.rho * Vr;
    double Vrdphi = Dphi.Vr - Dphi.rho * Vr;
    double Vrdtheta = Dtheta.Vr - Dtheta.rho * Vr;

    double Vtheta = U.Vtheta / U.rho;
    double Vthetadr = Dr.Vtheta - Dr.rho * Vtheta;
    double Vthetadphi = Dphi.Vtheta - Dphi.rho * Vtheta;
    double Vthetadtheta = Dtheta.Vtheta - Dtheta.rho * Vtheta;

    double Vphi = U.Vphi / U.rho;
    double Vphidr = Dr.Vphi - Dr.rho * Vphi;
    double Vphidphi = Dphi.Vphi - Dphi.rho * Vphi;
    double Vphidtheta = Dtheta.Vphi - Dtheta.rho * Vphi;*/

    // rho * vr = rho/dr * vr + rho * vr/dr
    //vr/dr = (rho * vr - rho/dr * vr)/rho
    double p = gamma * (U.E - 0.5 * U.rho * (U.Vr * U.Vr + U.Vtheta * U.Vtheta + U.Vphi * U.Vphi)
                        - 0.5/mu * (U.Br * U.Br + U.Bphi * U.Bphi + U.Btheta * U.Btheta));
    double dpdr = gamma * (Dr.E - 0.5 * Dr.rho * (U.Vr * U.Vr + U.Vphi * U.Vphi + U.Vtheta * U.Vtheta)
                           - U.rho * (Dr.Vr * U.Vr + Dr.Vphi * U.Vphi + Dr.Vtheta * U.Vtheta)
                           - (Dr.Br * U.Br + Dr.Bphi * U.Bphi + Dr.Btheta * U.Btheta)/mu);
    double dpdphi = gamma * (Dphi.E - 1.0/2*Dphi.rho* (U.Vr*U.Vr + U.Vphi*U.Vphi + U.Vtheta*U.Vtheta)
                           - U.rho * (Dphi.Vr * U.Vr + Dphi.Vphi * U.Vphi + Dphi.Vtheta * U.Vtheta)
                           - (Dphi.Br * U.Br + Dphi.Bphi * U.Bphi + Dphi.Btheta * U.Btheta)/mu);
    double dpdtheta =gamma * (Dtheta.E - 0.5 * Dtheta.rho * (U.Vr * U.Vr + U.Vphi * U.Vphi + U.Vtheta * U.Vtheta)
                              - U.rho * (Dtheta.Vr * U.Vr + Dtheta.Vphi * U.Vphi + Dtheta.Vtheta * U.Vtheta)
                              - (Dtheta.Br * U.Br + Dtheta.Bphi * U.Bphi + Dtheta.Btheta * U.Btheta)/mu);

    double P = p + 0.5 / mu * (U.Br * U.Br + U.Bphi * U.Bphi + U.Btheta * U.Btheta);
    double dPdr = dpdr + (Dr.Br * U.Br + Dr.Bphi * U.Bphi + Dr.Btheta * U.Btheta) / mu;
    double dPdphi = dpdphi + (Dphi.Br * U.Br + Dphi.Bphi * U.Bphi + Dphi.Btheta * U.Btheta) / mu;
    double dPdtheta = dpdtheta + (Dtheta.Br * U.Br + Dtheta.Bphi * U.Bphi + Dtheta.Btheta * U.Btheta) / mu;

    /*double drhodt=   Dr.rho * U.Vr + U.rho * Dr.Vr +
                     Dphi.rho * U.Vphi + U.rho * Dphi.Vphi;*/

   //double drhodt = U.Vr*Dr.rho + U.rho*Dr.Vr + 2.0/r*U.rho*U.Vr;
    double drhodt = 1. /(r * r) * (2 * r * U.rho * U.Vr + r * r * (Dr.rho * U.Vr + U.rho * Dr.Vr))
                 + 1.0 /(r * std::sin(theta)) * (std::cos(theta) * U.rho * U.Vtheta + std::sin(theta) * (Dtheta.rho * U.Vtheta + U.rho * Dtheta.Vtheta))
                 + 1.0 / r * (U.rho * Dphi.Vphi + U.Vphi * Dphi.rho);
    /*double drhodt = 1./(r*r) * (2 * r * U.Vr + r*r * Dr.Vr)
                    + 1.0/(r * std::sin(theta))*(std::cos(theta)*U.Vtheta + std::sin(theta)*Dtheta.Vtheta)
                    + 1.0/r * Dphi.Vphi;*/

    double Vrdt = -((U.rho * U.Vtheta*U.Vtheta + U.rho*U.Vphi*U.Vphi) / r
                        - dPdr
                        - (U.rho*G*Ms) / (r * r)
                        + U.Br / mu * Dr.Br
                        + 1.0/(r * mu) * (Dtheta.Br * U.Btheta + U.Br * Dtheta.Btheta)
                        + U.Bphi / mu * Dphi.Br
                        - (U.Btheta * U.Btheta + U.Bphi * U.Bphi) / (mu * r))
                     + 1.0/(r*r) * ( 2*r*U.rho * U.Vr * U.Vr + r*r*(Dr.rho * U.Vr*U.Vr + 2 * U.rho * U.Vr * Dr.Vr)
                     + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.rho*U.Vr*U.Vtheta + std::sin(theta)*(Dtheta.rho*U.Vr*U.Vtheta + U.rho*(Dtheta.Vr*U.Vtheta + U.Vr*Dtheta.Vtheta)))
                     + 1.0/(r*std::sin(theta)) * (Dphi.rho*U.Vr*U.Vphi + (Dphi.Vr*U.Vphi + U.Vr*Dphi.Vphi)));

    double Vphidt = -((-U.rho * U.Vr*U.Vphi - U.rho*U.Vtheta*U.Vphi*ctg(theta)) / r
                       - dPdphi/r
                       + U.Br / mu * Dr.Bphi
                       + 1.0/(r * r * mu) * (2*r*U.Bphi*U.Br + r*r*(Dr.Bphi * U.Br + U.Bphi * Dr.Br))
                       + 1.0 / (mu*r) *(Dtheta.Bphi * U.Btheta + U.Bphi * Dtheta.Btheta)
                       + U.Bphi/(r*mu) * Dphi.Bphi
                       + (U.Br * U.Bphi + U.Btheta * U.Bphi*ctg(theta)) / (mu * r))
                       + 1.0/(r*r) * ( 2*r*U.rho * U.Vphi * U.Vr + r*r*(Dr.rho * U.Vphi*U.Vr +U.rho*( Dr.Vphi* U.Vr + U.Vphi * Dr.Vr))
                       + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.rho*U.Vphi*U.Vtheta + std::sin(theta)*(Dtheta.rho*U.Vphi*U.Vtheta + U.rho*(Dtheta.Vphi * U.Vtheta + U.Vphi*Dtheta.Vtheta)))
                       + 1.0/(r*std::sin(theta)) * (Dphi.rho*U.Vphi*U.Vphi + 2*U.rho*(Dphi.Vphi*U.Vphi)));

    double Vthetadt = -((-U.rho * U.Vr*U.Vtheta + U.rho*U.Vphi*U.Vphi*ctg(theta)) / r
                           - dPdtheta/r
                           + 1.0/(r * r * mu) * (2*r*U.Btheta*U.Br + r*r*(Dr.Btheta * U.Br + U.Btheta * Dr.Br))
                           + 1.0/(r * std::sin(theta) * mu) * (std::cos(theta)*U.Btheta*U.Btheta + std::sin(theta) * 2 * U.Btheta * Dtheta.Btheta)
                           + 1.0 / (r * std::sin(theta) * mu) * (Dphi.Btheta * U.Bphi + U.Btheta * Dphi.Bphi)
                           + U.Br*U.Btheta/(r*mu)
                           - (U.Bphi * U.Bphi *ctg(theta)) / (mu * r))
                         + 1.0/(r*r) * ( 2*r*U.rho * U.Vtheta * U.Vr + r*r*(Dr.rho * U.Vtheta * U.Vr + U.rho*(Dr.Vtheta* U.Vr + U.Vtheta * Dr.Vr))
                         + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.rho*U.Vtheta*U.Vtheta + std::sin(theta)*(Dtheta.rho*U.Vtheta*U.Vtheta + 2* U.rho*(Dtheta.Vtheta*U.Vtheta)))
                         + 1.0/(r*std::sin(theta)) * (Dphi.rho*U.Vtheta*U.Vphi + U.rho*(Dtheta.Vphi*U.Vtheta + U.Vphi*Dtheta.Vtheta)));


   /* double Edt = 1.0/(r*r) * (2*r*U.E*Vr + r*r*(Dr.E*Vr + U.E*Vrdr/U.rho))
                 + 1.0/(r*std::sin(theta)) * (std::cos(theta)*U.E*Vtheta + std::sin(theta)*(Dtheta.E*Vtheta + U.E*Vthetadtheta/U.rho))
                 + 1.0/(r*std::sin(theta)) * (Dphi.E*Vphi + U.E*Vphidphi/U.rho)
                 -(- 1.0/(r*r) * (2*r*P*Vr + r*r*(dPdr*Vr + P*Vrdr/U.rho))
                   - 1.0/(r*std::sin(theta)) * (std::cos(theta)*P*Vtheta + std::sin(theta)*(dPdtheta*Vtheta + P*Vthetadtheta/U.rho))
                   - 1.0/(r*std::sin(theta)) *(dPdphi*Vphi + P * Vphidphi/U.rho)
                   - U.rho*Vr*G*Ms/(r*r)
                   + 1.0/(r*r*mu) * (2*r*U.Br*(U.Br*Vr+U.Btheta*Vtheta+U.Bphi*Vphi) + r*r*(Dr.Br*(U.Br*Vr+U.Btheta*Vtheta+U.Bphi*Vphi)
                        + U.Br*((Dr.Br*Vr + U.Br*Vrdr/U.rho)+(Dr.Btheta*Vtheta + U.Btheta*Vthetadr/U.rho)+(Dr.Bphi*Vphi + U.Bphi*Vphidr/U.rho))))
                   + 1.0/(r*std::sin(theta)*mu) *(std::cos(theta)*U.Btheta*(U.Br*Vr + U.Btheta*Vtheta+U.Bphi*Vphi)
                        + std::sin(theta)*(Dtheta.Btheta*(U.Br*Vr+U.Btheta*Vtheta+U.Bphi*Vphi)
                        + U.Btheta*((Dtheta.Br*Vr + U.Br*Vrdtheta/U.rho)+(Dtheta.Btheta*Vtheta + U.Btheta*Vthetadtheta/U.rho)+(Dtheta.Bphi*Vphi + U.Bphi*Vphidtheta/U.rho))))
                   + 1.0/(r*std::sin(theta)*mu) * (Dphi.Bphi*(U.Br*Vr+U.Btheta*Vtheta+U.Bphi*Vphi)
                        + U.Bphi*((Dphi.Br*Vr + U.Br*Vrdphi/U.rho)+(Dphi.Btheta*Vtheta + U.Btheta*Vthetadphi/U.rho)+(Dphi.Bphi*Vphi + U.Bphi*Vphidphi/U.rho)))
                         );*/

    if(std::abs(drhodt)>1)
        int a=0;



    return Cell{drhodt,
                Vrdt,
                Vphidt,
                Vthetadt,
                0,
                0,
                0,
                0};
}

double maxSpeed(Cell U){
    double p = (gamma-1) * (U.E - 1.0/2*U.rho* (U.Vr * U.Vr + U.Vphi * U.Vphi + U.Vtheta * U.Vtheta)
                        - 1.0/2*(U.Br * U.Br + U.Bphi * U.Bphi + U.Btheta * U.Btheta));
    double P = p + 1.0/2*(U.Br * U.Br + U.Bphi * U.Bphi + U.Btheta * U.Btheta);

    double B_2 = U.Br * U.Br + U.Bphi * U.Bphi + U.Btheta * U.Btheta;
    double B = std::sqrt(B_2);

    double cmax =
            std::sqrt(U.Vr * U.Vr + U.Vphi * U.Vphi + U.Vtheta * U.Vtheta)
            + 0.5*((gamma * p + B_2)/U.rho
            + std::sqrt(((gamma + B)/U.rho) * ((gamma + B)/U.rho) - 4 * (gamma * U.Br * U.Br) / (U.rho * U.rho)));

    return cmax*DT/CELL_SIZE;
}


void CalculateFlux(std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> out, SphericalGrid& in)
{
    for(int theta = 0; theta < in.getSizeTheta(); theta++) {
        for (int r = 0; r < in.getSizeR(); r++) {
            for (int phi = 0; phi < in.getSizePhi(); phi++) // 2->n ??
            {
                auto rR_current = (in.getCell(r,phi,theta) - in.getCell(r-1,phi,theta)) / nonZeroDenom(in.getCell(r+1,phi,theta) - in.getCell(r,phi,theta));
                auto rR_prev = (in.getCell(r-1,phi,theta) - in.getCell(r-2,phi,theta)) / nonZeroDenom(in.getCell(r,phi,theta) - in.getCell(r-1,phi,theta));

                auto rphi_current = (in.getCell(r,phi,theta) - in.getCell(r,phi-1,theta)) / nonZeroDenom(in.getCell(r,phi+1,theta) - in.getCell(r,phi,theta));
                auto rphi_prev = (in.getCell(r,phi-1,theta) - in.getCell(r,phi-2,theta)) / nonZeroDenom(in.getCell(r,phi,theta) - in.getCell(r,phi-1,theta));

                auto rtheta_current = (in.getCell(r,phi,theta) - in.getCell(r,phi,theta-1)) / nonZeroDenom(in.getCell(r,phi,theta+1) - in.getCell(r,phi,theta));
                auto rtheta_prev = (in.getCell(r,phi,theta-1) - in.getCell(r,phi,theta-2)) / nonZeroDenom(in.getCell(r,phi,theta) - in.getCell(r,phi,theta-1));

                auto uL_R = in.getCell(r-1,phi,theta) + 0.5 * SlopeLim(rR_prev) * (in.getCell(r,phi,theta) - in.getCell(r-1,phi,theta));
                auto uR_R = in.getCell(r,phi,theta) - 0.5 * SlopeLim(rR_current) * (in.getCell(r+1,phi,theta) - in.getCell(r,phi,theta));

                auto uL_phi = in.getCell(r,phi-1,theta) + 0.5 * SlopeLim(rphi_prev) * (in.getCell(r,phi,theta) - in.getCell(r,phi-1,theta));
                auto uR_phi = in.getCell(r,phi,theta) - 0.5 * SlopeLim(rphi_current) * (in.getCell(r,phi+1,theta) - in.getCell(r,phi,theta));

                auto uL_theta = in.getCell(r,phi,theta-1) + 0.5 * SlopeLim(rtheta_prev) * (in.getCell(r,phi,theta) - in.getCell(r,phi,theta-1));
                auto uR_theta = in.getCell(r,phi,theta) - 0.5 * SlopeLim(rtheta_current) * (in.getCell(r,phi,theta+1) - in.getCell(r,phi,theta));

#if defined(USE_CONST_A)
                double aR = A_SPEED;
                double aphi = A_SPEED;
                double atheta = A_SPEED;
#else
                double aR = maxSpeed(0.5 * (uL_R + uR_R));
                double aphi = maxSpeed(0.5 * (uL_phi + uR_phi));
                double atheta = maxSpeed(0.5 * (uL_theta + uR_theta));
#endif


                //std::cout << a << std::endl;
                // F{i-0.5} = 0.5 * (F(uR{i-0.5}) + F(uL{i-0.5}) - a * (uR{i-0.5} - uL{i-0.5}))
                std::get<0>(out).getCellRef(r,phi,theta) = (0.5 *
                        (F(uR_R, Cell::zeros(),Cell::zeros(),
                           in.getCell(r,phi,theta),
                           in.getRFromIndex(r),
                           in.getPhiFromIndex(phi),
                           in.getThetaFromIndex(theta))
                        + F(uL_R, Cell::zeros(),Cell::zeros(),
                            in.getCell(r,phi,theta),
                            in.getRFromIndex(r),
                            in.getPhiFromIndex(phi),
                            in.getThetaFromIndex(theta))
                        - aR * (uR_R - uL_R)));
                std::get<1>(out).getCellRef(r,phi,theta) = (0.5 *
                        (F(Cell::zeros(), Cell::zeros(),uR_phi,
                            in.getCell(r,phi,theta),
                            in.getRFromIndex(r),
                            in.getPhiFromIndex(phi),
                            in.getThetaFromIndex(theta))
                        + F(Cell::zeros(), Cell::zeros(),uL_phi,
                             in.getCell(r,phi,theta),
                             in.getRFromIndex(r),
                             in.getPhiFromIndex(phi),
                             in.getThetaFromIndex(theta))
                        - aphi * (uR_phi - uL_phi)));
                std::get<2>(out).getCellRef(r,phi,theta) = (0.5 *
                         (F(Cell::zeros(), uR_theta,Cell::zeros(),
                            in.getCell(r,phi,theta),
                            in.getRFromIndex(r),
                            in.getPhiFromIndex(phi),
                            in.getThetaFromIndex(theta))
                         + F(Cell::zeros(), uL_theta,Cell::zeros(),
                             in.getCell(r,phi,theta),
                             in.getRFromIndex(r),
                             in.getPhiFromIndex(phi),
                             in.getThetaFromIndex(theta))
                         - atheta * (uR_theta - uL_theta)));

            }

        }
    }
}

Cell S(int x,int y,Cell val)
{
    if (x>3 && x<7 && y>32 && y< 58)
        return {0,0,0,0};
    else return {0,0,0,0};
}

void InitialConditions(SphericalGrid& grid) {
    for (int x = 0; x < grid.getSizeR(); x++) {
        for (int y = 0; y < grid.getSizePhi(); y++) {
            double rho=0.06/std::pow(grid.getRFromIndex(x),2);
            double vx=100;
            double vy=0;
            double vz=0;
            double Bx=0.000;
            double By=0.000;
            double Bz=0.000;
            double T =273;
            // for (int y = 45; y < 55; y++){
            //        for (int x = 20; x < 40; x++){
            if(x>20 && x< 40 && y>45 && y<55)
            {
                rho=0.06*std::cos((x-30.0)/10);
               // vx=100;
            }


            //nk=rho/m m=1.6733e-27
            //k=1.38044e-23
            double E = 2 * rho * m_div_k * T/(gamma-1)
                    + rho * (vx*vx + vy*vy + vz*vz)
                    + (Bx*Bx + By*By + Bz*Bz) /(2*mu);

            grid.getCellRef(x,y,0) = {rho, vx, vy, vz,Bx,By,Bz,E};

        }
    }
}


void ApplyBoundaryConditions(SphericalGrid& grid)
{
    /*for (int x = 20; x < 40; x++){
        for (int y = 45; y < 55; y++){
            double rho = 0.02;
            double Vr = 0;
            double Vphi = 0;
            double Vtheta = 0;
            double Br = 0.000;
            double Bphi = 0.000;
            double Btheta = 0.000;
            double T = 273;

            double E = 2 * rho * m_div_k * T/(gamma-1)
                       + rho * (Vr*Vr + Vphi*Vphi + Vtheta*Vtheta)
                       + (Br*Br + Bphi*Bphi + Btheta*Btheta) /(2*mu);

            grid.mesh[x][y].rho = rho;
            grid.mesh[x][y].E = E;
        }
    }
*/
    double rho =0.01;
   /* for (int y=0;y<grid.getSizeR();y++) {
        grid.mesh[0][y].rho=rho;
        grid.mesh[0][y].E=E;
        grid.mesh[grid.sizeX][y].rho=rho;
        grid.mesh[grid.sizeX][y].E=E;
    }*/
    for (int x=0;x<grid.getSizePhi();x++) {
        Cell& c = grid.getCellRef(1,x,0);
        c.rho=rho;
        c.E=10000;
        //c.Vr=0;
        //c.Vphi=0;
        //c.Vtheta=0;
    }
}



void RKIntegrator(SphericalGrid& grid, double dt)
{
    double dx = CELL_SIZE;
    SphericalGrid fluxR = SphericalGrid::copyGrid(grid);
    SphericalGrid fluxPhi = SphericalGrid::copyGrid(grid);
    SphericalGrid fluxTheta = SphericalGrid::copyGrid(grid);
    std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> flux(fluxR,fluxPhi,fluxTheta);
    CalculateFlux(flux,grid);
    SphericalGrid k = SphericalGrid::copyGrid(grid);
    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                k.getCellRef(r,phi,theta) = grid.getCell(r,phi,theta) - 0.5 * dt / dx *
                                                 ( std::get<0>(flux).getCell(r+1,phi,theta) - std::get<0>(flux).getCell(r,phi,theta)
                                                  +std::get<1>(flux).getCell(r,phi+1,theta) - std::get<1>(flux).getCell(r,phi,theta))
                                                  + 0.5 * dt * S(phi, r, grid.getCell(r,phi,theta));
            }
        }
    }
    ApplyBoundaryConditions(k);
    CalculateFlux(flux, k);
    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                grid.getCellRef(r,phi,theta) = k.getCell(r,phi,theta) - dt / dx *
                                                                        (std::get<0>(flux).getCell(r+1,phi,theta) - std::get<0>(flux).getCell(r,phi,theta)
                                                                        +std::get<1>(flux).getCell(r,phi+1,theta) - std::get<1>(flux).getCell(r,phi,theta))
                                                                    +  dt * S(phi, r, grid.getCell(r,phi,theta));
            }
        }
    }
    ApplyBoundaryConditions(grid);

}
