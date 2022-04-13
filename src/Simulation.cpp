#include <algorithm>
#include <iostream>
#include "Simulation.h"
#include "Constants.h"

Cell SlopeLim(Cell r)
{
    /* return {std::max(0.0,std::max(std::min(2*r.rho,1.0),std::min(r.rho,2.0))),
             std::max(0.0,std::max(std::min(2*r.vx,1.0),std::min(r.vx,2.0))),
             std::max(0.0,std::max(std::min(2*r.vy,1.0),std::min(r.vy,2.0))),
                 std::max(0.0,std::max(std::min(2*r.B,1.0),std::min(r.B,2.0)))};*/
    return {std::max(0.0, std::min(1.0, r.rho)),
            std::max(0.0, std::min(1.0, r.vx)),
            std::max(0.0, std::min(1.0, r.vy)),
            std::max(0.0, std::min(1.0, r.vz)),
            std::max(0.0, std::min(1.0, r.Bx)),
            std::max(0.0, std::min(1.0, r.By)),
            std::max(0.0, std::min(1.0, r.Bz)),
            std::max(0.0, std::min(1.0, r.E))};

    /* return Cell{std::max(0.0, 1.5 * (r.rho * r.rho + r.rho) / (r.rho * r.rho + r.rho + 1)),
                 std::max(0.0, 1.5 * (r.vx * r.vx + r.vx) / (r.vx * r.vx + r.vx + 1)),
                          std::max(0.0, 1.5 * (r.vy * r.vy + r.vy) / (r.vy * r.vy + r.vy + 1)),
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
            nonZeroDouble(denom.vx),
            nonZeroDouble(denom.vy),
            nonZeroDouble(denom.vz),
            nonZeroDouble(denom.Bx),
            nonZeroDouble(denom.By),
            nonZeroDouble(denom.Bz),
            nonZeroDouble(denom.E)};
}

Cell F(Cell DX,Cell DY,Cell U)
{
    Cell DZ=Cell::zeros();


    double p = (gamma-1) * (U.E - 1.0/2 * U.rho* (U.vx * U.vx + U.vy * U.vy + U.vz * U.vz)
                        - 1.0/(2*mu)*(U.Bx*U.Bx + U.By*U.By + U.Bz * U.Bz));

    double pdx =  (gamma-1) * (DX.E - DX.rho * (U.vx * U.vx + U.vy * U.vy + U.vz * U.vz) / 2
                           - U.rho * (DX.vx * U.vx + DX.vy * U.vy + DX.vz * U.vz)
                           - (DX.Bx * U.Bx + DX.By * U.By + DX.Bz * U.Bz)/mu) ;

    double pdy =  (gamma-1) * (DY.E - DY.rho* (U.vx*U.vx + U.vy*U.vy + U.vz*U.vz) / 2
                           - U.rho * (DY.vx * U.vx + DY.vy * U.vy + DY.vz * U.vz)
                           - (DY.Bx * U.Bx + DY.By * U.By + DY.Bz * U.Bz)/mu) ;

    double pdz =  (gamma-1) * (DZ.E - DZ.rho* (U.vx*U.vx + U.vy*U.vy + U.vz*U.vz) / 2
                           - U.rho * (DZ.vx * U.vx + DZ.vy * U.vy + DZ.vz * U.vz)
                           - (DZ.Bx * U.Bx + DZ.By * U.By + DZ.Bz * U.Bz)/mu) ;

    //if(std::abs(pdx>0.001))
        //std::cout<<pdx<<"\n";


    double P = p + (U.Bx * U.Bx + U.By * U.By + U.Bz * U.Bz) / (2 * mu);

    double Pdx = (pdx + (DX.Bx * U.Bx + DX.By * U.By + DX.Bz * U.Bz) / mu);
    double Pdy = (pdy + (DY.Bx * U.Bx + DY.By * U.By + DY.Bz * U.Bz) / mu);
    double Pdz = (pdz + (DZ.Bx * U.Bx + DZ.By * U.By + DZ.Bz * U.Bz) / mu);

    double drhodt = DX.rho * U.vx + DY.rho * U.vy + DZ.rho * U.vz
                    + U.rho * (DX.vx + DY.vy + DZ.vz);


    double dvxdt = (DX.rho*U.vx*U.vx + 2 * U.rho * DX.vx * U.vx
                   +DY.rho*U.vx*U.vy + U.rho * DY.vx * U.vy + U.rho * U.vx * DY.vy
                   +DZ.rho*U.vx*U.vz + U.rho * DZ.vx * U.vz + U.rho * U.vx * DZ.vz
                   +Pdx// * U.rho//*0.01
                   -2 * U.Bx*DX.Bx / mu
                   -(DY.Bx * U.By + U.Bx * DY.By) / mu
                   -(DZ.Bx * U.Bz + U.Bx * DZ.Bz) / mu
                   -drhodt * U.vx);//* U.rho;


    double dvydt = (DX.rho * U.vx * U.vy + U.rho * DX.vx * U.vy + U.rho * U.vx * DX.vy
                   +DY.rho * U.vy * U.vy + 2 * U.rho * DY.vy * U.vy
                   +DZ.rho * U.vz * U.vy + U.rho * DZ.vz * U.vy + U.rho * U.vz * DZ.vy
                   +Pdy// * U.rho//*0.01
                   -(DX.Bx * U.By + U.Bx * DX.By) / mu
                   -2 * U.By * DY.By / mu
                   -(DZ.Bz * U.By + U.Bz * DZ.By) / mu
                   -drhodt * U.vy);//*U.rho;

    double dvzdt = (DX.rho * U.vx * U.vz + U.rho * DX.vx * U.vz + U.rho * U.vx * DX.vz
                   +DY.rho * U.vy * U.vz + U.rho * DY.vy * U.vz + U.rho * U.vy * DY.vz
                   +DZ.rho * U.vz * U.vz + 2 * U.rho * DZ.vz * U.vz
                   +Pdz// * U.rho//*0.01
                   -(DX.Bx * U.Bz + U.Bx * DX.Bz) / mu
                   -(DY.By * U.Bz + U.By * DY.Bz) / mu
                   -2 * U.Bz * DZ.Bz / mu
                   -drhodt * U.vz);// * U.rho;

    double dBxdt =  DY.vx * U.By + DY.By * U.vx
                    - DY.vy * U.Bx - DY.Bx * U.vy
                    + DZ.vx * U.Bz + DZ.Bz * U.vx
                    - DZ.vz * U.Bx - DZ.Bx * U.vz;

    double dBydt =  - DX.vx * U.By - DX.By * U.vx
                    + DX.vy * U.Bx + DX.Bx * U.vy
                    + DZ.vy * U.Bz + DZ.Bz * U.vy
                    - DZ.vz * U.By - DZ.By * U.vz;

    double dBzdt = - DX.vx * U.Bz - DX.Bz * U.vx
                   + DX.vz * U.Bx + DX.Bx * U.vz
                   - DY.vy * U.Bz - DY.Bz * U.vy
                   + DY.vz * U.By + DY.By * U.vz;



    double dEdt = (DX.E * U.vx + DY.E * U.vy + DZ.E * U.vz) + U.E * (DX.vx + DY.vy + DZ.vz)    //nabla(EV)
                  +(Pdx * U.vx + Pdy * U.vy + Pdz * U.vz) + P * (DX.vx + DY.vy + DZ.vz)  ;       //nabla(PV)

                  + (DX.vx * U.Bx * U.Bx + 2 * DX.Bx * U.Bx * U.vx
                  + DX.vy * U.Bx * U.By + DX.Bx * U.By * U.vy + DX.By * U.Bx * U.vy             //DX
                  + DX.vz * U.Bx * U.Bz + DX.Bx * U.Bz * U.vz + DX.Bz * U.Bx * U.vz

                  + DY.By * U.Bx * U.vx + DY.Bx * U.By * U.vx  +DY.vx * U.By * U.Bx
                  + 2 * DY.By * U.By * U.vy + DY.vy * U.By * U.By
                  + DY.By * U.Bz * U.vz + DY.Bz * U.By * U.vz + DY.vz * U.By * U.Bz

                  + DZ.vx * U.Bz * U.Bx + DZ.Bz * U.Bx * U.vx + DZ.Bx * U.Bz * U.vx
                  + DZ.vy * U.Bz * U.By + DZ.Bz * U.By * U.vy + DZ.By * U.Bz * U.vy
                  + DZ.vz * U.Bz * U.Bz + 2 * DZ.Bz * U.Bz * U.vz)/mu;


    return Cell{drhodt,
                dvxdt ,
                dvydt ,
                dvzdt,
                dBxdt,
                dBydt,
                dBzdt,
                dEdt};
}

double maxSpeed(Cell U){
    double p = (gamma-1) * (U.E - 1.0/2*U.rho* (U.vx*U.vx + U.vy*U.vy + U.vz*U.vz)
                        - 1.0/2*(U.Bx*U.Bx + U.By*U.By + U.Bz*U.Bz));
    double P = p + 1.0/2*(U.Bx*U.Bx + U.By*U.By + U.Bz*U.Bz);

    double B_2 = U.Bx*U.Bx + U.By*U.By + U.Bz*U.Bz;
    double B = std::sqrt(B_2);

    double cmax =
            std::sqrt(U.vx*U.vx + U.vy*U.vy + U.vz*U.vz)
            + 0.5*((gamma * p + B_2)/U.rho
            + std::sqrt(((gamma + B)/U.rho) * ((gamma + B)/U.rho) - 4*(gamma*U.Bx*U.Bx)/(U.rho*U.rho)));

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
                std::get<0>(out).getCellRef(r,phi,theta) = (0.5 * (F(uR_R, Cell::zeros(), in.getCell(r,phi,theta)) + F(uL_R, Cell::zeros(), in.getCell(r,phi,theta)) -
                        aR * (uR_R - uL_R)));
                std::get<1>(out).getCellRef(r,phi,theta) = (0.5 * (F(uR_phi, Cell::zeros(), in.getCell(r,phi,theta)) + F(uL_phi, Cell::zeros(), in.getCell(r,phi,theta)) -
                                                                aphi * (uR_phi - uL_phi)));
                std::get<2>(out).getCellRef(r,phi,theta) = (0.5 * (F(uR_theta, Cell::zeros(), in.getCell(r,phi,theta)) + F(uL_theta, Cell::zeros(), in.getCell(r,phi,theta)) -
                                                                atheta * (uR_theta - uL_theta)));

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
            double rho=0.01;
            double vx=0;
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
                rho=0.08*std::cos((x-30.0)/10);
                vx=0;
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
            double vx = 0;
            double vy = 0;
            double vz = 0;
            double Bx = 0.000;
            double By = 0.000;
            double Bz = 0.000;
            double T = 273;

            double E = 2 * rho * m_div_k * T/(gamma-1)
                       + rho * (vx*vx + vy*vy + vz*vz)
                       + (Bx*Bx + By*By + Bz*Bz) /(2*mu);

            grid.mesh[x][y].rho = rho;
            grid.mesh[x][y].E = E;
        }
    }
*/
    /*for (int y=0;y<grid.sizeY;y++) {
        double rho = 0.01;
        double vx = 0;
        double vy = 0;
        double vz = 0;
        double Bx = 0.000;
        double By = 0.000;
        double Bz = 0.000;
        double T = 273;

        double E = 2 * rho * m_div_k * T/(gamma-1)
                   + rho * (vx*vx + vy*vy + vz*vz)
                   + (Bx*Bx + By*By + Bz*Bz) /(2*mu);

        grid.mesh[0][y].rho=rho;
        grid.mesh[0][y].E=E;
        grid.mesh[grid.sizeX][y].rho=rho;
        grid.mesh[grid.sizeX][y].E=E;
    }
    for (int x=0;x<grid.sizeX;x++) {
        double rho = 0.01;
        double vx = 0;
        double vy = 0;
        double vz = 0;
        double Bx = 0.000;
        double By = 0.000;
        double Bz = 0.000;
        double T = 273;

        double E = 2 * rho * m_div_k * T/(gamma-1)
                   + rho * (vx*vx + vy*vy + vz*vz)
                   + (Bx*Bx + By*By + Bz*Bz) /(2*mu);


        grid.mesh[x][0].rho=rho;
        grid.mesh[x][0].E=E;
        grid.mesh[x][grid.sizeY].rho=rho;
        grid.mesh[x][grid.sizeY].E=E;

    }*/
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
                                                 (std::get<0>(flux).getCell(r,phi+1,theta) - std::get<0>(flux).getCell(r,phi,theta) +
                                                  std::get<1>(flux).getCell(r+1,phi,theta) - std::get<1>(flux).getCell(r,phi,theta))
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
                                                                        (std::get<0>(flux).getCell(r,phi+1,theta) - std::get<0>(flux).getCell(r,phi,theta) +
                                                                         std::get<1>(flux).getCell(r+1,phi,theta) - std::get<1>(flux).getCell(r,phi,theta))
                                                                    +  dt * S(phi, r, grid.getCell(r,phi,theta));
            }
        }
    }
    ApplyBoundaryConditions(grid);

}
