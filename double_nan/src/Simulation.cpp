#include <algorithm>
#include <iostream>
#include "Simulation.h"

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

    double gamma = 5.0/3 - 1;
    double p = gamma * (U.E - 1.0/2*U.rho* (U.vx*U.vx + U.vy*U.vy + U.vz*U.vz)
                        - 1.0/2*(U.Bx*U.Bx + U.By*U.By + U.Bz*U.Bz));

    double pdx =  gamma * (DX.E - 1.0/2*DX.rho* (U.vx*U.vx + U.vy*U.vy + U.vz*U.vz)
                           - U.rho*(DX.vx + DX.vy + DX.vz)
                           - (DX.Bx + DX.By + DX.Bz)) ;

    double pdy =  gamma * (DY.E - 1.0/2*DY.rho* (U.vx*U.vx + U.vy*U.vy + U.vz*U.vz)
                           - U.rho*(DY.vx + DY.vy + DY.vz)
                           - (DY.Bx + DY.By + DY.Bz)) ;

//dp/dt = p*Vx/dx + p/dx * Vx + p*Vy/dy + p/dy*Vy
    double drhodt = U.rho * DX.vx + DX.rho * U.vx + U.rho * DY.vy + DY.rho * U.vy;
    //p * Vx/dt = p/dx*Vx*Vx + 2*p*Vx/dx*Vx +
    //			p/dy*Vx*Vy + p*Vx/dy*Vy   + p*Vx*Vy/dy
    //			- P/dx
    //			-2*Bx*Bx/dx
    //			-Bx/dy*By - Bx*By/dy
    //			-p/dt * Vx
    double Pdx = (pdx + (DX.Bx + DX.By + DX.Bz)) / U.rho;
    //Pdx = 0; //10*DX.rho;

    double dvxdt = DX.rho*U.vx*U.vx + 2*U.rho*DX.vx*U.vx +
                   DY.rho*U.vx*U.vy + U.rho * DY.vx * U.vy + U.rho * U.vx * DY.vy
                   +Pdx
                   -2 * U.Bx*DX.Bx
                   -DY.Bx*U.By - U.Bx*DY.By
                   -drhodt * U.vx;
    //double pdvxdt = U.vx * DX.vx + U.vy * DY.vx + Pdx;
    //p * Vy/dt = p/dx*Vx*Vy + p*Vx/dx*Vy   + p*Vx*Vy/dx
    //			p/dy*Vy*Vy + 2*p*Vy/dy*Vy +
    //			-P/dy
    //			-Bx/dx*By - Bx*By/dx
    //			-2*By*By/dy
    //			-p/dt * Vy
    double Pdy = (pdy + (DY.Bx + DY.By + DY.Bz)) / U.rho;
    //Pdy = 0; //10*DY.rho;

    double dvydt = DX.rho*U.vx*U.vy + U.rho * DX.vx * U.vy + U.rho * U.vx * DX.vy +
                   DY.rho*U.vy*U.vy + 2*U.rho*DY.vy*U.vy +
                   +Pdy
                   -DX.Bx*U.By - U.Bx * DX.Bx
                   -2 * U.By * DY.By
                   -drhodt * U.vy;

    double dvzdt = DX.rho*U.vx*U.vz + U.rho * DX.vx * U.vz + U.rho * U.vx * DX.vz +
                   DY.rho*U.vy*U.vz + U.rho * DY.vy * U.vz + U.rho * U.vy * DY.vz +
                   -DX.By*U.Bz - U.Bx * DX.Bz
                   -DY.By*U.Bz - U.By * DY.Bz
                   -drhodt * U.vy;

    //B/dt = B*vx/dx + B/dx*vx + B*vy/dy + B/dy*vy
    double dBxdt =  U.Bx * DY.vy + DY.Bx * U.vy
                    -U.By * DY.vx - DY.By * U.vx;

    double dBydt =  U.By * DX.vx + DX.By * U.vx
                    -U.Bx * DX.vy - DX.Bx * U.vy;

    double dBzdt =  U.Bz * DX.vx + DX.Bz * U.vx
                    -U.Bx * DX.vz - DX.Bx * U.vz

                    +U.Bz * DY.vy + DY.Bz * U.vy
                    -U.By * DY.vz - DY.By * U.vz;

    //E*vx / dx = E/dx*vx + E*vx/dx
    //BxByVy/dx = (Bx)'*By*Vy + Bx*By'*Vy + Bx*By*Vy'
    // p* = 10*U.rho - (DX.Bx + DX.By + DX.Bz), p*dt = dPdt
    double P = p + 1.0/2*(U.Bx*U.Bx + U.By*U.By + U.Bz*U.Bz);

    double dExdt = DX.E*U.vx + U.E*DX.vx - Pdx*U.rho*U.vx + P*DX.vx
                   - 2*DX.Bx*U.vx - U.Bx*U.Bx*DX.vx
                   - DX.Bx*U.By*U.vy - U.Bx*DX.By*U.vy - U.Bx*U.By*DX.vy
                   - DX.Bx*U.Bz*U.vz - U.Bx*DX.Bz*U.vz - U.Bx*U.Bz*DX.vz;

    double dEydt = DY.E*U.vy + U.E*DY.vy - Pdy*U.rho*U.vy + P*DY.vy
                   - 2*DY.By*U.vy - U.By*U.By*DY.vy
                   - DY.Bx*U.By*U.vx - U.Bx*DY.By*U.vx - U.Bx*U.By*DY.vx
                   - DY.By*U.Bz*U.vz - U.By*DY.Bz*U.vz - U.By*U.Bz*DY.vz;

    double dEdt = dExdt + dEydt;
    dEdt = 0;

    if(p > 1000){
        std::cout << "P: " << P << std::endl;
        std::cout << "UE: " << U.E << std::endl;
        std::cout << "p: " << p << std::endl;
    }


    return Cell{drhodt,
                dvxdt ,
                dvydt ,
                dvzdt,
                dBxdt,
                dBydt,
                dBzdt,
                dEdt};
}


void CalculateFlux(Grid& out, Grid& in)
{
    for(int y = 0; y < in.sizeY; y++){
        unsigned int y_1= (y-1+in.sizeY)%(int)in.sizeY;
        unsigned int y1= (y+1+in.sizeY)%(int)in.sizeY;
        unsigned int y_2= (y-2+in.sizeY)%(int)in.sizeY;
        for (int x = 0; x < in.sizeX; x++) // 2->n ??
        {
            //r_i= (u{i} - u{i-1}) / (u{i+1}-u{i})   r_P = (P - E) / (W - P)
            unsigned int x_1= (x-1+in.sizeX)%(int)in.sizeX;
            unsigned int x1= (x+1+in.sizeX)%(int)in.sizeX;
            unsigned int x_2= (x-2+in.sizeX)%(int)in.sizeX;

            auto rx_current = (in.mesh[x][y]     - in.mesh[x_1][y]) / nonZeroDenom( in.mesh[x1][y]- in.mesh[x][y]);
            auto rx_prev =    (in.mesh[x_1][y] - in.mesh[x_2][y]) / nonZeroDenom(in.mesh[x][y] - in.mesh[x_1][y]);

            auto ry_current = (in.mesh[x][y]   - in.mesh[x][y_1]) / nonZeroDenom(in.mesh[x][y1] - in.mesh[x][y]);
            auto ry_prev =    (in.mesh[x][y_1] - in.mesh[x][y_2]) / nonZeroDenom(in.mesh[x][y] - in.mesh[x][y_1]);
            /* uL{i-0.5} = u{i-1} + 0.5 * phi(r{i-1}) * (u{i} - u{i - 1})
             uL{P} = E + 0.5 * phi * (P - E)
             uR{i-0.5} = u{i} - 0.5 * phi(r{i}) * (u{i+1} - u{i})
             uR{P} = P - 0.5 * phi * (W - P)
             phi - slope limiter*/

            auto uL_x = in.mesh[x_1][y] + 0.5 * SlopeLim(rx_prev) * (in.mesh[x][y] - in.mesh[x_1][y]);
            auto uR_x = in.mesh[x][y]   - 0.5 * SlopeLim(rx_current)   * (in.mesh[x1][y] - in.mesh[x][y]);

            auto uL_y = in.mesh[x][y_1] + 0.5 * SlopeLim(ry_prev) * (in.mesh[x][y] - in.mesh[x][y_1]);
            auto uR_y = in.mesh[x][y]   - 0.5 * SlopeLim(ry_current)   * (in.mesh[x][y1] - in.mesh[x][y]);

            //eigenvalues? -> max speed a
            double a = 500;
            // F{i-0.5} = 0.5 * (F(uR{i-0.5}) + F(uL{i-0.5}) - a * (uR{i-0.5} - uL{i-0.5}))
            out.mesh[x][y] =     (0.5 * (F(uR_x, Cell::zeros(),in.mesh[x][y]) + F(uL_x, Cell::zeros(),in.mesh[x][y]) - a * (uR_x - uL_x)));
            out.fluxMesh[x][y] = (0.5 * (F(Cell::zeros(), uR_y,in.mesh[x][y]) + F(Cell::zeros(), uL_y,in.mesh[x][y]) - a * (uR_y - uL_y)));

        }

    }

}

Cell S(int x,int y,Cell val)
{
    if (x>3 && x<7 && y>32 && y< 58)
        return {0,0,0,0};
    else return {0,0,0,0};
}

void ApplyBoundaryConditions(Grid& grid)
{
    for (int x=40;x<50;x++)
    {
        for (int y=15; y< 30;y++)
        {
            grid.mesh[x][y].rho=5;
            grid.mesh[x][y].vx=10;
            grid.mesh[x][y].vy=0;
            grid.mesh[x][y].Bz=0;
           // grid.mesh[x][y].E=10;
        }
    }
    for (int y=0;y<grid.sizeY;y++) {
        grid.mesh[0][y].rho=0.01;
        grid.mesh[0][y].vx=0;
        grid.mesh[0][y].vy=0;
        grid.mesh[0][y].vz=0;
        grid.mesh[0][y].Bx=0;
        grid.mesh[0][y].By=0;
        grid.mesh[0][y].Bz=0;
        grid.mesh[0][y].E = 10;
        grid.mesh[grid.sizeX-1][y].rho=0.01;
        grid.mesh[grid.sizeX-1][y].vx=0;
        grid.mesh[grid.sizeX-1][y].vy=0;
        grid.mesh[grid.sizeX-1][y].vz=0;
        grid.mesh[grid.sizeX-1][y].Bx=0;
        grid.mesh[grid.sizeX-1][y].By=0;
        grid.mesh[grid.sizeX-1][y].Bz=0;
        grid.mesh[grid.sizeX-1][y].E = 10;
    }
}



void RKIntegrator(Grid& grid, double dt)
{
    double dx = 1;
    Grid flux(grid.sizeX,grid.sizeY);
    flux.fluxMeshInit();
    CalculateFlux(flux,grid);
    Grid k(grid.sizeX,grid.sizeY);
    for (int y = 0; y < grid.sizeY; y++)
    {
        unsigned int y1= (y+1+grid.sizeY)%(int)grid.sizeY;
        for (int x = 0; x < grid.sizeX; x++)
        {
            unsigned int x1= (x+1+grid.sizeX)%(int)grid.sizeX;
            k.mesh[x][y] = grid.mesh[x][y] - 0.5 * dt / dx * (flux.mesh[x1][y]+ flux.fluxMesh[x][y1] - flux.mesh[x][y] - flux.fluxMesh[x][y]) +0.5*dt * S(x,y,grid.mesh[x][y]);
        }
    }
    ApplyBoundaryConditions(k);
    CalculateFlux(flux, k);
    for (int y = 0; y < grid.sizeY; y++)
    {
        unsigned int y1= (y+1+grid.sizeY)%(int)grid.sizeY;
        for (int x = 0; x < grid.sizeX; x++)
        {
            unsigned int x1= (x+1+grid.sizeX)%(int)grid.sizeX;
            grid.mesh[x][y] = grid.mesh[x][y] - dt / dx * (flux.mesh[x1][y]+ flux.fluxMesh[x][y1] - flux.mesh[x][y] - flux.fluxMesh[x][y]) + dt * S(x,y,grid.mesh[x][y]);
        }
    }
    ApplyBoundaryConditions(grid);

}
