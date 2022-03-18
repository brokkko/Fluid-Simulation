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
            std::max(0.0, std::min(1.0, r.B))};

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

//one vastd::absble !!
Cell nonZeroDenom(Cell denom)
{
    return {nonZeroDouble(denom.rho),nonZeroDouble(denom.vx),nonZeroDouble(denom.vy),nonZeroDouble(denom.B)};
}

Cell F(Cell DX,Cell DY,Cell U)
{
    // double dvaldt=   dudx.rho * U.vx + U.rho * dudx.vx +
    //  dudy.rho * U.vy + U.rho * dudy.vy;
//dp/dt = p*Vx/dx + p/dx * Vx + p*Vy/dy + p/dy*Vy
    double drhodt = U.rho * DX.vx + DX.rho * U.vx + U.rho * DY.vy + DY.rho * U.vy;
    //p * Vx/dt = p/dx*Vx*Vx + 2*p*Vx/dx*Vx +
    //			p/dy*Vx*Vy + p*Vx/dy*Vy   + p*Vx*Vy/dy
    //			- P/dx
    //			-2*Bx*Bx/dx
    //			-Bx/dy*By - Bx*By/dy
    //			-p/dt * Vx
    double Pdx = 10 * DX.rho/U.rho;
    double dvxdt = DX.rho*U.vx*U.vx + 2*U.rho*DX.vx*U.vx +
                   DY.rho*U.vx*U.vy + U.rho * DY.vx * U.vy + U.rho * U.vx * DY.vy
                   +Pdx
                   -drhodt * U.vx;
    - (U.vx*U.B)+ (U.vy*U.B);
    //double pdvxdt = U.vx * DX.vx + U.vy * DY.vx + Pdx;
    //p * Vy/dt = p/dx*Vx*Vy + p*Vx/dx*Vy   + p*Vx*Vy/dx
    //			p/dy*Vy*Vy + 2*p*Vy/dy*Vy +
    //			-P/dy
    //			-Bx/dx*By - Bx*By/dx
    //			-2*By*By/dy
    //			-p/dt * Vy
    double Pdy = 10 * DY.rho/U.rho;
    double dvydt = DX.rho*U.vx*U.vy + U.rho * DX.vx * U.vy + U.rho * U.vx * DX.vy +
                   DY.rho*U.vy*U.vy + 2*U.rho*DY.vy*U.vy +
                   +Pdy
                   -drhodt * U.vy;
    + (U.vy*U.B)-(U.vx*U.B);
    //double pdvydt =U.rho *  (U.vx * DX.vy + U.vy * DY.vy) + Pdy;
    //pdvxdt = (U.vx*DX.vx + U.vy * DY.vx) + DX.rho*10 /*+ U.B*(dudx.B + dudy.B)*/;
    //pdvydt = (U.vx*DX.vy + U.vy * DY.vy) + DY.rho*10/* + U.B*(dudx.B + dudy.B)*/;
    return Cell{drhodt, dvxdt ,dvydt ,0};
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
            double a = 300;
            // F{i-0.5} = 0.5 * (F(uR{i-0.5}) + F(uL{i-0.5}) - a * (uR{i-0.5} - uL{i-0.5}))
            out.mesh[x][y] =     (0.5 * (F(uR_x, {0,0,0,0},in.mesh[x][y]) + F(uL_x, {0,0,0,0},in.mesh[x][y]) - a * (uR_x - uL_x)));
            out.fluxMesh[x][y] = (0.5 * (F({0,0,0,0}, uR_y,in.mesh[x][y]) + F({0,0,0,0}, uL_y,in.mesh[x][y]) - a * (uR_y - uL_y)));

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
    for (int x=3;x<7;x++)
    {
        for (int y=32; y< 58;y++)
        {
            grid.mesh[x][y].rho=10;
            grid.mesh[x][y].vx=10;
            grid.mesh[x][y].vy=0;
        }
    }
    for (int y=0;y<grid.sizeY;y++) {
        grid.mesh[0][y].rho=0.01;
        grid.mesh[0][y].vx=0;
        grid.mesh[0][y].vy=0;
        grid.mesh[grid.sizeX-1][y].rho=0.01;
        grid.mesh[grid.sizeX-1][y].vx=0;
        grid.mesh[grid.sizeX-1][y].vy=0;
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
