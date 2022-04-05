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
    Cell DZ=Cell::zeros();

    double gamma = 5.0/3 - 1;
    double mu =1.2566e-8;

    double p = gamma * (U.E - 1.0/2 * U.rho* (U.vx * U.vx + U.vy * U.vy + U.vz * U.vz)
                        - 1.0/(2*mu)*(U.Bx*U.Bx + U.By*U.By + U.Bz * U.Bz));

    double pdx =  gamma * (DX.E - DX.rho * (U.vx * U.vx + U.vy * U.vy + U.vz * U.vz) / 2
                           - U.rho * (DX.vx * U.vx + DX.vy * U.vy + DX.vz * U.vz)
                           - (DX.Bx * U.Bx + DX.By * U.By + DX.Bz * U.Bz)/mu) ;

    double pdy =  gamma * (DY.E - DY.rho* (U.vx*U.vx + U.vy*U.vy + U.vz*U.vz) / 2
                           - U.rho * (DY.vx * U.vx + DY.vy * U.vy + DY.vz * U.vz)
                           - (DY.Bx * U.Bx + DY.By * U.By + DY.Bz * U.Bz)/mu) ;

    double pdz =  gamma * (DZ.E - DZ.rho* (U.vx*U.vx + U.vy*U.vy + U.vz*U.vz) / 2
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
                   +Pdx
                   -2 * U.Bx*DX.Bx / mu
                   -(DY.Bx * U.By + U.Bx * DY.By) / mu
                   -(DZ.Bx * U.Bz + U.Bx * DZ.Bz) / mu
                   -drhodt * U.vx);// / U.rho;


    double dvydt = (DX.rho * U.vx * U.vy + U.rho * DX.vx * U.vy + U.rho * U.vx * DX.vy
                   +DY.rho * U.vy * U.vy + 2 * U.rho * DY.vy * U.vy
                   +DZ.rho * U.vz * U.vy + U.rho * DZ.vz * U.vy + U.rho * U.vz * DZ.vy
                   +Pdy
                   -(DX.Bx * U.By + U.Bx * DX.By) / mu
                   -2 * U.By * DY.By / mu
                   -(DZ.Bz * U.By + U.Bz * DZ.By) / mu
                   -drhodt * U.vy);///U.rho;

    double dvzdt = (DX.rho * U.vx * U.vz + U.rho * DX.vx * U.vz + U.rho * U.vx * DX.vz
                   +DY.rho * U.vy * U.vz + U.rho * DY.vy * U.vz + U.rho * U.vy * DY.vz
                   +DZ.rho * U.vz * U.vz + 2 * U.rho * DZ.vz * U.vz
                   +Pdz
                   -(DX.Bx * U.Bz + U.Bx * DX.Bz) / mu
                   -(DY.By * U.Bz + U.By * DY.Bz) / mu
                   -2 * U.Bz * DZ.Bz / mu
                   -drhodt * U.vz);// / U.rho;

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



    //double dExdt = DX.E*U.vx + U.E*DX.vx + Pdx*U.rho*U.vx + P*DX.vx
                   - 2*DX.Bx*U.vx - U.Bx*U.Bx*DX.vx
                   - DX.Bx*U.By*U.vy - U.Bx*DX.By*U.vy - U.Bx*U.By*DX.vy
                   - DX.Bx*U.Bz*U.vz - U.Bx*DX.Bz*U.vz - U.Bx*U.Bz*DX.vz;

    //double dEydt = DY.E*U.vy + U.E*DY.vy + Pdy*U.rho*U.vy + P*DY.vy
                   - 2*DY.By*U.vy - U.By*U.By*DY.vy
                   - DY.Bx*U.By*U.vx - U.Bx*DY.By*U.vx - U.Bx*U.By*DY.vx
                   - DY.By*U.Bz*U.vz - U.By*DY.Bz*U.vz - U.By*U.Bz*DY.vz;

    ///double dEdt = dExdt + dEydt;

    double dEdt = -(DX.E * U.vx + DY.E * U.vy + DZ.E * U.vz) + U.E * (DX.vx + DY.vy + DZ.vz)    //nabla(EV)
                  -(Pdx * U.vx + Pdy * U.vy + Pdz * U.vz) + P * (DX.vx + DY.vy + DZ.vz)         //nabla(PV)

                  + (DX.vx * U.Bx * U.Bx + 2 * DX.Bx * U.Bx * U.vx
                  + DX.vy * U.Bx * U.By + DX.Bx * U.By * U.vy + DX.By * U.Bx * U.vy             //DX
                  + DX.vz * U.Bx * U.Bz + DX.Bx * U.Bz * U.vz + DX.Bz * U.Bx * U.vz

                  + DY.By * U.Bx * U.vx + DY.Bx * U.By * U.vx  +DY.vx * U.By * U.Bx
                  + 2 * DY.By * U.By * U.vy + DY.vy * U.By * U.By
                  + DY.By * U.Bz * U.vz + DY.Bz * U.By * U.vz + DY.vz * U.By * U.Bz

                  + DZ.vx * U.Bz * U.Bx + DZ.Bz * U.Bx * U.vx + DZ.Bx * U.Bz * U.vx
                  + DZ.vy * U.Bz * U.By + DZ.Bz * U.By * U.vy + DZ.By * U.Bz * U.vy
                  + DZ.vz * U.Bz * U.Bz + 2 * DZ.Bz * U.Bz * U.vz)/mu;


    double dExdt = DX.E*U.vx + U.E*DX.vx + Pdx*U.rho*U.vx + P*DX.vx
                   - 2*DX.Bx*U.vx - U.Bx*U.Bx*DX.vx
                   - DX.Bx*U.By*U.vy - U.Bx*DX.By*U.vy - U.Bx*U.By*DX.vy
                   - DX.Bx*U.Bz*U.vz - U.Bx*DX.Bz*U.vz - U.Bx*U.Bz*DX.vz;

    double dEydt = DY.E*U.vy + U.E*DY.vy + Pdy*U.rho*U.vy + P*DY.vy
                   - 2*DY.By*U.vy - U.By*U.By*DY.vy
                   - DY.Bx*U.By*U.vx - U.Bx*DY.By*U.vx - U.Bx*U.By*DY.vx
                   - DY.By*U.Bz*U.vz - U.By*DY.Bz*U.vz - U.By*U.Bz*DY.vz;

    dEdt = dExdt + dEydt;

    return Cell{drhodt,
                dvxdt ,
                dvydt ,
                dvzdt,
                0,
                0,
                0,
                dEdt};
}

double maxSpeed(Cell U){
    double gamma = 5.0/3 - 1;
    double p = gamma * (U.E - 1.0/2*U.rho* (U.vx*U.vx + U.vy*U.vy + U.vz*U.vz)
                        - 1.0/2*(U.Bx*U.Bx + U.By*U.By + U.Bz*U.Bz));
    double P = p + 1.0/2*(U.Bx*U.Bx + U.By*U.By + U.Bz*U.Bz);

    double B_2 = U.Bx*U.Bx + U.By*U.By + U.Bz*U.Bz;
    double B = std::sqrt(B_2);

    double cmax =
            std::sqrt(std::abs(U.vx)*std::abs(U.vx) + std::abs(U.vy)*std::abs(U.vy) + std::abs(U.vz)*std::abs(U.vz))
            + 1.0/2*((gamma*p + B_2)/U.rho
            + std::sqrt(((gamma + B)/U.rho) * ((gamma + B)/U.rho) - 4*(gamma*U.Bx*U.Bx)/U.rho*U.rho));

    return cmax;
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
            a = maxSpeed(in.mesh[x][y]);
             a = 1000;
            //std::cout << a << std::endl;
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

void InitialConditions(Grid& grid) {
    for (int x = 0; x < grid.sizeX; x++) {
        for (int y = 0; y < grid.sizeY; y++) {
            double rho=0.01;
            double vx=0;
            double vy=0;
            double vz=0;
            double Bx=0;
            double By=0;
            double Bz=0;
            double T =273;
            double m_div_k=8249.805773;
            double mu =1.2566e-8;
            // for (int y = 45; y < 55; y++){
            //        for (int x = 20; x < 40; x++){
            if(x>20 && x< 40 && y>45 && y<55)
            {
                rho=0.02;
                vx=0;
            }


            //nk=rho/m m=1.6733e-27
            //k=1.38044e-23
            double E = 3 * rho * m_div_k * T
                    + rho * (vx*vx + vy*vy + vz*vz)
                    + (Bx*Bx + By*By + Bz*Bz) /(2*mu);

            grid.mesh[x][y] = {rho, vx, vy, vz,Bx,By,Bz,E};

        }
    }
}


void ApplyBoundaryConditions(Grid& grid)
{
   /* double rho=0.01;
    double vx=0;
    double vy=0;
    double vz=0;
    double Bx=0;
    double By=0;
    double Bz=0;
    double T =273;
    double m_div_k=8249.805773;
    double mu =1.2566e-8;
    //nk=rho/m m=1.6733e-27
    //k=1.38044e-23
    double E = 3 * rho * m_div_k * T
               + rho * (vx*vx + vy*vy + vz*vz)
               + (Bx*Bx + By*By + Bz*Bz) /(2*mu);
    for (int y=0;y<grid.sizeY;y++) {
        grid.mesh[0][y].rho=rho;
        grid.mesh[0][y].vx=0;
        grid.mesh[0][y].vy=0;
        grid.mesh[0][y].vz=0;
        grid.mesh[0][y].Bx=0;
        grid.mesh[0][y].By=0;
        grid.mesh[0][y].Bz=0;
        grid.mesh[0][y].E = E;
        grid.mesh[grid.sizeX-1][y].rho=rho;
        grid.mesh[grid.sizeX-1][y].vx=0;
        grid.mesh[grid.sizeX-1][y].vy=0;
        grid.mesh[grid.sizeX-1][y].vz=0;
        grid.mesh[grid.sizeX-1][y].Bx=0;
        grid.mesh[grid.sizeX-1][y].By=0;
        grid.mesh[grid.sizeX-1][y].Bz=0;
        grid.mesh[grid.sizeX-1][y].E = E;
    }*/
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
