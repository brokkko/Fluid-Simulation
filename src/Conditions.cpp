#include "Conditions.h"
#include "Constants.h"

void InitialConditions(SphericalGrid& grid) {
    for(int th=0;th<grid.getSizeTheta();th++) {
        for (int x = 0; x < grid.getSizeR(); x++) {
            for (int y = 0; y < grid.getSizePhi(); y++) {
                double rho = small_rho / (grid.getRFromIndex(x) * grid.getRFromIndex(x)) *
                             (grid.getRFromIndex(0) * grid.getRFromIndex(0));
                double vx = 0;
                double vy = 0;
                double vz = 0;
                double Bx = 0;
                double By = 0;
                double Bz = 0;
                double T = 10;
                double E = 2 * rho * m_div_k * T / (gamma - 1)
                           + rho * (vx * vx + vy * vy + vz * vz) / 2
                           + (Bx * Bx + By * By + Bz * Bz) / 2;//(2);
                double p = (gamma - 1) *
                           (E - 0.5 * rho * (vx * vx + vy * vy + vz * vz) - 0.5 * (Bx * Bx + By * By + Bz * Bz));
                Cell &c = grid.getCellRef(x, y, th);
                c.p.rho = rho;
                c.p.Vr = vx;
                c.p.Vph = vy;
                c.p.Vth = vz;
                c.p.Br = Bx;
                c.p.Bph = By;
                c.p.Bth = Bz;
                c.p.P = p;

            }
        }
    }
    grid.UpdateCons();

}


void ApplyBoundaryConditions(SphericalGrid& grid, double t, double* dens,double* vels, double *temperature, double *magneticField){
    int r = (int) ((t/SOLAR_ROTATION)*180);
    for(int th=0;th<grid.getSizeTheta();th++) {
        for (int x = 0; x < grid.getSizePhi(); x++) {
            int row = (int)(180 + r - (double)x*180/grid.getSizePhi()) % 180;
            int col = (int)((double)th*60/grid.getSizeTheta()) % 60;
            double vx = vels[60*row + col];
            double vy = 0;
            double vz = 0;
            double Bx = magneticField[60*row + col];
            double By = 0.000;
            double Bz = 0.000;
            double T = temperature[60*row + col];
            double rho = dens[60*row + col];

            double E = 2 * rho * m_div_k * T / (gamma - 1)
                       + rho * (vx * vx + vy * vy + vz * vz) / 2
                       + (Bx * Bx + By * By + Bz * Bz) / 2;
            double p = (gamma - 1) * (E - 0.5 * rho * (vx * vx + vy * vy + vz * vz) - 0.5 * (Bx * Bx + By * By + Bz * Bz));
            Cell &c = grid.getCellRef(0, x, th);
            c.p.rho = rho;
            c.p.Vr = vx;
            c.p.Vph = vy;
            c.p.Vth = vz;
            c.p.Br = Bx;
            c.p.Bph = By;
            c.p.Bth = Bz;
            c.p.P = p;
        }
    }
    grid.UpdateCons();
}