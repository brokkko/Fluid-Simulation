#include "Conditions.h"
#include "Constants.h"

double* getDensity(){
    double *data = nullptr; int dataSize = 0;
    ReadData reader(DATA);
    reader.readData("D", &data, &dataSize);
    return data;
}

double* getVelocity(){
    double *data = nullptr; int dataSize = 0;
    ReadData reader(DATA);
    reader.readData("V1", &data, &dataSize);
    return data;
}

double* getTemperature(){
    double *data = nullptr; int dataSize = 0;
    ReadData reader(DATA);
    reader.readData("T", &data, &dataSize);
    return data;
}

double* getMagneticField(){
    double *data = nullptr; int dataSize = 0;
    ReadData reader(DATA);
    reader.readData("B1", &data, &dataSize);
    return data;
}



void InitialConditions(SphericalGrid& grid, double* dens,double* vels, double *temperature, double *magneticField ) {
    for(int th=0;th<grid.getSizeTheta();th++) {
        for (int x = 0; x < grid.getSizeR(); x++) {
            for (int y = 0; y < grid.getSizePhi(); y++) {
                int row = (int)(180 - (double)x*180/grid.getSizePhi()) % 180;
                int col = (int)((double)th*60/grid.getSizeTheta()) % 60;

                double rho = small_rho/*dens[60*row + col] */ / (grid.getRFromIndex(x) * grid.getRFromIndex(x)) *
                             (grid.getRFromIndex(0) * grid.getRFromIndex(0));
                double vx =0; vels[60*row + col];
                double vy = 0;
                double vz = 0;
                double Bx =0; magneticField[60*row + col]/ (grid.getRFromIndex(x) * grid.getRFromIndex(x)) *
                            (grid.getRFromIndex(0) * grid.getRFromIndex(0));;;
                double By = 0;
                double Bz = 0;
                double T =10; temperature[60*row + col]/ (grid.getRFromIndex(x) * grid.getRFromIndex(x)) *
                           (grid.getRFromIndex(0) * grid.getRFromIndex(0));

                if(th==5 && y> 20 && y<30 && x> 5 && x <15 )
                {
                   // rho *=20;
                    //vx=100000;
                    //T=100;
                }


                double p = 2 * rho * m_div_k * T;
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
    int r = (int) ((t/SOLAR_ROTATION)*180)%180;
    for(int th=0;th<grid.getSizeTheta();th++) {
        for (int x = 0; x < grid.getSizePhi(); x++) {
            int row = (int)(180 -r + (double)x*180/grid.getSizePhi()) % 180;
            int col = (int)((double)th*60/grid.getSizeTheta()) % 60;
            double vx = vels[60*row + col];
            double vy = 0;
            double vz = 0;
            double Bx =0;// magneticField[60*row + col];
            double By = 0.000;
            double Bz = 0;magneticField[60*row + col];
            double T = temperature[60*row + col];
            double rho =dens[60*row + col];

            double p = 2 * rho * m_div_k * T;
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