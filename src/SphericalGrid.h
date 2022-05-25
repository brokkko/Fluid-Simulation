//
// Created by alex on 12.04.2022.
//

#ifndef FLUID_SIMULATION_SPHERICALGRID_H
#define FLUID_SIMULATION_SPHERICALGRID_H
#include "Cell.h"
#include <cmath>

class SphericalGrid {
    Cell* mesh;
    unsigned int sizeR;
    unsigned int sizePhi;
    unsigned int sizeTheta;
    double maxPolarAngle;
    double minRadius;
    double maxRadius;

public:
    SphericalGrid(unsigned int sizeR,unsigned int sizePhi,unsigned int sizeTheta,double minRadius,double maxRadius, double maxPolarAngle);
    Cell getCell(int R,int Phi,int Theta);
    Cell& getCellRef(int R,int Phi,int Theta);
    unsigned int getSizeR() const {return sizeR;}
    unsigned int getSizePhi() const {return sizePhi;}
    unsigned int getSizeTheta() const {return sizeTheta;}
    double getRFromIndex(unsigned int R) const {return minRadius + (maxRadius-minRadius)/sizeR * (R+0.5);}
    double getPhiFromIndex(unsigned int Phi) const {return ((double)Phi+0.5)/sizePhi*2*M_PI;}
    double getThetaFromIndex(unsigned int Theta) const {return M_PI_2-maxPolarAngle + (2 * Theta + 1)*(maxPolarAngle/sizeTheta);}
    void UpdatePrim();
    void UpdateCons();
    static SphericalGrid copyGrid(SphericalGrid& grid);
    void Fill(double v);
    ~SphericalGrid();
    double getMaxRadius() const {return maxRadius;}
    double getMinRadius() const {return minRadius;}
};


#endif //FLUID_SIMULATION_SPHERICALGRID_H
