//
// Created by alex on 12.04.2022.
//

#include "SphericalGrid.h"
#include "Constants.h"

SphericalGrid::SphericalGrid(unsigned int sizeR,unsigned int sizePhi,unsigned int sizeTheta,double minRadius,double maxRadius, double maxPolarAngle)
{
    this->sizeR=sizeR;
    this->sizeTheta=sizeTheta;
    this->sizePhi = sizePhi;
    mesh = new Cell[sizeR*sizePhi*sizeTheta];
    this->minRadius = minRadius;
    this->maxRadius = maxRadius;
    this->maxPolarAngle = maxPolarAngle;
}

Cell SphericalGrid::getCell(int R,int Phi,int Theta)
{
    R=clamp(R,0,(int)sizeR-1);
    //Phi=clamp(Phi,0,(int)sizePhi-1);
    Phi=(Phi+sizePhi)%(int)sizePhi;
    Theta=clamp(Theta,0,(int)sizeTheta-1);
    return mesh[Phi + R*sizePhi + Theta*sizePhi*sizeR];
}
Cell& SphericalGrid::getCellRef(int R,int Phi,int Theta)
{
    return mesh[Phi + R*sizePhi + Theta*sizePhi*sizeR];
}
void SphericalGrid::Fill(double v) {
    for (int i = 0; i < sizeR * sizePhi * sizeTheta; i++) {
    mesh[i] = Cell::zeros();
    mesh[i].rho=double(i)/(sizeR * sizePhi * sizeTheta)*v;
    }
}

SphericalGrid SphericalGrid::copyGrid(SphericalGrid& grid)
{
    return SphericalGrid(grid.sizeR,grid.sizePhi,grid.sizeTheta,grid.minRadius,grid.maxRadius,grid.maxPolarAngle);

}
SphericalGrid::~SphericalGrid()
{
    delete[] mesh;
}
