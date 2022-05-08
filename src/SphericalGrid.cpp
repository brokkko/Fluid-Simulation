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

    for (int theta = 0; theta < sizeTheta; theta++) {
        for (int r = 0; r < sizeR; r++) {
            for (int phi = 0; phi < sizePhi; phi++) {

                double _r = minRadius + (maxRadius-minRadius)/sizeR * (r+0.5);
                double dr =(maxRadius-minRadius)/sizeR;
                double dphi = 2*M_PI/sizePhi;
                double dtheta = 2*maxPolarAngle/sizePhi;
                double _phi = ((double)phi+0.5)/sizePhi*2*M_PI;
                double _theta = M_PI_2-maxPolarAngle+((double)theta+0.5)/sizePhi*2*maxPolarAngle;
                double vol = _r*_r*std::sin(_theta)*dr*dphi*dtheta;
                mesh[phi + r*sizePhi + theta*sizePhi*sizeR] = Cell(vol,_r,_phi,_theta);
            }
        }
    }

}

void SphericalGrid::UpdatePrim()
{
    for(int i=0;i<sizeR*sizePhi*sizeTheta;i++)
    {
        mesh[i].UpdatePrim();
    }
}
void SphericalGrid::UpdateCons()
{
    for(int i=0;i<sizeR*sizePhi*sizeTheta;i++)
    {
        mesh[i].UpdateCons();
    }
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
    mesh[i] = Cell();
    mesh[i].p_rho= double(i) / (sizeR * sizePhi * sizeTheta) * v;
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
