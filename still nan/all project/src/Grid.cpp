#include "Grid.h"

Grid::Grid(unsigned int sizeX, unsigned int sizeY) {
    this->sizeX = sizeX;
    this->sizeY = sizeY;
    mesh = new Cell*[sizeX+1];
    for (int i=0; i<sizeX+1; i++){
        mesh[i] = new Cell[sizeY+1];
    }
    for (int i=0; i<sizeX+1; i++){
        mesh[i][0] = {0,0,0,0,0,0,0,0};
    }
    for (int i=0; i<sizeY+1; i++){
        mesh[0][i] = {0,0,0,0,0,0,0,0};
    }
    fluxMesh= nullptr;
}

void Grid::fluxMeshInit() {
    fluxMesh = new Cell*[sizeX+1];
    for (int i=0; i<sizeX+1; i++){
        fluxMesh[i] = new Cell[sizeY+1];
    }
    for (int i=0; i<sizeX+1; i++){
        fluxMesh[i][0] = {0,0,0,0,0, 0,0,0};
    }
    for (int i=0; i<sizeY+1; i++){
        fluxMesh[0][i] = {0,0,0,0,0,0,0,0};
    }
}

void Grid::Fill(double v) {
    for (int i = 0; i < sizeX; i++) {
        for (int j = 0; j < sizeY; j++) {
            mesh[i][j] = {v, 0, 0,0 ,0,0,0,200};
        }
    }
}

Grid::~Grid() {
    for (int i = 0; i < sizeX + 1; i++) {
        delete[] mesh[i];
    }
    delete[] mesh;
    if (fluxMesh) {
        for (int i = 0; i < sizeX + 1; i++) {
            delete[] fluxMesh[i];
        }
        delete[] fluxMesh;
    }
}
