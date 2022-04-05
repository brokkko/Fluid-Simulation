//
// Created by brokkko on 05.03.2022.
//

#include <sstream>
#include "ReadData.h"

ReadData::ReadData(const char *filePath) {
    n2 = 0, n3 = 0;
    nc_open(filePath, NC_NOWRITE, &this->fileID);
    getDataSize();
}

void ReadData::getDataSize() {
    int dimensionID;
    nc_inq_dimid(this->fileID, "n2", &dimensionID);
    nc_inq_dimlen(this->fileID, dimensionID, &this->n2);

    nc_inq_dimid(this->fileID, "n3", &dimensionID);
    nc_inq_dimlen(this->fileID, dimensionID, &this->n3);

    std::cout << n2 << " " << n3 << std::endl;
}

void ReadData::getData(const char *name, double *dataArray) const {
    int id;
    nc_inq_varid(this->fileID, name, &id);
    nc_get_var_double(this->fileID, id, dataArray);
}

void ReadData::readData(const char *name, double ** dataArray, int* dataArraySize) {
    *dataArray = new double [n3*n2];
    this->getData(name, *dataArray);
    *dataArraySize = n2*n3;
}

ReadData::~ReadData() {
    nc_close(this->fileID);
}




