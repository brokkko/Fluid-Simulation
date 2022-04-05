//
// Created by brokkko on 05.03.2022.
//

#ifndef INC_1_READDATA_H
#define INC_1_READDATA_H

#include <iostream>
#include <netcdf.h>

class ReadData {
private:
    int fileID;
    size_t n2, n3;
    void getDataSize();
    void getData(const char *name, double *dataArray) const;

public:
    ReadData(const char *filePath);
    void readData(const char *name, double ** dataArray, int* dataArraySize);
    ~ReadData();
};


#endif //INC_1_READDATA_H
