
#ifndef FLUID_SIMULATION_VECTOR_H
#define FLUID_SIMULATION_VECTOR_H
#include <vector>

struct Vector {
    double r;
    double ph;
    double th;

    friend Vector operator +(Vector a, Vector b);
    friend Vector operator -(Vector a, Vector b);
    friend Vector operator *(double a, Vector b);
    friend Vector operator /(Vector a, double b);
    friend Vector operator *(Vector a, Vector b);
};


#endif //FLUID_SIMULATION_VECTOR_H
