
#include "Vector.h"
Vector operator+(Vector a, Vector b)
{
    Vector c;
    c.r = a.r + b.r;
    c.ph = a.ph + b.ph;
    c.th = a.th + b.th;
    return c;
}

Vector operator-(Vector a, Vector b)
{
    Vector c;
    c.r = a.r - b.r;
    c.ph = a.ph - b.ph;
    c.th = a.th - b.th;
    return c;
}

Vector operator*(double a, Vector b)
{
    Vector c;
    c.r = a * b.r;
    c.ph = a * b.ph;
    c.th = a * b.th;
    return c;
}

Vector operator/(Vector a, double b)
{
    Vector c;
    c.r = a.r / b;
    c.ph = a.ph / b;
    c.th = a.th / b;
    return c;
}

Vector operator*(Vector a, Vector b)
{
    Vector c;
    c.r = a.r * b.r;
    c.ph = a.ph * b.ph;
    c.th = a.th * b.th;
    return c;
}