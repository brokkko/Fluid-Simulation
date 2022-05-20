#include "Cell.h"
#include "Constants.h"

Cell& Cell::zeros() {
    p.zeros();
    c.zeros();
    return *this;
}

PrimitiveVector& PrimitiveVector::zeros() {
    auto p = reinterpret_cast<double*>(this);
    for (int i=0;i<8;i++)
        p[i] = 0;
    return *this;
}
ConservativeVector& ConservativeVector::zeros() {
    auto p = reinterpret_cast<double*>(this);
    for (int i=0;i<8;i++)
        p[i] = 0;
    return *this;
}

PrimitiveVector PrimitiveVector::rotate(double phi,double theta)
{
    PrimitiveVector res = *this;
    res.Vr = Vr*std::cos(phi) - Vph*std::sin(phi);
    res.Vph =Vr*std::sin(phi) + Vph*std::cos(phi);

    res.Br = Br*std::cos(phi) - Bph*std::sin(phi);
    res.Bph =Br*std::sin(phi) + Bph*std::cos(phi);

    PrimitiveVector res2 = res;
    res2.Vr = res.Vr*std::cos(theta) - res.Vth*std::sin(theta);
    res2.Vth =res.Vr*std::sin(theta) + res.Vth*std::cos(theta);

    res2.Br = res.Br*std::cos(theta) - res.Bth*std::sin(theta);
    res2.Bth =res.Br*std::sin(theta) + res.Bth*std::cos(theta);

    return res2;
}
ConservativeVector ConservativeVector::rotate(double phi,double theta)
{
    ConservativeVector res = *this;
    res.Mr = Mr*std::cos(phi) - Mph*std::sin(phi);
    res.Mph =Mr*std::sin(phi) + Mph*std::cos(phi);

    res.Br = Br*std::cos(phi) - Bph*std::sin(phi);
    res.Bph =Br*std::sin(phi) + Bph*std::cos(phi);

    ConservativeVector res2 = res;
    res2.Mr =  res.Mr*std::cos(theta) - res.Mth*std::sin(theta);
    res2.Mth = res.Mr*std::sin(theta) + res.Mth*std::cos(theta);

    res2.Br =  res.Br*std::cos(theta) - res.Bth*std::sin(theta);
    res2.Bth = res.Br*std::sin(theta) + res.Bth*std::cos(theta);

    return res2;
}


int Cell::UpdatePrim()
{
    int res=0;
    if(c.m<0)
    {
#if defined(PRINT_NEG)
        std::cout<<"neg Rho\n";
#endif
        c.m=small_rho;//*volume;
        res=1;
    }
    p.rho = c.m;///volume;

    p.Vr=c.Mr/p.rho;///volume;
    p.Vph=c.Mph/p.rho;///volume;
    p.Vth=c.Mth/p.rho;///volume;
    p.Br=c.Br;///volume;
    p.Bph=c.Bph;///volume;
    p.Bth=c.Bth;///volume;
    double V2 = p.Vr * p.Vr + p.Vph * p.Vph + p.Vth * p.Vth;
    double M2 = c.Mr * c.Mr + c.Mph * c.Mph + c.Mth * c.Mth;
    double B2 = p.Br * p.Br + p.Bph * p.Bph + p.Bth * p.Bth;
    double cB2 = c.Br * c.Br + c.Bph * c.Bph + c.Bth * c.Bth;

    if(c.E<0) {
#if defined(PRINT_NEG)
        std::cout<<"neg E\n";
#endif
        c.E = small_P/(gamma-1) + 0.5*(M2/c.m + cB2/mu);
        res=2;
    }
    p.P= (gamma - 1) * (c.E/*/volume*/ - 0.5 * p.rho * V2 - 0.5 * B2/mu);// + 0.5 / mu * B2;
    if(p.P<0) {
#if defined(PRINT_NEG)
        std::cout << "neg P\n";
#endif
        p.P = small_P;
        c.E = (p.P / (gamma - 1) + 0.5 * p.rho * V2 + 0.5 * B2/mu);// * volume;
        res=3;
    }
    return res;
}
void Cell::UpdateCons()
{
    c.m = p.rho ;//* volume;
    c.Mr = p.rho * p.Vr ;//* volume;
    c.Mph = p.rho * p.Vph ;//* volume;
    c.Mth = p.rho * p.Vth ;//* volume;
    c.Br = p.Br ;//* volume;
    c.Bph = p.Bph ;//* volume;
    c.Bth = p.Bth ;//* volume;
    double V2 = p.Vr * p.Vr + p.Vph * p.Vph + p.Vth * p.Vth;
    //double M2 = c.Mr * c.Mr + c.Mph * c.Mph + c.Mth * c.Mth;
    double B2 = p.Br * p.Br + p.Bph * p.Bph + p.Bth * p.Bth;
    //double cB2 = c.Br * c.Br + c.Bph * c.Bph + c.Bth * c.Bth;
    //double prs = p.P - 0.5 * B2 / mu;
    c.E = (p.P / (gamma - 1) + 0.5 * p.rho * V2 + 0.5 * B2/mu) ;//* volume;

}


Cell Cell::operator+(const Cell right) const {
    Cell res = *this;
    res.p = res.p + right.p;
    res.c = res.c + right.c;
    return res;
}

template<class T>
T add(T l,T r)
{
    auto right_p= reinterpret_cast<const double*>(&r);
    T res = l;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<8;i++)
        res_p[i] += right_p[i];
    return res;
}

template<class T>
T sub(T l,T r)
{
    auto right_p= reinterpret_cast<const double*>(&r);
    T res = l;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<8;i++)
        res_p[i] -= right_p[i];
    return res;
}
template<class T>
T mul(T l,T r)
{
    auto right_p= reinterpret_cast<const double*>(&r);
    T res = l;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<8;i++)
        res_p[i] *= right_p[i];
    return res;
}

template<class T>
T div(T l,T r)
{
    auto right_p= reinterpret_cast<const double*>(&r);
    T res = l;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<8;i++)
        res_p[i] /= right_p[i];
    return res;
}

template<class T>
T div(T l,double r)
{
    T res = l;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<8;i++)
        res_p[i] /= r;
    return res;
}

template<class T>
T mul(T l,double r)
{
    T res = l;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<8;i++)
        res_p[i] *= r;
    return res;
}

PrimitiveVector PrimitiveVector::operator+(const PrimitiveVector right) const
{
    return add(*this,right);
}
ConservativeVector ConservativeVector::operator+(const ConservativeVector right) const
{
    return add(*this,right);
}
PrimitiveVector PrimitiveVector::operator-(const PrimitiveVector right) const
{
    return sub(*this,right);
}
ConservativeVector ConservativeVector::operator-(const ConservativeVector right) const
{
    return sub(*this,right);
}
PrimitiveVector PrimitiveVector::operator*(const PrimitiveVector right) const
{
    return mul(*this,right);
}
ConservativeVector ConservativeVector::operator*(const ConservativeVector right) const
{
    return mul(*this,right);
}
PrimitiveVector PrimitiveVector::operator/(const PrimitiveVector right) const
{
    return div(*this,right);
}
ConservativeVector ConservativeVector::operator/(const ConservativeVector right) const
{
    return div(*this,right);
}
PrimitiveVector PrimitiveVector::operator*(const double right) const
{
    return mul(*this,right);
}
ConservativeVector ConservativeVector::operator*(const double right) const
{
    return mul(*this,right);
}
PrimitiveVector PrimitiveVector::operator/(const double right) const
{
    return div(*this,right);
}
ConservativeVector ConservativeVector::operator/(const double right) const
{
    return div(*this,right);
}
PrimitiveVector operator*(const double l, const PrimitiveVector r) {
    return mul(r,l);
}
ConservativeVector operator*(const double l, const ConservativeVector r) {
    return mul(r,l);
}


Cell Cell::operator-(const Cell right) const {
    Cell res = *this;
    res.p = res.p - right.p;
    res.c = res.c - right.c;
    return res;
}

Cell Cell::operator*(const Cell right) const {
    Cell res = *this;
    res.p = res.p * right.p;
    res.c = res.c * right.c;
    return res;
}

Cell Cell::operator/(const Cell right) const {
    Cell res = *this;
    res.p = res.p / right.p;
    res.c = res.c / right.c;
    return res;
}

Cell Cell::operator*(const double right) const {
    Cell res = *this;
    res.p = res.p * right;
    res.c = res.c * right;
    return res;
}

Cell Cell::operator/(const double right) const {
    Cell res = *this;
    res.p = res.p / right;
    res.c = res.c / right;
    return res;
}

Cell operator*(const double l, const Cell r) {
    Cell res = r;
    res.p = l * r.p;
    res.c = l * r.c;
    return res;
}
