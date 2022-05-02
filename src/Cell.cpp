#include "Cell.h"
#include "Constants.h"

Cell& Cell::zeros() {
    auto p = reinterpret_cast<double*>(this);
    for (int i=0;i<16;i++)
        p[i] = 0;
    return *this;
}

void Cell::UpdatePrim()
{
    p_rho = c_rho/volume;
    p_Vr=c_Mr/p_rho/volume;
    p_Vph=c_Mph/p_rho/volume;
    p_Vth=c_Mth/p_rho/volume;
    p_Br=c_Br/volume;
    p_Bph=c_Bph/volume;
    p_Bth=c_Bth/volume;
    double V2 = p_Vr * p_Vr + p_Vph * p_Vph + p_Vth * p_Vth;
    double B2 = p_Br * p_Br + p_Bph * p_Bph + p_Bth * p_Bth;
    p_P= (gamma - 1) * (c_E/volume - 0.5 * p_rho * V2 - 0.5 / mu * B2) + 0.5 / mu * B2;
    if(p_P<0) {
        std::cout << "neg P\n";
    }
}
void Cell::UpdateCons()
{
    c_rho = p_rho * volume;
    c_Mr = p_rho * p_Vr * volume;
    c_Mph = p_rho * p_Vph * volume;
    c_Mth = p_rho * p_Vth * volume;
    c_Br = p_Br * volume;
    c_Bph = p_Bph * volume;
    c_Bth = p_Bth * volume;
    double V2 = p_Vr * p_Vr + p_Vph * p_Vph + p_Vth * p_Vth;
    double B2 = p_Br * p_Br + p_Bph * p_Bph + p_Bth * p_Bth;
    double p = p_P - 0.5 * B2 / mu;
    c_E = (p / (gamma - 1) + 0.5 * p_rho * V2 + 0.5*B2/mu) * volume;
    if(c_E<0) std::cout<<"neg E\n";
}


Cell Cell::operator+(const Cell right) const {
    auto right_p= reinterpret_cast<const double*>(&right);
    Cell res = *this;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<16;i++)
        res_p[i] += right_p[i];
    return res;
}

Cell Cell::operator-(const Cell right) const {
    auto right_p= reinterpret_cast<const double*>(&right);
    Cell res = *this;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<16;i++)
        res_p[i] -= right_p[i];
    return res;
}

Cell Cell::operator*(const Cell right) const {
    auto right_p= reinterpret_cast<const double*>(&right);
    Cell res = *this;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<16;i++)
        res_p[i] *= right_p[i];
    return res;
}

Cell Cell::operator/(const Cell right) const {
    auto right_p= reinterpret_cast<const double*>(&right);
    Cell res = *this;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<16;i++)
        res_p[i] /= right_p[i];
    return res;
}

Cell Cell::operator*(const double right) const {
    Cell res = *this;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<16;i++)
        res_p[i] *= right;
    return res;
}

Cell Cell::operator/(const double right) const {
    Cell res = *this;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<16;i++)
        res_p[i] /= right;
    return res;
}

Cell operator*(const double l, const Cell r) {
    auto right_p= reinterpret_cast<const double*>(&r);
    Cell res = r;
    auto res_p = reinterpret_cast<double*>(&res);
    for (int i=0;i<16;i++)
        res_p[i] = l * right_p[i];
    return res;
}
