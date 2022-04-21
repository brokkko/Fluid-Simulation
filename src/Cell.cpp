#include "Cell.h"

Cell Cell::zeros() {
    return Cell{ 0,0,0,0,0,0,0,0};
}

Cell Cell::operator+(const Cell r) const {
    return Cell{ rho + r.rho ,
                 Vr + r.Vr,
                 Vphi + r.Vphi,
                 Vtheta + r.Vphi,
                 Br + r.Br,
                 Bphi + r.Bphi,
                 Btheta + r.Btheta,
                 E+r.E,};
}

Cell Cell::operator-(const Cell r) const {
    return Cell{ rho - r.rho ,
                 Vr - r.Vr,
                 Vphi - r.Vphi,
                 Vtheta - r.Vphi,
                 Br - r.Br,
                 Bphi - r.Bphi,
                 Btheta - r.Btheta,
                 E - r.E,};
}

Cell Cell::operator*(const Cell r) const {
    return Cell{ rho * r.rho ,
                 Vr * r.Vr,
                 Vphi * r.Vphi,
                 Vtheta * r.Vphi,
                 Br * r.Br,
                 Bphi * r.Bphi,
                 Btheta * r.Btheta,
                 E * r.E,};
}

Cell Cell::operator/(const Cell r) const {
    return Cell{ rho / r.rho ,
                 Vr / r.Vr,
                 Vphi / r.Vphi,
                 Vtheta / r.Vphi,
                 Br / r.Br,
                 Bphi / r.Bphi,
                 Btheta / r.Btheta,
                 E / r.E,};
}

Cell Cell::operator*(const double r) const {
    return Cell{ rho * r ,
                 Vr * r,
                 Vphi * r,
                 Vtheta * r,
                 Br * r,
                 Bphi * r,
                 Btheta * r,
                 E * r};
}

Cell Cell::operator/(const double r) const {
    return Cell{ rho / r ,
                 Vr / r,
                 Vphi / r,
                 Vtheta / r,
                 Br / r,
                 Bphi / r,
                 Btheta / r,
                 E / r};
}

Cell operator*(const double l, const Cell r) {
    return Cell{ l*r.rho ,
                 l*r.Vr ,
                 l*r.Vphi ,
                 l*r.Vtheta ,
                 l*r.Br ,
                 l*r.Bphi ,
                 l*r.Btheta ,
                 l*r.E};
}
