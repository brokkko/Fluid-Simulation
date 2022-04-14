#include "Cell.h"

Cell Cell::zeros() {
    return Cell{ 0,0,0,0,0,0,0,0};
}

Cell Cell::operator+(const Cell r) const {
    return Cell{ rho + r.rho ,
                 rhoVr + r.rhoVr,
                 rhoVphi + r.rhoVphi,
                 rhoVtheta + r.rhoVphi,
                 Br + r.Br,
                 Bphi + r.Bphi,
                 Btheta + r.Btheta,
                 E+r.E,};
}

Cell Cell::operator-(const Cell r) const {
    return Cell{ rho - r.rho ,
                 rhoVr - r.rhoVr,
                 rhoVphi - r.rhoVphi,
                 rhoVtheta - r.rhoVphi,
                 Br - r.Br,
                 Bphi - r.Bphi,
                 Btheta - r.Btheta,
                 E - r.E,};
}

Cell Cell::operator*(const Cell r) const {
    return Cell{ rho * r.rho ,
                 rhoVr * r.rhoVr,
                 rhoVphi * r.rhoVphi,
                 rhoVtheta * r.rhoVphi,
                 Br * r.Br,
                 Bphi * r.Bphi,
                 Btheta * r.Btheta,
                 E * r.E,};
}

Cell Cell::operator/(const Cell r) const {
    return Cell{ rho / r.rho ,
                 rhoVr / r.rhoVr,
                 rhoVphi / r.rhoVphi,
                 rhoVtheta / r.rhoVphi,
                 Br / r.Br,
                 Bphi / r.Bphi,
                 Btheta / r.Btheta,
                 E / r.E,};
}

Cell Cell::operator*(const double r) const {
    return Cell{ rho * r ,
                 rhoVr * r,
                 rhoVphi * r,
                 rhoVtheta * r,
                 Br * r,
                 Bphi * r,
                 Btheta * r,
                 E * r};
}

Cell Cell::operator/(const double r) const {
    return Cell{ rho / r ,
                 rhoVr / r,
                 rhoVphi / r,
                 rhoVtheta / r,
                 Br / r,
                 Bphi / r,
                 Btheta / r,
                 E / r};
}

Cell operator*(const double l, const Cell r) {
    return Cell{ l*r.rho ,
                 l*r.rhoVr ,
                 l*r.rhoVphi ,
                 l*r.rhoVtheta ,
                 l*r.Br ,
                 l*r.Bphi ,
                 l*r.Btheta ,
                 l*r.E};
}
