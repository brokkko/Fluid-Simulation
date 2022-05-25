#include "Derivative.h"
#include "Constants.h"
#include "Vector.h"

Cell F(Cell Dr,Cell Dtheta, Cell Dphi, Cell U)
{
    double r =U.r;
    double phi = U.phi;
    double theta = U.theta;


    double drhodt =   Dr.p.rho*U.p.Vr + U.p.rho*Dr.p.Vr
                    + Dphi.p.rho*U.p.Vph + U.p.rho*Dphi.p.Vph
                    + Dtheta.p.rho*U.p.Vth + U.p.rho*Dtheta.p.Vth;

    Vector V_1 = {Dr.p.Vr * U.p.Vr + Dphi.p.Vr * U.p.Vph + Dtheta.p.Vr * U.p.Vth,
                 Dr.p.Vph * U.p.Vr + Dphi.p.Vph * U.p.Vph + Dtheta.p.Vph * U.p.Vth,
                 Dr.p.Vth * U.p.Vr + Dphi.p.Vth * U.p.Vph + Dtheta.p.Vth * U.p.Vth};

    Vector p_2 = { Dr.p.P,
                  Dphi.p.P,
                  Dtheta.p.P};

    Vector Phi = {G*Ms/(r*r), 0, 0};

    Vector A = { Dphi.p.Bth - Dtheta.p.Bph,
                -Dr.p.Bth + Dtheta.p.Br,
                Dr.p.Bph - Dphi.p.Br};

    Vector B_3 = { A.ph*U.p.Bth - A.th*U.p.Bph,
                  -A.r*U.p.Bth - A.th*U.p.Br,
                  A.r*U.p.Bph - A.th*U.p.Br};

    Vector V = V_1 - p_2/(-U.p.rho) + Phi - 1.0/(U.p.rho*mu)*B_3;

    // ------- B --------
    Vector DBdr = { (Dr.p.Vph*U.p.Bth + U.p.Vph*Dr.p.Bth) - (Dr.p.Vth*U.p.Bph + U.p.Vth*Dr.p.Bph),
                   - (Dr.p.Vr*U.p.Bth + U.p.Vr*Dr.p.Bth) + (Dr.p.Vth*U.p.Br + U.p.Vth*Dr.p.Br),
                   (Dr.p.Vr*U.p.Bph + U.p.Vr*Dr.p.Bph) - (Dr.p.Vph*U.p.Br + U.p.Vph*Dr.p.Br)};
    Vector DBdph = { (Dphi.p.Vph*U.p.Bth + U.p.Vph*Dphi.p.Bth) - (Dphi.p.Vth*U.p.Bph + U.p.Vth*Dphi.p.Bph),
                     - (Dphi.p.Vr*U.p.Bth + U.p.Vr*Dphi.p.Bth) + (Dphi.p.Vth*U.p.Br + U.p.Vth*Dphi.p.Br),
                     (Dphi.p.Vr*U.p.Bph + U.p.Vr*Dphi.p.Bph) - (Dphi.p.Vph*U.p.Br + U.p.Vph*Dphi.p.Br)};
    Vector DBdth = { (Dtheta.p.Vph*U.p.Bth + U.p.Vph*Dtheta.p.Bth) - (Dtheta.p.Vth*U.p.Bph + U.p.Vth*Dtheta.p.Bph),
                     - (Dtheta.p.Vr*U.p.Bth + U.p.Vr*Dtheta.p.Bth) + (Dtheta.p.Vth*U.p.Br + U.p.Vth*Dtheta.p.Br),
                     (Dtheta.p.Vr*U.p.Bph + U.p.Vr*Dtheta.p.Bph) - (Dtheta.p.Vph*U.p.Br + U.p.Vph*Dtheta.p.Br)};

    Vector B = { DBdph.th - DBdth.ph,
                DBdth.r - DBdr.th,
                DBdr.ph - DBdph.r};

    double p = U.p.Vr*Dr.p.P + U.p.Vph*Dphi.p.P + U.p.Vth*Dtheta.p.P;
    double v = Dr.p.Vr + Dphi.p.Vph + Dtheta.p.Vth;
    double P = p - (-gamma*U.p.P*v);

    Cell res = U;
    res.p.rho = drhodt;
    res.p.Vr = V.r;
    res.p.Vph = V.ph;
    res.p.Vth = V.th;
    res.p.Br  = B.r;
    res.p.Bph = B.ph;
    res.p.Bth = B.th;
    res.p.P = P;
    return res;
}
