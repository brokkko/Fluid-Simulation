#include "Derivative.h"
#include "Constants.h"

Cell F(Cell Dr,Cell Dtheta, Cell Dphi, Cell U)
{
    double r =U.r;
    double phi = U.phi;
    double theta = U.theta;


    double drhodt = 1. / (r * r) * (2 * r * U.p.rho * U.p.Vr + r * r * (Dr.p.rho * U.p.Vr + U.p.rho * Dr.p.Vr))
                    + 1.0 / (r * std::sin(theta)) * (std::cos(theta) * U.p.rho * U.p.Vth + std::sin(theta) * (Dtheta.p.rho * U.p.Vth + U.p.rho * Dtheta.p.Vth))
                    + 1.0 / (r * std::sin(theta)) * (U.p.rho * Dphi.p.Vph + U.p.Vph * Dphi.p.rho);

    drhodt = 1.0/(r*r) * (2*r*U.p.rho*U.p.Vr + r*r*(Dr.p.rho*U.p.Vr + U.p.rho*Dr.p.Vr))
             + 1.0/(r * std::sin(theta)) * (std::cos(theta)*U.p.rho*U.p.Vth + std::sin(theta)*(Dtheta.p.rho*U.p.Vth + U.p.rho*Dtheta.p.Vth))
             + 1.0/(r * std::sin(theta)) * (std::cos(theta)*U.p.rho*U.p.Vph + std::sin(theta)*(Dphi.p.rho*U.p.Vph + U.p.rho*Dphi.p.Vph))
             ;


   /* double dp.Vrdt = (Dr.p.rho*U.p.Vr*U.p.Vr + 2 * U.p.rho * Dr.p.Vr * U.p.Vr
                    +Dphi.p.rho*U.p.Vr*U.p.Vph + U.p.rho * Dphi.p.Vr * U.p.Vph + U.p.rho * U.p.Vr * Dphi.p.Vph
                    +Dtheta.p.rho*U.p.Vr*U.p.Vth + U.p.rho * Dtheta.p.Vr * U.p.Vth + U.p.rho * U.p.Vr * Dtheta.p.Vth
                    +Dr.p.P);*/



    double Mrdt = (-((U.p.rho * U.p.Vth * U.p.Vth + U.p.rho * U.p.Vph * U.p.Vph) / r
                     - Dr.p.P
                     - (U.p.rho * G * Ms) / (r * r)
                     + U.p.Br / mu * Dr.p.Br
                     + 1.0/(r * mu) * (Dtheta.p.Br * U.p.Bth + U.p.Br * Dtheta.p.Bth)
                     + U.p.Bph / mu * Dphi.p.Br
                     - (U.p.Bth * U.p.Bth + U.p.Bph * U.p.Bph) / (mu * r))
                   + 1.0/(r*r) * (2 * r * U.p.rho * U.p.Vr * U.p.Vr + r * r * (Dr.p.rho * U.p.Vr * U.p.Vr + 2 * U.p.rho * U.p.Vr * Dr.p.Vr))
                   + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.p.rho * U.p.Vr * U.p.Vth + std::sin(theta) * (Dtheta.p.rho * U.p.Vr * U.p.Vth + U.p.rho * (Dtheta.p.Vr * U.p.Vth + U.p.Vr * Dtheta.p.Vth)))
    + 1.0/(r*std::sin(theta)) * (Dphi.p.rho * U.p.Vr * U.p.Vph + (Dphi.p.Vr * U.p.Vph + U.p.Vr * Dphi.p.Vph)));

    double p_Mr = U.p.rho*U.p.Vr;
    double p_Mr_dph = Dphi.p.rho*U.p.Vr + U.p.rho*Dphi.p.Vr;
    double p_Mr_dth = Dtheta.p.rho*U.p.Vr + U.p.rho*Dtheta.p.Vr;
    double p_Mr_dr = Dr.p.rho*U.p.Vr + U.p.rho*Dr.p.Vr;

    Mrdt =  1.0/(r*r) * (2*r*(p_Mr*U.p.Vr - U.p.Br*U.p.Br) + r*r*((p_Mr_dr*U.p.Vr + p_Mr*Dr.p.Vr) - (Dr.p.Br*U.p.Br + U.p.Br*Dr.p.Br)))
            + 1.0/(r * std::sin(theta)) * (std::cos(theta)*(p_Mr*U.p.Vth - U.p.Br*U.p.Bth) + std::sin(theta)*((p_Mr_dth*U.p.Vth + p_Mr*Dtheta.p.Vth) - (Dtheta.p.Br*U.p.Bth + U.p.Br*Dtheta.p.Bth)))
            + 1.0/(r * std::sin(theta)) * (std::cos(theta)*(p_Mr*U.p.Vph - U.p.Br*U.p.Bph) + std::sin(theta)*((p_Mr_dph*U.p.Vph + p_Mr*Dphi.p.Vph) - (Dphi.p.Br*U.p.Bph + U.p.Br*Dphi.p.Bph)))
            + (Dr.p.P)
            - (U.p.rho * U.p.Vth*U.p.Vth - U.p.Bth*U.p.Bth)/U.p.rho
            - (U.p.rho * U.p.Vph*U.p.Vph - U.p.Bph*U.p.Bph)/U.p.rho;

    double Mphidt = (-((-U.p.rho * U.p.Vr * U.p.Vph - U.p.rho * U.p.Vth * U.p.Vph * ctg(theta)) / r
                       - Dphi.p.P/r
                       + U.p.Br / mu * Dr.p.Bph
                       + 1.0/(r * r * mu) * (2 * r * U.p.Bph * U.p.Br + r * r * (Dr.p.Bph * U.p.Br + U.p.Bph * Dr.p.Br))
                       + 1.0 / (mu*r) *(Dtheta.p.Bph * U.p.Bth + U.p.Bph * Dtheta.p.Bth)
                       + U.p.Bph / (r * mu) * Dphi.p.Bph
                       + (U.p.Br * U.p.Bph + U.p.Bth * U.p.Bph * ctg(theta)) / (mu * r))
                     + 1.0/(r*r) * (2 * r * U.p.rho * U.p.Vph * U.p.Vr + r * r * (Dr.p.rho * U.p.Vph * U.p.Vr + U.p.rho * (Dr.p.Vph * U.p.Vr + U.p.Vph * Dr.p.Vr)))
                     + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.p.rho * U.p.Vph * U.p.Vth + std::sin(theta) * (Dtheta.p.rho * U.p.Vph * U.p.Vth + U.p.rho * (Dtheta.p.Vph * U.p.Vth + U.p.Vph * Dtheta.p.Vth)))
                     + 1.0/(r*std::sin(theta)) * (Dphi.p.rho * U.p.Vph * U.p.Vph + 2 * U.p.rho * (Dphi.p.Vph * U.p.Vph)));

    double p_Mph = U.p.rho*U.p.Vph;
    double p_Mph_dph = Dphi.p.rho*U.p.Vph + U.p.rho*Dphi.p.Vph;
    double p_Mph_dth = Dtheta.p.rho*U.p.Vph + U.p.rho*Dtheta.p.Vph;
    double p_Mph_dr = Dr.p.rho*U.p.Vph + U.p.rho*Dr.p.Vph;

    Mphidt = 1.0/(r*r*r) * (3*r*r*(p_Mph*U.p.Vr - U.p.Bph*U.p.Br) + r*r*r*((p_Mph_dr*U.p.Vr + p_Mph*Dr.p.Vr) - (Dr.p.Bph*U.p.Br + U.p.Bph*Dr.p.Br)))
             + 1.0/(r * std::sin(theta) * std::sin(theta)) * (2*std::sin(theta)*std::cos(theta)*(p_Mph*U.p.Vth - U.p.Bph*U.p.Bth) + std::sin(theta)*std::sin(theta)*((p_Mph_dth*U.p.Vth + p_Mph*Dtheta.p.Vth) - (Dtheta.p.Bph*U.p.Bth + U.p.Bph*Dtheta.p.Bth)))
             + 1.0/(r * std::sin(theta)) * (std::cos(theta)*(p_Mph*U.p.Vph - U.p.Bph*U.p.Bph) + std::sin(theta)*((p_Mph_dph*U.p.Vph + p_Mph*Dphi.p.Vph) - (Dphi.p.Bph*U.p.Bph + U.p.Bph*Dphi.p.Bph)))
             + 1.0/(r * std::sin(theta)) * (Dphi.p.P);

    double Mthetadt = (-((-U.p.rho * U.p.Vr * U.p.Vth + U.p.rho * U.p.Vph * U.p.Vph * ctg(theta)) / r
                         - Dtheta.p.P/r
                         + 1.0/(r * r * mu) * (2 * r * U.p.Bth * U.p.Br + r * r * (Dr.p.Bth * U.p.Br + U.p.Bth * Dr.p.Br))
                         + 1.0/(r * std::sin(theta) * mu) * (std::cos(theta) * U.p.Bth * U.p.Bth + std::sin(theta) * 2 * U.p.Bth * Dtheta.p.Bth)
                         + 1.0 / (r * std::sin(theta) * mu) * (Dphi.p.Bth * U.p.Bph + U.p.Bth * Dphi.p.Bph)
                         + U.p.Br * U.p.Bth / (r * mu)
                         - (U.p.Bph * U.p.Bph * ctg(theta)) / (mu * r))
                       + 1.0/(r*r) * (2 * r * U.p.rho * U.p.Vth * U.p.Vr + r * r * (Dr.p.rho * U.p.Vth * U.p.Vr + U.p.rho * (Dr.p.Vth * U.p.Vr + U.p.Vth * Dr.p.Vr)))
                       + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.p.rho * U.p.Vth * U.p.Vth + std::sin(theta) * (Dtheta.p.rho * U.p.Vth * U.p.Vth + 2 * U.p.rho * (Dtheta.p.Vth * U.p.Vth)))
                       + 1.0/(r*std::sin(theta)) * (Dphi.p.rho * U.p.Vth * U.p.Vph + U.p.rho * (Dtheta.p.Vph * U.p.Vth + U.p.Vph * Dtheta.p.Vth)));

    double p_Mth = U.p.rho*U.p.Vth;
    double p_Mth_dph = Dphi.p.rho*U.p.Vth + U.p.rho*Dphi.p.Vth;
    double p_Mth_dth = Dtheta.p.rho*U.p.Vth + U.p.rho*Dtheta.p.Vth;
    double p_Mth_dr = Dr.p.rho*U.p.Vth + U.p.rho*Dr.p.Vth;

    Mthetadt = 1.0/(r*r) * (2*r*(p_Mth*U.p.Vr - U.p.Bth*U.p.Br) + r*r*((p_Mth_dr*U.p.Vr + p_Mth*Dr.p.Vr) - (Dr.p.Bth*U.p.Br + U.p.Bth*Dr.p.Br)))
               + 1.0/(r * std::sin(theta)) * (std::cos(theta)*(p_Mth*U.p.Vth - U.p.Bth*U.p.Bth) + std::sin(theta)*((p_Mth_dth*U.p.Vth + p_Mth*Dtheta.p.Vth) - (Dtheta.p.Bth*U.p.Bth + U.p.Bth*Dtheta.p.Bth)))
               + 1.0/(r * std::sin(theta)) * (std::cos(theta)*(p_Mth*U.p.Vph - U.p.Bth*U.p.Bph) + std::sin(theta)*((p_Mth_dph*U.p.Vph + p_Mth*Dphi.p.Vph) - (Dphi.p.Bth*U.p.Bph + U.p.Bth*Dphi.p.Bph)))
               + (Dtheta.p.P)
               - (U.p.rho * U.p.Vth*U.p.Vr - U.p.Bth*U.p.Br)/U.p.rho
               - ctg(theta)*(U.p.rho * U.p.Vph*U.p.Vph - U.p.Bph*U.p.Bph)/U.p.rho;

    double dPdt = gamma*U.p.P*Dr.p.Vr + U.p.Vr*Dr.p.P + gamma*U.p.P*Dphi.p.Vr + U.p.Vr*Dphi.p.P + gamma*U.p.P*Dtheta.p.Vr + U.p.Vr*Dtheta.p.P;
    double dPdt2 = gamma*U.p.P*(Dr.p.Vr+Dphi.p.Vph+Dtheta.p.Vth) + U.p.Vr*Dr.p.P + U.p.Vph*Dphi.p.P + U.p.Vth*Dtheta.p.P;


    double dPdt3 = gamma*U.p.P *((1.0/(r*r) * (2 * r * U.p.Vr + r*r*Dr.p.Vr))
                                 + 1.0/(r*std::sin(theta)) * (Dtheta.p.Vth * std::sin(theta) + U.p.Vth*std::cos(theta))
                                 + 1.0/(r*std::sin(theta)) * Dphi.p.Vph)
                   + U.p.Vr /(r * r) * (2 * r * U.p.P + r*r*(Dr.p.P))
                   + U.p.Vth/(r * std::sin(theta)) * (std::cos(theta) * U.p.P + std::sin(theta) *(Dtheta.p.P))
                   + U.p.Vph / r * (Dphi.p.P);

    double dVphiDt = ( U.p.Vr * Dr.p.Vph + U.p.Vph * Dphi.p.Vph + (1/U.p.rho) * Dphi.p.P);

    Cell res = U;
    res.p.rho = drhodt;
    res.p.Vr = 0;// (Mrdt - drhodt*U.p.Vr )/U.p.rho;
    //res.p.Vr =  (dp.Vrdt - drhodt*U.p.Vr )/U.p.rho;
    res.p.Vph =0;// (Mphidt - drhodt*U.p.Vph)/U.p.rho;
    //res.p.Vph = dVphiDt;
    res.p.Vth = 0; //(Mthetadt - drhodt*U.p.Vth)/U.p.rho;
    res.p.Vth = 0;
    res.p.Br  = 0;
    res.p.Bph = 0;
    res.p.Bth = 0;
    res.p.P = dPdt3;
    return res;
}
