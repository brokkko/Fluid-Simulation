#include "Derivative.h"
#include "Constants.h"

Cell F(Cell Dr,Cell Dtheta, Cell Dphi, Cell U)
{
    double r =U.r;
    double phi = U.phi;
    double theta = U.theta;


    double drhodt = 1. / (r * r) * (2 * r * U.p_rho * U.p_Vr + r * r * (Dr.p_rho * U.p_Vr + U.p_rho * Dr.p_Vr))
                    + 1.0 / (r * std::sin(theta)) * (std::cos(theta) * U.p_rho * U.p_Vth + std::sin(theta) * (Dtheta.p_rho * U.p_Vth + U.p_rho * Dtheta.p_Vth))
                    + 1.0 / (r * std::sin(theta)) * (U.p_rho * Dphi.p_Vph + U.p_Vph * Dphi.p_rho);

    drhodt = 1.0/(r*r) * (2*r*U.p_rho*U.p_Vr + r*r*(Dr.p_rho*U.p_Vr + U.p_rho*Dr.p_Vr))
             + 1.0/(r * std::sin(theta)) * (std::cos(theta)*U.p_rho*U.p_Vth + std::sin(theta)*(Dtheta.p_rho*U.p_Vth + U.p_rho*Dtheta.p_Vth))
             + 1.0/(r * std::sin(theta)) * (std::cos(theta)*U.p_rho*U.p_Vph + std::sin(theta)*(Dphi.p_rho*U.p_Vph + U.p_rho*Dphi.p_Vph))
             ;


   /* double dp_Vrdt = (Dr.p_rho*U.p_Vr*U.p_Vr + 2 * U.p_rho * Dr.p_Vr * U.p_Vr
                    +Dphi.p_rho*U.p_Vr*U.p_Vph + U.p_rho * Dphi.p_Vr * U.p_Vph + U.p_rho * U.p_Vr * Dphi.p_Vph
                    +Dtheta.p_rho*U.p_Vr*U.p_Vth + U.p_rho * Dtheta.p_Vr * U.p_Vth + U.p_rho * U.p_Vr * Dtheta.p_Vth
                    +Dr.p_P);*/



    double Mrdt = (-((U.p_rho * U.p_Vth * U.p_Vth + U.p_rho * U.p_Vph * U.p_Vph) / r
                     - Dr.p_P
                     - (U.p_rho * G * Ms) / (r * r)
                     + U.p_Br / mu * Dr.p_Br
                     + 1.0/(r * mu) * (Dtheta.p_Br * U.p_Bth + U.p_Br * Dtheta.p_Bth)
                     + U.p_Bph / mu * Dphi.p_Br
                     - (U.p_Bth * U.p_Bth + U.p_Bph * U.p_Bph) / (mu * r))
                   + 1.0/(r*r) * (2 * r * U.p_rho * U.p_Vr * U.p_Vr + r * r * (Dr.p_rho * U.p_Vr * U.p_Vr + 2 * U.p_rho * U.p_Vr * Dr.p_Vr))
                   + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.p_rho * U.p_Vr * U.p_Vth + std::sin(theta) * (Dtheta.p_rho * U.p_Vr * U.p_Vth + U.p_rho * (Dtheta.p_Vr * U.p_Vth + U.p_Vr * Dtheta.p_Vth)))
    + 1.0/(r*std::sin(theta)) * (Dphi.p_rho * U.p_Vr * U.p_Vph + (Dphi.p_Vr * U.p_Vph + U.p_Vr * Dphi.p_Vph)));

    double p_Mr = U.p_rho*U.p_Vr;
    double p_Mr_dph = Dphi.p_rho*U.p_Vr + U.p_rho*Dphi.p_Vr;
    double p_Mr_dth = Dtheta.p_rho*U.p_Vr + U.p_rho*Dtheta.p_Vr;
    double p_Mr_dr = Dr.p_rho*U.p_Vr + U.p_rho*Dr.p_Vr;

    Mrdt =  1.0/(r*r) * (2*r*(p_Mr*U.p_Vr - U.p_Br*U.p_Br) + r*r*((p_Mr_dr*U.p_Vr + p_Mr*Dr.p_Vr) - (Dr.p_Br*U.p_Br + U.p_Br*Dr.p_Br)))
            + 1.0/(r * std::sin(theta)) * (std::cos(theta)*(p_Mr*U.p_Vth - U.p_Br*U.p_Bth) + std::sin(theta)*((p_Mr_dth*U.p_Vth + p_Mr*Dtheta.p_Vth) - (Dtheta.p_Br*U.p_Bth + U.p_Br*Dtheta.p_Bth)))
            + 1.0/(r * std::sin(theta)) * (std::cos(theta)*(p_Mr*U.p_Vph - U.p_Br*U.p_Bph) + std::sin(theta)*((p_Mr_dph*U.p_Vph + p_Mr*Dphi.p_Vph) - (Dphi.p_Br*U.p_Bph + U.p_Br*Dphi.p_Bph)))
            + (Dr.p_P)
            - (U.p_rho * U.p_Vth*U.p_Vth - U.p_Bth*U.p_Bth)/U.p_rho
            - (U.p_rho * U.p_Vph*U.p_Vph - U.p_Bph*U.p_Bph)/U.p_rho;

    double Mphidt = (-((-U.p_rho * U.p_Vr * U.p_Vph - U.p_rho * U.p_Vth * U.p_Vph * ctg(theta)) / r
                       - Dphi.p_P/r
                       + U.p_Br / mu * Dr.p_Bph
                       + 1.0/(r * r * mu) * (2 * r * U.p_Bph * U.p_Br + r * r * (Dr.p_Bph * U.p_Br + U.p_Bph * Dr.p_Br))
                       + 1.0 / (mu*r) *(Dtheta.p_Bph * U.p_Bth + U.p_Bph * Dtheta.p_Bth)
                       + U.p_Bph / (r * mu) * Dphi.p_Bph
                       + (U.p_Br * U.p_Bph + U.p_Bth * U.p_Bph * ctg(theta)) / (mu * r))
                     + 1.0/(r*r) * (2 * r * U.p_rho * U.p_Vph * U.p_Vr + r * r * (Dr.p_rho * U.p_Vph * U.p_Vr + U.p_rho * (Dr.p_Vph * U.p_Vr + U.p_Vph * Dr.p_Vr)))
                     + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.p_rho * U.p_Vph * U.p_Vth + std::sin(theta) * (Dtheta.p_rho * U.p_Vph * U.p_Vth + U.p_rho * (Dtheta.p_Vph * U.p_Vth + U.p_Vph * Dtheta.p_Vth)))
                     + 1.0/(r*std::sin(theta)) * (Dphi.p_rho * U.p_Vph * U.p_Vph + 2 * U.p_rho * (Dphi.p_Vph * U.p_Vph)));

    double p_Mph = U.p_rho*U.p_Vph;
    double p_Mph_dph = Dphi.p_rho*U.p_Vph + U.p_rho*Dphi.p_Vph;
    double p_Mph_dth = Dtheta.p_rho*U.p_Vph + U.p_rho*Dtheta.p_Vph;
    double p_Mph_dr = Dr.p_rho*U.p_Vph + U.p_rho*Dr.p_Vph;

    Mphidt = 1.0/(r*r*r) * (3*r*r*(p_Mph*U.p_Vr - U.p_Bph*U.p_Br) + r*r*r*((p_Mph_dr*U.p_Vr + p_Mph*Dr.p_Vr) - (Dr.p_Bph*U.p_Br + U.p_Bph*Dr.p_Br)))
             + 1.0/(r * std::sin(theta) * std::sin(theta)) * (2*std::sin(theta)*std::cos(theta)*(p_Mph*U.p_Vth - U.p_Bph*U.p_Bth) + std::sin(theta)*std::sin(theta)*((p_Mph_dth*U.p_Vth + p_Mph*Dtheta.p_Vth) - (Dtheta.p_Bph*U.p_Bth + U.p_Bph*Dtheta.p_Bth)))
             + 1.0/(r * std::sin(theta)) * (std::cos(theta)*(p_Mph*U.p_Vph - U.p_Bph*U.p_Bph) + std::sin(theta)*((p_Mph_dph*U.p_Vph + p_Mph*Dphi.p_Vph) - (Dphi.p_Bph*U.p_Bph + U.p_Bph*Dphi.p_Bph)))
             + 1.0/(r * std::sin(theta)) * (Dphi.p_P);

    double Mthetadt = (-((-U.p_rho * U.p_Vr * U.p_Vth + U.p_rho * U.p_Vph * U.p_Vph * ctg(theta)) / r
                         - Dtheta.p_P/r
                         + 1.0/(r * r * mu) * (2 * r * U.p_Bth * U.p_Br + r * r * (Dr.p_Bth * U.p_Br + U.p_Bth * Dr.p_Br))
                         + 1.0/(r * std::sin(theta) * mu) * (std::cos(theta) * U.p_Bth * U.p_Bth + std::sin(theta) * 2 * U.p_Bth * Dtheta.p_Bth)
                         + 1.0 / (r * std::sin(theta) * mu) * (Dphi.p_Bth * U.p_Bph + U.p_Bth * Dphi.p_Bph)
                         + U.p_Br * U.p_Bth / (r * mu)
                         - (U.p_Bph * U.p_Bph * ctg(theta)) / (mu * r))
                       + 1.0/(r*r) * (2 * r * U.p_rho * U.p_Vth * U.p_Vr + r * r * (Dr.p_rho * U.p_Vth * U.p_Vr + U.p_rho * (Dr.p_Vth * U.p_Vr + U.p_Vth * Dr.p_Vr)))
                       + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.p_rho * U.p_Vth * U.p_Vth + std::sin(theta) * (Dtheta.p_rho * U.p_Vth * U.p_Vth + 2 * U.p_rho * (Dtheta.p_Vth * U.p_Vth)))
                       + 1.0/(r*std::sin(theta)) * (Dphi.p_rho * U.p_Vth * U.p_Vph + U.p_rho * (Dtheta.p_Vph * U.p_Vth + U.p_Vph * Dtheta.p_Vth)));

    double p_Mth = U.p_rho*U.p_Vth;
    double p_Mth_dph = Dphi.p_rho*U.p_Vth + U.p_rho*Dphi.p_Vth;
    double p_Mth_dth = Dtheta.p_rho*U.p_Vth + U.p_rho*Dtheta.p_Vth;
    double p_Mth_dr = Dr.p_rho*U.p_Vth + U.p_rho*Dr.p_Vth;

    Mthetadt = 1.0/(r*r) * (2*r*(p_Mth*U.p_Vr - U.p_Bth*U.p_Br) + r*r*((p_Mth_dr*U.p_Vr + p_Mth*Dr.p_Vr) - (Dr.p_Bth*U.p_Br + U.p_Bth*Dr.p_Br)))
               + 1.0/(r * std::sin(theta)) * (std::cos(theta)*(p_Mth*U.p_Vth - U.p_Bth*U.p_Bth) + std::sin(theta)*((p_Mth_dth*U.p_Vth + p_Mth*Dtheta.p_Vth) - (Dtheta.p_Bth*U.p_Bth + U.p_Bth*Dtheta.p_Bth)))
               + 1.0/(r * std::sin(theta)) * (std::cos(theta)*(p_Mth*U.p_Vph - U.p_Bth*U.p_Bph) + std::sin(theta)*((p_Mth_dph*U.p_Vph + p_Mth*Dphi.p_Vph) - (Dphi.p_Bth*U.p_Bph + U.p_Bth*Dphi.p_Bph)))
               + (Dtheta.p_P)
               - (U.p_rho * U.p_Vth*U.p_Vr - U.p_Bth*U.p_Br)/U.p_rho
               - ctg(theta)*(U.p_rho * U.p_Vph*U.p_Vph - U.p_Bph*U.p_Bph)/U.p_rho;

    double dPdt = gamma*U.p_P*Dr.p_Vr + U.p_Vr*Dr.p_P + gamma*U.p_P*Dphi.p_Vr + U.p_Vr*Dphi.p_P + gamma*U.p_P*Dtheta.p_Vr + U.p_Vr*Dtheta.p_P;
    double dPdt2 = gamma*U.p_P*(Dr.p_Vr+Dphi.p_Vph+Dtheta.p_Vth) + U.p_Vr*Dr.p_P + U.p_Vph*Dphi.p_P + U.p_Vth*Dtheta.p_P;


    double dPdt3 = gamma*U.p_P *((1.0/(r*r) * (2 * r * U.p_Vr + r*r*Dr.p_Vr))
                                 + 1.0/(r*std::sin(theta)) * (Dtheta.p_Vth * std::sin(theta) + U.p_Vth*std::cos(theta))
                                 + 1.0/(r*std::sin(theta)) * Dphi.p_Vph)
                   + U.p_Vr /(r * r) * (2 * r * U.p_P + r*r*(Dr.p_P))
                   + U.p_Vth/(r * std::sin(theta)) * (std::cos(theta) * U.p_P + std::sin(theta) *(Dtheta.p_P))
                   + U.p_Vph / r * (Dphi.p_P);

    double dVphiDt = ( U.p_Vr * Dr.p_Vph + U.p_Vph * Dphi.p_Vph + (1/U.p_rho) * Dphi.p_P);

    Cell res = U;
    res.p_rho = 0;//drhodt;
    res.p_Vr =  0 ; //(Mrdt - drhodt*U.p_Vr )/U.p_rho;
    //res.p_Vr =  (dp_Vrdt - drhodt*U.p_Vr )/U.p_rho;
    res.p_Vph = 0; // (Mphidt - drhodt*U.p_Vph)/U.p_rho;
    //res.p_Vph = dVphiDt;
    res.p_Vth = 0; //(Mthetadt - drhodt*U.p_Vth)/U.p_rho;
    res.p_Vth = 0;
    res.p_Br  = 0;
    res.p_Bph = 0;
    res.p_Bth = 0;
    res.p_P = 0; //dPdt3;
    return res;
}
