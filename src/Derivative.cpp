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


//    double dp_Vrdt = (Dr.p_rho*U.p_Vr*U.p_Vr + 2 * U.p_rho * Dr.p_Vr * U.p_Vr
//                    +Dphi.p_rho*U.p_Vr*U.p_Vph + U.p_rho * Dphi.p_Vr * U.p_Vph + U.p_rho * U.p_Vr * Dphi.p_Vph
//                    +Dtheta.p_rho*U.p_Vr*U.p_Vth + U.p_rho * Dtheta.p_Vr * U.p_Vth + U.p_rho * U.p_Vr * Dtheta.p_Vth
//                    +Dr.p_P);



    double Mrdt = (-((U.p_rho * U.p_Vth * U.p_Vth + U.p_rho * U.p_Vph * U.p_Vph) / r
                     - Dr.p_P
                     - (U.p_rho * G * Ms) / (r * r)
                     + U.p_Br / mu * Dr.p_Br
                     + 1.0/(r * mu) * (Dtheta.p_Br * U.p_Bth + U.p_Br * Dtheta.p_Bth)
                     + U.p_Bph / mu * Dphi.p_Br
                     - (U.p_Bth * U.p_Bth + U.p_Bph * U.p_Bph) / (mu * r))
                   + 1.0/(r*r) * (2 * r * U.p_rho * U.p_Vr * U.p_Vr + r * r * (Dr.p_rho * U.p_Vr * U.p_Vr + 2 * U.p_rho * U.p_Vr * Dr.p_Vr))
                   + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.p_rho * U.p_Vr * U.p_Vth + std::sin(theta) * (Dtheta.p_rho * U.p_Vr * U.p_Vth + U.p_rho * (Dtheta.p_Vr * U.p_Vth + U.p_Vr * Dtheta.p_Vth)))
    );//+ 1.0/(r*std::sin(theta)) * (Dphi.p_rho * U.p_Vr * U.p_Vph + (Dphi.p_Vr * U.p_Vph + U.p_Vr * Dphi.p_Vph)));

    double Mphidt = (-((-U.p_rho * U.p_Vr * U.p_Vph - U.p_rho * U.p_Vth * U.p_Vph * ctg(theta))*0 / r
                       - Dphi.p_P/r
                       + U.p_Br / mu * Dr.p_Bph
                       + 1.0/(r * r * mu) * (2 * r * U.p_Bph * U.p_Br + r * r * (Dr.p_Bph * U.p_Br + U.p_Bph * Dr.p_Br))
                       + 1.0 / (mu*r) *(Dtheta.p_Bph * U.p_Bth + U.p_Bph * Dtheta.p_Bth)
                       + U.p_Bph / (r * mu) * Dphi.p_Bph
                       + (U.p_Br * U.p_Bph + U.p_Bth * U.p_Bph * ctg(theta)) / (mu * r))
                     + 1.0/(r*r) * (2 * r * U.p_rho * U.p_Vph * U.p_Vr + r * r * (Dr.p_rho * U.p_Vph * U.p_Vr + U.p_rho * (Dr.p_Vph * U.p_Vr + U.p_Vph * Dr.p_Vr)))
                     + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.p_rho * U.p_Vph * U.p_Vth + std::sin(theta) * (Dtheta.p_rho * U.p_Vph * U.p_Vth + U.p_rho * (Dtheta.p_Vph * U.p_Vth + U.p_Vph * Dtheta.p_Vth)))
                     + 1.0/(r*std::sin(theta)) * (Dphi.p_rho * U.p_Vph * U.p_Vph + 2 * U.p_rho * (Dphi.p_Vph * U.p_Vph)));

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

    double dPdt = gamma*U.p_P*Dr.p_Vr + U.p_Vr*Dr.p_P + gamma*U.p_P*Dphi.p_Vr + U.p_Vr*Dphi.p_P + gamma*U.p_P*Dtheta.p_Vr + U.p_Vr*Dtheta.p_P;
    double dPdt2 = gamma*U.p_P*(Dr.p_Vr+Dphi.p_Vph+Dtheta.p_Vth) + U.p_Vr*Dr.p_P + U.p_Vph*Dphi.p_P + U.p_Vth*Dtheta.p_P;


    double dPdt3 = gamma*U.p_P *((1.0/(r*r) * (2 * r * U.p_Vr + r*r*Dr.p_Vr))
                                 + 1.0/(r*std::sin(theta)) * (Dtheta.p_Vth * std::sin(theta) + U.p_Vth*std::cos(theta))
                                 + 1.0/(r*std::sin(theta)) * Dphi.p_Vph)
                   + U.p_Vr /(r * r) * (2 * r * U.p_P + r*r*(Dr.p_P))
                   + U.p_Vth/(r * std::sin(theta)) * (std::cos(theta) * U.p_P + std::sin(theta) *(Dtheta.p_P))
                   + U.p_Vph / r * (Dphi.p_P);

    double dVphiDt = ( U.p_Vr * Dr.p_Vph + U.p_Vph * Dphi.p_Vph + (1/U.p_rho) * Dphi.p_P)/r;

    Cell res = U;
    res.p_rho = drhodt;
    res.p_Vr =  (Mrdt - drhodt*U.p_Vr )/U.p_rho;
    //res.p_Vr =  (dp_Vrdt - drhodt*U.p_Vr )/U.p_rho;
    //res.p_Vph = (Mphidt - drhodt*U.p_Vph)/U.p_rho;
    res.p_Vph = dVphiDt;
    res.p_Vth = (Mthetadt - drhodt*U.p_Vth)/U.p_rho;
    res.p_Vth = 0;
    res.p_Br  = 0;
    res.p_Bph = 0;
    res.p_Bth = 0;
    res.p_P = dPdt3;
    return res;
}
