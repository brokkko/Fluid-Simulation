#include <algorithm>
#include <iostream>
#include <tuple>
#include "Simulation.h"
#include "Constants.h"

Cell SlopeLim(Cell r)
{
    /* return {std::max(0.0,std::max(std::min(2*r.p_rho,1.0),std::min(r.p_rho,2.0))),
             std::max(0.0,std::max(std::min(2*r.p_Vr,1.0),std::min(r.p_Vr,2.0))),
             std::max(0.0,std::max(std::min(2*r.p_Vph,1.0),std::min(r.p_Vph,2.0))),
                 std::max(0.0,std::max(std::min(2*r.B,1.0),std::min(r.B,2.0)))};*/
    Cell res = r;
    auto res_p = reinterpret_cast<double*>(&res);
    auto r_p = reinterpret_cast<double*>(&r);
    for(int i=0;i<16;i++)
        res_p[i] = std::max(0.0, std::min(1.0, r_p[i]));
    return res;


   /* return {std::max(0.0, std::min(1.0, r.p_rho)),
            std::max(0.0, std::min(1.0, r.p_Vr)),
            std::max(0.0, std::min(1.0, r.p_Vph)),
            std::max(0.0, std::min(1.0, r.p_Vth)),
            std::max(0.0, std::min(1.0, r.p_Br)),
            std::max(0.0, std::min(1.0, r.p_Bph)),
            std::max(0.0, std::min(1.0, r.p_Bth)),
            std::max(0.0, std::min(1.0, r.c_E))};*/

    /* return Cell{std::max(0.0, 1.5 * (r.p_rho * r.p_rho + r.p_rho) / (r.p_rho * r.p_rho + r.p_rho + 1)),
                 std::max(0.0, 1.5 * (r.p_Vr * r.p_Vr + r.p_Vr) / (r.p_Vr * r.p_Vr + r.p_Vr + 1)),
                          std::max(0.0, 1.5 * (r.p_Vph * r.p_Vph + r.p_Vph) / (r.p_Vph * r.p_Vph + r.p_Vph + 1)),
                          std::max(0.0, 1.5 * (r.B * r.B + r.B) / (r.B * r.B + r.B + 1))};*/



}

double nonZeroDouble(double val)
{
    if (std::abs(val) < __DBL_EPSILON__)
    {
        if (val < 0)
            return -__DBL_EPSILON__;
        else
            return __DBL_EPSILON__;
    }
    return val;
}

Cell nonZeroDenom(Cell denom)
{
    Cell res = denom;
    auto res_p = reinterpret_cast<double*>(&res);
    auto d_p = reinterpret_cast<double*>(&denom);
    for(int i=0;i<16;i++)
        res_p[i] = nonZeroDouble(d_p[i]);
    return res;
}

Cell F(Cell Dr,Cell Dtheta, Cell Dphi, Cell U)
{
    double r =U.r;
    double phi = U.phi;
    double theta = U.theta;
    //theta=M_PI_2;

    // p_rho * vr = p_rho/dr * vr + p_rho * vr/dr
    //vr/dr = (p_rho * vr - p_rho/dr * vr)/p_rho
    double p = gamma * (U.c_E - 0.5 * U.p_rho * (U.p_Vr * U.p_Vr + U.p_Vth * U.p_Vth + U.p_Vph * U.p_Vph)
                        - 0.5/mu * (U.p_Br * U.p_Br + U.p_Bph * U.p_Bph + U.p_Bth * U.p_Bth));
    double dpdr = gamma * (Dr.c_E - 0.5 * Dr.p_rho * (U.p_Vr * U.p_Vr + U.p_Vph * U.p_Vph + U.p_Vth * U.p_Vth)
                           - U.p_rho * (Dr.p_Vr * U.p_Vr + Dr.p_Vph * U.p_Vph + Dr.p_Vth * U.p_Vth)
                           - (Dr.p_Br * U.p_Br + Dr.p_Bph * U.p_Bph + Dr.p_Bth * U.p_Bth) / mu);
    double dpdphi = gamma * (Dphi.c_E - 1.0 / 2 * Dphi.p_rho * (U.p_Vr * U.p_Vr + U.p_Vph * U.p_Vph + U.p_Vth * U.p_Vth)
                             - U.p_rho * (Dphi.p_Vr * U.p_Vr + Dphi.p_Vph * U.p_Vph + Dphi.p_Vth * U.p_Vth)
                             - (Dphi.p_Br * U.p_Br + Dphi.p_Bph * U.p_Bph + Dphi.p_Bth * U.p_Bth) / mu);
    double dpdtheta =gamma * (Dtheta.c_E - 0.5 * Dtheta.p_rho * (U.p_Vr * U.p_Vr + U.p_Vph * U.p_Vph + U.p_Vth * U.p_Vth)
                              - U.p_rho * (Dtheta.p_Vr * U.p_Vr + Dtheta.p_Vph * U.p_Vph + Dtheta.p_Vth * U.p_Vth)
                              - (Dtheta.p_Br * U.p_Br + Dtheta.p_Bph * U.p_Bph + Dtheta.p_Bth * U.p_Bth) / mu);

    double P = p + 0.5 / mu * (U.p_Br * U.p_Br + U.p_Bph * U.p_Bph + U.p_Bth * U.p_Bth);
    double dPdr = dpdr + (Dr.p_Br * U.p_Br + Dr.p_Bph * U.p_Bph + Dr.p_Bth * U.p_Bth) / mu;
    double dPdphi = dpdphi + (Dphi.p_Br * U.p_Br + Dphi.p_Bph * U.p_Bph + Dphi.p_Bth * U.p_Bth) / mu;
    double dPdtheta = dpdtheta + (Dtheta.p_Br * U.p_Br + Dtheta.p_Bph * U.p_Bph + Dtheta.p_Bth * U.p_Bth) / mu;

    /*double drhodt=   Dr.p_rho * U.p_Vr + U.p_rho * Dr.p_Vr +
                     Dphi.p_rho * U.p_Vph + U.p_rho * Dphi.p_Vph;*/

   //double drhodt = U.p_Vr*Dr.p_rho + U.p_rho*Dr.p_Vr + 2.0/r*U.p_rho*U.p_Vr;
    double drhodt = 1. /(r * r) * (2 * r * U.p_rho * U.p_Vr + r * r * (Dr.p_rho * U.p_Vr + U.p_rho * Dr.p_Vr))
                 + 1.0 /(r * std::sin(theta)) * (std::cos(theta) * U.p_rho * U.p_Vth + std::sin(theta) * (Dtheta.p_rho * U.p_Vth + U.p_rho * Dtheta.p_Vth))
                 + 1.0 / r * (U.p_rho * Dphi.p_Vph + U.p_Vph * Dphi.p_rho);
    /*double drhodt = 1./(r*r) * (2 * r * U.p_Vr + r*r * Dr.p_Vr)
                    + 1.0/(r * std::sin(theta))*(std::cos(theta)*U.p_Vth + std::sin(theta)*Dtheta.p_Vth)
                    + 1.0/r * Dphi.p_Vph;*/

    double Mrdt = (-((U.p_rho * U.p_Vth * U.p_Vth + U.p_rho * U.p_Vph * U.p_Vph) / r
                     - Dr.p_P
                     - (U.p_rho * G * Ms) / (r * r)
                     + U.p_Br / mu * Dr.p_Br
                     + 1.0/(r * mu) * (Dtheta.p_Br * U.p_Bth + U.p_Br * Dtheta.p_Bth)
                     + U.p_Bph / mu * Dphi.p_Br
                     - (U.p_Bth * U.p_Bth + U.p_Bph * U.p_Bph) / (mu * r))
                   + 1.0/(r*r) * (2 * r * U.p_rho * U.p_Vr * U.p_Vr + r * r * (Dr.p_rho * U.p_Vr * U.p_Vr + 2 * U.p_rho * U.p_Vr * Dr.p_Vr)
                                  + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.p_rho * U.p_Vr * U.p_Vth + std::sin(theta) * (Dtheta.p_rho * U.p_Vr * U.p_Vth + U.p_rho * (Dtheta.p_Vr * U.p_Vth + U.p_Vr * Dtheta.p_Vth)))
                                  + 1.0/(r*std::sin(theta)) * (Dphi.p_rho * U.p_Vr * U.p_Vph + (Dphi.p_Vr * U.p_Vph + U.p_Vr * Dphi.p_Vph))));

    double Mphidt = (-((-U.p_rho * U.p_Vr * U.p_Vph - U.p_rho * U.p_Vth * U.p_Vph * ctg(theta)) / r
                       - Dphi.p_P/r
                       + U.p_Br / mu * Dr.p_Bph
                       + 1.0/(r * r * mu) * (2 * r * U.p_Bph * U.p_Br + r * r * (Dr.p_Bph * U.p_Br + U.p_Bph * Dr.p_Br))
                       + 1.0 / (mu*r) *(Dtheta.p_Bph * U.p_Bth + U.p_Bph * Dtheta.p_Bth)
                       + U.p_Bph / (r * mu) * Dphi.p_Bph
                       + (U.p_Br * U.p_Bph + U.p_Bth * U.p_Bph * ctg(theta)) / (mu * r))
                     + 1.0/(r*r) * (2 * r * U.p_rho * U.p_Vph * U.p_Vr + r * r * (Dr.p_rho * U.p_Vph * U.p_Vr + U.p_rho * (Dr.p_Vph * U.p_Vr + U.p_Vph * Dr.p_Vr))
                                    + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.p_rho * U.p_Vph * U.p_Vth + std::sin(theta) * (Dtheta.p_rho * U.p_Vph * U.p_Vth + U.p_rho * (Dtheta.p_Vph * U.p_Vth + U.p_Vph * Dtheta.p_Vth)))
                                    + 1.0/(r*std::sin(theta)) * (Dphi.p_rho * U.p_Vph * U.p_Vph + 2 * U.p_rho * (Dphi.p_Vph * U.p_Vph))));

    double Mthetadt = (-((-U.p_rho * U.p_Vr * U.p_Vth + U.p_rho * U.p_Vph * U.p_Vph * ctg(theta)) / r
                         - Dtheta.p_P/r
                         + 1.0/(r * r * mu) * (2 * r * U.p_Bth * U.p_Br + r * r * (Dr.p_Bth * U.p_Br + U.p_Bth * Dr.p_Br))
                         + 1.0/(r * std::sin(theta) * mu) * (std::cos(theta) * U.p_Bth * U.p_Bth + std::sin(theta) * 2 * U.p_Bth * Dtheta.p_Bth)
                         + 1.0 / (r * std::sin(theta) * mu) * (Dphi.p_Bth * U.p_Bph + U.p_Bth * Dphi.p_Bph)
                         + U.p_Br * U.p_Bth / (r * mu)
                         - (U.p_Bph * U.p_Bph * ctg(theta)) / (mu * r))
                       + 1.0/(r*r) * (2 * r * U.p_rho * U.p_Vth * U.p_Vr + r * r * (Dr.p_rho * U.p_Vth * U.p_Vr + U.p_rho * (Dr.p_Vth * U.p_Vr + U.p_Vth * Dr.p_Vr))
                                      + 1.0/(r*std::sin(theta))*(std::cos(theta) * U.p_rho * U.p_Vth * U.p_Vth + std::sin(theta) * (Dtheta.p_rho * U.p_Vth * U.p_Vth + 2 * U.p_rho * (Dtheta.p_Vth * U.p_Vth)))
                                      + 1.0/(r*std::sin(theta)) * (Dphi.p_rho * U.p_Vth * U.p_Vph + U.p_rho * (Dtheta.p_Vph * U.p_Vth + U.p_Vph * Dtheta.p_Vth))));


    double Edt = 1.0/(r*r) * (2 * r * U.c_E/U.volume * U.p_Vr + r * r * (Dr.c_E * U.p_Vr + U.c_E * Dr.p_Vr))
                 + 1.0/(r*std::sin(theta)) * (std::cos(theta) * U.c_E * U.p_Vth + std::sin(theta) * (Dtheta.c_E * U.p_Vth + U.c_E * Dtheta.p_Vth))
                 + 1.0/(r*std::sin(theta)) * (Dphi.c_E * U.p_Vph + U.c_E * Dphi.p_Vph)
                 -(- 1.0/(r*r) * (2*r*P*U.p_Vr + r * r * (dPdr * U.p_Vr + P * Dr.p_Vr))
                   - 1.0/(r*std::sin(theta)) * (std::cos(theta)*P*U.p_Vth + std::sin(theta) * (dPdtheta * U.p_Vth + P * Dtheta.p_Vth))
                   - 1.0/(r*std::sin(theta)) *(dPdphi*U.p_Vph + P * Dphi.p_Vph)
                   - U.p_rho * U.p_Vr * G * Ms / (r * r)
                   + 1.0/(r*r*mu) * (2 * r * U.p_Br * (U.p_Br * U.p_Vr + U.p_Bth * U.p_Vth + U.p_Bph * U.p_Vph) + r * r * (Dr.p_Br * (U.p_Br * U.p_Vr + U.p_Bth * U.p_Vth + U.p_Bph * U.p_Vph)
                                                                                                                           + U.p_Br * ((Dr.p_Br * U.p_Vr + U.p_Br * Dr.p_Vr) + (Dr.p_Bth * U.p_Vth + U.p_Bth * Dr.p_Vth) + (Dr.p_Bph * U.p_Vph + U.p_Bph * Dr.p_Vph))))
                   + 1.0/(r*std::sin(theta)*mu) *(std::cos(theta) * U.p_Bth * (U.p_Br * U.p_Vr + U.p_Bth * U.p_Vth + U.p_Bph * U.p_Vph)
                        + std::sin(theta)*(Dtheta.p_Bth * (U.p_Br * U.p_Vr + U.p_Bth * U.p_Vth + U.p_Bph * U.p_Vph)
                        + U.p_Bth * ((Dtheta.p_Br * U.p_Vr + U.p_Br * Dtheta.p_Vr) + (Dtheta.p_Bth * U.p_Vth + U.p_Bth * Dtheta.p_Vth) + (Dtheta.p_Bph * U.p_Vph + U.p_Bph * Dtheta.p_Vph))))
                   + 1.0/(r*std::sin(theta)*mu) * (Dphi.p_Bph * (U.p_Br * U.p_Vr + U.p_Bth * U.p_Vth + U.p_Bph * U.p_Vph)
                        + U.p_Bph * ((Dphi.p_Br * U.p_Vr + U.p_Br * Dphi.p_Vr) + (Dphi.p_Bth * U.p_Vth + U.p_Bth * Dphi.p_Vth) + (Dphi.p_Bph * U.p_Vph + U.p_Bph * Dphi.p_Vph))))
                         ;

    if(std::abs(Edt)>10000)
        int a=0;

    Cell res = U;
    res.p_rho = drhodt;
    res.p_Vr =  (Mrdt - drhodt*U.p_Vr )/U.p_rho;
    res.p_Vph = (Mrdt - drhodt*U.p_Vph)/U.p_rho;
    res.p_Vth = (Mrdt - drhodt*U.p_Vth)/U.p_rho;
    res.p_Br = 0;
    res.p_Bph = 0;
    res.p_Bth = 0;
    res.p_P = 0;
    return res;
}

Cell FluxR(Cell U)
{
    double vB = U.p_Vr * U.p_Br + U.p_Vph * U.p_Bph + U.p_Vth * U.p_Bth;
    Cell res  = U;
    res.c_rho = U.c_Mr;
    res.c_Mr  = U.c_Mr  * U.p_Vr - U.p_Br * U.p_Br + U.p_P;
    res.c_Mph = U.c_Mph * U.p_Vr - U.p_Br * U.p_Bph;
    res.c_Mth = U.c_Mth * U.p_Vr - U.p_Br * U.p_Bth;

    res.c_Br = 0;
    res.c_Bph = 0;
    res.c_Bth = 0;
    res.c_E = (U.c_E + U.p_P)*U.p_Vr - U.p_Br*vB;
    return res;
}

Cell FluxPhi(Cell U)
{
    double vB = U.p_Vr * U.p_Br + U.p_Vph * U.p_Bph + U.p_Vth * U.p_Bth;
    Cell res  = U;
    res.c_rho = U.c_Mph;
    res.c_Mr  = U.c_Mr *  U.p_Vph - U.p_Bph * U.p_Br + U.p_P;
    res.c_Mph = U.c_Mph * U.p_Vph - U.p_Bph * U.p_Bph;
    res.c_Mth = U.c_Mth * U.p_Vph - U.p_Bph * U.p_Bth;

    res.c_Br = 0;
    res.c_Bph = 0;
    res.c_Bth = 0;
    res.c_E = (U.c_E + U.p_P)*U.p_Vph - U.p_Bph*vB;
    return res;
}
Cell FluxTheta(Cell U)
{
    double vB = U.p_Vr * U.p_Br + U.p_Vph * U.p_Bph + U.p_Vth * U.p_Bth;
    Cell res  = U;
    res.c_rho = U.c_Mth;
    res.c_Mr  = U.c_Mr  * U.p_Vth - U.p_Bth * U.p_Br + U.p_P;
    res.c_Mph = U.c_Mph * U.p_Vth - U.p_Bth * U.p_Bph;
    res.c_Mth = U.c_Mth * U.p_Vth - U.p_Bth * U.p_Bth;

    res.c_Br = 0;
    res.c_Bph = 0;
    res.c_Bth = 0;
    res.c_E = (U.c_E + U.p_P)*U.p_Vth - U.p_Bth*vB;
    return res;
}

double maxSpeed(Cell U){
    double p = gamma * (U.c_E - 0.5 * U.p_rho * (U.p_Vr * U.p_Vr + U.p_Vth * U.p_Vth + U.p_Vph * U.p_Vph)
                        - 0.5/mu * (U.p_Br * U.p_Br + U.p_Bph * U.p_Bph + U.p_Bth * U.p_Bth));
    double P = p + 0.5 / mu * (U.p_Br * U.p_Br + U.p_Bph * U.p_Bph + U.p_Bth * U.p_Bth);

    double B_2 = U.p_Br * U.p_Br + U.p_Bph * U.p_Bph + U.p_Bth * U.p_Bth;
    double B = std::sqrt(B_2);

    double cmax =
            std::sqrt(U.p_Vr * U.p_Vr + U.p_Vph * U.p_Vph + U.p_Vth * U.p_Vth)
            + 0.5*((gamma * p + B_2)/U.p_rho
            + std::sqrt(((gamma + B)/U.p_rho) * ((gamma + B) / U.p_rho) - 4 * (gamma * U.p_Br * U.p_Br) / (U.p_rho * U.p_rho)));

    return cmax*DT/CELL_SIZE;
}


void CalculateFlux(std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> out, SphericalGrid& in)
{
    for(int theta = 0; theta < in.getSizeTheta(); theta++) {
        for (int r = 0; r < in.getSizeR(); r++) {
            for (int phi = 0; phi < in.getSizePhi(); phi++) // 2->n ??
            {
                auto rR_current = (in.getCell(r,phi,theta) - in.getCell(r-1,phi,theta)) / nonZeroDenom(in.getCell(r+1,phi,theta) - in.getCell(r,phi,theta));
                auto rR_prev = (in.getCell(r-1,phi,theta) - in.getCell(r-2,phi,theta)) / nonZeroDenom(in.getCell(r,phi,theta) - in.getCell(r-1,phi,theta));

                auto rphi_current = (in.getCell(r,phi,theta) - in.getCell(r,phi-1,theta)) / nonZeroDenom(in.getCell(r,phi+1,theta) - in.getCell(r,phi,theta));
                auto rphi_prev = (in.getCell(r,phi-1,theta) - in.getCell(r,phi-2,theta)) / nonZeroDenom(in.getCell(r,phi,theta) - in.getCell(r,phi-1,theta));

                auto rtheta_current = (in.getCell(r,phi,theta) - in.getCell(r,phi,theta-1)) / nonZeroDenom(in.getCell(r,phi,theta+1) - in.getCell(r,phi,theta));
                auto rtheta_prev = (in.getCell(r,phi,theta-1) - in.getCell(r,phi,theta-2)) / nonZeroDenom(in.getCell(r,phi,theta) - in.getCell(r,phi,theta-1));

                auto uL_r = in.getCell(r-1,phi,theta) + 0.5 * SlopeLim(rR_prev) * (in.getCell(r,phi,theta) - in.getCell(r-1,phi,theta));
                auto uR_r = in.getCell(r,phi,theta) - 0.5 * SlopeLim(rR_current) * (in.getCell(r+1,phi,theta) - in.getCell(r,phi,theta));

                auto uL_ph = in.getCell(r,phi-1,theta) + 0.5 * SlopeLim(rphi_prev) * (in.getCell(r,phi,theta) - in.getCell(r,phi-1,theta));
                auto uR_ph = in.getCell(r,phi,theta) - 0.5 * SlopeLim(rphi_current) * (in.getCell(r,phi+1,theta) - in.getCell(r,phi,theta));

                auto uL_th = in.getCell(r,phi,theta-1) + 0.5 * SlopeLim(rtheta_prev) * (in.getCell(r,phi,theta) - in.getCell(r,phi,theta-1));
                auto uR_th = in.getCell(r,phi,theta) - 0.5 * SlopeLim(rtheta_current) * (in.getCell(r,phi,theta+1) - in.getCell(r,phi,theta));

#if defined(USE_CONST_A)
                double aR = A_SPEED;
                double aphi = A_SPEED;
                double atheta = A_SPEED;
#else
                double aR = maxSpeed(0.5 * (uL_R + uR_R));
                double aphi = maxSpeed(0.5 * (uL_phi + uR_phi));
                double atheta = maxSpeed(0.5 * (uL_theta + uR_theta));
#endif

                //std::cout << a << std::endl;
                // F{i-0.5} = 0.5 * (F(uR{i-0.5}) + F(uL{i-0.5}) - a * (uR{i-0.5} - uL{i-0.5}))
                std::get<0>(out).getCellRef(r,phi,theta) =
                        0.5*(FluxR(uR_r)+FluxR(uL_r))-aR*(uR_r-uL_r);
                std::get<1>(out).getCellRef(r,phi,theta) =
                        0.5*(FluxPhi(uR_ph)+FluxPhi(uL_ph))-aR*(uR_ph-uL_ph);
                std::get<2>(out).getCellRef(r,phi,theta) =
                        0.5*(FluxTheta(uR_th)+FluxTheta(uL_th))-aR*(uR_th-uL_th);
            }

        }
    }
}

Cell S(int x,int y,Cell val)
{
    if (x>3 && x<7 && y>32 && y< 58)
        return {0,0,0,0};
    else return {0,0,0,0};
}

void InitialConditions(SphericalGrid& grid) {
    for (int x = 0; x < grid.getSizeR(); x++) {
        for (int y = 0; y < grid.getSizePhi(); y++) {
            double rho=0.1;
            double vx=0;
            double vy=0;
            double vz=0;
            double Bx=0.000;
            double By=0.000;
            double Bz=0.000;
            double T =273;
            // for (int y = 45; y < 55; y++){
            //        for (int x = 20; x < 40; x++){
            if(x>20 && x< 40 && y>45 && y<55)
            {
                rho=0.36*std::cos((x-30.0)/10);
                 vx=1000000;
            }


            //nk=p_rho/m m=1.6733e-27
            //k=1.38044e-23
            double E = 2 * rho * m_div_k * T/(gamma-1)
                       + rho * (vx*vx + vy*vy + vz*vz)
                       + (Bx*Bx + By*By + Bz*Bz) /(2*mu);

            Cell& c =grid.getCellRef(x,y,0);
            c.p_rho = rho;
            c.p_Vr = vx;
            c.p_Vph = vy;
            c.p_Vth = vz;
            c.p_Br = Bx;
            c.p_Bph = By;
            c.p_Bth = Bz;
            c.p_P = 100000;
        }
    }
    grid.UpdateCons();
}


void ApplyBoundaryConditions(SphericalGrid& grid,double t)
{
   /* for (int x=0;x<grid.getSizePhi();x++) {

        double vx=0;
        double vy=0;
        double vz=0;
        double Bx=0.000;
        double By=0.000;
        double Bz=0.000;
        double T =273;
        double rho=0.06/std::pow(grid.getRFromIndex(1),2);
        if(x>40 && x< 45 || x>65 && x< 70)
        {
            rho=0.1;
            vx=1000000;
        }
        double E = 2 * rho * m_div_k * T/(gamma-1)
                   + rho * (vx*vx + vy*vy + vz*vz)
                   + (Bx*Bx + By*By + Bz*Bz) /(2*mu);


        Cell& c = grid.getCellRef(1,x,0);
        c.p_rho=rho;
        c.c_E=E;
        c.p_Vr=vx;
        c.p_Vth=vz;
        c.p_Vph=vy;
    }
    grid.UpdateCons();*/
}



void RKIntegrator(SphericalGrid& grid, double dt,double& t)
{
    t+=dt;

    double dr = grid.getRFromIndex(1)-grid.getRFromIndex(0);
    double dphi = (grid.getPhiFromIndex(1)-grid.getPhiFromIndex(0));
    double dtheta = grid.getThetaFromIndex(1)-grid.getThetaFromIndex(0);

    SphericalGrid fluxR = SphericalGrid::copyGrid(grid);
    SphericalGrid fluxPhi = SphericalGrid::copyGrid(grid);
    SphericalGrid fluxTheta = SphericalGrid::copyGrid(grid);
    std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> flux(fluxR,fluxPhi,fluxTheta);

    grid.UpdatePrim();

    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                Cell& c=grid.getCellRef(r,phi,theta);
                Cell Dr = (grid.getCell(r+1,phi,theta)-grid.getCell(r-1,phi,theta))/(2*dr);
                Cell Dphi = (grid.getCell(r,phi+1,theta)-grid.getCell(r,phi-1,theta))/(2*dphi * grid.getRFromIndex(r));
                Cell Dtheta = (grid.getCell(r,phi,theta+1)-grid.getCell(r,phi,theta-1))/(2*dtheta * grid.getRFromIndex(r));
                c = c -0.5 * dt * F(Dr,Dtheta,Dphi,c);
            }
        }
    }

    grid.UpdateCons();

    SphericalGrid k = SphericalGrid::copyGrid(grid);
    CalculateFlux(flux, grid);

    for (int theta = 0; theta < grid.getSizeTheta(); theta++) {
        for (int r = 0; r < grid.getSizeR(); r++) {
            for (int phi = 0; phi < grid.getSizePhi(); phi++) {
                Cell& c=grid.getCellRef(r,phi,theta);
                c = c - dt * (   (fluxR.getCell(r+1,phi,theta)     - fluxR.getCell(r,phi,theta)) /dr     +
                        (fluxPhi.getCell(r,phi+1,theta)   - fluxPhi.getCell(r,phi,theta))/(dphi * grid.getRFromIndex(r))   +
                        (fluxTheta.getCell(r,phi,theta+1) - fluxTheta.getCell(r,phi,theta))/(dtheta * grid.getRFromIndex(r)));
            }
        }
    }
    ApplyBoundaryConditions(grid,t);

}
