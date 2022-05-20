#include "fluxes.h"
#include "Constants.h"


Cell FluxR(Cell U)
{
    double B2 = U.p.Br * U.p.Br + U.p.Bph * U.p.Bph + U.p.Bth * U.p.Bth;
    double Pt = U.p.P + B2/2;
    double vB = U.p.Vr * U.p.Br + U.p.Vph * U.p.Bph + U.p.Vth * U.p.Bth;
    Cell res  = U;
    res.c.m =    U.p.Vr * U.p.rho;
    res.c.Mr  = (U.p.Vr * U.p.Vr  * U.p.rho - (U.p.Br * U.p.Br)/mu /*+ Pt*/);
    res.c.Mph = (U.p.Vr * U.p.Vph * U.p.rho - (U.p.Br * U.p.Bph)/mu);
    res.c.Mth = (U.p.Vr * U.p.Vth * U.p.rho - (U.p.Br * U.p.Bth)/mu);

    res.c.Br = 0.0;
    res.c.Bph = U.p.Vr*U.p.Bph - U.p.Br * U.p.Vph;
    res.c.Bth = U.p.Vr*U.p.Bth - U.p.Br * U.p.Vth;
    res.c.E = ((U.c.E/U.volume + Pt)*U.p.Vr - (U.p.Br*vB)/mu);
    return res;
}

Cell FluxPhi(Cell U)
{
    double B2 = U.p.Br * U.p.Br + U.p.Bph * U.p.Bph + U.p.Bth * U.p.Bth;
    double Pt = U.p.P + B2/2;
    double vB = U.p.Vr * U.p.Br + U.p.Vph * U.p.Bph + U.p.Vth * U.p.Bth;
    Cell res  = U;
    res.c.m = U.p.Vph * U.p.rho;
    res.c.Mr  = (U.p.Vr  * U.p.Vph * U.p.rho  - (U.p.Bph * U.p.Br)/mu);
    res.c.Mph =(U.p.Vph * U.p.Vph * U.p.rho  - (U.p.Bph * U.p.Bph)/mu /*+ Pt*/);
    res.c.Mth = (U.p.Vph * U.p.Vth * U.p.rho  - (U.p.Bph * U.p.Bth)/mu);

    //n - ph, t - r b - th
    res.c.Br  = U.p.Vph * U.p.Br  - U.p.Bph * U.p.Vr;
    res.c.Bph = 0.0;
    res.c.Bth = U.p.Vph * U.p.Bth - U.p.Bph * U.p.Vth;
    res.c.E = ((U.c.E/U.volume + Pt) * U.p.Vph - (U.p.Bph * vB)/mu);
    //res.c.E = (U.c.E + U.p.P)*U.p.Vph - U.p.Bph*vB;
    return res;
}


Cell FluxTheta(Cell U)
{
    double B2 = U.p.Br * U.p.Br + U.p.Bph * U.p.Bph + U.p.Bth * U.p.Bth;
    double Pt = U.p.P + B2/2;
    double vB = U.p.Vr * U.p.Br + U.p.Vph * U.p.Bph + U.p.Vth * U.p.Bth;
    Cell res  = U;
    res.c.m =  U.p.Vth * U.p.rho;
    res.c.Mr  = U.p.Vr  * U.p.Vth * U.p.rho - (U.p.Bth * U.p.Br )/mu;
    res.c.Mph = U.p.Vph * U.p.Vth * U.p.rho - (U.p.Bth * U.p.Bph)/mu;
    res.c.Mth = (U.p.Vth * U.p.Vth * U.p.rho - (U.p.Bth * U.p.Bth)/mu  /* + Pt*/);

    //n - th t - r b -ph
    res.c.Br  = U.p.Vth * U.p.Br  - U.p.Bth * U.p.Vr;
    res.c.Bph = U.p.Vth * U.p.Bph - U.p.Bth * U.p.Vph;
    res.c.Bth = 0.0;
    res.c.E = (U.c.E/U.volume + Pt)*U.p.Vth - (U.p.Bth * vB)/mu;
    return res;
}

Cell DiffusiveTerm(Cell R, Cell L)
{
    Cell res  = R;

    res.c.m = (R.p.rho - L.p.rho);
    res.c.Mr  = (R.p.rho * R.p.Vr - L.p.rho * L.p.Vr);
    res.c.Mph = (R.p.rho * R.p.Vph - L.p.rho * L.p.Vph);
    res.c.Mth = (R.p.rho * R.p.Vth - L.p.rho * L.p.Vth);

    res.c.Br = (R.p.Br - L.p.Br);
    res.c.Bph = (R.p.Bph - L.p.Bph);
    res.c.Bth = (R.p.Bth - L.p.Bth);
    res.c.E = (R.c.E/R.volume - L.c.E/L.volume);

    return res;

}


void CalculateFlux(std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> out, SphericalGrid& in, std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> grad)
{
    //std::cout << "---------------------------------------------------" << std::endl;
    for(int theta = 0; theta < in.getSizeTheta(); theta++) {
        for (int r = 0; r < in.getSizeR(); r++) {
            for (int phi = 0; phi < in.getSizePhi(); phi++) // 2->n ??
            {
                auto Dr       = std::get<T_R>(grad).getCell(r,phi,theta);
                auto Dr_1     = std::get<T_R>(grad).getCell(r-1,phi,theta);

                auto Dphi     = std::get<T_PHI>(grad).getCell(r,phi,theta);
                auto Dphi_1   = std::get<T_PHI>(grad).getCell(r,phi-1,theta);

                auto Dtheta   = std::get<T_THETA>(grad).getCell(r,phi,theta);
                auto Dtheta_1 = std::get<T_THETA>(grad).getCell(r,phi,theta-1);

                auto uL_r = in.getCell(r-1,phi,theta)  + 0.5 * Dr_1;
                auto uR_r = in.getCell(r,phi,theta)    - 0.5 * Dr;

                auto uL_ph = in.getCell(r,phi-1,theta) + 0.5 * Dphi_1;
                auto uR_ph = in.getCell(r,phi,theta)   - 0.5 * Dphi;

                auto uL_th = in.getCell(r,phi,theta-1) + 0.5 * Dtheta_1;
                auto uR_th = in.getCell(r,phi,theta)   - 0.5 * Dtheta;

                double dphi=M_PI/in.getSizePhi();
                double dtheta =(in.getThetaFromIndex(1)-in.getThetaFromIndex(0))/2;

                uL_ph.p = uL_ph.p.rotate(dphi,0);
                uR_ph.p = uR_ph.p.rotate(-dphi,0);

                uL_th.p = uL_th.p.rotate(0,dtheta);
                uR_th.p = uR_th.p.rotate(0,-dtheta);

                uL_r.UpdateCons();
                uR_r.UpdateCons();
                uL_ph.UpdateCons();
                uR_ph.UpdateCons();
                uL_th.UpdateCons();
                uR_th.UpdateCons();

                double C_uL_r = std::sqrt(gamma*uL_r.p.P/uL_r.p.rho) + std::abs(uL_r.p.Vr);
                double C_uR_r = std::sqrt(gamma*uR_r.p.P/uR_r.p.rho) + std::abs(uR_r.p.Vr);
                double C_r = std::max(C_uL_r, C_uR_r);

                double C_uL_ph = std::sqrt(gamma*uL_ph.p.P/uL_ph.p.rho) + std::abs(uL_ph.p.Vph);
                double C_uR_ph = std::sqrt(gamma*uR_ph.p.P/uR_ph.p.rho) + std::abs(uR_ph.p.Vph);
                double C_ph = std::max(C_uL_ph, C_uR_ph);

                double C_uL_th = std::sqrt(gamma*uL_th.p.P/uL_th.p.rho) + std::abs(uL_th.p.Vth);
                double C_uR_th = std::sqrt(gamma*uR_th.p.P/uR_th.p.rho) + std::abs(uR_th.p.Vth);
                double C_th = std::max(C_uL_th, C_uR_th);

               // maxSpeed(c);
#if defined(USE_CONST_A)
                double aR = A_SPEED;
                double aphi = A_SPEED;
                double atheta = A_SPEED;
#else
                double aR = C_r;
                double aphi = C_ph;
                double atheta = C_th;

#endif

                //std::cout << a << std::endl;
                // F{i-0.5} = 0.5 * (F(uR{i-0.5}) + F(uL{i-0.5}) - a * (uR{i-0.5} - uL{i-0.5}))
                std::get<T_R>(out).getCellRef(r,phi,theta) =
                        0.5 * (FluxR(uR_r)+FluxR(uL_r) - aR * DiffusiveTerm(uR_r, uL_r)) * uR_r.Sr;
                std::get<T_PHI>(out).getCellRef(r,phi,theta) =
                        0.5 * (FluxPhi(uR_ph)+FluxPhi(uL_ph) - aphi * DiffusiveTerm(uR_ph, uL_ph)) * uR_r.Sph;
                std::get<T_THETA>(out).getCellRef(r,phi,theta) =
                        0.5 * (FluxTheta(uR_th)+FluxTheta(uL_th) - atheta * DiffusiveTerm(uR_th, uL_th)) * uR_r.Sth;
            }

        }
    }
}


