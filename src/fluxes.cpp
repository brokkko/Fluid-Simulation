#include "fluxes.h"
#include "Constants.h"




Cell FluxR(Cell U)
{
    double vB = U.p_Vr * U.p_Br + U.p_Vph * U.p_Bph + U.p_Vth * U.p_Bth;
    Cell res  = U;

    res.c_rho =  U.p_Vr * U.p_rho * U.volume;
    res.c_Mr  = (U.p_Vr * U.p_Vr * U.p_rho  - U.p_Br * U.p_Br + U.p_P) * U.volume;
    res.c_Mph = (U.p_Vph * U.p_Vr * U.p_rho - U.p_Br * U.p_Bph) * U.volume;
    res.c_Mth = (U.p_Vth * U.p_Vr * U.p_rho - U.p_Br * U.p_Bth) * U.volume;

    res.c_Br = 0;
    res.c_Bph = 0;
    res.c_Bth = 0;
    res.c_E = ((U.c_E/U.volume + U.p_P)*U.p_Vr - U.p_Br*vB)*U.volume;
    return res;
}

Cell FluxPhi(Cell U)
{
    double vB = U.p_Vr * U.p_Br + U.p_Vph * U.p_Bph + U.p_Vth * U.p_Bth;
    Cell res  = U;
    res.c_rho = U.p_Vph * U.p_rho * U.volume;
    res.c_Mr  = (U.p_Vr  * U.p_rho * U.p_Vph - U.p_Bph * U.p_Br)* U.volume ;
    res.c_Mph = (U.p_Vph * U.p_rho * U.p_Vph - U.p_Bph * U.p_Bph + U.p_P)* U.volume;
    res.c_Mth = (U.p_Vth * U.p_rho * U.p_Vph - U.p_Bph * U.p_Bth)* U.volume;


    res.c_Br = 0;
    res.c_Bph = 0;
    res.c_Bth = 0;
    res.c_E = ((U.c_E/U.volume + U.p_P)*U.p_Vph - U.p_Bph*vB)*U.volume;
    //res.c_E = (U.c_E + U.p_P)*U.p_Vph - U.p_Bph*vB;
    return res;
}


//TODO: Multiply by volume
Cell FluxTheta(Cell U)
{
    double vB = U.p_Vr * U.p_Br + U.p_Vph * U.p_Bph + U.p_Vth * U.p_Bth;
    Cell res  = U;
    res.c_rho = U.c_Mth;
    res.c_Mr  = U.c_Mr  * U.p_Vth - U.p_Bth * U.p_Br ;
    res.c_Mph = U.c_Mph * U.p_Vth - U.p_Bth * U.p_Bph;
    res.c_Mth = U.c_Mth * U.p_Vth - U.p_Bth * U.p_Bth + U.p_P;

    res.c_Br = 0;
    res.c_Bph = 0;
    res.c_Bth = 0;
    res.c_E = (U.c_E + U.p_P)*U.p_Vth - U.p_Bth*vB;
    return res;
}


void CalculateFlux(std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> out, SphericalGrid& in, std::tuple<SphericalGrid&,SphericalGrid&,SphericalGrid&> grad)
{
    for(int theta = 0; theta < in.getSizeTheta(); theta++) {
        for (int r = 0; r < in.getSizeR(); r++) {
            for (int phi = 0; phi < in.getSizePhi(); phi++) // 2->n ??
            {
               /* auto Dr       = std::get<T_R>(grad).getCell(r,phi,theta);
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
                auto uR_th = in.getCell(r,phi,theta)   - 0.5 * Dtheta;*/

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
                std::get<T_R>(out).getCellRef(r,phi,theta) =
                        0.5*(FluxR(uR_r)+FluxR(uL_r))-aR*(uR_r-uL_r);
                std::get<T_PHI>(out).getCellRef(r,phi,theta) =
                        0.5*(FluxPhi(uR_ph)+FluxPhi(uL_ph))-aphi*(uR_ph-uL_ph);
                std::get<T_THETA>(out).getCellRef(r,phi,theta) =
                        0.5*(FluxTheta(uR_th)+FluxTheta(uL_th))-atheta*(uR_th-uL_th);
            }

        }
    }
}


