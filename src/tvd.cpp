#include "tvd.h"

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
        //res_p[i] = std::max(0.0, std::min(1.0, r_p[i]));
        res_p[i]=std::max(0.0,std::max(std::min(2*r_p[i],1.0),std::min(r_p[i],2.0)));
        //res_p[i]=std::max(0.0, 1.5 * (r_p[i] * r_p[i] + r_p[i]) / (r_p[i] * r_p[i] + r_p[i] + 1));
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