#include "Derivative.h"
#include "Constants.h"

Cell F(Cell Dr,Cell Dtheta, Cell Dphi, Cell U)
{
    Cell res = U;
    res.p_rho= Dr.p_rho * U.p_Vr + Dphi.p_rho * U.p_Vph
               + U.p_rho * (Dr.p_Vr + Dphi.p_Vph);
    res.p_Vr=0;
    res.p_Vph=0;
    res.p_Vth=0;
    res.p_Br  = 0;
    res.p_Bph = 0;
    res.p_Bth = 0;
    res.p_P = 0;
    return res;
}
