#include "bridge.hpp"
#include "yocto/math/opt/bracket.hpp"

double Bridge:: ProfileOfAlpha(const double alpha)
{
    return profile(alpha, __theta, __zeta, NULL);
}

double Bridge:: find_alpha(const double theta, const double zeta)
{
    __theta = theta;
    __zeta  = zeta;
    Function &F = fn_of_alpha;
    const double alpha_lo = resolution;
    const double value_lo = F(alpha_lo);

    const double alpha_up = numeric<double>::pi-resolution;
    const double value_up = F(alpha_up);

    Triplet alpha = { alpha_lo,0,alpha_up };
    Triplet value = { value_lo,0,value_up };

    bracket<double>::inside(F,alpha,value);
    std::cerr << "alpha=" << alpha << std::endl;
    std::cerr << "value=" << value << std::endl;

    if(value.b>0)
    {
        return -1; // no minium => no bridge !
    }

    return -1;
}
