#include "bridge.hpp"

#include "yocto/math/opt/bracket.hpp"
#include "yocto/math/opt/minimize.hpp"

double Bridge:: ProfileOfTheta(const double theta)
{
    return profile(__alpha, theta, __zeta, NULL);
}


#define HAS_LOWER 0x01
#define HAS_UPPER 0x02
#define HAS_BOTH  (HAS_LOWER|HAS_UPPER)

double Bridge:: find_theta(const double alpha, const double zeta, ios::ocstream *fp)
{
    assert(zeta>=-2);
    assert(alpha>0);
    assert(alpha<numeric<double>::pi);

    // setup
    __alpha     = alpha;
    __zeta      = zeta;
    Function &F = prof_theta;

    // initialize
    int    flags    = 0;
    double theta_lo = resolution;
    double value_lo = F(theta_lo);
    if(value_lo>0)
    {
        flags |= HAS_LOWER;
    }

    double theta_hi = numeric<double>::pi - resolution;
    double value_hi = F(theta_hi);
    if(value_hi>0)
    {
        flags |= HAS_UPPER;
    }


    if(fp)
    {
        (*fp)("%g %g %g %g\n", alpha, zeta, value_lo, value_hi);
    }

    if(0==flags)
    {
        throw exception("Invalid settings\n");
    }

    // finding a minimum
    Triplet Theta = { theta_lo, 0, theta_hi };
    Triplet Value = { value_lo, 0, value_hi };
    bracket<double>::inside(F, Theta, Value);
    optimize1D<double>::run(F, Theta, Value, resolution);

    if(Value.b<=0 && HAS_BOTH==flags)
    {
        std::cerr << "HAS_BOTH for alpha=" << Rad2Deg(alpha) << ", zeta=" << zeta << std::endl;
    }


    return 0;
}
