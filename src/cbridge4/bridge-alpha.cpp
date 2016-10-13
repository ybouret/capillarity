#include "bridge.hpp"
#include "yocto/math/opt/bracket.hpp"
#include "yocto/math/opt/minimize.hpp"

#define HAS_LOWER 0x01
#define HAS_UPPER 0x02
#define HAS_BOTH  (HAS_LOWER|HAS_UPPER)


size_t Bridge:: find_alpha(const double theta, const double zeta, double *alphas, bool &isFlat)
{
    assert(alphas!=NULL);
    assert(theta>0);
    assert(theta<numeric<double>::pi);

    assert(zeta>=-2);

    // check where we are
    const double half_res      = resolution*0.5;
    Triplet      critical_zeta = { CriticalZetaOfTheta(theta-half_res), CriticalZetaOfTheta(theta), CriticalZetaOfTheta(theta+half_res) };
    critical_zeta.sort();
    const double zeta_lower = critical_zeta.a;
    const double zeta_upper = critical_zeta.c;

    if(zeta>=zeta_lower && zeta<=zeta_upper )
    {
        alphas[0] =  numeric<double>::pi-theta;
        return 1;
    }


    // setup
    __theta     = theta;
    __zeta      = zeta;
    Function &F = prof_alpha;
    isFlat      = false;

    // find possibilities
    int          flags    = 0;
    const double alpha_lo = resolution;
    const double value_lo = F(alpha_lo);
    if(value_lo>0)
    {
        flags |= HAS_LOWER;
    }

    const double alpha_up = numeric<double>::pi-resolution;
    const double value_up = F(alpha_up);
    if(value_up>0)
    {
        flags |= HAS_UPPER;
    }

    Triplet Alpha = { alpha_lo, 0, alpha_up };
    Triplet Value = { value_lo, 0, value_up };

    bracket<double>::inside(F, Alpha, Value);
    optimize1D<double>::run(F, Alpha, Value, resolution);

    //std::cerr << "min=" << Value.b << std::endl;
    if(Value.b>0)
    {
        // no possible bridge
        return 0;
    }

    //std::cerr << "flags=" << flags << std::endl;
    //std::cerr << "zeta =" << zeta << " / critical=" << critical_zeta.b << std::endl;
    switch(flags)
    {
        case HAS_LOWER: alphas[0] = find_lower(alpha_lo,Alpha.b,F,resolution); return 1;
        case HAS_UPPER: alphas[0] = find_upper(Alpha.b,alpha_up,F,resolution); return 1;
        case HAS_BOTH:
        {
            alphas[0] = find_lower(alpha_lo,Alpha.b,F,resolution);
            alphas[1] = find_upper(Alpha.b,alpha_up,F,resolution);
            return 2;
        }

        default:
            break;
    }

    throw exception("find_alpha: invalid flags=%d", flags);

}
