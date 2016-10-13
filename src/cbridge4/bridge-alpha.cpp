#include "bridge.hpp"
#include "yocto/math/opt/bracket.hpp"
#include "yocto/math/opt/minimize.hpp"

#define HAS_LOWER 0x01
#define HAS_UPPER 0x02
#define HAS_BOTH  (HAS_LOWER|HAS_UPPER)


double Bridge:: find_alpha(const double theta, const double zeta, bool &isFlat)
{
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
        isFlat = true;
        return numeric<double>::pi-theta;
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

    std::cerr << "flags=" << flags << std::endl;
    std::cerr << "zeta =" << zeta << " / critical=" << critical_zeta.b << std::endl;
    switch(flags)
    {
        case HAS_LOWER: return find_lower(alpha_lo,Alpha.b,F,resolution);
        case HAS_UPPER: return find_upper(Alpha.b,alpha_up,F,resolution);
        case HAS_BOTH:
        {
            const double lower     = find_lower(alpha_lo,Alpha.b,F,resolution);
            const double upper     = find_upper(Alpha.b,alpha_up,F,resolution);
            bool         localFlat = false;
            const double th_lo     = find_theta(lower, zeta, localFlat);
            const double th_hi     = find_theta(upper, zeta, localFlat);
            std::cerr << "th_lo=" << Rad2Deg(th_lo) << std::endl;
            std::cerr << "th_hi=" << Rad2Deg(th_hi) << std::endl;
            const double delta_lo  = Fabs(th_lo-theta);
            const double delta_hi  = Fabs(th_hi-theta);
            if(delta_hi<delta_lo)
            {
                return upper;
            }
            else
            {
                if(delta_lo<delta_hi)
                {
                    return lower;
                }
                else
                {
                    return 0.5*(lower+upper);
                }
            }
        }

        default:
            break;
    }

    throw exception("find_alpha: invalid flags=%d", flags);

}
