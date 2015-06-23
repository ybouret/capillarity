#include "bridge.hpp"
#include "yocto/math/fcn/zfind.hpp"
#include "yocto/math/trigconv.hpp"

double Bridge:: FindAlpha(const double height,
                          const double theta)
{
    const double      omega_max = 180-theta;
    zfunction<double> zfn( lens->omega, omega_max );
    zfind<double>     solve(ATOL);
    double alpha_max = solve(zfn.call,0,180);
    std::cerr << "alpha_max=" << alpha_max << std::endl;

    double alpha_min = alpha_max;
    while(true)
    {
        alpha_min = alpha_max/2;
        std::cerr << "probe alpha_min=" << alpha_min << std::endl;
        if(alpha_min<ATOL)
        {
            return -1;
        }

        if(FinalRadius(height, theta, alpha_min))
        {
            break;
        }
        alpha_max = alpha_min;
    }
    std::cerr << "bracketed alpha: " << alpha_min << " => " << alpha_max << std::endl;

    // alpha_min: valid, alpha_max: invalid
    while(alpha_max-alpha_min>ATOL)
    {
        const double alpha_mid = 0.5*(alpha_min+alpha_max);
        if(FinalRadius(height, theta, alpha_mid))
        {
            alpha_min = alpha_mid;
        }
        else
        {
            alpha_max = alpha_mid;
        }
    }
    const double alpha = alpha_min;
    return alpha;
}
