#include "bridge.hpp"
#include "yocto/math/fcn/zfind.hpp"
#include "yocto/math/trigconv.hpp"

#include "yocto/ios/ocstream.hpp"

void  Bridge:: ScanAlpha(const double height, const double theta)
{
    ios::ocstream fp( vformat("alpha-theta%g.dat",theta),false);
    for(int alpha=1;alpha<=180;++alpha)
    {
        double res = 0;
        if( FinalRadius(height,theta,alpha) )
        {
            res = 1.0 / Y[1];
        }
        fp("%g %g\n", double(alpha), res);
    }
}


double Bridge:: FindAlpha(const double height,
                          const double theta) throw()
{
    const double      omega_max = 180-theta;
    zfunction<double> zfn( lens->omega, Deg2Rad(omega_max) );
    zfind<double>     solve(ATOL);
    double alpha_max = Rad2Deg(solve(zfn.call,0,numeric<double>::pi));
    //std::cerr << "alpha_max=" << alpha_max << std::endl;

    double alpha_min = alpha_max;
    while(true)
    {
        alpha_min = alpha_max/2;
        //std::cerr << "\talpha_min=" << alpha_min << std::endl;
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
    //std::cerr << "\t\tfinal: " << Y[1] << std::endl;

    //std::cerr << "bracketed alpha: " << alpha_min << " => " << alpha_max << std::endl;

    // alpha_min: valid, alpha_max: invalid
    while(alpha_max-alpha_min>ATOL)
    {
        const double alpha_mid = 0.5*(alpha_min+alpha_max);
        if(FinalRadius(height, theta, alpha_mid))
        {
            //std::cerr << "Good: " << Y[1] << std::endl;
            alpha_min = alpha_mid;
        }
        else
        {
            //std::cerr << "Bad : " << Y[1] << std::endl;
            alpha_max = alpha_mid;
        }
    }
    const double alpha = alpha_min;
    return alpha;
}
