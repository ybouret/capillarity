#include "bridge.hpp"
#include "yocto/math/trigconv.hpp"

double Bridge:: FindTheta(const double height, const double alpha_deg) throw()
{

    const double alpha     = Deg2Rad(alpha_deg);
    const double omega     = lens->omega(alpha);
    const double omega_deg = Rad2Deg(omega);
    double       theta_top = 180.0 - omega_deg;
    if(theta_top<=0)
    {
        return -1;
    }

    double theta_bot = 0;
    if( !FinalRadius(height, theta_bot, alpha_deg) )
    {
        return -1;
    }

    // assuming theta_bot is OK, theta_top is BAD
    while(theta_top - theta_bot>ATOL)
    {
        const double theta_mid = 0.5*(theta_bot+theta_top);
        if(FinalRadius(height, theta_mid, alpha_deg))
        {
            theta_bot = theta_mid;
        }
        else
        {
            theta_top = theta_mid;
        }
    }

    return theta_bot;
}
