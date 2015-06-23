#include "bridge.hpp"
#include "yocto/math/fcn/zfind.hpp"
#include "yocto/math/trigconv.hpp"

double Bridge::FindHmax(const double theta)
{

    double h_lo = 0;
    if( FindAlpha(h_lo, theta) < 0 )
    {
        return -1;
    }

    // assume h_lo is valid
    double h_up = h_lo;
    while(true)
    {
        h_up = h_lo + capillary_length;
        if(FindAlpha(h_up,theta)<0)
        {
            break;
        }
        h_lo = h_up;
    }
    std::cerr << "bracket h_max: " << h_lo << " => " << h_up << std::endl;

    while(h_up-h_lo>HTOL)
    {
        const double h_mid = 0.5*(h_lo+h_up);
        if(FindAlpha(h_mid, theta)<0)
        {
            h_up = h_mid;
        }
        else
        {
            h_lo = h_mid;
        }
    }

    return h_lo;
}