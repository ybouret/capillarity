#include "bridge.hpp"
#include "yocto/math/opt/bracket.hpp"

double Bridge:: ProfileOfTheta(const double theta)
{
    return profile(__alpha, theta, __Xi, NULL );
}

#define HAS_TOP   0x01
#define HAS_BOT   0x02

double Bridge:: find_theta(const double alpha,
                           const double Xi,
                           bool        *is_flat)
{
    __alpha        = alpha;
    __Xi           = Xi;
    Function &F    = fn_of_theta;
    int       flag = 0;

    if(is_flat) *is_flat = false;

    const double theta_lo = resolution;
    const double value_lo = F(theta_lo);
    if(value_lo>0)
    {
        flag |= HAS_BOT;
    }


    const double theta_up = numeric<double>::pi-resolution;
    const double value_up = F(theta_up);

    if(value_up>0)
    {
        flag |= HAS_TOP;
    }

    if(!flag)
    {
        throw exception("Bridge.find_theta: invalid setting!");
    }
    
    //__________________________________________________________________________
    //
    //
    // bracket potential minimun
    //
    //__________________________________________________________________________
    Triplet theta = { theta_lo,0,theta_up };
    Triplet value = { value_lo,0,value_up };

    bracket<double>::inside(F,theta,value);
    optimize1D<double>::run(F,theta,value,resolution);

    //__________________________________________________________________________
    //
    //
    // computing critical Xi to know what to do
    //
    //__________________________________________________________________________
    Triplet critical_Xi = { CriticalAlphaXi(alpha+resolution/2), CriticalAlphaXi(alpha), CriticalAlphaXi(alpha-resolution/2) };
    critical_Xi.sort(); // just to be sure
    const double upper_Xi = critical_Xi.c;
    const double lower_Xi = critical_Xi.a;
    assert(upper_Xi>=lower_Xi);

    //__________________________________________________________________________
    //
    //
    // check cases
    //
    //__________________________________________________________________________

    if(Xi>upper_Xi)
    {
        //______________________________________________________________________
        //
        // check sanity
        //______________________________________________________________________
        if(0==(flag&HAS_TOP))
        {
            throw exception("Bridge.find_theta: should have a TOP value!");
        }

        //______________________________________________________________________
        //
        // do we intercept ?
        //______________________________________________________________________
        if(value.b>0)
        {
            return -1; // too high !
        }

        return  __find_top(F, theta.b, theta_up, resolution);
    }
    else
    {
        if(Xi<lower_Xi)
        {

            //______________________________________________________________________
            //
            // check sanity
            //______________________________________________________________________
            if(0==(flag&HAS_BOT))
            {
                throw exception("Bridge.find_theta: should have a BOTTOM value!");
            }

            if(value.b>0)
            {
                return -1; // too low, shouldn't happen !
            }

            return __find_bot(F,theta_lo,theta.b,resolution);
        }
        else
        {
            // planar !
            if(is_flat) *is_flat = true;
            return numeric<double>::pi-alpha;
        }
    }
    

    return -1;
}
