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

double Bridge:: DeltaOfShift(const double shift)
{
    __success = true;
    bool is_flat = false;
    const double Xi    = __zeta+shift;
    const double theta = find_theta(__alpha,Xi,&is_flat);

    if(is_flat) return 0;

    if(theta<0)
    {
        __success = false;
        return 0;
    }
    
    (void) profile(__alpha,theta,Xi,NULL,true);
    double usink = compute_user_sink(__alpha,theta,Xi);
    return Xi-usink;
}

double Bridge:: find_theta_v2(double alpha, const double zeta, bool *is_flat)
{
    std::cerr << "zeta="  << zeta << std::endl;
    __alpha     = alpha;
    __zeta      = zeta;
    Function &F = delta_of_shift;

    double shift = 0;
    double dcurr = F(shift);
   // int    count = 0;
    std::cerr << "delta0=" << dcurr << std::endl;
    while(true)
    {
        shift -= dcurr;
        const double delta = F(shift);
        std::cerr << "shift=" << shift << std::endl;
        std::cerr << "delta=" << delta << std::endl;

        if(Fabs(delta)>=Fabs(dcurr))
        {
            break;
        }
        dcurr = delta;
    }

    {
        ios::wcstream fp("theta2.dat");
        const double Xi    = __zeta+shift;
        std::cerr << "Xi=" << Xi << std::endl;
        const double theta = find_theta(__alpha,Xi);
        (void) profile(__alpha,theta,Xi,NULL,true);
        double usink = compute_user_sink(__alpha,theta,Xi);
        (void)usink;
        profile(alpha, theta, Xi, &fp,false,-shift);
        SaveLens("shiftlens.dat",zeta);
    }

    return -1;
}
