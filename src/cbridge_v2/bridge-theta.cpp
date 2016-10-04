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
    __success            = true;
    bool         is_flat = false;
    const double Xi      = __zeta+shift;
    const double theta   = find_theta(__alpha,Xi,&is_flat);

    if(is_flat) return 0; // the very special case

    if(theta<0)
    {
        __success = false;
        return 0;
    }

    // recompute, storing coordinates
    (void) profile(__alpha,theta,Xi,NULL,true);

    // compute the sinking coefficient
    double usink = compute_user_sink(__alpha,theta,Xi);
    std::cerr << "usink=" << usink << std::endl;

    Bridge::SaveLens("lens.dat", Xi);
    {
        ios::wcstream fp("v2.dat");
        for(size_t i=1;i<=radii.size();++i)
        {
            fp("%g %g\n", radii[i], heights[i]);
        }
        fp << "\n";
        for(double xx=-1;xx<=1;xx+=0.1)
        {
            fp("%g %g\n", xx, usink);
        }
    }

    return (Xi-usink)-__zeta;
}

#include "yocto/math/fcn/zfind.hpp"

double Bridge:: find_theta_v2(double alpha, const double zeta, double &shift0, bool *global_is_flat)
{
    std::cerr << "zeta="  << zeta << std::endl;
    __alpha     = alpha;
    __zeta      = zeta;
    Function &F = delta_of_shift;

    double shift_prev = 0;
    double delta_prev = F(shift_prev);
    if(!__success)
    {
        std::cerr << "Cannot Initialize Search..." << std::endl;
        return -1;
    }
    const double dstep = -delta_prev * 0.1;

    while(true)
    {
        const double shift_next = shift_prev + dstep;
        const double delta_next = F(shift_next);
        std::cerr << "shift_prev=" << shift_prev << ", delta_prev=" << delta_prev << std::endl;
        std::cerr << "shift_next=" << shift_next << ", delta_next=" << delta_next << std::endl;

        if(!__success)
        {
            std::cerr << "Cannot find theta for shift=" << shift_next << std::endl;
            return -1;
        }
        if(delta_next*delta_prev<=0)
        {
            Triplet shift = { shift_prev, 0, shift_next };
            Triplet delta = { delta_prev, 0, delta_next };
            zfind<double> solve( shift_control * 0.1 );
            shift0 = solve.run(F,shift,delta);
            std::cerr << "shift0=" << shift0 << std::endl;
            const double Xi      = zeta+shift0;
            bool         is_flat = false;
            const double theta   = find_theta(alpha,Xi,&is_flat);
            if(global_is_flat) *global_is_flat = is_flat;
            std::cerr << "theta=" << Rad2Deg(theta) << std::endl;
            return theta;
        }
        shift_prev = shift_next;
        delta_prev = delta_next;
    }



#if 0
    double shift = 0;
    double delta = F(shift);
    if(!__success)
    {
        std::cerr << "Cannot Initialize Search..." << std::endl;
        return -1;
    }

    std::cerr << "shift=" << shift << ", delta=" << delta << std::endl;

    while(true)
    {
        shift -= delta;
        const double delta_new  = F(shift);
        if(!__success)
        {
            std::cerr << "Cannot find theta for shift=" << shift << std::endl;
            return -1;
        }
        std::cerr << "delta=" << delta_new << std::endl;
        if(Fabs(delta_new)>=Fabs(delta))
        {
            std::cerr << "STOP" << std::endl;
            break;
        }
        delta = delta_new;
        std::cerr << "shift=" << shift << ", delta=" << delta << std::endl;

    }

    const double Xi      = zeta+shift;
    bool         is_flat = false;
    const double theta   = find_theta(alpha,Xi,&is_flat);
    if(global_is_flat) *global_is_flat = is_flat;
    std::cerr << "theta=" << Rad2Deg(theta) << std::endl;
    return theta;
#endif

}
