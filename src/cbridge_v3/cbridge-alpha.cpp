#include "cbridge.hpp"
#include "yocto/math/opt/bracket.hpp"

double Bridge:: ProfileOfAlpha(const double alpha)
{
    return profile(alpha, __theta, __zeta, NULL);
}



//! assuming F(p_lo)=0, F(p_up)>0
double Bridge:: __find_top( Function &F, double p_lo, double p_up, const double res)
{
    assert(p_up>=p_lo);
    while(p_up-p_lo>res)
    {
        const double p_mid = 0.5*(p_up+p_lo);
        const double f_mid = F(p_mid);
        if(f_mid<=0)
        {
            p_lo = p_mid;
        }
        else
        {
            p_up = p_mid;
        }
    }
    return p_lo;
}

//! assuming F(p_lo)>0, F(p_up)=0

double Bridge::__find_bot( Function &F, double p_lo, double p_up, const double res)
{
    assert(p_up>=p_lo);
    while(p_up-p_lo>res)
    {
        const double p_mid = 0.5*(p_up+p_lo);
        const double f_mid = F(p_mid);
        if(f_mid<=0)
        {
            p_up = p_mid;
        }
        else
        {
            p_lo = p_mid;
        }
    }
    return p_up;

}

#include "yocto/fs/local-fs.hpp"
#include "yocto/ios/ocstream.hpp"

#define HAS_TOP   0x01
#define HAS_BOT   0x02
#define OUT_ALPHA 1

double Bridge:: find_alpha(const double theta,
                           const double zeta,
                           bool        *is_flat)
{
    __theta =  theta;
    __zeta  =  zeta;
    if(is_flat) *is_flat = false;

    int       flag = 0;
    Function &F    = prof_alpha;
    const double alpha_lo = resolution;
    const double value_lo = F(alpha_lo);
    if(value_lo>0)
    {
        flag |= HAS_BOT;
    }


    const double alpha_up = numeric<double>::pi-resolution;
    const double value_up = F(alpha_up);

    if(value_up>0)
    {
        flag |= HAS_TOP;
    }

    if(!flag)
    {
        throw exception("Bridge.find_alpha: invalid setting!");
    }

    //__________________________________________________________________________
    //
    //
    // bracket potential minimun
    //
    //__________________________________________________________________________
    Triplet alpha = { alpha_lo,0,alpha_up };
    Triplet value = { value_lo,0,value_up };

    bracket<double>::inside(F,alpha,value);
    optimize1D<double>::run(F,alpha,value,resolution);


    //__________________________________________________________________________
    //
    //
    // computing critical zeta to know what to do
    //
    //__________________________________________________________________________
    Triplet critical_zeta = { CriticalZetaOfTheta(theta+resolution/2), CriticalZetaOfTheta(theta), CriticalZetaOfTheta(theta-resolution/2) };
    critical_zeta.sort();
    const double upper_zeta = critical_zeta.c;
    const double lower_zeta = critical_zeta.a;
    assert(upper_zeta>=lower_zeta);

    //__________________________________________________________________________
    //
    //
    // check cases
    //
    //__________________________________________________________________________
    if(zeta>upper_zeta)
    {
        //______________________________________________________________________
        //
        // check sanity
        //______________________________________________________________________
        if(0==(flag&HAS_TOP))
        {
            throw exception("Bridge.find_alpha: should have a TOP value!");
        }

        //______________________________________________________________________
        //
        // do we intercept ?
        //______________________________________________________________________
        if(value.b>0)
        {
            return -1; // too high !
        }

        return  __find_top(F, alpha.b, alpha_up, resolution);
    }
    else
    {
        if(zeta<lower_zeta)
        {

            //______________________________________________________________________
            //
            // check sanity
            //______________________________________________________________________
            if(0==(flag&HAS_BOT))
            {
                throw exception("Bridge.find_alpha: should have a BOTTOM value!");
            }

            if(value.b>0)
            {
                return -1; // too low, shouldn't happen !
            }

            return __find_bot(F,alpha_lo,alpha.b,resolution);
        }
        else
        {
            // planar !
            if(is_flat) *is_flat = true;
            return numeric<double>::pi-theta;
        }
    }
    
    
    
    return -1;
}

double Bridge:: find_zeta_max(const double theta)
{
    // a good value
    double zeta_lo = 0.0;
    if( find_alpha(theta, zeta_lo, NULL) < 0 )
    {
        throw exception("Cannot find alpha for zeta=0");
    }

    // a bad balue
    double zeta_up = 1.0;
    while(find_alpha(theta,zeta_up,NULL)>=0)
    {
        zeta_lo = zeta_up;
        zeta_up += 1.0;
    }
    std::cerr << "Looking up between " << zeta_lo << " and " << zeta_up << std::endl;
    const double zeta_ftol = numeric<double>::sqrt_ftol;
    while( Fabs(zeta_up-zeta_lo) >= zeta_ftol * ( Fabs(zeta_lo) + Fabs(zeta_up) ) )
    {
        const double zeta_mid = 0.5*(zeta_lo+zeta_up);
        if(find_alpha(theta,zeta_mid,NULL) < 0 )
        {
            zeta_up = zeta_mid;
        }
        else
        {
            zeta_lo = zeta_mid;
        }
    }

    // highest good value
    return zeta_lo;
}
