#include "bridge.hpp"

#include "yocto/math/opt/bracket.hpp"
#include "yocto/math/opt/minimize.hpp"

double Bridge:: ProfileOfTheta(const double theta)
{
    return profile(__alpha, theta, __zeta, NULL);
}

double Bridge:: ProfileOfAlpha(const double alpha)
{
    return profile(alpha, __theta, __zeta, NULL);
}


#define HAS_LOWER 0x01
#define HAS_UPPER 0x02
#define HAS_BOTH  (HAS_LOWER|HAS_UPPER)

//! assuming F(lo)>0, F(hi)<=0
double Bridge:: find_lower(double lo,double hi, Function &F, const double resolution)
{
    assert(hi>=lo);
    while(hi-lo>resolution)
    {
        const double mid = 0.5*(hi+lo);
        if(F(mid)<=0)
        {
            hi = mid;
        }
        else
        {
            lo = mid;
        }
    }
    return hi;
}

//! assuming F(lo)<=0, F(hi)>0
double Bridge:: find_upper(double lo,double hi, Function &F, const double resolution)
{
    assert(hi>=lo);
    while(hi-lo>resolution)
    {
        const double mid = 0.5*(hi+lo);
        if(F(mid)<=0)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }
    return lo;
}

double Bridge:: find_theta(const double alpha, const double zeta, bool &isFlat)
{
    assert(zeta>=-2);
    assert(alpha>0);
    assert(alpha<numeric<double>::pi);

    //__________________________________________________________________________
    //
    //
    // computing critical zeta to know what to do
    //
    //__________________________________________________________________________
    const double half_res      = 0.5*resolution;
    Triplet      critical_zeta = { CriticalZetaOfAlpha(alpha+half_res), CriticalZetaOfAlpha(alpha), CriticalZetaOfAlpha(alpha-half_res) };
    critical_zeta.sort();
    const double upper_zeta = critical_zeta.c;
    const double lower_zeta = critical_zeta.a;
    assert(upper_zeta>=lower_zeta);

    if(zeta>lower_zeta&&zeta<upper_zeta)
    {
        isFlat = true;
        return numeric<double>::pi-alpha;
    }

    // setup
    __alpha     = alpha;
    __zeta      = zeta;
    Function &F = prof_theta;
    isFlat      = false;

    // initialize
    int    flags    = 0;
    double theta_lo = resolution;
    double value_lo = F(theta_lo);
    if(value_lo>0)
    {
        flags |= HAS_LOWER;
    }

    double theta_hi = numeric<double>::pi - resolution;
    double value_hi = F(theta_hi);
    if(value_hi>0)
    {
        flags |= HAS_UPPER;
    }
    
    if(0==flags)
    {
        throw exception("Invalid settings\n");
    }

    // finding a minimum
    Triplet Theta = { theta_lo, 0, theta_hi };
    Triplet Value = { value_lo, 0, value_hi };
    bracket<double>::inside(F, Theta, Value);
    optimize1D<double>::run(F, Theta, Value, resolution);




    if(zeta>upper_zeta)
    {
        if(0==(flags&HAS_UPPER))
        {
            throw exception("find_theta invalid setting: no valid upper value...");
        }
        if(Value.b>0)
        {
            return 0; // no bridge
        }
        else
        {
            return find_upper(Theta.b, theta_hi, F, resolution);
        }
    }
    else
    {
        if(zeta<lower_zeta)
        {
            if(0==(flags&HAS_LOWER))
            {
                throw exception("find_theta invalid setting: no valid upper value...");
            }

            if(Value.b>0)
            {
                return 0; // no bridge
            }
            else
            {
                return find_lower(theta_lo,Theta.b,F,resolution);
            }
        }
        else
        {
            // shouldn't happen
            isFlat = true;
            return numeric<double>::pi-alpha;
        }
    }
    
}

double Bridge:: find_zeta_max(const double theta)
{
    bool   isFlat  = false;
    double zeta_lo = 0;
    double alphas[2];

    if( find_alpha(theta,zeta_lo, alphas, isFlat) <= 0 )
    {
        throw exception("find_zeta_max: no solution for zeta=0");
    }
    //std::cerr << "alpha(" << Rad2Deg(theta) << "," << zeta_lo << ")=" << alpha << std::endl;
    double zeta_hi = 1.0;
    while( find_alpha(theta,zeta_hi,alphas,isFlat) > 0 )
    {
        zeta_hi += 0.5;
    }

    const double ztol = numeric<double>::sqrt_ftol;

    while( Fabs(zeta_hi-zeta_lo) > ztol* ( Fabs(zeta_hi)+Fabs(zeta_lo) ) )
    {
        const double zeta_mid = 0.5*(zeta_lo+zeta_hi);
        if( find_alpha(theta,zeta_mid,alphas,isFlat) >0)
        {
            zeta_lo = zeta_mid;
            //std::cerr << "alpha(" << Rad2Deg(theta) << "," << zeta_lo << ")=" << alpha << std::endl;
        }
        else
        {
            zeta_hi = zeta_mid;
        }
    }

    
    return zeta_lo;
}

