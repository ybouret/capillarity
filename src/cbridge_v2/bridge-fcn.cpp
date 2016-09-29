#include "bridge.hpp"
#include "yocto/math/opt/bracket.hpp"

double Bridge:: ProfileOfAlpha(const double alpha)
{
    return profile(alpha, __theta, __zeta, NULL);
}


//! assuming F(p_lo)=0, F(p_up)>0
static inline
double __find_top( Function &F, double p_lo, double p_up, const double res)
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

static inline
double __find_bot( Function &F, double p_lo, double p_up, const double res)
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

#define HAS_TOP 0x01
#define HAS_BOT 0x02
#define OUT_ALPHA 1

double Bridge:: find_alpha(const double theta, const double zeta)
{
    __theta = theta;
    __zeta  = zeta;


    int       flag = 0;
    Function &F    = fn_of_alpha;
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

    Triplet alpha = { alpha_lo,0,alpha_up };
    Triplet value = { value_lo,0,value_up };

    bracket<double>::inside(F,alpha,value);
    optimize1D<double>::run(F, alpha, value, resolution);
    
    if(value.b>0)
    {
        return -1; // no minium => no bridge !
    }


#if 1 == OUT_ALPHA
    vfs &fs = local_fs::instance();
    try { fs.remove_file("alpha_top.dat"); } catch(...) {}
    try { fs.remove_file("alpha_bot.dat"); } catch(...) {}
#endif

    double alpha_top = -1;
    double arate_top  = 0;
    if(0!=(flag&HAS_TOP))
    {
        alpha_top = __find_top(F, alpha.b, alpha_up, resolution);
        (void)F(alpha_top);
        //arate_top = Fabs(reduced_rate(param));
        //arate_top = Fabs(sin(param[BRIDGE_A]));
        arate_top = param[BRIDGE_U];

#if 1 == OUT_ALPHA
        ios::wcstream fp("alpha_top.dat");
        profile(alpha_top, theta, zeta, &fp);
        std::cerr << "alpha_top=" << Rad2Deg(alpha_top) << std::endl;
        std::cerr << "arate_top=" << arate_top << std::endl;
#endif
    }

    double alpha_bot = -1;
    double arate_bot =  0;
    if(0!=(flag&HAS_BOT))
    {
        alpha_bot = __find_bot(F,alpha_lo,alpha.b,resolution);
        (void)F(alpha_bot);
        //arate_bot = Fabs(reduced_rate(param));
        //arate_bot = Fabs(sin(param[BRIDGE_A]));
        arate_bot = param[BRIDGE_U];
        
#if 1 == OUT_ALPHA
        ios::wcstream fp("alpha_bot.dat");
        profile(alpha_bot, theta, zeta, &fp);
        std::cerr << "alpha_bot=" << Rad2Deg(alpha_bot) << std::endl;
        std::cerr << "arate_bot=" << arate_bot << std::endl;
#endif
    }

    switch(flag)
    {
        case HAS_BOT: return alpha_bot;
        case HAS_TOP: return alpha_top;
        case HAS_TOP|HAS_BOT:
            if(arate_top>arate_bot)
            {
                return alpha_top;
            }
            else
            {
                if(arate_bot>arate_top)
                {
                    return alpha_bot;
                }
                else
                {
                    return 0.5*(alpha_bot+alpha_top);
                }
            }
        default:
            break;
    }

    return -1;
}

#include "yocto/sort/quick.hpp"
#include "yocto/math/dat/linear.hpp"

double  Bridge:: get_user_rise(const double alpha, const double theta, const double zeta)
{

    const size_t n = heights.size(); assert(n==radii.size());
    std::cerr << "#heights=" << n << std::endl;
    co_qsort(heights,radii);
    slices.make(n);
    volume.make(n);
    const double phi0 = alpha+theta-numeric<double>::pi;
    if(phi0<=0)
    {
        // compute each slice area
        for(size_t i=1;i<=n;++i)
        {
            const double zz = heights[i];
            const double rr = radii[i];
            slices[i] = rr*rr;
            if(zz>zeta)
            {
                const double ri2 = max_of<double>(0,1.0-Square(1.0+zeta-zz));
                slices[i] -= ri2;
            }
        }
        slices[n] = 0;
        volume[1] = 0;
    }
    else
    {
        return 0;
    }

    for(size_t i=2;i<=n;++i)
    {
        volume[i] = volume[i-1] + 0.5*(heights[i]-heights[i-1]) * (slices[i]+slices[i-1]);
    }
    const double vhalf = 0.5*volume[n];
    const double urise = linear(vhalf,volume,heights);
    std::cerr << "urise=" << urise << "/zeta=" << zeta << std::endl;
    return urise;
}

