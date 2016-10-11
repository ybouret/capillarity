#include "cbridge.hpp"

double Bridge:: DeltaOfZeta(const double zeta)
{
    __success    = true;
    double theta = find_theta(__alpha,zeta,NULL);
    if(theta<0)
    {
        __success = false;
        return -1;
    }
    const double shift = compute_shift(__alpha,theta,zeta);
    //std::cerr << "shift=" << shift << std::endl;
    const double xi    = zeta - shift;
    //std::cerr << "xi   =" << xi << std::endl;
    return xi - __xi;
}

//#include "yocto/math/fcn/zfind.hpp"

double Bridge:: find_theta_by_xi(const double alpha, const double xi, bool *global_flat, double &shift)
{
    __alpha = alpha;
    __xi    = xi;
    if(global_flat) *global_flat = false;
    shift  = 0;

    Function &F = zprof_zeta;


    double zeta_a = xi;
    double jump_a = F(zeta_a);

    if(!__success)
    {
        std::cerr << "couldn't initialize zeta search" << std::endl;
        return -1;
    }

    static const double shrink = 0.8;
    double zeta_b = zeta_a;
    double jump_b = jump_a;
    while(true)
    {
        zeta_b = zeta_a - shrink * jump_a;
        jump_b = F(zeta_b);
        if(!__success)
        {
            std::cerr << "looking for zeta failure level-1" << std::endl;
            return -1;
        }
        if(jump_a * jump_b<=0) break;
        zeta_a = zeta_b;
        jump_a = jump_b;
    }

    //std::cerr << "Bracketed..." << std::endl;
    //std::cerr << "zeta  =" << zeta_a << " " << zeta_b << std::endl;
    //std::cerr << "jump  =" << jump_a << " " << jump_b << std::endl;

    const double dz = max_of<double>( numeric<double>::sqrt_ftol, odeint.eps );
    while(Fabs(zeta_a-zeta_b)>dz)
    {
        const double zeta_m = 0.5 * (zeta_a+zeta_b);
        const double jump_m = F(zeta_m);
        if(!__success)
        {
            std::cerr << "looking for theta failure level-2" << std::endl;
            return -1;
        }
        if(jump_m*jump_a>0)
        {
            zeta_a = zeta_m;
            jump_a = jump_m;
        }
        else
        {
            zeta_b = zeta_m;
            jump_b = jump_m;
        }
    }
    bool         is_flat = false;
    const double zeta    = 0.5 * (zeta_a+zeta_b);
    double       theta   = find_theta(alpha,zeta,&is_flat);
    if(global_flat) *global_flat = is_flat;
    shift                = compute_shift(alpha,theta,zeta);

    return theta;
}
