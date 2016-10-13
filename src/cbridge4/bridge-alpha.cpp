#include "bridge.hpp"
#include "yocto/math/opt/bracket.hpp"
#include "yocto/math/opt/minimize.hpp"

double Bridge:: find_alpha(const double theta, const double zeta, bool &isFlat)
{
    assert(theta>0);
    assert(theta<numeric<double>::pi);

    assert(zeta>=-2);

    const double half_res      = resolution*0.5;
    Triplet      critical_zeta = { CriticalZetaOfTheta(theta-half_res), CriticalZetaOfTheta(theta), CriticalZetaOfTheta(theta+half_res) };
    critical_zeta.sort();
    const double zeta_lower = critical_zeta.a;
    const double zeta_upper = critical_zeta.c;

    if(zeta>=zeta_lower && zeta<=zeta_upper )
    {
        isFlat = true;
        return numeric<double>::pi-theta;
    }


    // setup
    __theta     = theta;
    __zeta      = zeta;
    Function &F = prof_alpha;
    isFlat      = false;

    // find possibilities
    const double alpha_lo = resolution;
    const double value_lo = F(alpha_lo);

    const double alpha_up = numeric<double>::pi-resolution;
    const double value_up = F(alpha_up);

    Triplet Alpha = { alpha_lo, 0, alpha_up };
    Triplet Value = { value_lo, 0, value_up };

    bracket<double>::inside(F, Alpha, Value);
    optimize1D<double>::run(F, Alpha, Value, resolution);

    std::cerr << "min=" << Value.b << std::endl;
    if(Value.b>0)
    {
        // no possible bridge
        return 0;
    }

    const double lower_alpha = (value_lo<=0) ? alpha_lo : find_lower(alpha_lo, Alpha.b,F,resolution);
    const double upper_alpha = (value_up<=0) ? alpha_up : find_upper(Alpha.b,alpha_up,F,resolution);

    std::cerr << "Possible alpha between " << Rad2Deg(lower_alpha) << " and " << Rad2Deg(upper_alpha) << std::endl;



#if 0
    double zeta_ini   = CriticalZetaOfTheta(theta);
    double alpha_ini  = numeric<double>::pi-theta;
    double theta_ini  = theta;
    std::cerr << "zeta_ini =" << zeta_ini << std::endl;
    std::cerr << "alpha_ini=" << Rad2Deg(alpha_ini) << std::endl;
    std::cerr << "theta_ini=" << Rad2Deg(theta_ini) << std::endl;
    int count=0;
    for(;;)
    {
        if(++count>3) break;
        std::cerr << std::endl << "count=" << count << std::endl;

        const double zeta_step  = zeta  - zeta_ini;
        const double alpha_step = theta - theta_ini;
        std::cerr << "zeta_step=" << zeta_step << ", alpha_step=" << Rad2Deg(alpha_step) << std::endl;
        double       scan       = 1.0;
    FORWARD:
        std::cerr << "scan=" << scan << std::endl;
        double zeta_new  = zeta_ini  + scan * zeta_step;
        double alpha_new = alpha_ini + scan * alpha_step;
        if(alpha_new<=0||alpha_new>=numeric<double>::pi)
        {
            if( Fabs(alpha_new-alpha_ini) <= 0 )
            {
                std::cerr << "stuck in alpha" << std::endl;
                return 0;
            }
            scan *= 0.5;
            goto FORWARD;
        }

        double theta_new = find_theta(alpha_new, zeta_new, isFlat);
        if( theta_new <= 0 )
        {
            if( Fabs(zeta_new-zeta_ini) <= 0 && Fabs(alpha_new-alpha_ini) <= 0 )
            {
                std::cerr << "stuck in step" << std::endl;
                return 0;
            }
            scan *= 0.5;
            goto FORWARD;
        }

        std::cerr << "zeta_new =" << zeta_new << std::endl;
        std::cerr << "alpha_new=" << Rad2Deg(alpha_new) << std::endl;
        std::cerr << "theta_new=" << Rad2Deg(theta_new) << std::endl;

        alpha_ini = alpha_new;
        theta_ini = theta_new;

    }
#endif
    
    
    return 0;
}
