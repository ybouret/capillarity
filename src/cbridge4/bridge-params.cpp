#include "bridge.hpp"

////////////////////////////////////////////////////////////////////////////////
range_t Bridge:: find_zeta_range( const double alpha)
{
    bool         isFlat = false;
    const double zeta0  = CriticalZetaOfAlpha(alpha);
    range_t      zr;
    double      &zeta_min  = zr.vmin;
    double      &zeta_max  = zr.vmax;
    const double theta0    = find_theta(alpha, zeta0, isFlat);

    if(  theta0 <= 0 )
    {
        throw exception("No theta for critical zeta");
    }

    zeta_min  = zeta_max  =  zeta0;

    const double ztol = (numeric<double>::sqrt_ftol);

    {
        double zeta_top = zeta0;
        double zeta_bot = -2; // bad
        while( Fabs(zeta_top-zeta_bot) > ztol )
        {
            const double zeta_mid  = 0.5*(zeta_top+zeta_bot);
            const double theta_mid = find_theta(alpha,zeta_mid,isFlat);
            if(theta_mid<=0)
            {
                // invalid
                zeta_bot = zeta_mid;
            }
            else
            {
                // good
                zeta_top  = zeta_mid;
            }
        }
        zeta_min = zeta_top;
    }
   // std::cerr << "zeta_min  = " << zeta_min << std::endl;
    //std::cerr << "theta_max = " << theta_max << " => " << Rad2Deg(theta_max) << std::endl;

    {
        double zeta_bot = zeta0; // good
        double zeta_top = 0;     // is it bad ?
        while( find_theta(alpha, zeta_top,isFlat) > 0 )
        {
            zeta_top += 1;
        }

        while( Fabs(zeta_top-zeta_bot) > ztol )
        {
            const double zeta_mid  = 0.5*(zeta_top+zeta_bot);
            const double theta_mid = find_theta(alpha,zeta_mid,isFlat);
            if(theta_mid<=0)
            {
                // invalid
                zeta_top = zeta_mid;
            }
            else
            {
                // good
                zeta_bot  = zeta_mid;
            }
        }
        zeta_max = zeta_bot;
    }
    //std::cerr << "zeta_max  = " << zeta_max << std::endl;
    //std::cerr << "theta_min = " << theta_min<< " => " << Rad2Deg(theta_min) << std::endl;

    return zr;
}


