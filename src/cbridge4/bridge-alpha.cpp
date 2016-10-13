#include "bridge.hpp"

double Bridge:: find_alpha(const double theta, const double zeta, bool &isFlat)
{
    assert(theta>0);
    assert(theta<numeric<double>::pi);

    if(false)
    {
        ios::wcstream fp("alpha_try.dat");
        for(double alpha_deg = 1; alpha_deg <= 90; alpha_deg += 1 )
        {
            (std::cerr << ".").flush();
            const double alpha = Deg2Rad(alpha_deg);
            const double angle = find_theta(alpha, zeta, isFlat);
            fp("%g %g %g\n", alpha_deg, Rad2Deg(angle), Rad2Deg(theta) );
        }
        std::cerr << std::endl;
    }

    const double zeta0   = CriticalZetaOfTheta(theta);
    const double alpha0  = numeric<double>::pi-theta;
    const double theta0  = find_theta(alpha0,zeta0,isFlat);
    std::cerr << "zeta0 =" << zeta0  << std::endl;
    std::cerr << "theta0=" << Rad2Deg(theta0) << std::endl;

    return 0;
}
