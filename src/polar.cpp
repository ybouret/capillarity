#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"

#include "yocto/math/types.hpp"
#include "yocto/math/point2d.hpp"
#include "yocto/math/trigconv.hpp"

using namespace yocto;
using namespace math;

static inline double Func(double phi)
{
    const double arg = (2*phi/numeric<double>::pi);

    return 0.9+0.2*arg*arg;

}

YOCTO_PROGRAM_START()
{
    const double a = 10.0;
    const double b = 6.0;

    {
        ios::wcstream fp("ell.dat");
        ios::wcstream fp2("ell2.dat");
        const double phi = Deg2Rad(0.0);
        for(double theta=-1.6;theta<=1.6;theta+=0.01)
        {
            const double rho  = a*b/Hypotenuse(a*cos(theta),b*sin(theta));
            const double rho2 = a*b/Hypotenuse(a*cos(theta),b*sin(theta))*Func(theta);

            {
                const double X   =  rho * sin(theta);
                const double Y   = -rho * cos(theta);
                const double x   = cos(phi)*X-sin(phi)*Y;
                const double y   = sin(phi)*X+cos(phi)*Y;
                fp("%g %g %g %g\n", x, y, theta, rho );
            }

            {
                const double X   =  rho2 * sin(theta);
                const double Y   = -rho2 * cos(theta);
                const double x   = cos(phi)*X-sin(phi)*Y;
                const double y   = sin(phi)*X+cos(phi)*Y;
                fp2("%g %g %g %g\n", x, y, theta, rho2 );
            }
        }
    }
    
}
YOCTO_PROGRAM_END()
