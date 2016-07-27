#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"

#include "yocto/math/types.hpp"
#include "yocto/math/point2d.hpp"
#include "yocto/math/trigconv.hpp"

using namespace yocto;
using namespace math;

YOCTO_PROGRAM_START()
{
    const double a = 10.0;
    const double b = 6.0;

    {
        ios::wcstream fp("ell.dat");
        const double phi = Deg2Rad(20.0);
        for(double theta=-1.6;theta<=1.6;theta+=0.01)
        {
            const double rho = a*b/Hypotenuse(a*cos(theta),b*sin(theta));
            const double X   =  rho * sin(theta);
            const double Y   = -rho * cos(theta);
            const double x   = cos(phi)*X-sin(phi)*Y;
            const double y   = sin(phi)*X+cos(phi)*Y;
            fp("%g %g %g %g\n", x, y, theta, rho );
        }
    }

}
YOCTO_PROGRAM_END()
