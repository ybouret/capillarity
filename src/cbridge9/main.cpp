#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/sequence/lw-array.hpp"

YOCTO_PROGRAM_START()
{
    double __coefs[] = {63.037518,-2.5192661,10.626694,4.5368727};
    lw_array<double> coefs(__coefs,sizeof(__coefs)/sizeof(__coefs[0]));

    Derivative drvs;
    Lens lens(0.58205858,coefs,drvs);

    {
        ios::wcstream fp("lens.dat");
        for(double ad=0;ad<=180;ad+=0.1)
        {
            const double alpha = Deg2Rad(double(ad));
            const double rho   = lens.R(alpha);
            const double xx    = rho * sin(alpha);
            const double yy    = lens.R0 - rho * cos(alpha);
            fp("%g %g\n", xx, yy);
        }
    }

}
YOCTO_PROGRAM_END()

