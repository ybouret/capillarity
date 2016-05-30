#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"

YOCTO_PROGRAM_START()
{
    Bridge B(1e-5);

    double alpha_deg = 20;
    double theta_deg = 100;
    double zeta      = 0.0;

    {
        ios::wcstream fp("lens.dat");
        for(double ad=-180;ad<=180;ad+=0.1)
        {
            const double a = Deg2Rad(ad);
            fp("%g %g\n",sin(a),zeta+(1.0-cos(a)));
        }
    }

    {
        ios::wcstream fp("prof.dat");
        B.profile( Deg2Rad(alpha_deg), Deg2Rad(theta_deg), zeta, &fp);
    }

    
}
YOCTO_PROGRAM_END()
