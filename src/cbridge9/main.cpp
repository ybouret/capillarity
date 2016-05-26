#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/sequence/lw-array.hpp"
#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    SharedDerivative drvs( new Derivative() );
    if(argc>1)
    {
        shared_ptr<Lens> lens( Lens::load(argv[1],drvs) );

        {
            ios::wcstream fp("lens.dat");
            for(double ad=-180;ad<=180;ad+=0.1)
            {
                const double alpha = Deg2Rad(double(ad));
                const double rho   = lens->R(alpha);
                const double xx    = rho * sin(alpha);
                const double yy    = lens->R0 - rho * cos(alpha);
                fp("%g %g\n", xx, yy+0);
            }
        }

        {
            ios::wcstream fp("omega.dat");
            for(double ad=0;ad<=180;ad+=0.01)
            {
                const double alpha = Deg2Rad(ad);
                const double omega = lens->omega(alpha);
                fp("%g %g %g\n", ad, Rad2Deg(omega), (*drvs)(lens->R,alpha,1e-4));
            }
        }

        Bridge bridge(0.01);
        bridge.capillary_length = 2.7;






    }
}
YOCTO_PROGRAM_END()

