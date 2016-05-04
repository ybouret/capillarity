#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/sequence/lw-array.hpp"

YOCTO_PROGRAM_START()
{
#if 0
    double __coefs[] = {63.037518,-2.5192661,10.626694,4.5368727};
    lw_array<double> coefs(__coefs,sizeof(__coefs)/sizeof(__coefs[0]));

    Lens lens(0.58205858,coefs,drvs);
#endif
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
                fp("%g %g\n", xx, yy);
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

        Bridge bridge;

        bridge.compute_profile(*lens, Deg2Rad(40.0), Deg2Rad(100.0), 0.0);

    }
}
YOCTO_PROGRAM_END()

