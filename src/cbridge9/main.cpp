#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"

YOCTO_PROGRAM_START()
{

    const double beta = 0.1;
    const double rho_beta = 0.9;
    const double rho_beta_prime  = 0.1;
    const double R  = 1;

    LensExtend ll(beta,rho_beta,rho_beta_prime,R);

    {
        ios::wcstream fp("rho.dat");
        for(double alpha=beta;alpha<=numeric<double>::pi;alpha+=0.01)
        {
            const double rho = ll.compute(alpha);
            const double x   = rho*sin(alpha);
            const double y   = R-rho*cos(alpha);
            fp("%g %g %g %g\n",x,y,alpha,rho);
        }
    }
    
}
YOCTO_PROGRAM_END()

