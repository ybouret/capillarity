#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    Bridge B(0.01,1e-5);

    double zeta = 0;
    int    iarg = 0;


    if(argc>++iarg)
    {
        B.mu   = strconv::to<double>(argv[iarg],"mu");
    }
    std::cerr << "mu=" << B.mu << std::endl;

#if 0
    if(argc>++iarg)
    {
        zeta = strconv::to<double>(argv[iarg],"zeta");
    }
    std::cerr << "zeta=" << zeta << std::endl;
#endif


    double theta_deg = 100;
    if(argc>++iarg)
    {
        theta_deg = strconv::to<double>(argv[iarg],"theta");
    }
    const double theta = Deg2Rad(theta_deg);
    std::cerr << "theta=" << theta << std::endl;

#if 0
    {
        ios::wcstream fp("lens.dat");
        for(double ad=-180;ad<=180;ad+=0.1)
        {
            const double a = Deg2Rad(ad);
            fp("%g %g\n",sin(a),zeta+(1.0-cos(a)));
        }
    }
#endif
    

    (void) B.find_alpha(theta,zeta);
    
    
    
}
YOCTO_PROGRAM_END()
