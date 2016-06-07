#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    Bridge B(0.001,1e-5,0.1,1.0/100);

    double zeta = 0;
    int    iarg = 0;


    if(argc>++iarg)
    {
        B.mu   = strconv::to<double>(argv[iarg],"mu");
    }
    std::cerr << "mu=" << B.mu << std::endl;

#if 1
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
    std::cerr << "theta=" << theta_deg << std::endl;


    {
        ios::wcstream fp("lens.dat");
        for(double ad=-180;ad<=180;ad+=0.1)
        {
            const double a = Deg2Rad(ad);
            fp("%g %g\n",sin(a),zeta+(1.0-cos(a)));
        }
    }

#if 1
    {
        ios::wcstream pp("prof.dat");
        ios::wcstream ap("ans.dat");
        for(double a_deg=0.1;a_deg<=179;a_deg += 0.1)
        {
            const double alpha = Deg2Rad(a_deg);
            const double ans   = B.profile(alpha, theta, zeta, &pp);
            pp << "\n";
            ap("%g %g %g\n", a_deg, ans, sin(B.param[BRIDGE_A]) );
        }
    }
#endif

    const double alpha = B.find_alpha(theta,zeta);
    std::cerr << "alpha=" << Rad2Deg(alpha) << std::endl;
    
    
}
YOCTO_PROGRAM_END()
