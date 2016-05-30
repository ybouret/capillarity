#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    Bridge B(0.01,1e-5);

    double alpha_deg = 20;
    double theta_deg = 100;
    double zeta      = 0.0;
    
    if(argc>1)
    {
        alpha_deg = strconv::to<double>(argv[1],"alpha");
    }

    if(argc>2)
    {
        theta_deg = strconv::to<double>(argv[2],"theta");
    }

    if(argc>3)
    {
        zeta = strconv::to<double>(argv[3],"zeta");
    }

    if(argc>4)
    {
         B.mu = strconv::to<double>(argv[4],"mu");
    }




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
        ios::wcstream ap("ans.dat");
        ios::wcstream fp("prof.dat");
        for(theta_deg=2;theta_deg<175;theta_deg+=2)
        {
            const double ans = B.profile( Deg2Rad(alpha_deg), Deg2Rad(theta_deg), zeta, &fp);
            ap("%g %g\n", theta_deg, ans);
            fp << "\n";
        }

    }
#else
    {
        ios::wcstream ap("ans.dat");
        ios::wcstream fp("prof.dat");
        for(alpha_deg=0.1;alpha_deg<90;alpha_deg+=0.1)
        {
            const double ans = B.profile( Deg2Rad(alpha_deg), Deg2Rad(theta_deg), zeta, &fp);
            ap("%g %g %g %g\n", alpha_deg, ans, B.param[BRIDGE_U], sin(B.param[BRIDGE_A]));
            fp << "\n";
        }

    }
#endif

    
}
YOCTO_PROGRAM_END()
