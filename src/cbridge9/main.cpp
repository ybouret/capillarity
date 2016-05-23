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
        double alpha_deg = 20;
        double theta_deg = 110;
        double h         = 0;
        if(argc>2) alpha_deg = strconv::to<double>(argv[2],"alpha");
        if(argc>3) theta_deg = strconv::to<double>(argv[3],"theta");
        if(argc>4) h         = strconv::to<double>(argv[4],"height");


        {
            ios::wcstream fp("lens.dat");
            for(double ad=-180;ad<=180;ad+=0.1)
            {
                const double alpha = Deg2Rad(double(ad));
                const double rho   = lens->R(alpha);
                const double xx    = rho * sin(alpha);
                const double yy    = lens->R0 - rho * cos(alpha);
                fp("%g %g\n", xx, yy+h);
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
        bridge.capillary_length = 2.7;

        const double alpha_hor = lens->find_horizontal_alpha( Deg2Rad(theta_deg) );
        std::cerr << "alpha_hor(" << theta_deg << ")=" << Rad2Deg(alpha_hor) << std::endl;
        {
            ios::wcstream pp("prof_alpha_hor.dat");
            const double z_alpha_hor = bridge.compute_profile(*lens, alpha_hor, Deg2Rad(theta_deg), h, &pp);
            std::cerr << "z_alpha_hor=" << z_alpha_hor << std::endl;
        }

        return 0;

        ios::wcstream pp("prof.dat");
        ios::wcstream pa("ans.dat");
#if 0
        for(alpha_deg=0.5;alpha_deg<=90;alpha_deg+=0.5)
        {
            const double ans = bridge.compute_profile(*lens, Deg2Rad(alpha_deg), Deg2Rad(theta_deg), h, &pp);
            pp << "\n";
            pa("%g %g\n", alpha_deg, ans);
        }
#else
        for(theta_deg=10;theta_deg<=170;theta_deg+=1)
        {
            const double ans = bridge.compute_profile(*lens, Deg2Rad(alpha_deg), Deg2Rad(theta_deg), h, &pp);
            pp << "\n";
            pa("%g %g\n", theta_deg, ans);
        }
#endif



    }
}
YOCTO_PROGRAM_END()

