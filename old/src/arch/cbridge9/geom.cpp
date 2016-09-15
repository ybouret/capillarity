#include "yocto/program.hpp"
#include "bridge.hpp"

#include "yocto/string/conv.hpp"
#include "yocto/ios/ocstream.hpp"

YOCTO_PROGRAM_START()
{
    if(argc<=2)
    {
        throw exception("usage: %s R capillary_length",program);
    }

    const double R = strconv::to<double>(argv[1],"R");

    SharedDerivative drvs( new Derivative() );
    shared_ptr<Lens> lens( Lens::sphere(R,drvs) );

    Bridge bridge(0.01);
    bridge.capillary_length = strconv::to<double>(argv[2],"capillary_length");

    const string hmax_name = vformat("hmax_R=%g_L=%g.dat", R, bridge.capillary_length);
    const string prof_name = vformat("geom_R=%g_L=%g.dat", R, bridge.capillary_length);

    ios::ocstream::overwrite(hmax_name);
    ios::ocstream::overwrite(prof_name);

    for(double theta_deg = 5; theta_deg <= 175; theta_deg += 5 )
    {
        std::cerr << "theta=" << theta_deg << std::endl;
        const double theta = Deg2Rad( double(theta_deg) );
        const double hmax  = bridge.FindHMAX(*lens, theta);
        {
            ios::acstream fp(hmax_name);
            fp("%g %g %g\n",theta_deg,hmax,hmax/bridge.capillary_length);
        }

        {
            const size_t N = 100;
            for(size_t i=0;i<=N;++i)
            {
                const double h     = (i*hmax)/N;
                const double alpha = bridge.FindAlpha(*lens, theta, h);
                if(alpha<0)
                {
                    throw exception("unexpected alpha value...");
                }
                ios::acstream fp(prof_name);
                const double surf  = numeric<double>::pi * Square(lens->R(alpha) * sin(alpha));
                const double surf0 = numeric<double>::pi * Square(lens->R0);
                fp("%g %g %g %g\n", h, surf, surf/surf0, Rad2Deg(alpha));
            }
            {
                ios::acstream fp(prof_name);
                fp << "\n";
            }
        }
    }
    
    
    
}
YOCTO_PROGRAM_END()
