#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/sequence/lw-array.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/math/round.hpp"

YOCTO_PROGRAM_START()
{
    if(argc<=2)
    {
        throw exception("usage: %s lens_file.prm capillary_length [theta...]",program);
    }

    const string lens_name = argv[1];

    SharedDerivative drvs( new Derivative() );
    shared_ptr<Lens> lens( Lens::load(lens_name,drvs) );

    const string base_name = vfs::get_base_name(lens_name);

    {
        string output_lens = base_name;
        vfs::change_extension(output_lens, "lens.dat");
        std::cerr << "saving into " << output_lens << std::endl;
        ios::wcstream fp(output_lens);
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
        string output_omega = base_name;
        vfs::change_extension(output_omega, "omega.dat");
        ios::wcstream fp(output_omega);

        string output_surf = base_name;
        vfs::change_extension(output_surf,"surf.dat");
        ios::wcstream sp(output_surf);

        for(double ad=0;ad<=180;ad+=0.01)
        {
            const double alpha = Deg2Rad(ad);
            const double omega = lens->omega(alpha);
            fp("%g %g %g\n", ad, Rad2Deg(omega), (*drvs)(lens->R,alpha,1e-4));
            const double rr     = lens->R(alpha)*sin(alpha);
            sp("%g %g\n", ad, numeric<double>::pi * rr*rr);
        }
    }


    Bridge bridge(0.01);
    bridge.capillary_length = strconv::to<double>(argv[2],"capillary_length");

    string output_abaccus  = base_name;
    const string out_id = vformat("abaccus%g.dat",bridge.capillary_length);
    vfs::change_extension(output_abaccus, out_id);
    ios::ocstream::overwrite(output_abaccus);


    
    for(int i=3;i<argc;++i)
    {
        const double theta_deg = strconv::to<double>(argv[i],"theta");
        const double theta     = Deg2Rad(theta_deg);
        double       dh        = log_round_floor(bridge.capillary_length/50);
        int          count     = 0;
        ios::acstream fp(output_abaccus);

        while(true)
        {

            const double h     = count * dh;
            const double alpha = bridge.FindAlpha(*lens, theta, h);
            if(alpha<0)
            {
                break;
            }
            const double alpha_deg = Rad2Deg(alpha);
            const double surf      = numeric<double>::pi * Square( lens->R(alpha) * sin(alpha) );
            fp("%g %g %g\n", h, surf, alpha_deg);
            ++count;
        }
        fp << "\n";
    }
    
    
    
}
YOCTO_PROGRAM_END()
