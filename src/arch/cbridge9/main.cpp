#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/sequence/lw-array.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/io/data-set.hpp"

YOCTO_PROGRAM_START()
{

    if(argc<=2)
    {
        throw exception("usage: %s lens_file.prm capillary_lenght [FILES]",program);
    }



    SharedDerivative drvs( new Derivative() );
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
        ios::wcstream sp("surface.dat");
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

    for(int i=3;i<argc;++i)
    {
        const string input_name = argv[i];
        std::cerr << "processing " << input_name << std::endl;
        vector<double> h;
        vector<double> A;
        {
            data_set<double> ds;
            ds.use(1, h);
            ds.use(2, A);
            ios::icstream fp(input_name);
            ds.load(fp);
        }
        const size_t n = h.size();
        std::cerr << "loaded #data=" << n << std::endl;
        string output_name = vfs::get_base_name(input_name);
        vfs::change_extension(output_name, "xtr.dat");
        std::cerr << "saving into " << output_name << std::endl;

        {
            {
                ios::wcstream fp(output_name);
                fp << "#h A alpha theta\n";
            }
            
            for(size_t i=1;i<=n;++i)
            {
                const double surf  = A[i];
                const double alpha = lens->find_alpha(A[i]);
                const double alpha_deg = Rad2Deg(alpha);
                std::cerr << "A=" << surf << " => alpha=" << alpha_deg << std::endl;
                const double theta     = bridge.FindTheta(*lens, alpha, h[i]);
                if(theta<0)
                {
                    throw exception("couldn't find theta for h=%g", h[i]);
                }

                const double theta_deg = Rad2Deg(theta);
                std::cerr << "\ttheta=" << theta_deg << std::endl;
                ios::acstream fp(output_name);
                fp("%g %g %g %g\n", h[i], surf, alpha_deg, theta_deg);
            }
        }
        
    }
    
}
YOCTO_PROGRAM_END()

