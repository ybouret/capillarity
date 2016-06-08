#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    Bridge B(0.001,  // angular search
             1e-5,   // integrator ftol
             0.1,    // max angular turn in degrees
             1.0/100 // max speed rate in degrees
             );
    
    int    iarg = 0;
    
    
    if(argc>++iarg)
    {
        B.mu   = strconv::to<double>(argv[iarg],"mu");
    }
    std::cerr << "mu=" << B.mu << std::endl;
    
    
    
    //__________________________________________________________________________
    //
    // find zeta_max
    //__________________________________________________________________________
    const size_t N=101;
    vector<double> zeta(N);
    vector<double> surf(N);
    vector<double> sfit(N);

    const string heur_name = vformat("heur%g.dat",B.mu);
    ios::ocstream::overwrite(heur_name);
    
    for(int t_deg = 40; t_deg <= 170; t_deg+=10)
    {
        std::cerr << "theta=" << t_deg << std::endl; std::cerr.flush();
        const double theta     = Deg2Rad( double(t_deg) );
        const double zeta_max  = B.compute_zeta_max( theta);
        std::cerr << "\tzeta_max=" << zeta_max << std::endl;
        const double zeta_min  = -0.1;
        const double zeta_half = zeta_max/2;
        for(size_t i=1;i<=N;++i)
        {
            const double ztmp = zeta_min + (i-1)*(zeta_half-zeta_min)/double(N-1);
            zeta[i] = ztmp;
            const double alpha = B.find_alpha(theta,ztmp);
            surf[i]            = numeric<double>::pi * Square( sin(alpha) );
            ios::acstream hp(heur_name);
            hp("%.15g %.15g\n", zeta[i], surf[i]);
            std::cerr << "\tzeta=" << zeta[i] << " => " << surf[i] << std::endl;
        }

        {
            ios::acstream hp(heur_name);
            hp << "\n";
        }


    }
    
}
YOCTO_PROGRAM_END()
