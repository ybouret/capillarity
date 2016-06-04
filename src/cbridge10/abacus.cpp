#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    Bridge B(0.01,1e-5);
    
    int    iarg = 0;
    
    
    if(argc>++iarg)
    {
        B.mu   = strconv::to<double>(argv[iarg],"mu");
    }
    std::cerr << "mu=" << B.mu << std::endl;
    
    
    
    const string zeta_file = vformat("zeta_max%g.dat",B.mu);
    ios::ocstream::overwrite(zeta_file);
    const string abac_file = vformat("abacus%g.dat",B.mu);
    ios::ocstream::overwrite(abac_file);
    //__________________________________________________________________________
    //
    // find zeta_max
    //__________________________________________________________________________
    const double d_zeta = 0.002;
    for(int t_deg = 40; t_deg <= 170; t_deg+=10)
    {
        std::cerr << "theta=" << t_deg << std::endl; std::cerr.flush();
        const double theta    = Deg2Rad( double(t_deg) );
        const double zeta_max = B.compute_zeta_max( theta);
        std::cerr << "zeta_max=" << zeta_max << std::endl;
        {
            ios::acstream zp(zeta_file);
            zp("%.15g %.15g\n", double(t_deg), zeta_max);
        }
        {
            ios::acstream fp(abac_file);
            fp("#theta=%g\n", double(t_deg));
        }
        const double zeta_min = -0.1;
        const double N = ceil((zeta_max-zeta_min)/d_zeta);
        std::cerr << "N=" << N << std::endl;
        for(size_t i=0;i<=N;++i)
        {
            const double zeta = (i<=0) ? zeta_min : ( (i>=N)  ? zeta_max : zeta_min + (i*(zeta_max-zeta_min)/double(N)) );
            const double alpha = B.find_alpha(theta,zeta);
            const double surf  = numeric<double>::pi * Square( sin(alpha) );
            std::cerr << "\tzeta=" << zeta << " => " << alpha << " rad => " << Rad2Deg(alpha) << " deg" << std::endl;
            ios::acstream fp(abac_file);
            fp("%.15g %.15g %.15g %.15g\n",zeta,surf,Rad2Deg(alpha),alpha);
        }
        {
            ios::acstream fp(abac_file);
            fp << "\n";
        }
    }
    
}
YOCTO_PROGRAM_END()
