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
    
    //__________________________________________________________________________
    //
    // find zeta_max
    //__________________________________________________________________________
    for(int t_deg = 40; t_deg <= 170; t_deg+=10)
    {
        std::cerr << "theta=" << t_deg << std::endl; std::cerr.flush();
        const double theta = Deg2Rad( double(t_deg) );
        const double zeta_max = B.compute_zeta_max( theta);
        std::cerr << "zeta_max=" << zeta_max << std::endl;
        {
            ios::acstream zp(zeta_file);
            zp("%.15g %.15g\n", double(t_deg), zeta_max);
        }
    }
    
}
YOCTO_PROGRAM_END()
