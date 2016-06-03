#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/container/matrix.hpp"

YOCTO_PROGRAM_START()
{
    Bridge B(0.01,1e-5);
    
    const int    t_max = 175;
    const int    t_min = 25;
    const int    d_t   = 10;
    const size_t n     = (t_max-t_min)/d_t+1;
    
   
    double       mu[] = { 1, 10, 20, 30, 40, 50 };
    const size_t m = sizeof(mu)/sizeof(mu[0]);
    
    
    const string zfile = "zmax.dat";
    const string rfile = "ratio.dat";
    ios::ocstream::overwrite(zfile);
    ios::ocstream::overwrite(rfile);

    matrix<double> zmax(m,n);
    
    for(size_t j=1;j<=n;++j)
    {
        const double t_deg = t_min + (j-1)*d_t;
        const double theta = Deg2Rad(t_deg);
        std::cerr << "theta=" << t_deg << std::endl;
        {
            ios::acstream zp(zfile);
            zp("%.15g", t_deg);
            ios::acstream rp(rfile);
            rp("%.15g", t_deg);
        }
        for(size_t i=1;i<=m;++i)
        {
            B.mu    = mu[i-1];
            std::cerr << "\tmu=" << B.mu << std::endl;
            const double zeta_max = B.compute_zeta_max(theta);
            zmax[i][j] = zeta_max;
            
            ios::acstream zp(zfile);
            zp(" %.15g", zeta_max);
            
            ios::acstream rp(rfile);
            rp(" %15.g",zeta_max/zmax[i][1]);
        }
        
        {
            ios::acstream zp(zfile);
            zp << "\n";
            ios::acstream rp(rfile);
            rp << "\n";
        }
       
    }
    
    
    
    
    
}
YOCTO_PROGRAM_END()
