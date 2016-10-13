#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/sort/quick.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"
    const double hmin = Lua::Config::Get<lua_Number>(L, "hmin");
    const string filename = "abacus4.dat";

    ios::ocstream::overwrite(filename);
    const double zeta_min = hmin/B.R0;
    const size_t N        = 40;

    Vector h;
    Vector A;
    double alphas[2];
    bool   isFlat = false;
    
    for(double theta_deg=120; theta_deg <= 160; theta_deg += 10 )
    {
        std::cerr << "THETA=" << theta_deg << std::endl;
        std::cerr << "\tlooking for zeta_max..." << std::endl;
        const double theta    = Deg2Rad(theta_deg);
        const double zeta_max = B.find_zeta_max(theta);
        std::cerr << "zeta_max=" << zeta_max << std::endl;
        const double ampli = zeta_max - zeta_min;
        for(size_t i=0;i<=N;++i)
        {
            double   zeta = zeta_min;
            if(i>=N) zeta = zeta_max; else zeta += (i*ampli)/double(N);
            const size_t n = B.find_alpha(theta, zeta, alphas, isFlat);
            for(size_t j=0;j<n;++j)
            {
                const double alpha = alphas[j];
                h.push_back( zeta * B.R0 );
                A.push_back( numeric<double>::pi * Square( B.R0 * sin(alpha) ) );
            }
        }
        co_qsort(A,h);
        ios::acstream fp(filename);
        for(size_t i=1;i<=A.size();++i)
        {
            fp("%g %g\n", h[i], A[i]);
        }
        
    }

}
YOCTO_PROGRAM_END()

