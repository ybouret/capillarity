#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/sort/quick.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"
    const double hmin = Lua::Config::Get<lua_Number>(L, "hmin");
    const string filename  = "abacus4.dat";

    ios::ocstream::overwrite(filename);

    const double zeta_min = hmin/B.R0;
    const size_t N        = 200;

    Vector h;
    Vector A;
    double alphas[2];
    bool   isFlat = false;
    
    for(double theta_deg=110; theta_deg <= 175; theta_deg += 5 )
    {
        h.free();
        A.free();
        std::cerr << "THETA=" << theta_deg << std::endl;
        std::cerr << "\tlooking for zeta_max..." << std::endl;
        const double theta    = Deg2Rad(theta_deg);
        const double zeta_max = B.find_zeta_max(theta);
        const double ampli = zeta_max - zeta_min;
        std::cerr << "\tzeta_max=" << zeta_max << std::endl;

        std::cerr.flush();

        for(size_t i=0;i<=N;++i)
        {
            fprintf(stderr,"\t\t %10.2f%%   \r", (100.0*double(i))/double(N));
            double   zeta = zeta_min;
            if(i>=N) zeta = zeta_max; else zeta += (i*ampli)/double(N);
            const size_t n = B.find_alpha(theta, zeta, alphas, isFlat);
            for(size_t j=0;j<n;++j)
            {
                h.push_back( zeta * B.R0 );
                A.push_back( numeric<double>::pi * Square( B.R0 * sin(alphas[j]) ));
            }
        }
        fprintf(stderr,"\n");
        fflush(stderr);

        co_qsort(A,h);
        ios::acstream fp(filename);
        for(size_t i=1;i<=A.size();++i)
        {
            fp("%.15g %.15g\n",h[i],A[i]);
        }
        fp << "\n";

        //ios::ocstream::echo(filename, "\n");
        //ios::ocstream::echo(filename2,"\n");
    }

}
YOCTO_PROGRAM_END()

