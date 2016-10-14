#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/sort/quick.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"
    const double hmin = Lua::Config::Get<lua_Number>(L, "hmin");
    const string filename  = "abacus4.dat";
    const string filename2 = "abacus4b.dat";

    ios::ocstream::overwrite(filename);
    ios::ocstream::overwrite(filename2);

    const double zeta_min = hmin/B.R0;
    const size_t N        = 100;

    //Vector h;
    //Vector A;
    double alphas[2];
    bool   isFlat = false;
    
    for(double theta_deg=110; theta_deg <= 170; theta_deg += 10 )
    {
        std::cerr << "THETA=" << theta_deg << std::endl;
        std::cerr << "\tlooking for zeta_max..." << std::endl;
        const double theta    = Deg2Rad(theta_deg);
        const double zeta_max = B.find_zeta_max(theta);
        const double ampli = zeta_max - zeta_min;
        std::cerr << "\tzeta_max=" << zeta_max << std::endl;
        
        for(size_t i=0;i<=N;++i)
        {
            double   zeta = zeta_min;
            if(i>=N) zeta = zeta_max; else zeta += (i*ampli)/double(N);
            const size_t n = B.find_alpha(theta, zeta, alphas, isFlat);
            const double h = zeta*B.R0;

            switch(n)
            {
                case 0: break;
                case 1: ios::ocstream::echo(filename, "%g %g\n", h, numeric<double>::pi * Square( B.R0 * sin(alphas[0]) )); break;
                case 2:
                    ios::ocstream::echo(filename,  "%g %g\n", h, numeric<double>::pi * Square( B.R0 * sin(alphas[1]) ));
                    if(zeta>0)
                    {
                        ios::ocstream::echo(filename2, "%g %g\n", h, numeric<double>::pi * Square( B.R0 * sin(alphas[0]) ));
                    }
                    break;

            }
        }
        ios::ocstream::echo(filename, "\n");
        ios::ocstream::echo(filename2,"\n");
    }

}
YOCTO_PROGRAM_END()

