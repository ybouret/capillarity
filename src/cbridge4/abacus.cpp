#include "bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"
    const double h_min = Lua::Config::Get<lua_Number>(L, "h_min");
    ios::ocstream::overwrite("abacus.dat");
    ios::ocstream::overwrite("h_max.dat");
    const double zeta_min = h_min/B.R0;
    const size_t N        = 20;
    bool         isFlat   = false;

    for(double theta_deg=100; theta_deg <= 170; theta_deg += 10 )
    {
        std::cerr << "THETA=" << theta_deg << std::endl;
        std::cerr << "\tlooking for zeta_max..." << std::endl;
        const double theta    = Deg2Rad(theta_deg);
        const double zeta_max = B.find_zeta_max(theta);
        std::cerr << "\tzeta_max=" << zeta_max << std::endl;
        ios::ocstream::echo("h_max.dat","%g %g\n", theta_deg, zeta_max * B.R0);

        const double amplitude = zeta_max - zeta_min;
        for(size_t i=0;i<=N;++i)
        {
            double zeta = zeta_min;
            if(i>=N) zeta=zeta_max; else zeta+=(i*amplitude)/double(N);
            const double alpha = B.find_alpha(theta,zeta,isFlat);
            const double h     = zeta * B.R0;
            const double A     = numeric<double>::pi * Square( B.R0 * sin(alpha) );

            ios::acstream fp("abacus.dat");
            fp("%.15g %.15g %.15g\n",h,A,Rad2Deg(alpha));
            std::cerr.flush();
            fprintf(stderr,"\tzeta=%10.4f    \r",zeta);
        }
        fprintf(stderr,"\n");
        fflush(stderr);
        
        ios::ocstream::echo("abacus.dat","\n");

    }

}
YOCTO_PROGRAM_END()

