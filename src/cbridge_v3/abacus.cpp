#include "cbridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

    const double hmin = Lua::Config::Get<lua_Number>(L,"hmin");

    const size_t N    = 50;
    const string filename = "abacus.dat";
    ios::ocstream::overwrite(filename);

    for(int theta_deg = 60; theta_deg <= 160; theta_deg += 20 )
    {
        const double theta = Deg2Rad(double(theta_deg));
        std::cerr << std::endl;
        std::cerr << "theta=" << theta_deg << std::endl;

        const double zeta_max = B.find_zeta_max(theta);
        std::cerr << "\tzeta_max=" << zeta_max << std::endl;

        const double zeta_min = hmin/B.R0;
        for(size_t i=0;i<=N;++i)
        {
            const double zeta = zeta_min + ((zeta_max-zeta_min)*i)/double(N);
            const double alpha = B.find_alpha(theta, zeta, NULL);
            const double h     = zeta * B.R0;
            const double A     = numeric<double>::pi * Square( B.R0 * sin(alpha) );
            ios::acstream fp(filename);
            fp("%g %g\n",h,A);
        }
        {
            ios::acstream fp(filename); fp << "\n";
        }

    }

}
YOCTO_PROGRAM_END()
