#include "bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"
    const double hmin = Lua::Config::Get<lua_Number>(L, "hmin");
    const string filename = "abacus4.dat";

    ios::ocstream::overwrite(filename);
    const double zeta_min = hmin/B.R0;
    const size_t N        = 20;

    Vector h;
    Vector A;
    for(double theta_deg=100; theta_deg <= 170; theta_deg += 10 )
    {
        std::cerr << "THETA=" << theta_deg << std::endl;
        std::cerr << "\tlooking for zeta_max..." << std::endl;
        const double theta    = Deg2Rad(theta_deg);
        const double zeta_max = B.find_zeta_max(theta);
        std::cerr << "zeta_max=" << zeta_max << std::endl;
        ios::ocstream::echo(filename, "%g %g\n", theta_deg, zeta_max);

    }

}
YOCTO_PROGRAM_END()

