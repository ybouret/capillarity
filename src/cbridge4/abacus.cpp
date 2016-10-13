#include "bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"
    //const double alpha_deg = Lua::Config::Get<lua_Number>(L, "alpha" );
    //const double alpha     = Deg2Rad(alpha_deg);

    ios::ocstream::overwrite("abacus_alpha.dat");

    {
        ios::acstream fp("abacus_alpha.dat");
        fp("0  0  0 180 180\n");
    }

    for(double alpha_deg=1;alpha_deg<=30; alpha_deg += 1)
    {
        std::cerr << "\tALPHA=" << alpha_deg << std::endl;
        const double  alpha = Deg2Rad(alpha_deg);
        const range_t zr    = B.find_zeta_range(alpha);
        ios::acstream fp("abacus_alpha.dat");
        fp("%.15g %.15g %.15g\n",alpha_deg,zr.vmin,zr.vmax);
    }

}
YOCTO_PROGRAM_END()

