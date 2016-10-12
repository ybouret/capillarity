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
        range_t       tr;
        const range_t zr    = B.find_zeta_range(alpha,tr);
        ios::acstream fp("abacus_alpha.dat");
        fp("%g %g %g %g %g\n",alpha_deg,zr.vmin,zr.vmax,Rad2Deg(tr.vmin),Rad2Deg(tr.vmax));
    }

}
YOCTO_PROGRAM_END()

