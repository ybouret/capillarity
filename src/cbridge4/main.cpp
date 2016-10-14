#include "application.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#define DISCARD_BRIDGE
#include "main-core.cpp"



#if 0
    const string filename = "zrange.dat";
    ios::ocstream::overwrite(filename);
    ios::ocstream::echo(filename,"0 0 0 180 180\n");
    for(double alpha_deg=1;alpha_deg<=30;++alpha_deg)
    {
        std::cerr << std::endl << "ALPHA=" << alpha_deg << std::endl;
        const double  alpha = Deg2Rad(alpha_deg);
        range_t       tr;
        const range_t zr    = B.find_zeta_range(alpha,tr);
        std::cerr << "theta: " << Rad2Deg(tr.vmin) << " to " << Rad2Deg(tr.vmax) << std::endl;
        ios::ocstream::echo(filename,"%g %g %g %g %g\n", alpha_deg, zr.vmin, zr.vmax, Rad2Deg(tr.vmin), Rad2Deg(tr.vmax));
    }
#endif

    Application  app(L);
    const string filename = Lua::Config::Get<string>(L,"file");
    app.load(filename);

    app.process();
}
YOCTO_PROGRAM_END()
