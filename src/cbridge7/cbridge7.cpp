#include "yocto/program.hpp"
#include "bridge.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/code/rand.hpp"
YOCTO_PROGRAM_START()
{
    alea_init();
    if(argc<=2)
        throw exception("need radius clength");
    const double R = strconv::to<double>(argv[1],"R");
    const double C = strconv::to<double>(argv[2],"C");

    const Lens::Pointer lens( new SphericalLens(R));
    Bridge bridge(lens, C);
    bridge.Tests();

}
YOCTO_PROGRAM_END()
