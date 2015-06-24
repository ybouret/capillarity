#include "yocto/program.hpp"
#include "bridge.hpp"
#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    if(argc<=3)
        throw exception("need radius clength theta");
    const double R = strconv::to<double>(argv[1],"R");
    const double C = strconv::to<double>(argv[2],"C");
    
    const Lens::Pointer lens( new SphericalLens(R));
    Bridge              bridge(lens, C);

}
YOCTO_PROGRAM_END()
