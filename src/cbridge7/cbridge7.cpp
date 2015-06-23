#include "yocto/program.hpp"
#include "bridge.hpp"
#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    const Lens::Pointer lens( new SphericalLens(70) );
    Bridge bridge(lens, 70);

    for(int i=1;i<argc;++i)
    {
        const double theta = strconv::to<double>(argv[i],"theta");
        std::cerr << "alpha(h=0,theta=" << theta << ")=" << bridge.FindAlpha(0,theta) << std::endl;

    }

}
YOCTO_PROGRAM_END()
