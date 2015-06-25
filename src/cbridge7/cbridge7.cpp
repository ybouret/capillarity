#include "yocto/program.hpp"
#include "bridge.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/code/rand.hpp"
#include "datafile.hpp"

YOCTO_PROGRAM_START()
{
    alea_init();
    if(argc<=2)
        throw exception("need radius clength [datafiles]");
    const double R = strconv::to<double>(argv[1],"R");
    const double C = strconv::to<double>(argv[2],"C");

    Lens::Pointer lens( new SphericalLens(R));
    Bridge        bridge(lens, C);
    //bridge.Tests();
    for(int i=3;i<argc;++i)
    {
        const string filename = argv[i];
        std::cerr << "|_Loading " << filename << std::endl;
        DataFile     data(filename);
        
    }

}
YOCTO_PROGRAM_END()
