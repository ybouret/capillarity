#include "yocto/program.hpp"
#include "bridge.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/code/rand.hpp"
#include "datafile.hpp"
#include "yocto/fs/local-fs.hpp"

YOCTO_PROGRAM_START()
{
    alea_init();
    if(argc<=2)
        throw exception("need radius clength [datafiles]");
    const double R = strconv::to<double>(argv[1],"R");
    const double C = strconv::to<double>(argv[2],"C");

    Lens::Pointer lens( new SphericalLens(R));
    Bridge        bridge(lens, C);
    for(int i=3;i<argc;++i)
    {
        const string filename = argv[i];
        std::cerr << "|_Loading " << filename << std::endl;
        DataFile     data(filename,1,0);

#if 1
        const string b_name   = vfs::get_base_name(filename);
        const string savename = b_name + ".dat";
        std::cerr << "|_Saving to: " << savename << std::endl;
        bridge.Process(data,savename);
#endif

#if 0
        Optimizer Opt(bridge,data);
        vector<double> p(1,0.0);
        Opt.Run(p);
#endif
        
    }

}
YOCTO_PROGRAM_END()
