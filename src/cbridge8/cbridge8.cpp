#include "vlens.hpp"
#include "yocto/program.hpp"
#include "yocto/ptr/auto.hpp"


YOCTO_PROGRAM_START()
{
    auto_ptr<Lens> lens( new SphericalLens(70) );
    lens->initialize();
    std::cerr << "max angle=" << lens->max_surf_angle << std::endl;
    std::cerr << "max value=" << lens->max_surf_value << std::endl;

    if(argc>1)
    {
        const string filename = argv[1];
        lens.reset( VLens::ReadFrom(filename) );
        lens->initialize();
        std::cerr << "max angle=" << lens->max_surf_angle << std::endl;
        std::cerr << "max value=" << lens->max_surf_value << std::endl;
    }
}
YOCTO_PROGRAM_END()

