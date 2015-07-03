#include "vlens.hpp"
#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ptr/auto.hpp"

#include "yocto/threading/simd.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/code/rand.hpp"

using namespace threading;

YOCTO_PROGRAM_START()
{
    alea_init();

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

    SIMD   simd;
    Bridge bridge(*lens,2.7);
    simd.create<Bridge,const Lens&,double>(*lens,2.7);

#if 0
    ios::ocstream::overwrite("hmax.dat");
    for(int theta=5;theta<=175;theta+=5)
    {
        //bridge.ScanAlpha(0, theta);
        std::cerr << "theta=" << theta << std::endl;
        const double hmax = bridge.FindHmax(theta);
        ios::acstream fp("hmax.dat");
        fp("%d %g\n", theta, hmax);
    }
#endif
    
    //bridge.Tests();
    context::kernel kTest( cfunctor(Bridge::CallTests) );
    simd(kTest);

}
YOCTO_PROGRAM_END()

