#include "vlens.hpp"
#include "bridge.hpp"
#include "datafile.hpp"
#include "yocto/program.hpp"
#include "yocto/ptr/auto.hpp"

#include "yocto/threading/simd.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/code/rand.hpp"
#include "yocto/string/conv.hpp"

using namespace threading;

YOCTO_PROGRAM_START()
{
    alea_init();

    if(argc<=2)
    {
        throw exception("usage: %s lens.file capillary_length [files]", program);
    }
    const string   lens_filename = argv[1];
    const double   clength       = strconv::to<double>(argv[2],"capillary_length");
    auto_ptr<Lens> lens( VLens::ReadFrom(lens_filename) );
    auto_ptr<Lens> sphr( new SphericalLens(lens->rho(0)));
    SIMD   simd_shape;
    SIMD   simd_round;
    Bridge bridge(*lens,clength);
    simd_shape.create<Bridge,const Lens&,double>(*lens,clength);
    simd_round.create<Bridge,const Lens&,double>(*sphr,clength);

    const size_t nt = simd_shape.size; assert(simd_shape.size==simd_round.size);
    context::kernel kProc( cfunctor(Bridge::CallProcess) );
    
    for(int i=3;i<argc;++i)
    {
        const string filename = argv[i];
        DataFile     data(filename,1,0);

        // send args to all bridges
        for(size_t j=0;j<nt;++j)
        {
            simd_shape[j].as<Bridge>().args = &data;
            simd_round[j].as<Bridge>().args = &data;
        }

        simd_round(kProc);


        const string saveroot = vfs::get_base_name(filename);

        {
            const string savename = saveroot+"_round.dat";
            data.Save(savename);
        }

        simd_shape(kProc);
        {
            const string savename = saveroot+"_shape.dat";
            data.Save(savename);
        }

    }

#if 0
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
#endif


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

    //bridge.Tests();
    context::kernel kTest( cfunctor(Bridge::CallTests) );
    simd(kTest);
#endif


}
YOCTO_PROGRAM_END()

