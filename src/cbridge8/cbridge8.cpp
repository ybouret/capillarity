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
    SIMD::Kernel kProc( cfunctor(Bridge::CallProcess) );
    
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



}
YOCTO_PROGRAM_END()

