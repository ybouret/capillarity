#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/types.hpp"
#include "yocto/code/rand.hpp"

using namespace yocto;

YOCTO_PROGRAM_START()
{
    alea_init();
    const size_t N    = 200 + alea_lt(200);
    const double Xmin = -1;
    const double Xmax =  1;
    ios::ocstream fp("lens.dat",false);
    const double xc = 0.037;
    for(size_t i=0;i<=N;++i)
    {
        const double x =  Xmin + (i * (Xmax-Xmin) )/N;// + (alea<double>()-0.5) * 0.1;
        const double y = 0.8*((1.0 - cos(1.4*x) + (alea<double>()-0.5) * 0.01));
        //const double y = -sqrt(9-x*x);
        fp("%g %g\n",x+xc,y);
    }
}
YOCTO_PROGRAM_END()

