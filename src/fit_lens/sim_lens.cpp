#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/types.hpp"
#include "yocto/code/rand.hpp"

using namespace yocto;

YOCTO_PROGRAM_START()
{
    alea_init();
    const size_t N = 50 + alea_lt(100);
    const double Xmin = -1;
    const double Xmax = 1;
    ios::ocstream fp("lens.dat",false);
    for(size_t i=0;i<=N;++i)
    {
        const double x = Xmin + (i * (Xmax-Xmin) )/N;// + (alea<double>()-0.5) * 0.1;
        const double y = 0.8*((1.0 - cos(1.4*x) + (alea<double>()-0.5) * 0.02));
        fp("%g %g\n",x,y);
    }
}
YOCTO_PROGRAM_END()

