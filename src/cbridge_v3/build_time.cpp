#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"


using namespace yocto;
using namespace math;


YOCTO_PROGRAM_START()
{
    if(argc<=1)
    {
        throw exception("usage: %s datafile",program);
    }

    const string filename = argv[1];

    vector<double> h;
    vector<double> A;
    {
        data_set<double> ds;
        ds.use(1, A);
        ds.use(2, h);
        ios::icstream fp( filename );
        ds.load(fp);
    }
    const size_t n = h.size();
    std::cerr << "#data=" << n << std::endl;
}
YOCTO_PROGRAM_END()
