#include "yocto/program.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/string/conv.hpp"

using namespace yocto;
using namespace math;

YOCTO_PROGRAM_START()
{
    if(argc<=1)
        throw exception("usage: %s datafile", program);

    size_t skip = 0;
    if(argc>2) skip = strconv::to<size_t>(argv[2],"skip");

    vector<double> depth;
    vector<double> height;
    {
        ios::icstream fp(argv[1]);
        data_set<double> ds;
        ds.use(1,depth);
        ds.use(2,height);
        ds.load(fp);
    }
    assert( depth.size() == height.size() );

    if(depth.size() <= 0 )
    {
        throw exception("no data");
    }

    const double H0 = height[1];
    std::cerr << "H0=" << H0 << std::endl;

    if(skip>=depth.size())
    {
        throw exception("skipping too many data");
    }

    while(skip>0)
    {
        depth.pop_front();
        height.pop_front();
        --skip;
    }

    const size_t N = depth.size();
    std::cerr << "N=" << N << std::endl;
    for(size_t i=1;i<=N;++i)
    {
        std::cerr << depth[i] << " " << height[i] << std::endl;
    }

}
YOCTO_PROGRAM_END()
