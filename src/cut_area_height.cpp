#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/string/conv.hpp"

using namespace yocto;
using namespace math;

YOCTO_PROGRAM_START()
{
    if(argc<=2)
    {
        throw exception("usage: %s area_height.dat h_cut",program);
    }

    const string filename = argv[1];
    const double h_cut    = strconv::to<double>(argv[2],"h_cut");
    vector<double> area;
    vector<double> height;

    {
        data_set<double> ds;
        ds.use(1,area);
        ds.use(2,height);
        ios::icstream fp(filename);
        ds.load(fp);
    }
    const size_t N = area.size();
    std::cerr << "#data=" << N << std::endl;

    const string enfoncement_id = "enfoncement.dat";
    const string tirage_id      = "tirage.dat";

    ios::wcstream enfoncement(enfoncement_id);
    ios::wcstream tirage(tirage_id);

    for(size_t i=1;i<=N;++i)
    {
        const double h = height[i];
        const double a = area[i];
        if(h<=h_cut)
        {
            // enfoncement
            enfoncement("%.15g %.15g\n",a,h);
        }
        else
        {
            // tirage
            tirage("%.15g %.15g\n",a,h);
        }
    }


}
YOCTO_PROGRAM_END()
