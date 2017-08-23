#include "yocto/program.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/fit/glsf-spec.hpp"

using namespace yocto;
using namespace math;

static const double rescale = 1e-3; // microns to millimeters

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
        height[i] -= H0;
        height[i] *= rescale;
        depth[i]  *= rescale;
        std::cerr << depth[i] << " " << height[i] << std::endl;
    }

    vector<double> hfit(N);

    GLS<double>::Samples samples(1);
    GLS<double>::Sample &sample = samples.append(depth, height, hfit);

    vector<double> aorg(3);
    vector<double> aerr(aorg.size());
    vector<bool>   used(aorg.size(),true);

    if(!_GLS::Polynomial<double>::Start(sample,aorg))
    {
        throw exception("unable to initialize fit!!!!");
    }

    samples.prepare(aorg.size());

    _GLS::Polynomial<double> poly;
    GLS<double>::Function    fn(poly);
    if( ! samples.fit_with(fn, aorg, used, aerr) )
    {
        throw exception("unable to fully fit...");
    }
    GLS<double>::display(std::cerr, aorg, aerr);

    _GLS::Polynomial<double>::GnuPlot(std::cerr,aorg) << std::endl;
    {
        ios::wcstream fp("parabole.txt");
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g %g\n", depth[i], height[i], hfit[i]);
        }
    }


}
YOCTO_PROGRAM_END()
