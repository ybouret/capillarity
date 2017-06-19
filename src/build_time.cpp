#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/types.hpp"
#include "yocto/fs/vfs.hpp"

using namespace yocto;
using namespace math;


static double resolution = 1e9;

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
    vector<unit_t> ih(n);
    for(size_t i=1;i<=n;++i)
    {
        ih[i] = floor( h[i] * resolution + 0.5 );
        //std::cerr << h[i] << " => " << ih[i] << std::endl;
    }

    vector<double> t(n);
    for(size_t i=2;i<=n;++i)
    {
        const unit_t h_prev = ih[i-1];
        const unit_t h_curr = ih[i];
        if(h_curr<h_prev)
        {
            t[i] = t[i-1] + (h_prev-h_curr);
        }
        else
        {
            t[i] = t[i-1] + (h_curr-h_prev);
        }
    }
    for(size_t i=1;i<=n;++i)
    {
        t[i] /= resolution;
    }

    string       root = vfs::get_base_name(filename);
    const string rext = vfs::get_extension(root);
    vfs::change_extension(root,"time.");
    root += rext;
    std::cerr << "saving into " << root << std::endl;
    ios::wcstream fp(root);
    for(size_t i=1;i<=n;++i)
    {
        fp("%.15g %.15g %.15g\n",A[i],h[i],t[i]);
    }


}
YOCTO_PROGRAM_END()
