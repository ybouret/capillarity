#include "datafile.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/math/fcn/zfind.hpp"

using namespace math;

DataFile:: ~DataFile() throw()
{
}

DataFile:: DataFile( const string &filename) :
h(),
S(),
N(0)
{
    {
        data_set<double> ds;
        ds.use(1, h);
        ds.use(2, S);
        ios::icstream fp(filename);
        ds.load(fp);
    }
    (size_t &)N = h.size();
    alpha.make(N,0.0);
    theta.make(N,0.0);

#if 0
    zfind<double>     solve( Bridge::ATOL );
    zfunction<double> zfn( lens->surface, 0 );

    for(size_t i=1;i<=N;++i)
    {
        const double Si = S[i];
        if(Si<0||Si>lens->max_surface)
        {
            throw exception("Invalid Data: maximum surface exceeded!");
        }
        zfn.target = Si;
        triplet<double> aa = {0,0,lens->max_alpha};
        triplet<double> ss = {zfn.call(aa.a),0,zfn.call(aa.c)};
        if(ss.a*ss.b>0)
            throw exception("Unexpected Failure in height computation!");
    }
#endif

}
