#include "datafile.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/math/fcn/zfind.hpp"

using namespace math;

DataFile:: ~DataFile() throw()
{
}

DataFile:: DataFile(const string &filename,
                    const double speed,
                    const double start) :
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
    t.make(N,0.0);
    alpha.make(N,0.0);
    theta.make(N,0.0);

    if(h[N]>=h[1])
    {
        // increasing
        for(size_t i=1;i<=N;++i)
        {
            t[i] = start + h[i]/speed;
        }
    }
    else
    {
        // decreasing
        for(size_t i=1;i<=N;++i)
        {
            t[i] = start + (h[N]-h[i])/speed;
        }
    }
}

#if 0
#include "yocto/ios/ocstream.hpp"

void DataFile:: save(const string &filename) const
{
    ios::wcstream fp(filename);
    fp << "#h S alpha theta\n";
    for(size_t i=1;i<=N;++i)
    {
        fp("%g %g %g %g\n", h[i], S[i], alpha[i], theta[i] );
    }
}
#endif
