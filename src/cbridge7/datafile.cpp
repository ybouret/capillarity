#include "datafile.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"

DataFile:: ~DataFile() throw()
{
}

DataFile:: DataFile( const string &filename, Lens::Pointer &lens ) :
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

    for(size_t i=1;i<=N;++i)
    {
        
    }

}
