#include "bridge.hpp"
#include "datafile.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/trigconv.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/math/opt/cgrad.hpp"

Optimizer::Optimizer(Bridge &B, DataFile &D):
bridge(B),
data(D),
field(this, & Optimizer::Compute )
{

}

Optimizer:: ~Optimizer() throw() {}

double Optimizer:: Compute( const array_t &p )
{
    switch(p.size())
    {
        case 2:
            bridge.h_speed = p[2];

        case 1:
            bridge.h_shift = p[1];
            break;

        default:
            throw exception("invalid #variable=%u in Optimizer::Compute", unsigned(p.size()));
    }
    return bridge.Extract(data);
}

bool Optimizer:: Callback(const array_t &p)
{
    std::cerr << "params=" << p << std::endl;
    return true;
}

void Optimizer:: Run( array_t &p )
{
    const size_t n = p.size();
    switch (n) {
        case 2:

        case 1:
            break;

        default:
            throw exception("invalid #variable=%u in Optimizer::Run", unsigned(n));

    }

    vector<double> dp(n,0.0);
    double &dh = dp[1];
    for(size_t i=1;i<data.N;++i)
    {
        dh += data.h[i+1]-data.h[i];
    }
    dh /= (data.N-1);
    dh /= 10;
    std::cerr << "dh=" << dh << std::endl;

    cgrad<double> cg;
    cgrad<double>::callback cb(this, & Optimizer::Callback);
    cg.run(field, p, dp, 1e-4,&cb);


}
