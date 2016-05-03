#include "bridge.hpp"
#include "datafile.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/trigconv.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/math/opt/cgrad.hpp"

Optimizer::Optimizer(Bridge &B, DataFile &D):
bridge(B),
data(D),
F(this, & Optimizer::Function),
G(this, & Optimizer::Gradient)
{

}

Optimizer:: ~Optimizer() throw() {}

double Optimizer:: Function( const array_t &p )
{
    switch(p.size())
    {
        case 2:
            //bridge.h_speed = p[2];

        case 1:
            //bridge.h_shift = p[1];
            break;

        default:
            throw exception("invalid #variable=%u in Optimizer::Compute", unsigned(p.size()));
    }
    return bridge.Extract(data);
}

void Optimizer:: Gradient(array_t &g, const array_t &p)
{
    assert(g.size()==p.size());
    const size_t n = p.size();


    vector<double> q(n,0);
    for(size_t i=1;i<=n;++i)
    {
        g[i] = 0;
        q[i] = p[i];
    }


    switch(n)
    {
        case 1:
        {
            const double f0 = Function(p);
            q[1] = p[1] + dh;
            const double f1 = Function(q);
            g[1] = (f1-f0)/dh;
        }
            break;

        default:
            throw exception("No Gradient for #variable=%u", unsigned(n));
    }
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

    dh = 0;
    for(size_t i=1;i<data.N;++i)
    {
        dh += Fabs(data.h[i+1]-data.h[i]);
    }
    dh /= (data.N-1);
    dh /= 100;
    std::cerr << "dh=" << dh << std::endl;

    cgrad<double>::callback cb(this, & Optimizer::Callback);
    if(cgrad<double>::optimize(F, G, p, 1e-4, &cb))
    {
        (void)F(p);
        bridge.Process(data, "output.dat");
    }
    else
    {
        throw exception("Cannot optimize...");
    }
    
}
