#include "application.hpp"

#define __CORE 0x01
#define __SUBS 0x02

Application:: Application( lua_State *L ) :
Bridge(L), threading::par_server(true),
h(),
A(),
t(),
zeta(),
alpha(),
theta()
{
    threading::executor &self = *this;
    for(size_t i=0;i<self.num_threads();++i)
    {
        self[i].build<Bridge,lua_State *>(L);
    }

    mgr.enroll(h,__CORE);
    mgr.enroll(A,__CORE);
    mgr.enroll(t,__CORE);

    mgr.enroll(zeta,__SUBS);
    mgr.enroll(alpha,__SUBS);
    mgr.enroll(theta,__SUBS);

}

Application:: ~Application() throw()
{
}

#include "yocto/math/io/data-set.hpp"

void Application:: load( const string &filename )
{
    mgr.release_all();
    // loading
    {
        data_set<double> ds;
        ds.use(1, A);
        ds.use(2, h);
        ds.use(3, t);
        ios::icstream fp(filename);
        ds.load(fp);
    }

    // precomputing
    const size_t n = h.size();
    mgr.make_all(__SUBS,n);
    for(size_t i=1;i<=n;++i)
    {
        zeta[i] = h[i]/R0;
        const double aa = A[i] / A0;
        if(aa<=0||aa>1)
            throw exception("Invalid Area=%g",aa);
        alpha[i] =  asin( sqrt(aa) );
        theta[i] = 0;
    }
}
