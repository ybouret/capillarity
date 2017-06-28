#include "application.hpp"

#define __CORE 0x01
#define __SUBS 0x02

Application:: Application( Lua::State &L ) :
Bridge(L), threading::par_server(false),
h(),
A(),
tv(),
h_evap(),
h_corr(),
zeta(),
alpha(),
theta(),

coef_evap(0),
coef_push(0),
coef_pull(0)
{
    threading::executor &self = *this;
    for(size_t i=0;i<self.num_threads();++i)
    {
        self[i].build<Bridge,Lua::State&>(L);
    }

    mgr.enroll(h,__CORE);
    mgr.enroll(A,__CORE);

    mgr.enroll(zeta,__SUBS);
    mgr.enroll(alpha,__SUBS);
    mgr.enroll(theta,__SUBS);
    mgr.enroll(tv,__SUBS);
    mgr.enroll(h_evap,__SUBS);

}

Application:: ~Application() throw()
{
}


static const double h_scaling = 1e9;

void Application:: build_tv()
{
    const size_t n = h.size();
    std::cerr << "#data=" << n << std::endl;
    vector<unit_t> ih(n);
    for(size_t i=1;i<=n;++i)
    {
        ih[i] = floor( h[i] * h_scaling + 0.5 );
    }
    tv[1] = 0;
    for(size_t i=2;i<=n;++i)
    {
        const unit_t h_prev = ih[i-1];
        const unit_t h_curr = ih[i];
        if(h_curr<h_prev)
        {
            tv[i] = tv[i-1] + (h_prev-h_curr);
        }
        else
        {
            tv[i] = tv[i-1] + (h_curr-h_prev);
        }
    }
    for(size_t i=1;i<=n;++i)
    {
        tv[i] /= h_scaling;
    }

}

void Application:: correct_h()
{
    // evaporation
    const size_t n = h.size();
    for(size_t i=n;i>0;--i)
    {
        h_evap[i] = h[i] + coef_evap * tv[i];
    }

    // the pull/push
    for(size_t i=n;i>0;--i)
    {
        if(h[i]>h_evap[i])
        {
            // push
            h_corr[i] = h_evap[i] - coef_push * (h[i]-h_evap[i]);
        }
        else
        {
            // pull
            h_corr[i] = h_evap[i] + coef_push * (h_evap[i]-h[i]);
        }
    }
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
        ios::icstream fp(filename);
        ds.load(fp);
    }


    // precomputing
    const size_t n = h.size();
    mgr.make_all(__SUBS,n);

    build_tv();

    correct_h();

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
