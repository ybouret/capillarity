#include "application.hpp"

#define __CORE 0x01
#define __SUBS 0x02

Application:: Application( Lua::State &L ) :
Bridge(L), threading::par_server(false),
h(),
A(),
t(),
v(),
h_evap(),
h_corr(),
zeta(),
alpha(),
theta(),
main_rate(L.Get<lua_Number>("main_rate")),
coef_evap(L.Get<lua_Number>("coef_evap")),
coef_push(L.Get<lua_Number>("coef_push")),
coef_pull(L.Get<lua_Number>("coef_pull"))

{
    threading::executor &self = *this;
    for(size_t i=0;i<self.num_threads();++i)
    {
        self[i].build<Bridge,Lua::State&>(L);
    }

    mgr.enroll(h, __CORE);
    mgr.enroll(A, __CORE);

    mgr.enroll(zeta,   __SUBS);
    mgr.enroll(alpha,  __SUBS);
    mgr.enroll(theta,  __SUBS);
    mgr.enroll(t,      __SUBS);
    mgr.enroll(h_evap, __SUBS);
    mgr.enroll(h_corr, __SUBS);
    mgr.enroll(v,      __SUBS);
}

Application:: ~Application() throw()
{
}



#include "yocto/math/stat/int-hist.hpp"
#include "yocto/math/stat/descr.hpp"
#include "yocto/math/stat/dist.hpp"

#include "yocto/sort/remove-if.hpp"

static const double time_resolution = 1000;

static inline double get_xtime(const double tmx)
{
    return Floor( tmx * time_resolution + 0.5 );
}

void Application:: build_time()
{
    const size_t   n = h.size();

    vector<long>   dt(n-1);
    vector<long>   bins;
    vector<size_t> hist;

    //__________________________________________________________________________
    //
    // getting the dt histogram
    //__________________________________________________________________________
    for(size_t i=1;i<n;++i)
    {
        dt[i] = long( get_xtime( Fabs(h[i+1]-h[i])/main_rate ) );
    }

    i_histogram(bins,hist,dt);

    //__________________________________________________________________________
    //
    // taking the mode
    //__________________________________________________________________________
    const long dt_mode = histogram_mode(bins,hist);
    std::cerr << "dt_mode=" << dt_mode << std::endl;

    //__________________________________________________________________________
    //
    // checking for a lower speed, taking significative decreases
    //__________________________________________________________________________
    double lower_rate_factor = 1;
    {
        vector<double> dcoef(n-1,as_capacity);
        for(size_t i=1;i<n;++i)
        {
            const double dh = Fabs(h[i+1]-h[i]);
            if(dh>0)
            {
                double       local_rate = (dh/dt_mode) * time_resolution;
                if(local_rate<main_rate)
                {
                    const double d = main_rate/local_rate;
                    if(RInt(d)>1)
                    {
                        dcoef.push_back(d);
                    }
                }
            }
        }
        if(dcoef.size()>0)
        {
            std::cerr << "dcoef=" << dcoef << std::endl;
            double dmu=0,dsig=0;
            compute_average_and_stddev(dmu,dsig,dcoef);
            std::cerr << "dmu=" << dmu << ", dsig=" << dsig << std::endl;
            const double dlo = Floor(dmu);
            const double dup = Ceil(dmu);
            const double plo = gaussian<double>::pdf(dlo, dmu, dsig);
            const double pup = gaussian<double>::pdf(dup, dmu, dsig);
            std::cerr << "dlo=" << dlo << "@" << plo << std::endl;
            std::cerr << "dup=" << dup << "@" << pup << std::endl;
            if(plo>=pup)
            {
                lower_rate_factor = dlo;
            }
            else
            {
                lower_rate_factor = dup;
            }
        }
    }

    std::cerr << "lower_rate_factor=" << lower_rate_factor << std::endl;
    const double lower_rate = main_rate / lower_rate_factor;


    //__________________________________________________________________________
    //
    // reconstructing dt and speed
    //__________________________________________________________________________
    for(size_t i=1;i<n;++i)
    {
        const double dh = Fabs(h[i+1]-h[i]);
        v[i]  = main_rate;
        dt[i] = dt_mode;
        if(dh>0)
        {
            double local_rate = (dh/dt_mode) * time_resolution;
            if(local_rate<main_rate)
            {
                if( RInt(main_rate/local_rate) > 1 )
                {
                    v[i] = lower_rate;
                }
            }
            else
            {
                const double m  = local_rate/main_rate;
                const long   im = long(RInt(m));
                if(im>1)
                {
                    // increase time
                    dt[i] = im * dt_mode;
                }
            }
        }
    }


    //__________________________________________________________________________
    //
    // reconstructing time
    //__________________________________________________________________________
    const double h0  = h[1];
    const double dh0 = h[2] - h[1];
    const double v0  = (dh0 >= 0) ? v[1] : -v[1];
    std::cerr << "h0=" << h0 << std::endl;
    std::cerr << "v0=" << v0 << std::endl;

    t[1] = h0/v0;
    for(size_t i=2;i<=n;++i)
    {
        t[i] = t[i-1] + dt[i-1]/time_resolution;
    }
    v[n] = v[n-1];



    {
        ios::wcstream fp("thv.dat");
        fp << "#t h v\n";
        fp("0 0 %g\n", double(v[1]));
        for(size_t i=1;i<=n;++i)
        {
            fp("%g %g %g\n", double(t[i]), double(h[i]), double(v[i]) );
        }
    }



    exit(0);
}

void Application:: correct_h()
{

    // evaporation
    const size_t n = h.size();
    for(size_t i=n;i>0;--i)
    {
        h_evap[i] = h[i] + coef_evap * t[i];
        h_corr[i] = h_evap[i];
    }

    for(size_t i=n;i>0;--i)
    {
        if(h_evap[i]>0)
        {
            // pull !
            h_corr[i] = h_evap[i] + coef_pull * h_evap[i];
        }
        else
        {
            // push
            h_corr[i] = h_evap[i] + coef_push * h_evap[i];
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

    build_time();

    //build_tv();

    //correct_h();

    for(size_t i=1;i<=n;++i)
    {
        zeta[i] = h_corr[i]/R0;
        const double aa = A[i] / A0;
        if(aa<=0||aa>1)
            throw exception("Invalid Area=%g",aa);
        alpha[i] = asin( sqrt(aa) );
        theta[i] = 0;
    }
}
