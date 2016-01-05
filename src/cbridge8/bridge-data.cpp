#include "bridge.hpp"
#include "datafile.hpp"
#include "yocto/math/fcn/zfind.hpp"

void Bridge::Process(const Context &ctx) throw()
{
    assert(args);
    DataFile &data = *(DataFile *)args;
    size_t i  = 1;
    size_t length = data.N;

    ctx.split(i,length);

#if 0
    {
        scoped_lock guard(ctx.access);
        std::cerr << "Process from " << win.start << " to " << win.final << std::endl;
    }
#endif

    zfind<double>     solve( ATOL );
    zfunction<double> zfn(lens->surface,0);
    for(;length>0;++i,--length)
    {
        const double Si = data.S[i];
        if(Si<=0||Si>lens->max_surf_value)
        {
            data.alpha[i] = -1;
            data.theta[i] = -1;
            continue;
        }
        zfn.target = Si;
        triplet<double> aa = { 0,0,lens->max_surf_angle };
        triplet<double> ss = { zfn.call(aa.a),0,zfn.call(lens->max_surf_angle) };
        if(ss.a*ss.c>0)
        {
            data.alpha[i] = -1;
            data.theta[i] = -1;
            continue;
        }
        const double alpha = solve.run(zfn.call,aa,ss);
        data.alpha[i] = alpha;
        data.theta[i] = FindTheta(data.h[i], alpha);
    }
}

void Bridge::CallProcess(Context &ctx) throw()
{
    //Bridge &bridge = ctx.as<Bridge>();
    //bridge.Process(ctx);
}