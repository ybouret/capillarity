#include "bridge.hpp"
#include "yocto/threading/window.hpp"
#include "yocto/math/fcn/zfind.hpp"
#include "yocto/exception.hpp"
#include "datafile.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/math/trigconv.hpp"

void Bridge:: ExtractMT(const threading::context &ctx,
                        DataFile                 &data) throw()
{
    lastResult = 0;

    zfind<double>     solve(  ATOL );
    zfunction<double> zfn( lens->surface, 0 );

    const threading::window W(ctx,data.N,1);
    for(size_t i=W.start;i<=W.final;++i)
    {
        // get the surface
        const double Si = data.S[i];
        if(Si<0||Si>lens->max_surface)
        {
            lastResult = 1;
            return;
            //throw exception("Invalid Data: maximum surface exceeded!");
        }

        // find the corresponding polar angle
        zfn.target = Si;
        triplet<double> aa = {0,0,lens->max_alpha};
        triplet<double> ss = {zfn.call(aa.a),0,zfn.call(aa.c)};
        if(ss.a*ss.b>0)
        {
            lastResult = 2;
            return;
            //throw exception("Unexpected Failure in height computation!");
        }

        data.alpha[i] = Rad2Deg(solve.run(zfn.call,aa,ss));

        // then find the corresponding contact angle
        data.theta[i] = FindTheta(data.h[i]+data.t[i]*ev_rate, data.alpha[i]);
    }


}