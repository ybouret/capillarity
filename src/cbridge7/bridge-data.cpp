#include "bridge.hpp"
#include "datafile.hpp"
#include "yocto/math/fcn/zfind.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/trigconv.hpp"

#include "yocto/ios/ocstream.hpp"

void Bridge:: Process(DataFile &data, const string &savename)
{
    {
        ios::wcstream fp(savename);
        fp << "#h S alpha theta\n";
    }
    
    zfind<double>     solve(  ATOL );
    zfunction<double> zfn( lens->surface, 0 );

    for(size_t i=1;i<=data.N;++i)
    {
        const double Si = data.S[i];
        if(Si<0||Si>lens->max_surface)
        {
            throw exception("Invalid Data: maximum surface exceeded!");
        }
        zfn.target = Si;
        triplet<double> aa = {0,0,lens->max_alpha};
        triplet<double> ss = {zfn.call(aa.a),0,zfn.call(aa.c)};
        if(ss.a*ss.b>0)
            throw exception("Unexpected Failure in height computation!");
        data.alpha[i] = Rad2Deg(solve.run(zfn.call,aa,ss));
        data.theta[i] = FindTheta(data.h[i], data.alpha[i]);
        {
            ios::acstream fp(savename);
            fp("%g %g %g %g\n", data.h[i], data.S[i], data.alpha[i], data.theta[i] );
        }
        (std::cerr << ".").flush();
    }
    std::cerr << std::endl;

}
