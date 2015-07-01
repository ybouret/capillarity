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
        // get the surface
        const double Si = data.S[i];
        if(Si<0||Si>lens->max_surface)
        {
            throw exception("Invalid Data: maximum surface exceeded!");
        }

        // find the corresponding polar angle
        zfn.target = Si;
        triplet<double> aa = {0,0,lens->max_alpha};
        triplet<double> ss = {zfn.call(aa.a),0,zfn.call(aa.c)};
        if(ss.a*ss.b>0)
            throw exception("Unexpected Failure in height computation!");
        data.alpha[i] = Rad2Deg(solve.run(zfn.call,aa,ss));

        // then find the corresponding contact angle
        data.theta[i] = FindTheta(data.h[i]+h_shift+h_speed*i, data.alpha[i]);
        {
            ios::acstream fp(savename);
            fp("%g %g %g %g\n", data.h[i], data.S[i], data.alpha[i], data.theta[i] );
        }
        (std::cerr << ".").flush();
    }
    std::cerr << std::endl;

}

#include "yocto/code/utils.hpp"
double Bridge:: Extract( DataFile &data )
{
    zfind<double>     solve(  ATOL );
    zfunction<double> zfn( lens->surface, 0 );


    double Smax = 0;
    for(size_t i=1;i<=data.N;++i)
    {
        // get the surface
        const double Si = data.S[i];
        Smax = max_of(Si,Smax);
        if(Si<0||Si>lens->max_surface)
        {
            throw exception("Invalid Data: maximum surface exceeded!");
        }

        // find the corresponding polar angle
        zfn.target = Si;
        triplet<double> aa = {0,0,lens->max_alpha};
        triplet<double> ss = {zfn.call(aa.a),0,zfn.call(aa.c)};
        if(ss.a*ss.b>0)
            throw exception("Unexpected Failure in height computation!");
        data.alpha[i] = Rad2Deg(solve.run(zfn.call,aa,ss));

        // then find the corresponding contact angle
        data.theta[i] = FindTheta(data.h[i]+h_shift+h_speed*i, data.alpha[i]);
    }

    double err = 0;
    for(size_t i=2;i<data.N;++i)
    {
        const double tmp = (data.theta[i+1]-data.theta[i-1])*(data.S[i+1]-data.S[i-1])/Smax;
        err += tmp*tmp;
    }
    std::cerr << "\terr=" << err << std::endl;
    return err;

}