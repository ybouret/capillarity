#include "par_bridge.hpp"

ParBridge:: ~ParBridge() throw()
{
}

ParBridge:: ParBridge(lua_State *L) :
Crew(true),
area(),
height(),
alpha(),
zeta(),
theta(),
Theta(),
rate( Lua::Config::Get<lua_Number>(L,"rate" ) )
{
    for(size_t i=0;i<size;++i)
    {
        Context &ctx = (*this)[i];
        ctx.build<Bridge,lua_State *,double>(L,1);
    }
}

#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"

void ParBridge:: extract_from(lua_State *L)
{
    const string filename = Lua::Config::Get<string>(L, "file");


    {
        data_set<double> ds;
        ds.use(1, area);
        ds.use(2, height);
        ios::icstream fp(filename);
        ds.load(fp);
    }

    const size_t n0 = area.size();

    alpha.free();
    alpha.ensure(n0);
    zeta.free();
    zeta.ensure(n0);

    Bridge &bridge = (*this)[0].as<Bridge>();

    const double a0 = numeric<double>::pi * Square( bridge.R0 );
    for(size_t i=1;i<=n0;++i)
    {
        const double a = area[i];
        if(a>a0)
        {
            throw exception("area is too big!");
        }
        alpha.push_back( Asin(Sqrt(a/a0) ) );
        zeta.push_back( height[i]/bridge.R0 );
    }

    const size_t n = alpha.size();
    std::cerr << "computed #alpha=" << n << std::endl;
    theta.make(n);
    Theta.make(n);

    std::cerr << "Launching Extraction..." << std::endl;
    threading::kernel K(this, & ParBridge:: run );
    (*this)(K);

    {
        ios::wcstream fp("output.dat");
        for(size_t i=1;i<=alpha.size();++i)
        {
            fp("%g %g %g %g %g\n", zeta[i]*bridge.R0, a0 * Square( sin(alpha[i]) ), Rad2Deg(alpha[i]), Rad2Deg(theta[i]), Rad2Deg(Theta[i]) );
        }
    }


}

void ParBridge:: run(Context &ctx)
{
    Bridge &B     = ctx.as<Bridge>();
    size_t i      = 1;
    size_t length = alpha.size();
    ctx.split(i,length);
    double shift = 0;
    for(;length>0;++i,--length)
    {
        B.change_curv(1);
        const double zz = zeta[i]+(rate*i)/B.R0;
        theta[i] = B.find_theta(alpha[i],zz);
        B.change_curv(1);
        Theta[i] = B.find_theta_v2(alpha[i],zz,shift);
        //Theta[i] = B.find_theta(alpha[i],zz);
    }
}
