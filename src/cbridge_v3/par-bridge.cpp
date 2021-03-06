#include "par-bridge.hpp"
#include "yocto/math/io/data-set.hpp"

ParBridge:: ~ParBridge() throw()
{
    
}

#define __CORE 0x01
#define __SUBS 0x02

ParBridge:: ParBridge(lua_State *L) : Bridge(L), Crew(true),
zeta(),
surf(),
times(),
alpha(),
theta(),
theta2(),
rate( Lua::Config::Get<lua_Number>(L,"rate") ),
mgr()
{
    Crew &self = *this;
    for(size_t i=0;i<self.size;++i)
    {
        self[i].build<Bridge,lua_State *>(L);
    }

    mgr.enroll(zeta,__CORE);
    mgr.enroll(surf,__CORE);
    mgr.enroll(times,__CORE);

    mgr.enroll(alpha,__SUBS);
    mgr.enroll(theta,__SUBS);
    mgr.enroll(theta2,__SUBS);

    std::cerr << "rate=" << rate << std::endl;
}

void ParBridge:: load( const string &filename )
{
    mgr.free_all();
    {
        data_set<double> ds;
        ds.use(1,surf);
        ds.use(2,zeta);
        ds.use(3,times);
        ios::icstream fp(filename);
        ds.load(fp);
    }
    const double n = zeta.size();
    mgr.make_all(__SUBS,n);

    for(size_t i=n;i>0;--i)
    {
        zeta[i] /= R0;
        surf[i] /= A0;
        if(surf[i]>1||surf[i]<=0) throw exception("area invalid!");
        alpha[i] = Asin( Sqrt(surf[i] ) );
    }

}

void ParBridge:: find_theta()
{
    Crew &self = *this;
    threading::kernel k(this, & ParBridge::FindTheta );
    self.run(k);
    for(size_t i=0;i<self.size;++i)
    {
        self[i].as<Bridge>().change_curv(1);
    }
}

void ParBridge:: FindTheta( Context &ctx )
{
    size_t i      = 1;
    size_t length = zeta.size();
    ctx.split(i,length);
    Bridge &B     = ctx.as<Bridge>();
    for(;length>0;++i,--length)
    {
        theta[i]  = B.find_theta(alpha[i],zeta[i], NULL);
        theta2[i] = B.find_theta(alpha[i],zeta[i]+(rate*times[i])/R0,NULL);

#if 0
        //double shift = 0;
        Theta[i] = B.find_theta_by_xi(alpha[i],zeta[i], NULL, shift);

        B.change_curv(2);
        theta2[i] = B.find_theta(alpha[i],zeta[i], NULL);
        Theta2[i] = B.find_theta_by_xi(alpha[i],zeta[i], NULL, shift);
#endif

    }
}

