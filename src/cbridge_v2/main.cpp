#include "bridge.hpp"

#include "yocto/program.hpp"
#include "yocto/lua/lua-state.hpp"

#include "yocto/ios/ocstream.hpp"


YOCTO_PROGRAM_START()
{
    Lua::State VM;
    lua_State *L = VM();

    Lua::Config::DoString(L,"ftol=1e-5;");
    Lua::Config::DoString(L,"angle_control=5;");
    Lua::Config::DoString(L,"shift_control=0.1;");
    Lua::Config::DoString(L,"R0=2;");
    Lua::Config::DoString(L,"lambda=1;");
    Lua::Config::DoString(L,"theta=150;");
    Lua::Config::DoString(L,"zeta=0.2;");

    for(int i=1;i<argc;++i)
    {
        Lua::Config::DoString(L,argv[i]);
    }

    Bridge B(L);

    //B.compute_mu(82,2.7);
    ios::wcstream rp("result.dat");
    ios::wcstream fp("profile.dat");

    const double theta_deg = Lua::Config::Get<lua_Number>(L,"theta");
    const double theta     = Deg2Rad(theta_deg);
    const double zeta      = Lua::Config::Get<lua_Number>(L,"zeta");

    B.SaveLens("lens.dat", zeta);

    for(double alpha_deg=0.1;alpha_deg<=30;alpha_deg+=0.1)
    {
        const double ans = B.profile(Deg2Rad(alpha_deg), theta, zeta, &fp);
        fp << "\n";
        //rp("%.15g %.15g %.15g\n", alpha_deg, ans, (B.param[BRIDGE_A]) );
        rp("%.15g %.15g %.15g\n", alpha_deg, ans, B.reduced_rate(B.param) );

    }


}
YOCTO_PROGRAM_END()

