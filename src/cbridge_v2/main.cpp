#include "bridge.hpp"

#include "yocto/program.hpp"
#include "yocto/lua/lua-state.hpp"

#include "yocto/ios/ocstream.hpp"


YOCTO_PROGRAM_START()
{
    Lua::State VM;
    lua_State *L = VM();

    Lua::Config::DoString(L,"ftol=1e-5;");
    Lua::Config::DoString(L,"angle_control=0.01;");
    Lua::Config::DoString(L,"shift_control=0.01;");


    Bridge B(L);

    B.compute_mu(82,2.7);

    ios::wcstream rp("result.dat");
    ios::wcstream fp("profile.dat");
    for(int alpha_deg=5;alpha_deg<=20;alpha_deg+=1)
    {
        const double ans = B.profile(Deg2Rad(double(alpha_deg)), Deg2Rad(100.0), 0.0, &fp);
        fp << "\n";
        rp("%.15g %.15g\n", double(alpha_deg), ans);
    }


}
YOCTO_PROGRAM_END()

