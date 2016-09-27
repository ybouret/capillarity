#include "bridge.hpp"

#include "yocto/program.hpp"
#include "yocto/lua/lua-state.hpp"

#include "yocto/ios/ocstream.hpp"


YOCTO_PROGRAM_START()
{
    Lua::State VM;
    lua_State *L = VM();

    Lua::Config::DoString(L,"ftol=1e-5;");
    Lua::Config::DoString(L,"angle_control=0.1;");
    Lua::Config::DoString(L,"shift_control=0.01;");


    Bridge B(L);

    B.compute_mu(82,2.7);
    
    ios::wcstream fp("profile.dat");
    for(int alpha_deg=10;alpha_deg<=40;alpha_deg+=1)
    {
        B.profile(Deg2Rad(double(alpha_deg)), Deg2Rad(100.0), 0.0, &fp);
        fp << "\n";
    }


}
YOCTO_PROGRAM_END()

