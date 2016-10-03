#include "bridge.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"

#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
    Lua::State VM;
    lua_State *L = VM();

    Lua::Config::DoString(L,"ftol=1e-7;");
    Lua::Config::DoString(L,"angle_control=1;");
    Lua::Config::DoString(L,"shift_control=0.1;");
    Lua::Config::DoString(L,"R0=80;");
    Lua::Config::DoString(L,"lambda=2.72;");
    Lua::Config::DoString(L,"resolution=1e-3");

    for(int i=1;i<argc;++i)
    {
        Lua::Config::DoString(L,argv[i]);
    }

    //Bridge::curvature_coeff = 1;
    const double theta = Lua::Config::Get<lua_Number>(L,"theta");
    Bridge B(L);

    



}
YOCTO_PROGRAM_END()
