#include "bridge.hpp"

#include "yocto/program.hpp"
#include "yocto/lua/lua-state.hpp"

#include "yocto/ios/ocstream.hpp"

YOCTO_PROGRAM_START()
{
#include "find_core.cpp"

    const double alpha_deg = Lua::Config::Get<lua_Number>(L,"alpha");
    const double alpha     = Deg2Rad(alpha_deg);
    std::cerr << "alpha=" << alpha_deg << " (" << alpha << " rad)" << std::endl;

    const double zeta = Lua::Config::Get<lua_Number>(L,"zeta");
    double shift = 0;
    bool   is_flat = 0;
    B.find_theta_v2(alpha, zeta, shift, &is_flat);
    


}
YOCTO_PROGRAM_END()
