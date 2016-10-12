#include "bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

    double       zeta      = Lua::Config::Get<lua_Number>(L, "zeta" );
    double       theta_deg = Lua::Config::Get<lua_Number>(L,"theta");
    double       theta     = Deg2Rad(theta_deg);
    std::cerr << "theta=" << theta_deg << std::endl;

    B.SaveLens("lens.dat", zeta);

    bool isFlat = false;
    
    B.find_alpha(theta,zeta,isFlat);


}
YOCTO_PROGRAM_END()
