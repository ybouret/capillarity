#include "cbridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

    const double theta_deg = Lua::Config::Get<lua_Number>(L,"theta");
    const double theta     = Deg2Rad(theta_deg);

    std::cerr << "theta=" << theta_deg << std::endl;

    const double zeta_max = B.find_zeta_max(theta);
    


}
YOCTO_PROGRAM_END()
