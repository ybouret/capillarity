#include "bridge.hpp"

#include "yocto/program.hpp"
#include "yocto/lua/lua-state.hpp"

#include "yocto/ios/ocstream.hpp"

YOCTO_PROGRAM_START()
{
#define DISCARD_XI 1
#include "find_core.cpp"

    const double alpha_deg = Lua::Config::Get<lua_Number>(L,"alpha");
    const double alpha     = Deg2Rad(alpha_deg);
    std::cerr << "alpha=" << alpha_deg << " (" << alpha << " rad)" << std::endl;

    const double zeta = Lua::Config::Get<lua_Number>(L,"zeta");
    double shift = 0;
    bool   is_flat = 0;
    const double theta_opt = B.find_theta_v2(alpha, zeta, shift, &is_flat);
    B.SaveLens("lens.dat",zeta);
    if(theta_opt>0)
    {
        const double Xi = zeta+shift;
        std::cerr << "theta_opt=" << Rad2Deg(theta_opt) << std::endl;
        std::cerr << "shift=" << shift << std::endl;
        if(!is_flat)
        {
            ios::wcstream fp("theta_opt.dat");
            (void)B.profile(alpha, theta_opt, Xi, &fp,false,-shift);
        }
        else
        {
            ios::wcstream fp("theta_opt.dat");
            B.compute_start(alpha, theta_opt, Xi);
            for(double xx=0;xx<=1.4;xx+=0.1)
            {
                fp("%g %g\n", xx+B.param[BRIDGE_U], B.param[BRIDGE_V]-shift);
            }
        }

    }
    else
    {
        ios::ocstream::overwrite("theta_opt.dat");
        std::cerr << "NO BRIDGE!!!!" << std::endl;
    }



}
YOCTO_PROGRAM_END()
