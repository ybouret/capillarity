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

    {
        ios::wcstream fp("profile_theta.dat");
        ios::wcstream rp("results_theta.dat");

        for(double theta_deg=1;theta_deg<=179;theta_deg+=1)
        {
            const double ans = B.profile(alpha, Deg2Rad(theta_deg), Xi, &fp);
            rp("%g %g\n",theta_deg,ans);
            fp << "\n";
        }

    }
    B.SaveLens("lens.dat", Xi);

    bool   is_flat   = false;
    double theta_opt = B.find_theta(alpha, Xi, &is_flat );
    if(theta_opt>0)
    {
        std::cerr << "theta_opt=" << Rad2Deg(theta_opt) << std::endl;
        double usink = 0;
        if(!is_flat)
        {
            ios::wcstream fp("theta_opt.dat");
            (void)B.profile(alpha, theta_opt, Xi, &fp, true);

            usink = B.compute_user_sink(alpha,theta_opt,Xi);
            std::cerr << "usink=" << usink << ",zeta=" << Xi-usink << std::endl;
            fp << "\n";
            for(double xx=-1;xx<=1;xx+=0.1)
            {
                fp("%g %g 0\n",xx,usink);
            }
        }
        else
        {
            ios::wcstream fp("theta_opt.dat");
            B.compute_start(alpha, theta_opt, Xi);
            for(double xx=0;xx<=1.4;xx+=0.1)
            {
                fp("%g %g\n", xx+B.param[BRIDGE_U], B.param[BRIDGE_V]);
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
