#include "bridge.hpp"

#include "yocto/program.hpp"
#include "yocto/lua/lua-state.hpp"

#include "yocto/ios/ocstream.hpp"

YOCTO_PROGRAM_START()
{
#include "find_core.cpp"
    const double theta_deg = Lua::Config::Get<lua_Number>(L,"theta");
    const double theta     = Deg2Rad(theta_deg);
    std::cerr << "theta=" << theta_deg << std::endl;

    // make some profiles
    {
        ios::wcstream fp("profile_alpha.dat");
        ios::wcstream rp("results_alph.dat");
        for(double alpha_deg=0.1;alpha_deg<=90;alpha_deg+=0.1)
        {
            const double ans = B.profile(Deg2Rad(alpha_deg), theta, Xi, &fp);
            fp << "\n";
            rp("%.15g %.15g\n",alpha_deg,ans);
        }
    }

    B.SaveLens("lens.dat", Xi);

    bool         is_flat   = false;
    const double alpha_opt = B.find_alpha(theta,Xi,&is_flat);
    if(alpha_opt>0)
    {
        double usink = 0;
        std::cerr << "alpha_opt=" << Rad2Deg(alpha_opt) << std::endl;
        if(!is_flat)
        {
            ios::wcstream fp("alpha_opt.dat");
            (void)B.profile(alpha_opt, theta, Xi, &fp, true);

            usink = B.compute_user_sink(alpha_opt,theta,Xi);
            std::cerr << "usink=" << usink << ",zeta=" << Xi-usink << std::endl;
            fp << "\n";
            for(double xx=-1;xx<=1;xx+=0.1)
            {
                fp("%g %g 0\n",xx,usink);
            }
        }
        else
        {
            ios::wcstream fp("alpha_opt.dat");
            B.compute_start(alpha_opt, theta, Xi);
            for(double xx=0;xx<=1.4;xx+=0.1)
            {
                fp("%g %g\n", xx+B.param[BRIDGE_U], B.param[BRIDGE_V]);
            }
        }

        if(false)
        {
            B.SaveLens("real_lens.dat",Xi-usink);
            ios::wcstream fp("real_alpha.dat");
            (void)B.profile(alpha_opt, theta, Xi, &fp, true,-usink);
        }
    }
    else
    {
        ios::ocstream::overwrite("alpha_opt.dat");
        std::cerr << "\t\tNO BRIDGE!" << std::endl;
    }

    
}
YOCTO_PROGRAM_END()
