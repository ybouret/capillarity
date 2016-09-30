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
    Lua::Config::DoString(L,"R0=80;");
    Lua::Config::DoString(L,"lambda=2.72;");
    Lua::Config::DoString(L,"theta=150;");
    Lua::Config::DoString(L,"Xi=0.0;");
    Lua::Config::DoString(L,"resolution=0.01");
    
    for(int i=1;i<argc;++i)
    {
        Lua::Config::DoString(L,argv[i]);
    }

    Bridge B(L);




    const double theta_deg = Lua::Config::Get<lua_Number>(L,"theta");
    const double theta     = Deg2Rad(theta_deg);
    std::cerr << "theta=" << theta_deg << std::endl;

    const double Xi        = Lua::Config::Get<lua_Number>(L,"Xi");
    double       alpha_min = 0;
    double       zeta_max  = 0;
    const double Xi_max    = B.find_Xi_max(theta,alpha_min,zeta_max);
    std::cerr << "Xi_max   =" << Xi_max    << std::endl;
    std::cerr << "alpha_min=" << Rad2Deg(alpha_min) << std::endl;
    std::cerr << "zeta_max =" << zeta_max << ", " << zeta_max * B.R0 << " mm" << std::endl;
    {
        ios::wcstream fp("profmax.dat");
        (void)B.profile(alpha_min, theta, Xi_max, &fp, true);
        const double urise = B.compute_user_rise(alpha_min,theta,Xi_max);
        fp << "\n";
        fp("-1 %g 0\n",urise);
        fp("0  %g 0\n",urise);
        fp("%g %g 0\n",B.radii.front(),urise);
    }
    B.SaveLens("lensmax.dat", Xi_max);

    

    {
        ios::wcstream fp("profile.dat");
        ios::wcstream rp("results.dat");
        for(double alpha_deg=1;alpha_deg<=90;alpha_deg+=0.2)
        {
            const double ans = B.profile(Deg2Rad(alpha_deg), theta, Xi, &fp);
            fp << "\n";
            rp("%g %g %g\n",alpha_deg,ans,B.param[BRIDGE_U]);
        }
    }

    B.SaveLens("lens.dat", Xi);

    const double alpha_opt = B.find_alpha(theta,Xi);
    if(alpha_opt>0)
    {
        double urise = 0;
        std::cerr << "alpha_opt=" << Rad2Deg(alpha_opt) << std::endl;
        {
            ios::wcstream fp("alpha_opt.dat");
            (void)B.profile(alpha_opt, theta, Xi, &fp, true);

            urise = B.compute_user_rise(alpha_opt,theta,Xi);
            std::cerr << "urise=" << urise << ",zeta=" << Xi-urise << std::endl;
            fp << "\n";
            for(double xx=-1;xx<=1;xx+=0.1)
            {
                fp("%g %g 0\n",xx,urise);
            }
        }
        if(false)
        {
            B.SaveLens("real_lens.dat",Xi-urise);
            ios::wcstream fp("real_alpha.dat");
            (void)B.profile(alpha_opt, theta, Xi, &fp, true,-urise);
        }
    }

}
YOCTO_PROGRAM_END()

