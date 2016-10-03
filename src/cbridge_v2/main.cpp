#include "bridge.hpp"

#include "yocto/program.hpp"
#include "yocto/lua/lua-state.hpp"

#include "yocto/ios/ocstream.hpp"

YOCTO_PROGRAM_START()
{
    Lua::State VM;
    lua_State *L = VM();

    Lua::Config::DoString(L,"ftol=1e-7;");
    Lua::Config::DoString(L,"angle_control=1;");
    Lua::Config::DoString(L,"shift_control=0.1;");
    Lua::Config::DoString(L,"R0=80;");
    Lua::Config::DoString(L,"lambda=2.72;");
    Lua::Config::DoString(L,"theta=120;");
    Lua::Config::DoString(L,"Xi=0.0;");
    Lua::Config::DoString(L,"resolution=1e-3");
    
    for(int i=1;i<argc;++i)
    {
        Lua::Config::DoString(L,argv[i]);
    }

    Bridge B(L);


    const double zeta      = Lua::Config::Get<lua_Number>(L,"zeta");
    const double alpha_deg = Lua::Config::Get<lua_Number>(L,"alpha");
    const double alpha     = Deg2Rad(alpha_deg);
    std::cerr << "alpha=" << alpha << ", zeta=" << zeta << std::endl;
    B.find_theta_v2(alpha,zeta);


#if 0
    const double Xi        = Lua::Config::Get<lua_Number>(L,"Xi");
    std::cerr << "Xi=" << Xi << std::endl;

    const double alpha_deg = Lua::Config::Get<lua_Number>(L,"alpha");
    const double alpha     = Deg2Rad(alpha_deg);
    std::cerr << "alpha=" << alpha_deg << " (" << alpha << " rad)" << std::endl;

    {
        ios::wcstream fp("profile.dat");
        ios::wcstream rp("results.dat");

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
#endif


#if 0
    const double theta_deg = Lua::Config::Get<lua_Number>(L,"theta");
    const double theta     = Deg2Rad(theta_deg);
    std::cerr << "theta=" << theta_deg << std::endl;

#if 0
    double       alpha_min = 0;
    double       zeta_max  = 0;
    const double Xi_max    = B.find_Xi_max(theta,alpha_min,zeta_max);
    std::cerr << "Xi_max   =" << Xi_max    << std::endl;
    std::cerr << "alpha_min=" << Rad2Deg(alpha_min) << std::endl;
    std::cerr << "zeta_max =" << zeta_max << ", " << zeta_max * B.R0 << " mm" << std::endl;
    {
        ios::wcstream fp("profmax.dat");
        (void)B.profile(alpha_min, theta, Xi_max, &fp, true);
        const double usink = B.compute_user_sink(alpha_min,theta,Xi_max);
        fp << "\n";
        fp("-1 %g 0\n",usink);
        fp("0  %g 0\n",usink);
        fp("%g %g 0\n",B.radii.front(),usink);
    }
    B.SaveLens("lensmax.dat", Xi_max);
#endif
    

    {
        ios::wcstream fp("profile.dat");
        ios::wcstream rp("results.dat");
        for(double alpha_deg=0.1;alpha_deg<=90;alpha_deg+=0.1)
        {
            const double ans = B.profile(Deg2Rad(alpha_deg), theta, Xi, &fp);
            fp << "\n";
            //rp("%g %g %g\n",alpha_deg,ans,sin(B.param[BRIDGE_A]));
            //rp("%g %g %g\n",alpha_deg,ans,B.reduced_curv(B.param));
            rp("%.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",alpha_deg,ans,B.param[BRIDGE_U],B.param[BRIDGE_V],B.param[BRIDGE_A],B.reduced_rate(B.param),B.reduced_curv(B.param));
        }
    }

    B.SaveLens("lens.dat", Xi);

    bool is_flat = false;
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
        std::cerr << "\t\tNO BRIDGE!" << std::endl;
    }
#endif

}
YOCTO_PROGRAM_END()

