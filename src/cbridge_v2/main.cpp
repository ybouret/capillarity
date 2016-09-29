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
    Lua::Config::DoString(L,"R0=2;");
    Lua::Config::DoString(L,"lambda=1;");
    Lua::Config::DoString(L,"theta=150;");
    Lua::Config::DoString(L,"zeta=0.0;");

    for(int i=1;i<argc;++i)
    {
        Lua::Config::DoString(L,argv[i]);
    }

    Bridge B(L);

    //B.compute_mu(82,2.7);
    ios::wcstream rp("result.dat");
    ios::wcstream fp("profile.dat");

    const double theta_deg = Lua::Config::Get<lua_Number>(L,"theta");
    const double theta     = Deg2Rad(theta_deg);
    const double zeta      = Lua::Config::Get<lua_Number>(L,"zeta");

    std::cerr << "theta=" << theta_deg << std::endl;
    std::cerr << "zeta="  << zeta  << std::endl;
    B.SaveLens("lens.dat", zeta);

    for(double alpha_deg=0.2;alpha_deg<=50;alpha_deg+=0.2)
    {
        const double ans = B.profile(Deg2Rad(alpha_deg), theta, zeta, &fp);
        fp << "\n";
        rp("%.15g %.15g %.15g\n", alpha_deg, ans, B.param[BRIDGE_U] );
        //rp("%.15g %.15g %.15g\n", alpha_deg, ans, B.reduced_rate(B.param) );
    }

    double alpha_opt = B.find_alpha(theta,zeta);
    if(alpha_opt>0)
    {
        std::cerr << "alpha_opt=" << Rad2Deg(alpha_opt) << std::endl;
        {
            ios::wcstream ap("alpha_opt.dat");
            (void)B.profile(alpha_opt,theta,zeta,&ap,true);
            const double urise = B.get_user_rise(alpha_opt,theta,zeta);
            ap << "\n";
            ap("%.15g %.15g 0\n",  -1.0,       urise );
            ap("%.15g %.15g 0\n",  B.radii[1], urise );

        }

        if(false)
        {
            ios::wcstream sp("slice.dat");
            for(size_t i=1;i<=B.heights.size();++i)
            {
                sp("%.15g %.15g %.15g %.15g\n",B.heights[i],B.radii[i],B.slices[i],B.volume[i]);
            }
        }

    }


    if(false)
    {
        ios::wcstream ap("a_of_theta.dat");
        for(double th=1;th<=180;++th)
        {
            const double alpha = B.find_alpha(Deg2Rad(th),zeta);
            if(alpha>0)
            {
                ap("%.15g %.15g\n", th,Rad2Deg(alpha));
            }
        }
    }


}
YOCTO_PROGRAM_END()

