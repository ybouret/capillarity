#include "bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

#if 0
    double       zeta      = Lua::Config::Get<lua_Number>(L, "zeta" );
    double       theta_deg = Lua::Config::Get<lua_Number>(L,"theta");
    double       theta     = Deg2Rad(theta_deg);
    std::cerr << "theta=" << theta_deg << std::endl;

    B.SaveLens("lens.dat", zeta);

    {
        ios::wcstream fp("profile.dat");
        ios::wcstream rp("results.dat");
        for(double alpha_deg=1;alpha_deg<=90; alpha_deg += 0.2)
        {
            const double ans = B.profile(Deg2Rad(alpha_deg), theta, zeta, &fp);
            rp("%g %g\n",alpha_deg,ans);
            fp << "\n";
        }
    }

    bool isFlat = false;
    {
        ios::wcstream fp("alpha_opt.dat");

        const double alpha = B.find_alpha(theta,zeta,isFlat);
        std::cerr << "ALPHA=" << Rad2Deg(alpha) << std::endl;
        if(alpha<=0)
        {
            std::cerr << "No Bridge" << std::endl;
            fp("%g %g\n", B.start_u, B.start_v);
        }
        else
        {
            if(isFlat)
            {
                for(double xx=0;xx<=1.0;xx+=0.1)
                {
                    fp("%g %g\n", B.start_u+xx, B.start_v);
                }
            }
            else
            {
                (void) B.profile(alpha, theta,zeta,&fp);
            }
        }
        
    }
#endif

}
YOCTO_PROGRAM_END()
