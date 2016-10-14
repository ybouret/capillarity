#include "bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

    double       zeta      = Lua::Config::Get<lua_Number>(L, "zeta" );
    double       alpha_deg = Lua::Config::Get<lua_Number>(L,"alpha");
    double       alpha     = Deg2Rad(alpha_deg);
    std::cerr << "alpha=" << alpha_deg << std::endl;

    B.SaveLens("lens.dat", zeta);


    {
        ios::wcstream fp("profile.dat");
        ios::wcstream rp("results.dat");

        for(double theta_deg = 1; theta_deg <= 179; theta_deg += 1)
        {
            const double theta = Deg2Rad(theta_deg);
            const double ans   = B.profile(alpha,theta,zeta,&fp,false);
            fp << "\n";
            rp("%g %g %g\n", theta_deg, ans, sin(B.param[BRIDGE_A]));
        }

    }


    bool isFlat = false;
    {
        ios::wcstream fp("theta_opt.dat");

        const double theta = B.find_theta(alpha,zeta,isFlat);
        std::cerr << "THETA=" << Rad2Deg(theta) << std::endl;
        if(theta<=0)
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
            const double shift = B.compute_shift(alpha, theta, zeta);
            std::cerr << "shift=" << shift << std::endl;
        }
        
    }

}
YOCTO_PROGRAM_END()
