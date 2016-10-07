#include "cbridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

    const double zeta      = Lua::Config::Get<lua_Number>(L,"zeta");
    const double theta_deg = Lua::Config::Get<lua_Number>(L,"theta");
    const double theta     = Deg2Rad(theta_deg);

    std::cerr << "theta=" << theta_deg << std::endl;
    std::cerr << "zeta =" << zeta      << std::endl;

    Bridge::SaveLens("lens.dat", zeta);

    if(true)
    {
        ios::wcstream fp("profile.dat");
        ios::wcstream rp("results.dat");

        for(double alpha_deg=1;alpha_deg<=90; alpha_deg += 1)
        {
            const double alpha = Deg2Rad(alpha_deg);
            const double ans   = B.profile(alpha,theta,zeta, &fp);
            fp << "\n";
            rp("%.15g %.15g\n",alpha_deg,ans);
        }
    }

    bool is_flat = false;
    const double alpha_opt = B.find_alpha(theta,zeta,&is_flat);
    if(alpha_opt>=0)
    {
        std::cerr << "ALPHA=" << Rad2Deg(alpha_opt) << std::endl;
        B.compute_start(alpha_opt,theta,zeta);
        ios::wcstream fp("alpha_opt.dat");
        if(is_flat)
        {
            for(double xx=0;xx<=1;xx+=0.1)
            {
                fp("%g %g\n", B.start_u+xx, B.start_v );
            }
        }
        else
        {
            (void)B.profile(alpha_opt, theta, zeta, &fp);
        }
    }
    else
    {
        std::cerr << "No Bridge!" << std::endl;
        ios::ocstream::overwrite("alpha_opt.dat");
    }


}
YOCTO_PROGRAM_END()
