#include "cbridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

    const double zeta      = Lua::Config::Get<lua_Number>(L,"zeta");
    const double alpha_deg = Lua::Config::Get<lua_Number>(L,"alpha");
    const double alpha     = Deg2Rad(alpha_deg);

    std::cerr << "alpha=" << alpha_deg << std::endl;
    std::cerr << "zeta =" << zeta  << std::endl;

    Bridge::SaveLens("lens.dat", zeta);

    if(true)
    {
        ios::wcstream fp("profile.dat");
        ios::wcstream rp("results.dat");

        for(double theta_deg=1;theta_deg<=179; theta_deg += 1)
        {
            const double theta = Deg2Rad(theta_deg);
            const double ans   = B.profile(alpha,theta,zeta, &fp);
            fp << "\n";
            rp("%.15g %.15g\n",theta_deg,ans);
        }
    }

    bool         is_flat   = false;
    const double theta_opt = B.find_theta(alpha,zeta,&is_flat);
    if(theta_opt>=0)
    {
        std::cerr << "THETA=" << Rad2Deg(theta_opt) << std::endl;
        B.compute_start(alpha,theta_opt,zeta);
        ios::wcstream fp("theta_opt.dat");
        if(is_flat)
        {
            for(double xx=0;xx<=1;xx+=0.1)
            {
                fp("%g %g\n", B.start_u+xx, B.start_v );
            }
        }
        else
        {
            (void)B.profile(alpha, theta_opt, zeta, &fp);
        }
    }
    else
    {
        std::cerr << "No Bridge!" << std::endl;
        ios::ocstream::overwrite("theta_opt.dat");
    }

}
YOCTO_PROGRAM_END()
