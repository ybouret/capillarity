#include "bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

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
    
    B.find_alpha(theta,zeta,isFlat);


}
YOCTO_PROGRAM_END()
