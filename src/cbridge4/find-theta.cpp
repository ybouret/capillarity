#include "bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

    const double zeta      = Lua::Config::Get<lua_Number>(L, "zeta" );
    const double alpha_deg = Lua::Config::Get<lua_Number>(L,"alpha");
    const double alpha     = Deg2Rad(alpha_deg);
    std::cerr << "alpha=" << alpha_deg << std::endl;

    B.SaveLens("lens.dat", zeta);

    {
        ios::wcstream fp("profile.dat");
        ios::wcstream rp("results.dat");

        for(double theta_deg = 1; theta_deg <= 179; theta_deg += 1)
        {
            const double theta = Deg2Rad(theta_deg);
            const double ans = B.profile(alpha,theta,zeta,&fp,false);
            fp << "\n";
            rp("%g %g\n", theta_deg, ans );
        }
    }

}
YOCTO_PROGRAM_END()
