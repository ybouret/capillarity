#include "cbridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

    const double xi        = Lua::Config::Get<lua_Number>(L,"xi");
    const double alpha_deg = Lua::Config::Get<lua_Number>(L,"alpha");
    const double alpha     = Deg2Rad(alpha_deg);

    std::cerr << "alpha=" << alpha_deg << std::endl;
    std::cerr << "xi   =" << xi        << std::endl;

    Bridge::SaveLens("lens.dat", xi);
    ios::ocstream::overwrite("theta_xi.dat");

    bool         is_flat  = false;
    double       shift    = 0;
    const double theta_xi = B.find_theta_by_xi(alpha,xi,&is_flat,shift);

    std::cerr << "theta_xi=" << theta_xi << std::endl;


    if(theta_xi>=0)
    {
        const double zeta = xi+shift;
        std::cerr << "shift=" << shift << std::endl;
        std::cerr << "zeta =" << zeta  << std::endl;
        std::cerr << "THETA_Xi=" << Rad2Deg(theta_xi) << std::endl;

        ios::wcstream fp("theta_xi.dat");
        if(is_flat)
        {
            B.compute_start(alpha,theta_xi,zeta);
            for(double xx=0;xx<=1;xx+=0.1)
            {
                //fp("%g %g\n", B.start_u+xx, B.start_v );
            }
        }
        else
        {
            (void)B.profile(alpha, theta_xi, zeta, NULL, true);
            for(size_t i=1;i<=B.heights.size();++i)
            {
                fp("%.15g %.15g\n", B.radii[i], B.heights[i]-shift);
            }
        }
    }
    else
    {
        std::cerr << "No Bridge!" << std::endl;
    }

}
YOCTO_PROGRAM_END()
