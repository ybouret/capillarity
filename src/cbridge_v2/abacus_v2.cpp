#include "bridge.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"

#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
    Lua::State VM;
    lua_State *L = VM();

    Lua::Config::DoString(L,"ftol=1e-7;");
    Lua::Config::DoString(L,"angle_control=1;");
    Lua::Config::DoString(L,"shift_control=0.01;");
    Lua::Config::DoString(L,"R0=80;");
    Lua::Config::DoString(L,"lambda=2.72;");
    Lua::Config::DoString(L,"resolution=1e-3");

    for(int i=1;i<argc;++i)
    {
        Lua::Config::DoString(L,argv[i]);
    }

    //Bridge::curvature_coeff = 1;

    Bridge B(L);
    ios::ocstream::overwrite("xi_max.dat");
    const string filename = vformat("abacus_m%u.dat", unsigned(Bridge::curvature_coeff));
    ios::ocstream::overwrite(filename);

    const size_t N      = 80;
    for(int theta_deg = 60; theta_deg <= 170; theta_deg += 10 )
    {
        std::cerr << std::endl;
        std::cerr << "theta=" << theta_deg << std::endl;
        const double theta = Deg2Rad(double(theta_deg));
        double alpha_top = 0;
        double zeta_max  = 0;
        const double Xi_max = B.find_Xi_max(theta, alpha_top, zeta_max);
        {
            ios::acstream fp("xi_max.dat");
            fp("%.15g %.15g %.15g\n", double(theta_deg), Xi_max, zeta_max);
        }

        const double Xi_min = -Xi_max/2;
        //const double Xi_min = 0;
        for(size_t i=0;i<=N;++i)
        {
            const double Xi    = Xi_min + double(i)*(Xi_max-Xi_min)/N;
            const double alpha = B.find_alpha(theta,Xi);
            //std::cerr << "Xi=" << Xi << ", alpha=" << Rad2Deg(alpha) << std::endl;

            if(alpha<0) throw exception("invalid Xi=%g for theta=%g\n", Xi, double(theta_deg) );
            (void) B.profile(alpha, theta, Xi, NULL, true);
            const double usink = B.compute_user_sink(alpha, theta, Xi);
            const double zeta  = Xi - usink;
            ios::acstream fp(filename);
            fp("%.15g %.15g %.15g\n", Xi, Square(sin(alpha)), zeta );
        }

        {
            ios::acstream fp(filename);
            fp << "\n";
        }

    }

}
YOCTO_PROGRAM_END()

