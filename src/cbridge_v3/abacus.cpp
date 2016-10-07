#include "cbridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

    const double hmin = Lua::Config::Get<lua_Number>(L,"hmin");

    const size_t N    = 100;
    const string filename = "abacus.dat";
    ios::ocstream::overwrite(filename);

    //const double R3 = Cube(B.R0);
    for(int theta_deg = 60; theta_deg <= 170; theta_deg += 10 )
    {
        const double theta = Deg2Rad(double(theta_deg));
        std::cerr << std::endl;
        std::cerr << "theta=" << theta_deg << std::endl;

        const double zeta_max = B.find_zeta_max(theta);
        std::cerr << "\tzeta_max=" << zeta_max << std::endl;

        const double zeta_min = hmin/B.R0;

        for(size_t i=0;i<=N;++i)
        {
            //const double zeta  = (i<=0) ? zeta_min : zeta_min + ((zeta_max-zeta_min)*i)/double(N);
            double zeta        = 0;
            if(i<=0)
            {
                zeta = zeta_min;
            }
            else
            {
                if(i>=N)
                {
                    zeta = zeta_max;
                }
                else
                {
                    zeta =  zeta_min + ((zeta_max-zeta_min)*i)/double(N);
                }
            }
            const double alpha = B.find_alpha(theta, zeta, NULL);
            std::cerr.flush();
            fprintf(stderr,"\t\talpha=%10.6f   \r",Rad2Deg(alpha));
            fflush(stderr);
            (void)B.profile(alpha,theta,zeta,NULL);
            const double h     = zeta * B.R0;
            const double A     = numeric<double>::pi * Square( B.R0 * sin(alpha) );
            const double V     = B.last_space();// * R3;
            ios::acstream fp(filename);

#if 0
            double cvol = 0;
            if(zeta>=0)
            {
                cvol = Bridge::CapVolume(B.start_v-zeta);
            }
            else
            {
                cvol = Bridge::CapVolume(B.start_v-zeta)-Bridge::CapVolume(-zeta);
            }
#endif
            fp("%.15g %.15g %.15g %.15g %.15g %.15g\n",h,A,V, B.param[BRIDGE_Q], B.param[BRIDGE_q], B.last_cylinder_space());
        }
        std::cerr << std::endl;
        
        {
            ios::acstream fp(filename); fp << "\n";
        }

    }

}
YOCTO_PROGRAM_END()
