#include "bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"
    ios::ocstream::overwrite("abacus_alpha.dat");


    ios::ocstream::overwrite("h_max.dat");

    for(double theta_deg=100; theta_deg <= 170; theta_deg += 10 )
    {
        std::cerr << "THETA=" << theta_deg << std::endl;
        std::cerr << "\tlooking for zeta_max..." << std::endl;
        const double theta    = Deg2Rad(theta_deg);
        const double zeta_max = B.find_zeta_max(theta);
        std::cerr << "\tzeta_max=" << zeta_max << std::endl;
        ios::ocstream::echo("h_max.dat","%g %g\n", theta_deg, zeta_max * B.R0);

        

    }

}
YOCTO_PROGRAM_END()

