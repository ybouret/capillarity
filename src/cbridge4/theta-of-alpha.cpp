#include "bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"

    double       zeta      = L.Get<double>("zeta");

    B.SaveLens("lens.dat", zeta);

    ios::wcstream fp("profile.dat");
    ios::ocstream::overwrite("results.dat");
    bool isFlat;
    for(double alpha_deg=0.1;alpha_deg<=90; alpha_deg += 0.1)
    {
        const double alpha = Deg2Rad(alpha_deg);
        const double theta = B.find_theta(alpha,zeta,isFlat);
        std::cerr << "alpha=" << alpha_deg << ", theta=" << Rad2Deg(theta) << std::endl;
        if(theta>0)
        {
            B.profile(alpha, theta, zeta,&fp);
            fp << "\n";
        }
        ios::ocstream::echo("results.dat","%g %g\n", alpha_deg, Rad2Deg(theta) );
    }


}
YOCTO_PROGRAM_END()
