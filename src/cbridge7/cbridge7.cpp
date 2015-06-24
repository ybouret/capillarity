#include "yocto/program.hpp"
#include "bridge.hpp"
#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    if(argc<=3)
        throw exception("need radius clength theta");
    const double R = strconv::to<double>(argv[1],"R");
    const double C = strconv::to<double>(argv[2],"C");

    const Lens::Pointer lens( new SphericalLens(R));
    Bridge bridge(lens, C);


    const double theta  = strconv::to<double>(argv[3],"theta");
    const double alpha0 = bridge.FindAlpha(0,theta);
    std::cerr << "alpha(h=0,theta=" << theta << ")=" << alpha0 << std::endl;
    lens->Output(0);
    std::cerr << "out=" << bridge.FinalRadius(0, theta, alpha0, true) << std::endl;
    std::cerr << "fin=" << bridge.Y[1] << std::endl;
    const double hmax = bridge.FindHmax(theta);
    std::cerr << "h_max(theta=" << theta << ")=" << hmax << std::endl;
    

}
YOCTO_PROGRAM_END()
