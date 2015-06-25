#include "yocto/program.hpp"
#include "bridge.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/code/rand.hpp"
YOCTO_PROGRAM_START()
{
    alea_init();
    if(argc<=3)
        throw exception("need radius clength theta");
    const double R = strconv::to<double>(argv[1],"R");
    const double C = strconv::to<double>(argv[2],"C");

    const Lens::Pointer lens( new SphericalLens(R));
    Bridge bridge(lens, C);


    const double theta  = strconv::to<double>(argv[3],"theta");
    std::cerr << "Looking For Hmax, theta=" << theta << std::endl;
    const double hmax   = bridge.FindHmax(theta);
    std::cerr << "\thmax=" << hmax << std::endl;

    for(size_t i=1;i<=10;++i)
    {
        const double h = hmax * alea<double>();
        const double a = bridge.FindAlpha(h, theta);
        std::cerr << "alpha(h=" << h << ")=" << a << std::endl;
        const double t = bridge.FindTheta(h, a);
        std::cerr << "\trecovered theta=" << t << std::endl;
    }



#if 0
    const double alpha0 = bridge.FindAlpha(0,theta);
    std::cerr << "alpha(h=0,theta=" << theta << ")=" << alpha0 << std::endl;
    lens->Output(0);
    std::cerr << "out=" << bridge.FinalRadius(0, theta, alpha0, true) << std::endl;
    std::cerr << "fin=" << bridge.Y[1] << std::endl;
    std::cerr << "h_max(theta=" << theta << ")=" << hmax << std::endl;
#endif

}
YOCTO_PROGRAM_END()
