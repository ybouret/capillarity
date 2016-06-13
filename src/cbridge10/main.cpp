#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/threading/vpu.hpp"
#include "yocto/math/round.hpp"

YOCTO_PROGRAM_START()
{

    if(argc<=3)
    {
        throw exception("usage: %s R0 capillary_length theta", program);
    }
    const double R0                = strconv::to<double>(argv[1],"R0");
    const double capillary_length  = strconv::to<double>(argv[2],"mu");
    const double theta_deg         = strconv::to<double>(argv[3],"theta");
    const double theta             = Deg2Rad(theta_deg);
    const double R02               = R0*R0;

    std::cerr << "-- allocating resources" << std::endl;
    threading::processing_unit<DefaultBridge> cpu(new threading::crew(true));
    for(size_t i=0;i<cpu.cores;++i)
    {
        cpu.push_back().set_mu(R0,capillary_length);
    }
    Bridge &B = cpu[0];

    std::cerr << "-- computing zeta_max" << std::endl;
    const double zeta_max = B.compute_zeta_max(theta);
    std::cerr << "zeta_max=" << zeta_max << std::endl;

    std::cerr << "-- computing simulated range..." << std::endl;
    const double zeta_top = 0.75 * zeta_max;
    const double h_top    = R0 * zeta_top;
    const double dh       = log_round_floor(h_top/50.0);
    std::cerr << "h_top=" << h_top << std::endl;
    std::cerr << "dh   =" << dh    << std::endl;
    const size_t N        = 1+size_t(floor(h_top/dh));
    const size_t N1       = N-1;
    std::cerr << "h:0.." << N1*dh << std::endl;
    std::cerr << "in " << N << " points..." << std::endl;

    std::cerr << "-- building points..." << std::endl;
    vector<double> zeta(N);
    vector<double> height(N);
    vector<double> surf(N);

    for(size_t i=1;i<=N;++i)
    {
        height[i] = (i-1)*dh;
        zeta[i]   = height[i]/R0;
        surf[i]   = 0;
    }

    std::cerr << "-- compute points..." << std::endl;
    cpu.compile<double,double>();
    for(size_t i=0;i<cpu.cores;++i)
    {
        cpu[i].current_theta = theta;
    }
    int choice = BRIDGE_ZETA_TO_SURF;

    cpu.call(surf,zeta,&choice);
    for(size_t i=1;i<=N;++i)
    {
        surf[i] *= R02;
    }

    {
        ios::wcstream fp("simulated.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g\n", height[i], surf[i]);
        }
    }


}
YOCTO_PROGRAM_END()
