#include "application.hpp"
#include "yocto/program.hpp"
#include "yocto/threading/crew.hpp"
#include "yocto/seem/evaluator.hpp"

YOCTO_PROGRAM_START()
{
    if(argc<3)
    {
        throw exception("usage: %s R0 capillary_length rate",program);
    }
    Seem::Evaluator SEEM;
    const string    expr_R0      = argv[1];
    const string    expr_cap_len = argv[2];
    const string    expr_rate    = argv[3];

    const double    R0 = SEEM.run(expr_R0);
    const double    cap_len = SEEM.run(expr_cap_len);
    threading::crew par(true);

    Application     app;
    Bridge         &bridge = app.initialize(par,R0,cap_len);
    app.rate = SEEM.run(expr_rate);
    std::cerr << "mu=" << bridge.mu << std::endl;

    for(int argi = 4; argi < argc; ++argi )
    {
        const string filename = argv[argi];
        std::cerr << "-- loading " << filename << std::endl;
        app.load_from(filename);
        std::cerr << "-- building reduced variables" << std::endl;
        app.build_reduced_variables();
        std::cerr << "-- computing theta" << std::endl;
        app.compute_theta_using(par);

        {
            ios::wcstream fp("theta.dat");
            for(size_t i=1;i<=app.h.size();++i)
            {
                fp("%.15g %.15g %.15g\n", app.A[i], app.h[i], Rad2Deg(app.theta[i]));
            }
        }
    }
}
YOCTO_PROGRAM_END()

