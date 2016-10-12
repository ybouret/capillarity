#include "setup.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/fit/glsf-spec.hpp"
#include "yocto/container/utils.hpp"
#include "yocto/math/stat/descr.hpp"
#include "yocto/seem/evaluator.hpp"

YOCTO_PROGRAM_START()
{
    if(argc<=3)
    {
        throw exception("usage: %s R0 capillary_length datafiles", program);
    }

    Seem::Evaluator seem;
    const string    expr_R0     = argv[1];
    const string    expr_caplen = argv[2];
    const double    seem_R0     = seem.run(expr_R0);
    const double    seem_caplen = seem.run(expr_caplen);

    std::cerr << "R0=" << seem_R0 << std::endl;
    std::cerr << "caplen=" << seem_caplen << std::endl;


    Application app(new threading::crew(true),
                    seem_R0,
                    seem_caplen
                    );


    Setup          &setup = app.setup;
    vector<double> &zeta  = app.zeta;
    vector<double> &alpha = app.alpha;
    vector<double> &theta = app.theta;
    vector<double> &dzeta = app.dzeta;
    vector<double> &znew  = app.znew;
    vector<double> &zfit  = app.zfit;
    vector<double>  coef(2);

    for(int i=3;i<=3;++i)
    {
        const string filename = argv[i];
        const string rootname = vfs::get_base_name(filename);

        //______________________________________________________________________
        //
        // load data
        //______________________________________________________________________
        vector<double> height;
        vector<double> surface;
        vector<double> surffit;

        {
            data_set<double> ds;
            ds.use(2, height);
            ds.use(1, surface);
            ios::icstream fp(filename);
            ds.load(fp);
        }
        size_t N0 = height.size();
        surffit.make(N0);
        std::cerr << "#data=" << N0 << std::endl;

        if(false)
        {
            ios::wcstream fp("output.dat");
            for(size_t i=1;i<=N0;++i)
            {
                fp("%.15g %.15g\n", height[i], surface[i]);
            }
        }

        zeta.free();
        alpha.free();
        theta.free();
        dzeta.free();
        znew.free();
        vector<double> area(N0,as_capacity);

        const double   S0 = setup.S0;
        for(size_t i=1;i<=N0;++i)
        {
            const double s = surface[i];
#if 0
            if(area_max>area_min)
            {
                if(s<area_min||s>area_max)
                {
                    continue;
                }
            }
#endif
            zeta.push_back( height[i]/setup.R0 );
            const double ss = s / S0;
            if(ss>1)
                throw exception("surface is too high");
            alpha.push_back( asin( sqrt(ss) ) );
            area.push_back(s);
        }
        const size_t N = zeta.size();
        std::cerr << "using #data=" << N << std::endl;
        dzeta.make(N);
        for(size_t i=1;i<=N;++i)
        {
            dzeta[i] = zeta[i] - zeta[1];
        }
        theta.make(N);
        znew.make(N);
        zfit.make(N);

        //______________________________________________________________________
        //
        // extracting data
        //______________________________________________________________________
        std::cerr << "-- extracting theta using #cores=" << app.cores << std::endl;

        // initialize search
        vector<double> zeta0(N);
        tao::set(zeta0,zeta);
        if( !app.extract_theta() )
        {
            std::cerr << "COULDN'T EXTRACT THETA!" << std::endl;
            continue;
        }

        {
            string outname = rootname;
            vfs::change_extension(outname, "theta_raw.dat");
            ios::wcstream fp(outname);
            fp("#area theta h\n");
            for(size_t i=1;i<=N;++i)
            {
                const double t = Rad2Deg(theta[i]);
                fp("%.15g %.15g %.15g\n", area[i], t, setup.R0 * zeta0[i]  );
            }
        }

        
    }

}
YOCTO_PROGRAM_END()
