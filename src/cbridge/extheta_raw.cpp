#include "setup.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/fit/glsf-spec.hpp"
#include "yocto/container/utils.hpp"
#include "yocto/math/stat/descr.hpp"

YOCTO_PROGRAM_START()
{
    if(argc<=3)
    {
        throw exception("usage: %s R0 capillary_length datafiles", program);
    }


    Application app(new threading::crew(true),
                    strconv::to<double>(argv[1],"R0"),
                    strconv::to<double>(argv[2],"capillary_length")
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
            vfs::change_extension(outname, "theta.dat");
            ios::wcstream fp(outname);
            fp("#H theta alpha\n");
            for(size_t i=1;i<=N;++i)
            {
                const double t = Rad2Deg(theta[i]);
                fp("%.15g %.15g %.15g\n", setup.R0 * zeta0[i] , t , Rad2Deg(alpha[i]) );
            }
        }

        
    }

}
YOCTO_PROGRAM_END()
