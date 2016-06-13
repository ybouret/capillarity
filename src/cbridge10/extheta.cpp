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
    if(argc<=4)
    {
        throw exception("usage: %s R0 capillary_length area_min  area_max [datafiles]", program);
    }


    Application app(new threading::crew(true),
                    strconv::to<double>(argv[1],"R0"),
                    strconv::to<double>(argv[2],"capillary_length")
                    );

    const double area_min = strconv::to<double>(argv[3]);
    const double area_max = strconv::to<double>(argv[4]);

    Setup  &setup = app.setup;
    vector<double> &zeta  = app.zeta;
    vector<double> &alpha = app.alpha;
    vector<double> &theta = app.theta;
    vector<double> &dzeta = app.dzeta;
    vector<double> &znew  = app.znew;
    vector<double> &zfit  = app.zfit;
    vector<double>  coef(2);

    for(int i=5;i<argc;++i)
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
            ds.use(1, height);
            ds.use(2, surface);
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

#if 0
        //______________________________________________________________________
        //
        // autocut
        //______________________________________________________________________
        const Direction        hdir = (height[1] > height[N0]) ? Pushing : Pulling;
        switch (hdir)
        {
            case Pushing:
                std::cerr << "Pushing/Enfoncement" << std::endl;
                break;

            case Pulling:
                std::cerr << "Pulling/Tirage"   << std::endl;
                break;
        }

        const double CutOff = parted.split(hdir, height, surface, surffit);
        {
            string outname = rootname;
            vfs::change_extension(outname, "cut.dat");
            ios::wcstream fp(outname);
            for(size_t i=1;i<=N0;++i)
            {
                fp("%g %g %g\n", height[i], surface[i], surffit[i]);
            }
        }

        //______________________________________________________________________
        //
        // build data to extract
        //______________________________________________________________________
        vector<double> zeta(N0,as_capacity);  //!< reduced height
        vector<double> alpha(N0,as_capacity); //!< angle
        const size_t   N = setup.isolate(hdir, zeta, alpha, height, surface, CutOff);
        if(N<=0)
        {
            throw exception("No data remaining");
        }
        std::cerr << "using #data=" << N << std::endl;
#endif


        zeta.free();
        alpha.free();
        theta.free();
        dzeta.free();
        znew.free();

        const double   S0 = setup.S0;

        for(size_t i=1;i<=N0;++i)
        {
            const double s = surface[i];
            if(area_max>area_min)
            {
                if(s<area_min||s>area_max)
                {
                    continue;
                }
            }
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

        {
            string outname = rootname;
            vfs::change_extension(outname, "alpha.dat");
            ios::wcstream fp(outname);
            for(size_t i=1;i<=N;++i)
            {
                fp("%g %g %g\n", zeta[i] * setup.R0, numeric<double>::pi * Square( setup.R0 * sin(alpha[i]) ), Rad2Deg(alpha[i]) );
            }
        }


        const double zeta_max_exp = find_max_of(zeta);
        const double zeta_min_exp = find_min_of(zeta);
        std::cerr << "zeta_max_exp=" << zeta_max_exp << std::endl;
        std::cerr << "zeta_min_exp=" << zeta_min_exp << std::endl;

        //______________________________________________________________________
        //
        // extracting data
        //______________________________________________________________________
        std::cerr << "-- extracting theta using #cores=" << app.cores << std::endl;
        app.extract_theta();



        {
            string outname = rootname;
            vfs::change_extension(outname, "theta.dat");
            ios::wcstream fp(outname);

            for(size_t i=1;i<=N;++i)
            {
                const double t = round(10.0*Rad2Deg(theta[i]))/10.0;

                fp("%g %g %g\n", setup.R0 * zeta[i] , t , Rad2Deg(alpha[i]) );
            }
        }

        for(size_t iter=1;iter<=5;++iter)
        {
            double theta_ave = 0;
            compute_average(theta_ave, theta);
            std::cerr << "theta_ave=" << Rad2Deg(theta_ave) << std::endl;
            std::cerr << "-- inversion..." << std::endl;
            app.compile<double,double>();
            app.call(znew,alpha, &theta_ave);

            for(size_t i=1;i<=N;++i)
            {
                znew[i] -= zeta[i];
            }

            GLS<double>::Sample sample(dzeta,znew,zfit);
            _GLS::Polynomial<double>::Start(sample,coef);
            std::cerr << "coef=" << coef << std::endl;
            for(size_t i=1;i<=N;++i)
            {
                znew[i] = zeta[i] + _GLS::Polynomial<double>::Eval(dzeta[i],coef);
            }
            
            {
                ios::wcstream fp("znew.dat");
                for(size_t i=1;i<=N;++i)
                {
                    fp("%g %g %g\n", double(i), zeta[i], znew[i] );
                }
            }

            for(size_t i=1;i<=N;++i)
            {
                zeta[i] = znew[i];
            }
            app.extract_theta();
        }
        
        
    }
    
    
    
}
YOCTO_PROGRAM_END()
