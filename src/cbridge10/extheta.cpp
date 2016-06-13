#include "setup.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/fit/glsf-spec.hpp"
#include "yocto/container/utils.hpp"

YOCTO_PROGRAM_START()
{
    if(argc<=2)
    {
        throw exception("usage: %s R0 capillary_length [datafiles]", program);
    }

    Setup  setup( strconv::to<double>(argv[1],"R0"), strconv::to<double>(argv[2],"capillary_length"));
    Parted parted;
    GLS<double>::Function cut_pull( &parted, & Parted::pull );
    GLS<double>::Function cut_push( &parted, & Parted::push );
    GLS<double>::Callback cut_cb( &parted, & Parted::callback);

    for(int i=3;i<argc;++i)
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
        {
            ios::wcstream fp("output.dat");
            for(size_t i=1;i<=N0;++i)
            {
                fp("%.15g %.15g\n", height[i], surface[i]);
            }
        }

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

        vector<double> theta(N);
        double         aveth = 0;

        //______________________________________________________________________
        //
        // extracting data
        //______________________________________________________________________
        aveth=0;
        for(size_t i=1;i<=N;++i)
        {
            theta[i] = setup.compute_theta(alpha[i], zeta[i]);
            aveth   += theta[i];
        }
        aveth /= N;


        {
            string outname = rootname;
            vfs::change_extension(outname, "theta.dat");
            ios::wcstream fp(outname);

            for(size_t i=1;i<=N;++i)
            {
                fp("%g %g %g\n", zeta[i] ,  Rad2Deg(theta[i]), Rad2Deg(alpha[i]) );
            }
        }

#if 0
        {
            vector<double> coord(N);
            vector<double> dzfit(N);
            vector<double> dzeta(N);
            GLS<double>::Sample S(coord,dzeta,dzfit);
            vector<double> zparam(2);

            for(size_t i=1;i<=N;++i) coord[i] = i;

            ios::wcstream lp("corr.dat");

            string outname = rootname;
            vfs::change_extension(outname, "corr.dat");
            ios::wcstream fp(outname);
            for(size_t i=1;i<=N;++i)
            {
                fp("%.15g %.15g 0\n", coord[i],zeta[i]);
            }
            fp << "\n";


            double t0 = aveth;
            double H0 = setup.rebuild(t0, zparam, coord, alpha, zeta, dzeta, dzfit);
            std::cerr << "t0=" << Rad2Deg(t0) << std::endl;
            std::cerr << "H0=" << H0 << std::endl;
            lp("%g %g %g %g\n", Rad2Deg(t0), H0, zparam[2], zparam[1]);

            double t0deg = Rad2Deg(t0);
            double tp    = Deg2Rad( ceil(t0deg) + 1 );
            double Hp    = setup.rebuild(tp, zparam, coord, alpha, zeta, dzeta, dzfit);
            std::cerr << "tp=" << Rad2Deg(tp) << std::endl;
            std::cerr << "Hp=" << Hp << std::endl;

            if(Hp<H0)
            {
                lp("%g %g %g %g\n", Rad2Deg(tp), Hp, zparam[2], zparam[1]);
                for(int i=2;i<=90;i+=5)
                {
                    tp = Deg2Rad( ceil(t0deg)+i );
                    if(Rad2Deg(tp)>=175) break;
                    Hp = setup.rebuild(tp, zparam, coord, alpha, zeta, dzeta, dzfit);
                    std::cerr << "tp=" << Rad2Deg(tp) << std::endl;
                    std::cerr << "Hp=" << Hp << std::endl;
                    lp("%g %g %g %g\n", Rad2Deg(tp), Hp, zparam[2], zparam[1]);
                }
            }


            double tm = Deg2Rad( floor(Rad2Deg(aveth)) - 1 );
            double Hm = setup.rebuild(tm, zparam, coord, alpha, zeta, dzeta, dzfit);
            std::cerr << "tm=" << Rad2Deg(tm) << std::endl;
            std::cerr << "Hm=" << Hm << std::endl;

            if(Hm<H0)
            {
                lp("%g %g %g %g\n", Rad2Deg(tm), Hm, zparam[2], zparam[1]);
                for(int i=2;i<=90;i+=5)
                {
                    tm = Deg2Rad( floor(t0deg)-i );
                    if(Rad2Deg(tm)<=5) break;
                    Hm = setup.rebuild(tm, zparam, coord, alpha, zeta, dzeta, dzfit);
                    std::cerr << "tm=" << Rad2Deg(tm) << std::endl;
                    std::cerr << "Hm=" << Hm << std::endl;
                    lp("%g %g %g %g\n", Rad2Deg(tm), Hm, zparam[2], zparam[1]);
                }

            }

        }
#endif
        
        
    }
    
    
    
}
YOCTO_PROGRAM_END()
