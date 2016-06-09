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
        GLS<double>::Function *part = 0;
        switch (hdir)
        {
            case Pushing:
                std::cerr << "Pushing/Enfoncement" << std::endl;
                part = &cut_push;
                break;
                
            case Pulling:
                std::cerr << "Pulling/Tirage"   << std::endl;
                part = &cut_pull;
                break;
        }

        surffit.make(N0);
        GLS<double>::Samples       samples;
        const GLS<double>::Sample &sample = samples.append(height,surface,surffit);

        const size_t   nvar = 3;
        vector<double> aorg(nvar);
        vector<double> aerr(nvar);
        vector<bool>   used(nvar,true);
        {
            vector<double> p(2);
            if(!_GLS::Polynomial<double>::Start(sample, p))
            {
                throw exception("couldn't guess initial parameters");
            }
            std::cerr << "p=" << p << std::endl;
            aorg[1] = p[1];
            aorg[2] = p[2];
            aorg[3] = 0.5*(height[1]+height[N0]);
        }
        samples.prepare(nvar);


        if( !samples.fit_with(*part, aorg, used, aerr, &cut_cb) )
        {
            throw exception("couldn't part data...");
        }
        GLS<double>::display(std::cerr, aorg, aerr);

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
        const double R0 = setup.R0;
        const double S0 = numeric<double>::pi * R0*R0;
        const double CutOff = aorg[3];
        vector<double> zeta(N0,as_capacity); //!< reduced height
        vector<double> alpha(N0,as_capacity); //!< angle
        for(size_t i=1;i<=N0;++i)
        {
            const double h_i = height[i];
            const double hh  = h_i/R0;
            const double ss  = surface[i]/S0;
            if(ss>=1)
            {
                throw exception("invalid surface %g", surface[i]);
            }
            
            const double aa = asin( sqrt(ss) );
            switch(hdir)
            {
                case Pushing:
                    if(h_i<=CutOff)
                    {
                        zeta.push_back(hh);
                        alpha.push_back(aa);
                    }
                    break;

                case Pulling:
                    if(h_i>=CutOff)
                    {
                        zeta.push_back(hh);
                        alpha.push_back(aa);
                    }
                    break;
            }
        }

        const size_t N = zeta.size();
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
                fp("%g %g %g\n", zeta[i] * R0, S0 * Square( sin(alpha[i]) ), Rad2Deg(alpha[i]) );
            }
        }

        const double zeta_max_exp = find_max_of(zeta);
        const double zeta_min_exp = find_min_of(zeta);
        std::cerr << "zeta_max_exp=" << zeta_max_exp << std::endl;
        std::cerr << "zeta_min_exp=" << zeta_min_exp << std::endl;

        vector<double> theta(N);
        for(size_t i=1;i<=N;++i)
        {
            theta[i] = setup.compute_theta(alpha[i], zeta[i]);
        }

        {
            string outname = rootname;
            vfs::change_extension(outname, "theta.dat");
            ios::wcstream fp(outname);

            for(size_t i=1;i<=N;++i)
            {
                fp("%g %g %g\n", zeta[i] ,  Rad2Deg(theta[i]), Rad2Deg(alpha[i]) );
            }
        }


        const double rate = -0.00010;
        for(size_t i=1;i<=N;++i)
        {
            theta[i] = setup.compute_theta(alpha[i], zeta[i]+rate*(i-1));
        }

        {
            string outname = rootname;
            vfs::change_extension(outname, "theta0.dat");
            ios::wcstream fp(outname);

            for(size_t i=1;i<=N;++i)
            {
                fp("%g %g %g\n", zeta[i]+rate*(i-1),  Rad2Deg(theta[i]), Rad2Deg(alpha[i]) );
            }
        }




    }
    
    
    
}
YOCTO_PROGRAM_END()
