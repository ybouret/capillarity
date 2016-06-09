#include "setup.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/fit/glsf-spec.hpp"

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
        vector<double> height;
        vector<double> surface;
        vector<double> surffit;

        {
            data_set<double> ds;
            ds.use(1, height);
            ds.use(2, surface);
            ios::icstream fp(argv[i]);
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
            ios::wcstream fp("sample0.dat");
            for(size_t i=1;i<=N0;++i)
            {
                fp("%g %g %g\n", height[i], surface[i], surffit[i]);
            }
        }


    }
    
    
    
}
YOCTO_PROGRAM_END()
