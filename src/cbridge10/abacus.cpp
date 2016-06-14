#include "setup.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/container/matrix.hpp"

YOCTO_PROGRAM_START()
{

    if(argc<=2)
    {
        throw exception("usage: %s R0 capillary_length",program);
    }


    Setups app( new threading::crew(true),
               strconv::to<double>(argv[1],"R0"),
               strconv::to<double>(argv[2],"capillary_length") );

    Setup       &setup    = app[0];

    const string filename = vformat("abacusR0=%g_L=%g.dat", setup.R0, setup.capillary_length );
    const string fullname = vformat("abacusR0=%g_L=%g_full.dat",setup.R0, setup.capillary_length );
    const string areaname = vformat("athetaR0=%g_L=%g.dat", setup.R0, setup.capillary_length );

    ios::ocstream::overwrite(filename);
    ios::ocstream::overwrite(fullname);
    ios::ocstream::overwrite(areaname);

    const size_t N        = 50+1;

    const int th_deg_min = 5;
    const int th_deg_max = 175;
    const int th_deg_inc = 5;
    const int n_theta    = 1+(th_deg_max-th_deg_min)/th_deg_inc;

    matrix<double> A(2*n_theta,N);
    app.compile<double,double>();
    for(int th_deg = th_deg_min, iz=1,ia=2,count=1; th_deg <= th_deg_max; th_deg += th_deg_inc,iz+=2,ia+=2,++count)
    {
        std::cerr << "-- theta=" << th_deg << std::endl;
        const double  theta = Deg2Rad(double(th_deg));
        array<double> &zeta = A[iz];
        array<double> &area = A[ia];
        std::cerr << "\tcomputing zeta_max..." << std::endl;
        const double zeta_max = setup.bridge.compute_zeta_max(theta);
        std::cerr << "\tzeta_max=" << zeta_max << std::endl;

        zeta[1] = 0.0;
        for(size_t i=2;i<N;++i)
        {
            zeta[i] = (i*zeta_max)/double(N-1);
        }
        zeta[N] = zeta_max;
        double param = -theta;
        app.call(area,zeta,&param);

        for(size_t i=1;i<=N;++i)
        {
            zeta[i] *= setup.R0;
            area[i] *= setup.R02;
        }

        {
            ios::acstream fp(areaname);
            fp("%.15g %d\n", area[1], th_deg);
        }

        {
            ios::acstream fp(fullname);
            for(size_t i=1;i<=N;++i)
            {
                fp("%.15g %.15g\n", zeta[i], area[i]);
            }
            fp << '\n';
        }

        {
            ios::wcstream fp(filename);
            fp << "#";
            for(int j=1;j<=count;++j)
            {
                const int th_tmp = th_deg_min + (j-1) * th_deg_inc;
                if(j>1) fp << ' ';
                fp("h%d A%d",th_tmp,th_tmp);
            }
            fp << "\n";
            for(size_t i=1;i<=N;++i)
            {
                for(size_t j=1,jz=1,ja=2;j<=count;++j,jz+=2,ja+=2)
                {
                    if(j>1) fp << ' ';
                    fp("%.15g %.15g", A[jz][i], A[ja][i]);
                }
                fp << '\n';
            }
        }
    }


    
}
YOCTO_PROGRAM_END()
