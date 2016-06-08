#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/fit/glsf-spec.hpp"

YOCTO_PROGRAM_START()
{
    Bridge B(0.001,  // angular search
             1e-5,   // integrator ftol
             0.1,    // max angular turn in degrees
             1.0/100 // max speed rate in degrees
             );
    
    int    iarg = 0;
    
    
    if(argc>++iarg)
    {
        B.mu   = strconv::to<double>(argv[iarg],"mu");
    }
    std::cerr << "mu=" << B.mu << std::endl;
    
    
    
    //__________________________________________________________________________
    //
    // find zeta_max
    //__________________________________________________________________________
    const size_t N=101;
    vector<double> zeta(N);
    vector<double> surf(N);
    vector<double> sfit(N);

    GLS<double>::Function      poly = _GLS::Create<double,_GLS::Polynomial>();
    GLS<double>::Samples       samples(1);
    (void) samples.append(zeta, surf, sfit);
    
    const size_t   m=3;
    vector<double> aorg(m);
    vector<double> aerr(m);
    vector<bool>   used(m,true);
    samples.prepare(m);

    
    const string heur_name = vformat("heur%g.dat",B.mu);
    ios::ocstream::overwrite(heur_name);
    
    const string sfit_name = vformat("sfit%g.dat",B.mu);
    ios::ocstream::overwrite(sfit_name);

    const string coef_name = vformat("coef%g.dat",B.mu);
    ios::ocstream::overwrite(coef_name);

    
    for(int t_deg = 40; t_deg <= 170; t_deg+=10)
    {
        std::cerr << "theta=" << t_deg << std::endl; std::cerr.flush();
        const double theta     = Deg2Rad( double(t_deg) );
        const double zeta_max  = B.compute_zeta_max( theta);
        std::cerr << "\tzeta_max=" << zeta_max << std::endl;
        const double zeta_half = zeta_max/2;
        const double zeta_min  = -0.1;
        for(size_t i=1;i<=N;++i)
        {
            const double ztmp = zeta_min + (i-1)*(zeta_half-zeta_min)/double(N-1);
            zeta[i] = ztmp;
            const double alpha = B.find_alpha(theta,ztmp);
            surf[i]            = numeric<double>::pi * Square( sin(alpha) );
            ios::acstream hp(heur_name);
            hp("%.15g %.15g\n", zeta[i], surf[i]);
            std::cerr << "\tzeta=" << zeta[i] << " => " << surf[i] << std::endl;
        }
        {
            ios::acstream hp(heur_name);
            hp << "\n";
        }
        

        if(samples.fit_with(poly, aorg, used, aerr))
        {
            GLS<double>::display(std::cerr,aorg, aerr);
        }
        
        {
            ios::acstream fp(sfit_name);
            for(size_t i=0;i<=N;++i)
            {
                const double ztmp = zeta_min + (i*(zeta_max-zeta_min))/double(N);
                const double stmp = poly(ztmp,aorg);
                fp("%.15g %.15g\n",ztmp,stmp);
            }
            fp << "\n";
        }
        
        {
            ios::acstream fp(coef_name);
            fp("%.15g",double(t_deg));
            for(size_t i=1;i<=m;++i)
            {
                fp(" %.15g",aorg[i]);
            }
            fp << "\n";
        }
        
        
       

    }
    
}
YOCTO_PROGRAM_END()
