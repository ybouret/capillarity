#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/fit/glsf-spec.hpp"
#include "yocto/math/fcn/drvs.hpp"

YOCTO_PROGRAM_START()
{
    Bridge B(0.001,  // angular search
             1e-5,   // integrator ftol
             0.1,    // max angular turn in degrees
             1.0/100 // max speed rate in degrees
             );


    const size_t   N = 51;
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

    const double dzeta = 0.005;
    for(int iarg=1;iarg<argc;++iarg)
    {
        B.mu   = strconv::to<double>(argv[iarg],"mu");

        std::cerr << "mu=" << B.mu << std::endl;

        //derivative<double> drvs;

        const string param_name = vformat("param%g.dat",B.mu);
        ios::ocstream::overwrite(param_name);

        const string heur_name = vformat("heur%g.dat",B.mu);
        ios::ocstream::overwrite(heur_name);

        const string sfit_name = vformat("sfit%g.dat",B.mu);
        ios::ocstream::overwrite(sfit_name);


        const string delta_name = vformat("delta%g.dat",B.mu);
        ios::ocstream::overwrite(delta_name);

        const string coef_name = vformat("coef%g.dat",B.mu);
        ios::ocstream::overwrite(coef_name);


        for(int t_deg = 40; t_deg <= 170; t_deg+=10)
        {
            const double theta  = Deg2Rad(double(t_deg));
            std::cerr << "theta=" << t_deg << std::endl;

            std::cerr << "\t-- computing surface@zeta=0" << std::endl;
            const double surf0  = B.surface0(theta);
            std::cerr << "\t   surface0=" << surf0 << std::endl;

            std::cerr << "\t-- computing slope@zeta=0" << std::endl;
            B.current_theta = theta;
            const double slope  = (surf0-B.surface_of_zeta(-dzeta))/(dzeta);
            std::cerr << "\t-- slope approx " << slope << std::endl;

            std::cerr << "\t-- computing zeta_max" << std::endl;
            const double zeta_max = B.compute_zeta_max(theta);
            std::cerr << "\t   zeta_max=" << zeta_max << std::endl;

            std::cerr << "\t-- computing surf_min" << std::endl;
            const double surf_min = B.find_surface(theta,zeta_max);
            std::cerr << "\t   surf_min=" << surf_min << std::endl;

            {
                ios::acstream fp(param_name);
                fp("%.15g %.15g %.15g %.15g %.15g\n", double(t_deg), surf0, slope, zeta_max, surf_min);
            }


            const double zeta_min = 0;
            const double zeta_top = zeta_max/2;
            for(size_t i=1;i<=N;++i)
            {
                const double ztmp = (i<=1) ? zeta_min : ( (i>=N) ? zeta_top : zeta_min + (i-1)*(zeta_top-zeta_min)/double(N-1) );
                const double stmp = B.find_surface(theta,ztmp);
                std::cerr << "\t\t" << ztmp << " " << stmp << std::endl;
                { ios::acstream fp(heur_name); fp("%.15g %.15g\n", ztmp, stmp); }
                zeta[i] = ztmp;
                surf[i] = stmp;
            }
            { ios::acstream fp(heur_name); fp << "\n"; }

            {
                ios::acstream fp(sfit_name);
                for(size_t i=1;i<=N;++i)
                {
                    fp("%.15g %.15g\n", zeta[i], surf0+slope*zeta[i]);
                }
                fp << "\n";
            }

            {
                ios::acstream fp(delta_name);
                const double delta_min = surf_min - (surf0+slope*zeta_max);

                for(size_t i=1;i<=N;++i)
                {
                    fp("%.15g %.15g\n", zeta[i]/zeta_max, (surf[i]-(surf0+slope*zeta[i]))/delta_min );
                }
                fp << "\n";
            }

#if 1
            aorg[1] = surf0;
            used[1] = false;

            aorg[3] = 0;
            used[3] = true;


            if(samples.fit_with(poly, aorg, used, aerr))
            {
                GLS<double>::display(std::cerr,aorg, aerr);
                ios::acstream fp(sfit_name);
                for(size_t i=1;i<=N;++i)
                {
                    fp("%.15g %.15g\n", zeta[i], sfit[i]);
                }
                fp << "\n";
            }
            else
            {
                throw exception("unexpected fit failure...");
            }

            {
                ios::acstream fp(coef_name);
                fp("%.15g", double(t_deg));
                for(size_t i=1;i<=m;++i)
                {
                    fp(" %.15g",aorg[i]);
                }
                fp("\n");
            }
#endif


        }

    }

    return 0;


#if 0
    //__________________________________________________________________________
    //
    // find zeta_max
    //__________________________________________________________________________
    const size_t N=21;
    vector<double> zeta(N);
    vector<double> surf(N);
    vector<double> sfit(N);

    GLS<double>::Function      poly = _GLS::Create<double,_GLS::Polynomial>();
    GLS<double>::Samples       samples(1);
    (void) samples.append(zeta, surf, sfit);

    const size_t   m=2;
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
        //const double zeta_min  = -0.1;
        const double zeta_min = -zeta_half;

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


        used[1] = false;
        aorg[1] = numeric<double>::pi * Square( sin(B.find_alpha(theta,0.0)) );

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
#endif
    
}
YOCTO_PROGRAM_END()
