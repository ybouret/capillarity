#include "application.hpp"
#include "yocto/program.hpp"
#include "yocto/sys/wtime.hpp"

YOCTO_PROGRAM_START()
{
#define DISCARD_BRIDGE
#include "main-core.cpp"



#if 0
    const string filename = "zrange.dat";
    ios::ocstream::overwrite(filename);
    ios::ocstream::echo(filename,"0 0 0 180 180\n");
    for(double alpha_deg=1;alpha_deg<=30;++alpha_deg)
    {
        std::cerr << std::endl << "ALPHA=" << alpha_deg << std::endl;
        const double  alpha = Deg2Rad(alpha_deg);
        range_t       tr;
        const range_t zr    = B.find_zeta_range(alpha,tr);
        std::cerr << "theta: " << Rad2Deg(tr.vmin) << " to " << Rad2Deg(tr.vmax) << std::endl;
        ios::ocstream::echo(filename,"%g %g %g %g %g\n", alpha_deg, zr.vmin, zr.vmax, Rad2Deg(tr.vmin), Rad2Deg(tr.vmax));
    }
#endif


    Application  app(L);
    const string dirName  = argv[1];

    app.load_v2(dirName);
    wtime chrono;
    chrono.start();
    app.process();
    const double ell = chrono.query();
    std::cerr << "Done in " << ell << " seconds" << std::endl;
    
#if 0
    const string filename = L.Get<string>("file");
    app.load(filename);
    wtime chrono;
    chrono.start();
    app.process();
    const double ell = chrono.query();
    std::cerr << "Done in " << ell << " seconds" << std::endl;

    {
        string output = filename;
        string extension = "cbridge.dat";
        // changing extension if necessary
        extension = vformat("cbridge%g.dat",app.percent);
        vfs::change_extension(output, extension);
        std::cerr << "Saving into " << output << std::endl;


        {
            ios::wcstream fp(output);

            fp("#h A alpha theta h_corr t\r\n");
            for(size_t i=1;i<=app.h.size();++i)
            {
                fp("%g %g %g %g %g %g\r\n", app.h[i], app.A[i], Rad2Deg(app.alpha[i]), Rad2Deg(app.theta[i]), app.h_corr[i], app.t[i]);
            }
        }

#if 0
        const string outdir = vfs::get_file_dir(filename);
        const string outlog = outdir + "coeffs.txt";
        std::cerr << "Saving into " << outlog << std::endl;
        {
            ios::wcstream fp(outlog);
            fp("scan_rate_in_mm_per_second %g\r\n", app.main_rate);
            fp("eval_rate_in_mm_per_second %g\r\n", app.evap_rate);
            fp("coef_push                  %g\r\n", app.coef_push);
            //fp("coef_pull                  %g\r\n", app.coef_pull);
        }
#endif
    }
#endif
    
}
YOCTO_PROGRAM_END()
