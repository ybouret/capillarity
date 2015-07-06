#include "yocto/program.hpp"
#include "bridge.hpp"
#include "vlens.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/fit/lsf.hpp"
#include "yocto/fs/local-fs.hpp"

#include "yocto/threading/simd.hpp"

using namespace threading;

static const double HRES = 0.01;
class Slope
{
public:
    const double S0;
    inline Slope(const double surf0) throw() : S0(surf0)
    {
    }
    inline ~Slope() throw()
    {
    }

    double Eval(double h, const array_t &a)
    {
        assert(a.size()>0);
        double v = S0 - a[1] * h;
        if(a.size()>1)
        {
            v -= a[2] * h*h;
        }
        return v;
    }

};

YOCTO_PROGRAM_START()
{
    // reading params
    if(argc<=2)
        throw exception("need lens.file clength");
    const string   lens_name = argv[1];
    auto_ptr<Lens> full_lens( VLens::ReadFrom(lens_name) );
    const double R = full_lens->rho(0);
    const double C = strconv::to<double>(argv[2],"C");
    auto_ptr<Lens> spherical( new SphericalLens(R) );
    
    // preparing output
    string root = vfs::get_base_name(lens_name);
    vfs::remove_extension(root);
    std::cerr << "root=" << root << std::endl;
    std::cerr << "R   =" << R    << std::endl;
    std::cerr << "C   =" << C    << std::endl;
    
    
    const string f_name = vformat("abacus_%s_C%g.dat", root.c_str(), C);
    const string s_name = vformat("abacus_%s_C%g_R%g.dat", root.c_str(), C, R);
    const string h_name = vformat("height_%s_C%g.dat", root.c_str(),C);
    
    ios::ocstream::overwrite(f_name);
    ios::ocstream::overwrite(s_name);
    ios::ocstream::overwrite(h_name);

    SIMD simd(2,0);
    simd[0].build<Bridge,const Lens &, double>(*full_lens,C);
    simd[1].build<Bridge,const Lens &, double>(*spherical,C);

    Bridge & f_bridge  = simd[0].as<Bridge>();
    Bridge & s_bridge  = simd[1].as<Bridge>();
    
    
    static const char   wheel[] = "-\\|/";
    static const size_t nwheel  = sizeof(wheel)/sizeof(wheel[0])-1;
    size_t              iwheel  = 0;
    
    context::kernel kHmax( cfunctor(Bridge::CallHmax) );
    
    bool skip=false;
    for(int theta = 30; theta < 180; theta += 5 )
    {
        std::cerr << "theta=" << theta << std::endl;
        
        // find max height to integrate
        f_bridge.param1 = s_bridge.param1 = theta;
        simd(kHmax);
        
        //const double f_hmax = f_bridge.FindHmax(theta);
        //const double s_hmax = s_bridge.FindHmax(theta);
        const double f_hmax = f_bridge.result;
        const double s_hmax = s_bridge.result;
        
        {
            ios::acstream fp(h_name);
            fp("%d %g %g\n",theta,f_hmax,s_hmax);
        }
        std::cerr << "\theight_max=" << f_hmax << ", " << s_hmax << std::endl;
        
        // draw the abacus
        const size_t NH = 2+size_t(Ceil(max_of(f_hmax,s_hmax)/HRES));
        if(skip)
        {
            {
                ios::acstream fp(f_name);
                fp << "\n";
            }
            {
                ios::acstream fp(f_name);
                fp << "\n";
            }
        }
        else
        {
            skip=true;
        }
        
        for(size_t i=0;i<=NH;++i)
        {
            (std::cerr << "[" << wheel[++iwheel%nwheel] << "]  \r").flush();
            const double f_height = clamp<double>(0,(i*f_hmax)/double(NH),f_hmax);
            const double s_height = clamp<double>(0,(i*s_hmax)/double(NH),s_hmax);

        }
        std::cerr << std::endl;

    }
    
    
    
#if 0
    auto_ptr<Lens> lens( new SphericalLens(R));
    Bridge         bridge(*lens, C);

    const string a_name = vformat("abacus_R%g_C%g.dat",R,C);
    const string h_name = vformat("height_R%g_C%g.dat",R,C);
    const string f_name = vformat("curves_R%g_C%g.dat",R,C);
    const string s_name = vformat("slopes_R%g_C%g.dat",R,C);
    ios::ocstream::overwrite(a_name);
    ios::ocstream::overwrite(h_name);
    ios::ocstream::overwrite(f_name);
    ios::ocstream::overwrite(s_name);

    bool skip = false;

    static const char   wheel[] = "-\\|/";
    static const size_t nwheel  = sizeof(wheel)/sizeof(wheel[0])-1;
    size_t              iwheel  = 0;

    for(int theta = 30; theta < 180; theta += 10 )
    {
        std::cerr << "theta=" << theta << std::endl;

        // find max height to integrate
        const double hmax = bridge.FindHmax(theta);
        {
            ios::acstream fp(h_name);
            fp("%d %g\n",theta,hmax);
        }
        std::cerr << "\theight_max=" << hmax << std::endl;

        // draw the abacus
        const size_t NH = 2+size_t(Ceil(hmax/HRES));
        if(skip)
        {
            {
                ios::acstream fp(a_name);
                fp << "\n";
            }
            {
                ios::acstream fp(f_name);
                fp << "\n";
            }
        }
        else
        {
            skip=true;
        }

        {
            ios::acstream fp(a_name);
            fp("#h S(theta=%d)\n", theta);
        }



        const size_t   NF = NH/2;
        vector<double> fit_h(NF,as_capacity);
        vector<double> fit_s(NF,as_capacity);
        vector<double> fit_v(NF,as_capacity);

        for(size_t i=0;i<=NH;++i)
        {
            const double height     = clamp<double>(0,(i*hmax)/double(NH),hmax);
            const double alpha_deg  = bridge.FindAlpha(height, theta);
            const double alpha      = Deg2Rad(alpha_deg);
            const double S          = lens->surface(alpha);
            {
                ios::acstream fp(a_name);
                fp("%g %g\n", height, S);
            }
            (std::cerr << "\t[" << wheel[++iwheel%nwheel] << "]  \r").flush();
            if(i<NF)
            {
                fit_h.push_back(height);
                fit_s.push_back(S);
                fit_v.push_back(0);
            }
        }
        assert(NF==fit_h.size());
        assert(NF==fit_s.size());
        assert(NF==fit_v.size());

        std::cerr << std::endl;
        LeastSquares<double>::Samples samples;
        samples.append(fit_h,fit_s,fit_v);

        const size_t   nf = 1;
        vector<double> aorg(nf,0);
        vector<double> aerr(nf,0);
        vector<bool>   used(nf,true);
        samples.prepare(nf);

        Slope slope(fit_s[1]);
        LeastSquares<double>::Function FF( &slope, & Slope::Eval );
        LeastSquares<double>           Fit;

        aorg[1] = (fit_s[1]-fit_s[NF])/(fit_h[NF]-fit_h[1]);

        if(!Fit(samples, FF, aorg, used, aerr, 0))
        {
            throw exception("unable to fit...");
        }

        std::cerr << "aorg=" << aorg << std::endl;
        std::cerr << "aerr=" << aerr << std::endl;
        std::cerr << "R" << nf << "=" << samples.corr() << std::endl;
        std::cerr << slope.S0 << "-(" << aorg[1] << ")*x";
        if(nf>1) {
            std::cerr << "-(" << aorg[2] << ")*x**2";
        }
        std::cerr << std::endl;
        {
            ios::acstream fp(f_name);
            fp("#h  Fit(theta=%d,%g-(%g)*x-(%g)x**2)\n", theta, slope.S0, aorg[1], (nf>1?aorg[2]:0.0) );
            for(size_t i=1;i<=NF;++i)
            {
                fp("%g %g\n", fit_h[i], fit_v[i]);
            }
        }

        {
            ios::acstream fp(s_name);
            fp("%d %g %g",theta,slope.S0,aorg[1]);
            if(nf>1)
            {
                fp(" %g", aorg[2]);
            }
            fp << "\n";
        }
        
    }
#endif
    
}
YOCTO_PROGRAM_END()
