#include "yocto/program.hpp"
#include "bridge.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/fit/lsf.hpp"

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
    if(argc<=2)
        throw exception("need radius clength");
    const double R = strconv::to<double>(argv[1],"R");
    const double C = strconv::to<double>(argv[2],"C");

    Lens::Pointer lens( new SphericalLens(R));
    Bridge        bridge(lens, C);

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
    
    
}
YOCTO_PROGRAM_END()
