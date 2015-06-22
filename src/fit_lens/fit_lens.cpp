#include "yocto/math/io/data-set.hpp"
#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/opt/cgrad.hpp"
#include "yocto/math/fit/lsf.hpp"
#include "yocto/code/ipower.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/math/kernel/crout.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/sort/quick.hpp"

using namespace yocto;
using namespace math;

class Lens
{
public:
    vector<double>         X;
    vector<double>         Y;
    vector<double>         F;
    vector<double>         alpha;
    vector<double>         rho;
    const size_t           N;
    const double           Xpos;
    const double           Ypos;
    double                 mu;
    mutable vector<double> stk;

    explicit Lens(const string &filename, const double ratio) :
    X(),
    Y(),
    F(),
    alpha(),
    rho(),
    N(0),
    Xpos(0),
    Ypos(0),
    mu(10),
    stk()
    {
        data_set<double> ds;
        ds.use(1, X);
        ds.use(2, Y);
        ios::icstream fp(filename);
        ds.load(fp);
        (size_t &)N = X.size();
        if(N<2)
            throw exception("not enough data points");
        F.    make(N,0);
        alpha.make(N,0);
        rho.  make(N,0);
        stk.  make(N,0);
        std::cerr << "ratio=" << ratio << std::endl;
        for(size_t i=1;i<=N;++i)
        {
            X[i] *= ratio;
            Y[i] *= ratio;
        }
    }

    inline void BuildWith(const double Xc,const double Yc)
    {
        for(size_t i=1;i<=N;++i)
        {
            const double dx = X[i]-Xc;
            const double dy = Yc - Y[i];
            rho[i]    = Hypotenuse(dx, dy);
            alpha[i]  = 2 * atan(dx/(dy+rho[i]));
        }
        (double &)Xpos      = Xc;
        (double &)Ypos      = Yc;
    }

    inline void SavePolar() const
    {
        ios::ocstream fp( "polar.dat", false);
        for(size_t i=1;i<=N;++i)
        {
            const double a = alpha[i];
            const double s = sin(a);
            const double c = cos(a);
            const double x = Xpos + rho[i] * s;
            const double y = Ypos - rho[i] * c;
            fp("%g %g %g %g\n",x,y,a,rho[i]);
        }
    }

    inline void SaveRadii( ) const
    {
        ios::ocstream fp( "radii.dat", false);
        for(size_t i=1;i<=N;++i)
        {
            const double a = alpha[i];
            const double s = sin(a);
            const double c = cos(a);
            const double x = Xpos + rho[i] * s;
            const double y = Ypos - rho[i] * c;
            fp("%g %g\n",x,y);
            fp("%g %g\n\n",Xpos,Ypos);
        }

    }


    void SaveProfile(const array<double> &params) const
    {
        ios::ocstream fp("profile.dat",false);
        //const double S         = params[4];
        const double alpha_min = -3.14;
        const double alpha_max =  3.14;
        const int    NA        = 500;
        for(int i=0;i<=NA;++i)
        {
            const double a = alpha_min + i*(alpha_max-alpha_min)/NA;
            const double s = sin(a);
            const double c = cos(a);
            const double r = RhoEx(a,params);
            const double x = Xpos + r * s;
            const double y = Ypos - r * c;
            fp("%g %g\n",x,y);
        }
    }


    inline double Energy() const
    {
        double res = 0;
        double ave = 0;
        for(size_t i=1;i<=N;++i)
        {
            ave += rho[i];
        }
        ave/=N;
        for(size_t i=1;i<=N;++i)
        {
            stk[i] = Square(rho[i]-ave);
        }
        quicksort(stk);
        for(size_t i=1;i<=N;++i)
        {
            res += stk[i];
        }
        return sqrt(res/N);
    }


    double H(const array<double> &q)
    {
        const double Xc = q[1];
        const double Yc = q[2];
        BuildWith(Xc, Yc);
        return Energy();
    }


    bool CB(const array<double> &q) const
    {
        std::cerr << "q=" << q << std::endl;
        return true;
    }

    virtual ~Lens() throw()
    {

    }

    inline
    double G(const double x,const array<double> &Q) const
    {
        const double a=Q[1];
        const double b=Q[2];
        const double K=Q[3];
        return
        a
        + b              * ipower(x,2)
        + (6*K-3*b-6*a)  * ipower(x,4)
        + (-8*K+3*b+8*a) * ipower(x,6)
        + (3*K-b-3*a)    * ipower(x,8);
    }

    //! fit/value function
    double Rho(double alpha, const array<double> &Q) const
    {
        const double S = Q[4];
        return G(alpha*S,Q);
    }

    double RhoEx(double alpha, const array<double> &Q) const
    {
        const double K = Q[3];
        const double S = Q[4];
        const double alpha_max = 1.0/S;
        return Fabs(alpha) >= alpha_max ? K : Rho(alpha,Q);
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Lens);
};


#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    const char *program = vfs::get_base_name(argv[0]);
    if(argc<=3)
        throw exception("%s: need pixels mm lens_pixels",program);

    const double pixels   = strconv::to<double>(argv[1],"pixels");
    const double mm       = strconv::to<double>(argv[2],"mm");
    const string filename = argv[3];

    Lens lens(filename,mm/pixels);
    const size_t nvar = 2;    //Xc,Yc
    vector<double> q(nvar,0);
    q[1] = lens.X[lens.N/2];
    q[2] = 1;
    vector<double> dq(nvar,1e-4);

    numeric<double>::scalar_field H(  &lens, & Lens::H );
    cgrad<double>::callback       cb( &lens, & Lens::CB);
    cgrad<double> CG;

    if(CG.run(H,q,dq,1e-4,&cb))
    {
        std::cerr << "Polar Fit SUCCESS" << std::endl;
        const double Xc = q[1];
        const double Yc = q[2];
        lens.BuildWith(Xc,Yc);
        lens.SavePolar();
        lens.SaveRadii();

        LeastSquares<double>::Samples samples;
        samples.append(lens.alpha,lens.rho,lens.F);

        size_t         nf   = 4;
        vector<double> aorg(nf,0);
        vector<double> aerr(nf,0);
        vector<bool>   used(nf,true);
        LeastSquares<double>::Function FF( &lens, & Lens::Rho );
        LeastSquares<double>           Fit;

        samples.prepare(nf);

        aorg[1] = lens.rho[lens.N/2];
        aorg[2] = 0;
        aorg[3] = 2*aorg[1];
        aorg[4] = 1;



        std::cerr << "Fitting..." << std::endl;
        Fit.verbose = true;
        if(Fit(samples, FF, aorg, used, aerr, 0))
        {
            for( size_t i=1; i <= nf; ++i )
            {
                std::cerr << "a[" << i << "]=" << aorg[i] << " +/- " << aerr[i] << " (" << (100.0*aerr[i]/Fabs(aorg[i])) << "%)" << std::endl;
            }
            std::cerr << "R" << nf << "=" << samples.corr() << std::endl;
            lens.SaveProfile(aorg);
            //lens.SaveOmega(aorg);
        }
        else
        {
            throw exception("Couldn't fit...");
        }

    }
    else
    {
        throw exception("Couldn't Fit Polar Shape...");
    }

    
    
}
YOCTO_PROGRAM_END()


