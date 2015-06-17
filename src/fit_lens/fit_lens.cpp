#include "yocto/math/io/data-set.hpp"
#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/opt/cgrad.hpp"
#include "yocto/math/fit/lsf.hpp"
#include "yocto/code/ipower.hpp"

using namespace yocto;
using namespace math;

class Lens
{
public:
    vector<double> X;
    vector<double> Y;
    vector<double> F;
    vector<double> alpha;
    vector<double> rho;
    const size_t   N;
    explicit Lens(const string &filename) :
    X(),
    Y(),
    F(),
    alpha(),
    rho(),
    N(0)
    {
        data_set<double> ds;
        ds.use(1, X);
        ds.use(2, Y);
        ios::icstream fp(filename);
        ds.load(fp);
        (size_t &)N = X.size();
        F.make(N,0);
        alpha.make(N,0);
        rho.make(N,0);
    }

    void BuildWith(const double Xc,const double Yc)
    {
        //assert(R>0);
        for(size_t i=1;i<=N;++i)
        {
            const double dx = X[i]-Xc;
            const double dy = Yc - Y[i];
            rho[i]    = Hypotenuse(dx, dy);
            alpha[i] =  2 * atan(dx/(dy+rho[i]));
#if 0
            ios::ocstream fp( vformat("polar%g.dat",R),false);
            for(size_t i=1;i<=N;++i)
            {
                const double a = alpha[i];
                const double s = sin(a);
                const double c = cos(a);
                const double x = xc+rho[i] * s;
                const double y = R - rho[i] * c;
                fp("%g %g %g %g\n",x,y,a,rho[i]);
            }
#endif
        }
    }

    void SavePolar(const double Xc, const double Yc) const
    {
        ios::ocstream fp( "polar.dat", false);
        for(size_t i=1;i<=N;++i)
        {
            const double a = alpha[i];
            const double s = sin(a);
            const double c = cos(a);
            const double x = Xc + rho[i] * s;
            const double y = Yc - rho[i] * c;
            fp("%g %g %g %g\n",x,y,a,rho[i]);
        }

    }

    inline double Energy()
    {
        double res = 0;
        double ave = 0;
        for(size_t i=1;i<=N;++i)
        {
            ave += rho[i];
            //res += Square(rho[i]);
        }
        ave/=N;
        for(size_t i=1;i<=N;++i)
        {
            res += Square(rho[i]-ave);
        }
        return sqrt(res/N);
    }


    double H( const array<double> &q)
    {
        const double Xc = q[1];
        const double Yc = q[2];
        BuildWith(Xc, Yc);
        return Energy();
    }


    bool CB(const array<double> &q)
    {
        std::cerr << "q=" << q << std::endl;
        return true;
    }

    virtual ~Lens() throw()
    {

    }

    double Eval(double angle, const array<double> &a )
    {
        assert(a.size()>0);
        double ans = a[1];
        for(size_t i=2;i<=a.size();++i)
        {
            ans += a[i] * ipower(angle,2*(i-1));
        }
        return ans;
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Lens);
};

YOCTO_PROGRAM_START()
{
    for(int i=1;i<argc;++i)
    {
        const string filename = argv[i];
        Lens lens(filename);
        const size_t nvar = 2;
        vector<double> q(nvar,0);
        q[1] = lens.X[lens.N/2];
        q[2] = 1;
        vector<double> dq(nvar,1e-4);

        numeric<double>::scalar_field H(  &lens, & Lens::H );
        cgrad<double>::callback       cb( &lens, & Lens::CB);
        cgrad<double> CG;

        if(CG.run(H,q,dq,1e-4,&cb))
        {
            std::cerr << "SUCCESS" << std::endl;
            lens.SavePolar(q[1],q[2]);

            LeastSquares<double>::Samples samples;
            samples.append(lens.alpha,lens.rho,lens.F);

            size_t         nf   = 5;
            vector<double> aorg(nf,0);
            vector<double> aerr(nf,0);
            vector<bool>   used(nf,true);
            LeastSquares<double>::Function FF( &lens, & Lens::Eval );
            LeastSquares<double>           Fit;

            samples.prepare(nf);
            std::cerr << "Fitting..." << std::endl;

            if(Fit(samples, FF, aorg, used, aerr, 0))
            {
                for( size_t i=1; i <= nf; ++i )
                {
                    std::cerr << "a[" << i << "]=" << aorg[i] << " +/- " << aerr[i] << std::endl;
                }
                std::cerr << "R" << nf << "=" << samples.corr() << std::endl;
            }
            std::cerr << aorg[1];
            for(size_t i=2;i<=nf;++i)
            {
                std::cerr << "+(" << aorg[i] << ")*(x**" << (2*(i-1)) << ")";
            }
            std::cerr << std::endl;
        }



    }
}
YOCTO_PROGRAM_END()


