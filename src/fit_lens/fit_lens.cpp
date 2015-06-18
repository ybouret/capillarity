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
        if(N<2)
            throw exception("not enough data points");
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
    
    void SaveRadii(const double Xc, const double Yc) const
    {
        ios::ocstream fp( "radii.dat", false);
        for(size_t i=1;i<=N;++i)
        {
            const double a = alpha[i];
            const double s = sin(a);
            const double c = cos(a);
            const double x = Xc + rho[i] * s;
            const double y = Yc - rho[i] * c;
            fp("%g %g\n",x,y);
            fp("%g %g\n\n",Xc,Yc);
        }

    }
    
    void SaveProfile(const double Xc,
                     const double Yc,
                     const array<double> &params) const
    {
        ios::ocstream fp("profile.dat",false);
        const double alpha_min = alpha[1];
        const double alpha_max = alpha[N];
        const int    NA        = 50;
        for(int i=0;i<=NA;++i)
        {
            const double a = alpha_min + i*(alpha_max-alpha_min)/NA;
            const double s = sin(a);
            const double c = cos(a);
            const double r = Rho(a,params);
            const double x = Xc + r * s;
            const double y = Yc - r * c;
            fp("%g %g %g %g\n",x,y,Rad2Deg(a),Curvature(a,params));
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


    bool CB(const array<double> &q) const
    {
        std::cerr << "q=" << q << std::endl;
        return true;
    }

    virtual ~Lens() throw()
    {

    }

    double Rho(double angle, const array<double> &a ) const
    {
        assert(a.size()>0);
        double ans = a[1];
        for(size_t i=2;i<=a.size();++i)
        {
            ans += a[i] * ipower(angle,2*(i-1));
        }
        return ans;
    }

    double RhoPrime(double angle, const array<double> &a) const
    {
        double ans = 0;
        for(size_t i=2;i<=a.size();++i)
        {
            const size_t p = 2*(i-1);
            ans += a[i] * p * ipower(angle,p-1);
        }
        return ans;
    }
    
    double RhoSecond(double angle, const array<double> &a) const
    {
        double ans = 0;
        for(size_t i=2;i<=a.size();++i)
        {
            const size_t p = 2*(i-1);
            ans += a[i] * p * (p-1) * ipower(angle,p-2);
        }
        return ans;
    }

    double Curvature(double angle, const array<double> &a) const
    {
        const double r0 = Rho(angle,a);
        const double r1 = RhoPrime(angle,a);
        const double r2 = RhoSecond(angle,a);
        const double c_den = sqrt(r0*r0+r1*r1);
        const double c_num = r0*r0 + 2*r1*r1 - r0*r2;
        return c_num/c_den;
    }
    
    void Continuity(const array<double> &a)
    {
        const double ac = 0.5*(-alpha[1]+alpha[N]);
        const double rp = RhoPrime(ac, a);
        const double rc = Rho(ac,a);
        
        std::cerr << "alpha_c=" << Rad2Deg(ac) << std::endl;
        std::cerr << "rho_c  =" << rc << std::endl;
        std::cerr << "drho_c =" << rp << std::endl;
        
        matrix<double> M(3,3);
        vector<double> X(3,0.0);
        X[2]    = rp;
        
        M[1][1] = ipower(ac,2);
        M[1][2] = ipower(ac,4);
        M[1][3] = ipower(ac,6);
        
        M[2][1] = 2 * ipower(ac,1);
        M[2][2] = 4 * ipower(ac,3);
        M[2][3] = 6 * ipower(ac,5);

        const double af = numeric<double>::pi;
        M[3][1] = 2 * ipower(af,1);
        M[3][2] = 4 * ipower(af,3);
        M[3][3] = 6 * ipower(af,5);
        
        if( !crout<double>::build(M) )
        {
            throw exception("invalid continuity...");
        }
        crout<double>::solve(M, X);
        std::cerr << "X=" << X << std::endl;
        std::cerr << rc
        << "+(" << X[1] << ")*x**2"
        << "+(" << X[2] << ")*x**4"
        << "+(" << X[3] << ")*x**6"
        << std::endl;

        
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
            const double Xc = q[1];
            const double Yc = q[2];
            lens.SavePolar(Xc,Yc);
            lens.SaveRadii(Xc,Yc);
            
            LeastSquares<double>::Samples samples;
            samples.append(lens.alpha,lens.rho,lens.F);

            size_t         nf   = 3;
            vector<double> aorg(nf,0);
            vector<double> aerr(nf,0);
            vector<bool>   used(nf,true);
            LeastSquares<double>::Function FF( &lens, & Lens::Rho );
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
            lens.SaveProfile(Xc, Yc, aorg);
            //lens.Continuity(aorg);
        }



    }
}
YOCTO_PROGRAM_END()


