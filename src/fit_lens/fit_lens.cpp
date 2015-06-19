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
    const double   alpha_cut;
    const double   Xpos;
    const double   Ypos;

    explicit Lens(const string &filename) :
    X(),
    Y(),
    F(),
    alpha(),
    rho(),
    N(0),
    alpha_cut(0),
    Xpos(0),
    Ypos(0)
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
        }
        (double &)alpha_cut = 0.5*(Fabs(alpha[1])+Fabs(alpha[N]));
        (double &)Xpos      = Xc;
        (double &)Ypos      = Yc;
    }

    void SavePolar() const
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

    void SaveRadii( ) const
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
        const double alpha_min = alpha[1];
        const double alpha_max = alpha[N];
        const int    NA        = 50;
        for(int i=0;i<=NA;++i)
        {
            const double a = alpha_min + i*(alpha_max-alpha_min)/NA;
            const double s = sin(a);
            const double c = cos(a);
            const double r = Rho(a,params);
            const double x = Xpos + r * s;
            const double y = Ypos - r * c;
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

    double RhoEx( const double angle, const array<double> &a)
    {
        if( Fabs(angle) <= alpha_cut )
        {
            return Rho(angle,a);
        }
        else
        {
            return Rho(alpha_cut,a);
        }
    }

    void Continuity(const array<double> &a)
    {

        const double rp = RhoPrime(alpha_cut, a);
        const double rc = Rho(alpha_cut,a);

        std::cerr << "alpha_cut=" << Rad2Deg(alpha_cut) << std::endl;
        std::cerr << "rho_c  =" << rc << std::endl;
        std::cerr << "drho_c =" << rp << std::endl;

        ios::ocstream fp("extra.dat",false);
        const int NA = 100;
        const double alpha_max = Deg2Rad(90.0);

        for(int i=-NA;i<=NA;++i)
        {
            const double angle = (i*alpha_max)/double(NA);
            const double r     = RhoEx(angle,a);
            fp("%g %g %g %g\n",Xpos+r*sin(angle),Ypos-r*cos(angle),angle,r);
        }

    }


    double GetOmega(const double angle, const array<double> &a)
    {
        if(Fabs(angle)<=alpha_cut)
        {
            const  double r = Rho(angle,a);
            const  double rp = RhoPrime(angle,a);
            return angle - asin( rp/Hypotenuse(r, rp) );
        }
        else
        {
            return angle;
        }
    }

    void SaveOmega(const array<double> &a)
    {
        ios::ocstream fp("omega.dat",false);
        const int    NA        = 100;
        const double alpha_max = Deg2Rad(90.0);

        for(int i=0;i<=NA;++i)
        {
            const double angle = (i*alpha_max)/double(NA);
            const double omega = GetOmega(angle,a);
            fp("%g %g\n", angle, omega);
        }

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
            lens.BuildWith(Xc,Yc);
            lens.SavePolar();
            lens.SaveRadii();

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
            lens.SaveProfile(aorg);
            lens.Continuity(aorg);
            lens.SaveOmega(aorg);
        }
        
        
        
    }
}
YOCTO_PROGRAM_END()


