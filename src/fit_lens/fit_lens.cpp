#include "yocto/math/io/data-set.hpp"
#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/opt/cgrad.hpp"

using namespace yocto;
using namespace math;

class Lens
{
public:
    vector<double> X;
    vector<double> Y;
    vector<double> alpha;
    vector<double> rho;
    const size_t   N;
    explicit Lens(const string &filename) :
    X(),
    Y(),
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
        alpha.make(N,0);
        rho.make(N,0);
    }

    void BuildWith(const double xc,const double R)
    {
        //assert(R>0);
        for(size_t i=1;i<=N;++i)
        {
            const double dx = X[i]-xc;
            const double dy = R - Y[i];
            rho[i]   = Hypotenuse(dx, dy);
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

    void SavePolar(const double xc, const double R) const
    {
        ios::ocstream fp( "polar.dat", false);
        for(size_t i=1;i<=N;++i)
        {
            const double a = alpha[i];
            const double s = sin(a);
            const double c = cos(a);
            const double x = xc+ rho[i] * s;
            const double y = R - rho[i] * c;
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
        return sqrt(res);
    }


    double H( const array<double> &q)
    {
        const double xc = q[1];
        const double R  = q[2];
        BuildWith(xc, R);
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
        q[1] = 0;
        q[2] = 1;
        vector<double> dq(nvar,1e-4);

        //lens.BuildWith(0,1);
        //lens.BuildWith(0,2);
        numeric<double>::scalar_field F(  &lens, & Lens::H );
        cgrad<double>::callback       cb( &lens, & Lens::CB);
        cgrad<double> CG;

        if(CG.run(F,q,dq,1e-4,&cb))
        {
            std::cerr << "SUCCESS" << std::endl;
            lens.SavePolar(0, q[1]);
        }



    }
}
YOCTO_PROGRAM_END()


