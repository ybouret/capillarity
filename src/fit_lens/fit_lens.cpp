#include "yocto/math/io/data-set.hpp"
#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/types.hpp"

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
        assert(R>0);
        for(size_t i=1;i<=N;++i)
        {
            const double dx = X[i]-xc;
            const double dy = R - Y[i];
            rho[i]   = Hypotenuse(dx, dy);
            alpha[i] =  2 * atan(dx/(dy+rho[i]));
            ios::ocstream fp( vformat("polar%g.dat",R),false);
            for(size_t i=1;i<=N;++i)
            {
                const double a = alpha[i];
                const double s = sin(a);
                const double c = cos(a);
                const double x = rho[i] * s;
                const double y = R - rho[i] * c;
                fp("%g %g %g %g\n",x,y,a,rho[i]);
            }
        }
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
        lens.BuildWith(0,1);
        lens.BuildWith(0,2);
    }
}
YOCTO_PROGRAM_END()


