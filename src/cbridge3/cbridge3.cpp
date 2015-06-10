#include "yocto/sequence/vector.hpp"
#include "yocto/sequence/many-arrays.hpp"
#include "yocto/program.hpp"
#include "yocto/math/types.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/math/kernel/tao.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fcn/zfind.hpp"
#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/math/trigconv.hpp"

using namespace yocto;
using namespace math;

static const size_t   NVAR = 3;
typedef array<double> array_t;
typedef ode::Field<double>::Equation equation;
typedef ode::Field<double>::Callback callback;

#if 0
static double       ATOL = 1e-5;
static double       BTOL = 1e-4;
static double       dY   = 1e-4;
static double       AMIN = 0.001;
#endif


class Bridge
{
public:
    ode::driverCK<double>::type odeint;

    double R;
    double kappa;
    vector<double> V;

    Bridge() : odeint(1e-6),
    R(100),
    kappa(3),
    V(NVAR,0)
    {
        odeint.start(NVAR);
    }


    ~Bridge() throw()
    {
    }

    void Equation( array_t &dVds, const double s, const array_t &V)
    {
        const double r = V[1];
        const double z = V[2];
        const double u = V[3];

        const double cu = Cos(u);
        const double su = Sin(u);

        dVds[1] = cu;
        dVds[2] = su;

        dVds[3] = -(su/r - 2.0 * kappa * kappa * z);


    }

    inline void OutputBridge() const
    {
        const size_t  NB = 100;
        ios::ocstream fp("bridge.dat",false);
        for(size_t i=0;i<=NB;++i)
        {
            const double angle = (numeric<double>::two_pi*i)/NB;
            fp("%g %g\n",R*sin(angle),R*(1.0-cos(angle)));
        }
    }

    void Build( const double alpha, const double theta )
    {
        std::cerr << "alpha=" << alpha << std::endl;
        std::cerr << "theta=" << theta << std::endl;

        const double aa = Deg2Rad(alpha);
        const double u0 = Deg2Rad((alpha+theta)-180.0);
        const double r0 = R * Sin(aa);
        const double z0 = R * (1.0-Cos(aa));

        V[1] = r0;
        V[2] = z0;
        V[3] = u0;

        std::cerr << "dr=" << cos(u0) << std::endl;
        std::cerr << "dz=" << sin(u0) << std::endl;

        equation eq(this, & Bridge::Equation);

        double ds = 1e-3;
        double h  = 1e-8;
        ios::ocstream fp("profile.dat",false);
        fp("%g %g %g\n",V[1],V[2],V[3]);
        for(size_t i=0;;++i)
        {
            const double s0 = i*ds;
            const double s1 = (i+1)*ds;
            odeint(eq,V,s0,s1,h,NULL);
            fp("%g %g %g\n",V[1],V[2],V[3]);
            if(s1>=R)
                break;
        }



    }



private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
};

#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    if(argc<=2)
    {
        throw exception("need alpha theta");
    }
    const double alpha = strconv::to<double>(argv[1],"alpha");
    const double theta = strconv::to<double>(argv[2],"theta");
    Bridge B;
    B.Build(alpha,theta);
    B.OutputBridge();
}
YOCTO_PROGRAM_END()
