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
#include "yocto/math/types.hpp"

using namespace yocto;
using namespace math;

static const size_t   NVAR = 4;
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

    Bridge() : odeint(1e-7),
    R(100),
    kappa(1.0/R),
    V(NVAR,0)
    {
        odeint.start(NVAR);
    }


    ~Bridge() throw()
    {
    }

    void Equation( array_t &dVds, const double s, const array_t &V)
    {
        const double r     = V[1];
        const double z     = V[2];
        const double speed = Hypotenuse(V[3],V[4]);
        const double drds  = V[3]/speed;
        const double dzds  = V[4]/speed;
        dVds[1] = drds;
        dVds[2] = dzds;

        const double phi = kappa*kappa * z + dzds/r;
        dVds[3] =  dzds * phi;
        dVds[4] = -drds * phi;

    }

    void Legalize( array_t &V, double)
    {
        const double speed = Hypotenuse(V[3],V[4]);
        V[3] /= speed;
        V[4] /= speed;
    }

    inline void OutputBridge() const
    {
        const size_t  NB = 1024;
        ios::ocstream fp("bridge.dat",false);
        for(size_t i=0;i<=NB;++i)
        {
            const double angle = (numeric<double>::two_pi*i)/NB;
            fp("%g %g\n",R*sin(angle),R*(1.0-cos(angle)));
        }
    }

    void Build( const double theta, const double alpha )
    {
        std::cerr << "theta=" << theta << std::endl;
        std::cerr << "alpha=" << alpha << std::endl;

        const double aa = Deg2Rad(alpha);
        const double u0 = Deg2Rad((alpha+theta)-180.0);
        const double r0 = R * Sin(aa);
        const double z0 = R * (1.0-Cos(aa)); // + h

        V[1] = r0;
        V[2] = z0;
        V[3] = Cos(u0); // drds
        V[4] = Sin(u0); // dzds

        std::cerr << "dr=" << cos(u0) << std::endl;
        std::cerr << "dz=" << sin(u0) << std::endl;

        equation eq(this, & Bridge::Equation);
        callback cb(this, & Bridge::Legalize);
        
        double ds = 1e-2;
        double h  = 1e-8;
        ios::ocstream fp("profile.dat",false);
        fp("%g %g\n",V[1],V[2]);
        for(size_t i=0;;++i)
        {
            const double s0 = i*ds;
            const double s1 = (i+1)*ds;
            odeint(eq,V,s0,s1,h,&cb);
            fp("%g %g\n",V[1],V[2]);
            if(s1>=R*3)
                break;
        }



    }



private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
};

#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
#if 1
    if(argc<=2)
    {
        throw exception("need theta alpha");
    }
    const double theta = strconv::to<double>(argv[1],"theta");
    const double alpha = strconv::to<double>(argv[2],"alpha");
    Bridge B;
    B.Build(theta,alpha);
    B.OutputBridge();
#endif

#if 0
    if(argc<=1)
    {
        throw exception("need theta");
    }
    const double theta = strconv::to<double>(argv[1],"theta");
    Bridge B;
    const string fn = vformat("phase%g.dat",theta);
    ios::ocstream::overwrite(fn);
    for(int alpha=1;alpha<=179;++alpha)
    {
        B.Build(theta,alpha);
        ios::ocstream fp(fn,true);
        fp("%d %g %g\n",alpha,B.V[1],B.V[2]);
    }
#endif

}
YOCTO_PROGRAM_END()
