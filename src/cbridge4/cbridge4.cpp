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

static const size_t   NVAR = 5;
typedef array<double> array_t;
typedef ode::Field<double>::Equation equation;
typedef ode::Field<double>::Callback callback;

#define I_R  1
#define I_Z  2
#define I_U  3
#define I_DR 4
#define I_DZ 5

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

    void Equation( array_t &dVds, const double , const array_t &V)
    {
        const double r = V[I_R];
        const double z = V[I_Z];
        const double u = V[I_U];

        const double drds = Cos(u);
        const double dzds = Sin(u);
        const double phi  = kappa*kappa*z + dzds/r;
        const double duds = -phi;

        dVds[I_R]  = drds;
        dVds[I_Z]  = dzds;
        dVds[I_U]  = duds;
        dVds[I_DR] = -duds*dzds;
        dVds[I_DZ] =  duds*drds;
    }

    void Legalize( array_t &V, double)
    {

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

        V[I_R]  = r0;
        V[I_Z]  = z0;
        V[I_U]  = u0;
        V[I_DR] = Cos(u0); // drds
        V[I_DZ] = Sin(u0); // dzds

        std::cerr << "dr=" << V[I_DR] << std::endl;
        std::cerr << "dz=" << V[I_DZ] << std::endl;

        equation eq(this, & Bridge::Equation);
        //callback cb(this, & Bridge::Legalize);
        
        double ds = 1e-2;
        double h  = 1e-8;
        ios::ocstream fp("profile.dat",false);
        fp("%g %g 0 %g\n",V[I_R],V[I_Z],V[I_U]);
        for(size_t i=0;;++i)
        {
            const double s0 = i*ds;
            const double s1 = (i+1)*ds;
            odeint(eq,V,s0,s1,h,NULL);
            fp("%g %g %g %g\n",V[I_R],V[I_Z],s1,V[I_U]);
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
