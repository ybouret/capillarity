#include "yocto/program.hpp"
#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace math;




class diffBridge
{
public:
    double kappa;

    inline  diffBridge() throw() : kappa(2) {}
    inline ~diffBridge() throw() {}

    inline void computeZ( array<double> &dYdz, double z, const array<double> &Y )
    {
        const double r    = Y[1];
        const double drdz = Y[2];
        const double speed = sqrt(1.0+drdz*drdz);

        dYdz[1] = drdz;
        dYdz[2] = speed*speed*speed*(-kappa*kappa*z+1.0/(r*speed));
    }

    inline void computeS( array<double> &dYds, double , const array<double> &Y )
    {
        const double r   = Y[1];
        const double z   = Y[2];
        const double phi = Y[3];

        const double C = cos(phi);
        const double S = sin(phi);

        dYds[1] = C; // dr/ds
        dYds[2] = S; // dz/ds
        dYds[3] = -(kappa*kappa*z+S/r); // dphi/ds

    }

    inline void computeQ( array<double> &dYds, double s, const array<double> &Y )
    {
        const double r    = Y[1];
        const double z    = Y[2];
        const double drds = Y[3];
        const double dzds = Y[4];

        const double speed    = Hypotenuse(drds, dzds);
        const double speed_sq = speed*speed;
        const double f        = speed*speed_sq*(kappa*kappa*z + dzds/r/speed);

        dYds[1] =  drds;
        dYds[2] =  dzds;
        dYds[3] =  dzds * f / speed_sq;
        dYds[4] = -drds * f / speed_sq;


    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(diffBridge);
};



static inline void save_bridge(const double h, const double R)
{
    ios::wcstream fp("bridge.dat");
    const size_t N  = 1000;
    for(size_t i=0;i<=N;++i)
    {
        const double a = Deg2Rad( (i*360.0)/N );
        fp("%g %g\n", R*sin(a), h + R*(1.0-cos(a)) );
    }
}

YOCTO_PROGRAM_START()
{
    ode::driverCK<double>::type odeint(1e-7);

    odeint.start(2);

    const double R     = 1;
    const double h     = 0.01;
    const double theta = Deg2Rad(50.0);
    const double alpha = Deg2Rad(30.0);
    const double z0    = h + R*(1.0-cos(alpha));
    const double r0    = R*sin(alpha);

    save_bridge(h,R);

    vector<double> Y(2);

    Y[1] = r0;
    Y[2] = cos(theta+alpha)/sin(theta+alpha);

    const double dz   = 1e-3;
    double       dzz  = dz/10;


    std::cerr << "generating profile from Z equation..." << std::endl;

    diffBridge bridge;
    ode::Field<double>::Equation eqZ( &bridge, &diffBridge::computeZ);
    try
    {
        ios::wcstream fp("zprofile.dat");
        fp("%g %g\n", r0, z0);

        for(size_t i=0;i<=400;++i)
        {
            const double zstart = z0 - i*dz;
            const double zfinal = z0 - (i+1)*dz;
            odeint(eqZ,Y,zstart,zfinal,dzz,NULL);
            fp("%g %g\n", Y[1], zfinal);
        }
    }
    catch(const exception &e)
    {
        std::cerr << "couldnt integrate zProfile" << std::endl;
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }

    std::cerr << "generation profile from S equation..." << std::endl;
    const double ds  = 1e-3;
    double       dss = ds/10;
    ode::Field<double>::Equation eqS( &bridge, &diffBridge::computeS);

    Y.make(3);
    odeint.start(3);
    const double phi0=theta+alpha-numeric<double>::pi;
    Y[1] = r0;
    Y[2] = z0;
    Y[3] = phi0;
    {
        ios::wcstream fp("sprofile.dat");
        fp("%g %g %g\n", r0, z0, Y[3]);
        for(size_t i=0;i<=800;++i)
        {
            const double start = i*ds;
            const double final = (i+1)*ds;
            odeint(eqS,Y,start,final,dss,NULL);
            fp("%g %g %g\n", Y[1], Y[2], Y[3]);
        }
    }

    std::cerr << "generation profile from Q equation" << std::endl;
    Y.make(4);
    odeint.start(4);
    Y[1] = r0;
    Y[2] = z0;
    Y[3] = cos(phi0);
    Y[4] = sin(phi0);
    ode::Field<double>::Equation eqQ( &bridge, &diffBridge::computeQ);

    {
        ios::wcstream fp("qprofile.dat");
        fp("%g %g\n", r0, z0);
        for(size_t i=0;i<=800;++i)
        {
            const double start = i*ds;
            const double final = (i+1)*ds;
            odeint(eqQ,Y,start,final,dss,NULL);
            fp("%g %g\n", Y[1], Y[2]);
        }
    }

    
}
YOCTO_PROGRAM_END()
