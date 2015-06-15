#include "yocto/sequence/vector.hpp"
#include "yocto/sequence/many-arrays.hpp"
#include "yocto/program.hpp"
#include "yocto/math/trigconv.hpp"

#include "yocto/code/utils.hpp"
#include "yocto/math/kernel/tao.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/ptr/auto.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/math/fcn/zfind.hpp"

using namespace yocto;
using namespace math;

typedef array<double> array_t;
static const size_t NVAR      = 2;
static double       ATOL      = 1e-5;  // tolerance on alpha
static double       BTOL      = 1e-5;  // tolerance on beta
static size_t       NUM_STEPS = 10000; // used #steps

class Bridge
{
public:

    double K;      //!< kappa * R
    many_arrays<double,memory::global> arrays;
    array_t &U;
    array_t &k1;
    array_t &k2;
    array_t &k3;
    array_t &k4;
    array_t &V;
    array_t &U0;
    mutable double last_umin;

    Bridge() :
    K(1),
    arrays(7),
    U(  arrays.next_array() ),
    k1( arrays.next_array() ),
    k2( arrays.next_array() ),
    k3( arrays.next_array() ),
    k4( arrays.next_array() ),
    V(  arrays.next_array() ),
    U0( arrays.next_array() ),
    last_umin(0)
    {
        arrays.allocate(NVAR);
    }

    virtual ~Bridge() throw()
    {
    }

    inline void Evaluate(array<double>       &dudy,
                         const double         y,
                         const array<double> &u) throw()
    {
        const double r     = u[1];
        const double drdy  = u[2];
        const double speed = Hypotenuse(1.0,drdy);
        const double accel = ((K*K*y)+1.0/r/speed)*speed*speed*speed;
        dudy[1] = drdy;
        dudy[2] = accel;
    }

    inline void RK4(const double y_ini,
                    const double y_end) throw()
    {
        const double h     = y_end - y_ini;
        const double half  = h*0.5;
        const double y_mid = y_ini + half;
        Evaluate(k1,y_ini,U);

        tao::setprobe(V, U, half, k1);
        Evaluate(k2,y_mid,V);

        tao::setprobe(V,U,half,k2);
        Evaluate(k3,y_mid,V);

        tao::setprobe(V,U,h,k3);
        Evaluate(k4,y_end,V);

        for(size_t i=U.size();i>0;--i)
        {
            U[i] += h * (k1[i]+k2[i]+k2[i]+k3[i]+k3[i]+k4[i]) / 6.0;
        }
    }


    inline void OutputBridge(const double beta) const
    {
        const size_t NB = 1024;
        ios::ocstream fp("bridge.dat",false);
        for(size_t i=0;i<=NB;++i)
        {
            const double angle = (numeric<double>::two_pi * i)/NB;
            fp("%g %g\n", Sin(angle), beta + (1.0-Cos(angle)));
        }
    }


    bool FinalRadius(const double beta,
                     const double theta,
                     const double alpha,
                     const bool   save=false)
    {

        //std::cerr << "beta  = " << beta  << std::endl;
        //std::cerr << "theta = " << theta << std::endl;
        //std::cerr << "alpha = " << alpha << std::endl;
        // initial conditions
        if(theta+alpha>=180)
        {
            return false;
        }


        const double aa  = Deg2Rad(alpha);
        const double r_c = Sin(aa);
        const double y_c = beta + (1.0-Cos(aa));
        const double rot = Deg2Rad(alpha+theta);
        const double u_c = Cos(rot)/Sin(rot);
        const size_t ny  = NUM_STEPS;
        const double y_b = beta+1.0;

        U[1] = r_c; // initial radius
        U[2] = u_c; // initial slope drdy

        last_umin    = U[1];
        bool success = true;
        auto_ptr<ios::ocstream> fp;
        if(save)
        {
            fp.reset( new ios::ocstream("profile.dat", false) );
            (*fp)("%g %g\n", r_c, y_c);
        }


        double reg[3] = { U[1],0,0 }; // curvature detector
        size_t num    = 1;
        for(size_t i=ny;i>0;--i)
        {
            //__________________________________________________________________
            //
            // forward
            //__________________________________________________________________
            const double y_ini = (y_c*i)/ny;
            const double y     = (y_c*(i-1))/ny;
            RK4(y_ini,y);

            //__________________________________________________________________
            //
            // test valid radius
            //__________________________________________________________________
            const double u = U[1];
            if( isnan(u) || isinf(u) || u <= 0)
            {
                //std::cerr << "invalid u=" << u << std::endl;
                success = false;
                break;
            }

            //__________________________________________________________________
            //
            // test valid position
            //__________________________________________________________________
            const double dist_sq = Square(y-y_b) + Square(u);
            if(dist_sq<1)
            {
                //std::cerr << "inside bridge" << std::endl;
                success = false;
                break;
            }

            //__________________________________________________________________
            //
            // test valid curvature
            //__________________________________________________________________
            if(num<3)
            {
                reg[num++] = U[1];
            }
            else
            {
                reg[0] = reg[1];
                reg[1] = reg[2];
                reg[2] = U[1];
                if((reg[1]+reg[1])>= (reg[0]+reg[2]) )
                {
                    success = false;
                    break;
                }
            }

            //__________________________________________________________________
            //
            // may save
            //__________________________________________________________________
            if(save&&u<4)
            {
                (*fp)("%g %g\n",u,y);
            }
            
            if(u<last_umin)
            {
                last_umin = u;
            }
        }
        return success;
        
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
};

YOCTO_PROGRAM_START()
{
    
}
YOCTO_PROGRAM_END()
