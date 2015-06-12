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

using namespace yocto;
using namespace math;

typedef array<double> array_t;
static const size_t NVAR      = 2;
static double       ATOL      = 1e-4;
static size_t       NUM_STEPS = 10000;

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

    Bridge() :
    K(1),
    arrays(7),
    U(  arrays.next_array() ),
    k1( arrays.next_array() ),
    k2( arrays.next_array() ),
    k3( arrays.next_array() ),
    k4( arrays.next_array() ),
    V(  arrays.next_array() ),
    U0( arrays.next_array() )
    {
        arrays.allocate(NVAR);
    }

    ~Bridge() throw()
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
                     const double alpha, bool save=false)
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


        bool success = true;
        auto_ptr<ios::ocstream> fp;
        if(save)
        {
            fp.reset( new ios::ocstream("profile.dat", false) );
            (*fp)("%g %g\n", r_c, y_c);
        }

        for(size_t i=ny;i>0;--i)
        {
            // forward
            const double y_ini = (y_c*i)/ny;
            const double y_end = (y_c*(i-1))/ny;
            RK4(y_ini,y_end);

            // test
            const double u = U[1];
            if( isnan(u) || isinf(u) || u <= 0)
            {
                //std::cerr << "invalid u=" << u << std::endl;
                success = false;
                break;
            }
            const double y       = y_end;
            const double dist_sq = Square(y-y_b) + Square(u);
            if(dist_sq<1)
            {
                //std::cerr << "inside bridge" << std::endl;
                success = false;
                break;
            }
            if(save&&u<4)
            {
                (*fp)("%g %g\n",u,y);
            }

        }
        return success;

    }



    void ScanAlpha(double beta,double theta)
    {
        ios::ocstream fp("inv.dat",false);
        for(int alpha=1;alpha<=179;++alpha)
        {
            double value = 0;
            if(FinalRadius(beta,theta,alpha))
            {
                value = 1.0/U[1];
            }
            fp("%g %g\n", double(alpha), value);
        }
    }

    inline double FindAlpha(double beta,double theta)
    {
        std::cerr << "FindAlpha(beta=" << beta << ",theta=" << theta << ")" << std::endl;

        if(theta<=0||theta>=180)
        {
            return -1;
        }

        // bracketing alpha
        double hi = 180 - theta; // assuming false
        double lo = hi/2;
        while( !FinalRadius(beta,theta,lo) )
        {
            hi  = lo;
            lo /= 2;
            if(lo<ATOL)
            {
                std::cerr << "Impossible..." << std::endl;
                return -1;
            }
        }

        //std::cerr << "BracketAlpha: " << lo << "->" << hi << std::endl;
        while(hi-lo>ATOL)
        {
            const double mid = (lo+hi)*0.5;
            if(FinalRadius(beta,theta,mid))
            {
                // good
                lo = mid;
            }
            else
            {
                // bad
                hi = mid;
            }

        }
        // last good
        const double alpha = lo;
        //FinalRadius(beta,theta,alpha);
        return alpha;
    }


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
};



YOCTO_PROGRAM_START()
{
    vfs &fs = local_fs::instance();
    fs.try_remove_file("profile.dat");
    
    if(argc<=2)
    {
        throw exception("need beta theta");
    }
    const double beta  = strconv::to<double>(argv[1],"beta");
    const double theta = strconv::to<double>(argv[2],"theta");
    
    Bridge B;
    B.OutputBridge(beta);
    B.ScanAlpha(beta, theta);
    const double alpha = B.FindAlpha(beta,theta);
    if(alpha>=0)
    {
    (void)B.FinalRadius(beta,theta,alpha,true);
    }



}
YOCTO_PROGRAM_END()
