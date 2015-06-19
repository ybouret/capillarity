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
    
    //__________________________________________________________________________
    //
    //
    // The Equation
    //
    //__________________________________________________________________________
    inline void Evaluate(array<double>       &dudy,
                         const double         y,
                         const array<double> &u) throw()
    {
        const double r     = u[1];
        const double drdy  = u[2];
        const double speed = Hypotenuse(1.0,drdy);
        const double Ksq   = K*K;
        const double accel = ((Ksq)*y+1.0/r/speed)*speed*speed*speed;
        dudy[1] = drdy;
        dudy[2] = accel;
    }
    
    
    //__________________________________________________________________________
    //
    //
    // Integrator
    //
    //__________________________________________________________________________
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
    
    
    //__________________________________________________________________________
    //
    //
    // Saving the bridge
    //
    //__________________________________________________________________________
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
    
    
    //__________________________________________________________________________
    //
    //
    // Can We Find a Final Radius
    //
    //__________________________________________________________________________
    bool FinalRadius(const double beta,
                     const double theta,
                     const double alpha,
                     const bool   save=false)
    {
        
        //______________________________________________________________________
        //
        // initial conditions
        //______________________________________________________________________
        
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
    
    //__________________________________________________________________________
    //
    //
    // Try to find Alpha for a given beta and theta
    //
    //__________________________________________________________________________
    inline double FindAlpha(double beta,double theta)
    {
        //std::cerr << "FindAlpha(beta=" << beta << ",theta=" << theta << ")" << std::endl;
        
        if(theta<0||theta>=180)
        {
            return -1;
        }
        
        //______________________________________________________________________
        //
        // bracketing alpha
        //______________________________________________________________________
        double hi = 180 - theta; // assuming false
        double lo = hi/2;
        while( !FinalRadius(beta,theta,lo) )
        {
            hi  = lo;
            lo /= 2;
            if(lo<ATOL)
            {
                //std::cerr << "Impossible..." << std::endl;
                return -1;
            }
        }
        
        //______________________________________________________________________
        //
        // bissecting
        //______________________________________________________________________
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
        
        //______________________________________________________________________
        //
        // last good
        //______________________________________________________________________
        const double alpha = lo;
        return alpha;
    }
    
    //__________________________________________________________________________
    //
    //
    // Try to find theta for a given beta and alpha
    //
    //__________________________________________________________________________
    inline double FindTheta(const double beta, const double alpha)
    {
        if(alpha<=0)
            return -1;
        
        //______________________________________________________________________
        //
        // assuming a good value
        //______________________________________________________________________
        double theta_lo  = 0;
        if( !FinalRadius(beta, theta_lo, alpha))
        {
            return -1;
        }
        
        //______________________________________________________________________
        //
        // assuming bad value
        //______________________________________________________________________
        double theta_up = 180.0 - alpha;
        
        // bracket between good and bad
        while(theta_up-theta_lo>ATOL)
        {
            const double theta_mid = 0.5 * ( theta_up + theta_lo );
            if( FinalRadius(beta,theta_mid,alpha) )
            {
                theta_lo = theta_mid;
            }
            else
            {
                theta_up = theta_mid;
            }
            
        }
        const double theta = theta_lo;
        (void) FinalRadius(beta,theta,alpha);
        return theta;
    }
    
    //__________________________________________________________________________
    //
    //
    // Scan Alpha
    //
    //__________________________________________________________________________
    void ScanAlpha(double beta,double theta)
    {
        ios::ocstream fp(vformat("alpha_ofK%g_beta%g_theta%g.dat",K,beta,theta),false);
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
    
    //__________________________________________________________________________
    //
    //
    // Scan Theta
    //
    //__________________________________________________________________________
    void ScanTheta(double beta, double alpha)
    {
        ios::ocstream fp( vformat("theta_ofK%g_beta%g_alpha%g.dat",K,beta,alpha),false);
        for(int theta=0;theta<=179;++theta)
        {
            double value = 0;
            if(FinalRadius(beta,theta,alpha))
            {
                value = 1.0/U[1];
            }
            fp("%g %g\n", double(theta), value);
        }
        
    }
    
    //__________________________________________________________________________
    //
    //
    // Finding BetaMax
    //
    //__________________________________________________________________________
    inline double FindBetaMax(const double theta)
    {
        std::cerr << "FindBetaMax(theta=" << theta << ")" << std::endl;
        // check that a zero height is valid
        double beta_lo = 0;
        if(FindAlpha(beta_lo,theta)<0)
        {
            return -1;
        }
        
        // find an invalid beta
        double d_beta  = 1.0/K;
        double beta_up = d_beta;
        while(FindAlpha(beta_up,theta)>=0)
        {
            beta_lo = beta_up;
            beta_up += d_beta;
        }
        
        // beta_lo: valid, beta_up: invalid
        while(beta_up-beta_lo>BTOL)
        {
            const double beta_mid = 0.5 * (beta_lo+beta_up);
            if(FindAlpha(beta_mid,theta)<0)
            {
                beta_up = beta_mid;
            }
            else
            {
                beta_lo = beta_mid;
            }
        }
        const double beta  = beta_lo;
        const double alpha = FindAlpha(beta,theta);
        (void) FinalRadius(beta,theta,alpha,theta<=5);
        return beta;
    }
    
    
    //__________________________________________________________________________
    //
    //
    // Simulate a curve
    //
    //__________________________________________________________________________
    void GenerateCurve(const double theta, vector<double> &Beta, vector<double> &Alpha)
    {
        Beta. free();
        Alpha.free();
        
        const double beta_max = FindBetaMax(theta);
        std::cerr << "beta_max=" << beta_max << std::endl;
        if(beta_max<0)
        {
            return;
        }
        const double d_beta   = beta_max / 50;
        for(double beta = 0; beta < beta_max; beta += d_beta)
        {
            std::cerr << "beta= " << beta << "          \r"; std::cerr.flush();
            const double alpha = FindAlpha(beta, theta);
            if(alpha>=0)
            {
                Beta.  push_back(beta);
                Alpha. push_back(alpha);
            }
            else
            {
                break;
            }
        }
        std::cerr << std::endl;
        
    }
    
    //__________________________________________________________________________
    //
    //
    // Inverse a Curve
    //
    //__________________________________________________________________________
    void InverseCurve( array<double> &Theta, const array<double> &Beta, const array<double> &Alpha)
    {
        const size_t n = Theta.size(); assert(Beta.size()==n); assert(Alpha.size()==n);
        for(size_t i=1;i<=n;++i)
        {
            Theta[i] = FindTheta(Beta[i], Alpha[i]);
            std::cerr << "Theta(beta=" << Beta[i] <<",alpha=" << Alpha[i] << ",K=" << K << ")=" << Theta[i] << std::endl;
        }
        
    }
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
};


class DataFile
{
public:
    explicit DataFile(const string &filename,
                      const double R) :
    Height(),
    Surface(),
    Count(0),
    beta(),
    alpha(),
    theta()
    {
        data_set<double> ds;
        ds.use(1, Height);
        ds.use(2, Surface);
        ios::icstream fp(filename);
        ds.load(fp);
        (size_t &)Count = Height.size();
        beta.make(Count,0);
        alpha.make(Count,0);
        theta.make(Count,0);
        std::cerr << "Count=" << Count << std::endl;
        
        const double area = numeric<double>::pi*R*R;
        std::cerr << "Normalizing Area=" << area << " mm^2" << std::endl;
        for(size_t j=Count;j>0;--j)
        {
            beta[j] = Height[j]/R;
            const double sa = Sqrt(clamp<double>(0,Surface[j]/area,1.0));
            alpha[j] = Rad2Deg(Asin(sa));
        }
        
    }
    
    
    
    
    virtual ~DataFile() throw()
    {
    }
    
    vector<double> Height;
    vector<double> Surface;
    const size_t   Count;
    vector<double> beta;
    vector<double> alpha;
    vector<double> theta;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(DataFile);
};


YOCTO_PROGRAM_START()
{
#if 0
    if(argc<=2)
    {
        throw exception("Need K theta");
    }
    
    const double K     = strconv::to<double>(argv[1],"K");
    const double theta = strconv::to<double>(argv[2],"theta");
    
    Bridge B;
    B.K = K;
    
    vector<double> Beta;
    vector<double> Alpha;
    B.GenerateCurve(theta, Beta, Alpha);
    const string fn = vformat("alpha_K%g_theta%g.dat",K,theta);
    {
        ios::ocstream fp(fn,false);
        for(size_t i=1;i<=Beta.size();++i)
        {
            fp("%g %g\n", Beta[i], Alpha[i]);
        }
    }
    std::cerr << "Inversing..." << std::endl;
    vector<double> Theta(Beta.size(),0);
    B.InverseCurve(Theta,Beta, Alpha);
#endif
    
#if 1
    if(argc<=2)
    {
        throw exception("need R[mm], 1/kappa [mm]");
    }
    const double R     = strconv::to<double>(argv[1],"R");
    const double kappa = 1.0/strconv::to<double>(argv[2],"1/kappa");
    const double K     = R*kappa;
    
    Bridge B; B.K = K;
    
    for(int i=3;i<argc;++i)
    {
        const string filename = argv[i];
        DataFile     datafile(filename,R);
        string bname = vfs::get_base_name(filename);
        bname += ".dat";
        {
            ios::ocstream fp(bname,false);
            fp << "#h S beta alpha theta\n";
        }

        {
            std::cerr << "Processing " << vfs::get_base_name(filename);
            for(size_t j=1;j<=datafile.Count;++j)
            {
                const double alpha = datafile.alpha[j];
                const double beta  = datafile.beta[j];
                datafile.theta[j]  = B.FindTheta(beta,alpha);
                std::cerr << "\th=" << datafile.Height[j] << " => theta=" << datafile.theta[j] << std::endl;
                ios::ocstream fp(bname,true);
                fp("%g %g %g %g %g\n",datafile.Height[j], datafile.Surface[j], datafile.beta[j], datafile.alpha[j], datafile.theta[j]);
            }
            
        }
    }
#endif
    
    
}
YOCTO_PROGRAM_END()
