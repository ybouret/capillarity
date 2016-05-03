#include "yocto/program.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sort/quick.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/stat/descr.hpp"
#include "yocto/code/ipower.hpp"
#include "yocto/math/trigconv.hpp"

using namespace yocto;
using namespace math;

#define NVAR   3

#define I_XC   1
#define I_YC   2
#define I_R0   3

class Shape
{
public:

    vector<double> X;     //!< X coordinates of the lens profile
    vector<double> Y;     //!< Y coordinates of the lens profile
    vector<double> rho;   //!< radius
    vector<double> alpha; //!< alpha
    vector<double> alpha_sq; //
    vector<double> rho_v2;   //!< rho for alpha_sq;
    const size_t   N; //!< #points
    double         Xg;
    double         scaling;



    inline Shape(const string &filename, const double px2mm ) :
    X(),
    Y(),
    rho(),
    alpha(),
    alpha_sq(),
    rho_v2(),
    N(0),
    scaling(0)
    {
        std::cerr << "// loading data from '" << filename << "'" << std::endl;
        ios::icstream fp(filename);
        {
            data_set<double> ds;
            ds.use(1, X);
            ds.use(2, Y);
            ds.load(fp);
        }
        (size_t &)N = X.size();
        rho.make(N);
        alpha.make(N);
        alpha_sq.make(N);
        rho_v2.make(N);
        std::cerr << "// found " << N << " coordinates" << std::endl;
        for(size_t i=N;i>0;--i)
        {
            X[i] *= px2mm;
            Y[i] *= px2mm;
        }
        co_qsort(X,Y);
        {
            double Yg=0;
            for(size_t i=N;i>0;--i)
            {
                Xg += X[i];
                Yg += Y[i];
            }
            Xg/=N;
            Yg/=N;
            for(size_t i=N;i>0;--i)
            {
                scaling += Square(X[i]-Xg) + Square(Y[i]-Yg);
            }
            scaling=Sqrt(scaling/N);
        }
        std::cerr << "scaling=" << scaling << std::endl;

    }

    inline ~Shape() throw()
    {
    }

    void buildPolar(const double center_x,const double center_y)
    {
        for(size_t i=N;i>0;--i)
        {
            const double dx = X[i] - center_x;
            const double dy = center_y - Y[i];
            const double rr = Hypotenuse(dx, dy);
            const double at = atan(dx/(dy+rr));
            rho[i]      = rr;
            alpha[i]    = at+at;
            rho_v2[i]   = rr;
            alpha_sq[i] = alpha[i]*alpha[i];
        }
        co_qsort(alpha_sq,rho_v2);
    }

    //! return variance of estimated radius
    double center_energy(const array<double> &param)
    {
        assert(param.size()>=2);
        const double center_x = param[I_XC];
        const double center_y = param[I_YC];
        buildPolar(center_x, center_y);
        double ave=0;
        double sig=0;
        compute_average_and_stddev(ave,sig,rho);
        return sig/scaling;
    }

    void savePolar(const string &filename) const
    {
        ios::wcstream fp(filename);
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g %g %g %g\n", X[i], Y[i], alpha[i], rho[i], alpha[i] * alpha[i]);
        }
    }

    //! return function to minimize
    double profile_energy(const array<double> &param)
    {
        assert(param.size()>=NVAR);
        const double center_x = param[I_XC];
        const double center_y = param[I_YC];
        buildPolar(center_x, center_y);

        double H = 0;
        for(size_t i=N;i>0;--i)
        {
            const double rr = rho[i];
            const double rf = profile(alpha[i],param);
            H += Square(rr-rf);
        }
        H = Sqrt(H/N)/scaling;
        //H /= (N*scaling*scaling);
        return H;
    }

    bool profile_callback(const array<double> &param)
    {
        const double H = profile_energy(param);
        std::cerr << "H=" << H << std::endl;
        return true;
    }

    inline double profile(const double angle, const array<double> &param) const
    {
        double        ans  = param[I_R0];
        const  size_t nvar = param.size();
        const double  a2   = angle*angle;
        for(size_t i=NVAR+1,j=1;i<=nvar;++i,++j)
        {
            ans += param[i] * ipower(a2,j);
        }
        return ans;
    }

    void saveApprox(const string &filename, const array<double> &param) const
    {
        ios::wcstream fp(filename);
        assert(param.size()>=NVAR);
        const double center_x = param[I_XC];
        const double center_y = param[I_YC];

        for(size_t i=1;i<=N;++i)
        {
            const double rr = profile(alpha[i],param);
            const double xx = center_x + rr * sin(alpha[i]);
            const double yy = center_y - rr * cos(alpha[i]);
            fp("%g %g %g %g %g\n", alpha[i], rho[i], rr, xx, yy);
        }
    }

    void saveRadii(const string &filename, const array<double> &param) const
    {
        ios::wcstream fp(filename);
        assert(param.size()>=NVAR);
        const double center_x = param[I_XC];
        const double center_y = param[I_YC];
        for(size_t i=1;i<=N;++i)
        {
            const double rr = profile(alpha[i],param);
            const double xx = center_x + rr * sin(alpha[i]);
            const double yy = center_y - rr * cos(alpha[i]);
            fp("%g %g\n", center_x, center_y);
            fp("%g %g\n\n", xx, yy);
        }

    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Shape);
};


#include "yocto/math/opt/cgrad.hpp"

YOCTO_PROGRAM_START()
{
    if(argc<=3)
        throw exception("%s: need pixels mm lens_pixels",program);

    const double pixels   = strconv::to<double>(argv[1],"pixels");
    const double mm       = strconv::to<double>(argv[2],"mm");
    const string filename = argv[3];

    Shape shape(filename,mm/pixels);
    cgrad<double> CG;
    const double ftol = 1e-5;

    //__________________________________________________________________________
    //
    // guess center
    //__________________________________________________________________________
    const size_t ndof = 3;
    const size_t nvar = NVAR + ndof;
    vector<double> param(nvar);
    vector<bool>   param_used(nvar,false);
    vector<double> param_scal(nvar,1e-4);

    {
        numeric<double>::scalar_field F( &shape, & Shape::center_energy);
        param[I_XC] = shape.Xg;
        param[I_YC] = shape.scaling;
        param_used[1] = param_used[2] = true;
        if(!CG.run(F,param,param_used,param_scal, ftol))
        {
            throw exception("cannot estimate center!");
        }
        (void)F(param);
    }
    std::cerr << "center_approx=" << param << std::endl;
    shape.savePolar("shape0.dat");

    //__________________________________________________________________________
    //
    // guess other parameters
    //__________________________________________________________________________
    param[I_R0] = shape.rho_v2[1];
    param_used.make(nvar,true);
    {
        numeric<double>::scalar_field F( &shape, & Shape::profile_energy);
        cgrad<double>::callback CB(&shape, & Shape::profile_callback);

        param_used.make(nvar,true);
        //for(size_t i=NVAR+1;i<=nvar;++i) param_used[i]=false;
        if(!CG.run(F,param,param_used,param_scal, ftol, &CB))
        {
            throw exception("cannot estimate parameters!");
        }
        (void)F(param);

        std::cerr << "param_approx=" << param << std::endl;
        shape.saveApprox("shape1.dat",param);
    }
    const double beta = (-shape.alpha.front()+shape.alpha.back())/2;
    std::cerr << "beta=" << beta << " (" <<  Rad2Deg(beta) << " deg)" << std::endl;
    shape.saveRadii("radii.dat", param);
    
}
YOCTO_PROGRAM_END()
