
#include "yocto/math/io/data-set.hpp"
#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/opt/cgrad.hpp"
#include "yocto/math/fit/glsf.hpp"
#include "yocto/code/ipower.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/math/core/lu.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/sort/quick.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/math/fcn/drvs.hpp"

using namespace yocto;
using namespace math;

static inline double G(double X,const double b)
{
    const double c   = 6.0-3*b;
    const double d   = 3*b-8.0;
    const double e   = 3.0-b;
    const double X2  = X*X;
    const double X4  = X2*X2;
    const double X6  = X2*X4;
    const double X8  = X4*X4;
    return b*X2+c*X4+d*X6+e*X8;
}

static inline double GS(const double alpha,const double b, const double beta)
{
    if(Fabs(alpha)>=beta)
    {
        return 1.0;
    }
    else
    {
        return G(alpha/beta,b);
    }
}

static inline double RhoFit(const double alpha, const array<double> &params)
{
    const double a    = params[1];
    const double b    = params[2];
    const double K    = params[3];
    const double beta = params[4];

    return a+ (K-a)*GS(alpha,b,beta);
}


class Lens
{
public:
    vector<double>            X;
    vector<double>            Y;
    vector<double>            F;
    vector<double>            alpha; // in degree
    vector<double>            rho;
    const size_t              N;
    const double              Xpos;
    const double              Ypos;
    mutable vector<double>    stk;
    vector<double>            params;
    numeric<double>::function Rho;
    numeric<double>::function RhoRad;
    derivative<double>        drvs;

    explicit Lens(const string &filename, const double ratio) :
    X(),
    Y(),
    F(),
    alpha(),
    rho(),
    N(0),
    Xpos(0),
    Ypos(0),
    stk(),
    params(),
    Rho(   this, & Lens::Rho__   ),
    RhoRad(this, & Lens::RhoRad__),
    drvs()
    {
        data_set<double> ds;
        ds.use(1, X);
        ds.use(2, Y);
        ios::icstream fp(filename);
        ds.load(fp);
        (size_t &)N = X.size();
        if(N<2)
            throw exception("not enough data points");
        F.    make(N,0);
        alpha.make(N,0);
        rho.  make(N,0);
        stk.  make(N,0);
        std::cerr << "ratio=" << ratio << std::endl;
        for(size_t i=1;i<=N;++i)
        {
            X[i] *= ratio;
            Y[i] *= ratio;
        }
    }



    virtual ~Lens() throw()
    {
    }

    inline void BuildWith(const double Xc,const double Yc) throw()
    {
        for(size_t i=1;i<=N;++i)
        {
            const double dx = X[i]-Xc;
            const double dy = Yc - Y[i];
            rho[i]    = Hypotenuse(dx, dy);
            alpha[i]  = Rad2Deg(2.0*atan(dx/(dy+rho[i])));
        }
        (double &)Xpos      = Xc;
        (double &)Ypos      = Yc;
    }

    inline double Energy() const
    {
        double res = 0;
        double ave = 0;
        for(size_t i=1;i<=N;++i)
        {
            ave += rho[i];
        }
        ave/=N;
        for(size_t i=1;i<=N;++i)
        {
            stk[i] = Square(rho[i]-ave);
        }
        quicksort(stk);
        for(size_t i=1;i<=N;++i)
        {
            res += stk[i];
        }
        return sqrt(res/N);
    }


    double H(const array<double> &q) throw()
    {
        const double Xc = q[1];
        const double Yc = q[2];
        BuildWith(Xc, Yc);
        return Energy();
    }

    bool CB(const array<double> &q) const
    {
        std::cerr << "q=" << q << std::endl;
        return true;
    }


    inline void SavePolar() const
    {
        ios::ocstream fp( "polar.dat", false);
        for(size_t i=1;i<=N;++i)
        {
            const double a = Deg2Rad(alpha[i]);
            const double s = sin(a);
            const double c = cos(a);
            const double x = Xpos + rho[i] * s;
            const double y = Ypos - rho[i] * c;
            fp("%g %g %g %g\n",x,y,alpha[i],rho[i]);
        }
    }

    inline void SaveRadii( ) const
    {
        ios::ocstream fp( "radii.dat", false);
        for(size_t i=1;i<=N;++i)
        {
            const double a = Deg2Rad(alpha[i]);
            const double s = sin(a);
            const double c = cos(a);
            const double x = Xpos + rho[i] * s;
            const double y = Ypos - rho[i] * c;
            fp("%g %g\n",x,y);
            fp("%g %g\n\n",Xpos,Ypos);
        }

    }

    inline double Rho__(double alpha)
    {
        return RhoFit(alpha, params);
    }

    inline double RhoRad__(double alpha)
    {
        return RhoFit(Rad2Deg(alpha),params);
    }


    inline double Omega(double alpha)
    {
        const double r0 = Rho(alpha);
        const double r1 = drvs(RhoRad,Deg2Rad(alpha),1e-4);
        return alpha - Rad2Deg( Asin(r1/Hypotenuse(r0, r1)) );
    }

    inline void SaveProfile( )
    {
        ios::ocstream fp( "profile.dat", false);
        for(double angle=-180;angle<=180;angle+=0.2)
        {
            const double a = Deg2Rad(angle);
            const double s = sin(a);
            const double c = cos(a);
            const double r = Rho(angle);
            const double x = Xpos + r * s;
            const double y = Ypos - r * c;
            fp("%g %g %g %g\n",x,y, angle, r);
            //fp("%g %g\n\n",Xpos,Ypos);
        }

    }

    inline void SaveOmega()
    {
        ios::ocstream fp( "omega.dat", false);
        for(double angle=0;angle<=180;angle+=0.1)
        {
            fp("%g %g\n", angle, Omega(angle));
        }
    }



private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Lens);

};


YOCTO_PROGRAM_START()
{
    //const char *program = vfs::get_base_name(argv[0]);
    if(argc<=3)
        throw exception("%s: need pixels mm lens_pixels",program);

    const double pixels   = strconv::to<double>(argv[1],"pixels");
    const double mm       = strconv::to<double>(argv[2],"mm");
    const string filename = argv[3];

    Lens lens(filename,mm/pixels);
    {
        const size_t nvar = 2;    //Xc,Yc
        vector<double> q(nvar,0);
        q[1] = lens.X[lens.N/2];
        q[2] = 1;
        vector<double> dq(nvar,1e-4);

        numeric<double>::scalar_field H(  &lens, & Lens::H );
        cgrad<double>::callback       cb( &lens, & Lens::CB);
        cgrad<double> CG;
        vector<bool> used(q.size(),true);

        if(!CG.run(H,q,used,dq,1e-4,&cb))
        {
            throw exception("Couldn't find the center of the lens...");
        }

        (void)H(q);
    }
    lens.SavePolar();
    lens.SaveRadii();

    {
        const size_t   nvar = 4;
        vector<double> &aorg = lens.params;
        aorg.make(nvar,0.0);
        vector<double> aerr(nvar,0);
        vector<bool>   used(nvar,true);

        const array<double> &rho   = lens.rho;
        const array<double> &alpha = lens.alpha;
        const size_t         N     = lens.N;
        //a
        aorg[1] = rho[N/2];

        //b
        aorg[2] = 4*(rho[N/4]-2* rho[N/2]+ rho[(3*N)/4]);

        //K
        aorg[3] = max_of(rho[1],rho[N]);

        //beta
        aorg[4] = 1.1 * max_of(alpha[1],alpha[N]);

        {
            ios::wcstream fp("sample.dat");
            for(double angle=1.2*lens.alpha.front();angle<=1.2*lens.alpha.back();angle+=0.1)
            {
                fp("%g %g\n", angle, RhoFit(angle,aorg));
            }
        }

        std::cerr << "Fitting..." << std::endl;
        GLS<double>::Samples samples;
        samples.append(lens.alpha,lens.rho,lens.F);
        samples.prepare(nvar);

        GLS<double>::Function FF( cfunctor2(RhoFit) );


        if(!samples.fit_with(FF, aorg, used, aerr, 0))
        {
            throw exception("Couldn't fit..");
        }
        //Fit.display(std::cerr,aorg,aerr);
        //std::cerr << "R" << nvar << "=" << samples.corr() << std::endl;
        //lens.SaveProfile(aorg);
        //lens.SaveOmega(aorg);
        lens.SaveProfile();
        lens.SaveOmega();
    }

    {
        ios::wcstream fp("lens.prm");
        for(size_t i=1;i<=lens.params.size();++i)
        {
            fp("%.8g\n", lens.params[i] );
        }
    }
    
    
}
YOCTO_PROGRAM_END()
