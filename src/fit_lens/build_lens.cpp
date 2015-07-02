
#include "yocto/math/io/data-set.hpp"
#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/opt/cgrad.hpp"
#include "yocto/math/fit/lsf.hpp"
#include "yocto/code/ipower.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/math/kernel/crout.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/sort/quick.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/trigconv.hpp"

using namespace yocto;
using namespace math;

class Lens
{
public:
    vector<double>         X;
    vector<double>         Y;
    vector<double>         F;
    vector<double>         alpha; // in degree
    vector<double>         rho;
    const size_t           N;
    const double           Xpos;
    const double           Ypos;
    mutable vector<double> stk;

    explicit Lens(const string &filename, const double ratio) :
    X(),
    Y(),
    F(),
    alpha(),
    rho(),
    N(0),
    Xpos(0),
    Ypos(0),
    stk()
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


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Lens);

};


YOCTO_PROGRAM_START()
{
    const char *program = vfs::get_base_name(argv[0]);
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

        if(!CG.run(H,q,dq,1e-4,&cb))
        {
            throw exception("Couldn't find the center of the lens...");
        }
        
        (void)H(q);
    }
    lens.SavePolar();
    lens.SaveRadii();

}
YOCTO_PROGRAM_END()
