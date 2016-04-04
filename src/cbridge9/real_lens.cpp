#include "yocto/math/opt/cgrad.hpp"
#include "yocto/program.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sort/quick.hpp"
#include "yocto/code/ipower.hpp"
#include "yocto/math/fit/glsf.hpp"
#include "yocto/math/trigconv.hpp"

using namespace yocto;
using namespace math;

#define NEXTRA   3

#define INDEX_R0 1
#define INDEX_Xc 2
#define INDEX_Yc 3

namespace
{

    class Shape
    {
    public:
        const double px2mm;
        const size_t ndof;
        const size_t nvar;
        const size_t N;
        vector<double> X;
        vector<double> Y;
        vector<double> rho;
        vector<double> alpha;
        vector<double> rho_f;
        double       xmiddle;
        double       scaling;
        numeric<double>::scalar_field energy;
        cgrad<double>::callback       cb;

        inline Shape(const string &filename,
                     const double user_px2mm,
                     const size_t user_ndof) :
        px2mm(user_px2mm),
        ndof(user_ndof),
        nvar(ndof+NEXTRA),
        N(0),
        X(),
        Y(),
        rho(),
        alpha(),
        xmiddle(0),
        scaling(0),
        energy(this, & Shape::ComputeH ),
        cb( this, & Shape::CB)
        {
            ios::icstream fp(filename);

            //__________________________________________________________________
            //
            // Loading Pixels Coordinates
            //__________________________________________________________________
            data_set<double> ds;
            ds.use(1, X);
            ds.use(2, Y);
            ds.load(fp);

            //__________________________________________________________________
            //
            // prepare real size data
            //__________________________________________________________________
            (size_t &)N = X.size();
            rho.make(N);
            alpha.make(N);
            rho_f.make(N);
            std::cerr << "#coord=" << N << std::endl;
            for(size_t i=N;i>0;--i)
            {
                X[i] *= px2mm;
                Y[i] *= px2mm;
            }
            co_qsort(X,Y);

            //__________________________________________________________________
            //
            // find the scaling for CG
            //__________________________________________________________________
            scaling = 0;
            xmiddle = 0;
            for(size_t i=1;i<N;++i)
            {
                xmiddle += X[i];
                for(size_t j=i+1;j<=N;++j)
                {
                    const double tmp = Hypotenuse(X[j]-X[i], Y[j]-Y[i]);
                    if(tmp>scaling) scaling=tmp;
                }
            }
            xmiddle/=N;
            std::cerr << "xmiddle=" << xmiddle << std::endl;
            std::cerr << "scaling=" << scaling << std::endl;

        }

        inline ~Shape() throw()
        {
        }

        void ComputePolar(const array<double> &params )
        {
            assert(params.size()>=nvar);
            //const double R0 = params[m+1];
            //const double K  = params[m+2];
            const double Xc = params[ndof+INDEX_Xc];
            const double Yc = params[ndof+INDEX_Yc];
            for(size_t i=N;i>0;--i)
            {
                const double dx = X[i] - Xc;
                const double dy = Yc   - Y[i];
                const double rr = Hypotenuse(dx,dy);
                rho[i]   = rr;
                alpha[i] = 2*atan(dx/(dy+rr));
            }
        }

        double ComputeRho(const double alpha, const array<double> &params) const
        {
            assert(params.size()>=nvar);
            const double R0 = params[ndof+INDEX_R0];
            //const double Xc = params[m+INDEX_Xc];
            //const double Yc = params[m+INDEX_Yc];
            return R0 + GetP( Square(alpha/numeric<double>::pi), params);
        }



        // get the reduced Polynomial
        double GetP(const double X, const array<double> &params) const
        {
            assert(params.size()>=ndof);
            const size_t M   = ndof+2;
            const size_t Mm1 = M-1;
            double U = 0;
            double V = 0;
            for(size_t i=1;i<=ndof;++i)
            {
                U -= params[i];
                V -= i * params[i];
            }
            const double CMm1 = M*U-V;
            const double CM   = V-Mm1*U;
            double ans = 0;
            for(size_t i=1;i<=ndof;++i)
            {
                ans += params[i]  * ipower(X,i);
            }
            ans += CMm1 * ipower(X,Mm1);
            ans += CM   * ipower(X,M);
            return ans;
        }

        double ComputeH(const array<double> &params) const
        {
            //const double R0 = params[ndof+INDEX_R0];
            //const double K  = params[ndof+INDEX_K];
            const double Xc = params[ndof+INDEX_Xc];
            const double Yc = params[ndof+INDEX_Yc];

            double H = 0;
            for(size_t i=N;i>0;--i)
            {
                const double dx = X[i] - Xc;
                const double dy = Yc - Y[i];
                const double rr = Hypotenuse(dx, dy);
                const double aa = 2*atan(dx/(rr+dy));
                const double rf = ComputeRho(aa,params);
                const double dr = (rr-rf)/scaling;
                H += dr*dr;
            }
            return H/N;
        }


        void SaveShape(const string &filename, const array<double> &params) const
        {
            ios::wcstream fp(filename);
            const double Xc = params[ndof+INDEX_Xc];
            const double Yc = params[ndof+INDEX_Yc];
            for(size_t i=1;i<=N;++i)
            {
                const double dx = X[i] - Xc;
                const double dy = Yc - Y[i];
                const double rr = Hypotenuse(dx, dy);
                const double aa = 2*atan(dx/(rr+dy));
                const double rf = ComputeRho(aa,params);
                const double x  = Xc + rf * sin(aa);
                const double y  = Yc - rf * cos(aa);
                fp("%g %g %g %g\n", X[i], Y[i], x, y);
            }
        }

        void SavePolar(const string &filename, const array<double> &params) const
        {
            ios::wcstream fp(filename);
            const double Xc = params[ndof+INDEX_Xc];
            const double Yc = params[ndof+INDEX_Yc];
            for(size_t i=1;i<=N;++i)
            {
                const double dx = X[i] - Xc;
                const double dy = Yc - Y[i];
                const double rr = Hypotenuse(dx, dy);
                const double aa = 2*atan(dx/(rr+dy));
                const double rf = ComputeRho(aa,params);
                //const double x  = Xc + rf * sin(aa);
                //const double y  = Yc - rf * cos(aa);
                fp("%g %g %g\n", aa, rr, rf );
            }

        }

        bool CB( const array<double> &params)
        {
            std::cerr << "params=" << params << std::endl;
            return true;
        }

        void SaveFull( const string &filename, const array<double> &params) const
        {
            ios::wcstream fp(filename);
            const double Xc = params[ndof+INDEX_Xc];
            const double Yc = params[ndof+INDEX_Yc];
            for(int ai=-180;ai<=180;++ai)
            {
                const double aa = Deg2Rad( double(ai) );
                const double rf = ComputeRho(aa,params);
                const double x  = Xc + rf * sin(aa);
                const double y  = Yc - rf * cos(aa);
                fp("%g %g %g %g\n",x,y, aa, rf);
            }
        }


    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(Shape);
    };

}


YOCTO_PROGRAM_START()
{
    if(argc<=3)
        throw exception("%s: need pixels mm lens_pixels",program);

    const double pixels   = strconv::to<double>(argv[1],"pixels");
    const double mm       = strconv::to<double>(argv[2],"mm");
    const string filename = argv[3];

    Shape shape(filename,mm/pixels,2);
    const size_t ndof = shape.ndof;
    const size_t nvar = ndof + NEXTRA;
    vector<double> params(nvar,0);
    vector<double> dparam(nvar,1e-4);

    params[ndof+INDEX_Xc] = shape.xmiddle;
    params[ndof+INDEX_Yc] = shape.scaling;
    params[ndof+INDEX_R0] = shape.scaling;

    vector<bool> used(nvar,false);

    // fist part
    used[ndof+INDEX_Xc] = true;
    used[ndof+INDEX_Yc] = true;
    used[ndof+INDEX_R0] = true;

    cgrad<double> CG;

    std::cerr << "average R..." << std::endl;
    if(!CG.run(shape.energy, params, used, dparam, 1e-5, & shape.cb))
    {
        throw exception("cannot find average radius!");
    }

    shape.SaveShape("shape0.dat",params);
    shape.SavePolar("polar0.dat",params);

    if(ndof>0)
    {
        std::cerr << "refining..." << std::endl;
        shape.ComputePolar(params);
        used.make(nvar, false);
        for(size_t i=1;i<=ndof;++i) used[i] = true;
        used[ndof+INDEX_R0] = true;

        GLS<double>::Samples samples;
        samples.append(shape.alpha, shape.rho, shape.rho_f);
        samples.prepare(nvar);

        GLS<double>::Function F( &shape, &Shape::ComputeRho);
        vector<double> aerr(nvar,0);

        if(!samples.fit_with(F,params, used, aerr))
        {
            throw exception("Cannot match initial DOF...");
        }

        GLS<double>::display(std::cerr, params, aerr);

        std::cerr << "final full fit..." << std::endl;
        used.make(nvar,true);
        if(!CG.run(shape.energy, params, used, dparam, 1e-5, & shape.cb))
        {
            throw exception("cannot find full radius!");
        }
        shape.SaveShape("shape1.dat",params);
        shape.SavePolar("polar1.dat",params);


    }

    shape.SaveFull("shape_full.dat", params);


}
YOCTO_PROGRAM_END()
