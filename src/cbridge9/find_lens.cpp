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

using namespace yocto;
using namespace math;

namespace
{


    class LensShape
    {
    public:
        const double   px2mm;
        const size_t   N;
        vector<double> X;
        vector<double> Y;
        vector<double> rho;
        vector<double> alpha;
        vector<double> rho_f;
        double         xmiddle;
        double         scaling;
        double         Xc;
        double         Yc;
        numeric<double>::scalar_field energy;
        cgrad<double>::callback       cb;
        GLS<double>::Function         ppoly;

        inline LensShape(const string &filename, const double user_px2mm) :
        px2mm(user_px2mm),
        N(0),
        X(), Y(), rho(),
        alpha(),
        rho_f(),
        xmiddle(0),
        scaling(0),
        Xc(0),
        Yc(0),
        energy(this, & LensShape::H),
        cb(    this, & LensShape::CB),
        ppoly( this, & LensShape::ComputeLSF)
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

        inline ~LensShape() throw()
        {
        }

        void ComputePolar(const array<double> &a)
        {
            assert(a.size()>=2);
            Xc = a[1];
            Yc = a[2];
            for(size_t i=N;i>0;--i)
            {
                const double dx = X[i]-Xc;
                const double dy = Yc-Y[i];
                const double rr = Hypotenuse(dx,dy);
                rho[i]   = rr;
                alpha[i] = 2.0 * atan(dx/(rr+dy));
            }
        }

        void SavePolar(const array<double> &a) const
        {
            ios::wcstream fp("polar.dat");
            for(size_t i=1;i<=N;++i)
            {
                fp("%g %g %g\n", alpha[i], rho[i], computeR(alpha[i],a) );
            }
        }

        void SaveShape() const
        {
            ios::wcstream fp("shape.dat");
            for(size_t i=1;i<=N;++i)
            {
                //fp("%g %g\n",Xc,Yc);
                fp("%g %g\n",X[i],Y[i]);
                //fp("\n");
            }
        }

        void SaveFitted(const array<double> &a) const
        {
            ios::wcstream fp("fitted.dat");
            for(size_t i=1;i<=N;++i)
            {
                const double R = computeR(alpha[i], a);
                const double x = Xc + R * sin(alpha[i]);
                const double y = Yc - R * cos(alpha[i]);
                fp("%g %g\n",x,y);
            }

        }

        double computeR(const double alpha, const array<double> &a) const
        {
            assert(a.size()>=3);
            const size_t nvar = a.size();
            double RR = a[3]*scaling;
            for(size_t j=4;j<=nvar;++j)
            {
                const size_t p = 2+2*(j-4);
                RR += a[j] * ipower(alpha, p);
            }
            return RR;
        }


        double H(const array<double> &a)
        {
            assert(a.size()>=3);
            // compute rho/alpha
            ComputePolar(a);

            // deduce energy
            double HH = 0;
            for(size_t i=N;i>0;--i)
            {
                double rr  = computeR(alpha[i],a);
                double dr  = (rho[i] - rr)/scaling;

                HH += dr*dr;
            }
            return HH/N;
        }

        bool CB(const array<double> &a) const
        {
            std::cerr << "a=" << a << std::endl;
            return true;
        }

        double ComputeLSF(const double aa, const array<double> &aorg)
        {
            assert(aorg.size()>0);
            double RR = scaling * aorg[1];
            for(size_t i=2;i<=aorg.size();++i)
            {
                const size_t p = 2+2*(i-2);
                RR += aorg[i] * ipower(aa, p);
            }
            return RR;
        }

    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(LensShape);
    };

}

YOCTO_PROGRAM_START()
{
    if(argc<=3)
        throw exception("%s: need pixels mm lens_pixels",program);

    const double pixels   = strconv::to<double>(argv[1],"pixels");
    const double mm       = strconv::to<double>(argv[2],"mm");
    const string filename = argv[3];

    LensShape shape(filename,mm/pixels);

    vector<double> params(6);
    vector<double> dparam(params.size(),1e-4);
    double &Xc = params[1];
    double &Yc = params[2];
    double &RR = params[3];
    Xc = shape.xmiddle;
    Yc = shape.scaling;
    RR = 1;




    (void) shape.energy(params);
    std::cerr << "params=" << params << std::endl;
    {
        vector<double> aorg(params.size()-2);
        const size_t   nvar = aorg.size();
        vector<bool>   used(nvar,true);
        vector<double> aerr(nvar);
        GLS<double>::Samples samples;
        samples.append(shape.alpha, shape.rho, shape.rho_f);
        samples.prepare(nvar);
        if(!samples.fit_with(shape.ppoly, aorg, used, aerr))
        {
            throw exception("Couldn't guess initial polynomial...");
        }
        GLS<double>::display(std::cerr, aorg, aerr);
        for(size_t i=1;i<=nvar;++i)
        {
            params[i+2] = aorg[i];
        }
    }
    (void) shape.energy(params);
    shape.SavePolar(params);
    shape.SaveShape();
    shape.SaveFitted(params);
    //return 0;


    cgrad<double> CG;
    vector<bool> used(params.size(),true);
    CG.run(shape.energy,
           params,
           used,
           dparam,
           1e-5,
           //0
           & shape.cb
           );
    (void) shape.energy(params);
    shape.SavePolar(params);
    shape.SaveShape();
    shape.SaveFitted(params);
    std::cerr << std::endl;
    std::cerr << "params=" << params << std::endl;
}
YOCTO_PROGRAM_END()
