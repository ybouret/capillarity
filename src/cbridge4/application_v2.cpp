#include "application.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/lang/pattern/matching.hpp"
#include "yocto/math/fit/glsf-spec.hpp"

static inline void __load(const string &filename,
                          Vector &A,
                          Vector &h )
{
    std::cerr << "-- Loading from " << filename << std::endl;
    data_set<double> ds;
    ds.use(1,A);
    ds.use(2,h);
    ios::icstream fp(filename);
    ds.load(fp);
}

class DoublePoly
{
public:
    vector<double> Alo;
    vector<double> Aup;

    double         Xlo;
    double         Ylo;
    double         Slo;
    double         Xm;
    double         Ym;
    double         Xup;
    double         Yup;
    double         Sup;

    explicit DoublePoly() :
    Alo(4),
    Aup(4)
    {
    }

    double eval(double x)
    {
        if(x<=Xm)
        {
            return _GLS::Polynomial<double>::Eval(x,Alo);
        }
        else
        {
            return _GLS::Polynomial<double>::Eval(x,Aup);
        }
    }

    void compute()
    {
        matrix<double> Mlo(4);
        matrix<double> Mup(4);

        {
            {
                for(size_t i=1;i<=4;++i)
                {
                    Mup[1][i] = Mlo[1][i] = ipower(Xm,i-1);
                    Mlo[2][i] = ipower(Xlo,i-1);
                    Mup[2][i] = ipower(Xup,i-1);
                }
                Alo[1] = Ym;
                Alo[2] = Ylo;
                Aup[1] = Ym;
                Aup[2] = Yup;
            }

            {
                Mup[3][1] = Mlo[3][1] = 0;
                Mup[4][1] = Mlo[4][1] = 0;
                for(size_t i=2;i<=4;++i)
                {
                    Mup[3][i] = Mlo[3][i] = (i-1)*ipower(Xm,i-2);
                    Mlo[4][i] = (i-1)*ipower(Xlo,i-2);
                    Mup[4][i] = (i-1)*ipower(Xup,i-2);
                }
                Alo[3] = 0;
                Aup[3] = 0;
                Alo[4] = Slo;
                Aup[4] = Sup;
            }
        }

        if(! LU<double>::build(Mlo) ) throw exception("Singular Cubic@Lower");
        if(! LU<double>::build(Mup) ) throw exception("Singular Cubic@Upper");

        LU<double>::solve(Mlo,Alo);
        LU<double>::solve(Mup,Aup);



    }

    virtual ~DoublePoly() throw()
    {
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(DoublePoly);
};

void Application:: load_v2(const string &dirName)
{
    local_fs &fs = local_fs::instance();
    //vfs::as_directory((string &)dirName);

    string enfoncement;
    string plateau;
    string tirage;

    {
        Lang::Matching sm_enf = "enfoncement";
        Lang::Matching sm_pla = "plateau";
        Lang::Matching sm_tir = "tirage";
        auto_ptr<vfs::scanner> scan( fs.new_scanner(dirName) );
        for( const vfs::entry *ep = scan->next(); ep; ep = scan->next() )
        {
            if(ep->is_directory()) continue;
            if(sm_enf.partly_matches(ep->base_name))
                enfoncement = ep->path;
            if(sm_pla.partly_matches(ep->base_name))
                plateau     = ep->path;
            if(sm_tir.partly_matches(ep->base_name))
                tirage      = ep->path;
        }
    }
    if(enfoncement.size()<=0) throw exception("enfoncement not found in %s", *dirName);
    if(plateau.size()    <=0) throw exception("plateau     not found in %s", *dirName);
    if(tirage.size()<=0     ) throw exception("tirage      not found in %s", *dirName);



    Vector h1,A1,h2,A2,h3,A3;

    __load(enfoncement,A1,h1);
    __load(plateau,A2,h2);
    __load(tirage,A3,h3);

    const size_t n1 = h1.size();
    const size_t n2 = h2.size();
    const size_t n3 = h3.size();

    offset1 = 0;
    offset2 = offset1 + n1;
    offset3 = offset2 + n2;

    const size_t n  = n1+n2+n3;
    mgr.make_all(__SUBS|__CORE,n);

    // be sure to get data
    {
        size_t j=1;
        for(size_t i=1;i<=n1;++i,++j)
        {
            h[j] = h1[i];
            A[j] = A1[i];
        }

        for(size_t i=1;i<=n2;++i,++j)
        {
            h[j] = h2[i];
            A[j] = A2[i];
        }

        for(size_t i=1;i<=n3;++i,++j)
        {
            h[j] = h3[i];
            A[j] = A3[i];
        }

        assert(1+n==j);

        for(size_t i=1;i<=n;++i)
        {
            h_corr[i] = h[i];
        }
    }



    // now change water level
    for(size_t i=1;i<=n1;++i)
    {
        h_corr[i] = h[i] + (shift+(p1/100.0)*h[i]);
    }

    //__________________________________________________________________________
    //
    //
    //
    //__________________________________________________________________________
    DoublePoly DP;
    DP.Xlo = h[n1];
    DP.Ylo = (h_corr[n1] - h[n1]); // water height in mm
    DP.Slo = p2/100.0;

    DP.Xm  = 0.417;
    DP.Ym  = -0.18; // water height in mm

    DP.Xup = 2994e-3;
    DP.Yup = 28e-3;
    DP.Sup = p3/100.0;

    DP.compute();

    {
        ios::wcstream fp("para.dat");

        for(double hh=DP.Xlo;hh<=DP.Xup;hh += 0.01)
        {
            fp("%g %g\n", hh, DP.eval(hh));
        }

    }

#if 0
    const double Y1 = (h_corr[n1] - h[n1]); // water height
    const double X1 = h[n1];

    const double Xm =  0.417;
    const double Ym = -180e-3;

    const double X3 = 2994e-3;
    const double Y3 = 28e-3;

    matrix<double> M(6);
    vector<double> A(6);

    std::cerr << "X1=" << X1 << std::endl;
    std::cerr << "Y1=" << Y1 << std::endl;

    std::cerr << "Xm=" << Xm << std::endl;
    std::cerr << "Ym=" << Ym << std::endl;

    std::cerr << "X3=" << X3 << std::endl;
    std::cerr << "Y3=" << Y3 << std::endl;

    // points
    {
        array<double> &row = M[1];
        for(size_t i=1;i<=6;++i) row[i] = ipower(X1,i-1);
        A[1] = Y1;
    }

    {
        array<double> &row = M[2];
        for(size_t i=1;i<=6;++i) row[i] = ipower(Xm,i-1);
        A[2] = Ym;
    }

    {
        array<double> &row = M[3];
        for(size_t i=1;i<=6;++i) row[i] = ipower(X3,i-1);
        A[3] = Y3;
    }

    //derivatives
    {
        array<double> &row = M[4];
        row[1]=0;
        for(size_t i=2;i<=6;++i) row[i] = (i-1)*ipower(X1,i-2);
        A[4] = p2/100.0;
    }

    {
        array<double> &row = M[5];
        row[1]=0;
        for(size_t i=2;i<=6;++i) row[i] = (i-1)*ipower(Xm,i-2);
        A[5] = 0;
    }

    {
        array<double> &row = M[6];
        row[1]=0;
        for(size_t i=2;i<=6;++i) row[i] = (i-1)*ipower(X3,i-2);
        A[6] = p3/100.0;
    }


    std::cerr << "M=" << M << std::endl;
    std::cerr << "A=" << A << std::endl;

    if( !LU<double>::build(M) )
    {
        throw exception("Singular Water Height Extrapolation!");
    }

    LU<double>::solve(M,A);
    std::cerr << "A=" << A << std::endl;
    {
        ios::wcstream fp("para.dat");

        for(double hh=X1;hh<=X3;hh += 0.01)
        {
            fp("%g %g\n", hh, _GLS::Polynomial<double>::Eval(hh,A) );
        }

    }
#endif

#if 0
    const double shift = 0;
    //const double __h0  = h1[1];
    const double b1    = 1.0 + (p1/100.0);
    const double b2    = 1.0 + (p2/100.0);
    const double b3    = 1.0 + (p3/100.0);

    size_t j=0;
    for(size_t i=1;i<=n1;++i)
    {
        ++j;
        h[j]      = h1[i];
        //h_corr[j] = shift + b1*(h1[i] - __h0);
        h_corr[j] = shift + b1*(h1[i]);
        A[j]      = A1[i];
    }

    const double   H1 = h_corr[j];
    const double __h1 = h[j];
    for(size_t i=1;i<=n2;++i)
    {
        ++j;
        h[j]      = h2[i];
        h_corr[j] = H1 + b2*(h2[i]-__h1);
        A[j]      = A2[i];
    }

    const double   H2 = h_corr[j];
    const double __h2 = h[j];
    for(size_t i=1;i<=n3;++i)
    {
        ++j;
        h[j]      = h3[i];
        h_corr[j] = H2 + b3*(h3[i]-__h2);
        A[j]      = A3[i];
    }

    {
        ios::wcstream fp("sample.dat");
        for(size_t i=1;i<=n;++i)
        {
            fp("%g %g %g\n", h[i], h_corr[i], A[i]);
        }
    }
#endif

    pre_process();



}
