#include "yocto/math/core/lu.hpp"
#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/code/ipower.hpp"

using namespace yocto;
using namespace math;

typedef matrix<double> Matrix;
typedef vector<double> Vector;

static double Ym = -180;
static double hm =  510;

static double h1 = -1538;
static double Y1 = 88;

static double h3 =  2994;
static double Y3 =  28;
static double sigma1 = -0.17;
static double sigma3 =  0.09;

YOCTO_PROGRAM_START()
{
    Matrix M(4);
    const double X1 = (h1-hm);
    const double X3 = (h3-hm);

    {
        array<double> &R1 = M[1];
        for(size_t i=1;i<=4;++i) R1[i] = ipower(X1,i+1);
    }

    {
        array<double> &R2 = M[2];
        for(size_t i=1;i<=4;++i) R2[i] = ipower(X3,i+1);
    }

    {
        array<double> &R3 = M[3];
        for(size_t i=1;i<=4;++i) R3[i] = (i+1)*ipower(X1,i);
    }

    {
        array<double> &R4 = M[4];
        for(size_t i=1;i<=4;++i) R4[i] = (i+1)*ipower(X3,i);
    }

    Vector rhs(4);
    rhs[1] = Y1-Ym;
    rhs[2] = Y3-Ym;
    rhs[3] = sigma1;
    rhs[4] = sigma3;


    std::cerr << "M=" << M << std::endl;
    std::cerr << "rhs=" << rhs << std::endl;
    if(!LU<double>::build(M))
    {
        throw exception("Singular System!");
    }
    LU<double>::solve(M,rhs);
    std::cerr << "rhs=" << rhs << std::endl;
    std::cerr.flush();

    const string YM = vformat("(%g)",Ym);
    const string XX = vformat("(x-%g)",hm);
    const string AA = vformat("(%g)",rhs[1]);
    const string BB = vformat("(%g)",rhs[2]);
    const string CC = vformat("(%g)",rhs[3]);
    const string DD = vformat("(%g)",rhs[4]);

    std::cerr << YM;
    std::cerr << "+" << AA << "*" << XX << "*" << XX;
    std::cerr << "+" << BB << "*" << XX << "*" << XX << "*" << XX;
    std::cerr << "+" << CC << "*" << XX << "*" << XX << "*" << XX << "*" << XX;
    std::cerr << "+" << DD << "*" << XX << "*" << XX << "*" << XX << "*" << XX << "*" << XX;

    std::cerr << std::endl;


}
YOCTO_PROGRAM_END()


