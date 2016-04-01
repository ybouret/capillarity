#ifndef LENS_INCLUDED
#define LENS_INCLUDED 1

#include "yocto/math/fcn/drvs.hpp"

using namespace yocto;
using namespace math;

class LensExtend
{
public:
    const double beta;
    const double rho_beta;
    const double rho_beta_prime;
    const double rho_beta_second;
    const double R;
    const double U;
    const double V;
    const double W;

    LensExtend(const double _beta,
               const double _rho_beta,
               const double _rho_beta_prime,
               const double _rho_beta_second,
               const double _R) throw();
    ~LensExtend() throw();

    double compute(double alpha) const throw();

    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(LensExtend);
};


#endif
