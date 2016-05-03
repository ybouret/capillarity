#ifndef LENS_INCLUDED
#define LENS_INCLUDED 1

#include "yocto/math/fcn/drvs.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/trigconv.hpp"

using namespace yocto;
using namespace math;

typedef numeric<double>::function Function;
typedef derivative<double>        Derivative;

class Lens
{
public:
    const double   beta; //!< upper limit of fitted value
    const double   pimb; //!< pi - beta
    vector<double> coef; //!< polynomial coefficient of alpha<2>
    const double   rho_beta; //!< fitted@beta
    const double   rho_beta_prime; //!< derivative of fitted@beta
    Function       rho_fitted;
    explicit Lens(const double         user_beta,
                  const array<double> &user_coef,
                  Derivative          &drvs );

    virtual ~Lens() throw();

    double compute_fitted(const double alpha);


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Lens);
};

#endif
