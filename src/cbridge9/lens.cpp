#include "lens.hpp"
#include "yocto/code/ipower.hpp"
#include "yocto/math/core/tao.hpp"

Lens:: ~Lens() throw()
{}


Lens:: Lens(const double         user_beta,
            const array<double> &user_coef,
            Derivative          &drvs) :
beta(user_beta),
pimb(numeric<double>::pi-beta),
coef(user_coef.size()),
rho_beta(0),
rho_beta_prime(0),
rho_fitted(this, &Lens::compute_fitted)
{
    tao::set(coef,user_coef);
    (double&)rho_beta       = compute_fitted(beta);
    (double&)rho_beta_prime = drvs(rho_fitted,beta,1e-4);
}


double Lens:: compute_fitted(const double alpha)
{
    double       rho = 0;
    const double a2  = alpha*alpha;
    for(size_t i=coef.size(),p=i-1;i>0;--i,--p)
    {
        rho += coef[i] * ipower(a2,p);
    }
    return rho;
}