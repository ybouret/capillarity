#include "lens.hpp"

LensExtend:: LensExtend(const double _beta,
                        const double _rho_beta,
                        const double _rho_beta_prime,
                        const double _R) throw():
beta(_beta),
pimb(numeric<double>::pi-beta),
rho_beta(_rho_beta),
rho_beta_prime(_rho_beta_prime),
R(_R),
U(0),
V(0)
{
    static const double pi = numeric<double>::pi;
    assert(beta<pi);
    const double rho1 = pimb*rho_beta_prime;
    const double rho0 = rho_beta;
    {
        const double UU = 3*R - 2*rho1 -3*rho0;
        (double &)U = UU;
    }

    {
        const double VV = -2*R + rho1 + 2*rho0;
        (double &)V = VV;
    }
}


LensExtend:: ~LensExtend() throw()
{
}


double LensExtend:: compute(double alpha) const throw()
{
    const double del = alpha - beta;
    const double X   = del/pimb;
    return rho_beta + del * rho_beta_prime + X*X*(U+X*V);
}