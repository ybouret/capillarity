#include "lens.hpp"

LensExtend:: LensExtend(const double _beta,
                        const double _rho_beta,
                        const double _rho_beta_prime,
                        const double _rho_beta_second,
                        const double _R) throw():
beta(_beta),
rho_beta(_rho_beta),
rho_beta_prime(_rho_beta_prime),
rho_beta_second(_rho_beta_second),
R(_R),
U(0),
V(0),
W(0)
{
    static const double pi = numeric<double>::pi;
    assert(beta>0);
    assert(beta<pi);
    const double del  = pi-beta;
    const double del2 = del*del;
    const double del3 = del2*del;
    const double del4 = del2*del2;
    const double del5 = del2*del3;
    const double rho0 = rho_beta;
    const double rho1 = del * rho_beta_prime;
    const double rho2 = del2 * rho_beta_second;

    {
        const double UU   = (20*R - 3*rho2 - 12*rho1 - 20*rho0)/del3;
        (double &)U = UU;
    }

    {
        const double VV   = (-30*R + 3*rho2 + 16*rho1 + 30 * rho0)/del4;
        (double &)V = VV;
    }

    {
        const double WW = (12*R - rho2 - 6*rho1 -12*rho0)/del5;
        (double &)W = WW;
    }

}


LensExtend:: ~LensExtend() throw()
{
}


double LensExtend:: compute(double alpha) const throw()
{
    const double del = alpha - beta;
    const double del2 = del*del;
    const double del3 = del2*del;
    const double del4 = del2*del2;
    const double del5 = del2*del3;
    return rho_beta
    + del * rho_beta_prime
    + (del2*rho_beta_second+del3*U+del4*V+del5*W)*0.5;
}