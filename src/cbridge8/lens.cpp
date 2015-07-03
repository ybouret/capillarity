#include "lens.hpp"
#include "yocto/math/trigconv.hpp"

Lens:: ~Lens() throw() {}

Lens:: Lens() :
drvs(),
rho(    this, & Lens::rho_     ),
rho_rad(this, & Lens::rho_rad_ ),
omega(  this, & Lens:: omega_  ),
surface(this, & Lens::surface_ ),
negsurf(this, & Lens::negsurf_ ),
max_surf_angle(0),
max_surf_value(0)
{

}

double Lens:: rho_(double alpha) throw()
{
    return ComputeRho(alpha);
}

double Lens:: rho_rad_(double alpha_rad) throw()
{
    return rho_( Rad2Deg(alpha_rad) );
}

double  Lens:: drho(double alpha) throw()
{
    return drvs(rho_rad,Deg2Rad(alpha),1e-4);
}

double Lens:: omega_(double alpha) throw()
{
    const double r0 = rho_(alpha);
    const double r1 = drho(alpha);
    return alpha - Rad2Deg( Asin(r0/Hypotenuse(r0,r1) ) );
}

double Lens:: surface_(double alpha) throw()
{
    const double r = rho_(alpha);
    const double s = Sin( Deg2Rad(alpha) );
    return numeric<double>::pi * Square(r*s);
}

double Lens:: negsurf_(double alpha) throw()
{
    return -surface_(alpha);
}

#include "yocto/math/opt/bracket.hpp"
#include "yocto/math/opt/minimize.hpp"
#include "yocto/exception.hpp"

void Lens:: initialize()
{
    triplet<double> A = { 0, 0, 180 };
    triplet<double> S = { 0, 0, 0   };
    if( !bracket<double>::inside(negsurf, A, S) )
    {
        throw exception("Lens::initialize: can't bracket maximum");
    }
    minimize<double>(negsurf, A, S, 1e-4);
    (double &)max_surf_angle = A.b;
    (double &)max_surf_value = surface(max_surf_angle);
}



SphericalLens:: SphericalLens(const double usr_radius) : Lens(),
radius(usr_radius)
{

}

SphericalLens:: ~SphericalLens() throw()
{

}

double SphericalLens:: ComputeRho(const double alpha) throw()
{
    return radius;
}

Lens * SphericalLens:: clone() const
{
    return new SphericalLens(this->radius);
}

