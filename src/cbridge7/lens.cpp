#include "lens.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/trigconv.hpp"

Lens::Lens() :
drvs(),
scaling(LENS_DRVS_DEFAULT_SCALING),
radius(this, &Lens::Radius__),
dradius(this, &Lens::dRadius__),
omega(this, &Lens::Omega__),
surface(this,&Lens::Surface__)
{
}

Lens:: ~Lens() throw()
{
    
}


double Lens:: Radius__(const double alpha) const
{
    return max_of<double>(0,GetRadius(alpha));
}


double Lens:: dRadius__(const double alpha)
{
    return drvs(radius,alpha,scaling);
}

double Lens:: Omega__(const double alpha)
{
    const double r0 = radius(alpha);
    const double r1 = dradius(alpha);
    return alpha - Rad2Deg(Asin(r1/Hypotenuse(r0, r1)));
}

double Lens:: Surface__(const double alpha) 
{
    return numeric<double>::pi * Square( radius(alpha) * Sin(alpha) );
}


SphericalLens:: SphericalLens(const double RR) :
rho(RR)
{
    if(rho<=0) throw exception("Spherical Lens: negative radius");
}

SphericalLens:: ~SphericalLens() throw() {}

double SphericalLens:: GetRadius(const double alpha) const
{
    return rho;
}