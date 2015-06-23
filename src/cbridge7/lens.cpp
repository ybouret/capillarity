#include "lens.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/exception.hpp"

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
    return alpha - (Asin(r1/Hypotenuse(r0, r1)));
}

double Lens:: Surface__(const double alpha)
{
    return numeric<double>::pi * Square( radius(alpha) * Sin(alpha) );
}

V2D Lens:: profile(const double alpha) const
{
    const double r     = Radius__(alpha);
    return V2D(r*Sin(alpha),Radius__(0)-r*Cos(alpha));
}

#include "yocto/ios/ocstream.hpp"

void Lens:: Output(const double height) const
{
    ios::ocstream fp("lens.dat",false);
    const size_t NL = 200;
    for(size_t i=0;i<=NL;++i)
    {
        const double alpha = (i*numeric<double>::pi)/NL;
        const V2D    q     = profile(alpha);
        fp("%g %g\n", q.x, q.y+height);
    }
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

