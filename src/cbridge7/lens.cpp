#include "lens.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/opt/bracket.hpp"
#include "yocto/math/opt/minimize.hpp"
#include "yocto/math/trigconv.hpp"

Lens::Lens() :
drvs(),
scaling(LENS_DRVS_DEFAULT_SCALING),
radius(this, &Lens::Radius__),
dradius(this, &Lens::dRadius__),
omega(this, &Lens::Omega__),
surface(this,&Lens::Surface__),
negsurf(this,&Lens::NegSurf__),
max_surface(0),
max_alpha(0),
max_alpha_deg(0)
{


}

void Lens:: Initialize()
{
    triplet<double> a = { 0, 0, numeric<double>::pi };
    triplet<double> s = { 0,0,0 };
    if( !bracket<double>::inside(negsurf, a, s) )
    {
        throw exception("Invalid Lens Profile!");
    }
    minimize<double>(negsurf, a, s, 0);
    (double&)max_alpha     = a.b;
    (double&)max_alpha_deg = Rad2Deg(max_alpha);
    (double&)max_surface   = surface(max_alpha);
    std::cerr << "LensMaxAlpha   = " << max_alpha_deg << std::endl;
    std::cerr << "LensMaxSurface = " << max_surface   << std::endl;
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

double Lens:: NegSurf__(const double alpha)
{
    return -Surface__(alpha);
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

