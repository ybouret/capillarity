#include "lens.hpp"

Lens::Lens() :
drvs(),
scaling(LENS_DRVS_DEFAULT_SCALING),
radius(this, &Lens::Radius__)
{
}

Lens:: ~Lens() throw()
{
    
}


double Lens:: Radius__(const double alpha) const
{
    return GetRadius(alpha);
}