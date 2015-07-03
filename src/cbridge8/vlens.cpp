#include "vlens.hpp"
#include "yocto/exception.hpp"

VLens:: ~VLens() throw()
{
}

VLens:: VLens(const array<double> &usr_params ) :
Lens(),
params(nparam,0.0)
{
    if(usr_params.size()<nparam)
        throw exception("Missing Parameters for VLens");
    for(size_t i=1;i<=nparam;++i)
    {
        params[i] = usr_params[i];
    }
}

double VLens:: ComputeRho(const double alpha) throw()
{
    return RhoFit(alpha, params);
}

Lens * VLens:: clone() const
{
    return new VLens(this->params);
}