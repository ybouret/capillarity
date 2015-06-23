#ifndef LENS_INCLUDED
#define LENS_INCLUDED 1

#include "yocto/math/fcn/drvs.hpp"

using namespace yocto;
using namespace math;

#define LENS_DRVS_DEFAULT_SCALING 1e-4
typedef numeric<double>::function function_type;

class Lens
{
public:
    derivative<double> drvs;
    double             scaling; //!< scaling, set to LENS_DRVS_DEFAULT_SCALING
    function_type      radius;
    
    virtual ~Lens() throw();


protected:
    explicit Lens();
    virtual double GetRadius(const double alpha) const = 0;

    double Radius__(const double alpha) const;

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Lens);
};

#endif
