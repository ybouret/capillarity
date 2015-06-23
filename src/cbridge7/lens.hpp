#ifndef LENS_INCLUDED
#define LENS_INCLUDED 1

#include "yocto/math/fcn/drvs.hpp"
#include "yocto/counted-object.hpp"
#include "yocto/ptr/arc.hpp"

using namespace yocto;
using namespace math;

#define LENS_DRVS_DEFAULT_SCALING 1e-4
typedef numeric<double>::function function_type;

//! warning all angle are in degrees
class Lens : public counted_object
{
public:
    typedef arc_ptr<Lens> Pointer;
    
    derivative<double> drvs;
    double             scaling; //!< angle scaling, set to LENS_DRVS_DEFAULT_SCALING
    function_type      radius;  //!< radius(alpha)
    function_type      dradius; //!< d_radius/d_alpha
    function_type      omega;   //!< omega(alpha)
    function_type      surface; //!< surface(alpha)
    virtual ~Lens() throw();


protected:
    explicit Lens();
    virtual double GetRadius(const double alpha) const = 0;

    double Radius__(const double alpha) const;
    double dRadius__(const double alpha);
    double Omega__(const double alpha);
    double Surface__(const double alpha);

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Lens);
};


class SphericalLens : public Lens
{
public:
    explicit SphericalLens(const double RR);
    virtual ~SphericalLens() throw();
    const double rho;

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(SphericalLens);
    virtual double GetRadius(const double alpha) const;
};

#endif
