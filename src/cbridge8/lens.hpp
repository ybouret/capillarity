#ifndef LENS_INCLUDED
#define LENS_INCLUDED 1

#include "yocto/math/fcn/drvs.hpp"

using namespace yocto;
using namespace math;

typedef numeric<double>::function function_type;

class Lens : public object
{
public:
    derivative<double> drvs;
    function_type      rho;
    function_type      rho_rad;
    function_type      omega;
    function_type      surface;
    function_type      negsurf; //!< to find max surface
    const double       max_surf_angle;
    const double       max_surf_value;
    virtual ~Lens() throw();


    virtual  double ComputeRho(const double alpha) throw() = 0; //!< alpha in degree
    virtual  Lens  *clone() const = 0;
    double   rho_(double alpha) throw();
    double   rho_rad_(double alpha_rad) throw();
    double   drho(double alpha) throw();
    double   omega_(double alpha) throw();
    double   surface_(double alpha) throw();
    double   negsurf_(double alpha) throw();

    void initialize();

protected:
    explicit Lens();

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Lens);

};

class SphericalLens : public Lens
{
public:
    explicit SphericalLens(const double usr_radius);
    virtual ~SphericalLens() throw();

    const double radius;
    virtual  double ComputeRho(const double alpha) throw();
    virtual  Lens  *clone() const;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(SphericalLens);
};

#endif

