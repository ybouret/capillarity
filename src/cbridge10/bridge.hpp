#ifndef BRIDGE_INCLUDED
#define BRIDGE_INCLUDED 1

#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/ios/ostream.hpp"
#include "yocto/math/trigconv.hpp"

using namespace yocto;
using namespace math;


#define BRIDGE_U 1
#define BRIDGE_V 2
#define BRIDGE_A 3
#define BRIDGE_N 3

typedef numeric<double>::function    Function;
typedef ode::Field<double>::Equation Equation;
typedef ode::Field<double>::Callback Callback;

class Bridge
{
public:

    explicit Bridge(const double delta_degrees, const double ftol);
    virtual ~Bridge() throw();
    
    const size_t                nvar;
    double                      mu;
    ode::driverCK<double>::type odeint;
    vector<double>              pprev;
    vector<double>              param;
    vector<double>              prate;
    bool                        flag;
    Equation                    eq;
    Callback                    cb;
    const double                delta;

    double profile( const double alpha, const double theta, const double zeta, ios::ostream *fp );

    double compute_zeta_max( const double theta );

    double find_alpha( const double theta, const double zeta );

    double goodness(const double u, const double v, const double phi) const throw();

    void   compute_rates(array<double> &dY, const array<double> &Y) throw();
    void   compute_start(const double alpha, const double theta, const double zeta);
    

private:
    double   v_center;
    double   mu2;
    double   current_theta;
    double   current_zeta;
    Function fn_of_alpha;

    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
    void __Eq( array<double> &dYdt, double, const array<double> &Y);
    void __Cb( array<double> &Y, double );
    double __profile_of_alpha(const double theta);

};

#endif
