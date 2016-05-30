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

    explicit Bridge(const double ftol);
    virtual ~Bridge() throw();
    
    const size_t                nvar;
    double                      mu;
    ode::driverCK<double>::type odeint;
    vector<double>              pprev;
    vector<double>              param;
    bool                        flag;
    Equation                    eq;
    Callback                    cb;

    double profile( const double alpha, const double theta, const double zeta, ios::ostream *fp );


private:
    double v_center;
    double mu2;
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
    void __Eq( array<double> &dYdt, double, const array<double> &Y);
    void __Cb( array<double> &Y, double );
};

#endif
