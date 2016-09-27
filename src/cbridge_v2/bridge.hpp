#ifndef BRIDGE_INCLUDED
#define BRIDGE_INCLUDED 1

#include "yocto/lua/lua-config.hpp"
#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/ios/ostream.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/math/opt/minimize.hpp"
#include "yocto/math/core/tao.hpp"

//#include "yocto/threading/crew.hpp"


using namespace yocto;
using namespace math;

#define BRIDGE_U 1
#define BRIDGE_V 2
#define BRIDGE_A 3
#define BRIDGE_N 3


typedef numeric<double>::function    Function;
typedef ode::Field<double>::Equation Equation;
typedef ode::Field<double>::Callback Callback;
typedef ode::driverCK<double>::type  DEsolver;
typedef vector<double> Vector;

class Bridge
{
public:
    virtual ~Bridge() throw();
    explicit Bridge( lua_State *L );

    const size_t nvar;     //!< BRIDGE_N
    DEsolver     odeint;   //!< adaptive solve, initialized with L->ftol (should be 1e-4-1e-7)
    Equation     profEq;   //!< differential equation
    bool         status;   //!< no error during integration

    double         mu;         //!< scaling
    double         mu2;        //!< mu^2...
    Vector         pprev;      //!< starting point for integration
    Vector         param;      //!< current point for intergation
    double         center_v;       //!< center of lens pos, 1.0+zeta
    const double   angle_control;  //!< angle control in radians
    const double   shift_control;  //!< max speed = tau max

    static size_t  curvature_coeff;

    //! compute mu and mu2
    void compute_mu(const double R0, const double capillary_length);

    //! initialize param and center_v
    void compute_start(const double alpha, const double theta, const double zeta);

    //! compute a profile...
    double profile( const double alpha, const double theta, const double zeta, ios::ostream *fp );

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
    void ProfileEq( array<double> &dYds, double, const array<double> &Y );
    
};


#endif

