#ifndef BRIDGE_INCLUDED
#define BRIDGE_INCLUDED 1

#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/ios/ostream.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/math/opt/minimize.hpp"

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

    explicit Bridge(const double delta_degrees,
                    const double ftol,
                    const double actrl_degrees,
                    const double sctrl);
    virtual ~Bridge() throw();
    
    const size_t                nvar;
    double                      mu;
    ode::driverCK<double>::type odeint; //!< use ftol as parameter
    vector<double>              pprev;
    vector<double>              param;
    bool                        flag;
    Equation                    eq;
    Callback                    cb;
    const double                delta; //!< search resolution in radians
    const double                angle_control;  //!< angle control in radians
    const double                shift_control; //!< max speed = tau max
    Function                    surface0;      //!< find_surface(theta,0.0)
    Function                    surface_of_zeta; //!< for a current_theta

    double profile( const double alpha, const double theta, const double zeta, ios::ostream *fp );
    double compute_zeta_max( const double theta );
    double find_alpha( const double theta, const double zeta );
    double find_surface( const double theta, const double zeta);

    void   compute_rates(array<double> &dY, const array<double> &Y) throw();
    void   compute_start(const double alpha, const double theta, const double zeta);
    double find_theta( const double alpha, const double zeta );


    double   v_center;
    double   mu2;
    double   current_theta;
    double   current_zeta;
    double   current_alpha;
    Function fn_of_alpha;
    Function fn_of_theta;
    optimize1D<double>::event check0;

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
    void __Eq( array<double> &dYdt, double, const array<double> &Y);
    void __Cb( array<double> &Y, double );
    double __profile_of_alpha(const double alpha);
    double __profile_of_theta(const double theta);
    double __surfac0_of_theta(const double theta);
    double __surface_of_zeta(const double zeta);
    
    bool   __check0(const triplet<double> &X, const triplet<double> &Y);
};

#endif
