#ifndef BRIDGE_INCLUDED
#define BRIDGE_INCLUDED 1

#include "yocto/lua/lua-config.hpp"
#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/math/opt/minimize.hpp"
#include "yocto/math/core/tao.hpp"
#include "yocto/math/round.hpp"
#include "yocto/container/tuple.hpp"

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
typedef vector<double>               Vector;
typedef triplet<double>              Triplet;




class Bridge
{
public:


    virtual ~Bridge() throw();
    explicit Bridge( lua_State *L );

    const size_t   nvar;     //!< BRIDGE_N
    DEsolver       odeint;   //!< adaptive solve, initialized with L->ftol (should be 1e-4-1e-7)
    Equation       profEq;   //!< differential equation
    bool           status;   //!< no error during integration
    const double   R0;       //!< radius
    const double   lambda;   //!< capillary length
    const double   mu2;            //!< mu^2...
    const double   mu;             //!< scaling
    Vector         pprev;          //!< starting point for integration
    Vector         param;          //!< current point for intergation
    double         center_v;       //!< center of lens pos, 1.0+zeta
    const double   resolution;     //!< angle resolution: search between resolution and PI-resolution
    const double   angle_control;  //!< angle control in radians (converted from L->angle_control)
    const double   shift_control;  //!< max speed = tau max
    Function       fn_of_alpha;    //!< ProfileOfAlpha
    Function       fn_of_theta;
    static size_t  curvature_coeff;


    //! initialize param and center_v
    void compute_start(const double alpha,
                       const double theta,
                       const double Xi);

    //! compute a profile...
    /**
     - return GetValue(v0,final_v) or exactly zero
     - the last param are always compute
     */
    double profile(const double alpha,
                   const double theta,
                   const double Xi,
                   ios::ostream *fp,
                   const bool    store_data=false,
                   const double  shift=0.0);

    //! debug, get dphi/dtau
    double reduced_rate( const array<double> &arr ) const throw();

    //! debug, get d^2phi/dtau^2
    double reduced_curv( const array<double> &arr ) const throw();


    //! this is the function to get a continuous level indicator
    /**
     - if (v0>=0) return Fabs(min_of(v0,v));
     - if (v0<0) return  Fabs(max_of(v0,v));
     */
    static double GetValue(const double v0, const double v) throw();

    //! save a drawing of the lens
    static void SaveLens(const string &filename, const double shift );

    static double CriticalThetaXi(const double theta)  throw();
    static double CriticalAlphaXi(const double alpha)  throw();

    //! finding alpha for a given theta and and Xi
    double find_alpha(const double theta, const double Xi, bool *is_flat=0);

    //! finding Xi_max(theta)
    double find_Xi_max(const double theta,double &alpha_min,double &zeta_max);

    //! finding theta for a given theta and Xi
    double find_theta(const double alpha, const double Xi, bool *id_flat=0);

    


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
    void   ProfileEq( array<double> &dYds, double, const array<double> &Y );
    double ProfileOfAlpha(const double alpha);
    double ProfileOfTheta(const double theta);
    double __alpha;
    double __theta;
    double __Xi;
    static double __find_top( Function &F, double p_lo, double p_up, const double res);
    static double __find_bot( Function &F, double p_lo, double p_up, const double res);

public:
    Vector heights;
    Vector radii;
    Vector slices;
    Vector volume;
    
    double   compute_user_sink(const double alpha, const double theta, const double Xi);

};


#endif

