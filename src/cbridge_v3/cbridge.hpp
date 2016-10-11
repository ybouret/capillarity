#ifndef BRIDGE_INCLUDED
#define BRIDGE_INCLUDED 1

#include "yocto/lua/lua-config.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/ios/icstream.hpp"
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
#define BRIDGE_Q 4
#define BRIDGE_q 5
#define BRIDGE_N 5


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
    explicit Bridge( lua_State *L);

    const size_t   nvar;     //!< BRIDGE_N
    DEsolver       odeint;   //!< adaptive solve, initialized with L->ftol (should be 1e-4-1e-7)
    Equation       profEq;   //!< differential equation
    bool           status;   //!< no error during integration
    const double   R0;       //!< radius
    const double   A0;       //!< pi * R0^2
    const double   lambda;   //!< capillary length
    const double   mu2;            //!< mu^2...
    const double   mu;             //!< scaling
    Vector         pprev;          //!< starting point for integration
    Vector         param;          //!< current point for intergation
    double         center_v;       //!< center of lens pos, 1.0+zeta
    double         start_u;        //!< initial u position
    double         start_v;        //!< initial v position
    const double   resolution;     //!< angle resolution: search between resolution and PI-resolution
    const double   angle_control;  //!< angle control in radians (converted from L->angle_control)
    const double   shift_control;  //!< max speed = tau max
    Function       prof_alpha;
    Function       prof_theta;
    Function       zprof_zeta;
    
    //! re-compute mu2 and mu
    void change_curv(const double curvature_coeff);

    //! initialize param, center_v, start_v and start_u
    void compute_start(const double alpha,
                       const double theta,
                       const double zeta);

    //! compute a profile...
    /**
     - return GetValue(v0,final_v) or exactly zero
     - the last param are always compute
     */
    double profile(const double  alpha,
                   const double  theta,
                   const double  zeta,
                   ios::ostream *fp,
                   const bool    record = false);

    //! this is the function to get a continuous level indicator
    /**
     - if (v0>=0) return Fabs(min_of(v0,v));
     - if (v0<0) return  Fabs(max_of(v0,v));
     */
    static double GetValue(const double v0, const double v) throw();

    //! save a drawing of the lens
    static void SaveLens(const string &filename, const double shift );

    static double CriticalZetaOfTheta(const double theta)  throw();
    static double CriticalZetaOfAlpha(const double alpha)  throw();

    double find_alpha(const double theta,
                      const double zeta,
                      bool        *is_flat);

    double find_theta(const double alpha,
                      const double zeta,
                      bool        *is_flat);

    double find_zeta_max(const double theta);
    double last_space() const throw();

    static double CapVolume(const double hh);
    double last_cylinder_space() const throw();

    //! zeta must be set
    double compute_dVdu(const double u,const double v,const double phi) const throw();

    double find_theta_by_xi(const double alpha, const double zeta, bool *if_flat, double &shift);


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
    double __alpha; //!< for functions
    double __theta; //!< for functions
    double __zeta;  //!< for functions
    double __xi;    //!< for functions
    bool   __success; //!< for xi...
    
    void   ProfileEq( array<double> &dYds, double, const array<double> &Y );
    double ProfileOfAlpha(const double alpha);
    double ProfileOfTheta(const double theta);
    double  __find_top( Function &F, double p_lo, double p_up, const double res);
    double  __find_bot( Function &F, double p_lo, double p_up, const double res);
    double DeltaOfZeta(const double zeta);

public:
    size_t last_counts;
    Vector heights;
    Vector radii;
    Vector volumes;
    Vector angles;
    
    double compute_shift(const double alpha, const double theta, const double zeta);
};

#endif
