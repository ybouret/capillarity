#ifndef BRIDGE_INCLUDED
#define BRIDGE_INCLUDED 1

#include "yocto/lua/lua-state.hpp"
#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/math/opt/minimize.hpp"
#include "yocto/math/core/tao.hpp"
#include "yocto/math/round.hpp"
#include "yocto/container/tuple.hpp"
#include "yocto/container/manager.hpp"

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

YOCTO_PAIR_DECL(STANDARD,range_t,double,vmin,double,vmax);
inline range_t() throw() : vmin(0), vmax(0) {}
inline range_t & operator=( const range_t &r ) throw()
{
    vmin = r.vmin; vmax = r.vmax; return *this;
}
inline void check() throw() { if(vmax>vmin) cswap(vmin,vmax); }
YOCTO_PAIR_END();


class Bridge
{
public:


    virtual ~Bridge() throw();
    explicit Bridge( Lua::State &L);

    const size_t   nvar;     //!< BRIDGE_N
    DEsolver       odeint;   //!< adaptive solve, initialized with L->ftol (should be 1e-4-1e-7)
    Equation       profEq;   //!< differential equation
    bool           status;   //!< no error during integration
    const double   R0;       //!< radius
    const double   A0;       //!< pi * R0^2
    const double   lambda;   //!< capillary length
    const double   mu2;            //!< mu^2 = R0^2/lambda^2
    const double   mu;             //!< scaling
    Vector         pprev;          //!< starting point for integration
    Vector         param;          //!< current point for intergation
    double         center_v;       //!< center of lens pos, 1.0+zeta
    double         start_u;        //!< initial u position
    double         start_v;        //!< initial v position
    const double   resolution;     //!< angle resolution: search between resolution and PI-resolution
    const double   angle_control;  //!< angle control in radians (converted from L->angle_control)
    const double   shift_control;  //!< max speed = tau max
    Function       prof_theta;
    Function       prof_alpha;

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

    //! compute shift from volume splitting
    double compute_shift(const double alpha, const double theta, const double zeta);

    double find_theta(const double alpha, const double zeta, bool &isFlat);
    static double CriticalZetaOfAlpha(const double alpha)  throw();
    static double CriticalZetaOfTheta(const double theta)  throw();

    range_t find_zeta_range( const double alpha, range_t &tr);

    size_t find_alpha(const double theta, const double zeta, double *alphas, bool &isFlat);


    double find_zeta_max(const double theta);



private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
    double __alpha; //!< for functions
    double __theta; //!< for functions
    double __zeta;  //!< for functions

    void   ProfileEq( array<double> &dYds, double, const array<double> &Y );
    double ProfileOfTheta( const double theta );
    double ProfileOfAlpha( const double alpha );

    //! assuming F(lo)>0, F(hi)<=0
    static double find_lower(double lo,double hi, Function &F, const double resolution);

    //! assuming F(lo)<=0, F(hi)>0
    static double find_upper(double lo,double hi, Function &F, const double resolution);

public:
    size_t last_counts;
    Vector heights;
    Vector radii;
    Vector volumes;
    Vector angles;
    container_manager_of<Vector> mgr;
};

#endif
