#include "bridge.hpp"
#include "yocto/math/core/tao.hpp"

Bridge:: Bridge(const double ftol) :
nvar( BRIDGE_N ),
mu(1),
odeint(ftol),
pprev(nvar),
param(nvar),
flag(false),
eq( this, & Bridge::__Eq ),
cb( this, & Bridge::__Cb ),
v_center(0)
{
    odeint.start(nvar);
    
}

Bridge:: ~Bridge() throw()
{
}


void Bridge:: __Eq( array<double> &dYdt, double, const array<double> &Y)
{
    const double u   = Y[BRIDGE_U];
    const double v   = Y[BRIDGE_V];
    const double phi = Y[BRIDGE_A];
    const double C   = cos(phi);
    const double S   = sin(phi);

    if(u>0)
    {
        dYdt[BRIDGE_U] = C;
        dYdt[BRIDGE_V] = S;
        dYdt[BRIDGE_A] = v/(mu*mu) - S/u;
    }
    else
    {
        flag = false;
        tao::ld(dYdt,0);
    }
}


void Bridge:: __Cb( array<double> &Y, double )
{
    const double u = Y[BRIDGE_U];

    // negative radius detected
    if(u<=0)
    {
        std::cerr << "negative radius" << std::endl;
        flag = false;
        return;
    }

    const double v = Y[BRIDGE_V];
    const double dv = v_center - v;
    if( Hypotenuse(u, dv) < 1.0 )
    {
        std::cerr << "return in lens" << std::endl;
        flag = false;
        return;
    }
}


double Bridge:: profile( const double alpha, const double theta, const double zeta, ios::ostream *fp )
{
    //__________________________________________________________________________
    //
    // initialize geometry
    //__________________________________________________________________________

    assert(alpha>0);
    assert(alpha<numeric<double>::pi);
    flag = true;
    v_center = 1.0+zeta;
    param[BRIDGE_U] = sin(alpha);
    param[BRIDGE_V] = v_center - cos(alpha);
    param[BRIDGE_A] = alpha+theta-numeric<double>::pi;


    tao::set(pprev, param);

    //__________________________________________________________________________
    //
    // initialize step: TODO set as arguments!
    //__________________________________________________________________________
    const double dtau_lower = param[BRIDGE_U]/100.0;
    const double dtau_scale = 1.0/100;
    double dtau = min_of(dtau_lower,dtau_scale);
    double tctl = dtau/10.0;

    double tau = 0;
    if(fp) (*fp)("%g %g %g\n", param[BRIDGE_U],param[BRIDGE_V],param[BRIDGE_A]);
    while(true)
    {
        const double tau_next = tau + dtau;
        odeint(eq,param,tau,tau_next,tctl,&cb);

        if(!flag)
            break;

        //______________________________________________________________________
        //
        // would break...
        //______________________________________________________________________
        const double dr_prev = cos( pprev[BRIDGE_A] );
        const double dr_curr = cos( param[BRIDGE_A] );
        if(dr_prev>=0&&dr_curr<0)
        {
            std::cerr << "returning!" << std::endl;
            break;
        }

        tau = tau_next;
        tao::set(pprev,param);
        dtau = tctl;
        if(fp) (*fp)("%g %g %g\n", param[BRIDGE_U],param[BRIDGE_V],param[BRIDGE_A]);
        if(tau>10)
            break;
    }

    return 0;
}
