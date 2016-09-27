#include "bridge.hpp"

Bridge:: ~Bridge() throw()
{
}

Bridge:: Bridge( lua_State *L ) :
nvar(BRIDGE_N),
odeint( Lua::Config::Get<lua_Number>(L,"ftol") ),
profEq( this, & Bridge::ProfileEq ),
status(false),
mu(1),
mu2(1),
pprev(nvar),
param(nvar),
center_v(0),
angle_control( Lua::Config::Get<lua_Number>(L,"angle_control") ),
shift_control( Lua::Config::Get<lua_Number>(L,"shift_control") )
{
    odeint.start(nvar);
    std::cerr << "ftol="       << odeint.eps  << std::endl;
    std::cerr << "shift_control=" << shift_control << std::endl;
    std::cerr << "angle_control=" << angle_control << std::endl;
}

size_t Bridge:: curvature_coeff = 1;


void Bridge::ProfileEq(array<double> &dYds, double, const array<double> &Y)
{
    const double u   = Y[BRIDGE_U];
    const double v   = Y[BRIDGE_V];
    const double phi = Y[BRIDGE_A];
    const double C   = cos(phi);
    const double S   = sin(phi);

    if(u>0)
    {
        dYds[BRIDGE_U] = C;
        dYds[BRIDGE_V] = S;
        dYds[BRIDGE_A] = mu2 * v - S/u;
    }
    else
    {
        status  = false;
        tao::ld(dYds,0);
    }

}


void Bridge:: compute_mu(const double R0, const double capillary_length)
{
    mu2 = curvature_coeff * Square(R0) / Square(capillary_length);
    mu  = Sqrt(mu2);
}


void Bridge:: compute_start(const double alpha, const double theta, const double zeta)
{
    status          = true;
    center_v        = 1.0 + zeta;
    param[BRIDGE_U] = sin(alpha);
    param[BRIDGE_V] = zeta + (1.0-cos(alpha));
    param[BRIDGE_A] = alpha+theta-numeric<double>::pi;
}


static inline
double GetValue(const double v0, const double v)
{
    if(v0>=0)
    {
        return Fabs(min_of(v0,v));
    }
    else
    {
        return Fabs(max_of(v0,v));
    }
}


#include "yocto/math/round.hpp"

//! compute a profile
double Bridge:: profile(const double alpha,
                        const double theta,
                        const double zeta, ios::ostream *fp )
{
    const double SAFETY = 0.1;

    assert(alpha>0);
    assert(alpha<numeric<double>::pi);


    //__________________________________________________________________________
    //
    // initialize parameters
    //__________________________________________________________________________

    compute_start(alpha,theta,zeta);
    tao::set(pprev,param);
    const double v0 = param[BRIDGE_V];


    // loop
    double       tau      = 0;
    double       ctl      = SAFETY * shift_control;
    if(fp) (*fp)("%.15g %.15g %.15g\n", param[BRIDGE_U],param[BRIDGE_V], param[BRIDGE_A]);
    while(true)
    {
        //______________________________________________________________________
        //
        // initialize dtau according to position
        //______________________________________________________________________
        double       dtau         = min_of(shift_control,param[BRIDGE_U]*SAFETY);


        //______________________________________________________________________
        //
        // compute the angular rate and correct dtau
        //______________________________________________________________________
        const double angular_rate = Fabs(mu2 *param[BRIDGE_V] - sin( param[BRIDGE_A] ) / param[BRIDGE_U]);
        if( angular_rate * dtau > angle_control )
        {
            dtau = log_round_floor(angle_control/angular_rate);
        }

        //______________________________________________________________________
        //
        // take the step
        //______________________________________________________________________
        const double tau_next = tau + dtau;
        ctl = min_of(ctl,dtau);
        odeint(profEq,param,tau,tau_next,ctl,NULL);
        if(!status)
        {
            tao::set(param,pprev);
            const double u   = pprev[BRIDGE_U];
            const double v   = pprev[BRIDGE_V];
            const double phi = pprev[BRIDGE_V];


            if(fp) (*fp)("%.15g %.15g %.15g\n", u,v,phi);
            return GetValue(v0,v);
        }



        //______________________________________________________________________
        //
        // check if profile is done...
        //______________________________________________________________________


        //______________________________________________________________________
        //
        // prepare next step
        //______________________________________________________________________
        tau = tau_next;
        tao::set(pprev,param);
        std::cerr << "tau=" << tau << std::endl;
        if(fp) (*fp)("%.15g %.15g %.15g\n",param[BRIDGE_U],param[BRIDGE_V],param[BRIDGE_A]);
        if(tau>1)
        {
            break;
        }

    }

    // end of simulation
    {
        //const double u   = pprev[BRIDGE_U];
        const double v   = pprev[BRIDGE_V];
        //const double phi = pprev[BRIDGE_V];

        return GetValue(v0,v);
    }
}

