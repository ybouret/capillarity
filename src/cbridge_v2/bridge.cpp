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
    //
    // initialize parameters
    //
    //__________________________________________________________________________

    compute_start(alpha,theta,zeta);
    tao::set(pprev,param);
    const double v0 = param[BRIDGE_V];


    //__________________________________________________________________________
    //
    //
    // loop
    //
    //__________________________________________________________________________

    double       tau      = 0;
    double       ctl      = SAFETY * shift_control;
    if(fp) (*fp)("%.15g %.15g %.15g\n", param[BRIDGE_U],param[BRIDGE_V], param[BRIDGE_A]);
    while(true)
    {
        //______________________________________________________________________
        //
        //
        // initialize dtau according to position
        //
        //______________________________________________________________________
        double       dtau         = min_of(shift_control,param[BRIDGE_U]*SAFETY);


        //______________________________________________________________________
        //
        //
        // compute the angular rate and correct dtau
        //
        //______________________________________________________________________
        const double angular_rate = Fabs(mu2 *param[BRIDGE_V] - sin( param[BRIDGE_A] ) / param[BRIDGE_U]);
        if( angular_rate * dtau > angle_control )
        {
            dtau = log_round_floor(angle_control/angular_rate);
        }

        //______________________________________________________________________
        //
        //
        // take the step
        //
        //______________________________________________________________________
        const double tau_next = tau + dtau;
        ctl = min_of(ctl,dtau);
        odeint(profEq,param,tau,tau_next,ctl,NULL);
        if(!status)
        {
            // something went wrong
            tao::set(param,pprev);
            if(fp) (*fp)("%.15g %.15g %.15g\n",param[BRIDGE_U],param[BRIDGE_V],param[BRIDGE_A]);
            return GetValue(v0,param[BRIDGE_V]);
        }



        //______________________________________________________________________
        //
        //
        // check if profile is done...
        //
        //______________________________________________________________________

        //______________________________________________________________________
        //
        // check zero crossing
        //______________________________________________________________________
        if(true)
        {


            const double v_prev = pprev[BRIDGE_V];
            const double v_curr = param[BRIDGE_V];
            if(v_curr*v_prev<=0)
            {
                const double X  = clamp<double>(0,-v_prev/(v_curr-v_prev),1);
                param[BRIDGE_V] = 0;
                param[BRIDGE_U] = pprev[BRIDGE_U] + X*(param[BRIDGE_U]-pprev[BRIDGE_U]);
                param[BRIDGE_A] = pprev[BRIDGE_A] + X*(param[BRIDGE_A]-pprev[BRIDGE_A]);
                if(fp) (*fp)("%.15g %.15g %.15g\n",param[BRIDGE_U],param[BRIDGE_V],param[BRIDGE_A]);
                return 0;
            }
        }

        //______________________________________________________________________
        //
        // if u was increasing then decreases => invalid u extremum
        //______________________________________________________________________
        if(true)
        {
            const double du_prev = cos( pprev[BRIDGE_A] );
            const double du_curr = cos( param[BRIDGE_A] );
            if(du_prev>=0&&du_curr<0)
            {
                //std::cerr << "u-returning!" << std::endl;
                const double X   = clamp<double>(0.0,-du_prev/(du_curr-du_prev),1.0);
                const double v   = pprev[BRIDGE_V] + X * (param[BRIDGE_V]-pprev[BRIDGE_V]);
                const double u   = pprev[BRIDGE_U] + X * (param[BRIDGE_U]-pprev[BRIDGE_U]);
                const double phi = pprev[BRIDGE_A] + X * (param[BRIDGE_A]-pprev[BRIDGE_A]);
                param[BRIDGE_V] = v;
                param[BRIDGE_U] = u;
                param[BRIDGE_A] = phi;
                if(fp) (*fp)("%.15g %.15g %.15g\n", u,v,phi);
                return GetValue(v0,v);
            }
        }

        //______________________________________________________________________
        //
        // if dv is zero => v extremum
        //______________________________________________________________________
        if(true)
        {
            const double dv_prev = sin( pprev[BRIDGE_A] );
            const double dv_curr = sin( param[BRIDGE_A] );
            if(dv_prev*dv_curr<=0)
            {
                assert(Fabs(dv_prev)>0||Fabs(dv_curr)>0);
                const double X   = clamp<double>(0.0,-dv_prev/(dv_curr-dv_prev),1.0);
                const double v   = pprev[BRIDGE_V] + X * (param[BRIDGE_V]-pprev[BRIDGE_V]);
                const double u   = pprev[BRIDGE_U] + X * (param[BRIDGE_U]-pprev[BRIDGE_U]);
                const double phi = pprev[BRIDGE_A] + X * (param[BRIDGE_A]-pprev[BRIDGE_A]);
                param[BRIDGE_V] = v;
                param[BRIDGE_U] = u;
                param[BRIDGE_A] = phi;
                if(fp) (*fp)("%.15g %.15g %.15g\n", u,v,phi);
                return GetValue(v0,v);
            }
        }


        //______________________________________________________________________
        //
        // prepare next step
        //______________________________________________________________________
        tau = tau_next;
        tao::set(pprev,param);
        if(fp) (*fp)("%.15g %.15g %.15g\n",param[BRIDGE_U],param[BRIDGE_V],param[BRIDGE_A]);
        if(tau>1)
        {
            break;
        }

    }

    // end of simulation
    return GetValue(v0,param[BRIDGE_V]);

}

