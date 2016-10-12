#include "bridge.hpp"

void Bridge:: compute_start(const double alpha, const double theta, const double Xi)
{
    status          = true;
    center_v        = 1.0 + Xi;
    param[BRIDGE_U] = sin(alpha);
    param[BRIDGE_V] = Xi + (1.0-cos(alpha));
    param[BRIDGE_A] = alpha+theta-numeric<double>::pi;
    //param[BRIDGE_Q] = 0;
}

double Bridge:: GetValue(const double v0, const double v) throw()
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
        //dYds[BRIDGE_Q] = 0;
    }
    else
    {
        status  = false;
        tao::ld(dYds,0);
    }
    
}

double Bridge:: reduced_rate( const array<double> &Y ) const throw()
{
    assert(Y.size()>=nvar);
    const double u   = Y[BRIDGE_U];
    const double v   = Y[BRIDGE_V];
    const double phi = Y[BRIDGE_A];
    const double S   = sin(phi);

    return (mu2*v-S/u)/mu2;
}

double Bridge:: reduced_curv( const array<double> &Y ) const throw()
{
    assert(Y.size()>=nvar);
    const double u   = Y[BRIDGE_U];
    const double v   = Y[BRIDGE_V];
    const double phi = Y[BRIDGE_A];
    const double S   = sin(phi);
    const double C   = cos(phi);

    const double lhs = mu2 *(S-v*C/u);
    const double rhs = (C*S)/(u*u);
    return (lhs+rhs+rhs)/mu2;
}



#include "yocto/math/round.hpp"
#include "yocto/math/point2d.hpp"
typedef point2d<double> P2D;

//! compute a profile
double Bridge:: profile(const double  alpha,
                        const double  theta,
                        const double  Xi,
                        ios::ostream *fp,
                        const bool    store_data,
                        const double  shift)
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

    compute_start(alpha,theta,Xi);
    tao::set(pprev,param);
    const double v0 = param[BRIDGE_V];
    heights.free();
    radii.free();

    if(store_data)
    {
        heights.ensure(1000);
        radii.ensure(1000);
    }
    //__________________________________________________________________________
    //
    //
    // loop
    //
    //__________________________________________________________________________
#define SAVE_STATUS() do { \
if(fp) (*fp)("%.15g %.15g %.15g\n", param[BRIDGE_U],param[BRIDGE_V]+shift, param[BRIDGE_A]);\
if(store_data) { radii.push_back(param[BRIDGE_U]); heights.push_back(param[BRIDGE_V]); } \
} while(false)

    double       tau      = 0;
    double       ctl      = SAFETY * shift_control;
    SAVE_STATUS();

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
            SAVE_STATUS();
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
                assert(Fabs(v_curr)>0 || Fabs(v_prev)>0);
                const double X  = clamp<double>(0,-v_prev/(v_curr-v_prev),1);
                param[BRIDGE_U] = pprev[BRIDGE_U] + X*(param[BRIDGE_U]-pprev[BRIDGE_U]);
                param[BRIDGE_V] = 0;
                param[BRIDGE_A] = pprev[BRIDGE_A] + X*(param[BRIDGE_A]-pprev[BRIDGE_A]);
                SAVE_STATUS();
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
                const double u   = pprev[BRIDGE_U] + X * (param[BRIDGE_U]-pprev[BRIDGE_U]);
                const double v   = pprev[BRIDGE_V] + X * (param[BRIDGE_V]-pprev[BRIDGE_V]);
                const double phi = pprev[BRIDGE_A] + X * (param[BRIDGE_A]-pprev[BRIDGE_A]);
                param[BRIDGE_V] = v;
                param[BRIDGE_U] = u;
                param[BRIDGE_A] = phi;
                SAVE_STATUS();
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
                const double u   = pprev[BRIDGE_U] + X * (param[BRIDGE_U]-pprev[BRIDGE_U]);
                const double v   = pprev[BRIDGE_V] + X * (param[BRIDGE_V]-pprev[BRIDGE_V]);
                const double phi = pprev[BRIDGE_A] + X * (param[BRIDGE_A]-pprev[BRIDGE_A]);
                param[BRIDGE_V]  = v;
                param[BRIDGE_U]  = u;
                param[BRIDGE_A]  = phi;
                SAVE_STATUS();
                return GetValue(v0,v);
            }
        }

        //______________________________________________________________________
        //
        // check not into lens
        //______________________________________________________________________
        if(true)
        {
            const P2D B(param[BRIDGE_U],param[BRIDGE_V]);
            const P2D OB(B.x,B.y-center_v);
            const double OB2 = OB.norm2();
            if(OB2<=1)
            {
                const P2D A(pprev[BRIDGE_U],pprev[BRIDGE_V]);
                const P2D BA(B,A);
                const P2D OA(A.x,A.y-center_v);
                const double a = BA.norm2();
                const double b = OB*BA;
                const double c = OB2-1.0;

                // compute reduced discriminant
                const double sD   = sqrt(max_of(0.0,b*b - a*c));
                const double beta = clamp<double>(0,(-b + sD)/a,1);
                const double X   = 1.0 - beta;
                const double u   = pprev[BRIDGE_U] + X * (param[BRIDGE_U]-pprev[BRIDGE_U]);
                const double v   = pprev[BRIDGE_V] + X * (param[BRIDGE_V]-pprev[BRIDGE_V]);
                const double phi = pprev[BRIDGE_A] + X * (param[BRIDGE_A]-pprev[BRIDGE_A]);
                param[BRIDGE_U]  = u;
                param[BRIDGE_V]  = v;
                param[BRIDGE_A]  = phi;
                SAVE_STATUS();
                return GetValue(v0,v);
            }
        }


        //______________________________________________________________________
        //
        // prepare next step
        //______________________________________________________________________
        tau = tau_next;
        tao::set(pprev,param);
        SAVE_STATUS();
    }
    
    // end of simulation
    return GetValue(v0,param[BRIDGE_V]);
    
}

