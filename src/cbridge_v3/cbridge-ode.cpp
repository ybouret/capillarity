#include "cbridge.hpp"

void Bridge:: compute_start(const double alpha,
                            const double theta,
                            const double zeta)
{
    __zeta          = zeta;
    status          = true;
    center_v        = 1.0 + zeta;
    start_v         = zeta + (1.0-cos(alpha));
    start_u         = sin(alpha);
    param[BRIDGE_U] = start_u;
    param[BRIDGE_V] = start_v;
    param[BRIDGE_A] = alpha+theta-numeric<double>::pi;
    param[BRIDGE_Q] = 0;
    param[BRIDGE_q] = 0;
    last_counts = 0;
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

double Bridge:: last_space() const throw()
{
    return  (param[BRIDGE_Q] - last_cylinder_space());
}

double Bridge:: last_cylinder_space() const throw()
{
    if(__zeta>=0)
    {
        return  CapVolume(start_v - __zeta);
    }
    else
    {
        return CapVolume(start_v - __zeta) - CapVolume(-__zeta);
    }
}


void Bridge::ProfileEq(array<double> &dYds, double, const array<double> &Y)
{
    static const double fac = -numeric<double>::pi;
    const double u   = Y[BRIDGE_U];
    const double v   = Y[BRIDGE_V];
    const double phi = Y[BRIDGE_A];
    const double C   = cos(phi);
    const double S   = sin(phi);

    if(u>0)
    {
        dYds[BRIDGE_U]   = C;
        dYds[BRIDGE_V]   = S;
        dYds[BRIDGE_A]   = mu2 * v - S/u;
        const double ff  = fac * S;
        dYds[BRIDGE_Q]   = ff * (u*u);
        if(v>=__zeta)
        {
            dYds[BRIDGE_q]   = ff * max_of<double>(0,1.0-Square(1.0+__zeta-v));
        }
        else
        {
            dYds[BRIDGE_q] = 0;
        }
        //std::cerr << "dq(" << v << ")=" << dYds[BRIDGE_q] << std::endl;
    }
    else
    {
        status  = false;
        tao::ld(dYds,0);
    }

}

double Bridge:: compute_dVdu(const double u,const double v,const double phi) const throw()
{
    static const double fac = -numeric<double>::pi;

    const double S = sin(phi);
    const double ff = S * fac;
    double       dVdtau = ff * (u*u);
    if(v>=__zeta)
    {
        dVdtau -= ff * max_of<double>(0,1.0-Square(1.0+__zeta-v));
    }
    return dVdtau / cos(phi);
}

#include "yocto/math/round.hpp"
#include "yocto/math/point2d.hpp"
typedef point2d<double> P2D;

//! compute a profile
double Bridge:: profile(const double  alpha,
                        const double  theta,
                        const double  zeta,
                        ios::ostream *fp,
                        const bool    record)
{
    const double SAFETY = 0.1;

    assert(alpha>0);
    assert(alpha<numeric<double>::pi);
    mgr.free_all();

    //__________________________________________________________________________
    //
    //
    // initialize parameters
    //
    //__________________________________________________________________________

    compute_start(alpha,theta,zeta);
    tao::set(pprev,param);
    const double v0 = start_v;

    //__________________________________________________________________________
    //
    //
    // loop
    //
    //__________________________________________________________________________
#define SAVE_STATUS() do { \
++last_counts;\
if(fp) (*fp)("%.15g %.15g %.15g\n", param[BRIDGE_U],param[BRIDGE_V], param[BRIDGE_A]);\
if(record) {\
heights.push_back(param[BRIDGE_V]);                \
radii.push_back(param[BRIDGE_U]);                  \
volumes.push_back(param[BRIDGE_Q]-param[BRIDGE_q]);\
angles.push_back(param[BRIDGE_A]);                 \
}\
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
                tao::setbar(param,pprev,X,param);
                param[BRIDGE_V] = 0;
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
                tao::setbar(param,pprev,X,param);
                SAVE_STATUS();
                return GetValue(v0,param[BRIDGE_V]);
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
                tao::setbar(param,pprev,X,param);
                SAVE_STATUS();
                return GetValue(v0,param[BRIDGE_V]);
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
                const double X    = 1.0 - beta;
                tao::setbar(param,pprev,X,param);
                SAVE_STATUS();
                return GetValue(v0,param[BRIDGE_V]);
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

#include "yocto/sort/quick.hpp"
#include "yocto/math/dat/linear.hpp"

double Bridge:: compute_shift(const double alpha, const double theta, const double zeta)
{

    if( profile(alpha, theta, zeta, NULL, false) > 0)
    {
        throw exception("compute_shift: NO BRIDGE!");
    }
    mgr.ensure_all(last_counts);

    (void)profile(alpha, theta, zeta, NULL, true);
    if(heights.size()<=2)
    {
        return 0; // must be flat...
    }
    co_qsort(volumes,heights);
    if(false)
    {
        ios::wcstream fp("shift.dat");
        for(size_t i=1;i<=volumes.size();++i)
        {
            fp("%g %g\n", volumes[i], heights[i]);
        }
    }
    const double shift = linear((volumes.back()+volumes.front())*0.5, volumes, heights);

    //std::cerr << "volume=" << volumes.back() << std::endl;
    //std::cerr << "shift =" << shift << "/zeta=" << zeta << std::endl;
    return shift;
}

