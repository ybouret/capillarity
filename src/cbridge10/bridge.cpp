#include "bridge.hpp"
#include "yocto/math/core/tao.hpp"

Bridge:: Bridge(const double delta_degrees,
                const double ftol,
                const double actrl_degrees,
                const double sctrl) :
nvar( BRIDGE_N ),
mu(1),
odeint(ftol),
pprev(nvar),
param(nvar),
flag(false),
eq( this, & Bridge::__Eq ),
cb( this, & Bridge::__Cb ),
delta( Deg2Rad(delta_degrees) ),
angle_control( Deg2Rad(actrl_degrees) ),
shift_control( sctrl ),
surface0(this, & Bridge:: __surfac0_of_theta ),
surface_of_zeta(this, & Bridge:: __surface_of_zeta),
v_center(0),
mu2(0),
current_theta(0),
current_zeta(0),
current_alpha(0),
fn_of_alpha( this, & Bridge:: __profile_of_alpha ),
fn_of_theta( this, & Bridge:: __profile_of_theta ),
check0(this, & Bridge:: __check0)
{
    odeint.start(nvar);
    std::cerr << "odeint.eps=" << odeint.eps << std::endl;
    std::cerr << "angular_search=" << Rad2Deg(delta) << std::endl;
    std::cerr << "angle_control="  << Rad2Deg(angle_control) << std::endl;
    std::cerr << "shift_control="  << shift_control          << std::endl;
}


Bridge:: ~Bridge() throw()
{
}



void   Bridge:: compute_rates(array<double> &dYdt, const array<double> &Y)  throw()
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
        dYdt[BRIDGE_A] = mu2 * v - S/u;
    }
    else
    {
        flag  = false;
        tao::ld(dYdt,0);
    }

}

void Bridge:: __Eq( array<double> &dYdt, double, const array<double> &Y)
{
    compute_rates(dYdt,Y);
}

void Bridge:: __Cb( array<double> &Y, double )
{
    const double u = Y[BRIDGE_U];
    //const double v  = Y[BRIDGE_V];

    // negative radius detected
    if(u<=0)
    {
        std::cerr << "negative radius" << std::endl;
        flag  = false;
        return;
    }

}



void Bridge:: compute_start(const double alpha, const double theta, const double zeta)
{
    v_center = 1.0+zeta;
    mu2      = mu*mu;
    param[BRIDGE_U] = sin(alpha);
    param[BRIDGE_V] = v_center - cos(alpha);
    param[BRIDGE_A] = alpha+theta-numeric<double>::pi;
}

#include "yocto/math/point2d.hpp"
#include "yocto/math/round.hpp"


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

double Bridge:: profile(const double  alpha,
                        const double  theta,
                        const double  zeta,
                        ios::ostream *fp )
{
    //__________________________________________________________________________
    //
    // initialize geometry
    //__________________________________________________________________________
    static const double SAFETY = 10.0;
    assert(alpha>0);
    assert(alpha<numeric<double>::pi);
    flag     = true;
    compute_start(alpha, theta, zeta);
    //std::cerr << "Bridge.mu=" << mu << "|profile(alpha=" << Rad2Deg(alpha) << ",theta=" << Rad2Deg(theta) << ",zeta=" << zeta << ")" << std::endl;

    tao::set(pprev, param);

    //__________________________________________________________________________
    //
    // initialize step
    //__________________________________________________________________________
    const double v0          = param[BRIDGE_V];


    double tau = 0;
    if(fp) (*fp)("%.15g %.15g %.15g\n", param[BRIDGE_U],param[BRIDGE_V], param[BRIDGE_A]);
    while(true)
    {
        //______________________________________________________________________
        //
        // initialize dtau according to position
        //______________________________________________________________________
        double       dtau         = min_of(shift_control,param[BRIDGE_U]/SAFETY);


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
        double       tctl     = dtau/SAFETY;
        odeint(eq,param,tau,tau_next,tctl,&cb);

        if(!flag)
        {

            const double u   = pprev[BRIDGE_U];
            const double v   = pprev[BRIDGE_V];
            const double phi = pprev[BRIDGE_V];
            param[BRIDGE_V] = v;
            param[BRIDGE_U] = u;
            param[BRIDGE_A] = phi;

            if(fp) (*fp)("%.15g %.15g %.15g\n", u,v,phi);
            return GetValue(v0,v);
        }


        //______________________________________________________________________
        //
        // would break...
        //______________________________________________________________________
        if(true)
        {
            // if returning into lens
            const point2d<double> J(param[BRIDGE_U],param[BRIDGE_V]-v_center);
            const double Jnorm = Hypotenuse(J.x,J.y);

            if( Jnorm < 1.0 )
            {
                // assuming pprev is outside
                const point2d<double> I(pprev[BRIDGE_U],pprev[BRIDGE_V]-v_center);
                const point2d<double> IJ(J.x-I.x,param[BRIDGE_V]-pprev[BRIDGE_V]);
                const double a   = IJ.norm2();
                const double b   = I*IJ;
                const double c   = max_of(Hypotenuse(I.x, I.y),1.0)-1;
                const double DeltaPrime = max_of(b*b-a*c,0.0);
                const double X = clamp<double>(0, -(b+sqrt(DeltaPrime))/a, 1);

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

        // if u was increasing then decreases
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

        if(true)
        {
            // if dv is zero
            const double dv_prev = sin( pprev[BRIDGE_A] );
            const double dv_curr = sin( param[BRIDGE_A] );
            if(dv_prev*dv_curr<=0)
            {
                //std::cerr << "z-extremum" << std::endl;
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

        if(true)
        {
            // if v crosses 0
            const double v_prev = pprev[BRIDGE_V];
            const double v_curr = param[BRIDGE_V];
            if(v_prev*v_curr<=0)
            {
                //std::cerr << "z-crossing" << std::endl;
                assert(Fabs(v_prev)>0||Fabs(v_curr)>0);
                const double X = clamp<double>(0,-v_prev/(v_curr-v_prev),1.0);
                const double v   = 0;//pprev[BRIDGE_V] + X * (param[BRIDGE_V]-pprev[BRIDGE_V]);
                const double u   = pprev[BRIDGE_U] + X * (param[BRIDGE_U]-pprev[BRIDGE_U]);
                const double phi = pprev[BRIDGE_A] + X * (param[BRIDGE_A]-pprev[BRIDGE_A]);
                param[BRIDGE_V] = v;
                param[BRIDGE_U] = u;
                param[BRIDGE_A] = phi;

                if(fp) (*fp)("%.15g %.15g %.15g\n", u,v,phi);
                return GetValue(v0,v);
            }
        }


        tau = tau_next;
        tao::set(pprev,param);
        if(fp) (*fp)("%.15g %.15g %.15g\n", param[BRIDGE_U],param[BRIDGE_V], param[BRIDGE_A]);
        if(tau>10)
            break;
    }

    return GetValue(v0,param[BRIDGE_V]);
}


double Bridge:: __profile_of_alpha(const double alpha)
{
    return profile(alpha, current_theta, current_zeta, NULL);
}

double Bridge:: __profile_of_theta(const double theta)
{
    return profile(current_alpha, theta, current_zeta, NULL);
}


#include "yocto/math/opt/bracket.hpp"
#include "yocto/ios/ocstream.hpp"

#include "yocto/math/opt/bracket.hpp"

bool Bridge:: __check0(const triplet<double> &X, const triplet<double> &F)
{
    return F.b<=0;
}

double Bridge:: find_alpha(const double theta, const double zeta)
{
    //std::cerr << "theta=" << Rad2Deg(theta) << ", zeta=" << zeta << std::endl;

    //__________________________________________________________________________
    //
    // prepare functions
    //__________________________________________________________________________
    current_theta = theta;
    current_zeta  = zeta;
    Function &f   = fn_of_alpha;

#if 0
    {
        ios::wcstream fp("find-alpha.dat");
        ios::wcstream pp("prof-alpha.dat");
        for(double a_deg = 0.1; a_deg <= 179; a_deg += 0.2 )
        {
            const double alpha = Deg2Rad(a_deg);
            const double ans   = profile(alpha, theta, zeta, &pp);
            fp("%g %g %g\n", a_deg, ans, sin(param[BRIDGE_A]));
            pp << "\n";
        }
    }
#endif

    //-- take full interval
    const double alpha_lower = delta;
    const double alpha_upper = numeric<double>::pi - delta;
    const double value_lower = f(alpha_lower);
    const double value_upper = f(alpha_upper);

    if(value_upper<=0)
        throw exception("upper angle has negative value!");

    //-- best effort minimum
    triplet<double> X = { alpha_lower, 0, alpha_upper };
    triplet<double> F = { value_lower, 0, value_upper };
    bracket<double>::inside(f,X,F);

    if(F.b>0)
    {
        //std::cerr << "Need to optimize..." << std::endl;
        static const double xtol = log_round_floor(numeric<double>::sqrt_epsilon*numeric<double>::pi);
        optimize1D<double>::run_until(check0,f, X, F, xtol);
        //TODO: check on result !
    }

    const double alpha_optim = X.b;
    const double value_optim = F.b;
    //std::cerr << "alpha_optim=" << Rad2Deg(alpha_optim) << std::endl;
    //std::cerr << "value_optim=" << value_optim << std::endl;

    if(value_optim>0)
    {
        // no possible intersection
        return -1;
    }

    (void)f(alpha_optim);
    double slope_optim = sin( param[BRIDGE_A] ); // slope at intercept

    // find alpha_top: smallest not zero value
    double alpha_top = alpha_upper;
    double slope_top = slope_optim;
    {
        double alpha_tmp = alpha_optim;
        while(alpha_top-alpha_tmp>delta)
        {
            const double alpha_mid = clamp(alpha_tmp,0.5*(alpha_tmp+alpha_top),alpha_top);
            const double value_mid = f(alpha_mid);
            if(value_mid>0)
            {
                alpha_top = alpha_mid;
            }
            else
            {
                alpha_tmp = alpha_mid;
                slope_top = sin( param[BRIDGE_A] );
            }
        }
    }

    double alpha = alpha_top;


    //std::cerr << "alpha=" << Rad2Deg(alpha) << std::endl;
    if(zeta<0&&value_lower>0)
    {
        //std::cerr << "There exists a secondary value!" << std::endl;
        double alpha_bot = alpha_lower;
        double alpha_tmp = alpha_optim;
        double slope_bot = slope_optim;
        while(alpha_tmp-alpha_bot>delta)
        {
            const double alpha_mid = clamp(alpha_bot,0.5*(alpha_bot+alpha_tmp),alpha_tmp);
            const double value_mid = f(alpha_mid);
            if(value_mid>0)
            {
                alpha_bot = alpha_mid;
            }
            else
            {
                alpha_tmp = alpha_mid;
                slope_bot = sin( param[BRIDGE_A] );
            }
        }

        //std::cerr << "slope_bot=" << slope_bot << std::endl;
        //std::cerr << "slope_top=" << slope_top << std::endl;
        if(Fabs(slope_bot)<Fabs(slope_top))
        {
            alpha = alpha_bot;
        }

    }

#if 0
    {
        ios::wcstream fp("good-alpha.dat");
        (void) profile(alpha, theta, zeta, &fp);
    }
#endif

    return alpha;
}


double Bridge:: find_surface( const double theta, const double zeta)
{
    return numeric<double>::pi * Square( sin( find_alpha(theta, zeta) ) );
}


double Bridge:: __surfac0_of_theta(const double theta)
{
    return find_surface(theta, 0.0);
}

double Bridge:: __surface_of_zeta(const double zeta)
{
    return find_surface(current_theta,zeta);
}

double Bridge:: compute_zeta_max(const double theta)
{
    //std::cerr << "find zeta_max(theta=" << Rad2Deg(theta) << ")" << std::endl;
    double zeta_lo = 0.0;
    if( find_alpha(theta, zeta_lo) <= 0 )
    {
        throw exception("couldn't find bridge for zeta=%g", zeta_lo);
    }
    double zeta_hi = 1;
    while( find_alpha(theta,zeta_hi) > 0 )
    {
        zeta_hi += 0.5;
        //std::cerr << "zeta_hi=" << zeta_hi << std::endl;
    }
    
    while( zeta_hi - zeta_lo > odeint.eps )
    {
        const double zeta_mid  = clamp(zeta_lo,0.5*(zeta_lo+zeta_hi),zeta_hi);
        const double alpha_mid = find_alpha(theta,zeta_mid);
        //std::cerr << "alpha(" << zeta_mid << ")=" << alpha_mid << std::endl;
        if(alpha_mid<0)
        {
            zeta_hi = zeta_mid;
        }
        else
        {
            zeta_lo = zeta_mid;
        }
    }
    
    return zeta_lo;
}

double Bridge:: find_theta( const double alpha, const double zeta )
{
    current_alpha = alpha;
    current_zeta  = zeta;
    Function &f   = fn_of_theta;

    const double theta_hi = numeric<double>::pi - delta;
    const double value_hi = f(theta_hi);

#if 1
    {
        ios::wcstream pp("prof-theta.dat");
        ios::wcstream fp("find-theta.dat");
        for(double theta_deg=0.5;theta_deg<=175;theta_deg+=1)
        {
            const double ans = profile(alpha, Deg2Rad(theta_deg), zeta, &pp);
            pp << "\n";
            fp("%g %g\n", theta_deg,ans);
        }
    }
#endif


    if(value_hi<=0)
    {
        // no possible intercept...
        std::cerr << "no possible intercept level-1" << std::endl;
        return -1;
    }

    const double theta_lo = delta;
    const double value_lo = f(theta_lo);

    triplet<double> X = { theta_lo, 0, theta_hi };
    triplet<double> F = { value_lo, 0, value_hi };

    bracket<double>::inside(f,X,F);
    if(F.b>0 && !optimize1D<double>::run_until(check0, f, X, F, 0) )
    {
        // no possible intercept
        std::cerr << "no possible intercept level-2" << std::endl;
        return -1;
    }

    double theta_l = X.b;
    double theta_r = theta_hi;
    while( (theta_r-theta_l) > delta )
    {
        const double theta_m = clamp(theta_l,0.5*(theta_l+theta_r),theta_r);
        const double value_m = f(theta_m);
        if(value_m<=0)
        {
            theta_l = theta_m;
        }
        else
        {
            theta_r = theta_m;
        }
    }

    const double theta = theta_l;
    std::cerr << "theta=" << Rad2Deg(theta) << std::endl;
#if 1
    {
        ios::wcstream fp("good-theta.dat");
        (void) profile(alpha, theta, zeta, &fp);
    }
#endif

    return theta;
}

