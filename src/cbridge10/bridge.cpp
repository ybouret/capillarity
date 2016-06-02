#include "bridge.hpp"
#include "yocto/math/core/tao.hpp"

Bridge:: Bridge(const double delta_degrees, const double ftol) :
nvar( BRIDGE_N ),
mu(1),
odeint(ftol),
pprev(nvar),
param(nvar),
prate(nvar),
flag(false),
eq( this, & Bridge::__Eq ),
cb( this, & Bridge::__Cb ),
delta( Deg2Rad(delta_degrees) ),
v_center(0),
mu2(0),
current_theta(0),
current_zeta(0),
fn_of_alpha( this, & Bridge:: __profile_of_alpha )
{
    odeint.start(nvar);

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

/*
double Bridge:: angular_rate(const vector<double> &Y) const throw()
{
    return mu2 * Y[BRIDGE_V] - sin( Y[BRIDGE_A] ) / Y[BRIDGE_U];
}
*/

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

double Bridge:: goodness(const double u, const double v, const double phi) const throw()
{
    //return sign_of(v);
    return Fabs(v);
    //return ((mu2*v-sin(phi)/u)/mu2);
    //return v < 0 ? -1 : (0<v ? 1 : 0);
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


double Bridge:: profile( const double alpha, const double theta, const double zeta, ios::ostream *fp )
{
    //__________________________________________________________________________
    //
    // initialize geometry
    //__________________________________________________________________________

    assert(alpha>0);
    assert(alpha<numeric<double>::pi);
    flag     = true;
    compute_start(alpha, theta, zeta);


    tao::set(pprev, param);

    //__________________________________________________________________________
    //
    // initialize step: TODO set as arguments!
    //__________________________________________________________________________
    const double max_delta_a = Deg2Rad(0.1);
    const double max_delta_l = 1.0/100.0;

    

    double tau = 0;
    if(fp) (*fp)("%.15g %.15g %.15g\n", param[BRIDGE_U],param[BRIDGE_V], param[BRIDGE_A]);
    while(true)
    {
        // initialize to max delta_l
        double       dtau = min_of(max_delta_l,param[BRIDGE_U]/10);

        const double angular_rate = Fabs(mu2 *param[BRIDGE_V] - sin( param[BRIDGE_A] ) / param[BRIDGE_U]);
        if( angular_rate * dtau > max_delta_a )
        {
            dtau = log_round_floor(max_delta_a/angular_rate);
        }

        const double tau_next = tau + dtau;
        double tctl = dtau/10.0;
        odeint(eq,param,tau,tau_next,tctl,&cb);

        if(!flag)
        {

            const double u   = pprev[BRIDGE_U];
            const double v   = pprev[BRIDGE_V];
            const double phi = pprev[BRIDGE_V];
            if(fp) (*fp)("%.15g %.15g %.15g\n", u,v,phi);
            return goodness(u, v, phi);
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
                //std::cerr << "return in lens" << std::endl;
                const double v   = pprev[BRIDGE_V] + X * (param[BRIDGE_V]-pprev[BRIDGE_V]);
                const double u   = pprev[BRIDGE_U] + X * (param[BRIDGE_U]-pprev[BRIDGE_U]);
                const double phi = pprev[BRIDGE_A] + X * (param[BRIDGE_A]-pprev[BRIDGE_A]);
                if(fp) (*fp)("%.15g %.15g %.15g\n", u,v,phi);

                return goodness(u, v, phi);
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
                const double X = clamp<double>(0.0,-du_prev/(du_curr-du_prev),1.0);
                const double v   = pprev[BRIDGE_V] + X * (param[BRIDGE_V]-pprev[BRIDGE_V]);
                const double u   = pprev[BRIDGE_U] + X * (param[BRIDGE_U]-pprev[BRIDGE_U]);
                const double phi = pprev[BRIDGE_A] + X * (param[BRIDGE_A]-pprev[BRIDGE_A]);
                if(fp) (*fp)("%.15g %.15g %.15g\n", u,v,phi);
                return goodness(u, v, phi);
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
                const double X = clamp<double>(0.0,-dv_prev/(dv_curr-dv_prev),1.0);
                const double v   = pprev[BRIDGE_V] + X * (param[BRIDGE_V]-pprev[BRIDGE_V]);
                const double u   = pprev[BRIDGE_U] + X * (param[BRIDGE_U]-pprev[BRIDGE_U]);
                const double phi = pprev[BRIDGE_A] + X * (param[BRIDGE_A]-pprev[BRIDGE_A]);
                if(fp) (*fp)("%.15g %.15g %.15g\n", u,v,phi);
                return goodness(u, v, phi);
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
                if(fp) (*fp)("%.15g %.15g %.15g\n", u,v,phi);
                return 0;
                //return goodness(u, v, phi);
            }
        }


        tau = tau_next;
        tao::set(pprev,param);
        if(fp) (*fp)("%.15g %.15g %.15g\n", param[BRIDGE_U],param[BRIDGE_V], param[BRIDGE_A]);
        if(tau>10)
            break;
    }

    return goodness(param[BRIDGE_U], param[BRIDGE_V], param[BRIDGE_A]);

}

double Bridge:: __profile_of_alpha(const double alpha)
{
    return profile(alpha, current_theta, current_zeta, NULL);
}

#include "yocto/math/opt/bracket.hpp"
#include "yocto/ios/ocstream.hpp"

#include "yocto/math/opt/bracket.hpp"
#include "yocto/math/opt/minimize.hpp"

double Bridge:: find_alpha(const double theta, const double zeta)
{
    std::cerr << "theta=" << Rad2Deg(theta) << ", zeta=" << zeta << std::endl;

    //__________________________________________________________________________
    //
    // prepare functions
    //__________________________________________________________________________
    current_theta = theta;
    current_zeta  = zeta;
    Function &f   = fn_of_alpha;
    //const double alpha_max = numeric<double>::pi-theta;


#if 1
    {
        ios::wcstream fp("find-alpha.dat");
        ios::wcstream pp("prof-alpha.dat");
        const double a_deg_max = 170;//floor(Rad2Deg(alpha_max));
        for(double a_deg = 0.2; a_deg <= a_deg_max; a_deg += 0.2 )
        {
            fp("%.15g %.15g\n", a_deg, profile(Deg2Rad(a_deg), theta, zeta, &pp) );
            pp << "\n";
        }
    }

#endif

    double alpha_lo = delta;
    double value_lo = f(alpha_lo);
    double alpha_hi = numeric<double>::pi/2;
    double value_hi = f(alpha_hi);

    if(value_hi<=0) throw exception("invalid value@90deg");

    triplet<double> X = { alpha_lo, 0, alpha_hi };
    triplet<double> F = { value_lo, 0, value_hi };

    bracket<double>::inside(f, X, F);
    optimize(f, X, F, delta);
    double alpha_opt = X.b;
    double value_opt = F.b;

    std::cerr << "alpha_opt=" << Rad2Deg(alpha_opt) << std::endl;
    std::cerr << "value_opt=" << value_opt          << std::endl;

    if(value_opt>0)
    {
        std::cerr << "no intersection!" << std::endl;
        return -1;
    }

    // right value, for any zeta
    double alpha_r = alpha_opt;
    double diff_r  = value_hi;
    {
        double alpha_top = alpha_hi;
        while((alpha_top-alpha_r)>delta)
        {
            const double alpha_mid = clamp(alpha_r,0.5*(alpha_r+alpha_top),alpha_top);
            const double value_mid = f(alpha_mid);
            if(value_mid<=0)
            {
                alpha_r = alpha_mid;
            }
            else
            {
                alpha_top = alpha_mid;
                diff_r    = value_mid;
            }
        }
    }

    std::cerr << "alpha_r=" << Rad2Deg(alpha_r) << std::endl;
    std::cerr << "diff_r =" << diff_r           << std::endl;

    double alpha = alpha_r;

    if(zeta<0&&value_lo>0)
    {
        std::cerr << "looking for alpha_l..." << std::endl;
        double alpha_l = alpha_opt;
        double diff_l  = value_lo;
        {
            double alpha_bot = alpha_lo;
            while( (alpha_l-alpha_bot)>delta )
            {
                const double alpha_mid = clamp(alpha_bot,0.5*(alpha_l+alpha_bot),alpha_l);
                const double value_mid = f(alpha_mid);
                if(value_mid<=0)
                {
                    alpha_l = alpha_mid;
                }
                else
                {
                    alpha_bot = alpha_mid;
                    diff_l    = value_mid;
                }
            }
        }
        std::cerr << "alpha_l=" << Rad2Deg(alpha_l) << std::endl;
        std::cerr << "diff_l =" << diff_l           << std::endl;
        if(diff_l<diff_r)
        {
            alpha=alpha_l;
        }
    }


    {
        ios::wcstream ap("good-alpha.dat");
        (void)profile(alpha, theta, zeta, &ap);

    }

    return alpha;
    
}

double Bridge:: compute_zeta_max(const double theta)
{
    double zeta_lo = 0.0;
    if( find_alpha(theta, zeta_lo) <= 0 )
    {
        throw exception("couldn't find bridge for zeta=%g", zeta_lo);
    }

    double zeta_hi = 1;
    while( find_alpha(theta,zeta_lo) > 0 )
    {
        zeta_hi += 0.5;
    }

    return 0;
}
