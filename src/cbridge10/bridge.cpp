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
    return v;
    //return ((mu2*v-sin(phi)/u)/mu2);
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
    v_center = 1.0+zeta;
    mu2      = mu*mu;
    param[BRIDGE_U] = sin(alpha);
    param[BRIDGE_V] = v_center - cos(alpha);
    param[BRIDGE_A] = alpha+theta-numeric<double>::pi;


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
                return goodness(u, v, phi);
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

double Bridge:: find_alpha(const double theta, const double zeta)
{
    // prepare functions
    current_theta = theta;
    current_zeta  = zeta;
    //Function &F  = fn_of_alpha;

#if 1
    {
        ios::wcstream fp("find-alpha.dat");
        ios::wcstream pp("prof-alpha.dat");
        for(double a_deg = 0.5; a_deg <= 90; a_deg += 0.5 )
        {
            fp("%.15g %.15g\n", a_deg, profile(Deg2Rad(a_deg), theta, zeta, &pp) );
            pp << "\n";
        }
    }
#endif

    
    return 0;
}


