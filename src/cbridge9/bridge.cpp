#include "bridge.hpp"
#include "yocto/math/core/tao.hpp"

Bridge:: Bridge(const double delta_deg) :
nvar(3),
flag(true),
capillary_length(1),
current_height(0),
current_center(0),
current_alpha(0),
current_theta(0),
current_lens(NULL),
current_fp(0),
param(nvar),
pprev(nvar),
odeint(1e-7),
Eq( this, &Bridge::__Eq ),
Cb( this, &Bridge::__Cb ),
FnOfAlpha(this, &Bridge::__profile_of_alpha),
FnOfTheta(this, &Bridge::__profile_of_theta),
delta( Deg2Rad(clamp<double>(1e-3,delta_deg,1)) )
{
    odeint.start(nvar);
}

Bridge:: ~Bridge() throw()
{
}



void Bridge:: __Eq(Array &dQds, const double s, const Array &Q)
{
    const double r   = Q[BRIDGE_R];
    if(flag&&r>0)
    {
        const double z   = Q[BRIDGE_Z];
        const double phi = Q[BRIDGE_A];


        const double C = cos(phi);
        const double S = sin(phi);

        const double L2 = capillary_length*capillary_length*2;
        dQds[BRIDGE_R]  = C; // dr/ds
        dQds[BRIDGE_Z]  = S; // dz/ds
        dQds[BRIDGE_A]  = (z/L2-S/r); // dphi/ds

    }
    else
    {
        tao::ld(dQds,0);
    }
}


void Bridge:: __Cb(Array       &Q,
                   const double s)
{
    const double r = Q[BRIDGE_R];
    if(r<=0)
    {
        // invalid radius
        std::cerr << "r<=0" << std::endl;
        flag=false;
        return;
    }

    {
        const double z      = Q[BRIDGE_Z];
        const double dx     = r;
        const double dy     = current_center-z;
        const double rr     = Hypotenuse(dx, dy);
        const double alpha  = 2.0*atan(dx/(dy+rr));
        const double lens_r = current_lens->R(alpha);
        if(rr<lens_r)
        {
            // came back into bridge
            std::cerr << "rr=" << rr << "<" << lens_r << std::endl;
            flag = false;
            return;
        }
    }
}

#include "yocto/ios/ocstream.hpp"

double Bridge:: compute_profile(Lens          &lens,
                                const double   alpha,
                                const double   theta,
                                const double   height,
                                ios::ostream *fp)
{
    //std::cerr << "alpha =" << Rad2Deg(alpha) << ", theta =" << Rad2Deg(theta) << ",  height=" << height << std::endl;

    //__________________________________________________________________________
    //
    // initializing startup point
    //__________________________________________________________________________
    current_height = height;
    current_alpha  = alpha;
    current_theta  = theta;
    current_lens   = &lens;
    current_fp     = fp;

    //__________________________________________________________________________
    //
    // run
    //__________________________________________________________________________
    return __compute_profile();

}

double Bridge:: __compute_profile()
{
    flag = true;
    current_center = current_height+current_lens->R0;
    current_lens->starting_point(param, current_alpha, current_theta, current_height);

    const double s_max   = Fabs(param[BRIDGE_R]);
    const double s_cap   = capillary_length;
    double       ds      = min_of(s_max,s_cap)/100.0;
    const double ds_max  = s_cap/10.0;
    double       ds_ctrl = ds/10.0;
    double       s       = 0;
    size_t       iter    = 1;

    if(current_fp) (*current_fp)("%g %g %g %g\n",param[1],param[2],s,Rad2Deg(param[3]));


    //__________________________________________________________________________
    //
    // initialize loop
    //__________________________________________________________________________
    tao::set(pprev,param);

    const double z0  = param[BRIDGE_Z];
    const double dz0 = sin(param[BRIDGE_A]);
    if(z0*dz0<=0)
    {
        while(true)
        {
            double s_next = s+ds;
            odeint(Eq,param,s,s_next,ds_ctrl,&Cb);
            ds = min_of(ds_ctrl,ds_max);
            ++iter;
            if(current_fp) (*current_fp)("%g %g %g %g\n",param[1],param[2],s,Rad2Deg(param[3]));
            if(!flag) break;


            //______________________________________________________________________
            //
            // finding stop conditions pprev->param
            //______________________________________________________________________

            // crossing the line
            const double z_curr = param[BRIDGE_Z];

            if(z_curr*z0<=0)
            {
                //std::cerr << "crossed->stop" << std::endl;
                return 0;
            }


            // going backwards
            {
                const double drds_prev = cos( pprev[BRIDGE_A] );
                const double drds_curr = cos( param[BRIDGE_A] );
                if(drds_prev>0 && drds_curr<=0)
                {
                    //std::cerr << "going back->stop" << std::endl;
                    break;
                }
            }

            // dzds extremum
            {
                const double dzds_prev = sin( pprev[BRIDGE_A] );
                const double dzds_curr = sin( param[BRIDGE_A] );
                if(dzds_curr*dzds_prev<=0)
                {
                    //std::cerr << "z extremum->stop" << std::endl;
                    break;
                }
            }


            tao::set(pprev, param);
            s=s_next;
            if(s>=10000*capillary_length)
                break;
        }
    }


    const double zz = param[BRIDGE_Z];
    return fabs(zz);
    //return zz;
}


double Bridge:: __profile_of_alpha(const double alpha)
{
    current_alpha = alpha;
    return __compute_profile();
}

double Bridge:: __profile_of_theta( const double theta )
{
    current_theta = theta;
    return __compute_profile();
}

#include "yocto/math/opt/bracket.hpp"
#include "yocto/math/opt/minimize.hpp"

double Bridge:: FindTheta( Lens &lens, const double alpha, const double height )
{
    std::cerr << "alpha=" << Rad2Deg(alpha) << "; h=" << height << std::endl;
    current_lens   = &lens;
    current_alpha  = alpha;
    current_height = height;
    current_fp     = NULL;
    Function &F    = FnOfTheta;

#if 0
    {
        ios::wcstream fp("find-theta.dat");
        for(double theta_deg=1;theta_deg<=179;theta_deg+=0.5)
        {
            fp("%g %g\n", theta_deg, F( Deg2Rad(theta_deg) ) );
        }
    }
#endif

    const double th_lo = delta;
    const double fn_lo = F(th_lo);
    const double th_hi = numeric<double>::pi-delta;
    const double fn_hi  = F(th_hi);

    triplet<double> th = { th_lo, 0, th_hi };
    triplet<double> fn = { fn_lo, 0, fn_hi };
    bracket<double>::inside(F, th, fn);
    minimize(F, th, fn, 1e-5);
    //std::cerr << "th=" << th << std::endl;
    //std::cerr << "fn=" << fn << std::endl;

    const double th_md = th.b;
    const double fn_md = fn.b;
    //std::cerr << "th_md=" << Rad2Deg(th_md) << std::endl;
    //std::cerr << "fn_md=" << fn_md << std::endl;

    if(fn_md>0)
    {
        return -1; // no intersection
    }

    if(fn_hi>0)
    {

        double th_left  = th_md;
        double th_right = th_hi;
        while(th_right-th_left>delta)
        {
            const double th_middle = clamp(th_left, 0.5*(th_left+th_right), th_right);
            const double fn_middle = F(th_middle);
            if(fn_middle<=0)
            {
                th_left = th_middle;
            }
            else
            {
                th_right = th_middle;
            }
        }
        return 0.5*(th_left+th_right);

    }
    else
    {
        double th_left  = th_lo;
        double th_right = th_md;
        while(th_right-th_left>delta)
        {
            const double th_middle = clamp(th_left, 0.5*(th_left+th_right), th_right);
            const double fn_middle = F(th_middle);
            //std::cerr << Rad2Deg(th_middle) << " => " << fn_middle << std::endl;
            if(fn_middle<=0)
            {
                th_right = th_middle;
            }
            else
            {
                th_left = th_middle;
            }
        }
        return 0.5*(th_left+th_right);
    }

}

double Bridge:: FindAlpha( Lens &lens, const double theta, const double height )
{
    std::cerr << "theta=" << Rad2Deg(theta) << "; h=" << height << std::endl;
    current_lens   = &lens;
    current_theta  = theta;
    current_height = height;
    current_fp     = NULL;
    Function &F    = FnOfAlpha;

#if 0
    {
        ios::wcstream fp("find-alpha.dat");
        for(double alpha_deg=0.1;alpha_deg<=179;alpha_deg+=0.1)
        {
            fp("%g %g\n", alpha_deg, F( Deg2Rad(alpha_deg) ) );
        }
    }
#endif

    const double al_lo = delta;
    const double al_hi = numeric<double>::pi - delta;
    const double fn_lo = F(al_lo);
    const double fn_hi = F(al_hi);

    triplet<double> al = { al_lo, 0, al_hi };
    triplet<double> fn = { fn_lo, 0, fn_hi };
    bracket<double>::inside(F, al, fn);
    minimize(F, al, fn, 1e-5);
    //std::cerr << "al=" << al << std::endl;
    //std::cerr << "fn=" << fn << std::endl;
    const double al_mn = al.b;
    const double fn_mn = fn.b;

    if(fn_mn>0)
    {
        return -1; // not possible
    }


    double al_left  = al_mn;
    double al_right = al_hi;
    while(al_right-al_left>delta)
    {
        const double al_middle = clamp(al_left,0.5*(al_left+al_right),al_right);
        const double fn_middle = F(al_middle);
        if(fn_middle<=0)
        {
            al_left = al_middle;
        }
        else
        {
            al_right = al_middle;
        }
    }

    const double ans = 0.5*(al_left+al_right);
    lens.starting_point(pprev, ans, theta, height);
    const double z0 = pprev[BRIDGE_Z];
    if(z0>=0)
    {
        return ans;
    }
    else
    {
        al_left  = al_lo;
        al_right = al_mn;
        while(al_right-al_left>delta)
        {
            const double al_middle = clamp(al_left,0.5*(al_left+al_right),al_right);
            const double fn_middle = F(al_middle);
            if(fn_middle<=0)
            {
                al_right = al_middle;
            }
            else
            {
                al_left = al_middle;
            }
        }
        return 0.5*(al_left+al_right);
    }
}



