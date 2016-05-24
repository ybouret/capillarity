#include "bridge.hpp"
#include "yocto/math/core/tao.hpp"

Bridge:: Bridge() :
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
FnOfTheta(this, &Bridge::__profile_of_theta)
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
    std::cerr << "alpha =" << Rad2Deg(alpha) << std::endl;
    std::cerr << "theta =" << Rad2Deg(theta) << std::endl;
    std::cerr << "height=" << height << std::endl;

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
        }
    }
    
    
    const double zz = param[BRIDGE_Z];
    //return fabs(zz);
    return zz;
}


double Bridge:: __profile_of_alpha(const double alpha)
{
    current_alpha = alpha;
    return __compute_profile();
}

double Bridge:: __profile_of_theta( const double theta )
{
    std::cerr << "profile_of_theta=" << theta << std::endl;
    current_theta = theta;
    return __compute_profile();
}


double Bridge:: FindTheta( Lens &lens, const double alpha, const double height )
{
    std::cerr << "alpha=" << Rad2Deg(alpha) << "; h=" << height << std::endl;
    current_lens   = &lens;
    current_alpha  = alpha;
    current_height = height;
    current_fp     = NULL;
    Function &F    = FnOfTheta;



    return 0;

#if 0
    double theta_min = Deg2Rad(0.1);
    if(F(theta_min)>0)
    {
        return Deg2Rad(double(-1));
    }

    return Deg2Rad(double(1));
#endif

}



