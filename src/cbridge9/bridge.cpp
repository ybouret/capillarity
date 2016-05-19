#include "bridge.hpp"
#include "yocto/math/core/tao.hpp"

Bridge:: Bridge() :
nvar(3),
flag(true),
capillary_length(1),
current_height(0),
current_center(0),
current_lens(NULL),
param(nvar),
pprev(nvar),
odeint(1e-7),
Eq( this, &Bridge::__Eq ),
Cb( this, &Bridge::__Cb )
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
    flag           = true;
    current_height = height;
    current_center = height+lens.R0;
    current_lens   = &lens;

    lens.starting_point(param, alpha, theta, height);

    const double ds      = capillary_length/1000.0;
    double       ds_ctrl = ds/10.0;
    double       s       = 0;
    size_t       iter    = 1;

    if(fp) (*fp)("%g %g %g %g\n",param[1],param[2],s,Rad2Deg(param[3]));


    //__________________________________________________________________________
    //
    // initialize loop
    //__________________________________________________________________________
    tao::set(pprev,param);

    const double z0 = param[BRIDGE_Z];

    while(true)
    {
        double s_next = iter * ds;
        odeint(Eq,param,s,s_next,ds_ctrl,&Cb);
        s=s_next;
        ++iter;
        if(fp) (*fp)("%g %g %g %g\n",param[1],param[2],s,Rad2Deg(param[3]));
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
            //return 0;
        }


        // going backwards
        {
            const double drds_prev = cos( pprev[BRIDGE_A] );
            const double drds_curr = cos( param[BRIDGE_A] );
            if(drds_prev>0 && drds_curr<=0)
            {
                std::cerr << "going back->stop" << std::endl;
                return z_curr;
            }
        }

        // dzds extremum
        {
            const double dzds_prev = sin( pprev[BRIDGE_A] );
            const double dzds_curr = sin( param[BRIDGE_A] );
            if(dzds_curr*dzds_prev<=0)
            {
                std::cerr << "z extremum->stop" << std::endl;
                break;
            }
        }


        tao::set(pprev, param);

#if 0
        const double z1 = param[BRIDGE_Z];
        if(z0*z1<=0)
        {
            std::cerr << "\t\tcrossing" << std::endl;
            z0=z1;
            break;
        }
        z0 = z1;
#endif
    }
    
    const double zz = param[BRIDGE_Z];
    //return (zz<=0?-1:1);
    return zz;
}