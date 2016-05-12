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

        if(C<=0)
        {
            std::cerr << "going back!" << std::endl;
            flag = false;
        }
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

double Bridge:: compute_profile(Lens        &lens,
                                const double   alpha,
                                const double   theta,
                                const double   height,
                                ios::ostream *fp)
{
    std::cerr << "alpha =" << Rad2Deg(alpha) << std::endl;
    std::cerr << "theta =" << Rad2Deg(theta) << std::endl;
    std::cerr << "height=" << height << std::endl;

    flag           = true;
    current_height = height;
    current_center = height+lens.R0;
    current_lens   = &lens;

    lens.starting_point(param, alpha, theta, height);

    const double ds      = capillary_length/100.0;
    double       ds_ctrl = ds/10.0;
    double       s       = 0;
    size_t       iter    = 1;

    if(fp) (*fp)("%g %g %g %g\n",param[1],param[2],s,Rad2Deg(param[3]));

    double dz0 = sin( param[BRIDGE_A] );

    //double dzds_curr = sin( param[BRIDGE_A] );
    while(true)
    {
        double s_next = iter * ds;
        odeint(Eq,param,s,s_next,ds_ctrl,&Cb);
        s=s_next;
        ++iter;
        if(fp) (*fp)("%g %g %g %g\n",param[1],param[2],s,Rad2Deg(param[3]));
        //if(!flag||s>=2*lens.R0) break;
        if(!flag) break;

        const double dz1 = sin( param[BRIDGE_A] );
        if(dz0*dz1<0)
        {
            std::cerr << "\t\textremum" << std::endl;
            break;
        }
        dz0 = dz1;
    }

    return dz0;
}