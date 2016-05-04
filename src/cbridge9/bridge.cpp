#include "bridge.hpp"

Bridge:: Bridge() :
nvar(3),
flag(true),
capillary_length(1),
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
    const double z   = Q[BRIDGE_Z];
    const double phi = Q[BRIDGE_A];
    

    const double C = cos(phi);
    const double S = sin(phi);

    const double l2 = capillary_length*capillary_length;
    dQds[BRIDGE_R] = C; // dr/ds
    dQds[BRIDGE_Z] = S; // dz/ds
    dQds[BRIDGE_A] = -(z/l2+S/r); // dphi/ds

}


void Bridge:: __Cb(Array &Qtmp,const double s)
{
    const double r = Qtmp[BRIDGE_R];
    if(r<=0)
    {
        flag=false;
        return;
    }
}

bool Bridge:: compute_profile(Lens        &lens,
                              const double alpha,
                              const double theta,
                              const double height)
{
    flag = true;

    lens.starting_point(param, alpha, theta, height);

    const double ds      = capillary_length/100.0;
    double       ds_ctrl = ds/10.0;
    double       s       = 0;
    size_t       iter    = 1;
    while(true)
    {
        double s_next = iter * ds;
        odeint(Eq,param,s,s_next,ds_ctrl,NULL);
        s=s_next;
        ++iter;
        break;
    }

    return true;
}