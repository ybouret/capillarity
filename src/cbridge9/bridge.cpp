#include "bridge.hpp"

Bridge:: Bridge() :
nvar(3),
flag(true),
odeint(1e-7),
Eq( this, &Bridge::Eq ),
Cb( this, &Bridge::Cb )
{
    odeint.start(nvar);
}

Bridge:: ~Bridge() throw()
{
}


void Bridge:: initialize()
{
    flag = true;
    
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
