#include "bridge.hpp"

Bridge:: Bridge() :
nvar(3),
flag(true),
rk4(),
Q(rk4.YY),
Eq( this, &Bridge::Eq ),
Cb( this, &Bridge::Cb )
{
    rk4.allocate(nvar);
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
