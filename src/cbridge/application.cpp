#include "application.hpp"

Application:: ~Application() throw()
{
}

Application:: Application() :
R0(0),
cap_len(0),
A0(0)
{
    
}

Bridge & Application:: initialize(KernelExecutor &kExec,const double user_R0, const double user_capillary_length)
{
    R0      = user_R0;
    cap_len = user_capillary_length;
    A0      = numeric<double>::pi * Square(R0);
    for(size_t i=0;i<kExec.num_threads();++i)
    {
        Context &ctx    = kExec[i];
        Bridge  &bridge = ctx.make<DefaultBridge>();
        bridge.set_mu(R0,cap_len);
    }
    return kExec[0].as<Bridge>();
}


void Application:: build_reduced_variables()
{
    for(size_t i=h.size();i>0;--i)
    {
        zeta[i] = h[i]/R0;
        const double ss = Sqrt(A[i]/A0);
        if( Fabs(ss) > 1 )
        {
            throw exception("invalid area in build_reduced_variables");
        }
        alpha[i] = Asin(ss);
    }

}



void Application:: load_from(const string &filename)
{
    
}
