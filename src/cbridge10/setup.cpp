#include "setup.hpp"

Setup:: ~Setup() throw()
{
}

Setup:: Setup(const double  user_R0,
              const double user_capillary_length) :
bridge(BRIDGE_SEARCH_DEGREES,
       BRIDGE_INTEGRATOR_FTOL,
       BRIDGE_INTEGRATOR_DEGREES,
       BRIDGE_INTEGRATOR_LENGTH),
R0(user_R0),
capillary_length(user_capillary_length)
{
    assert(R0>0);
    assert(capillary_length>0);
    
}


double Setup:: compute_theta(const double alpha, const double zeta)
{
    bridge.mu = sqrt( Square(R0)/(2*Square(capillary_length) ));
    //std::cerr << "mu=" << bridge.mu << std::endl;

    return bridge.find_theta(alpha,zeta);
}
