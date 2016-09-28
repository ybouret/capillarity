#include "bridge.hpp"

double Bridge:: ProfileOfAlpha(const double alpha)
{
    return profile(alpha, __theta, __zeta, NULL);
}
