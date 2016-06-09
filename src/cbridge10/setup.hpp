#ifndef SETUP_INCLUDED
#define SETUP_INCLUDED 1

#include "bridge.hpp"
#include "yocto/string.hpp"
#include "yocto/math/fit/glsf.hpp"

#define BRIDGE_SEARCH_DEGREES     0.001
#define BRIDGE_INTEGRATOR_FTOL    1e-5
#define BRIDGE_INTEGRATOR_DEGREES 0.1
#define BRIDGE_INTEGRATOR_LENGTH  0.01

class Setup
{
public:

    explicit Setup(const double  user_R0,
                   const double  user_capillary_length);
    Bridge bridge;
    double R0;
    double capillary_length;

    virtual ~Setup() throw();

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Setup);
};

enum Direction
{
    Pulling, //!< increasing H
    Pushing, //!< decreasing H
};


class Parted
{
public:
    explicit Parted() throw() {}
    virtual ~Parted() throw() {}

    //! decreasing heigh
    double push( double x, const array<double> &aorg )
    {
        assert(aorg.size()>=3);
        const double intercept = aorg[1];
        const double slope     = aorg[2];
        const double cutoff    = aorg[3];
        const double surfmin   = intercept + slope * cutoff;
        return (x <= cutoff) ? intercept+slope*x : surfmin;
    }

    //! increasing h
    double pull( double x, const array<double> &aorg )
    {
        assert(aorg.size()>=3);
        const double intercept = aorg[1];
        const double slope     = aorg[2];
        const double cutoff    = aorg[3];
        const double surfmax   = intercept + slope * cutoff;
        return (x >= cutoff) ? intercept+slope*x : surfmax;
    }

    bool callback(const GLS<double>::Samples &, const array<double> &aorg )
    {
        std::cerr << "aorg=" << aorg << std::endl;
        return true;
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parted);
};

#endif

