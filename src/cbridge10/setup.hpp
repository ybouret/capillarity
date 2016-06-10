#ifndef SETUP_INCLUDED
#define SETUP_INCLUDED 1

#include "bridge.hpp"
#include "yocto/string.hpp"
#include "yocto/math/fit/glsf.hpp"

#define BRIDGE_SEARCH_DEGREES     0.001
#define BRIDGE_INTEGRATOR_FTOL    1e-5
#define BRIDGE_INTEGRATOR_DEGREES 0.1
#define BRIDGE_INTEGRATOR_LENGTH  0.01

enum Direction
{
    Pulling, //!< increasing H
    Pushing, //!< decreasing H
};


class Setup
{
public:

    explicit Setup(const double  user_R0,
                   const double  user_capillary_length);
    Bridge bridge;
    double R0;
    double capillary_length;

    virtual ~Setup() throw();

    double compute_theta(const double alpha, const double zeta);

    size_t isolate(const Direction      hdir,
                   vector<double>      &zeta,
                   vector<double>      &alpha,
                   const array<double> &height,
                   const array<double> &surface,
                   const double         CutOff);


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Setup);
};



class Parted
{
public:
    explicit Parted() throw();
    virtual ~Parted() throw();
    double push( double x, const array<double> &aorg );
    double pull( double x, const array<double> &aorg );
    bool   callback(const GLS<double>::Samples &, const array<double> &aorg );
    double split(const Direction hdir,const array<double> &height, const array<double> &surface, array<double> &surffit);


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parted);
};

#endif

