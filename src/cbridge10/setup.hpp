#ifndef SETUP_INCLUDED
#define SETUP_INCLUDED 1

#include "bridge.hpp"
#include "yocto/string.hpp"
#include "yocto/math/fit/glsf-spec.hpp"


enum Direction
{
    Pulling, //!< increasing H
    Pushing, //!< decreasing H
};

#define SETUP_EXTRACT_THETA 1

class Setup
{
public:

    explicit Setup(const double  user_R0,
                   const double  user_capillary_length);
    DefaultBridge bridge;
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

    double rebuild(const double         th,
                   array<double>       &zparam,
                   const array<double> &coord,
                   const array<double> &alpha,
                   const array<double> &zeta,
                   array<double>       &dzeta,
                   array<double>       &dzfit
                   );


    void run( threading::context &ctx, array<double> &target, const array<double> &source, const array<double> &second, void *args);
    

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Setup);
};


//! will provide a cutoff
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

