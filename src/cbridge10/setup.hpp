#ifndef SETUP_INCLUDED
#define SETUP_INCLUDED 1

#include "bridge.hpp"
#include "yocto/string.hpp"
#include "yocto/math/fit/glsf-spec.hpp"
#include "yocto/threading/vpu.hpp"


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
    const double R0;
    const double capillary_length;
    const double R02;
    const double S0;
    
    virtual ~Setup() throw();
    
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

    //! compute zeta from alpha if theta>=0, compute surf from zeta if theta<0
    void run( threading::context &ctx, array<double> &target, const array<double> &source, void *args);

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Setup);
};


class Setups : public threading::processing_unit<Setup>
{
public:
    explicit Setups(threading::kernel_executor *kxp,
                    const double                R0,
                    const double                capillary_length);

    virtual ~Setups() throw();


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Setups);
};

class Application : public Setups
{
public:

    explicit Application( threading::kernel_executor *kxp,
                         const double                 R0,
                         const double                 capillary_length);

    virtual ~Application() throw();

    Setup &setup;

    vector<double> zeta;
    vector<double> alpha;
    vector<double> theta;
    vector<double> dzeta;
    vector<double> znew;
    vector<double> zfit;
    
    void extract_theta();

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Application);
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

