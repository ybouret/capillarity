#include "setup.hpp"

////////////////////////////////////////////////////////////////////////////////
Setup:: ~Setup() throw()
{
}

Setup:: Setup(const double user_R0,
              const double user_capillary_length) :
bridge(),
R0(user_R0),
capillary_length(user_capillary_length),
R02( Square(R0) ),
S0( R02 * numeric<double>::pi )
{
    assert(R0>0);
    assert(capillary_length>0);
    bridge.set_mu(R0,capillary_length);
}


size_t Setup:: isolate(const Direction      hdir,
                       vector<double>      &zeta,
                       vector<double>      &alpha,
                       const array<double> &height,
                       const array<double> &surface,
                       const double         CutOff)
{
    const size_t N0 = height.size();
    zeta.free();
    alpha.free();
    const double   S0 = numeric<double>::pi * R0*R0;
    for(size_t i=1;i<=N0;++i)
    {
        const double h_i = height[i];
        const double hh  = h_i/R0;
        const double ss  = surface[i]/S0;
        if(ss>=1)
        {
            throw exception("invalid surface %g", surface[i]);
        }

        const double aa = asin( sqrt(ss) );
        switch(hdir)
        {
            case Pushing:
                if(h_i<=CutOff)
                {
                    zeta.push_back(hh);
                    alpha.push_back(aa);
                }
                break;

            case Pulling:
                if(h_i>=CutOff)
                {
                    zeta.push_back(hh);
                    alpha.push_back(aa);
                }
                break;
        }
    }
    return zeta.size();
}

double Setup:: rebuild(const double         th,
                       array<double>       &zparam,
                       const array<double> &coord,
                       const array<double> &alpha,
                       const array<double> &zeta,
                       array<double>       &dzeta,
                       array<double>       &dzfit
                       )
{
    const size_t N = coord.size();
    for(size_t i=1;i<=N;++i)
    {
        dzeta[i] = bridge.find_zeta(alpha[i], th) - zeta[i];
    }

    GLS<double>::Sample S(coord,dzeta,dzfit);
    _GLS::Polynomial<double>::Start(S,zparam);
    return Fabs(zparam[2]);
}

////////////////////////////////////////////////////////////////////////////////

Parted:: Parted() throw()
{}

Parted:: ~Parted() throw() {}


//! decreasing heigh
double Parted:: push( double x, const array<double> &aorg )
{
    assert(aorg.size()>=3);
    const double intercept = aorg[1];
    const double slope     = aorg[2];
    const double cutoff    = aorg[3];
    const double surfmin   = intercept + slope * cutoff;
    return (x <= cutoff) ? intercept+slope*x : surfmin;
}

//! increasing h
double Parted:: pull( double x, const array<double> &aorg )
{
    assert(aorg.size()>=3);
    const double intercept = aorg[1];
    const double slope     = aorg[2];
    const double cutoff    = aorg[3];
    const double surfmax   = intercept + slope * cutoff;
    return (x >= cutoff) ? intercept+slope*x : surfmax;
}

bool Parted:: callback(const GLS<double>::Samples &, const array<double> &aorg )
{
    std::cerr << "-- cutting at " << aorg[3] << std::endl;
    return true;
}

#include "yocto/math/fit/glsf-spec.hpp"

double Parted:: split( const Direction hdir, const array<double> &height, const array<double> &surface, array<double> &surffit)
{
    GLS<double>::Samples       samples;
    const GLS<double>::Sample &sample = samples.append(height,surface,surffit);

    const size_t   N0   = height.size();
    const size_t   nvar = 3;
    vector<double> aorg(nvar);
    vector<double> aerr(nvar);
    vector<bool>   used(nvar,true);
    {
        vector<double> p(2);
        if(!_GLS::Polynomial<double>::Start(sample, p))
        {
            throw exception("couldn't guess initial parameters");
        }
        std::cerr << "p=" << p << std::endl;
        aorg[1] = p[1];
        aorg[2] = p[2];
        aorg[3] = 0.5*(height[1]+height[N0]);
    }
    samples.prepare(nvar);

    GLS<double>::Callback  cb(this, & Parted::callback );
    GLS<double>::Function  *part = 0;
    GLS<double>::Function  _push(this, & Parted:: push );
    GLS<double>::Function  _pull(this, & Parted:: pull );

    switch(hdir)
    {
        case Pushing: part = & _push; break;
        case Pulling: part = & _pull; break;
    }

    if( !samples.fit_with(*part, aorg, used, aerr, &cb) )
    {
        throw exception("couldn't part data...");
    }
    GLS<double>::display(std::cerr, aorg, aerr);
    return aorg[3];
}

void Setup:: run(threading::context  &ctx,
                 array<double>       &target,
                 const array<double> &source,
                 const array<double> &second,
                 void                *args)
{
    assert(args);
    const int choice = *(int *)args;
    switch(choice)
    {
        case SETUP_EXTRACT_THETA:
        {
            array<double>       &theta = target;
            const array<double> &alpha = source;
            const array<double> &zeta  = second;
            size_t offset = 1;
            size_t length = theta.size();
            ctx.split(offset, length);
            for(size_t i=offset,count=length;count-->0;++i)
            {
                theta[i] = bridge.find_theta(alpha[i], zeta[i]);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
Setups:: ~Setups() throw()
{

}

Setups:: Setups(threading::kernel_executor *kxp,
                const double                R0,
                const double                capillary_length) :
threading::processing_unit<Setup>(kxp)
{
    for(size_t i=0;i<cores;++i)
    {
        append<double,double>(R0,capillary_length);
    }
}


////////////////////////////////////////////////////////////////////////////////

Application:: ~Application() throw() {}

Application:: Application( threading::kernel_executor *kxp,
                          const double                 user_R0,
                          const double                 user_capillary_length) :
Setups(kxp,user_R0,user_capillary_length),
setup( (*this)[0] )
{

}

void Application:: extract_theta()
{
    int choice = SETUP_EXTRACT_THETA;
    compile<double,double,double>();
    call(theta, alpha, zeta, &choice);
}

