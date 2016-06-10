#include "setup.hpp"

////////////////////////////////////////////////////////////////////////////////
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

    GLS<double>::Callback cb(this, & Parted::callback );
    GLS<double>::Function *part = 0;
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

