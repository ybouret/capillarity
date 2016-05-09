#ifndef BRIDGE_INCLUDED
#define BRIDGE_INCLUDED 1

#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "lens.hpp"


typedef ode::Field<double>::Array    Array;
typedef ode::Field<double>::Equation Equation;
typedef ode::Field<double>::Callback Callback;

class Bridge
{
public:
    const size_t                 nvar;
    bool                         flag;
    double                       capillary_length;
    double                       current_height;   //!< for computations
    double                       current_center;   //!< center position
    Lens                        *current_lens;     //!< for computations
    vector<double>               param;
    ode::driverCK<double>::type  odeint;
    Equation                     Eq;
    Callback                     Cb;

    explicit Bridge();
    virtual ~Bridge() throw();


    double compute_profile(Lens        &lens,
                           const double alpha,
                           const double theta,
                           const double height);



private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
    void __Eq(Array &,const double,const Array&);
    void __Cb(Array &,const double);
    
};


#endif
