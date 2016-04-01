#ifndef BRIDGE_INCLUDED
#define BRIDGE_INCLUDED 1

#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "lens.hpp"


#define BRIDGE_R   1
#define BRIDGE_Z   2
#define BRIDGE_PSI 3

typedef ode::Field<double>::Array    Array;
typedef ode::Field<double>::Equation Equation;
typedef ode::Field<double>::Callback Callback;

class Bridge
{
public:
    const size_t      nvar;
    bool              flag;
    ode::driverCK<double>::type  odeint;
    Equation          Eq;
    Callback          Cb;
    
    explicit Bridge();
    virtual ~Bridge() throw();


    void initialize();
    


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
    void __Eq(Array &,const double,const Array);
    void __Cb(Array &,const double);

};


#endif
