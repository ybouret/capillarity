#ifndef BRIDGE_INCLUDED
#define BRIDGE_INCLUDED 1

#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "lens.hpp"
#include "yocto/ios/ostream.hpp"


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
    double                       current_alpha;    //!< for computation
    double                       current_theta;    //!< for computation
    Lens                        *current_lens;     //!< for computations
    ios::ostream                *current_fp;       //!< during computation
    vector<double>               param;
    vector<double>               pprev;
    ode::driverCK<double>::type  odeint;
    Equation                     Eq;
    Callback                     Cb;
    Function                     FnOfAlpha;
    Function                     FnOfTheta;
    const double                 delta; //!< angle delta, in radiants

    explicit Bridge(const double delta_deg);
    virtual ~Bridge() throw();


    double compute_profile(Lens        &lens,
                           const double alpha,
                           const double theta,
                           const double height,
                           ios::ostream *fp);


    double FindTheta( Lens &lens, const double alpha, const double height );


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
    void __Eq(Array &,const double,const Array&);
    void __Cb(Array &,const double); // testing valid position
    double __profile_of_alpha( const double alpha );
    double __profile_of_theta( const double theta );

    double __compute_profile();

};


#endif
