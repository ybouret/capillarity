#ifndef BRIDGE_INCLUDED
#define BRIDGE_INCLUDED 1

#include "lens.hpp"
#include "yocto/ptr/auto.hpp"
#include "yocto/sequence/many-arrays.hpp"
#include "yocto/threading/crew.hpp"

using namespace threading;
typedef threading::context                 Context;
typedef array<double>                      array_t;
typedef many_arrays<double,memory::global> arrays_t;

class Bridge
{
public:
    static size_t NUM_STEPS;
    static double ATOL;
    static double HTOL;

    auto_ptr<Lens> lens;
    const double   capillary_length;
    arrays_t       arrays;
    array_t       &Y;
    array_t       &k1;
    array_t       &k2;
    array_t       &k3;
    array_t       &k4;
    array_t       &V;
    double        param1;
    double        param2;
    double        result;
    void         *args;   // datafile...

    explicit Bridge(const Lens &usr_lens, const double clength);
    virtual ~Bridge() throw();
    
    //! Y=[r drdz]
    void Evaluate( array_t &dYdz, double z, const array_t &Y ) throw();

    //! make a step
    void RK4(const double z_ini,const double z_end) throw();

    //! warning, angles are in degrees
    /**
     no throw if no save...
     */
    bool FinalRadius(const double height,
                     const double theta,
                     const double alpha,
                     const bool   save=false);


    //! warning, theta is in degree
    double FindAlpha(const double height,const double theta) throw();

    //! theta in degrees
    double FindHmax(const double theta) throw();

    //! scanning alpha for algo debug
    void ScanAlpha(const double height, const double theta);

    //! warning, alpha is in degree, return as well
    double FindTheta(const double height, const double alpha) throw();


    void Tests(const Context &) throw();
    static void CallTests( Context &) throw();
    
    //! result = FindHmax(param1)
    static void CallHmax( Context &) throw();

    //! result = FindAlpha(height=param1,theta=param2);
    static void CallAlpha( Context &) throw();

    //! process data (address in args)
    void Process(const Context &ctx) throw();
    static void CallProcess(Context &ctx) throw();
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);

};

#endif

