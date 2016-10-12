#ifndef APPLICATION_INCLUDED
#define APPLICATION_INCLUDED 1

#include "bridge.hpp"
#include "yocto/threading/kernel-executor.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/container/manager.hpp"

typedef threading::kernel             Kernel;
typedef threading::kernel_executor    KernelExecutor;
typedef threading::context            Context;
typedef vector<double>                Vector;
typedef container_manager_for<Vector> VecMgr;

class Application
{
public:
    explicit Application();
    virtual ~Application() throw();

    Bridge &initialize(KernelExecutor &kExec,const double user_R0,const double user_capillary_length);

    double         R0;
    double         cap_len;
    double         A0;      //!< pi * R0^2 in mm2
    double         rate;
    
    vector<double> A;       //!< area in mm2
    vector<double> h;       //!< height in mm
    vector<double> zeta;    //!< h/R0
    vector<double> alpha;   //!< alpha=asin( sqrt(A/A0) );
    vector<double> theta;
    VecMgr         vecs; //A,h,zeta,alpha
    VecMgr         subs; //!< zeta, alpha to setup memory

    void build_reduced_variables();
    void load_from(const string &filename);
    void compute_theta_using(KernelExecutor &kExec);


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Application);
    void compute_theta_kernel( Context & );

};

#endif
