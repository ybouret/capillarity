#include "application.hpp"

Application:: ~Application() throw()
{
}

Application:: Application() :
R0(0),
cap_len(0),
A0(0),
A(),h(),zeta(),alpha(),
vecs(),
subs()
{
    vecs.enroll(A);
    vecs.enroll(h);
    vecs.enroll(zeta);    subs.enroll(zeta);
    vecs.enroll(alpha);   subs.enroll(alpha);
}

Bridge & Application:: initialize(KernelExecutor &kExec,const double user_R0, const double user_capillary_length)
{
    std::cerr << "-- initialize Bridges for threads..." << std::endl;
    R0      = user_R0;
    cap_len = user_capillary_length;
    A0      = numeric<double>::pi * Square(R0);
    for(size_t i=0;i<kExec.num_threads();++i)
    {
        Context &ctx    = kExec[i];
        Bridge  &bridge = ctx.make<DefaultBridge>();
        bridge.set_mu(R0,cap_len);
    }
    std::cerr << "-- done" << std::endl;
    return kExec[0].as<DefaultBridge>();
}


void Application:: build_reduced_variables()
{
    for(size_t i=h.size();i>0;--i)
    {
        zeta[i] = h[i]/R0;
        const double ss = Sqrt(A[i]/A0);
        if( Fabs(ss) > 1 )
        {
            throw exception("invalid area in build_reduced_variables");
        }
        alpha[i] = Asin(ss);
    }

}

void Application:: compute_theta_using(KernelExecutor &kExec)
{
    Kernel K(this, & Application::compute_theta_kernel);
    kExec(K);
}

#define COMPUTE_THETA(I)


void Application:: compute_theta_kernel( Context &ctx )
{
    size_t offset = 1;
    size_t length = h.size();
    ctx.split(offset,length);
    YOCTO_LOOP_FUNC(length,COMPUTE_THETA,offset);
}

#include "yocto/math/io/data-set.hpp"

void Application:: load_from(const string &filename)
{
    vecs.free_all();
    ios::icstream fp(filename);
    data_set<double> ds;
    ds.use(1,A);
    ds.use(2,h);
    ds.load(fp);
    subs.make_all(A.size());
}
