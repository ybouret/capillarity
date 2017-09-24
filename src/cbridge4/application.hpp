#ifndef APPLICATION_INCLUDED
#define APPLICATION_INCLUDED 1


#include "yocto/threading/scheme/server.hpp"
#include "bridge.hpp"

typedef threading::context Context;
typedef threading::kernel  Job;

#define __CORE 0x01
#define __SUBS 0x02

class Application : public Bridge, public threading::par_server
{
public:
    virtual ~Application() throw();
    explicit Application( Lua::State &L );

    Vector h;       //! original recording
    Vector A;       //! original recording
    Vector t;       //!< rebuilt time
    Vector v;       //!< rebuilt rate
    Vector delta_t; //!< rebuilt dt
    Vector h_evap;  //!< with evaporation
    Vector h_corr;  //!< corrected with push/pull
    Vector zeta;
    Vector alpha;
    Vector theta;

    double main_rate; //!< main rate in mm/s (around 5-e3)
    double evap_rate; //!< evap rate in mm/s
    double percent;   //!< percent correction
    
    // during pull, corrected height w.r.t h_evap
    inline double parabole( const double x )
    {
        return (-0.0185)*(x)+(0.0345238)*(x*x);
        //return (-0.171752)+(-0.0185)*(x)+(0.0345238)*(x*x);
    }

    void load(const string &filename);

    void load_v2(const string &dirname);


    void build_time();


    void correct_h();

    void pre_process();


    void process()
    {
        flush();
        std::cerr << "enqueuing #" << h.size() << " jobs..." << std::endl;
        for(size_t i=1;i<=h.size();++i)
        {
            enqueue_copy(i,this, &Application::run1);
        }
        flush();
    }


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Application);
    container_manager_of<Vector> mgr;

    void run1(const size_t i, threading::context &ctx)
    {
        assert(i>0);
        assert(i<=theta.size());
        
        Bridge &B      = ctx.as<Bridge>();
        bool    isFlat = false;

        theta[i]  = B.find_theta(alpha[i], zeta[i], isFlat);
    }

};

#endif



