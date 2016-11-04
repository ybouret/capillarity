#ifndef APPLICATION_INCLUDED
#define APPLICATION_INCLUDED 1


#include "yocto/threading/scheme/server.hpp"
#include "bridge.hpp"

typedef threading::context Context;
typedef threading::kernel  Job;

class Application : public Bridge, public threading::par_server
{
public:
    virtual ~Application() throw();
    explicit Application( lua_State *L );

    Vector h;
    Vector A;
    Vector t;
    Vector zeta;
    Vector alpha;
    Vector theta;

    void load(const string &filename);

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



