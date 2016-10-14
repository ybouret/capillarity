#ifndef APPLICATION_INCLUDED
#define APPLICATION_INCLUDED 1


#include "yocto/threading/engine.hpp"
#include "bridge.hpp"

typedef threading::context Context;
typedef threading::job     Job;

class Application : public Bridge, public threading::engine
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
            //enqueue2(this, &Application::run1, i);
        }
        flush();
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Application);
    container_manager_of<Vector> mgr;

    void run1(size_t i, threading::context &ctx)
    {
        {
            YOCTO_LOCK(ctx.access);
            std::cerr << "i=" << i << std::endl;
            std::cerr << "theta.size()=" << theta.size() << std::endl;
        }
        Bridge &B = ctx.as<Bridge>();
        //assert(i>0);
        //assert(i<=theta.size());
        //theta[i] = B.mu;
    }

};

#endif



