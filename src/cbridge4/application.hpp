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


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Application);
};

#endif



