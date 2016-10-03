#ifndef PAR_BRIDGE_INCLUDED
#define PAR_BRIDGE_INCLUDED 1

#include "yocto/threading/crew.hpp"
#include "bridge.hpp"

typedef threading::context Context;
typedef threading::crew    Crew;

class ParBridge : public Crew
{
public:
    explicit ParBridge(lua_State *L);
    virtual ~ParBridge() throw();

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(ParBridge);
};

#endif
