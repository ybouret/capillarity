#include "par_bridge.hpp"

ParBridge:: ~ParBridge() throw()
{
}

ParBridge:: ParBridge(lua_State *L) :
Crew(true)
{
    for(size_t i=0;i<size;++i)
    {
        Context &ctx = (*this)[i];
        ctx.build<Bridge,lua_State *>(L);
    }
}
