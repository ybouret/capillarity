#include "par-bridge.hpp"

ParBridge:: ~ParBridge() throw()
{
    
}

ParBridge:: ParBridge(lua_State *L) : Crew(true)
{
    Crew &self = *this;
    for(size_t i=0;i<self.size;++i)
    {
        self[i].build<Bridge,lua_State *>(L);
    }
}
