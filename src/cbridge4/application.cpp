#include "application.hpp"

Application:: Application( lua_State *L ) : Bridge(L), threading::engine(true)
{
    threading::dispatcher &self = *this;
    for(size_t i=0;i<self.num_threads();++i)
    {
        self[i].build<Bridge,lua_State *>(L);
    }
}

Application:: ~Application() throw()
{
}
