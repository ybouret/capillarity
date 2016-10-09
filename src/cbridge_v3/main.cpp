#include "par-bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
#include "main-core.cpp"
    
    ParBridge Bridges(L);
    
}
YOCTO_PROGRAM_END()

