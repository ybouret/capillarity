#include "yocto/program.hpp"
#include "bridge.hpp"

YOCTO_PROGRAM_START()
{
    const Lens::Pointer lens( new SphericalLens(70) );
    
}
YOCTO_PROGRAM_END()
