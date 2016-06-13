#include "bridge.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    DefaultBridge B;

    const double R0                = strconv::to<double>(argv[1],"R0");
    const double capillary_length  = strconv::to<double>(argv[2],"mu");
    B.set_mu(R0, capillary_length);
    
    


}
YOCTO_PROGRAM_END()
