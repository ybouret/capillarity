#include "bridge.hpp"
#include "yocto/exception.hpp"

size_t Bridge::NUM_STEPS = 10000;
double Bridge::ATOL      = 1e-5;
double Bridge::HTOL      = 1e-5;

Bridge:: ~Bridge() throw()
{

}

Bridge:: Bridge(const Lens  &usr_lens,
                const double usr_clength ) :
lens( usr_lens.clone() ),
capillary_length(usr_clength),
arrays(10),
Y(  arrays.next_array() ),
k1( arrays.next_array() ),
k2( arrays.next_array() ),
k3( arrays.next_array() ),
k4( arrays.next_array() ),
V(  arrays.next_array() )
{
    if( capillary_length <= 0 )
    {
        throw exception("negative capillary length");
    }
    lens->initialize();

    std::cerr << "\t(*) Bridge Initialized" << std::endl;


}
