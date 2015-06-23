#include "bridge.hpp"
#include "yocto/exception.hpp"

size_t Bridge::NUM_STEPS = 10000;
double Bridge::ATOL      = 1e-5;
double Bridge::HTOL      = 1e-4;

Bridge:: ~Bridge() throw() {}

Bridge:: Bridge(const Lens::Pointer &usr_lens,const double usr_clength) :
lens(usr_lens),
capillary_length(usr_clength),
arrays(10),
Y(  arrays.next_array() ),
k1( arrays.next_array() ),
k2( arrays.next_array() ),
k3( arrays.next_array() ),
k4( arrays.next_array() ),
V(  arrays.next_array() )
{
    if(capillary_length<=0) throw exception("negative capillary length");
    arrays.allocate(2);
    assert(2==Y.size());
}

