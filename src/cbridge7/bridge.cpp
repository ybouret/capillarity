#include "bridge.hpp"
#include "yocto/exception.hpp"

size_t Bridge::NUM_STEPS = 10000;
double Bridge::ATOL      = 1e-5;
double Bridge::HTOL      = 1e-5;

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
    lens->Initialize();
}

#include "yocto/code/rand.hpp"

void Bridge:: Tests()
{
    alea_init();
    for(ptrdiff_t theta=15;theta<=165;theta+=15)
    {
        std::cerr << "Test theta=" << theta << std::endl;
        std::cerr << "|_Looking for HMAX" << std::endl;
        const double hmax = FindHmax(theta);
        std::cerr << "\t\tHMAX=" << hmax << std::endl;
        for(int i=1;i<=10;++i)
        {
            const double h = alea<double>() * hmax;
            std::cerr << " |_Testing h=" << h << std::endl;
            const double a = FindAlpha(h,theta);
            std::cerr << "  |_alpha=" << a << std::endl;
            const double t = FindTheta(h,a);
            if(RInt(t)!=theta)
            {
                throw exception("Wrong inversion!");
            }

        }
    }
}
