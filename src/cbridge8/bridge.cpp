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

    arrays.allocate(2);

    std::cerr << "\t(*) Bridge Initialized" << std::endl;


}

#include "yocto/code/rand.hpp"

void Bridge:: Tests()
{
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
            std::cerr << " |_alpha=" << a << std::endl;
            const double t = FindTheta(h,a);
            std::cerr << "   |_theta=" << t << std::endl;

            if(RInt(t)!=theta)
            {
                throw exception("Wrong inversion!");
            }

        }
    }
}

#include "yocto/threading/window.hpp"

void Bridge:: TestsMT(const threading::context &ctx ) throw()
{
    const threading::window win(ctx,35,1);
    {
        scoped_lock guard(ctx.access);
        Bridge *self = this;
        std::cerr << "Bridge@" << (void*)self << "/" << (void*)&ctx << " : " << win.start << " -> " << win.final << std::endl;
    }

    for(size_t i=win.start;i<=win.final;++i)
    {
        const int theta = i*5;
        {
            scoped_lock guard(ctx.access);
            std::cerr << "theta=" << theta << std::endl;
        }
        const double hmax = FindHmax(theta);
        for(int i=1;i<=10;++i)
        {
            const double h = alea<double>() * hmax;
            const double a = FindAlpha(h,theta);
            const double t = FindTheta(h,a);

            if(RInt(t)!=theta)
            {
                scoped_lock guard(ctx.access);
                std::cerr << "wrong inversion for theta=" << theta << ", h=" << h  << std::endl;
            }
            
        }

    }
}


void Bridge:: CallTests(threading::context &ctx) throw()
{
    Bridge &bridge = ctx.as<Bridge>();
    bridge.TestsMT(ctx);
}