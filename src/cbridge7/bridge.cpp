#include "bridge.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/kernel/tao.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/ptr/auto.hpp"
#include "yocto/ios/ocstream.hpp"

size_t Bridge::NUM_STEPS = 10000;

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
}

void Bridge::Evaluate(array_t &dYdz, double z, const array_t &Y) throw()
{
    const double r     = Y[1];
    const double drdz  = Y[2];
    const double speed = Hypotenuse(1.0,drdz);
    const double accel = (z/Square(capillary_length)+1.0/r/speed)*speed*speed*speed;
    dYdz[1] = drdz;
    dYdz[2] = accel;
}

void Bridge:: RK4(const double z_ini,const double z_end) throw()
{
    const double h     = z_end - z_ini;
    const double half  = h*0.5;
    const double z_mid = z_ini + half;
    Evaluate(k1,z_ini,Y);

    tao::setprobe(V, Y, half, k1);
    Evaluate(k2,z_mid,V);

    tao::setprobe(V,Y,half,k2);
    Evaluate(k3,z_mid,V);

    tao::setprobe(V,Y,h,k3);
    Evaluate(k4,z_end,V);

    for(size_t i=Y.size();i>0;--i)
    {
        Y[i] += h * (k1[i]+k2[i]+k2[i]+k3[i]+k3[i]+k4[i]) / 6.0;
    }

}


bool Bridge:: FinalRadius(const double height, const double theta_deg, const double alpha_deg, const bool save)
{
    static vfs & fs = local_fs::instance();
    if(save)
    {
        fs.try_remove_file("profile.dat");
    }

    // check possibility
    const double theta = Deg2Rad(theta_deg);
    const double alpha = Deg2Rad(alpha_deg);
    const double omega = lens->omega(alpha);
    const double angle = omega+theta;
    if(angle>=numeric<double>::pi)
    {
        return false;
    }


    // prepare initial conditions
    const double Rc   = lens->radius(0);
    const double R0   = lens->radius(alpha);
    const double z0   = height + Rc - R0 * Cos(alpha);
    const double r0   = R0 * Sin(alpha);
    const double drdz = Cos(angle)/Sin(angle);
    Y[1] = r0;
    Y[2] = drdz;

    auto_ptr<ios::ostream> fp;
    if(save)
    {
        fp.reset(new ios::ocstream("profile.dat",false) );
        (*fp)("%g %g\n", r0, z0 );
    }

    double reg[3] = { Y[1],0,0 }; // curvature detector
    size_t num    = 1;
    for(size_t i=NUM_STEPS;i>0;--i)
    {
        const double z_ini = (i*z0)/NUM_STEPS;
        const double z     = ((i-1)*z0)/NUM_STEPS;
        RK4(z_ini, z);

        // analyze
        const double r = Y[0];
        

        if(save && (r<=4*Rc))
        {
            (*fp)("%g %g\n", r, z );
        }
    }



    return false;
}