#include "bridge.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
// set up
//
////////////////////////////////////////////////////////////////////////////////

Bridge:: ~Bridge() throw()
{
}


#define GET_FROM_LUA(ID) Lua::Config::Get<lua_Number>(L,#ID)
#define INI_FROM_LUA(ID) ID(GET_FROM_LUA(ID))
Bridge:: Bridge( lua_State *L ) :
nvar(BRIDGE_N),
odeint( Lua::Config::Get<lua_Number>(L,"ftol") ),
profEq( this, & Bridge::ProfileEq ),
status(false),
INI_FROM_LUA(R0),
INI_FROM_LUA(lambda),
mu2((curvature_coeff * Square(R0)) / Square(lambda)),
mu(sqrt(mu2)),
pprev(nvar),
param(nvar),
center_v(0),
INI_FROM_LUA(resolution),
INI_FROM_LUA(angle_control),
INI_FROM_LUA(shift_control),
fn_of_alpha(this, & Bridge:: ProfileOfAlpha),
__alpha(0),
__theta(0),
__Xi(0)
{
    odeint.start(nvar);
    std::cerr << "ftol="       << odeint.eps  << std::endl;
    std::cerr << "resolution(deg)   =" << resolution << std::endl;
    (double&)resolution = Deg2Rad(resolution);
    std::cerr << "resolution(rad)   =" << resolution << std::endl;
    std::cerr << "shift_control(tau)=" << shift_control << std::endl;
    std::cerr << "angle_control(deg)=" << angle_control << std::endl;
    (double&)angle_control = Deg2Rad(angle_control);
    std::cerr << "angle_control(rad)=" << angle_control << std::endl;

    std::cerr << "R0    =" << R0 << std::endl;
    std::cerr << "lambda=" << lambda << std::endl;
    std::cerr << "mu    =" << mu << std::endl;



}

size_t Bridge:: curvature_coeff = 1;



#include "yocto/ios/ocstream.hpp"
void Bridge:: SaveLens(const string &filename, const double shift )
{
    ios::wcstream fp(filename);
    const double vc = 1.0 + shift;
    for(double alpha=0;alpha<=360;alpha+=0.1)
    {
        const double angle = Deg2Rad(alpha);
        fp("%.15g %.15g\n", sin(angle), vc - cos(angle));
    }
}




