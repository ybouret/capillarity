#include "par-bridge.hpp"
#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
    std::cerr << "-- Reading Parameters" << std::endl;
    Lua::State VM;
    lua_State *L = VM();

    // default parameters
    Lua::Config::DoString(L,"ftol=1e-6;");
    Lua::Config::DoString(L,"angle_control=0.5;");
    Lua::Config::DoString(L,"shift_control=0.001;");
    Lua::Config::DoString(L,"resolution=1e-4");

    // user's parameter
    for(int i=1;i<argc;++i)
    {
        Lua::Config::DoString(L,argv[i]);
    }

    std::cerr << "-- Loading Data" << std::endl;
    ParBridge    B(L);
    const string filename = Lua::Config::Get<string>(L,"file");

    B.load(filename);
    std::cerr << "-- Processing..." << std::endl;
    B.find_theta();

    string root = vfs::get_base_name(filename);
    vfs::change_extension(root,"theta.dat");
    std::cerr << "-- Saving into " << root << std::endl;
    
    {
        ios::wcstream fp(root);
        for(size_t i=1;i<=B.zeta.size();++i)
        {
            fp("%.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
               B.zeta[i]*B.R0,
               B.surf[i]*B.A0,
               B.times[i],
               B.zeta[i]*B.R0 + B.times[i] * B.rate,
               Rad2Deg(B.alpha[i]),
               Rad2Deg(B.theta[i]),
               Rad2Deg(B.theta2[i])
               );
        }
    }

}
YOCTO_PROGRAM_END()

