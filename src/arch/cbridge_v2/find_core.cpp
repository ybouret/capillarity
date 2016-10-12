

Lua::State VM;
lua_State *L = VM();

// default parameters
Lua::Config::DoString(L,"ftol=1e-7;");
Lua::Config::DoString(L,"angle_control=1;");
Lua::Config::DoString(L,"shift_control=0.01;");
Lua::Config::DoString(L,"R0=80;");
Lua::Config::DoString(L,"lambda=2.72;");
Lua::Config::DoString(L,"Xi=0.0;");
Lua::Config::DoString(L,"resolution=1e-3");

// modified
for(int i=1;i<argc;++i)
{
    Lua::Config::DoString(L,argv[i]);
}

// load bridge
Bridge B(L,1);

#if !defined(DISCARD_XI)
// load other parameters
const double Xi = Lua::Config::Get<lua_Number>(L,"Xi");
std::cerr << "Xi=" << Xi << std::endl;
#endif
