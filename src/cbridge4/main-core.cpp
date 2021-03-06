

Lua::State L;

// default parameters
L.DoString("ftol=1e-6;");
L.DoString("angle_control=0.5;");
L.DoString("shift_control=0.001;");
L.DoString("R0=82;");
L.DoString("lambda=2.72;");
L.DoString("resolution=1e-4");
L.DoString("evap_rate=0");
L.DoString("coef_push=0");
L.DoString("coef_pull=0");
L.DoString("main_rate=5e-3");



// modified
for(int i=2;i<argc;++i)
{
    L.DoString(argv[i]);
}

// load bridge
#if !defined(DISCARD_BRIDGE)
Bridge B(L);
#endif

