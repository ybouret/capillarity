#ifndef WINDOW_INCLUDED
#define WINDOW_INCLUDED

#include "yocto/threading/window.hpp"
#include "datafile.hpp"

class Window : public threading::window
{
public:
    explicit Window(const threading::context &ctx,  DataFile &df );
    virtual ~Window() throw();
    
    DataFile &data;

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Window);
};

#endif

