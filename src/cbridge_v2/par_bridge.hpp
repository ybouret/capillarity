#ifndef PAR_BRIDGE_INCLUDED
#define PAR_BRIDGE_INCLUDED 1

#include "yocto/threading/crew.hpp"
#include "bridge.hpp"

typedef threading::context Context;
typedef threading::crew    Crew;

class ParBridge : public Crew
{
public:
    explicit ParBridge(lua_State *L);
    virtual ~ParBridge() throw();

    void extract_from(lua_State *L);

    vector<double> area;   //!< mm^2
    vector<double> height; //!< mm

    vector<double> alpha;  //!< angle, radians
    vector<double> zeta;   //!< h/R0
    vector<double> theta;  //!< theta..



private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(ParBridge);
    void run( Context &ctx );

};

#endif
