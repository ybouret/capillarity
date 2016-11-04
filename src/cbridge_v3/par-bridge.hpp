#ifndef PAR_BRIDGE_INCLUDED
#define PAR_BRIDGE_INCLUDED 1

#include "cbridge.hpp"
#include "yocto/threading/scheme/simd.hpp"

typedef threading::context Context;
typedef threading::SIMD    Crew;
typedef container_manager_of<Vector> VMgr;

//! base class bridge to have a copy of global parameters
class ParBridge : public Bridge, public Crew
{
public:
    explicit ParBridge( lua_State *L );
    virtual ~ParBridge() throw();


    Vector zeta;   // h/R0
    Vector surf;   // area/A0
    Vector times;  // units ?
    Vector alpha;  // asin( sqrt(surf) )
    Vector theta;  // theta..
    Vector theta2; // with rate...
    
    void load( const string &filename );


    void find_theta(); //!< in parallel
    double rate;

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(ParBridge);
    VMgr mgr;

    void FindTheta( Context &ctx );


};

#endif


