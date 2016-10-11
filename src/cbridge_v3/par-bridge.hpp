#ifndef PAR_BRIDGE_INCLUDED
#define PAR_BRIDGE_INCLUDED 1

#include "cbridge.hpp"
#include "yocto/threading/crew.hpp"

typedef threading::context Context;
typedef threading::crew    Crew;
typedef container_manager_of<Vector> VMgr;

//! base class bridge to have a copy of global parameters
class ParBridge : public Bridge, public Crew
{
public:
    explicit ParBridge( lua_State *L );
    virtual ~ParBridge() throw();


    Vector zeta;   // h/R0
    Vector surf;   // area/A0
    Vector alpha;  // asin( sqrt(surf) )
    Vector theta;  // theta..
    Vector theta2; // theta
    Vector Theta;  //
    Vector Theta2;
    void load( const string &filename );


    void find_theta(); //!< in parallel
    double rate;

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(ParBridge);
    VMgr mgr;

    void FindTheta( Context &ctx );


};

#endif


