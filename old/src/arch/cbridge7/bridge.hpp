#include "lens.hpp"
#include "yocto/sequence/many-arrays.hpp"
#include "yocto/string.hpp"
#include "yocto/threading/context.hpp"

typedef array<double>                      array_t;
typedef many_arrays<double,memory::global> arrays_t;

class DataFile;

#define BRIDGE_ERR_MAX_SURFACE_EXCEEDED 1
#define BRIDGE_ERR_HEIGHT_NOT_FOUND     2

class Bridge
{
public:
    static size_t NUM_STEPS;
    static double ATOL;
    static double HTOL;

    Lens::Pointer lens;
    const double  capillary_length;
    arrays_t      arrays;
    array_t      &Y;
    array_t      &k1;
    array_t      &k2;
    array_t      &k3;
    array_t      &k4;
    array_t      &V;
    double        ev_rate;   //!< evaporation rate
    int           lastResult; //!< for ExtractMT

    explicit Bridge(const Lens::Pointer &usr_lens, const double usr_clength);
    virtual ~Bridge() throw();

    //! Y=[r drdz]
    void Evaluate( array_t &dYdz, double z, const array_t &Y ) throw();

    //! make a step
    void RK4(const double z_ini,const double z_end) throw();

    //! warning, angles are in degrees
    /**
     no throw if no save...
     */
    bool FinalRadius(const double height,
                     const double theta,
                     const double alpha,
                     const bool   save=false);


    //! warning, theta is in degree
    double FindAlpha(const double height,const double theta) throw();

    //! theta in degrees
    double FindHmax(const double theta) throw();

    //! scanning alpha for algo debug
    void ScanAlpha(const double height, const double theta);

    //! warning, alpha is in degree, return as well
    double FindTheta(const double height, const double alpha) throw();


    void Tests();

    void   Process( DataFile &data, const string &savename );
    double Extract( DataFile &data );
    
    void ExtractMT(const threading::context &ctx, DataFile &data) throw();



private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
};

class Optimizer
{
public:
    explicit Optimizer(Bridge &B, DataFile &D);
    virtual ~Optimizer() throw();

    Bridge   &bridge;
    DataFile &data;
    numeric<double>::scalar_field F;
    numeric<double>::vector_field G;
    double    dh;


    double    Function( const array_t &p );
    void      Gradient( array_t &g, const array_t &p);
    bool      Callback( const array_t &p );

    void      Run( array_t &p );
    
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Optimizer);
};
