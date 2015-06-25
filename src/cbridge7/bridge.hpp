#include "lens.hpp"
#include "yocto/sequence/many-arrays.hpp"

typedef array<double>                      array_t;
typedef many_arrays<double,memory::global> arrays_t;


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

    explicit Bridge(const Lens::Pointer &usr_lens, const double usr_clength);
    virtual ~Bridge() throw();

    //! Y=[r drdz]
    void Evaluate( array_t &dYdz, double z, const array_t &Y ) throw();

    //! make a step
    void RK4(const double z_ini,const double z_end) throw();

    //! warning, angles are in degrees
    bool FinalRadius(const double height,
                     const double theta,
                     const double alpha,
                     const bool   save=false);
    

    //! warning, theta is in degree
    double FindAlpha(const double height,const double theta);

    //! theta in degrees
    double FindHmax(const double theta);

    //! scanning alpha for algo debug
    void ScanAlpha(const double height, const double theta);

    //! warning, alpha is in degree, return as well
    double FindTheta(const double height, const double alpha);


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
};