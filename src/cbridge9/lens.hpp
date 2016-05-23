#ifndef LENS_INCLUDED
#define LENS_INCLUDED 1

#include "yocto/math/fcn/drvs.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/ios/imstream.hpp"
#include "yocto/ptr/shared.hpp"

using namespace yocto;
using namespace math;

typedef numeric<double>::function Function;
typedef derivative<double>        Derivative;
typedef shared_ptr<Derivative>    SharedDerivative;


#define BRIDGE_R   1
#define BRIDGE_Z   2
#define BRIDGE_A   3


class Lens : public object
{
public:
    const double     beta;         //!< upper limit of fitted value
    const double     pimb;         //!< pi - beta
    vector<double>   coef;         //!< polynomial coefficient of alpha<2>
    const double     R0;           //!< fitted@0
    const double     R_beta;       //!< fitted@beta
    const double     R_beta_prime; //!< derivative of fitted@beta
    const double     R_pi;         //!< R@pi
    const double     U;            //!< for extend
    const double     V;            //!< for extend
    Function         R_fitted;     //!< 0->beta
    Function         R;            //!< 0->beta->pi
    Function         omega;        //!< tangent angle
    SharedDerivative drvs;         //!< to compute derivatices
    double           alpha_h;      //!< 1e-4

    explicit Lens(const double            user_beta,
                  const array<double>    &user_coef,
                  const SharedDerivative &user_drvs );

    virtual ~Lens() throw();

    double compute_fitted(const double alpha); //!< fitted radius
    double compute_extend(const double alpha); //!< extended radius
    double compute_omega(const double alpha);  //!< tangent angle

    static Lens * load( ios::istream &fp,       const SharedDerivative &user_drvs );
    static Lens * load( const string &filename, const SharedDerivative &user_drvs );

    void starting_point( array<double> &param, const double alpha, const double theta, const double height);

    double find_alpha(const double surf);

private:
    double local_surface;
    double z_surface(const double alpha);
    
    //double __minus_surface(const double alpha);
    YOCTO_DISABLE_COPY_AND_ASSIGN(Lens);
    //void compute_max_surface();
};

#endif
