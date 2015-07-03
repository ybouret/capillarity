#ifndef VLENS_INCLUDED
#define VLENS_INCLUDED 1

#include "lens.hpp"
#include "yocto/sequence/vector.hpp"

class VLens : public Lens
{
public:
    static const size_t nparam = 4;
    vector<double>      params;

    explicit VLens( const array<double> &usr_params );
    virtual ~VLens() throw();
    
    virtual  double ComputeRho(const double alpha) throw();
    virtual  Lens  *clone() const;

    static inline double G(double X,const double b) throw()
    {
        const double c   = 6.0-3*b;
        const double d   = 3*b-8.0;
        const double e   = 3.0-b;
        const double X2  = X*X;
        const double X4  = X2*X2;
        const double X6  = X2*X4;
        const double X8  = X4*X4;
        return b*X2+c*X4+d*X6+e*X8;
    }

    static inline double GS(const double alpha,const double b, const double beta) throw()
    {
        if(Fabs(alpha)>=beta)
        {
            return 1.0;
        }
        else
        {
            return G(alpha/beta,b);
        }
    }

    static inline double RhoFit(const double alpha, const array<double> &Params)
    {
        const double a    = Params[1];
        const double b    = Params[2];
        const double K    = Params[3];
        const double beta = Params[4];
        
        return a+ (K-a)*GS(alpha,b,beta);
    }
    


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(VLens);

};


#endif
