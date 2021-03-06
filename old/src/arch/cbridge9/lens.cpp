#include "lens.hpp"
#include "yocto/code/ipower.hpp"
#include "yocto/math/core/tao.hpp"
#include "yocto/math/fcn/intg.hpp"

Lens:: ~Lens() throw()
{}


Lens:: Lens(const double            user_beta,
            const array<double>    &user_coef,
            const SharedDerivative &user_drvs) :
beta(user_beta),
pimb(numeric<double>::pi-beta),
coef(user_coef.size()),
R0(0),
R_beta(0),
R_beta_prime(0),
R_pi(0),
U(0),
V(0),
R_fitted(this, &Lens::compute_fitted),
R(       this, &Lens::compute_extend),
omega(   this, &Lens::compute_omega ),
drvs(user_drvs),
alpha_h(1e-4),
local_surface(0)
{
    integrator<double> intg;
    // TODO: check values...
    tao::set(coef,user_coef);
    (double&)R0           = compute_fitted(0);
    (double&)R_beta       = compute_fitted(beta);
    (double&)R_beta_prime = (*drvs)(R_fitted,beta,alpha_h);
    (double&)R_pi         = intg(0,beta,R_fitted,1e-5)/beta;

    std::cerr << "R0          = " << R0           << std::endl;
    std::cerr << "beta        = " << beta         << std::endl;
    std::cerr << "R_beta      = " << R_beta       << std::endl;
    std::cerr << "R_beta_prime= " << R_beta_prime << std::endl;
    std::cerr << "R_pi        = " << R_pi         << std::endl;

    const double dR = R_pi - R_beta;
    const double sp = pimb*R_beta_prime;
    (double &)U = sp - 4.0 * dR;
    (double &)V = 2.0*dR - sp;
    std::cerr << "U=" << U << "; V=" << V << std::endl;
    //compute_max_surface();
}


double Lens:: compute_fitted(const double alpha)
{
    double       rho = 0;
    const double a2  = alpha*alpha;
    for(size_t i=coef.size(),p=i-1;i>0;--i,--p)
    {
        rho += coef[i] * ipower(a2,p);
    }
    return rho;
}

double Lens:: compute_extend(const double alpha)
{
    if(fabs(alpha)<=beta)
    {
        return compute_fitted(alpha);
    }
    else
    {
        const double X = Square( (numeric<double>::pi-fabs(alpha)) / pimb );
        return R_pi + X*(U+V*X)*0.5;
    }
}

double Lens:: compute_omega(const double alpha)
{
    const double dr = (*drvs)(R,alpha,alpha_h);
    const double rr = R(alpha);
    return alpha - asin( dr/Hypotenuse(dr,rr) );
}


#include "yocto/exceptions.hpp"
#include "yocto/string/tokenizer.hpp"
#include "yocto/string/conv.hpp"

static inline bool is_sep(char C) throw()
{
    return C == ' ' || C == '\t' || C == ',';
}

Lens * Lens:: load( ios::istream &fp, const SharedDerivative &user_drvs )
{
    static const char fn[] = "Lens::load";
    string         line;
    vector<string> words(8,as_capacity);

    if( fp.read_line(line) <= 0 ) throw imported::exception(fn,"missing line containing beta");
    words.free();
    tokenizer::split(words, line, is_sep);
    if(words.size()<=0) throw imported::exception(fn,"missing beta value");
    const double __beta =  strconv::to<double>(words[1],"beta");

    line.clear();
    if(fp.read_line(line)<=0) throw imported::exception(fn,"missing line containing parameters");
    words.free();
    tokenizer::split(words, line, is_sep);
    vector<double> __coef( words.size(), as_capacity);
    for(size_t i=1;i<=words.size();++i)
    {
        __coef.push_back( strconv::to<double>(words[i],"coef") );
    }


    return new Lens(__beta,__coef,user_drvs);
}

#include "yocto/ios/icstream.hpp"

Lens * Lens:: load( const string &filename, const SharedDerivative &user_drvs )
{
    ios::icstream fp(filename);
    return load(fp,user_drvs);
}


void Lens:: starting_point(array<double> &param,
                           const double   alpha,
                           const double   theta,
                           const double   height)
{
    const double rr = R(alpha);
    param[BRIDGE_R] = rr * sin(alpha);
    param[BRIDGE_Z] = height + R0 - rr * cos(alpha);
    param[BRIDGE_A] = (omega(alpha)+theta) - numeric<double>::pi;
}



#if 0
#include "yocto/math/opt/bracket.hpp"
#include "yocto/math/opt/minimize.hpp"

void Lens:: compute_max_surface()
{
    Function F(this, & Lens::__minus_surface );
    triplet<double> Angle = { 0, 0, numeric<double>::pi };
    triplet<double> mSurf = { F(Angle.a), F(Angle.b), F(Angle.c) };
    if( ! bracket<double>::inside(F, Angle, mSurf) )
    {
        throw exception("couldn't bracket max surface !");
    }
    std::cerr << "Angle=" << Angle << std::endl;
    std::cerr << "mSurf=" << mSurf << std::endl;
    minimize(F,Angle, mSurf, 1e-5);
    (double &)max_angle   =  Angle.b;
    (double &)max_surface = -mSurf.b;
    std::cerr << "max_angle=" << max_angle << std::endl;
    std::cerr << "max_surface=" << max_surface << std::endl;
    exit(1);
}
#endif

double Lens:: z_surface(const double alpha)
{
    return numeric<double>::pi * Square( R(alpha)*sin(alpha) ) - local_surface;
}

#include "yocto/math/fcn/zfind.hpp"

double Lens:: find_alpha(const double surface)
{
    local_surface = surface;
    Function  F(this, & Lens::z_surface );
    triplet<double> a = { 0,0,0 };
    triplet<double> Z = { F(a.a),0,0};
    bool   found = false;
    for(int a_deg=1;a_deg<=180;++a_deg)
    {
        a.c = Deg2Rad(double(a_deg));
        Z.c = F(a.c);
        if(Z.a*Z.c<=0)
        {
            found = true;
            break;
        }
    }
    if(!found)
    {
        throw exception("couldn't find angle for surface=%g", surface);
    }
    zfind<double> solve(1e-5);
    return solve.run(F,a,Z);
}

double Lens:: find_horizontal_alpha(const double theta)
{
    // Phi=0=omega+theta-pi => omega = pi - theta
    zfunction<double> F(omega,numeric<double>::pi-theta);
    triplet<double>   a = { 0,0,numeric<double>::pi };
    triplet<double>   f = { F.call(a.a), 0, F.call(a.c) };
    if( f.a * f.c > 0 )
    {
        throw exception("couldn't find horizontal alpha for theta=%g", theta);
    }
    zfind<double> solve(1e-5);
    return solve.run(F.call,a,f);
}

double Lens:: find_horizontal_theta(const double alpha)
{
    // Phi=0=omega+theta-pi => theta = pi - omega(alpha)
    return numeric<double>::pi - omega(alpha);
}

Lens * Lens:: sphere(const double radius, const SharedDerivative &user_drvs )
{
    vector<double> user_coef(1);
    user_coef[1] = fabs(radius);
    return new Lens(1,user_coef,user_drvs);
}


