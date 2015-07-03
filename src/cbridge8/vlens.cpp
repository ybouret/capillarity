#include "vlens.hpp"
#include "yocto/exception.hpp"

VLens:: ~VLens() throw()
{
}

VLens:: VLens(const array<double> &usr_params ) :
Lens(),
params(nparam,0.0)
{
    if(usr_params.size()<nparam)
        throw exception("Missing Parameters for VLens");
    for(size_t i=1;i<=nparam;++i)
    {
        params[i] = usr_params[i];
    }
}

double VLens:: ComputeRho(const double alpha) throw()
{
    return RhoFit(alpha, params);
}

Lens * VLens:: clone() const
{
    return new VLens(this->params);
}

#include "yocto/ios/icstream.hpp"
#include "yocto/string/conv.hpp"
static inline bool is_bad(char C) throw()
{
    return C == ' ' || C == '\t';
}

VLens * VLens:: ReadFrom(const string &filename)
{
    ios::icstream fp(filename);
    vector<double> Params(nparam,as_capacity);
    string line;
    for(unsigned i=1;i<=nparam;++i)
    {
        line.clear();
        if(fp.read_line(line) < 0 )
        {
            throw exception("missing param#%u in %s", i, filename.c_str());
        }
        line.clean(is_bad);
        Params.push_back( strconv::to<double>(line,"parameter") );
    }
    return new VLens( Params );
}

