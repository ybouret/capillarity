#ifndef DATAFILE_INCLUDED
#define DATAFILE_INCLUDED 1

#include "yocto/sequence/vector.hpp"
#include "yocto/string.hpp"

using namespace yocto;
//using namespace math;


class DataFile
{
public:
    vector<double> h;
    vector<double> S;
    const size_t   N; //!< number of points
    vector<double> alpha;
    vector<double> theta;
    DataFile( const string &filename);
    ~DataFile() throw();

    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(DataFile);
};

#endif

