#ifndef DATAFILE_INCLUDED
#define DATAFILE_INCLUDED 1

#include "yocto/sequence/vector.hpp"
#include "yocto/string.hpp"

using namespace yocto;
//using namespace math;


// usual speed: 5 microns/s


class DataFile
{
public:
    vector<double> h;
    vector<double> S;
    const size_t   N; //!< number of points
    vector<double> t; //!< time, deduce from speed
    vector<double> alpha;
    vector<double> theta;

    //! build with t = start + h/speed or t=start+(h_end-h)/speed
    /**
     \param filename datafile
     \param speed    speed in mm/s
     \param start    start time in s
     */
    DataFile(const string &filename,
             const double speed,
             const double start);
    ~DataFile() throw();

    void Save(const string &filename) const;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(DataFile);
};

#endif

