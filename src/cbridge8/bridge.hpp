#ifndef BRIDGE_INCLUDED
#define BRIDGE_INCLUDED 1

#include "lens.hpp"
#include "yocto/ptr/auto.hpp"
#include "yocto/sequence/many-arrays.hpp"

typedef array<double>                      array_t;
typedef many_arrays<double,memory::global> arrays_t;

class Bridge
{
public:
    auto_ptr<Lens> lens;
    const double   capillary_length;
    arrays_t       arrays;
    array_t       &Y;
    array_t       &k1;
    array_t       &k2;
    array_t       &k3;
    array_t       &k4;
    array_t       &V;

    explicit Bridge(const Lens &usr_lens, const double clength);
    virtual ~Bridge() throw();
    

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);

};

#endif

