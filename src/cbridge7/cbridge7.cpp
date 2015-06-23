#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/sequence/many-arrays.hpp"
#include "yocto/math/trigconv.hpp"

#include "yocto/code/utils.hpp"
#include "yocto/math/kernel/tao.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/ptr/auto.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/math/fcn/zfind.hpp"

using namespace yocto;
using namespace math;

typedef array<double> array_t;


class Bridge
{
public:
    explicit Bridge()
    {

    }


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bridge);
};

YOCTO_PROGRAM_START()
{

}
YOCTO_PROGRAM_END()
