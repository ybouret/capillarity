#include "yocto/program.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/tiff.hpp"
#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/particles.hpp"
#include "yocto/gfx/ops/histogram.hpp"
//#include "yocto/gfx/ops/blur.hpp"
#include "yocto/gfx/ops/filter.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/gfx/rawpix.hpp"

#include "yocto/math/alg/shapes.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/stat/descr.hpp"

using namespace yocto;
using namespace gfx;

YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(tiff);

    imageIO          &IMG = image::instance();

    if(argv[1]>0)
    {
        const string  filename = argv[1];
        const pixmapf original( IMG.loadf(filename,0) );
        const unit_t  w = original.w;
        const unit_t  h = original.h;
        Filter<float> F(w,h);

        //______________________________________________________________________
        //
        // extracting foreground...
        //______________________________________________________________________
        


    }


}
YOCTO_PROGRAM_END()
