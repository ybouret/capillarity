#include "yocto/program.hpp"
#include "yocto/graphics/image/png.hpp"
#include "yocto/graphics/image/jpeg.hpp"
#include "yocto/graphics/image/tiff.hpp"
#include "yocto/graphics/ops/gradient.hpp"
#include "yocto/graphics/ops/blobs.hpp"
#include "yocto/graphics/ops/histogram.hpp"
#include "yocto/graphics/ops/blur.hpp"
#include "yocto/graphics/ops/mix.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/graphics/rawpix.hpp"

#include "yocto/math/alg/shapes.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/stat/descr.hpp"

using namespace yocto;
using namespace graphics;

YOCTO_PROGRAM_START()
{
    // image formats
    image               &IMG = image::instance();
    const image::format &PNG = IMG.declare( new png_format()  );
    IMG.declare( new jpeg_format() );
    IMG.declare( new tiff_format() );

    // parallel
    threading::engine server(false);
    xpatches          xps;

    if(argc<=1)
    {
        throw exception("usage: %s filename", program);
    }

    const string filename = argv[1];


    pixmapf           pxm( IMG.loadf(filename,NULL) );
    PNG.save("image.png", pxm, NULL);
    xpatch::create(xps, pxm, &server);

    const unit_t w = pxm.w;
    const unit_t h = pxm.h;


    pixmapf grd(w,h);
    {
        pixmapf tmp(w,h);
        compute_gradient(grd, tmp, pxm, xps, &server);
    }
    PNG.save("image_grad.png", grd, NULL);


}
YOCTO_PROGRAM_END()
