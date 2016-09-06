#include "yocto/program.hpp"

#include "yocto/gfx/image.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/tiff.hpp"

#include "yocto/gfx/ops/edges.hpp"


using namespace yocto;
using namespace gfx;
using namespace math;

YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(tiff);

    imageIO          &IMG = image::instance();
    //float             sig = 1.4f;
    if(argc>1)
    {
        const string  filename = argv[1];
        std::cerr << "-- Loading " << filename << std::endl;
        const pixmapf source( IMG.loadf(filename,0) );
        std::cerr << "-- Prepare Engine" << std::endl;
        xpatches xps(source,true);

        // hard copy
        pixmapf img(source);
        IMG.save("img.png",img,0);
        const unit_t w = source.w;
        const unit_t h = source.h;

        std::cerr << "Edges..." << std::endl;
        EdgeDetector ED(w,h);
        stencil_scharr_x Gx;
        stencil_scharr_y Gy;

        ED.build_from(img,Gx,Gy,xps);
        IMG.save("img-tags.png", ED.tags, ED.tags.colors, 0);
        std::cerr << "#edges=" << ED.edges.size() << std::endl;
    }

}
YOCTO_PROGRAM_END()

