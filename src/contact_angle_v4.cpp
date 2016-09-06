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
        const pixmap3 source( IMG.load3(filename,0) );
        std::cerr << "-- Prepare Engine" << std::endl;
        xpatches xps(source,true);

        // hard copy
        pixmapf img(source,RGB::to_float,source);
        IMG.save("img.png",img,0);
        const unit_t w = source.w;
        const unit_t h = source.h;

        // detect edges
        std::cerr << "Edges..." << std::endl;
        EdgeDetector     ED(w,h);
        stencil_sobel_x5 Gx;
        stencil_sobel_y5 Gy;

        ED.build_from(img,Gx,Gy,xps);
        IMG.save("img-tags.png", ED.tags, ED.tags.colors, 0);
        std::cerr << "#edges=" << ED.edges.size() << std::endl;

        pixmap3 tgt(source);
        for(size_t i=ED.edges.size();i>0;--i)
        {
            const particle &pa = *ED.edges[i];
            pa.mask(tgt, named_color::fetch(pa.tag+ED.tags.colors.shift), 255);
        }
        IMG.save("img-edges.png", tgt, 0);

    }

}
YOCTO_PROGRAM_END()

