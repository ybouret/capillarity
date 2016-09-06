#include "yocto/program.hpp"

#include "yocto/gfx/image.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/tiff.hpp"

#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/filter.hpp"
#include "yocto/string/conv.hpp"

using namespace yocto;
using namespace gfx;
using namespace math;

YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(tiff);

    imageIO          &IMG = image::instance();
    float             sig = 1.4f;
    if(argc>1)
    {
        const string  filename = argv[1];
        std::cerr << "-- Loading " << filename << std::endl;
        const pixmap3 source( IMG.load3(filename,0) );
        const unit_t w = source.w;
        const unit_t h = source.h;
        std::cerr << "-- Prepare Engine" << std::endl;
        xpatches xps(source,true);

        // hard copy
        const pixmapf img0(source,RGB::to_float,source);
        IMG.save("img0.png",img0,0);

        pixmapf       img(img0);
        if(argc>2)
        {
            sig = strconv::to_float(argv[2],"sig");
        }

        if(sig>0)
        {
            std::cerr << "-- Gaussian Filtering, sigma=" << sig << std::endl;
            const stencil_gauss   g5(2,sig);
            stencil::dispatcher   dsp(w,h);
            dsp(g5,img,img0,xps);
        }
        else
        {
            std::cerr << "-- Using Raw Image..." << std::endl;
        }
        IMG.save("img.png",img,0);


        // detect edges
        std::cerr << "-- Detecting Edges..." << std::endl;
        EdgeDetector     ED(w,h);
        stencil_scharr_x5 Gx;
        stencil_scharr_y5 Gy;

        ED.build_from(img,Gx,Gy,xps);
        ED.tags.colors.shift = YGFX_RED;
        IMG.save("img-grad.png", ED, 0 );
        IMG.save("img-tags.png", ED.tags, ED.tags.colors, 0);
        std::cerr << "#edges=" << ED.edges.size() << std::endl;


        pixmap3 tgt(source);
        for(size_t i=ED.edges.size();i>0;--i)
        {
            const particle &pa = *ED.edges[i];
            pa.mask(tgt, named_color::fetch(pa.tag+ED.tags.colors.shift), 255);
        }
        IMG.save("img-edges.png", tgt, 0);

        std::cerr << "-- Closing..." << std::endl;
        Filter<float> F(w,h);
        F.Close(ED,xps);
        IMG.save("img-close.png", ED, 0);

    }

}
YOCTO_PROGRAM_END()

