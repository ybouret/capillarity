#include "yocto/program.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/tiff.hpp"
#include "yocto/gfx/ops/stencil.hpp"
#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/filter.hpp"

using namespace yocto;
using namespace gfx;

YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(tiff);

    imageIO &IMG = image::instance();
    float    sig = 1.4f;

    if(argc>1)
    {
        //______________________________________________________________________
        //
        // loading original file, in colors
        //______________________________________________________________________
        const string  filename = argv[1];
        const pixmap3 origin( IMG.load3(filename,NULL) );
        IMG.save("img-origin.png",origin,NULL);

        //______________________________________________________________________
        //
        // making a grayscale picture
        //______________________________________________________________________
        pixmapf      img0(origin,RGB::to_float,origin);
        const unit_t h = origin.h;
        const unit_t w = origin.w;
        xpatches     xps(origin,false);
        std::cerr << "parallel=" << xps.size() << std::endl;
        pixmapf img(w,h);

        //______________________________________________________________________
        //
        // image, raw or filtered if sigma>0
        //______________________________________________________________________
        if(sig>0)
        {
            std::cerr << "-- Gauss " << sig << std::endl;
            stencil_gauss g5(2,sig);
            g5.apply(img,img0,xps);
        }
        else
        {
            img.copy(img0);
        }
        IMG.save("img-filter.png",img,NULL);

        //______________________________________________________________________
        //
        // edge detections
        //______________________________________________________________________
        std::cerr << "-- Edges Detection" << std::endl;
        EdgeDetector      ED(w,h);
        tagmap           &tags  = ED.tags;
        particles        &edges = ED.edges;
        stencil_scharr_x5 gx;
        stencil_scharr_y5 gy;
        ED.build_from(img,gx,gy,xps);
        IMG.save("img-grad.png", ED, 0 );
        IMG.save("img-tags.png",tags,tags.colors,NULL);
        std::cerr << "#edges=" << edges.size() << std::endl;

        //______________________________________________________________________
        //
        // copy edges over original picture
        //______________________________________________________________________
        pixmap3 tgt(origin);
        for(size_t i=edges.size();i>0;--i)
        {
            const particle &pa = *edges[i];
            pa.mask(tgt, named_color::fetch(pa.tag+tags.colors.shift), 255);
        }
        IMG.save("img-edges.png", tgt, 0);

        //______________________________________________________________________
        //
        // close shapes
        //______________________________________________________________________
        std::cerr << "-- Closing..." << std::endl;
        Filter<float> F(w,h);
        F.Close(ED,xps);
        IMG.save("img-close.png", ED, 0);

        //______________________________________________________________________
        //
        // detect final blobs
        //______________________________________________________________________
        std::cerr << "-- Final Blobs..." << std::endl;
        tags.build(ED,8);
        edges.load(tags);
        IMG.save("img-blobs.png", tags, tags.colors, 0);
        std::cerr << "#edges=" << edges.size() << std::endl;

        while( edges.size() > 2 ) edges.pop_back();
        if(edges.size()<=0)
            throw exception("not enough edges");

        tgt.copy(origin);
        for(size_t i=edges.size();i>0;--i)
        {
            const particle &pa = *edges[i];
            pa.mask(tgt, named_color::fetch(pa.tag+tags.colors.shift), 255);
            //std::cerr << "#points=" << pa.size << std::endl;
        }
        IMG.save("img-final.png", tgt, 0);
    }

}
YOCTO_PROGRAM_END()
