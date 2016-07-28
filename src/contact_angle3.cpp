#include "yocto/program.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/tiff.hpp"
#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/particles.hpp"
#include "yocto/gfx/ops/histogram.hpp"
#include "yocto/gfx/ops/filter.hpp"
#include "yocto/gfx/ops/transform.hpp"
#include "yocto/gfx/ops/blur.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/gfx/rawpix.hpp"

#include "yocto/math/alg/shapes2d.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/stat/descr.hpp"
#include "yocto/string/conv.hpp"

#include "yocto/math/alg/shapes2d.hpp"
#include "yocto/math/opt/cgrad.hpp"

using namespace yocto;
using namespace gfx;
using namespace math;

static inline
float saturate(const float x,const float p)
{
    return powf(x,p);
}

static inline
float invert(const float x)
{
    return 1.0f-x;
}

YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(tiff);

    imageIO          &IMG = image::instance();

    if(argv[1]>0)
    {
        // Load original image
        std::cerr << "-- Loading" << std::endl;
        const string  filename = argv[1];
        const pixmapf original( IMG.loadf(filename,0) );
        const unit_t  w = original.w;
        const unit_t  h = original.h;
        xpatches      xps(original, new threading::engine(4,0,true));
        transform     TR;
        Filter<float> F(w,h);
        Edges         edges(w,h);
        Histogram     H;

        // copy
        pixmapf img(original);
        IMG.save("img.png",img,NULL);
        std::cerr << "-- Edges" << std::endl;
        edges.build_from(img,xps);
        IMG.save("img_edges0.png",edges,NULL);

        std::cerr << "-- Histogram..." << std::endl;
        H.reset();
        H.update(edges,xps);
        const size_t ht = H.threshold();
        std::cerr << "ht=" << ht << std::endl;
        std::cerr << "-- Prepare" << std::endl;
        threshold::apply(img,ht,edges,threshold::keep_foreground);
        IMG.save("img_edges1.png",img,NULL);

        tagmap tags(w,h);
        tags.build(img,8);
        particles pa;
        pa.load(tags);
        pa.remove_shallow_with(tags);
        IMG.save("tags.png",tags,tags.colors,NULL);


    }

}
YOCTO_PROGRAM_END()