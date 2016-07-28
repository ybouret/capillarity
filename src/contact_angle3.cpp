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

        // histogram
        std::cerr << "-- Histogram..." << std::endl;
        H.reset();
        H.update(edges,xps);
        const size_t tt = H.threshold();
        std::cerr << "tt=" << tt << std::endl;

        {
            ios::wcstream fp("hist.dat");
            for(size_t i=0;i<H.bins;++i)
            {
                fp("%g %g\n", double(i), double(H.count[i]));
            }
        }

        threshold::apply(img,tt,edges,threshold::keep_foreground);
        IMG.save("img_edges1.png",img,0);
        threshold::apply(img,tt/2,edges,threshold::keep_foreground);
        IMG.save("img_edges2.png",img,0);


    }

}
YOCTO_PROGRAM_END()