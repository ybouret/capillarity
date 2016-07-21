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
#include "yocto/string/conv.hpp"

using namespace yocto;
using namespace gfx;

static unit_t y_limit = 0;
static inline bool is_bad_vertex(const vertex &v) throw()
{
    return v.y>=y_limit;
}


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
        y_limit         = h/2;
        xpatches      xps(original, new threading::engine(4,0,true) );
        Filter<float> F(w,h);
        size_t        nf = 0;
        if(argc>2)
        {
            nf = strconv::to<size_t>(argv[2],"#filters");
        }

        //______________________________________________________________________
        //
        //
        // extracting foreground...
        //
        //______________________________________________________________________
        pixmapf  img(original);
        IMG.save("img.png",img,NULL);

        std::cerr << "-- Filtering #" << nf << std::endl;
        for(size_t i=1;i<=nf;++i)
        {
            (std::cerr << ".").flush();
            F.Median(img,xps);
        }
        if(nf>0) std::cerr << std::endl;
        IMG.save("img_f.png",img,NULL);

        std::cerr << "-- Building Edges" << std::endl;
        Edges edges(w,h);
        edges.build_from(img,xps);
        IMG.save("edges.png",edges,NULL);


        std::cerr << "-- Extracting Foreground" << std::endl;
        pixmap1 edgesU(edges,gist::float2byte,edges);
        pixmap1 fg(w,h);
        separate(threshold::keep_foreground,fg,edgesU,xps);

        IMG.save("fg.png",fg,0);

        std::cerr << "-- Extracting Tags" << std::endl;
        tagmap tags(w,h);
        tags.build(fg,8);
        IMG.save("tags.png", tags, tags.colors, NULL);

        std::cerr << "-- Extract Particles and remove Shallows..." << std::endl;
        particles pa;
        pa.load(tags);
        std::cerr << "#particles=" << pa.size() << std::endl;
        pa.remove_shallow_with(tags);
        std::cerr << "#particles=" << pa.size() << std::endl;
        IMG.save("tags_lw.png", tags, tags.colors, NULL);

        std::cerr << "-- Removing Vertices" << std::endl;
        pa.reject_all_vertices_from(tags,is_bad_vertex);
        IMG.save("tags_ok.png",tags,tags.colors,NULL);
        std::cerr << "#particles=" << pa.size() << std::endl;

        if(pa.size()<2)
        {
            throw exception("Must have at least 2 particles...");
        }
        particle &p1 = *pa[1];
        particle &p2 = *pa[2];

        pixmapf wksp(w,h);
        p1.transfer(wksp,original);
        p2.transfer(wksp,original);
        IMG.save("wksp.png",wksp,NULL);


    }


}
YOCTO_PROGRAM_END()
